from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

CANONICAL_AA = "ACDEFGHIKLMNPQRSTVWY"
AA_BACKGROUND = {
    "X": 0.0,
    "S": 0.08621733200464735,
    "L": 0.07651555159603048,
    "P": 0.07203973199729977,
    "T": 0.06922675985837555,
    "A": 0.06405039516162363,
    "G": 0.0596202682782342,
    "V": 0.06135659146120252,
    "R": 0.05800586887338755,
    "D": 0.05276268719088534,
    "N": 0.051973971225069915,
    "Y": 0.0454190182724939,
    "H": 0.05010516212877539,
    "F": 0.04232959372517324,
    "I": 0.0396433672087032,
    "M": 0.033543274904825976,
    "K": 0.03183794929256655,
    "Q": 0.03255870021445727,
    "E": 0.03274491525034556,
    "W": 0.030877254212225442,
    "C": 0.009171607143677185,
}

CHI_SQUARE_THRESHOLD = 3.841
DEFAULT_MAX_COMBINATIONS = 1_000_000


@dataclass
class MannWhitneyResult:
    k: int
    total_kmers: int
    ad_elevated: int
    nc_elevated: int
    result_file: Path
    ad_file: Path
    nc_file: Path
    matrix_file: Path
    volcano_file: Path


@dataclass(frozen=True)
class CohortTilingResult:
    """Per-cohort counts needed for both matrix output and statistical testing."""
    raw_filtered_counts: Dict[str, Dict[str, int]]
    passed_kmers: set[str]
    active_patients: list[str]
    universe: Tuple[str, ...]
    totals: Dict[str, int]


ProgressCallback = Callable[[int, str], None]


def _group_label_for_filename(label: str, fallback: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in label.strip())
    return cleaned or fallback


def _with_optional_suffix(base: str, *suffixes: str) -> str:
    parts = [base]
    parts.extend(suffix.strip("_") for suffix in suffixes if suffix and suffix.strip("_"))
    return "_".join(parts)


def prebuild_kmers(
    kmer_length: int,
    wildcard_positions: Sequence[int] | None = None,
    alphabet: str = CANONICAL_AA,
    max_combinations: int = DEFAULT_MAX_COMBINATIONS,
) -> Tuple[str, ...]:
    """Return every ordered k-mer combination, fixing wildcard positions to X."""
    if kmer_length <= 0:
        raise ValueError("kmer_length must be positive.")
    if not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError("alphabet must contain unique amino-acid symbols.")

    wildcard_set = set(wildcard_positions or [])
    invalid = sorted(position for position in wildcard_set if not 0 <= position < kmer_length)
    if invalid:
        raise ValueError(f"Wildcard positions outside the k-mer: {invalid}")

    variable_count = kmer_length - len(wildcard_set)
    combination_count = len(alphabet) ** variable_count
    if combination_count > max_combinations:
        raise ValueError(
            f"Prebuilding {combination_count:,} k-mers exceeds the configured limit "
            f"of {max_combinations:,}. Reduce k, add wildcards, or explicitly raise "
            "max_combinations if enough memory is available."
        )

    kmers = []
    for residues in product(alphabet, repeat=variable_count):
        residue_iterator = iter(residues)
        kmer = [
            "X" if position in wildcard_set else next(residue_iterator)
            for position in range(kmer_length)
        ]
        kmers.append("".join(kmer))
    return tuple(kmers)


def calculate_expected_frequency(
    kmer: str,
    total_count: int,
    aa_background: Dict[str, float] = AA_BACKGROUND,
) -> float:
    """Calculate E = T * product(p(aa)), ignoring generated X wildcards."""
    expected = float(total_count)
    for amino_acid in kmer:
        if amino_acid == "X":
            continue
        frequency = aa_background.get(amino_acid)
        if frequency is None:
            return 0.0
        expected *= frequency
    return expected


def split_input_by_group(
    input_path: Path,
    positive_keyword: str = "AD",
    negative_keyword: str = "NC"
) -> Tuple[Path, Path]:
    """Split paired sequence-count columns into cohort files based on keywords.

    Matches columns starting with positive_keyword and negative_keyword.
    """
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()
    if suffix == ".xlsx":
        df = pd.read_excel(input_path)
    elif suffix == ".csv":
        df = pd.read_csv(input_path)
    else:
        raise ValueError("Unsupported file format; use .csv or .xlsx.")

    positive_keyword = positive_keyword.strip()
    negative_keyword = negative_keyword.strip()
    if not positive_keyword or not negative_keyword:
        raise ValueError("Positive and negative keywords are required.")

    positive_columns = []
    negative_columns = []

    def matches_keyword(header: str, keyword: str) -> bool:
        return header.lower().startswith(keyword.lower())

    i = 0
    while i < df.shape[1] - 1:
        header = str(df.columns[i]).strip()
        if not header or header.startswith("Unnamed:"):
            i += 1
            continue
        if matches_keyword(header, positive_keyword):
            positive_columns.extend(df.columns[i : i + 2])
            i += 2
        elif matches_keyword(header, negative_keyword):
            negative_columns.extend(df.columns[i : i + 2])
            i += 2
        else:
            i += 1

    if not positive_columns:
        raise ValueError(f"No columns matched the positive keyword '{positive_keyword}'.")
    if not negative_columns:
        raise ValueError(f"No columns matched the negative keyword '{negative_keyword}'.")

    df_pos = df[positive_columns].copy()
    df_neg = df[negative_columns].copy()

    for subset in (df_pos, df_neg):
        for idx in range(1, subset.shape[1], 2):
            subset.columns.values[idx] = "count"

    pos_label = _group_label_for_filename(positive_keyword, "positive")
    neg_label = _group_label_for_filename(negative_keyword, "negative")
    output_pos = input_path.parent / f"{input_path.stem}_{pos_label}{input_path.suffix}"
    output_neg = input_path.parent / f"{input_path.stem}_{neg_label}{input_path.suffix}"

    if suffix == ".xlsx":
        df_pos.to_excel(output_pos, index=False)
        df_neg.to_excel(output_neg, index=False)
    else:
        df_pos.to_csv(output_pos, index=False)
        df_neg.to_csv(output_neg, index=False)

    return output_pos, output_neg


def apply_chi_square_filter(
    raw_counts: Dict[str, int],
    universe: Iterable[str],
    total_count: int,
    aa_background: Dict[str, float] = AA_BACKGROUND,
    threshold: float = CHI_SQUARE_THRESHOLD,
) -> Tuple[Dict[str, int], set[str]]:
    """Apply per-patient chi-square calculation.

    Every universe row is evaluated. Passing zero rows are included in
    ``passed_kmers`` but sparse output omits zeros (matrix restores them).
    """
    filtered_raw: Dict[str, int] = {}
    passed_kmers: set[str] = set()
    if total_count <= 0:
        return filtered_raw, passed_kmers

    for kmer in universe:
        observed = raw_counts.get(kmer, 0)
        expected = calculate_expected_frequency(kmer, total_count, aa_background)
        if expected <= 0:
            continue
        chi_square = (observed - expected) ** 2 / expected
        if chi_square > threshold:
            passed_kmers.add(kmer)
            if observed:
                filtered_raw[kmer] = observed
    return filtered_raw, passed_kmers


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        return pd.read_excel(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError("Unsupported file format; use .csv or .xlsx.")


def _legacy_patient_name(column_name: object) -> str:
    """Match the shortened patient identifiers from column names."""
    stem = Path(str(column_name)).stem
    parts = stem.split("_")
    return "_".join(parts[:2]) if len(parts) > 1 else stem


def _iter_patient_raw_counts(
    path: Path,
    kmer_length: int,
    wildcard_positions: Sequence[int],
    universe_set: set[str],
    progress_callback: Optional[ProgressCallback] = None,
) -> Iterable[Tuple[str, Dict[str, int], int]]:
    """Tile samples while counting windows correctly."""
    data = _read_table(path)
    if data.shape[1] % 2:
        raise ValueError("Input must contain sequence/count column pairs.")

    patient_names: set[str] = set()
    total_patients = data.shape[1] // 2
    for col_idx, column in enumerate(range(0, data.shape[1], 2)):
        sequence_column = data.columns[column]
        count_column = data.columns[column + 1]
        patient = _legacy_patient_name(sequence_column)
        if patient in patient_names:
            raise ValueError(
                f"Multiple sample columns resolve to patient name {patient!r}. "
                "Remove longitudinal duplicates before running analysis."
            )
        patient_names.add(patient)
        counts: Dict[str, int] = defaultdict(int)
        total = 0

        for sequence, count in zip(data[sequence_column], data[count_column]):
            if pd.isna(sequence) or pd.isna(count) or not isinstance(sequence, str):
                continue
            try:
                count = int(count)
            except (TypeError, ValueError):
                continue

            window_count = max(0, len(sequence) - kmer_length + 1)
            total += window_count * count
            for start in range(window_count):
                kmer = list(sequence[start : start + kmer_length])
                for position in wildcard_positions:
                    kmer[position] = "X"
                kmer_string = "".join(kmer)
                if kmer_string in universe_set:
                    counts[kmer_string] += count

        if progress_callback:
            progress_callback(15 + int((col_idx / total_patients) * 15), f"Tiling patients ({col_idx + 1}/{total_patients})")

        yield patient, dict(counts), total


def tile_patient_file(
    path: Path,
    kmer_length: int = 4,
    wildcard_positions: Sequence[int] | None = None,
    aa_background: Dict[str, float] = AA_BACKGROUND,
    chi_square_threshold: float = CHI_SQUARE_THRESHOLD,
    alphabet: str = CANONICAL_AA,
    max_combinations: int = DEFAULT_MAX_COMBINATIONS,
    progress_callback: Optional[ProgressCallback] = None,
) -> CohortTilingResult:
    """Tile and product-filter every sample independently."""
    path = Path(path)
    wildcard_positions = sorted(set(wildcard_positions or []))
    universe = prebuild_kmers(
        kmer_length,
        wildcard_positions=wildcard_positions,
        alphabet=alphabet,
        max_combinations=max_combinations,
    )
    raw_filtered_counts: Dict[str, Dict[str, int]] = {}
    active_patients: list[str] = []
    totals: Dict[str, int] = {}
    passed_kmers: set[str] = set()

    for patient, counts, total in _iter_patient_raw_counts(
        path,
        kmer_length,
        wildcard_positions,
        set(universe),
        progress_callback=progress_callback,
    ):
        filtered_raw, patient_passed = apply_chi_square_filter(
            counts,
            universe,
            total,
            aa_background=aa_background,
            threshold=chi_square_threshold,
        )
        raw_filtered_counts[patient] = filtered_raw
        totals[patient] = total
        if patient_passed:
            active_patients.append(patient)
        passed_kmers.update(patient_passed)

    return CohortTilingResult(
        raw_filtered_counts=raw_filtered_counts,
        passed_kmers=passed_kmers,
        active_patients=active_patients,
        universe=universe,
        totals=totals,
    )


def build_kmer_matrix(
    patient_counts: Dict[str, Dict[str, int | float]],
    kmers: Iterable[str],
    output_path: Path | None = None,
    normalize: bool = False,
    normalization_totals: Dict[str, int] | None = None,
) -> pd.DataFrame:
    """Build a zero-filled matrix, optionally using pre-filter sample totals."""
    index = pd.Index(tuple(kmers), name="kmer")
    matrix = pd.DataFrame.from_dict(patient_counts).reindex(index).fillna(0.0)
    if normalize:
        if normalization_totals is None:
            raise ValueError("normalization_totals are required when normalize=True.")
        for patient in matrix.columns:
            total = normalization_totals.get(str(patient), 0)
            if total:
                matrix[patient] /= total
    if output_path is not None:
        matrix.to_csv(output_path)
    return matrix


def _filter_kmers_by_zero_percentage(
    matrix: pd.DataFrame,
    kmers: Sequence[str],
    max_zero_percentage: float,
) -> list[str]:
    """Return output k-mers at or below the allowed percentage of zero cells."""
    if not 0 <= max_zero_percentage <= 100:
        raise ValueError("max_zero_percentage must be between 0 and 100.")
    zero_percentages = matrix.eq(0).mean(axis=1).mul(100)
    return [
        kmer
        for kmer in kmers
        if float(zero_percentages.loc[kmer]) <= max_zero_percentage
    ]


def _write_volcano_plot(
    results: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a self-contained interactive HTML volcano plot."""
    rows = [
        [
            str(row.kmer),
            float(row.p_value),
            float(row.mean_rank_diff),
            float(row[3]),
        ]
        for row in results.itertuples(index=False)
    ]
    point_data = json.dumps(rows, separators=(",", ":"))
    html = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Interactive MWU volcano plot</title>
<style>
:root{font-family:Inter,system-ui,sans-serif;color:#17212b;background:#f4f7f9}*{box-sizing:border-box}
body{margin:0;padding:22px}.shell{max-width:1500px;margin:auto}.top{display:flex;gap:14px;align-items:end;flex-wrap:wrap;margin-bottom:14px}
h1{font-size:24px;margin:0 auto 0 0}.control{display:grid;gap:4px;font-size:12px;color:#52606d}.control input{height:36px;border:1px solid #b0bec5;border-radius:7px;padding:0 10px;background:white}.search-result{width:100%;min-height:18px;margin:-7px 0 7px;color:#455a64;font-size:12px}.search-result.error{color:#b71c1c}
button{height:36px;border:0;border-radius:7px;padding:0 14px;background:#263238;color:white;cursor:pointer}.stats{display:grid;grid-template-columns:repeat(4,minmax(120px,1fr));gap:10px;margin-bottom:12px}
.stat{background:white;border:1px solid #dce4e8;border-radius:9px;padding:10px 14px}.stat b{font-size:20px;display:block}.stat span{font-size:12px;color:#607d8b}
.plot{position:relative;height:min(72vh,820px);min-height:520px;background:white;border:1px solid #dce4e8;border-radius:10px;overflow:hidden}canvas{width:100%;height:100%;display:block;cursor:crosshair}
.tip{position:absolute;display:none;pointer-events:none;background:#17212bee;color:white;padding:9px 11px;border-radius:7px;font-size:12px;line-height:1.5;box-shadow:0 4px 16px #0004}
.help{margin:9px 3px;color:#607d8b;font-size:12px}.legend{display:flex;gap:15px;align-items:center}.dot{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:5px}.red{background:#d32f2f}.blue{background:#1976d2}.gray{background:#90a4ae}
@media(max-width:700px){.stats{grid-template-columns:1fr 1fr}}
</style></head><body><div class="shell">
<div class="top"><h1>Interactive MWU volcano plot</h1>
<label class="control">Search k-mer<input id="search" placeholder="e.g. ACDG" autocomplete="off"></label><button id="find">Find k-mer</button>
<label class="control">Q-value cutoff<input id="cutoff" type="number" min="0" max="1" step="0.001" value="0.05"></label>
<button id="reset">Reset view</button></div><div class="search-result" id="searchResult" aria-live="polite">Enter a complete k-mer to locate and mark it.</div>
<div class="stats"><div class="stat"><b id="total"></b><span>Tested k-mers</span></div><div class="stat"><b id="sig"></b><span>Q-significant</span></div><div class="stat"><b id="pos"></b><span>Positive significant</span></div><div class="stat"><b id="neg"></b><span>Negative significant</span></div></div>
<div class="plot"><canvas id="canvas"></canvas><div class="tip" id="tip"></div></div>
<div class="help"><span class="legend"><span><i class="dot red"></i>Positive significant</span><span><i class="dot blue"></i>Negative significant</span><span><i class="dot gray"></i>Not significant</span></span>Scroll to zoom · drag to pan · double-click or Reset view to restore · hover for k-mer details. Y-axis is −log10(p-value); significance uses the selected Q-value cutoff.</div>
</div><script>
const raw=__POINT_DATA__;
const points=raw.map(r=>({k:r[0],p:Math.max(r[1],Number.MIN_VALUE),x:r[2],q:r[3],y:-Math.log10(Math.max(r[1],Number.MIN_VALUE))}));
const canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d'),tip=document.getElementById('tip'),cutoffEl=document.getElementById('cutoff'),searchEl=document.getElementById('search'),searchResult=document.getElementById('searchResult');
const pad={l:72,r:24,t:24,b:58};let screenX=new Float32Array(points.length),screenY=new Float32Array(points.length),drag=null,hoverFrame=0,selected=-1;
const fullX=Math.max(1,points.reduce((m,p)=>Math.max(m,Math.abs(p.x)),0)),fullY=Math.max(1.5,points.reduce((m,p)=>Math.max(m,p.y),0))*1.05;let view={xmin:-fullX,xmax:fullX,ymin:0,ymax:fullY};
function cutoff(){const v=Number(cutoffEl.value);return Number.isFinite(v)?Math.min(1,Math.max(0,v)):0.05}
function resize(){const dpr=devicePixelRatio||1,r=canvas.getBoundingClientRect();canvas.width=Math.round(r.width*dpr);canvas.height=Math.round(r.height*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);draw()}
function dims(){const r=canvas.getBoundingClientRect();return{w:r.width,h:r.height,pw:r.width-pad.l-pad.r,ph:r.height-pad.t-pad.b}}
function sx(x,d){return pad.l+(x-view.xmin)/(view.xmax-view.xmin)*d.pw}function sy(y,d){return pad.t+d.ph-(y-view.ymin)/(view.ymax-view.ymin)*d.ph}
function fmt(v){return v<0.001?v.toExponential(3):v.toFixed(5)}
function draw(){const d=dims(),qcut=cutoff();ctx.clearRect(0,0,d.w,d.h);ctx.fillStyle='#fff';ctx.fillRect(0,0,d.w,d.h);ctx.font='12px system-ui';ctx.strokeStyle='#dce4e8';ctx.fillStyle='#52606d';ctx.lineWidth=1;
for(let i=0;i<=5;i++){const y=view.ymin+(view.ymax-view.ymin)*i/5,Y=sy(y,d);ctx.beginPath();ctx.moveTo(pad.l,Y);ctx.lineTo(d.w-pad.r,Y);ctx.stroke();ctx.fillText(y.toFixed(2),8,Y+4)}
for(let i=0;i<=6;i++){const x=view.xmin+(view.xmax-view.xmin)*i/6,X=sx(x,d);ctx.beginPath();ctx.moveTo(X,pad.t);ctx.lineTo(X,d.h-pad.b);ctx.stroke();ctx.textAlign='center';ctx.fillText(x.toFixed(2),X,d.h-pad.b+20)}ctx.textAlign='start';
const lineY=sy(-Math.log10(.05),d);if(lineY>=pad.t&&lineY<=d.h-pad.b){ctx.strokeStyle='#ef6c00';ctx.setLineDash([6,5]);ctx.beginPath();ctx.moveTo(pad.l,lineY);ctx.lineTo(d.w-pad.r,lineY);ctx.stroke();ctx.setLineDash([])}
let sig=0,pos=0,neg=0;for(let i=0;i<points.length;i++){const p=points[i],X=sx(p.x,d),Y=sy(p.y,d);screenX[i]=X;screenY[i]=Y;const s=p.q<=qcut;if(s){sig++;p.x>0?pos++:p.x<0&&neg++}if(X<pad.l||X>d.w-pad.r||Y<pad.t||Y>d.h-pad.b)continue;ctx.fillStyle=s?(p.x>0?'#d32f2f':p.x<0?'#1976d2':'#7b1fa2'):'#90a4ae';ctx.globalAlpha=s?.82:.42;ctx.fillRect(X-1.5,Y-1.5,3,3)}ctx.globalAlpha=1;
if(selected>=0){const p=points[selected],X=screenX[selected],Y=screenY[selected];if(X>=pad.l&&X<=d.w-pad.r&&Y>=pad.t&&Y<=d.h-pad.b){ctx.strokeStyle='#ffb300';ctx.fillStyle='#fff8e1';ctx.lineWidth=3;ctx.beginPath();ctx.arc(X,Y,9,0,Math.PI*2);ctx.fill();ctx.stroke();ctx.strokeStyle='#3e2723';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(X-13,Y);ctx.lineTo(X+13,Y);ctx.moveTo(X,Y-13);ctx.lineTo(X,Y+13);ctx.stroke();ctx.font='bold 13px system-ui';const tw=ctx.measureText(p.k).width;const lx=Math.min(d.w-pad.r-tw-12,Math.max(pad.l,X+12)),ly=Math.max(pad.t+19,Y-12);ctx.fillStyle='#fff8e1';ctx.strokeStyle='#ffb300';ctx.lineWidth=1;ctx.fillRect(lx,ly-17,tw+10,22);ctx.strokeRect(lx,ly-17,tw+10,22);ctx.fillStyle='#3e2723';ctx.textAlign='start';ctx.fillText(p.k,lx+5,ly)}}
ctx.fillStyle='#263238';ctx.font='14px system-ui';ctx.textAlign='center';ctx.fillText('Mean rank difference',pad.l+d.pw/2,d.h-12);ctx.save();ctx.translate(17,pad.t+d.ph/2);ctx.rotate(-Math.PI/2);ctx.fillText('−log10(p-value)',0,0);ctx.restore();
document.getElementById('total').textContent=points.length.toLocaleString();document.getElementById('sig').textContent=sig.toLocaleString();document.getElementById('pos').textContent=pos.toLocaleString();document.getElementById('neg').textContent=neg.toLocaleString()}
function reset(){view={xmin:-fullX,xmax:fullX,ymin:0,ymax:fullY};tip.style.display='none';draw()}
canvas.addEventListener('wheel',e=>{e.preventDefault();const d=dims(),r=canvas.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top,fx=(mx-pad.l)/d.pw,fy=1-(my-pad.t)/d.ph,z=e.deltaY>0?1.18:.84,cx=view.xmin+fx*(view.xmax-view.xmin),cy=view.ymin+fy*(view.ymax-view.ymin),nw=(view.xmax-view.xmin)*z,nh=(view.ymax-view.ymin)*z;view={xmin:cx-fx*nw,xmax:cx+(1-fx)*nw,ymin:Math.max(0,cy-fy*nh),ymax:Math.max(0,cy-fy*nh)+nh};draw()},{passive:false});
canvas.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,v:{...view}};canvas.setPointerCapture(e.pointerId)});canvas.addEventListener('pointermove',e=>{if(drag){const d=dims(),dx=(e.clientX-drag.x)/d.pw*(drag.v.xmax-drag.v.xmin),dy=(e.clientY-drag.y)/d.ph*(drag.v.ymax-drag.v.ymin);view={xmin:drag.v.xmin-dx,xmax:drag.v.xmax-dx,ymin:Math.max(0,drag.v.ymin+dy),ymax:Math.max(0,drag.v.ymin+dy)+(drag.v.ymax-drag.v.ymin)};draw();return}cancelAnimationFrame(hoverFrame);hoverFrame=requestAnimationFrame(()=>showTip(e))});canvas.addEventListener('pointerup',()=>drag=null);canvas.addEventListener('pointercancel',()=>drag=null);canvas.addEventListener('dblclick',reset);
function showTip(e){const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;let best=-1,dist=64;for(let i=0;i<points.length;i++){const dx=screenX[i]-x,dy=screenY[i]-y,dd=dx*dx+dy*dy;if(dd<dist){dist=dd;best=i}}if(best<0){tip.style.display='none';return}const p=points[best];tip.innerHTML=`<b>${p.k}</b><br>p-value: ${fmt(p.p)}<br>Q value: ${fmt(p.q)}<br>Mean rank difference: ${p.x.toFixed(4)}`;tip.style.display='block';tip.style.left=Math.min(r.width-190,x+14)+'px';tip.style.top=Math.max(8,y-72)+'px'}
function locateKmer(allowPartial=false){const key=searchEl.value.trim().toUpperCase();searchResult.classList.remove('error');if(!key){selected=-1;searchResult.textContent='Enter a complete k-mer to locate and mark it.';draw();return}let matches=[];for(let i=0;i<points.length;i++){if(points[i].k===key||(allowPartial&&points[i].k.includes(key)))matches.push(i)}if(matches.length!==1){selected=-1;searchResult.classList.add('error');searchResult.textContent=matches.length?`${matches.length.toLocaleString()} partial matches. Enter the complete k-mer.`:`No tested k-mer matches “${key}”.`;draw();return}selected=matches[0];const p=points[selected],xr=Math.max(.15,fullX*.08),yr=Math.max(.2,fullY*.08),bottom=Math.max(0,p.y-yr);view={xmin:p.x-xr,xmax:p.x+xr,ymin:bottom,ymax:bottom+2*yr};searchResult.textContent=`Located ${p.k} — p=${fmt(p.p)}, Q=${fmt(p.q)}, mean rank difference=${p.x.toFixed(4)}`;draw()}
searchEl.addEventListener('input',()=>{const key=searchEl.value.trim().toUpperCase();if(points.some(p=>p.k===key))locateKmer();else{selected=-1;searchResult.classList.remove('error');searchResult.textContent=key?'Keep typing, then press Enter or Find k-mer.':'Enter a complete k-mer to locate and mark it.';draw()}});searchEl.addEventListener('keydown',e=>{if(e.key==='Enter'){e.preventDefault();locateKmer(true)}});document.getElementById('find').addEventListener('click',()=>locateKmer(true));cutoffEl.addEventListener('input',draw);document.getElementById('reset').addEventListener('click',reset);new ResizeObserver(resize).observe(canvas.parentElement);resize();
</script></body></html>""".replace("__POINT_DATA__", point_data)
    output_path.write_text(html, encoding="utf-8")


def run_mannwhitney(
    positive: CohortTilingResult,
    negative: CohortTilingResult,
    output_prefix: Path,
    matrix_file: Path,
    positive_label: str = "AD",
    negative_label: str = "NC",
    tested_kmers: Sequence[str] | None = None,
    progress_callback: Optional[ProgressCallback] = None,
) -> MannWhitneyResult:
    """Test the union of passing k-mers using normalized filtered values.

    The directional AD/NC files apply significance cutoff to raw p-value.
    """
    tested_kmers = list(
        tested_kmers
        if tested_kmers is not None
        else sorted(positive.passed_kmers | negative.passed_kmers)
    )
    if not tested_kmers:
        raise ValueError("No k-mers remain for Mann-Whitney testing.")

    positive_patients = positive.active_patients
    negative_patients = negative.active_patients
    if not positive_patients or not negative_patients:
        raise ValueError("Both cohorts need at least one patient with a passing k-mer.")

    # Build matrices with raw filtered counts
    positive_matrix = build_kmer_matrix(
        {patient: positive.raw_filtered_counts[patient] for patient in positive_patients},
        tested_kmers,
    )
    negative_matrix = build_kmer_matrix(
        {patient: negative.raw_filtered_counts[patient] for patient in negative_patients},
        tested_kmers,
    )

    # Normalize for Mann-Whitney testing
    for patient in positive_patients:
        positive_matrix[patient] /= positive.totals[patient]
    for patient in negative_patients:
        negative_matrix[patient] /= negative.totals[patient]

    p_values = []
    mean_rank_differences = []

    for index, kmer in enumerate(tqdm(tested_kmers, desc="Running Mann-Whitney U tests"), start=1):
        positive_values = positive_matrix.loc[kmer].to_numpy()
        negative_values = negative_matrix.loc[kmer].to_numpy()

        combined = np.concatenate([positive_values, negative_values])
        ranks = rankdata(combined)
        mean_rank_diff = float(ranks[: len(positive_values)].mean() - ranks[len(positive_values) :].mean())

        try:
            _, p_value = mannwhitneyu(positive_values, negative_values, alternative="two-sided")
        except ValueError:
            p_value = 1.0

        p_values.append(p_value)
        mean_rank_differences.append(mean_rank_diff)

        total = len(tested_kmers)
        if progress_callback and (index == total or index % max(1, total // 20) == 0):
            progress_callback(60 + int((index / total) * 30), f"Running Mann-Whitney tests ({index}/{total})")

    fdr = multipletests(p_values, method="fdr_tsbky")[1]

    result_df = pd.DataFrame(
        {
            "kmer": tested_kmers,
            "p_value": p_values,
            "mean_rank_diff": mean_rank_differences,
            "Q value": fdr,
        }
    ).sort_values("p_value")

    result_file = output_prefix.with_suffix(".csv")
    volcano_file = output_prefix.with_name(f"{output_prefix.stem}_volcano.html")
    positive_file_label = _group_label_for_filename(positive_label, "positive")
    negative_file_label = _group_label_for_filename(negative_label, "negative")
    ad_file = output_prefix.with_name(f"{output_prefix.stem}_{positive_file_label}.csv")
    nc_file = output_prefix.with_name(f"{output_prefix.stem}_{negative_file_label}.csv")

    result_df.to_csv(result_file, index=False)
    _write_volcano_plot(result_df, volcano_file)
    result_df[result_df["mean_rank_diff"] > 0].to_csv(ad_file, index=False)
    result_df[result_df["mean_rank_diff"] < 0].to_csv(nc_file, index=False)

    return MannWhitneyResult(
        k=len(tested_kmers[0]) if tested_kmers else 0,
        total_kmers=len(result_df),
        ad_elevated=int((result_df["mean_rank_diff"] > 0).sum()),
        nc_elevated=int((result_df["mean_rank_diff"] < 0).sum()),
        result_file=result_file,
        ad_file=ad_file,
        nc_file=nc_file,
        matrix_file=matrix_file,
        volcano_file=volcano_file,
    )


def analyze_single_k(
    input_path: Path,
    pos_file: Path,
    neg_file: Path,
    k: int,
    wildcard_positions: Sequence[int] | None = None,
    normalize: bool = True,
    max_zero_percentage: float = 100.0,
    workdir: Optional[Path] = None,
    positive_label: str = "AD",
    negative_label: str = "NC",
    progress_callback: Optional[ProgressCallback] = None,
) -> MannWhitneyResult:
    """Analyze a single k value for both cohorts.

    This is the main entry point for analysis.
    """
    workdir = Path(workdir or input_path.parent)
    workdir.mkdir(parents=True, exist_ok=True)

    wildcard_positions = list(wildcard_positions or [])
    if not 0 <= max_zero_percentage <= 100:
        raise ValueError("max_zero_percentage must be between 0 and 100.")

    if progress_callback:
        progress_callback(10, "Tiling positive cohort")
    positive = tile_patient_file(
        pos_file,
        kmer_length=k,
        wildcard_positions=wildcard_positions,
        progress_callback=progress_callback,
    )

    if progress_callback:
        progress_callback(30, "Tiling negative cohort")
    negative = tile_patient_file(
        neg_file,
        kmer_length=k,
        wildcard_positions=wildcard_positions,
        progress_callback=progress_callback,
    )

    if progress_callback:
        progress_callback(50, "Building k-mer matrix")

    # Use union of passing kmers from both cohorts
    tested_kmers = sorted(positive.passed_kmers | negative.passed_kmers)

    # Build output matrix with all tested kmers
    wildcard_label = f"[{''.join(str(i) for i in wildcard_positions)}]" if wildcard_positions else ""
    matrix_file = workdir / f"{_with_optional_suffix(input_path.stem, 'matrix', f'{k}mers', wildcard_label)}.csv"

    # Build and save the matrix
    matrix = build_kmer_matrix(
        {**positive.raw_filtered_counts, **negative.raw_filtered_counts},
        tested_kmers,
        output_path=matrix_file,
        normalize=normalize,
        normalization_totals={**positive.totals, **negative.totals},
    )

    # Keep the saved matrix intact. Apply the zero-rate cutoff only to the
    # k-mers passed into Mann-Whitney and written to its result CSV files.
    output_kmers = _filter_kmers_by_zero_percentage(
        matrix,
        tested_kmers,
        max_zero_percentage,
    )
    if not output_kmers:
        raise ValueError(
            "No k-mers remain after applying the maximum zero-percentage cutoff."
        )

    output_prefix = workdir / _with_optional_suffix(input_path.stem, "U_test", f"{k}mers", wildcard_label)
    if progress_callback:
        progress_callback(60, "Running Mann-Whitney tests")

    result = run_mannwhitney(
        positive,
        negative,
        output_prefix=output_prefix,
        matrix_file=matrix_file,
        positive_label=positive_label,
        negative_label=negative_label,
        tested_kmers=output_kmers,
        progress_callback=progress_callback,
    )

    if progress_callback:
        progress_callback(95, "Writing output files")
    return result


def analyze_groups(
    input_path: Path,
    pos_file: Path,
    neg_file: Path,
    k_values: Iterable[int],
    wildcard_positions: Sequence[int] | None = None,
    normalize: bool = True,
    max_zero_percentage: float = 100.0,
    workdir: Optional[Path] = None,
) -> List[MannWhitneyResult]:
    """Run analysis for multiple k values."""
    workdir = Path(workdir or input_path.parent)
    workdir.mkdir(parents=True, exist_ok=True)

    wildcard_positions = list(wildcard_positions or [])
    results: List[MannWhitneyResult] = []

    for k in tqdm(list(k_values), desc="Processing k values"):
        results.append(
            analyze_single_k(
                input_path,
                pos_file,
                neg_file,
                k=k,
                wildcard_positions=wildcard_positions,
                normalize=normalize,
                max_zero_percentage=max_zero_percentage,
                workdir=workdir,
            )
        )

    return results


__all__ = [
    "AA_BACKGROUND",
    "CANONICAL_AA",
    "CHI_SQUARE_THRESHOLD",
    "CohortTilingResult",
    "MannWhitneyResult",
    "analyze_groups",
    "analyze_single_k",
    "apply_chi_square_filter",
    "build_kmer_matrix",
    "calculate_expected_frequency",
    "prebuild_kmers",
    "run_mannwhitney",
    "split_input_by_group",
    "tile_patient_file",
]
