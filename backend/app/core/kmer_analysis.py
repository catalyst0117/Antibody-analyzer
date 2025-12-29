from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

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


def split_input_by_group(input_path: Path) -> Tuple[Path, Path]:
    if input_path.suffix == ".xlsx":
        df = pd.read_excel(input_path)
    elif input_path.suffix == ".csv":
        df = pd.read_csv(input_path)
    else:
        raise ValueError("Unsupported file format.")

    ad_cols: List[str] = []
    nc_cols: List[str] = []

    i = 0
    while i < df.shape[1] - 1:
        header = str(df.columns[i]).strip()
        if not header or header.startswith("Unnamed:"):
            i += 1
            continue
        if header.startswith("AD_"):
            ad_cols.extend(df.columns[i : i + 2])
            i += 2
        elif header.startswith("NC_"):
            nc_cols.extend(df.columns[i : i + 2])
            i += 2
        else:
            i += 1

    df_ad = df[ad_cols].copy()
    df_nc = df[nc_cols].copy()

    for subset in (df_ad, df_nc):
        for idx in range(1, subset.shape[1], 2):
            subset.columns.values[idx] = "count"

    output_ad = input_path.parent / f"{input_path.stem}_AD{input_path.suffix}"
    output_nc = input_path.parent / f"{input_path.stem}_NC{input_path.suffix}"

    if input_path.suffix == ".xlsx":
        df_ad.to_excel(output_ad, index=False)
        df_nc.to_excel(output_nc, index=False)
    else:
        df_ad.to_csv(output_ad, index=False)
        df_nc.to_csv(output_nc, index=False)

    return output_ad, output_nc


def apply_chi_square_filter(
    kmer_dict: Dict[str, int],
    aa_background: Dict[str, float],
    threshold: float = CHI_SQUARE_THRESHOLD,
    normalize_expected: bool = True,
) -> Dict[str, int]:
    filtered: Dict[str, int] = {}
    total_count = sum(kmer_dict.values())
    if total_count == 0:
        return filtered

    for kmer, observed in kmer_dict.items():
        expected = 0.0
        for aa in kmer:
            expected += aa_background.get(aa, 0.0) * total_count
        if normalize_expected and kmer_dict:
            expected /= len(kmer_dict)
        if expected == 0:
            continue
        chi2 = (observed - expected) ** 2 / expected
        if chi2 >= threshold:
            filtered[kmer] = observed
    return filtered


def tile_patient_file(
    path: Path,
    kmer_length: int,
    wildcard_positions: Sequence[int] | None = None,
    apply_chi_square: bool = True,
    aa_background: Optional[Dict[str, float]] = None,
    normalize_expected: bool = True,
) -> Tuple[Dict[str, Dict[str, int]], set[str]]:
    wildcard_positions = list(wildcard_positions or [])
    if path.suffix == ".xlsx":
        df = pd.read_excel(path)
    elif path.suffix == ".csv":
        df = pd.read_csv(path)
    else:
        raise ValueError("Unsupported file format.")

    patient_dicts: Dict[str, Dict[str, int]] = {}
    total_kmer_dict: Dict[str, int] = defaultdict(int)

    for col in range(0, df.shape[1], 2):
        seq_col = df.columns[col]
        count_col = df.columns[col + 1]
        patient_name = str(seq_col)
        patient_dict: Dict[str, int] = defaultdict(int)

        for seq, count in zip(df[seq_col], df[count_col]):
            if pd.isna(seq) or pd.isna(count):
                continue
            if not isinstance(seq, str):
                continue
            try:
                count = int(count)
            except Exception:
                continue

            for i in range(len(seq) - kmer_length + 1):
                kmer = list(seq[i : i + kmer_length])
                for pos in wildcard_positions:
                    if 0 <= pos < len(kmer):
                        kmer[pos] = "X"
                kmer_str = "".join(kmer)
                patient_dict[kmer_str] += count
                total_kmer_dict[kmer_str] += count

        patient_dicts[patient_name] = dict(patient_dict)

    filtered_dict = dict(total_kmer_dict)
    if apply_chi_square:
        if aa_background is None:
            aa_background = AA_BACKGROUND
        filtered_dict = apply_chi_square_filter(
            filtered_dict,
            aa_background=aa_background,
            threshold=CHI_SQUARE_THRESHOLD,
            normalize_expected=normalize_expected,
        )

    return patient_dicts, set(filtered_dict.keys())


def build_kmer_matrix(
    patient_dicts: Dict[str, Dict[str, int]],
    kmers_filter: Optional[Iterable[str]] = None,
    normalize: bool = True,
) -> pd.DataFrame:
    kmers = set()
    for data in patient_dicts.values():
        kmers.update(data.keys())
    if kmers_filter is not None:
        kmers &= set(kmers_filter)
    kmers = sorted(kmers)

    matrix = pd.DataFrame(index=kmers)
    for patient, counts in patient_dicts.items():
        series = pd.Series(counts, dtype=float)
        if normalize:
            total = series.sum()
            if total:
                series = series / total
        matrix[patient] = series
    matrix.fillna(0.0, inplace=True)
    return matrix


def run_mannwhitney(
    matrix: pd.DataFrame,
    pos_cols: Sequence[str],
    neg_cols: Sequence[str],
    output_prefix: Path,
    matrix_file: Path,
) -> MannWhitneyResult:
    kmers = matrix.index
    pvals: List[float] = []
    mean_rank_diffs: List[float] = []

    for kmer in kmers:
        pos_vals = matrix.loc[kmer, pos_cols].to_numpy()
        neg_vals = matrix.loc[kmer, neg_cols].to_numpy()

        combined = np.concatenate([pos_vals, neg_vals])
        ranks = rankdata(combined)
        rank_pos = ranks[: len(pos_vals)]
        rank_neg = ranks[len(pos_vals) :]
        mean_rank_diff = float(rank_pos.mean() - rank_neg.mean())

        try:
            _, p_value = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")
        except ValueError:
            p_value = 1.0

        pvals.append(p_value)
        mean_rank_diffs.append(mean_rank_diff)

    fdr = multipletests(pvals, method="fdr_tsbky")[1]

    result_df = pd.DataFrame(
        {
            "kmer": kmers,
            "p_value": pvals,
            "mean_rank_diff": mean_rank_diffs,
            "fdr_corrected_p": fdr,
        }
    ).sort_values("p_value")

    result_file = output_prefix.with_suffix(".csv")
    ad_file = output_prefix.with_name(output_prefix.stem + "_AD.csv")
    nc_file = output_prefix.with_name(output_prefix.stem + "_NC.csv")

    result_df.to_csv(result_file, index=False)
    result_df[result_df["mean_rank_diff"] > 0].to_csv(ad_file, index=False)
    result_df[result_df["mean_rank_diff"] < 0].to_csv(nc_file, index=False)

    return MannWhitneyResult(
        k=len(kmers[0]) if kmers else 0,
        total_kmers=len(result_df),
        ad_elevated=int((result_df["mean_rank_diff"] > 0).sum()),
        nc_elevated=int((result_df["mean_rank_diff"] < 0).sum()),
        result_file=result_file,
        ad_file=ad_file,
        nc_file=nc_file,
        matrix_file=matrix_file,
    )


def analyze_groups(
    input_path: Path,
    pos_file: Path,
    neg_file: Path,
    k_values: Iterable[int],
    wildcard_positions: Sequence[int] | None = None,
    normalize: bool = True,
    workdir: Optional[Path] = None,
) -> List[MannWhitneyResult]:
    workdir = Path(workdir or input_path.parent)
    wildcard_positions = list(wildcard_positions or [])
    results: List[MannWhitneyResult] = []

    for k in k_values:
        pos_dicts, filtered_pos = tile_patient_file(
            pos_file,
            kmer_length=k,
            wildcard_positions=wildcard_positions,
        )
        neg_dicts, filtered_neg = tile_patient_file(
            neg_file,
            kmer_length=k,
            wildcard_positions=wildcard_positions,
        )

        filter_set = filtered_pos & filtered_neg if filtered_pos and filtered_neg else filtered_pos or filtered_neg
        matrix = build_kmer_matrix({**pos_dicts, **neg_dicts}, kmers_filter=filter_set, normalize=normalize)

        wildcard_label = ''.join(str(i) for i in wildcard_positions) or 'no_wildcards'
        matrix_file = workdir / f"{input_path.stem}_matrix_{k}mers_{wildcard_label}.csv"
        matrix.to_csv(matrix_file)

        pos_cols = list(pos_dicts.keys())
        neg_cols = list(neg_dicts.keys())

        output_prefix = workdir / f"{input_path.stem}_U_test_{k}mers_{wildcard_label}"
        result = run_mannwhitney(matrix, pos_cols, neg_cols, output_prefix=output_prefix, matrix_file=matrix_file)
        results.append(result)

    return results


__all__ = [
    "AA_BACKGROUND",
    "CHI_SQUARE_THRESHOLD",
    "MannWhitneyResult",
    "analyze_groups",
    "build_kmer_matrix",
    "run_mannwhitney",
    "split_input_by_group",
    "tile_patient_file",
]
