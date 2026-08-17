from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
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


def run_mannwhitney(
    positive: CohortTilingResult,
    negative: CohortTilingResult,
    output_prefix: Path,
    matrix_file: Path,
    positive_label: str = "AD",
    negative_label: str = "NC",
    progress_callback: Optional[ProgressCallback] = None,
) -> MannWhitneyResult:
    """Test the union of passing k-mers using normalized filtered values.

    The directional AD/NC files apply significance cutoff to raw p-value.
    """
    tested_kmers = sorted(positive.passed_kmers | negative.passed_kmers)
    if not tested_kmers:
        raise ValueError("No k-mers passed chi-square filtering in either cohort.")

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
    positive_file_label = _group_label_for_filename(positive_label, "positive")
    negative_file_label = _group_label_for_filename(negative_label, "negative")
    ad_file = output_prefix.with_name(f"{output_prefix.stem}_{positive_file_label}.csv")
    nc_file = output_prefix.with_name(f"{output_prefix.stem}_{negative_file_label}.csv")

    result_df.to_csv(result_file, index=False)
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
    )


def analyze_single_k(
    input_path: Path,
    pos_file: Path,
    neg_file: Path,
    k: int,
    wildcard_positions: Sequence[int] | None = None,
    normalize: bool = True,
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
