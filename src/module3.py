import json
import os
import re
import sys
from typing import Optional

import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

from LocalProteomeDB import LocalProteomeDB
from local_mapping import map_local_kmers


AA_PATTERN_WITH_WILDCARD = re.compile(r"^[A-Z]+$")


def _normalize_header_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).strip().lower())


def _normalize_kmer_value(val) -> str:
    if pd.isna(val):
        return ""
    return re.sub(r"\s+", "", str(val).strip().upper())


def _find_qvalue_column(df: pd.DataFrame) -> str:
    candidates = []
    for col in df.columns:
        norm = _normalize_header_name(col)
        if norm in {"qvalue", "qvalues", "qval", "qvals"}:
            candidates.append(col)
    if candidates:
        return candidates[0]
    raise ValueError(
        f"Could not find a Q-value column. Available columns: {list(df.columns)}"
    )


def _find_sequence_column(df: pd.DataFrame) -> str:
    """
    Find the sequence column robustly.

    Header typos like 'tetramer', 'tetramers', '6 mer', '6-mer', etc.
    should not break the run. We therefore score columns by both header hints
    and content consistency.
    """
    best_col = None
    best_score = float("-inf")

    for col in df.columns:
        series = df[col].dropna().astype(str)
        if series.empty:
            continue

        sample = series.head(200)
        normalized_values = [_normalize_kmer_value(v) for v in sample]
        alpha_values = [v for v in normalized_values if v and AA_PATTERN_WITH_WILDCARD.fullmatch(v)]
        if not alpha_values:
            continue

        lengths = [len(v) for v in alpha_values]
        mode_len = pd.Series(lengths).mode().iloc[0]
        consistent = sum(1 for x in lengths if x == mode_len)
        alpha_ratio = len(alpha_values) / max(len(sample), 1)
        consistent_ratio = consistent / max(len(alpha_values), 1)

        norm_header = _normalize_header_name(col)
        header_bonus = 0
        if any(token in norm_header for token in [
            "kmer", "kmers", "peptide", "peptides", "sequence", "seq", "motif", "motifs"
        ]):
            header_bonus += 2
        if any(token in norm_header for token in [
            "tetramer", "tetramers", "trimer", "trimers", "pentamer", "pentamers", "hexamer", "hexamers"
        ]):
            header_bonus += 2
        if re.search(r"\b\d+\s*[- ]?mer[s]?\b", str(col), flags=re.IGNORECASE):
            header_bonus += 2

        score = (alpha_ratio * 5) + (consistent_ratio * 5) + header_bonus
        if score > best_score:
            best_score = score
            best_col = col

    if best_col is None:
        raise ValueError(
            "Could not identify the k-mer sequence column from the input file. "
            f"Available columns: {list(df.columns)}"
        )

    return best_col


def _infer_single_k_from_series(series: pd.Series) -> int:
    normalized = [_normalize_kmer_value(v) for v in series.dropna()]
    kmers = [v for v in normalized if v]
    if not kmers:
        raise ValueError("The detected sequence column is empty.")

    lengths = sorted({len(v) for v in kmers})
    if len(lengths) != 1:
        raise ValueError(
            f"The detected sequence column contains mixed sequence lengths: {lengths}. "
            "Module 3 local mapping expects one k per run."
        )
    return lengths[0]


def extract_kmer_list(input_path, output_txt, top_n: Optional[int] = None):
    if input_path.endswith(".txt"):
        raise ValueError(
            f"extract_kmer_list() expected a CSV/XLS/XLSX table with headers, "
            f"but got a text file instead: {input_path}"
        )

    if input_path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(input_path, header=0)
    else:
        try:
            df = pd.read_csv(
                input_path,
                encoding="utf-8-sig",
                header=0,
                skip_blank_lines=True,
            )
        except Exception:
            df = pd.read_csv(
                input_path,
                sep=None,
                engine="python",
                encoding="utf-8-sig",
                header=0,
                skip_blank_lines=True,
            )

    df.columns = [str(c).strip() for c in df.columns]

    seq_col = _find_sequence_column(df)
    q_col = _find_qvalue_column(df)

    df = df[[seq_col, q_col]].copy()
    df[seq_col] = df[seq_col].map(_normalize_kmer_value)
    df[q_col] = pd.to_numeric(df[q_col], errors="coerce")
    df = df.dropna(subset=[seq_col, q_col])
    df = df[df[seq_col] != ""]

    inferred_k = _infer_single_k_from_series(df[seq_col])

    bad_rows = df[~df[seq_col].str.fullmatch(r"[A-Z]+")]
    if not bad_rows.empty:
        raise ValueError(
            f"Detected non-alphabetic sequence values in column '{seq_col}'. "
            f"Example values: {bad_rows[seq_col].head(5).tolist()}"
        )

    df = df[df[seq_col].str.len() == inferred_k]
    df = df.sort_values(q_col, ascending=True)

    if top_n is not None:
        df = df.head(top_n)

    with open(output_txt, "w") as f:
        for _, row in df.iterrows():
            f.write(f"{row[seq_col]}\t{row[q_col]}\n")

    print(
        f"✔ Wrote {len(df)} {inferred_k}-mers "
        f"(top {top_n if top_n is not None else 'ALL'}) to {output_txt} "
        f"using sequence column '{seq_col}' and q-value column '{q_col}'"
    )
    return output_txt, inferred_k, seq_col, q_col

def _count_unique_kmers_from_clean_txt(clean_txt_path: str) -> tuple[int, int]:
    kmers = set()
    lengths = set()
    with open(clean_txt_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            kmer = line.split("\t")[0].strip().upper()
            if kmer:
                kmers.add(kmer)
                lengths.add(len(kmer))

    if not kmers:
        raise ValueError(f"No k-mers found in clean txt file: {clean_txt_path}")
    if len(lengths) != 1:
        raise ValueError(
            f"Mixed k-mer lengths found in clean txt file {clean_txt_path}: {sorted(lengths)}"
        )
    return len(kmers), next(iter(lengths))


def add_random_expected_hits(dict_csv_path: str, clean_txt_path: str) -> None:
    """
    Adds 2 columns to DictFile*.csv:
      - ExpectedRandomHits
      - Hits/ExpectedRandomHits

    K = number of unique input k-mers
    universe_size = 20**k
    expected hits = max(L-k+1, 0) * (K / 20**k)

    Wildcards are intentionally not expanded here. The goal is to keep the same
    simple baseline as before, now generalized from fixed 4-mers to arbitrary k.
    """
    df = pd.read_csv(dict_csv_path)
    if df.empty:
        print(f"⚠ Dict file is empty, skipping: {dict_csv_path}")
        return

    length_col = None
    for c in df.columns:
        if str(c).strip().lower() in {
            "length of the protein", "protein length", "proteinlength", "length_aa", "length"
        }:
            length_col = c
            break
    if length_col is None:
        for c in df.columns:
            if "length" in str(c).lower():
                length_col = c
                break
    if length_col is None:
        raise ValueError(
            f"Could not find a protein length column in {dict_csv_path}. "
            f"Expected something like 'Length of the Protein'."
        )

    hits_col = None
    for c in df.columns:
        if str(c).strip().lower() in {"# of hits", "hits", "hit", "hit_count", "hitcount", "num_hits"}:
            hits_col = c
            break
    if hits_col is None:
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c != length_col]
        if not numeric_cols:
            raise ValueError(f"Could not find a '# of hits' numeric column in {dict_csv_path}.")
        hits_col = numeric_cols[0]

    K, k = _count_unique_kmers_from_clean_txt(clean_txt_path)
    universe_size = 20 ** k
    p = K / universe_size

    def _expected(L):
        if pd.isna(L):
            return None
        try:
            L = int(float(L))
        except Exception:
            return None
        windows = max(L - k + 1, 0)
        return windows * p

    expected_col = "ExpectedRandomHits"
    ratio_col = "Hits/ExpectedRandomHits"

    df[expected_col] = df[length_col].apply(_expected)

    def _ratio(row):
        exp = row[expected_col]
        obs = row[hits_col]
        if pd.isna(exp) or exp == 0:
            return None
        return float(obs) / float(exp)

    df[ratio_col] = df.apply(_ratio, axis=1)
    df.to_csv(dict_csv_path, index=False)
    print(f"✔ Added random-hit expectation columns to: {dict_csv_path}")


def run_mapping(
    pos_sig_path: str,
    neg_sig_path: str,
    fasta_path: str,
    out_dir: str,
    top_n: Optional[int] = None,
    wildcards: bool = False,
    q_cutoff: float = 0.01,
):
    """
    Local-proteome-only Module 3 runner.

    This version intentionally excludes SQL mode and old check-only helpers.
    """
    os.makedirs(out_dir, exist_ok=True)

    clean_pos = os.path.join(out_dir, "clean_pos.txt")
    clean_neg = os.path.join(out_dir, "clean_neg.txt")
    pos_out = os.path.join(out_dir, "DictFilePos.csv")
    neg_out = os.path.join(out_dir, "DictFileNeg.csv")
    man_pos = os.path.join(out_dir, "manhattan_pos.json")
    man_neg = os.path.join(out_dir, "manhattan_neg.json")

    clean_pos, pos_k, pos_seq_col, pos_q_col = extract_kmer_list(pos_sig_path, clean_pos, top_n=top_n)
    clean_neg, neg_k, neg_seq_col, neg_q_col = extract_kmer_list(neg_sig_path, clean_neg, top_n=top_n)

    if pos_k != neg_k:
        raise ValueError(
            f"Positive and negative files produced different k values: {pos_k} vs {neg_k}. "
            "Run Module 3 separately for each k."
        )

    proteome = LocalProteomeDB(fasta_path)

    ret_pos, mp1 = map_local_kmers(
        clean_pos, proteome, pos_out, wildcards=wildcards, q_cutoff=q_cutoff
    )
    ret_neg, mp2 = map_local_kmers(
        clean_neg, proteome, neg_out, wildcards=wildcards, q_cutoff=q_cutoff
    )

    add_random_expected_hits(pos_out, clean_pos)
    add_random_expected_hits(neg_out, clean_neg)

    with open(man_pos, "w") as f:
        json.dump(mp1, f)
    with open(man_neg, "w") as f:
        json.dump(mp2, f)

    print("\n=== COMPLETE ===")
    print(f"Positive input: sequence='{pos_seq_col}', q='{pos_q_col}', k={pos_k}")
    print(f"Negative input: sequence='{neg_seq_col}', q='{neg_q_col}', k={neg_k}")
    return ret_pos, ret_neg


if __name__ == "__main__":
    POS_SIG = "/Users/ciao/Downloads/Wildcard_test/SI_file1_U_test_5mers_3_AD.csv"
    NEG_SIG = "/Users/ciao/Downloads/Wildcard_test/SI_file1_U_test_5mers_3_NC.csv"
    FASTA_PATH = "/Users/ciao/Desktop/module3_test/UinProtHumanProteoCanonicalComplete.fasta.txt"
    OUT_DIR = "/Users/ciao/Downloads/Wildcard_test/"
    top_n = 50

    run_mapping(
        POS_SIG,
        NEG_SIG,
        FASTA_PATH,
        OUT_DIR,
        top_n=top_n,
        wildcards=True,
        q_cutoff=0.01,
    )

    print("\nAll done!\n")
