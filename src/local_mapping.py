import pandas as pd
from collections import defaultdict
from tqdm import tqdm


def compute_sequence_coverage(hits, protein_length):
    """
    hits: list of (kmer, pos, qv)
    protein_length: int

    Returns:
        coverage_ratio (float)
        covered_aa (int)
    """
    if not hits or protein_length == 0:
        return 0.0, 0

    intervals = []
    for kmer, pos, _ in hits:
        k = len(kmer)
        intervals.append((pos, pos + k))

    intervals.sort(key=lambda x: x[0])

    covered = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start > cur_end:
            covered += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    covered += cur_end - cur_start

    return covered / protein_length, covered


def map_local_kmers(
    kmer_file,
    db,
    out_path,
    heap_size=25,
    pos_diff=100,
    wildcards=False,
    q_cutoff=0.01,
):
    """
    Local FASTA proteome mapper, generalized to k-mer-aware input.

    Wildcard behavior:
      - if wildcards=False, everything is treated as exact string matching
      - if wildcards=True, any k-mer containing 'X' is matched as a wildcard
        pattern, where X can be any amino acid
      - wildcard matching is done by anchor-based scanning, not by expanding
        X into every possible amino-acid combination
    """

    data = {}
    total_count = 0
    kept_count = 0
    kept_lengths = set()

    with open(kmer_file) as f:
        for line_no, line in enumerate(f, start=1):
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue

            kmer, qv = parts
            kmer = str(kmer).strip().upper()
            if not kmer:
                continue

            qv = float(qv)
            total_count += 1

            if qv <= q_cutoff:
                data[kmer] = qv
                kept_count += 1
                kept_lengths.add(len(kmer))

    if not data:
        print(f"Loaded {total_count} kmers, kept 0 with q ≤ {q_cutoff}")
        df = pd.DataFrame(columns=[
            "Description",
            "Q score",
            "Length of the Protein",
            "# of hits",
            "Hits/Protein length",
            "Sequence coverage",
            "Covered amino acids",
            "K-mer, Position",
        ])
        df.to_csv(out_path, index=False)
        print("✔ Local mapping saved to:", out_path)
        return out_path, []

    if len(kept_lengths) != 1:
        raise ValueError(
            f"Input contains mixed k-mer lengths after q-value filtering: {sorted(kept_lengths)}. "
            f"Module 3 local mapping expects one k per run."
        )

    k = next(iter(kept_lengths))
    db.ensure_k(k)

    wildcard_count = sum(1 for x in data if "X" in x)
    print(f"Loaded {total_count} kmers, kept {kept_count} with q ≤ {q_cutoff}")
    print(f"Detected k={k} from input content")
    print(f"Mapping {len(data)} {k}-mers using local proteome...")
    if wildcards:
        print(f"Wildcard mode ON: {wildcard_count} pattern(s) contain X")

    protein_hits = defaultdict(list)
    for kmer, qv in tqdm(data.items(), desc=f"Mapping {k}-mers"):
        if wildcards and "X" in kmer:
            hits = db.search_pattern(kmer, k)
        else:
            hits = db.search_kmer(kmer, k)

        for pid, pos in hits:
            protein_hits[pid].append((kmer, pos, qv))

    rows = []
    manhattan = []

    for pid, hits in tqdm(protein_hits.items(), desc="Scoring proteins"):
        desc = db.get_protein_desc(pid)
        plen = db.get_protein_len(pid)
        hits_sorted = sorted(hits, key=lambda x: x[1])

        z = sum(1 / qv for (_, _, qv) in hits_sorted) / len(hits_sorted)
        seq_cov_ratio, covered_aa = compute_sequence_coverage(hits_sorted, plen)

        rows.append([
            desc,
            z,
            plen,
            len(hits_sorted),
            len(hits_sorted) / plen,
            seq_cov_ratio,
            covered_aa,
            str([(t, p) for (t, p, _) in hits_sorted]),
        ])

    df = pd.DataFrame(rows, columns=[
        "Description",
        "Q score",
        "Length of the Protein",
        "# of hits",
        "Hits/Protein length",
        "Sequence coverage",
        "Covered amino acids",
        "K-mer, Position",
    ])

    df.to_csv(out_path, index=False)
    print("✔ Local mapping saved to:", out_path)

    return out_path, manhattan


map_local_tetramers = map_local_kmers
