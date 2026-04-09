from collections import defaultdict
from tqdm import tqdm


class LocalProteomeDB:
    """
    Loads a FASTA proteome and builds a k-mer -> [(protein_id, pos)] index.

    Wildcard-aware behavior:
      - ordinary k-mers use the exact prebuilt k-mer index
      - patterns containing 'X' are matched without expanding all possibilities
      - wildcard matching uses the longest concrete segment as an anchor and
        verifies full-pattern compatibility only on candidate positions
    """

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path

        self.proteins = {}          # pid -> (desc, seq)
        self.kmer_indices = {}      # k -> defaultdict(list)
        self.index = None
        self.tet_index = None

        self._load_fasta()

        # compatibility helpers used elsewhere in the project
        self.sequences = {pid: seq for pid, (_, seq) in self.proteins.items()}
        self.lengths = {pid: len(seq) for pid, seq in self.sequences.items()}

    def _load_fasta(self):
        print("\nLoading FASTA proteome:", self.fasta_path)
        pid = None
        desc = None
        seq_chunks = []

        with open(self.fasta_path, "r", encoding="utf-8", errors="ignore") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if pid:
                        self.proteins[pid] = (desc, "".join(seq_chunks).upper())

                    raw = line[1:].strip()
                    tokens = raw.split()
                    id_part = tokens[0] if tokens else raw
                    parts = id_part.split("|")
                    accession = parts[1] if len(parts) > 1 else id_part
                    description = " ".join(tokens[1:]) if len(tokens) > 1 else ""

                    desc = f"{accession} {description}".strip()
                    pid = accession
                    seq_chunks = []
                else:
                    seq_chunks.append(line)

        if pid:
            self.proteins[pid] = (desc, "".join(seq_chunks).upper())

        print(f"✔ Loaded {len(self.proteins)} proteins from FASTA")

    def _build_kmer_index(self, k: int):
        if not isinstance(k, int) or k <= 0:
            raise ValueError(f"k must be a positive integer, got {k!r}")

        if k in self.kmer_indices:
            self.index = self.kmer_indices[k]
            if k == 4:
                self.tet_index = self.kmer_indices[k]
            return self.kmer_indices[k]

        print(f"Indexing {k}-mers across proteome...")
        kmer_index = defaultdict(list)

        for pid, (_, seq) in tqdm(self.proteins.items()):
            L = len(seq)
            if L < k:
                continue
            for i in range(L - k + 1):
                kmer = seq[i:i + k].upper()
                kmer_index[kmer].append((pid, i))

        self.kmer_indices[k] = kmer_index
        self.index = kmer_index
        if k == 4:
            self.tet_index = kmer_index

        print(f"✔ Indexed {len(kmer_index)} unique {k}-mers\n")
        return kmer_index

    def ensure_k(self, k: int):
        return self._build_kmer_index(k)

    def _pattern_matches_at(self, seq: str, start: int, pattern: str) -> bool:
        k = len(pattern)
        if start < 0 or start + k > len(seq):
            return False

        window = seq[start:start + k]
        for a, b in zip(pattern, window):
            if a != "X" and a != b:
                return False
        return True

    def _longest_concrete_segment(self, pattern: str):
        best = None
        best_start = None
        run_start = None

        for i, ch in enumerate(pattern):
            if ch != "X":
                if run_start is None:
                    run_start = i
            else:
                if run_start is not None:
                    seg = pattern[run_start:i]
                    if best is None or len(seg) > len(best):
                        best = seg
                        best_start = run_start
                    run_start = None

        if run_start is not None:
            seg = pattern[run_start:]
            if best is None or len(seg) > len(best):
                best = seg
                best_start = run_start

        return best, best_start

    def search_pattern(self, pattern: str, k: int | None = None):
        pattern = str(pattern).strip().upper()
        if not pattern:
            return []

        if k is None:
            k = len(pattern)
        if len(pattern) != k:
            raise ValueError(f"Pattern length {len(pattern)} does not match k={k}")

        if "X" not in pattern:
            return self.search_kmer(pattern, k)

        anchor, offset = self._longest_concrete_segment(pattern)
        if not anchor:
            # All wildcards: every window of length k is a match.
            hits = []
            for pid, seq in self.sequences.items():
                L = len(seq)
                if L < k:
                    continue
                hits.extend((pid, i) for i in range(L - k + 1))
            return hits

        hits = []
        for pid, seq in self.sequences.items():
            search_from = 0
            while True:
                found = seq.find(anchor, search_from)
                if found == -1:
                    break
                start = found - offset
                if self._pattern_matches_at(seq, start, pattern):
                    hits.append((pid, start))
                search_from = found + 1
        return hits

    def search_kmer(self, kmer: str, k: int | None = None):
        kmer = str(kmer).strip().upper()
        if not kmer:
            return []
        if k is None:
            k = len(kmer)
        index = self._build_kmer_index(k)
        return index.get(kmer, [])

    def search_tetramer(self, tetramer):
        return self.search_kmer(tetramer, 4)

    def get_protein_desc(self, pid):
        return self.proteins[pid][0]

    def get_protein_len(self, pid):
        return len(self.proteins[pid][1])
