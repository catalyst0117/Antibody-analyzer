from __future__ import annotations

import gzip
import re
import shutil
import uuid
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from Bio.Seq import Seq
import pandas as pd


VALID_BASES = {"A", "T", "C", "G"}


def _ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _normalize_sample_name(path: Path) -> str:
    stem = path.stem
    if path.suffix == ".gz":
        stem = Path(path.stem).stem
    return stem.split(".")[0]


@dataclass
class SampleSummary:
    sample_name: str
    total_sequences: int
    unique_sequences: int
    filtered_sequences: int
    output_columns: List[str] = field(default_factory=list)


@dataclass
class FastqResult:
    summary: List[SampleSummary]
    excel_path: Path
    filtered_outputs: Dict[str, Path]
    background_dump: Optional[Path]


class FastqProcessor:
    """High level wrapper around the FASTQ processing workflow used by module1."""

    def __init__(self, workdir: Path):
        self.workdir = _ensure_directory(workdir)
        self.pattern = re.compile(r"TTCTCACTCT.{0,36}")

    def _extract_peptides_from_fastq(self, fastq_path: Path) -> Iterable[str]:
        peptides = set()
        with fastq_path.open("r") as handle:
            for line in handle:
                if line.startswith("@") or line.startswith("+"):
                    continue
                match = self.pattern.search(line)
                if not match:
                    continue
                nt_seq = match.group(0)[10:]
                if len(nt_seq) != 36 or set(nt_seq) - VALID_BASES:
                    continue
                aa = str(Seq(nt_seq).translate()).replace("*", "X")
                if len(aa) >= 12:
                    peptides.add(aa)
        return peptides

    def _prepare_background(self, background_file: Optional[Path]) -> tuple[set[str], Optional[Path]]:
        if not background_file:
            return set(), None
        if not background_file.exists():
            raise FileNotFoundError(f"Background file {background_file} does not exist")

        peptides = set(self._extract_peptides_from_fastq(background_file))
        if not peptides:
            return set(), None

        dump_path = self.workdir / "processed_background_peptides.txt"
        with dump_path.open("w") as fh:
            for seq in sorted(peptides):
                fh.write(seq + "\n")
        return peptides, dump_path

    def _open_fastq(self, path: Path) -> Path:
        if path.suffix == ".gz":
            target = self.workdir / f"{uuid.uuid4().hex}_{path.stem}"
            with gzip.open(path, "rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            return target
        return path

    def process(
        self,
        fastq_files: Iterable[Path],
        background_file: Optional[Path] = None,
        output_name: str = "sequence_matrix.xlsx",
    ) -> FastqResult:
        fastq_files = list(fastq_files)
        if not fastq_files:
            raise ValueError("No FASTQ files provided")

        background_set, background_dump = self._prepare_background(background_file)

        sample_counters: Dict[str, Counter[str]] = {}
        filtered_outputs: Dict[str, Path] = {}
        summaries: List[SampleSummary] = []

        for file_path in fastq_files:
            sample_name = _normalize_sample_name(file_path)
            opened_path = self._open_fastq(file_path)
            translated: List[str] = []
            filtered: List[str] = []

            with opened_path.open("r") as handle:
                for line in handle:
                    match = self.pattern.search(line)
                    if not match:
                        continue
                    seq = match.group(0)[10:]
                    if len(seq) != 36 or set(seq) - VALID_BASES:
                        continue
                    aa = str(Seq(seq).translate()).replace("*", "X")
                    if len(aa) < 12:
                        continue
                    if aa in background_set:
                        filtered.append(aa)
                    else:
                        translated.append(aa)

            counter = Counter(translated)
            sample_counters[sample_name] = counter

            filtered_path = self.workdir / f"{sample_name}_filtered_out.txt"
            if filtered:
                with filtered_path.open("w") as fh:
                    for seq in sorted(set(filtered)):
                        fh.write(seq + "\n")
                filtered_outputs[sample_name] = filtered_path

            summaries.append(
                SampleSummary(
                    sample_name=sample_name,
                    total_sequences=len(translated) + len(filtered),
                    unique_sequences=len(counter),
                    filtered_sequences=len(filtered),
                )
            )

        df_parts: List[pd.DataFrame] = []
        for sample, counter in sample_counters.items():
            sorted_pairs = sorted(counter.items(), key=lambda item: -item[1])
            sequences = [seq for seq, _ in sorted_pairs]
            counts = [count for _, count in sorted_pairs]
            sample_df = pd.DataFrame({
                f"{sample}_Sequence": sequences,
                f"{sample}_Count": counts,
            })
            df_parts.append(sample_df)
            for summary in summaries:
                if summary.sample_name == sample:
                    summary.output_columns = list(sample_df.columns)
                    break

        max_len = max(len(df) for df in df_parts)
        aligned_parts = [df.reindex(range(max_len)) for df in df_parts]
        final_df = pd.concat(aligned_parts, axis=1)

        excel_path = self.workdir / (output_name if output_name.endswith(".xlsx") else f"{output_name}.xlsx")
        final_df.to_excel(excel_path, index=False)

        return FastqResult(
            summary=summaries,
            excel_path=excel_path,
            filtered_outputs=filtered_outputs,
            background_dump=background_dump,
        )

__all__ = ["FastqProcessor", "FastqResult", "SampleSummary"]
