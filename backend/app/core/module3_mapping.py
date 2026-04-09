from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
import re


@dataclass
class Module3Result:
    output_folder_name: str
    positive_mapping_file: Path
    negative_mapping_file: Path
    positive_clean_file: Path
    negative_clean_file: Path
    positive_manhattan_file: Path
    negative_manhattan_file: Path
    run_summary_file: Path


@lru_cache(maxsize=1)
def _load_module3_module():
    repo_root = Path(__file__).resolve().parents[3]
    module_path = repo_root / "src" / "module3.py"
    spec = importlib.util.spec_from_file_location("antibody_analyzer_module3", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load Proteome Mapping support from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sanitize_folder_name(name: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")
    return sanitized or "proteome_mapping"


def resolve_output_folder_name(
    *,
    requested_name: str | None,
    source_filename: str,
    top_n: int | None,
) -> str:
    if requested_name and requested_name.strip():
        return _sanitize_folder_name(requested_name)

    stem = Path(source_filename).stem
    stem = re.sub(
        r"(?i)(?:[_-]?)(positive_output|negative_output|ad_output|nc_output)$",
        "",
        stem,
    )
    stem = re.sub(r"(?i)(?:[_-]?)(positive|negative|ad|nc)$", "", stem)
    stem = stem.rstrip("_- ")
    if not stem:
        stem = "proteome_mapping"

    top_label = top_n if top_n is not None else "all"
    return _sanitize_folder_name(f"{stem}_top{top_label}")


def build_run_summary_text(
    *,
    output_folder_name: str,
    positive_input_filename: str,
    negative_input_filename: str,
    fasta_filename: str,
    top_n: int | None,
    wildcards: bool,
    q_cutoff: float,
    output_filenames: list[str],
) -> str:
    timestamp = datetime.now().isoformat(timespec="seconds")
    wildcard_label = "Enabled" if wildcards else "Disabled"
    top_label = str(top_n) if top_n is not None else "All rows"

    return "\n".join(
        [
            "Proteome Mapping Outputs",
            "",
            "These files are included in the downloadable ZIP bundle.",
            "",
            "Settings",
            f"top_n: {top_label}",
            f"Wildcard matching: {wildcard_label}",
            f"q_cutoff: {q_cutoff}",
            f"Output folder name: {output_folder_name}",
            f"Positive input filename: {positive_input_filename}",
            f"Negative input filename: {negative_input_filename}",
            f"FASTA/reference filename: {fasta_filename}",
            f"Timestamp: {timestamp}",
            f"Output files: {', '.join(output_filenames)}",
            "",
        ]
    )


def run_module3_mapping(
    *,
    positive_file: Path,
    negative_file: Path,
    fasta_file: Path,
    output_dir: Path,
    output_folder_name: str,
    top_n: int | None,
    wildcards: bool,
    q_cutoff: float,
) -> Module3Result:
    module3 = _load_module3_module()
    output_dir.mkdir(parents=True, exist_ok=True)

    module3.run_mapping(
        str(positive_file),
        str(negative_file),
        str(fasta_file),
        str(output_dir),
        top_n=top_n,
        wildcards=wildcards,
        q_cutoff=q_cutoff,
    )

    output_filenames = [
        "DictFilePos.csv",
        "DictFileNeg.csv",
        "clean_pos.txt",
        "clean_neg.txt",
        "manhattan_pos.json",
        "manhattan_neg.json",
    ]
    run_summary_file = output_dir / "run_summary.txt"
    run_summary_file.write_text(
        build_run_summary_text(
            output_folder_name=output_folder_name,
            positive_input_filename=positive_file.name,
            negative_input_filename=negative_file.name,
            fasta_filename=fasta_file.name,
            top_n=top_n,
            wildcards=wildcards,
            q_cutoff=q_cutoff,
            output_filenames=output_filenames,
        ),
        encoding="utf-8",
    )

    return Module3Result(
        output_folder_name=output_folder_name,
        positive_mapping_file=output_dir / "DictFilePos.csv",
        negative_mapping_file=output_dir / "DictFileNeg.csv",
        positive_clean_file=output_dir / "clean_pos.txt",
        negative_clean_file=output_dir / "clean_neg.txt",
        positive_manhattan_file=output_dir / "manhattan_pos.json",
        negative_manhattan_file=output_dir / "manhattan_neg.json",
        run_summary_file=run_summary_file,
    )


__all__ = [
    "Module3Result",
    "build_run_summary_text",
    "resolve_output_folder_name",
    "run_module3_mapping",
]
