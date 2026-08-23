from pathlib import Path

from backend.app.core.fastq_processing import _normalize_peptide, _normalize_sample_name


def test_x_is_replaced_only_in_peptide_content() -> None:
    assert _normalize_peptide("ACDXFGHIKLMN") == "ACDQFGHIKLMN"

    # File and sample names never pass through peptide normalization.
    input_path = Path("AD1.txt")
    assert input_path.name == "AD1.txt"
    assert _normalize_sample_name(input_path) == "AD1"


def test_x_in_filename_is_not_treated_as_a_peptide_residue() -> None:
    input_path = Path("ADX1.txt")
    assert input_path.name == "ADX1.txt"
    assert _normalize_sample_name(input_path) == "ADX1"
