from pathlib import Path

import pytest

from backend.app.core.kmer_analysis import (
    TETRAMER_LENGTH,
    calculate_expected_frequency,
    prebuild_kmers,
)
from backend.app.core.module3_mapping import run_module3_mapping


def test_prebuild_kmers_builds_exact_tetramers() -> None:
    kmers = prebuild_kmers(TETRAMER_LENGTH, alphabet="AC")

    assert len(kmers) == 16
    assert set(kmers) == {
        "AAAA", "AAAC", "AACA", "AACC",
        "ACAA", "ACAC", "ACCA", "ACCC",
        "CAAA", "CAAC", "CACA", "CACC",
        "CCAA", "CCAC", "CCCA", "CCCC",
    }


@pytest.mark.parametrize("k", [1, 3, 5, 8])
def test_prebuild_kmers_rejects_non_tetramer_lengths(k: int) -> None:
    with pytest.raises(ValueError, match="tetramers"):
        prebuild_kmers(k, alphabet="AC")


def test_prebuild_kmers_rejects_wildcard_positions() -> None:
    with pytest.raises(ValueError, match="Wildcard"):
        prebuild_kmers(TETRAMER_LENGTH, wildcard_positions=[1], alphabet="AC")


def test_expected_frequency_rejects_wildcard_sequences() -> None:
    with pytest.raises(ValueError, match="Wildcard"):
        calculate_expected_frequency("AXAA", total_count=10)


def test_proteome_mapping_rejects_wildcard_mode_before_file_processing(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Wildcard"):
        run_module3_mapping(
            positive_file=tmp_path / "positive.csv",
            negative_file=tmp_path / "negative.csv",
            fasta_file=tmp_path / "proteome.fasta",
            output_dir=tmp_path / "output",
            output_folder_name="output",
            top_n=None,
            wildcards=True,
            q_cutoff=0.01,
        )
