from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel


class FastqSampleSummary(BaseModel):
    sample_name: str
    total_sequences: int
    unique_sequences: int
    filtered_sequences: int
    output_columns: List[str]


class FastqResponse(BaseModel):
    result_id: str
    summary: List[FastqSampleSummary]
    excel_filename: str
    background_filename: Optional[str]
    filtered_files: List[str]


class KmerResultSummary(BaseModel):
    k: int
    total_kmers: int
    ad_elevated: int
    nc_elevated: int
    result_filename: str
    ad_filename: str
    nc_filename: str
    matrix_filename: str


class KmerResponse(BaseModel):
    result_id: str
    runs: List[KmerResultSummary]


__all__ = [
    "FastqSampleSummary",
    "FastqResponse",
    "KmerResultSummary",
    "KmerResponse",
]
