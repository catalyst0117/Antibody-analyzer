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
    volcano_filename: str


class KmerResponse(BaseModel):
    result_id: str
    runs: List[KmerResultSummary]


class KmerTaskCreatedResponse(BaseModel):
    task_id: str


class KmerTaskStatusResponse(BaseModel):
    task_id: str
    status: str
    progress: int
    message: str
    result: Optional[KmerResponse] = None
    error: Optional[str] = None


class Module3Response(BaseModel):
    result_id: str
    output_folder_name: str
    positive_mapping_filename: str
    negative_mapping_filename: str
    positive_clean_filename: str
    negative_clean_filename: str
    positive_manhattan_filename: str
    negative_manhattan_filename: str
    top_n: Optional[int]
    wildcards: bool
    q_cutoff: float


__all__ = [
    "FastqSampleSummary",
    "FastqResponse",
    "KmerResultSummary",
    "KmerResponse",
    "KmerTaskCreatedResponse",
    "KmerTaskStatusResponse",
    "Module3Response",
]
