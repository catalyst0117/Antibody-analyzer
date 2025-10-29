from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from app.core.fastq_processing import FastqProcessor
from app.core.kmer_analysis import analyze_groups, split_input_by_group
from app.models.responses import (
    FastqResponse,
    FastqSampleSummary,
    KmerResponse,
    KmerResultSummary,
)
from app.utils.result_store import ResultStore

router = APIRouter(prefix="/api")


def get_store() -> ResultStore:
    # Singleton pattern using function attribute
    if not hasattr(get_store, "_store"):
        get_store._store = ResultStore()
    return get_store._store  # type: ignore[attr-defined]


@router.post("/process-fastq", response_model=FastqResponse)
async def process_fastq(
    files: List[UploadFile] = File(...),
    background_file: UploadFile | None = File(None),
    output_name: str = Form("sequence_matrix.xlsx"),
    store: ResultStore = Depends(get_store),
) -> FastqResponse:
    if not files:
        raise HTTPException(status_code=400, detail="At least one FASTQ file must be provided.")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        saved_files: List[Path] = []
        for index, upload in enumerate(files):
            filename = upload.filename or f"sample_{index}.fastq"
            destination = tmp_path / filename
            destination.write_bytes(await upload.read())
            saved_files.append(destination)

        bg_path = None
        if background_file:
            bg_filename = background_file.filename or "background.fastq"
            bg_path = tmp_path / bg_filename
            bg_path.write_bytes(await background_file.read())

        processor = FastqProcessor(workdir=tmp_path / "outputs")
        result = processor.process(saved_files, background_file=bg_path, output_name=output_name)

        summary_payload = [
            FastqSampleSummary(
                sample_name=item.sample_name,
                total_sequences=item.total_sequences,
                unique_sequences=item.unique_sequences,
                filtered_sequences=item.filtered_sequences,
                output_columns=item.output_columns,
            )
            for item in result.summary
        ]

        files_to_archive = [result.excel_path]
        if result.background_dump:
            files_to_archive.append(result.background_dump)
        files_to_archive.extend(result.filtered_outputs.values())

        result_id = store.create_result(
            summary={
                "type": "fastq",
                "samples": [summary.dict() for summary in summary_payload],
            },
            files=files_to_archive,
        )

        return FastqResponse(
            result_id=result_id,
            summary=summary_payload,
            excel_filename=result.excel_path.name,
            background_filename=result.background_dump.name if result.background_dump else None,
            filtered_files=[path.name for path in result.filtered_outputs.values()],
        )


@router.post("/analyze-kmers", response_model=KmerResponse)
async def analyze_kmers(
    data_file: UploadFile = File(...),
    k_min: int = Form(4),
    k_max: int = Form(7),
    wildcard_positions: str = Form(""),
    normalize: bool = Form(True),
    store: ResultStore = Depends(get_store),
) -> KmerResponse:
    if k_min < 4 or k_max < k_min:
        raise HTTPException(status_code=400, detail="Invalid k-mer range.")

    wildcard_positions_list = [
        int(pos.strip())
        for pos in wildcard_positions.split(",")
        if pos.strip()
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        data_filename = data_file.filename or "patient_data.xlsx"
        data_path = tmp_path / data_filename
        data_path.write_bytes(await data_file.read())

        pos_file, neg_file = split_input_by_group(data_path)
        results = analyze_groups(
            data_path,
            pos_file,
            neg_file,
            k_values=range(k_min, k_max + 1),
            wildcard_positions=wildcard_positions_list,
            normalize=normalize,
            workdir=tmp_path / "outputs",
        )

        response_runs = [
            KmerResultSummary(
                k=item.k,
                total_kmers=item.total_kmers,
                ad_elevated=item.ad_elevated,
                nc_elevated=item.nc_elevated,
                result_filename=item.result_file.name,
                ad_filename=item.ad_file.name,
                nc_filename=item.nc_file.name,
                matrix_filename=item.matrix_file.name,
            )
            for item in results
        ]

        files_to_archive = [pos_file, neg_file]
        for item in results:
            files_to_archive.extend([item.result_file, item.ad_file, item.nc_file, item.matrix_file])

        result_id = store.create_result(
            summary={
                "type": "kmer",
                "runs": [run.dict() for run in response_runs],
            },
            files=files_to_archive,
        )

        return KmerResponse(result_id=result_id, runs=response_runs)


@router.get("/results/{result_id}")
async def get_result_summary(result_id: str, store: ResultStore = Depends(get_store)):
    summary = store.get_summary(result_id)
    if not summary:
        raise HTTPException(status_code=404, detail="Result not found")
    return summary


@router.get("/results/{result_id}/download")
async def download_result(result_id: str, store: ResultStore = Depends(get_store)):
    archive_path = store.get_archive_path(result_id)
    if not archive_path or not archive_path.exists():
        raise HTTPException(status_code=404, detail="Result not found")
    return FileResponse(
        path=archive_path,
        filename=archive_path.name,
        media_type="application/zip",
    )
