from __future__ import annotations

import tempfile
import threading
from pathlib import Path
from typing import List

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from app.core.fastq_processing import FastqProcessor
from app.core.kmer_analysis import analyze_single_k, split_input_by_group
from app.core.module3_mapping import resolve_output_folder_name, run_module3_mapping
from app.models.responses import (
    FastqResponse,
    FastqSampleSummary,
    KmerResponse,
    KmerResultSummary,
    KmerTaskCreatedResponse,
    KmerTaskStatusResponse,
    Module3Response,
)
from app.utils.kmer_task_store import KmerTaskStore
from app.utils.result_store import ResultStore
from typing import Optional

router = APIRouter(prefix="/api")


def get_store() -> ResultStore:
    # Singleton pattern using function attribute
    if not hasattr(get_store, "_store"):
        get_store._store = ResultStore()
    return get_store._store  # type: ignore[attr-defined]


def get_kmer_task_store() -> KmerTaskStore:
    if not hasattr(get_kmer_task_store, "_store"):
        get_kmer_task_store._store = KmerTaskStore()
    return get_kmer_task_store._store  # type: ignore[attr-defined]


def _run_kmer_task(
    *,
    task_id: str,
    data_bytes: bytes,
    data_filename: str,
    k: int,
    wildcard_positions_list: List[int],
    normalize: bool,
    archive_name: str,
    store: ResultStore,
    task_store: KmerTaskStore,
) -> None:
    try:
        task_store.update_task(task_id, status="running", progress=5, message="Preparing analysis")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_path = tmp_path / data_filename
            data_path.write_bytes(data_bytes)

            task_store.update_task(task_id, progress=8, message="Splitting AD and NC cohorts")
            pos_file, neg_file = split_input_by_group(data_path)
            result = analyze_single_k(
                data_path,
                pos_file,
                neg_file,
                k=k,
                wildcard_positions=wildcard_positions_list,
                normalize=normalize,
                workdir=tmp_path / "outputs",
                progress_callback=lambda progress, message: task_store.update_task(
                    task_id,
                    status="running",
                    progress=progress,
                    message=message,
                ),
            )

            response_run = KmerResultSummary(
                k=result.k,
                total_kmers=result.total_kmers,
                ad_elevated=result.ad_elevated,
                nc_elevated=result.nc_elevated,
                result_filename=result.result_file.name,
                ad_filename=result.ad_file.name,
                nc_filename=result.nc_file.name,
                matrix_filename=result.matrix_file.name,
            )

            files_to_archive = [pos_file, neg_file, result.result_file, result.ad_file, result.nc_file, result.matrix_file]
            task_store.update_task(task_id, progress=98, message="Bundling result files")
            result_id = store.create_result(
                summary={
                    "type": "kmer",
                    "runs": [response_run.model_dump()],
                },
                files=files_to_archive,
                download_name=archive_name or f"{Path(data_filename).stem}_k{k}_results",
            )

            task_store.complete_task(task_id, KmerResponse(result_id=result_id, runs=[response_run]))
    except Exception as exc:
        task_store.fail_task(task_id, str(exc))


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
            download_name=f"{Path(output_name).stem}_results",
        )

        return FastqResponse(
            result_id=result_id,
            summary=summary_payload,
            excel_filename=result.excel_path.name,
            background_filename=result.background_dump.name if result.background_dump else None,
            filtered_files=[path.name for path in result.filtered_outputs.values()],
        )




@router.post("/analyze-kmers", response_model=KmerTaskCreatedResponse)
async def analyze_kmers(
    data_file: UploadFile = File(...),
    k: int = Form(4),
    wildcard_positions: str = Form(""),
    normalize: bool = Form(True),
    archive_name: str = Form(""),
    store: ResultStore = Depends(get_store),
    task_store: KmerTaskStore = Depends(get_kmer_task_store),
) -> KmerTaskCreatedResponse:
    if k < 4:
        raise HTTPException(status_code=400, detail="k must be >= 4.")
    wildcard_positions_list = [
        int(pos.strip())
        for pos in wildcard_positions.split(",")
        if pos.strip()
    ]
    data_bytes = await data_file.read()
    data_filename = data_file.filename or "patient_data.xlsx"

    task_id = task_store.create_task()
    worker = threading.Thread(
        target=_run_kmer_task,
        kwargs={
            "task_id": task_id,
            "data_bytes": data_bytes,
            "data_filename": data_filename,
            "k": k,
            "wildcard_positions_list": wildcard_positions_list,
            "normalize": normalize,
            "archive_name": archive_name,
            "store": store,
            "task_store": task_store,
        },
        daemon=True,
    )
    worker.start()

    return KmerTaskCreatedResponse(task_id=task_id)


@router.get("/analyze-kmers/{task_id}", response_model=KmerTaskStatusResponse)
async def get_kmer_task_status(
    task_id: str,
    task_store: KmerTaskStore = Depends(get_kmer_task_store),
) -> KmerTaskStatusResponse:
    task = task_store.get_task(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    return KmerTaskStatusResponse(
        task_id=task.task_id,
        status=task.status,
        progress=task.progress,
        message=task.message,
        result=task.result,
        error=task.error,
    )


@router.post("/module3-map", response_model=Module3Response)
async def module3_map(
    positive_file: UploadFile = File(...),
    negative_file: UploadFile = File(...),
    proteome_fasta: UploadFile = File(...),
    top_n: Optional[int] = Form(None),
    wildcards: bool = Form(False),
    q_cutoff: float = Form(0.01),
    output_folder_name: str = Form(""),
    store: ResultStore = Depends(get_store),
) -> Module3Response:
    if top_n is not None and top_n <= 0:
        raise HTTPException(status_code=400, detail="top_n must be greater than 0.")
    if q_cutoff < 0 or q_cutoff > 1:
        raise HTTPException(status_code=400, detail="q_cutoff must be between 0 and 1.")

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            positive_path = tmp_path / (positive_file.filename or "positive_significance.csv")
            negative_path = tmp_path / (negative_file.filename or "negative_significance.csv")
            fasta_path = tmp_path / (proteome_fasta.filename or "proteome.fasta")
            resolved_output_folder_name = resolve_output_folder_name(
                requested_name=output_folder_name,
                source_filename=positive_path.name,
                top_n=top_n,
            )

            positive_path.write_bytes(await positive_file.read())
            negative_path.write_bytes(await negative_file.read())
            fasta_path.write_bytes(await proteome_fasta.read())

            result = run_module3_mapping(
                positive_file=positive_path,
                negative_file=negative_path,
                fasta_file=fasta_path,
                output_dir=tmp_path / resolved_output_folder_name,
                output_folder_name=resolved_output_folder_name,
                top_n=top_n,
                wildcards=wildcards,
                q_cutoff=q_cutoff,
            )

            files_to_archive = [
                result.positive_mapping_file,
                result.negative_mapping_file,
                result.positive_clean_file,
                result.negative_clean_file,
                result.positive_manhattan_file,
                result.negative_manhattan_file,
                result.run_summary_file,
            ]
            result_id = store.create_result(
                summary={
                    "type": "proteome_mapping",
                    "output_folder_name": result.output_folder_name,
                    "positive_mapping_filename": result.positive_mapping_file.name,
                    "negative_mapping_filename": result.negative_mapping_file.name,
                    "positive_clean_filename": result.positive_clean_file.name,
                    "negative_clean_filename": result.negative_clean_file.name,
                    "positive_manhattan_filename": result.positive_manhattan_file.name,
                    "negative_manhattan_filename": result.negative_manhattan_file.name,
                    "run_summary_filename": result.run_summary_file.name,
                    "top_n": top_n,
                    "wildcards": wildcards,
                    "q_cutoff": q_cutoff,
                },
                files=files_to_archive,
                download_name=result.output_folder_name,
                archive_root_name=result.output_folder_name,
            )

            return Module3Response(
                result_id=result_id,
                output_folder_name=result.output_folder_name,
                positive_mapping_filename=result.positive_mapping_file.name,
                negative_mapping_filename=result.negative_mapping_file.name,
                positive_clean_filename=result.positive_clean_file.name,
                negative_clean_filename=result.negative_clean_file.name,
                positive_manhattan_filename=result.positive_manhattan_file.name,
                negative_manhattan_filename=result.negative_manhattan_file.name,
                top_n=top_n,
                wildcards=wildcards,
                q_cutoff=q_cutoff,
            )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get("/results/{result_id}")
async def get_result_summary(result_id: str, store: ResultStore = Depends(get_store)):
    summary = store.get_summary(result_id)
    if not summary:
        raise HTTPException(status_code=404, detail="Result not found")
    return summary


@router.get("/results/{result_id}/download")
async def download_result(result_id: str, store: ResultStore = Depends(get_store)):
    archive_path = store.get_archive_path(result_id)
    download_name = store.get_download_name(result_id)
    if not archive_path or not archive_path.exists():
        raise HTTPException(status_code=404, detail="Result not found")
    return FileResponse(
        path=archive_path,
        filename=download_name or archive_path.name,
        media_type="application/zip",
    )
