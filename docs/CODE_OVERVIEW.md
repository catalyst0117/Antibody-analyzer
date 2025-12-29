# Code Overview: Antibody Analyzer

This document explains how the current backend, frontend, and core processing modules fit together so you can quickly orient yourself in the codebase.

## High-level architecture

```
Browser (React)  --->  FastAPI API  --->  Core processing modules (FASTQ + k-mer)
     |                     |                      |
     |                 ResultStore            Pandas / NumPy / SciPy
     |                     |
     |             ZIP archives for download
```

- **Frontend** (`frontend/`) provides the user interface for uploading files, selecting parameters, and downloading results.
- **Backend** (`backend/`) exposes two main workflows via HTTP endpoints: FASTQ processing and k-mer enrichment.
- **Core modules** (`backend/app/core/`) adapt the original desktop workflows (also mirrored in `src/`) into server-friendly functions.
- **Result storage** (`backend/app/utils/result_store.py`) archives outputs into ZIP files so the UI can download one bundle per run.

## Backend walkthrough (FastAPI)

### Entry point

- **`backend/app/main.py`** creates the FastAPI app, sets CORS, and includes the API router from `backend/app/api/routes.py`.

### API routes

File: **`backend/app/api/routes.py`**

- `POST /api/process-fastq`
  - Accepts multiple FASTQ / FASTQ.GZ files.
  - Optionally accepts a background FASTQ file to filter peptides.
  - Calls `FastqProcessor.process(...)`.
  - Stores outputs in `ResultStore` and returns a `FastqResponse`.

- `POST /api/analyze-kmers`
  - Accepts one merged AD/NC spreadsheet (XLSX or CSV).
  - Accepts `k_min`, `k_max`, and optional wildcard positions.
  - Calls `split_input_by_group(...)` and `analyze_groups(...)`.
  - Stores outputs in `ResultStore` and returns a `KmerResponse`.

- `GET /api/results/{result_id}`
  - Returns the summary JSON for a previously processed run.

- `GET /api/results/{result_id}/download`
  - Returns a ZIP archive containing all generated outputs for that run.

### Response models

File: **`backend/app/models/responses.py`**

- Defines `FastqResponse`, `KmerResponse`, and the per-sample summary objects returned by the API.

### Result storage

File: **`backend/app/utils/result_store.py`**

- `ResultStore.create_result(...)`
  - Copies outputs into a working directory.
  - Zips artifacts into a single archive.
  - Returns a unique `result_id` used by the API.
- `RESULT_ARCHIVE_DIR` environment variable lets you choose where archives are stored.

## FASTQ processing workflow

File: **`backend/app/core/fastq_processing.py`**

Key pieces:

- `FastqProcessor.process(...)`
  - Iterates each FASTQ file, searches reads for the motif `TTCTCACTCT.{0,36}`.
  - Extracts the 36-nt region after the motif and translates to amino acids.
  - Filters sequences that appear in the background FASTQ (if provided).
  - Builds a per-sample count matrix (sequence + count columns).
  - Writes a single Excel file combining all sample columns.
  - Writes optional filtered-out sequence logs per sample and background peptide dump.

- `SampleSummary` and `FastqResult`
  - Small dataclasses that describe totals and output paths.

## K-mer enrichment workflow

File: **`backend/app/core/kmer_analysis.py`**

Key pieces:

- `split_input_by_group(...)`
  - Reads the merged AD/NC spreadsheet.
  - Splits out two files: `*_AD` and `*_NC`, using column prefixes.

- `tile_patient_file(...)`
  - Slides a k-mer window over each peptide sequence.
  - Optionally applies wildcard positions (replacing positions with `X`).
  - Aggregates counts per patient and globally.
  - Uses a chi-square filter to drop non-enriched k-mers.

- `build_kmer_matrix(...)`
  - Builds a k-mer-by-patient matrix with optional normalization.

- `run_mannwhitney(...)`
  - Runs Mann–Whitney U tests comparing AD vs NC columns.
  - Writes:
    - `_U_test_*.csv` (p-values + mean rank differences)
    - `_AD.csv` (AD-elevated k-mers)
    - `_NC.csv` (NC-elevated k-mers)

- `analyze_groups(...)`
  - Orchestrates the tiling, filtering, matrix creation, and test run for each k value.

## Frontend walkthrough (React)

### Entry points

- **`frontend/src/main.tsx`** bootstraps React, React Router, and global CSS.
- **`frontend/src/App.tsx`** defines two routes:
  - `/fastq` (FASTQ Processing)
  - `/kmer` (K-mer Analysis)

### Pages

- **`frontend/src/pages/FastqProcessingPage.tsx`**
  - Uploads multiple FASTQ files and optional background file.
  - Posts to `/api/process-fastq` using `frontend/src/api/client.ts`.
  - Renders response summaries and download link.

- **`frontend/src/pages/KmerAnalysisPage.tsx`**
  - Uploads a merged AD/NC spreadsheet.
  - Allows setting `k_min`, `k_max`, wildcard positions, and normalization.
  - Posts to `/api/analyze-kmers` and lists each k-mer run summary.

### Shared UI components

- **`frontend/src/components/HeaderNav.tsx`** – top navigation for FASTQ vs k-mer screens.
- **`frontend/src/components/SectionCard.tsx`** – layout card wrapper.
- **`frontend/src/components/StatusBanner.tsx`** – success/error display.
- **`frontend/src/components/DownloadButton.tsx`** – triggers ZIP download for result IDs.

### Data access + state helpers

- **`frontend/src/api/client.ts`**
  - Central wrapper for API calls and base URL configuration.

- **`frontend/src/hooks/useAsyncTask.ts`**
  - Shared async state handler (`isLoading`, `error`, `data`).

## Original research code (reference)

The legacy desktop scripts live in **`src/`** and document the original algorithmic steps. The FastAPI service mirrors these workflows but wraps them for web requests.

## Typical end-to-end flow

1. User uploads data in the React UI.
2. React posts to the FastAPI endpoints in `backend/app/api/routes.py`.
3. FastAPI calls the core processing functions in `backend/app/core/`.
4. Outputs are archived in `ResultStore` and returned to the UI.
5. The UI downloads the ZIP bundle for analysis results.
