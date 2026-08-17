# Antibody Analyzer

Antibody Analyzer is a full-stack web app for antibody sequencing workflows. It wraps the research modules in `src/` with a FastAPI backend and a React/Vite frontend so users can upload datasets, monitor run progress, and download reproducible result bundles from the browser.

## Current workflows

- **FASTQ Processor**: upload one or more FASTQ or FASTQ.GZ files, optionally subtract one or more background files, and generate a ranked peptide matrix.
- **K-mer Analysis**: upload either a merged cohort spreadsheet with configurable positive/negative keywords or already-separated positive and negative cohort files, run Mann-Whitney U testing for one k value, track progress in the UI, and download the generated CSV bundle.
- **Proteome Mapping**: upload positive and negative significance files plus a local proteome FASTA, or open it from a completed K-mer Analysis run with the generated U-test files attached automatically.

## Highlights

- FastAPI backend with downloadable ZIP bundles for every completed run
- React single-page app with section-based navigation
- Session-level UI persistence for recent form state and latest results when switching sections
- Friendly download naming for generated archives
- Vercel-compatible Python entrypoint at `api/main.py`
- [Matrix comparison utility](#matrix-comparison-utility)

## Project structure

```text
api/            Deployment entrypoint for FastAPI on Vercel
backend/        FastAPI application and workflow orchestration
frontend/       React + Vite single-page application
src/            Original research modules kept alongside the app
requirements.txt
run_project_windows.bat
```

## Local development

### Prerequisites

- Python 3.10+
- Node.js 18+

### Backend

From the repository root:

```bash
python -m venv backend/.venv
source backend/.venv/bin/activate
pip install -r requirements.txt
cd backend
uvicorn app.main:app --reload
```

The API will be available at `http://localhost:8000/api`.

### Frontend

```bash
cd frontend
npm install
npm run dev
```

The development UI runs at `http://localhost:5173`.

By default the frontend talks to `http://localhost:8000/api`. To point it elsewhere, set `VITE_API_BASE_URL` in `frontend/.env`.

### Windows helper

For local development on Windows, you can use:

```bat
run_project_windows.bat
```

That script opens one terminal for the backend and one for the frontend.

## How each workflow behaves

### FASTQ Processor

Inputs:
- one or more sequencing files
- optional background files in FASTQ, FASTQ.GZ, or TXT format
- output spreadsheet name

Outputs:
- primary Excel matrix
- filtered peptide files
- optional background dump
- downloadable ZIP bundle from the UI

### K-mer Analysis

Inputs:
- merged cohort spreadsheet with user-selected positive and negative column keywords
- or already-separated positive and negative cohort files
- single k value
- optional wildcard positions
- optional maximum percentage of zero-valued samples for statistical output
- optional custom download bundle name
- normalize-counts toggle

Behavior:
- runs as a background task
- progress is exposed through the API and shown in the React UI
- builds the complete k-mer universe and applies product-based chi-square
  filtering independently to each sample
- always uses pre-filter sample totals for normalized Mann-Whitney inputs
- optionally writes the downloadable matrix as raw counts or normalized values
- can omit sparse k-mers from Mann-Whitney and result CSVs without changing the
  downloadable matrix
- generated U-test positive/negative files can be handed directly to Proteome Mapping
- preserves latest visible state when navigating away and back in the same browser session

Outputs:
- split cohort files
- Mann-Whitney summary CSV
- positive-elevated CSV
- negative-elevated CSV
- matrix CSV with the k-mer column included first
- interactive volcano plot HTML with hover details, zoom, pan, and an
  adjustable Q-value cutoff
- downloadable ZIP bundle

### Proteome Mapping

Inputs:
- positive significance file
- negative significance file
- proteome FASTA
- optional output folder name
- optional `top_n`
- wildcard toggle
- `q_cutoff`

Behavior:
- uses local FASTA mapping only
- can auto-attach positive and negative U-test CSVs from a completed K-mer Analysis run
- still allows users to manually replace any carried-over files
- generates a default output folder name when none is supplied
- keeps the latest result panel visible when revisiting the section in the same browser session

Outputs:
- positive mapping CSV
- negative mapping CSV
- cleaned positive and negative text outputs
- Manhattan JSON files
- `run_summary.txt`
- downloadable ZIP bundle rooted at the resolved output folder name

## API summary

- `POST /api/process-fastq`
  Processes FASTQ uploads and returns a `result_id` for bundle download.
- `POST /api/analyze-kmers`
  Starts a background k-mer task and returns a `task_id`.
- `GET /api/analyze-kmers/{task_id}`
  Returns task status, progress, and final result when complete.
- `POST /api/module3-map`
  Runs Proteome Mapping and returns filenames plus a `result_id`.
- `GET /api/results/{result_id}`
  Returns stored result metadata.
- `GET /api/results/{result_id}/files/{filename}`
  Downloads one generated file from a stored result bundle, used by the K-mer to Proteome Mapping handoff.
- `GET /api/results/{result_id}/download`
  Downloads the stored ZIP bundle for a completed run.

## Matrix comparison utility

The comparator accepts `.csv` and `.xlsx` matrices in any combination. For Excel
files, it reads the active worksheet and uses calculated cell values. Pass the
approved or legacy matrix first and the newly generated matrix second:

```bash
python3 -m tests.matrix_compare \
  path/to/expected_matrix.csv \
  path/to/generated_matrix.xlsx \
  --output-dir test-artifacts/module2-matrix
```

The command checks rows and sample columns by name, compares numeric values with a
small floating-point tolerance, and writes:

- `test-artifacts/module2-matrix/matrix_diff.csv`: `actual - expected` for every
  shared cell
- `test-artifacts/module2-matrix/matrix_diff.html`: a visual heatmap where blue is
  lower, red is higher, and gray marks missing rows or columns

Open `matrix_diff.html` in a browser to inspect the largest differences. The command
exits with status `0` when the matrices match and `1` when they differ, so status
`1` means the report was generated successfully but differences were found.

## Deployment notes

- The frontend is a Vite app and can be deployed separately or alongside the backend.
- The FastAPI deployment entrypoint for platforms like Vercel is [`api/main.py`](api/main.py), which re-exports the app from `backend/app/main.py`.
- If you deploy the frontend separately from the backend, set `VITE_API_BASE_URL` to the deployed API base URL.

## Verification

Recommended application checks:

```bash
cd frontend
npm run build
```

```bash
cd backend
source .venv/bin/activate
python -m uvicorn app.main:app --reload
```
