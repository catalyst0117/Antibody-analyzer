# Antibody Analyzer

Antibody Analyzer is a full-stack web app for antibody sequencing workflows. It wraps the research modules in `src/` with a FastAPI backend and a React/Vite frontend so users can upload datasets, monitor run progress, and download reproducible result bundles from the browser.

## Current workflows

- **FASTQ Processor**: upload one or more FASTQ or FASTQ.GZ files, optionally subtract a background library, and generate a ranked peptide matrix.
- **K-mer Analysis**: upload a merged AD/NC spreadsheet, run Mann-Whitney U testing for a single k value at a time, track progress in the UI, and download the generated CSV bundle.
- **Proteome Mapping**: upload positive and negative significance files plus a local proteome FASTA, map significant k-mers onto the proteome, and download a bundled output folder with a run summary.

## Highlights

- FastAPI backend with downloadable ZIP bundles for every completed run
- React single-page app with section-based navigation
- Session-level UI persistence for recent form state and latest results when switching sections
- Friendly download naming for generated archives
- Vercel-compatible Python entrypoint at `api/main.py`

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
- optional background FASTQ
- output spreadsheet name

Outputs:
- primary Excel matrix
- filtered peptide files
- optional background dump
- downloadable ZIP bundle from the UI

### K-mer Analysis

Inputs:
- merged AD/NC spreadsheet
- single k value
- optional wildcard positions
- optional custom download bundle name
- normalize-counts toggle

Behavior:
- runs as a background task
- progress is exposed through the API and shown in the React UI
- preserves latest visible state when navigating away and back in the same browser session

Outputs:
- split cohort files
- Mann-Whitney summary CSV
- AD-enriched CSV
- NC-enriched CSV
- matrix CSV
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
  Downloads the stored ZIP bundle for a completed run.

## Deployment notes

- The frontend is a Vite app and can be deployed separately or alongside the backend.
- The FastAPI deployment entrypoint for platforms like Vercel is [`api/main.py`](api/main.py), which re-exports the app from `backend/app/main.py`.
- If you deploy the frontend separately from the backend, set `VITE_API_BASE_URL` to the deployed API base URL.

## Verification

The repository does not currently include a formal automated test suite. Recommended checks:

```bash
cd frontend
npm run build
```

```bash
cd backend
source .venv/bin/activate
python -m uvicorn app.main:app --reload
```
