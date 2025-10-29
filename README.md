# Antibody Analyzer

A full-stack web application for processing antibody sequencing datasets. The backend exposes the original FASTQ processing and k-mer enrichment modules through a FastAPI service, while the frontend provides a responsive dashboard for uploading data, monitoring progress, and downloading result bundles.

## Features

- **FASTQ processing** – upload multiple FASTQ/FASTQ.GZ files, optionally subtract background reads, and generate ranked peptide matrices.
- **K-mer analysis** – split merged AD/NC spreadsheets, tile peptides across k-mer sizes, and run Mann–Whitney U tests with FDR correction.
- **Downloadable artifacts** – every analysis run stores Excel/CSV outputs in a reusable archive.
- **Responsive UI** – React interface optimized for desktop and tablet workflows.

## Project structure

```
backend/        FastAPI application wrapping the antibody modules
frontend/       React + Vite single-page application
src/            Original research modules (kept for reference)
```

## Getting started

### Backend

```bash
cd backend
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

The API is exposed at `http://localhost:8000/api`.

### Frontend

```bash
cd frontend
npm install
npm run dev
```

By default the UI expects the backend to run on `http://localhost:8000`. You can change the API target by defining `VITE_API_BASE_URL` in a `.env` file inside `frontend/`.

## Running the workflows

1. Open the frontend (`npm run dev`) and visit `http://localhost:5173`.
2. Navigate to **FASTQ Processor** to upload sequencing runs and optional background files.
3. Navigate to **K-mer Analysis** to upload the merged AD/NC spreadsheet and configure k-mer parameters.
4. After processing, download the generated archives from the UI; archives include all CSV/XLSX outputs alongside filtered peptide logs.

## Tests

The repository currently does not ship automated tests. Recommended manual checks:

- `npm run build` inside `frontend/`
- `uvicorn app.main:app --reload` inside `backend/`
