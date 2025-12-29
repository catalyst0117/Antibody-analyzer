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
docs/           Extended documentation and code overview
```

## Documentation

- [Code overview](docs/CODE_OVERVIEW.md) for a guided tour of the backend, frontend, and processing modules.

## Getting started

### Prerequisites

- Python 3.10+
- Node.js 18+ with npm
- (Optional) [UVicorn](https://www.uvicorn.org/) and [FastAPI](https://fastapi.tiangolo.com/) knowledge if you plan to tweak the server

### Quick start

1. **Clone the repository** (or download the source bundle) and open a terminal in the project root.
2. **Install backend dependencies** and start the API server:

   ```bash
   cd backend
   python -m venv .venv
   source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
   pip install -r requirements.txt
   uvicorn app.main:app --reload
   ```

   The API becomes available at `http://localhost:8000/api` (the interactive docs are at `http://localhost:8000/docs`).

3. **Install frontend dependencies** and start the development server in a second terminal:

   ```bash
   cd frontend
   npm install
   npm run dev
   ```

   The Vite dev server exposes the UI at `http://localhost:5173`. The app automatically proxies requests to `http://localhost:8000`.

### Configuration notes

- To point the frontend at a different backend URL, create a `.env` file inside `frontend/` containing `VITE_API_BASE_URL="https://your-api-host"`.
- The backend stores generated archives in a temporary directory. Set the environment variable `RESULT_ARCHIVE_DIR` before launching `uvicorn` to persist archives in a custom location.
- When deploying to production, run `npm run build` in `frontend/` and serve the generated `dist/` directory with your preferred web server.

## Running the workflows

1. Open the frontend (`npm run dev`) and visit `http://localhost:5173`.
2. Navigate to **FASTQ Processor** to upload sequencing runs and optional background files.
3. Navigate to **K-mer Analysis** to upload the merged AD/NC spreadsheet and configure k-mer parameters.
4. After processing, download the generated archives from the UI; archives include all CSV/XLSX outputs alongside filtered peptide logs.

## Tests

The repository currently does not ship automated tests. Recommended manual checks:

- `npm run build` inside `frontend/`
- `uvicorn app.main:app --reload` inside `backend/`
