import { FormEvent, useMemo, useState } from "react";

import { apiClient } from "../api/client";
import { FastqResponse, FastqSampleSummary } from "../api/types";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";

export function FastqProcessingPage() {
  const [files, setFiles] = useState<File[]>([]);
  const [backgroundFiles, setBackgroundFiles] = useState<File[]>([]);
  const [outputName, setOutputName] = useState("sequence_matrix.xlsx");
  const [result, setResult] = useState<FastqResponse | null>(null);

  const { execute, loading, error } = useAsyncTask(async (formData: FormData) => {
    const response = await apiClient.post<FastqResponse>("/process-fastq", formData, {
      headers: { "Content-Type": "multipart/form-data" },
    });
    return response.data;
  });

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    if (files.length === 0) {
      alert("Please add at least one FASTQ or FASTQ.GZ file.");
      return;
    }
    const formData = new FormData();
    files.forEach((file) => formData.append("files", file));
    backgroundFiles.forEach((file) => formData.append("background_files", file));
    formData.append("output_name", outputName);

    const response = await execute(formData);
    if (response) {
      setResult(response);
    }
  };

  const summaryRows = useMemo(() => result?.summary ?? [], [result]);
  const formatSelectedFiles = (selectedFiles: File[]) => {
    if (selectedFiles.length === 0) {
      return "No files selected";
    }
    if (selectedFiles.length === 1) {
      return selectedFiles[0].name;
    }
    return `${selectedFiles.length} files selected: ${selectedFiles.map((file) => file.name).join(", ")}`;
  };

  const renderUploadCard = (
    inputId: string,
    label: string,
    helperText: string,
    accept: string,
    selectedFiles: File[],
    onFilesChange: (files: File[]) => void,
  ) => (
    <div className="upload-card fastq-upload-card">
      <div className="upload-card__top">
        <div>
          <span className="upload-card__label">{label}</span>
          <p className="upload-card__helper">{helperText}</p>
        </div>
        <label className="upload-trigger" htmlFor={inputId}>
          Choose files
        </label>
      </div>
      <input
        id={inputId}
        className="sr-only-file-input"
        type="file"
        accept={accept}
        multiple
        onChange={(event) => {
          const targetFiles = event.target.files;
          onFilesChange(targetFiles ? Array.from(targetFiles) : []);
        }}
      />
      <div className={`upload-file-display ${selectedFiles.length > 0 ? "upload-file-display--selected" : ""}`}>
        {formatSelectedFiles(selectedFiles)}
      </div>
    </div>
  );

  return (
    <div className="page-grid">
      <SectionCard
        title="Process FASTQ libraries"
        description="Upload sequencing runs to generate ranked peptide matrices with optional background subtraction."
      >
        <form className="form-grid" onSubmit={handleSubmit}>
          {renderUploadCard(
            "fastq-files",
            "FASTQ files",
            "Supports plain or gzipped FASTQ files. Multiple uploads allowed.",
            ".fastq,.fq,.fastq.gz,.fq.gz",
            files,
            setFiles,
          )}

          {renderUploadCard(
            "fastq-background-files",
            "Background files (optional)",
            "Upload one or more FASTQ, FASTQ.GZ, or TXT background files. They will be combined before filtering.",
            ".fastq,.fq,.fastq.gz,.fq.gz,.txt",
            backgroundFiles,
            setBackgroundFiles,
          )}

          <label className="form-field">
            <span>Output name</span>
            <input
              type="text"
              value={outputName}
              onChange={(event) => setOutputName(event.target.value)}
              placeholder="sequence_matrix.xlsx"
            />
          </label>

          <div className="form-actions">
            <button type="submit" className="primary-button" disabled={loading}>
              {loading ? "Processing..." : "Run analysis"}
            </button>
            {result && <DownloadButton resultId={result.result_id} />}
          </div>
        </form>
        {error && <StatusBanner tone="error" title="Failed to process FASTQ files" message={error} />}
        {result && !error && (
          <StatusBanner
            tone="success"
            title="FASTQ analysis complete"
            message={`Primary matrix saved as ${result.excel_filename}`}
          />
        )}
      </SectionCard>

      {summaryRows.length > 0 && result && (
        <SectionCard
          title="Sample overview"
          description="Per-sample counts after background filtering. Columns match the Excel output."
          actions={<DownloadButton resultId={result.result_id} label="Download full bundle" />}
        >
          <div className="table-wrapper">
            <table>
              <thead>
                <tr>
                  <th>Sample</th>
                  <th>Total reads</th>
                  <th>Unique peptides</th>
                  <th>Filtered (background)</th>
                  <th>Output columns</th>
                </tr>
              </thead>
              <tbody>
                {summaryRows.map((row: FastqSampleSummary) => (
                  <tr key={row.sample_name}>
                    <td>{row.sample_name}</td>
                    <td>{row.total_sequences.toLocaleString()}</td>
                    <td>{row.unique_sequences.toLocaleString()}</td>
                    <td>{row.filtered_sequences.toLocaleString()}</td>
                    <td>{row.output_columns.join(", ")}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </SectionCard>
      )}
    </div>
  );
}
