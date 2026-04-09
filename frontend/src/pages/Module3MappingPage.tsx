import { FormEvent, useEffect, useMemo, useState } from "react";

import { apiClient } from "../api/client";
import { Module3Response } from "../api/types";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";

const PROTEOME_MAPPING_SESSION_KEY = "proteome-mapping-page-state";

type ProteomeMappingSessionState = {
  positiveFileName: string | null;
  negativeFileName: string | null;
  proteomeFileName: string | null;
  outputFolderName: string;
  topNValue: string;
  qCutoffValue: string;
  wildcards: boolean;
  result: Module3Response | null;
};

function loadProteomeMappingSessionState(): ProteomeMappingSessionState | null {
  if (typeof window === "undefined") {
    return null;
  }
  const raw = window.sessionStorage.getItem(PROTEOME_MAPPING_SESSION_KEY);
  if (!raw) {
    return null;
  }
  try {
    return JSON.parse(raw) as ProteomeMappingSessionState;
  } catch {
    return null;
  }
}

export function Module3MappingPage() {
  const savedState = loadProteomeMappingSessionState();
  const [positiveFile, setPositiveFile] = useState<File | null>(null);
  const [negativeFile, setNegativeFile] = useState<File | null>(null);
  const [proteomeFile, setProteomeFile] = useState<File | null>(null);
  const [positiveFileName, setPositiveFileName] = useState<string | null>(savedState?.positiveFileName ?? null);
  const [negativeFileName, setNegativeFileName] = useState<string | null>(savedState?.negativeFileName ?? null);
  const [proteomeFileName, setProteomeFileName] = useState<string | null>(savedState?.proteomeFileName ?? null);
  const [outputFolderName, setOutputFolderName] = useState(savedState?.outputFolderName ?? "");
  const [topNValue, setTopNValue] = useState(savedState?.topNValue ?? "");
  const [qCutoffValue, setQCutoffValue] = useState(savedState?.qCutoffValue ?? "0.01");
  const [wildcards, setWildcards] = useState(savedState?.wildcards ?? false);
  const [result, setResult] = useState<Module3Response | null>(savedState?.result ?? null);

  const { execute, loading, error } = useAsyncTask(async (formData: FormData) => {
    const response = await apiClient.post<Module3Response>("/module3-map", formData, {
      headers: { "Content-Type": "multipart/form-data" },
    });
    return response.data;
  });

  useEffect(() => {
    if (typeof window === "undefined") {
      return;
    }
    const nextState: ProteomeMappingSessionState = {
      outputFolderName,
      positiveFileName,
      negativeFileName,
      proteomeFileName,
      topNValue,
      qCutoffValue,
      wildcards,
      result,
    };
    window.sessionStorage.setItem(PROTEOME_MAPPING_SESSION_KEY, JSON.stringify(nextState));
  }, [negativeFileName, outputFolderName, positiveFileName, proteomeFileName, qCutoffValue, result, topNValue, wildcards]);

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();

    if (!positiveFile || !negativeFile || !proteomeFile) {
      alert("Please upload the positive file, negative file, and proteome FASTA.");
      return;
    }

    const parsedTopN = topNValue.trim() === "" ? null : Number.parseInt(topNValue, 10);
    if (parsedTopN !== null && (!Number.isInteger(parsedTopN) || parsedTopN <= 0)) {
      alert("top_n must be a positive integer.");
      return;
    }

    const parsedQCutoff = Number.parseFloat(qCutoffValue);
    if (Number.isNaN(parsedQCutoff) || parsedQCutoff < 0 || parsedQCutoff > 1) {
      alert("q_cutoff must be a number between 0 and 1.");
      return;
    }

    const formData = new FormData();
    formData.append("positive_file", positiveFile);
    formData.append("negative_file", negativeFile);
    formData.append("proteome_fasta", proteomeFile);
    formData.append("output_folder_name", outputFolderName);
    if (parsedTopN !== null) {
      formData.append("top_n", String(parsedTopN));
    }
    formData.append("wildcards", String(wildcards));
    formData.append("q_cutoff", String(parsedQCutoff));

    const response = await execute(formData);
    if (response) {
      setResult(response);
    }
  };

  const outputFiles = useMemo(() => {
    if (!result) {
      return [];
    }

      return [
        result.positive_mapping_filename,
        result.negative_mapping_filename,
        result.positive_clean_filename,
        result.negative_clean_filename,
        result.positive_manhattan_filename,
        result.negative_manhattan_filename,
        "run_summary.txt",
      ];
  }, [result]);

  return (
    <div className="page-grid proteome-page">
      <SectionCard
        title="Proteome Mapping"
        description="Map significant positive and negative k-mers onto a local proteome FASTA and download the bundled mapping outputs."
      >
        <form className="proteome-form" onSubmit={handleSubmit}>
          <section className="proteome-section">
            <div className="proteome-section__header">
              <h3>Input Files</h3>
              <p>Upload the significance tables and the reference FASTA needed for one mapping run.</p>
            </div>

            <div className="upload-card-grid">
              <div className="upload-card">
                <div className="upload-card__top">
                  <div>
                    <span className="upload-card__label">Positive significance file</span>
                    <p className="upload-card__helper">
                      Upload the positive significance output from Module 2 or another compatible table.
                    </p>
                  </div>
                  <label className="upload-trigger" htmlFor="proteome-positive-file">
                    Choose file
                  </label>
                </div>
                <input
                  id="proteome-positive-file"
                  className="sr-only-file-input"
                  type="file"
                  accept=".xlsx,.xls,.csv,.txt"
                  onChange={(event) => {
                    const files = event.target.files;
                    const selectedFile = files && files.length > 0 ? files[0] : null;
                    setPositiveFile(selectedFile);
                    setPositiveFileName(selectedFile?.name ?? null);
                  }}
                />
                <div className={`upload-file-display ${positiveFileName ? "upload-file-display--selected" : ""}`}>
                  {positiveFileName ?? "No file selected"}
                </div>
              </div>

              <div className="upload-card">
                <div className="upload-card__top">
                  <div>
                    <span className="upload-card__label">Negative significance file</span>
                    <p className="upload-card__helper">
                      Upload the matching negative significance file with the same k-mer length.
                    </p>
                  </div>
                  <label className="upload-trigger" htmlFor="proteome-negative-file">
                    Choose file
                  </label>
                </div>
                <input
                  id="proteome-negative-file"
                  className="sr-only-file-input"
                  type="file"
                  accept=".xlsx,.xls,.csv,.txt"
                  onChange={(event) => {
                    const files = event.target.files;
                    const selectedFile = files && files.length > 0 ? files[0] : null;
                    setNegativeFile(selectedFile);
                    setNegativeFileName(selectedFile?.name ?? null);
                  }}
                />
                <div className={`upload-file-display ${negativeFileName ? "upload-file-display--selected" : ""}`}>
                  {negativeFileName ?? "No file selected"}
                </div>
              </div>

              <div className="upload-card">
                <div className="upload-card__top">
                  <div>
                    <span className="upload-card__label">Proteome FASTA file</span>
                    <p className="upload-card__helper">
                      Local FASTA only. Proteome Mapping will build the required k-mer index from this file.
                    </p>
                  </div>
                  <label className="upload-trigger" htmlFor="proteome-fasta-file">
                    Choose file
                  </label>
                </div>
                <input
                  id="proteome-fasta-file"
                  className="sr-only-file-input"
                  type="file"
                  accept=".fasta,.fa,.faa,.txt"
                  onChange={(event) => {
                    const files = event.target.files;
                    const selectedFile = files && files.length > 0 ? files[0] : null;
                    setProteomeFile(selectedFile);
                    setProteomeFileName(selectedFile?.name ?? null);
                  }}
                />
                <div className={`upload-file-display ${proteomeFileName ? "upload-file-display--selected" : ""}`}>
                  {proteomeFileName ?? "No file selected"}
                </div>
              </div>
            </div>
          </section>

          <section className="proteome-section">
            <div className="proteome-section__header">
              <h3>Run Settings</h3>
              <p>Set the output naming and filtering options for this mapping job.</p>
            </div>

            <label className="form-field proteome-full-width-field">
              <span>Output folder name</span>
              <input
                type="text"
                value={outputFolderName}
                onChange={(event) => setOutputFolderName(event.target.value)}
                placeholder="Optional. Leave blank to auto-generate"
              />
              <small>
                Optional. If left blank, a dataset name will be generated from the uploaded significance file and top_n.
              </small>
            </label>

            <div className="settings-grid">
              <label className="form-field">
                <span>top_n</span>
                <input
                  type="text"
                  inputMode="numeric"
                  pattern="[0-9]*"
                  value={topNValue}
                  onChange={(event) => setTopNValue(event.target.value)}
                  placeholder="Leave empty for all"
                />
                <small>Optional. Keep only the top N lowest-q rows from each file.</small>
              </label>
              <label className="form-field">
                <span>q_cutoff</span>
                <input
                  type="text"
                  value={qCutoffValue}
                  onChange={(event) => setQCutoffValue(event.target.value)}
                  placeholder="0.01"
                />
                <small>Rows above this q-value are excluded from the proteome mapping step.</small>
              </label>
            </div>

            <label className="form-field form-checkbox settings-toggle">
              <input
                type="checkbox"
                checked={wildcards}
                onChange={(event) => setWildcards(event.target.checked)}
              />
              <span>Enable wildcard X matching</span>
            </label>
          </section>

          <div className="proteome-actions">
            <button type="submit" className="primary-button" disabled={loading}>
              {loading ? "Running Proteome Mapping..." : "Start Proteome Mapping"}
            </button>
            {result && <DownloadButton resultId={result.result_id} label="Download Proteome Mapping bundle" />}
          </div>
        </form>

        {error && <StatusBanner tone="error" title="Proteome Mapping failed" message={error} />}
        {result && !error && (
          <StatusBanner
            tone="success"
            title="Proteome Mapping complete"
            message="The positive and negative proteome mapping outputs are ready to download."
          />
        )}
      </SectionCard>

      {result && outputFiles.length > 0 && (
        <SectionCard
          title="Proteome Mapping Outputs"
          description="These files are included in the downloadable ZIP bundle."
          actions={<DownloadButton resultId={result.result_id} label="Download combined outputs" />}
        >
          <div className="table-wrapper">
            <table>
              <thead>
                <tr>
                  <th>Setting</th>
                  <th>Value</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>top_n</td>
                  <td>{result.top_n ?? "All rows"}</td>
                </tr>
                <tr>
                  <td>Output folder name</td>
                  <td>{result.output_folder_name}</td>
                </tr>
                <tr>
                  <td>Wildcard matching</td>
                  <td>{result.wildcards ? "Enabled" : "Disabled"}</td>
                </tr>
                <tr>
                  <td>q_cutoff</td>
                  <td>{result.q_cutoff}</td>
                </tr>
                <tr>
                  <td>Output files</td>
                  <td>{outputFiles.join(", ")}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </SectionCard>
      )}
    </div>
  );
}
