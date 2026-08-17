import { FormEvent, useEffect, useMemo, useState } from "react";
import { Link } from "react-router-dom";

import { apiClient } from "../api/client";
import {
  KmerResponse,
  KmerResultSummary,
  KmerTaskCreatedResponse,
  KmerTaskStatusResponse,
} from "../api/types";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";

const KMER_SESSION_KEY = "kmer-analysis-page-state";

type KmerInputMode = "merged" | "separate";

type KmerSessionState = {
  inputMode: KmerInputMode;
  dataFileName: string | null;
  positiveFileName: string | null;
  negativeFileName: string | null;
  positiveKeyword: string;
  negativeKeyword: string;
  kValue: string;
  archiveName: string;
  wildcards: string;
  normalize: boolean;
  result: KmerResponse | null;
  taskId: string | null;
  taskStatus: KmerTaskStatusResponse | null;
  loading: boolean;
  error: string | null;
};

function loadKmerSessionState(): KmerSessionState | null {
  if (typeof window === "undefined") {
    return null;
  }
  const raw = window.sessionStorage.getItem(KMER_SESSION_KEY);
  if (!raw) {
    return null;
  }
  try {
    return JSON.parse(raw) as KmerSessionState;
  } catch {
    return null;
  }
}

export function KmerAnalysisPage() {
  const savedState = loadKmerSessionState();
  const [inputMode, setInputMode] = useState<KmerInputMode>(savedState?.inputMode ?? "merged");
  const [dataFile, setDataFile] = useState<File | null>(null);
  const [positiveFile, setPositiveFile] = useState<File | null>(null);
  const [negativeFile, setNegativeFile] = useState<File | null>(null);
  const [dataFileName, setDataFileName] = useState<string | null>(savedState?.dataFileName ?? null);
  const [positiveFileName, setPositiveFileName] = useState<string | null>(savedState?.positiveFileName ?? null);
  const [negativeFileName, setNegativeFileName] = useState<string | null>(savedState?.negativeFileName ?? null);
  const [positiveKeyword, setPositiveKeyword] = useState(savedState?.positiveKeyword ?? "AD");
  const [negativeKeyword, setNegativeKeyword] = useState(savedState?.negativeKeyword ?? "NC");
  const [kValue, setKValue] = useState(savedState?.kValue ?? "4");
  const [archiveName, setArchiveName] = useState(savedState?.archiveName ?? "");
  const [wildcards, setWildcards] = useState(savedState?.wildcards ?? "");
  const [normalize, setNormalize] = useState(savedState?.normalize ?? true);
  const [result, setResult] = useState<KmerResponse | null>(savedState?.result ?? null);
  const [taskId, setTaskId] = useState<string | null>(savedState?.taskId ?? null);
  const [taskStatus, setTaskStatus] = useState<KmerTaskStatusResponse | null>(savedState?.taskStatus ?? null);
  const [loading, setLoading] = useState(savedState?.loading ?? false);
  const [error, setError] = useState<string | null>(savedState?.error ?? null);

  useEffect(() => {
    if (typeof window === "undefined") {
      return;
    }
    const nextState: KmerSessionState = {
      inputMode,
      kValue,
      dataFileName,
      positiveFileName,
      negativeFileName,
      positiveKeyword,
      negativeKeyword,
      archiveName,
      wildcards,
      normalize,
      result,
      taskId,
      taskStatus,
      loading,
      error,
    };
    window.sessionStorage.setItem(KMER_SESSION_KEY, JSON.stringify(nextState));
  }, [
    archiveName,
    dataFileName,
    error,
    inputMode,
    kValue,
    loading,
    negativeFileName,
    negativeKeyword,
    normalize,
    positiveFileName,
    positiveKeyword,
    result,
    taskId,
    taskStatus,
    wildcards,
  ]);

  useEffect(() => {
    if (!taskId) {
      return;
    }

    let cancelled = false;
    const pollStatus = async () => {
      try {
        const response = await apiClient.get<KmerTaskStatusResponse>(`/analyze-kmers/${taskId}`);
        if (cancelled) {
          return;
        }

        const nextStatus = response.data;
        setTaskStatus(nextStatus);

        if (nextStatus.status === "succeeded") {
          setLoading(false);
          setResult(nextStatus.result);
          return;
        }

        if (nextStatus.status === "failed") {
          setLoading(false);
          setError(nextStatus.error ?? nextStatus.message);
          return;
        }

        window.setTimeout(pollStatus, 700);
      } catch (pollError) {
        if (cancelled) {
          return;
        }
        const message = pollError instanceof Error ? pollError.message : "Unable to fetch progress";
        setLoading(false);
        setError(message);
      }
    };

    void pollStatus();

    return () => {
      cancelled = true;
    };
  }, [taskId]);

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    if (inputMode === "merged" && !dataFile) {
      alert("Upload the merged cohort spreadsheet first.");
      return;
    }
    if (inputMode === "merged" && (!positiveKeyword.trim() || !negativeKeyword.trim())) {
      alert("Enter both positive and negative column keywords.");
      return;
    }
    if (inputMode === "separate" && (!positiveFile || !negativeFile)) {
      alert("Upload both positive and negative cohort files.");
      return;
    }
    const parsedK = Number.parseInt(kValue, 10);
    if (!Number.isInteger(parsedK) || parsedK < 4 || parsedK > 10) {
      setError("k-mer size must be an integer between 4 and 10.");
      return;
    }

    const wildcardTokens = wildcards
      .split(",")
      .map((token) => token.trim())
      .filter(Boolean);
    if (wildcardTokens.some((token) => !/^\d+$/.test(token))) {
      setError("Wildcard positions must be comma-separated integers.");
      return;
    }
    const wildcardValues = wildcardTokens.map(Number);
    if (new Set(wildcardValues).size !== wildcardValues.length) {
      setError("Wildcard positions must be unique.");
      return;
    }
    if (wildcardValues.some((position) => position < 0 || position >= parsedK)) {
      setError(`Wildcard positions must be between 0 and ${parsedK - 1}.`);
      return;
    }
    const combinationCount = 20 ** (parsedK - wildcardValues.length);
    if (combinationCount > 1_000_000) {
      setError(
        `This configuration requires ${combinationCount.toLocaleString()} prebuilt k-mers. ` +
          "Reduce k or add wildcard positions to stay at or below 1,000,000.",
      );
      return;
    }

    const formData = new FormData();
    formData.append("input_mode", inputMode);
    if (inputMode === "merged" && dataFile) {
      formData.append("data_file", dataFile);
      formData.append("positive_keyword", positiveKeyword.trim());
      formData.append("negative_keyword", negativeKeyword.trim());
    }
    if (inputMode === "separate" && positiveFile && negativeFile) {
      formData.append("positive_file", positiveFile);
      formData.append("negative_file", negativeFile);
    }
    formData.append("k", String(parsedK));
    formData.append("archive_name", archiveName);
    formData.append("wildcard_positions", wildcards);
    formData.append("normalize", String(normalize));

    setError(null);
    setTaskId(null);
    setResult(null);
    setTaskStatus({
      task_id: "pending",
      status: "queued",
      progress: 0,
      message: "Uploading analysis input",
      result: null,
      error: null,
    });
    setLoading(true);

    try {
      const response = await apiClient.post<KmerTaskCreatedResponse>("/analyze-kmers", formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });
      setTaskId(response.data.task_id);
    } catch (submitError) {
      const message = submitError instanceof Error ? submitError.message : "Unable to start analysis";
      setLoading(false);
      setError(message);
    }
  };

  const runs = useMemo(() => result?.runs ?? [], [result]);
  const proteomeMappingState =
    result && runs[0]
      ? {
          resultId: result.result_id,
          positiveFilename: runs[0].ad_filename,
          negativeFilename: runs[0].nc_filename,
          inputMode,
          positiveKeyword,
          negativeKeyword,
        }
      : undefined;

  const renderUploadCard = (
    inputId: string,
    label: string,
    helperText: string,
    selectedFileName: string | null,
    onFileChange: (file: File | null) => void,
  ) => (
    <div className="upload-card kmer-upload-card">
      <div className="upload-card__top">
        <div>
          <span className="upload-card__label">{label}</span>
          <p className="upload-card__helper">{helperText}</p>
        </div>
        <label className="upload-trigger" htmlFor={inputId}>
          Choose file
        </label>
      </div>
      <input
        id={inputId}
        className="sr-only-file-input"
        type="file"
        accept=".xlsx,.csv"
        onChange={(event) => {
          const files = event.target.files;
          const selectedFile = files && files.length > 0 ? files[0] : null;
          onFileChange(selectedFile);
        }}
      />
      <div className={`upload-file-display ${selectedFileName ? "upload-file-display--selected" : ""}`}>
        {selectedFileName ?? "No file selected"}
      </div>
    </div>
  );

  return (
    <div className="page-grid">
      <SectionCard
        title="K-mer enrichment analysis"
        description="Build the complete k-mer universe, apply product-based chi-square filtering to each sample, and run Mann–Whitney U tests."
      >
        <form className="form-grid" onSubmit={handleSubmit}>
          <fieldset className="form-field input-mode-group">
            <legend>Input type</legend>
            <div className="input-mode-options">
              <label className={`input-mode-option ${inputMode === "merged" ? "input-mode-option--active" : ""}`}>
                <input
                  type="radio"
                  name="kmer-input-mode"
                  value="merged"
                  checked={inputMode === "merged"}
                  onChange={() => setInputMode("merged")}
                />
                <span>Merged cohort file</span>
                <small>Split columns using positive and negative keywords.</small>
              </label>
              <label className={`input-mode-option ${inputMode === "separate" ? "input-mode-option--active" : ""}`}>
                <input
                  type="radio"
                  name="kmer-input-mode"
                  value="separate"
                  checked={inputMode === "separate"}
                  onChange={() => setInputMode("separate")}
                />
                <span>Separate cohort files</span>
                <small>Upload positive and negative cohort files directly.</small>
              </label>
            </div>
          </fieldset>

          {inputMode === "merged" ? (
            <>
              {renderUploadCard(
                "kmer-merged-file",
                "Merged cohort file",
                "Original file should contain alternating cohort columns with sequence/count pairs.",
                dataFileName,
                (selectedFile) => {
                  setDataFile(selectedFile);
                  setDataFileName(selectedFile?.name ?? null);
                },
              )}

              <div className="settings-grid">
                <label className="form-field">
                  <span>Positive column keyword</span>
                  <input
                    type="text"
                    value={positiveKeyword}
                    onChange={(event) => setPositiveKeyword(event.target.value)}
                    placeholder="AD"
                  />
                  <small>Columns starting with this keyword are treated as the positive cohort.</small>
                </label>

                <label className="form-field">
                  <span>Negative column keyword</span>
                  <input
                    type="text"
                    value={negativeKeyword}
                    onChange={(event) => setNegativeKeyword(event.target.value)}
                    placeholder="NC"
                  />
                  <small>Columns starting with this keyword are treated as the negative cohort.</small>
                </label>
              </div>
            </>
          ) : (
            <div className="settings-grid">
              {renderUploadCard(
                "kmer-positive-file",
                "Positive cohort file",
                "Upload a file containing positive cohort sequence/count columns.",
                positiveFileName,
                (selectedFile) => {
                  setPositiveFile(selectedFile);
                  setPositiveFileName(selectedFile?.name ?? null);
                },
              )}

              {renderUploadCard(
                "kmer-negative-file",
                "Negative cohort file",
                "Upload a file containing negative cohort sequence/count columns.",
                negativeFileName,
                (selectedFile) => {
                  setNegativeFile(selectedFile);
                  setNegativeFileName(selectedFile?.name ?? null);
                },
              )}
            </div>
          )}

          <div className="form-field form-field--inline">
            <label>
              <span>k-mer size</span>
              <input
                type="text"
                inputMode="numeric"
                pattern="[0-9]*"
                value={kValue}
                onChange={(event) => setKValue(event.target.value)}
                placeholder="4"
              />
              <small>Enter one integer between 4 and 10.</small>
            </label>
          </div>

          <label className="form-field">
            <span>Download bundle name</span>
            <input
              type="text"
              value={archiveName}
              onChange={(event) => setArchiveName(event.target.value)}
              placeholder="e.g. patient_cohort_k4_results"
            />
            <small>Optional. Used as the ZIP filename when users download the result bundle.</small>
          </label>

          <label className="form-field">
            <span>Wildcard positions</span>
            <input
              type="text"
              placeholder="e.g. 1,3"
              value={wildcards}
              onChange={(event) => setWildcards(event.target.value)}
            />
            <small>Leave empty for strict k-mers. Use unique, comma-separated, zero-indexed positions.</small>
          </label>

          <label className="form-field form-checkbox">
            <input
              type="checkbox"
              checked={normalize}
              onChange={(event) => setNormalize(event.target.checked)}
            />
            <span>Normalize the matrix and statistical values to proportions</span>
            <small>Uses each sample's pre-filter tiled total. Turn this off to keep raw filtered counts.</small>
          </label>

          <div className="form-actions">
            <button type="submit" className="primary-button" disabled={loading}>
              {loading ? "Running analysis..." : "Start analysis"}
            </button>
            {result && <DownloadButton resultId={result.result_id} label="Download CSV bundle" />}
          </div>
        </form>

        {taskStatus && loading && (
          <div className="progress-panel" aria-live="polite">
            <div className="progress-panel__header">
              <strong>{taskStatus.message}</strong>
              <span>{taskStatus.progress}%</span>
            </div>
            <div className="progress-bar" role="progressbar" aria-valuenow={taskStatus.progress} aria-valuemin={0} aria-valuemax={100}>
              <div className="progress-bar__fill" style={{ width: `${taskStatus.progress}%` }} />
            </div>
          </div>
        )}

        {error && <StatusBanner tone="error" title="K-mer analysis failed" message={error} />}
        {result && !error && (
          <StatusBanner
            tone="success"
            title="K-mer analysis ready"
            message="Download CSV exports for the completed k-mer run."
          />
        )}
      </SectionCard>

      {runs.length > 0 && result && (
        <SectionCard
          title="Mann–Whitney summaries"
          description="Highlights the number of kmers enriched in each cohort for the selected k value."
          actions={
            <div className="section-card__actions-group">
              <Link to="/module3" state={proteomeMappingState} className="secondary-button">
                Open Proteome Mapping
              </Link>
              <DownloadButton resultId={result.result_id} label="Download combined outputs" />
            </div>
          }
        >
          <div className="table-wrapper">
            <table>
              <thead>
                <tr>
                  <th>k-mer size</th>
                  <th>Total kmers</th>
                  <th>Positive-elevated</th>
                  <th>Negative-elevated</th>
                  <th>Files</th>
                </tr>
              </thead>
              <tbody>
                {runs.map((run: KmerResultSummary) => (
                  <tr key={run.k}>
                    <td>{run.k}</td>
                    <td>{run.total_kmers.toLocaleString()}</td>
                    <td>{run.ad_elevated.toLocaleString()}</td>
                    <td>{run.nc_elevated.toLocaleString()}</td>
                    <td>
                      <ul className="file-list">
                        <li>{run.result_filename}</li>
                        <li>{run.ad_filename}</li>
                        <li>{run.nc_filename}</li>
                        <li>{run.matrix_filename}</li>
                      </ul>
                    </td>
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
