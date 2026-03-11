import { FormEvent, useEffect, useMemo, useState } from "react";

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

export function KmerAnalysisPage() {
  const [dataFile, setDataFile] = useState<File | null>(null);
  const [kValue, setKValue] = useState("4");
  const [archiveName, setArchiveName] = useState("");
  const [wildcards, setWildcards] = useState("");
  const [normalize, setNormalize] = useState(true);
  const [result, setResult] = useState<KmerResponse | null>(null);
  const [taskId, setTaskId] = useState<string | null>(null);
  const [taskStatus, setTaskStatus] = useState<KmerTaskStatusResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

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
    if (!dataFile) {
      alert("Upload the merged AD/NC spreadsheet first.");
      return;
    }
    const parsedK = Number.parseInt(kValue, 10);
    if (!Number.isInteger(parsedK) || parsedK < 4 || parsedK > 10) {
      setError("k-mer size must be an integer between 4 and 10.");
      return;
    }

    const formData = new FormData();
    formData.append("data_file", dataFile);
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

  return (
    <div className="page-grid">
      <SectionCard
        title="K-mer enrichment analysis"
        description="Split AD vs. NC cohorts, tile peptides, and run Mann–Whitney U tests for one k-mer size at a time."
      >
        <form className="form-grid" onSubmit={handleSubmit}>
          <label className="form-field">
            <span>Merged cohort file</span>
            <input
              type="file"
              accept=".xlsx,.csv"
              onChange={(event) => {
                const files = event.target.files;
                setDataFile(files && files.length > 0 ? files[0] : null);
              }}
            />
            <small>Original file should contain alternating AD_/NC_ columns with sequence/count pairs.</small>
          </label>

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
            <small>Leave empty for strict kmers. Use comma separated zero-index positions.</small>
          </label>

          <label className="form-field form-checkbox">
            <input
              type="checkbox"
              checked={normalize}
              onChange={(event) => setNormalize(event.target.checked)}
            />
            <span>Normalize counts to proportions before statistical testing</span>
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
          actions={<DownloadButton resultId={result.result_id} label="Download combined outputs" />}
        >
          <div className="table-wrapper">
            <table>
              <thead>
                <tr>
                  <th>k-mer size</th>
                  <th>Total kmers</th>
                  <th>AD-elevated</th>
                  <th>NC-elevated</th>
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
