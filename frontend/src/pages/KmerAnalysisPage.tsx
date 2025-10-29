import { FormEvent, useMemo, useState } from "react";

import { apiClient } from "../api/client";
import { KmerResponse, KmerResultSummary } from "../api/types";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";

export function KmerAnalysisPage() {
  const [dataFile, setDataFile] = useState<File | null>(null);
  const [kMin, setKMin] = useState(4);
  const [kMax, setKMax] = useState(7);
  const [wildcards, setWildcards] = useState("");
  const [normalize, setNormalize] = useState(true);
  const [result, setResult] = useState<KmerResponse | null>(null);

  const { execute, loading, error } = useAsyncTask(async (formData: FormData) => {
    const response = await apiClient.post<KmerResponse>("/analyze-kmers", formData, {
      headers: { "Content-Type": "multipart/form-data" },
    });
    return response.data;
  });

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    if (!dataFile) {
      alert("Upload the merged AD/NC spreadsheet first.");
      return;
    }
    const formData = new FormData();
    formData.append("data_file", dataFile);
    formData.append("k_min", String(kMin));
    formData.append("k_max", String(kMax));
    formData.append("wildcard_positions", wildcards);
    formData.append("normalize", String(normalize));

    const response = await execute(formData);
    if (response) {
      setResult(response);
    }
  };

  const runs = useMemo(() => result?.runs ?? [], [result]);

  return (
    <div className="page-grid">
      <SectionCard
        title="K-mer enrichment analysis"
        description="Split AD vs. NC cohorts, tile peptides, and run Mann–Whitney U tests across k-mer windows."
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
              <span>k-min</span>
              <input
                type="number"
                min={4}
                max={10}
                value={kMin}
                onChange={(event) => setKMin(Number(event.target.value))}
              />
            </label>
            <label>
              <span>k-max</span>
              <input
                type="number"
                min={kMin}
                max={10}
                value={kMax}
                onChange={(event) => setKMax(Number(event.target.value))}
              />
            </label>
          </div>

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
              {loading ? "Running tests..." : "Start analysis"}
            </button>
            {result && <DownloadButton resultId={result.result_id} label="Download CSV bundle" />}
          </div>
        </form>

        {error && <StatusBanner tone="error" title="K-mer analysis failed" message={error} />}
        {result && !error && (
          <StatusBanner
            tone="success"
            title="K-mer analysis ready"
            message="Download CSV exports for each window size."
          />
        )}
      </SectionCard>

      {runs.length > 0 && result && (
        <SectionCard
          title="Mann–Whitney summaries"
          description="Highlights the number of kmers enriched in each cohort for every k value."
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
