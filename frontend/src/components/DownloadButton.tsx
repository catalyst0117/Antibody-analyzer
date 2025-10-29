import { useState } from "react";

import { apiClient } from "../api/client";

type DownloadButtonProps = {
  resultId: string;
  label?: string;
};

export function DownloadButton({ resultId, label = "Download results" }: DownloadButtonProps) {
  const [downloading, setDownloading] = useState(false);

  const handleDownload = async () => {
    setDownloading(true);
    try {
      const response = await apiClient.get(`/results/${resultId}/download`, {
        responseType: "blob",
      });
      const blob = new Blob([response.data], { type: "application/zip" });
      const url = window.URL.createObjectURL(blob);
      const anchor = document.createElement("a");
      anchor.href = url;
      anchor.download = `antibody-results-${resultId}.zip`;
      anchor.click();
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error("Download failed", error);
    } finally {
      setDownloading(false);
    }
  };

  return (
    <button className="primary-button" onClick={handleDownload} disabled={downloading}>
      {downloading ? "Preparing archive..." : label}
    </button>
  );
}
