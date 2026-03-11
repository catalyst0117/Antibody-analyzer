import { useState } from "react";

import { apiClient } from "../api/client";

type DownloadButtonProps = {
  resultId: string;
  label?: string;
};

export function DownloadButton({ resultId, label = "Download results" }: DownloadButtonProps) {
  const [downloading, setDownloading] = useState(false);

  const getDownloadName = (contentDisposition: string | undefined, fallbackName: string) => {
    if (!contentDisposition) {
      return fallbackName;
    }

    const utf8Match = contentDisposition.match(/filename\*=UTF-8''(?<name>[^;]+)/i);
    if (utf8Match?.groups?.name) {
      return decodeURIComponent(utf8Match.groups.name);
    }

    const plainMatch = contentDisposition.match(/filename="?(?<name>[^"]+)"?/i);
    return plainMatch?.groups?.name ?? fallbackName;
  };

  const handleDownload = async () => {
    setDownloading(true);
    try {
      const response = await apiClient.get(`/results/${resultId}/download`, {
        responseType: "blob",
      });
      const blob = new Blob([response.data], { type: "application/zip" });
      const url = window.URL.createObjectURL(blob);
      const anchor = document.createElement("a");
      const disposition = response.headers["content-disposition"] as string | undefined;
      const downloadName = getDownloadName(disposition, `antibody-results-${resultId}.zip`);
      anchor.href = url;
      anchor.download = downloadName;
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
