import { Navigate, Route, Routes, useLocation } from "react-router-dom";

import { FastqProcessingPage } from "./pages/FastqProcessingPage";
import { KmerAnalysisPage } from "./pages/KmerAnalysisPage";
import { HeaderNav } from "./components/HeaderNav";

export default function App() {
  const location = useLocation();
  return (
    <div className="app-shell">
      <HeaderNav currentPath={location.pathname} />
      <main className="app-main">
        <Routes>
          <Route path="/fastq" element={<FastqProcessingPage />} />
          <Route path="/kmer" element={<KmerAnalysisPage />} />
          <Route path="*" element={<Navigate to="/fastq" replace />} />
        </Routes>
      </main>
      <footer className="app-footer">
        <p>
          Built for antibody sequencing workflows. Upload your datasets, monitor progress,
          and download reproducible reports.
        </p>
      </footer>
    </div>
  );
}
