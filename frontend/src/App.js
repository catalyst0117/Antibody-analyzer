import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { Navigate, Route, Routes, useLocation } from "react-router-dom";
import { FastqProcessingPage } from "./pages/FastqProcessingPage";
import { KmerAnalysisPage } from "./pages/KmerAnalysisPage";
import { Module3MappingPage } from "./pages/Module3MappingPage";
import { HeaderNav } from "./components/HeaderNav";
export default function App() {
    const location = useLocation();
    return (_jsxs("div", { className: "app-shell", children: [_jsx(HeaderNav, { currentPath: location.pathname }), _jsx("main", { className: "app-main", children: _jsxs(Routes, { children: [_jsx(Route, { path: "/fastq", element: _jsx(FastqProcessingPage, {}) }), _jsx(Route, { path: "/kmer", element: _jsx(KmerAnalysisPage, {}) }), _jsx(Route, { path: "/module3", element: _jsx(Module3MappingPage, {}) }), _jsx(Route, { path: "*", element: _jsx(Navigate, { to: "/fastq", replace: true }) })] }) }), _jsx("footer", { className: "app-footer", children: _jsx("p", { children: "Built for antibody sequencing workflows. Upload your datasets, monitor progress, and download reproducible reports." }) })] }));
}
