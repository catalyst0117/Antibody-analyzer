import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { useMemo, useState } from "react";
import { apiClient } from "../api/client";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";
export function KmerAnalysisPage() {
    const [dataFile, setDataFile] = useState(null);
    const [kMin, setKMin] = useState(4);
    const [kMax, setKMax] = useState(7);
    const [wildcards, setWildcards] = useState("");
    const [normalize, setNormalize] = useState(true);
    const [result, setResult] = useState(null);
    const [useRange, setUseRange] = useState(false);
    const { execute, loading, error } = useAsyncTask(async (formData) => {
        const response = await apiClient.post("/analyze-kmers", formData, {
            headers: { "Content-Type": "multipart/form-data" },
        });
        return response.data;
    });
    const handleSubmit = async (event) => {
        event.preventDefault();
        if (!dataFile) {
            alert("Upload the merged AD/NC spreadsheet first.");
            return;
        }
        const formData = new FormData();
        formData.append("data_file", dataFile);
        formData.append("k_min", String(kMin));
        if (useRange) {
            formData.append("k_max", String(kMax));
        }
        formData.append("wildcard_positions", wildcards);
        formData.append("normalize", String(normalize));
        const response = await execute(formData);
        if (response) {
            setResult(response);
        }
    };
    const runs = useMemo(() => result?.runs ?? [], [result]);
    return (_jsxs("div", { className: "page-grid", children: [_jsxs(SectionCard, { title: "K-mer enrichment analysis", description: "Split AD vs. NC cohorts, tile peptides, and run Mann\u2013Whitney U tests across k-mer windows.", children: [_jsxs("form", { className: "form-grid", onSubmit: handleSubmit, children: [_jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Merged cohort file" }), _jsx("input", { type: "file", accept: ".xlsx,.csv", onChange: (event) => {
                                            const files = event.target.files;
                                            setDataFile(files && files.length > 0 ? files[0] : null);
                                        } }), _jsx("small", { children: "Original file should contain alternating AD_/NC_ columns with sequence/count pairs." })] }), _jsx("div", { className: "form-field form-field--inline", children: _jsxs("label", { children: [_jsx("span", { children: "k-min" }), _jsx("input", { type: "number", min: 4, max: 10, value: kMin, onChange: (event) => setKMin(Number(event.target.value)) })] }) }), _jsxs("label", { className: "form-field form-checkbox", children: [_jsx("input", { type: "checkbox", checked: useRange, onChange: (event) => setUseRange(event.target.checked) }), _jsx("span", { children: "Run a range of k values" })] }), useRange && (_jsx("div", { className: "form-field form-field--inline", children: _jsxs("label", { children: [_jsx("span", { children: "k-max" }), _jsx("input", { type: "number", min: kMin, max: 10, value: kMax, onChange: (event) => setKMax(Number(event.target.value)) })] }) })), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Wildcard positions" }), _jsx("input", { type: "text", placeholder: "e.g. 1,3", value: wildcards, onChange: (event) => setWildcards(event.target.value) }), _jsx("small", { children: "Leave empty for strict kmers. Use comma separated zero-index positions." })] }), _jsxs("label", { className: "form-field form-checkbox", children: [_jsx("input", { type: "checkbox", checked: normalize, onChange: (event) => setNormalize(event.target.checked) }), _jsx("span", { children: "Normalize counts to proportions before statistical testing" })] }), _jsxs("div", { className: "form-actions", children: [_jsx("button", { type: "submit", className: "primary-button", disabled: loading, children: loading ? "Running tests..." : "Start analysis" }), result && _jsx(DownloadButton, { resultId: result.result_id, label: "Download CSV bundle" })] })] }), error && _jsx(StatusBanner, { tone: "error", title: "K-mer analysis failed", message: error }), result && !error && (_jsx(StatusBanner, { tone: "success", title: "K-mer analysis ready", message: "Download CSV exports for each window size." }))] }), runs.length > 0 && result && (_jsx(SectionCard, { title: "Mann\u2013Whitney summaries", description: "Highlights the number of kmers enriched in each cohort for every k value.", actions: _jsx(DownloadButton, { resultId: result.result_id, label: "Download combined outputs" }), children: _jsx("div", { className: "table-wrapper", children: _jsxs("table", { children: [_jsx("thead", { children: _jsxs("tr", { children: [_jsx("th", { children: "k-mer size" }), _jsx("th", { children: "Total kmers" }), _jsx("th", { children: "AD-elevated" }), _jsx("th", { children: "NC-elevated" }), _jsx("th", { children: "Files" })] }) }), _jsx("tbody", { children: runs.map((run) => (_jsxs("tr", { children: [_jsx("td", { children: run.k }), _jsx("td", { children: run.total_kmers.toLocaleString() }), _jsx("td", { children: run.ad_elevated.toLocaleString() }), _jsx("td", { children: run.nc_elevated.toLocaleString() }), _jsx("td", { children: _jsxs("ul", { className: "file-list", children: [_jsx("li", { children: run.result_filename }), _jsx("li", { children: run.ad_filename }), _jsx("li", { children: run.nc_filename }), _jsx("li", { children: run.matrix_filename })] }) })] }, run.k))) })] }) }) }))] }));
}
