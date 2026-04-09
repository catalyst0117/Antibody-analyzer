import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { useMemo, useState } from "react";
import { apiClient } from "../api/client";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";
export function FastqProcessingPage() {
    const [files, setFiles] = useState([]);
    const [background, setBackground] = useState(null);
    const [outputName, setOutputName] = useState("sequence_matrix.xlsx");
    const [result, setResult] = useState(null);
    const { execute, loading, error } = useAsyncTask(async (formData) => {
        const response = await apiClient.post("/process-fastq", formData, {
            headers: { "Content-Type": "multipart/form-data" },
        });
        return response.data;
    });
    const handleSubmit = async (event) => {
        event.preventDefault();
        if (files.length === 0) {
            alert("Please add at least one FASTQ or FASTQ.GZ file.");
            return;
        }
        const formData = new FormData();
        files.forEach((file) => formData.append("files", file));
        if (background) {
            formData.append("background_file", background);
        }
        formData.append("output_name", outputName);
        const response = await execute(formData);
        if (response) {
            setResult(response);
        }
    };
    const summaryRows = useMemo(() => result?.summary ?? [], [result]);
    return (_jsxs("div", { className: "page-grid", children: [_jsxs(SectionCard, { title: "Process FASTQ libraries", description: "Upload sequencing runs to generate ranked peptide matrices with optional background subtraction.", children: [_jsxs("form", { className: "form-grid", onSubmit: handleSubmit, children: [_jsxs("label", { className: "form-field", children: [_jsx("span", { children: "FASTQ files" }), _jsx("input", { className: "file-input", type: "file", accept: ".fastq,.fq,.fastq.gz,.fq.gz", multiple: true, onChange: (event) => {
                                            const targetFiles = event.target.files;
                                            setFiles(targetFiles ? Array.from(targetFiles) : []);
                                        } }), _jsx("small", { children: "Supports plain or gzipped FASTQ files. Multiple uploads allowed." })] }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Background FASTQ (optional)" }), _jsx("input", { className: "file-input", type: "file", accept: ".fastq,.fq,.fastq.gz,.fq.gz", onChange: (event) => {
                                            const fileList = event.target.files;
                                            setBackground(fileList && fileList.length > 0 ? fileList[0] : null);
                                        } }), _jsx("small", { children: "Sequences detected in the background file will be filtered out." })] }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Output name" }), _jsx("input", { type: "text", value: outputName, onChange: (event) => setOutputName(event.target.value), placeholder: "sequence_matrix.xlsx" })] }), _jsxs("div", { className: "form-actions", children: [_jsx("button", { type: "submit", className: "primary-button", disabled: loading, children: loading ? "Processing..." : "Run analysis" }), result && _jsx(DownloadButton, { resultId: result.result_id })] })] }), error && _jsx(StatusBanner, { tone: "error", title: "Failed to process FASTQ files", message: error }), result && !error && (_jsx(StatusBanner, { tone: "success", title: "FASTQ analysis complete", message: `Primary matrix saved as ${result.excel_filename}` }))] }), summaryRows.length > 0 && result && (_jsx(SectionCard, { title: "Sample overview", description: "Per-sample counts after background filtering. Columns match the Excel output.", actions: _jsx(DownloadButton, { resultId: result.result_id, label: "Download full bundle" }), children: _jsx("div", { className: "table-wrapper", children: _jsxs("table", { children: [_jsx("thead", { children: _jsxs("tr", { children: [_jsx("th", { children: "Sample" }), _jsx("th", { children: "Total reads" }), _jsx("th", { children: "Unique peptides" }), _jsx("th", { children: "Filtered (background)" }), _jsx("th", { children: "Output columns" })] }) }), _jsx("tbody", { children: summaryRows.map((row) => (_jsxs("tr", { children: [_jsx("td", { children: row.sample_name }), _jsx("td", { children: row.total_sequences.toLocaleString() }), _jsx("td", { children: row.unique_sequences.toLocaleString() }), _jsx("td", { children: row.filtered_sequences.toLocaleString() }), _jsx("td", { children: row.output_columns.join(", ") })] }, row.sample_name))) })] }) }) }))] }));
}
