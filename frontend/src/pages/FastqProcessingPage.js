import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { useMemo, useState } from "react";
import { apiClient } from "../api/client";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";
export function FastqProcessingPage() {
    const [files, setFiles] = useState([]);
    const [backgroundFiles, setBackgroundFiles] = useState([]);
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
        backgroundFiles.forEach((file) => formData.append("background_files", file));
        formData.append("output_name", outputName);
        const response = await execute(formData);
        if (response) {
            setResult(response);
        }
    };
    const summaryRows = useMemo(() => result?.summary ?? [], [result]);
    const formatSelectedFiles = (selectedFiles) => {
        if (selectedFiles.length === 0) {
            return "No files selected";
        }
        if (selectedFiles.length === 1) {
            return selectedFiles[0].name;
        }
        return `${selectedFiles.length} files selected: ${selectedFiles.map((file) => file.name).join(", ")}`;
    };
    const renderUploadCard = (inputId, label, helperText, accept, selectedFiles, onFilesChange) => (_jsxs("div", { className: "upload-card fastq-upload-card", children: [_jsxs("div", { className: "upload-card__top", children: [_jsxs("div", { children: [_jsx("span", { className: "upload-card__label", children: label }), _jsx("p", { className: "upload-card__helper", children: helperText })] }), _jsx("label", { className: "upload-trigger", htmlFor: inputId, children: "Choose files" })] }), _jsx("input", { id: inputId, className: "sr-only-file-input", type: "file", accept: accept, multiple: true, onChange: (event) => {
                    const targetFiles = event.target.files;
                    onFilesChange(targetFiles ? Array.from(targetFiles) : []);
                } }), _jsx("div", { className: `upload-file-display ${selectedFiles.length > 0 ? "upload-file-display--selected" : ""}`, children: formatSelectedFiles(selectedFiles) })] }));
    return (_jsxs("div", { className: "page-grid", children: [_jsxs(SectionCard, { title: "Process FASTQ libraries", description: "Upload sequencing runs to generate ranked peptide matrices with optional background subtraction.", children: [_jsxs("form", { className: "form-grid", onSubmit: handleSubmit, children: [renderUploadCard("fastq-files", "FASTQ files", "Supports plain or gzipped FASTQ files. Multiple uploads allowed.", ".fastq,.fq,.fastq.gz,.fq.gz", files, setFiles), renderUploadCard("fastq-background-files", "Background files (optional)", "Upload one or more FASTQ, FASTQ.GZ, or TXT background files. They will be combined before filtering.", ".fastq,.fq,.fastq.gz,.fq.gz,.txt", backgroundFiles, setBackgroundFiles), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Output name" }), _jsx("input", { type: "text", value: outputName, onChange: (event) => setOutputName(event.target.value), placeholder: "sequence_matrix.xlsx" })] }), _jsxs("div", { className: "form-actions", children: [_jsx("button", { type: "submit", className: "primary-button", disabled: loading, children: loading ? "Processing..." : "Run analysis" }), result && _jsx(DownloadButton, { resultId: result.result_id })] })] }), error && _jsx(StatusBanner, { tone: "error", title: "Failed to process FASTQ files", message: error }), result && !error && (_jsx(StatusBanner, { tone: "success", title: "FASTQ analysis complete", message: `Primary matrix saved as ${result.excel_filename}` }))] }), summaryRows.length > 0 && result && (_jsx(SectionCard, { title: "Sample overview", description: "Per-sample counts after background filtering. Columns match the Excel output.", actions: _jsx(DownloadButton, { resultId: result.result_id, label: "Download full bundle" }), children: _jsx("div", { className: "table-wrapper", children: _jsxs("table", { children: [_jsx("thead", { children: _jsxs("tr", { children: [_jsx("th", { children: "Sample" }), _jsx("th", { children: "Total reads" }), _jsx("th", { children: "Unique peptides" }), _jsx("th", { children: "Filtered (background)" }), _jsx("th", { children: "Output columns" })] }) }), _jsx("tbody", { children: summaryRows.map((row) => (_jsxs("tr", { children: [_jsx("td", { children: row.sample_name }), _jsx("td", { children: row.total_sequences.toLocaleString() }), _jsx("td", { children: row.unique_sequences.toLocaleString() }), _jsx("td", { children: row.filtered_sequences.toLocaleString() }), _jsx("td", { children: row.output_columns.join(", ") })] }, row.sample_name))) })] }) }) }))] }));
}
