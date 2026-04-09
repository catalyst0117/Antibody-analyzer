import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { useEffect, useMemo, useState } from "react";
import { apiClient } from "../api/client";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
import { useAsyncTask } from "../hooks/useAsyncTask";
const PROTEOME_MAPPING_SESSION_KEY = "proteome-mapping-page-state";
function loadProteomeMappingSessionState() {
    if (typeof window === "undefined") {
        return null;
    }
    const raw = window.sessionStorage.getItem(PROTEOME_MAPPING_SESSION_KEY);
    if (!raw) {
        return null;
    }
    try {
        return JSON.parse(raw);
    }
    catch {
        return null;
    }
}
export function Module3MappingPage() {
    const savedState = loadProteomeMappingSessionState();
    const [positiveFile, setPositiveFile] = useState(null);
    const [negativeFile, setNegativeFile] = useState(null);
    const [proteomeFile, setProteomeFile] = useState(null);
    const [positiveFileName, setPositiveFileName] = useState(savedState?.positiveFileName ?? null);
    const [negativeFileName, setNegativeFileName] = useState(savedState?.negativeFileName ?? null);
    const [proteomeFileName, setProteomeFileName] = useState(savedState?.proteomeFileName ?? null);
    const [outputFolderName, setOutputFolderName] = useState(savedState?.outputFolderName ?? "");
    const [topNValue, setTopNValue] = useState(savedState?.topNValue ?? "");
    const [qCutoffValue, setQCutoffValue] = useState(savedState?.qCutoffValue ?? "0.01");
    const [wildcards, setWildcards] = useState(savedState?.wildcards ?? false);
    const [result, setResult] = useState(savedState?.result ?? null);
    const { execute, loading, error } = useAsyncTask(async (formData) => {
        const response = await apiClient.post("/module3-map", formData, {
            headers: { "Content-Type": "multipart/form-data" },
        });
        return response.data;
    });
    useEffect(() => {
        if (typeof window === "undefined") {
            return;
        }
        const nextState = {
            outputFolderName,
            positiveFileName,
            negativeFileName,
            proteomeFileName,
            topNValue,
            qCutoffValue,
            wildcards,
            result,
        };
        window.sessionStorage.setItem(PROTEOME_MAPPING_SESSION_KEY, JSON.stringify(nextState));
    }, [negativeFileName, outputFolderName, positiveFileName, proteomeFileName, qCutoffValue, result, topNValue, wildcards]);
    const handleSubmit = async (event) => {
        event.preventDefault();
        if (!positiveFile || !negativeFile || !proteomeFile) {
            alert("Please upload the positive file, negative file, and proteome FASTA.");
            return;
        }
        const parsedTopN = topNValue.trim() === "" ? null : Number.parseInt(topNValue, 10);
        if (parsedTopN !== null && (!Number.isInteger(parsedTopN) || parsedTopN <= 0)) {
            alert("top_n must be a positive integer.");
            return;
        }
        const parsedQCutoff = Number.parseFloat(qCutoffValue);
        if (Number.isNaN(parsedQCutoff) || parsedQCutoff < 0 || parsedQCutoff > 1) {
            alert("q_cutoff must be a number between 0 and 1.");
            return;
        }
        const formData = new FormData();
        formData.append("positive_file", positiveFile);
        formData.append("negative_file", negativeFile);
        formData.append("proteome_fasta", proteomeFile);
        formData.append("output_folder_name", outputFolderName);
        if (parsedTopN !== null) {
            formData.append("top_n", String(parsedTopN));
        }
        formData.append("wildcards", String(wildcards));
        formData.append("q_cutoff", String(parsedQCutoff));
        const response = await execute(formData);
        if (response) {
            setResult(response);
        }
    };
    const outputFiles = useMemo(() => {
        if (!result) {
            return [];
        }
        return [
            result.positive_mapping_filename,
            result.negative_mapping_filename,
            result.positive_clean_filename,
            result.negative_clean_filename,
            result.positive_manhattan_filename,
            result.negative_manhattan_filename,
            "run_summary.txt",
        ];
    }, [result]);
    return (_jsxs("div", { className: "page-grid proteome-page", children: [_jsxs(SectionCard, { title: "Proteome Mapping", description: "Map significant positive and negative k-mers onto a local proteome FASTA and download the bundled mapping outputs.", children: [_jsxs("form", { className: "proteome-form", onSubmit: handleSubmit, children: [_jsxs("section", { className: "proteome-section", children: [_jsxs("div", { className: "proteome-section__header", children: [_jsx("h3", { children: "Input Files" }), _jsx("p", { children: "Upload the significance tables and the reference FASTA needed for one mapping run." })] }), _jsxs("div", { className: "upload-card-grid", children: [_jsxs("div", { className: "upload-card", children: [_jsxs("div", { className: "upload-card__top", children: [_jsxs("div", { children: [_jsx("span", { className: "upload-card__label", children: "Positive significance file" }), _jsx("p", { className: "upload-card__helper", children: "Upload the positive significance output from Module 2 or another compatible table." })] }), _jsx("label", { className: "upload-trigger", htmlFor: "proteome-positive-file", children: "Choose file" })] }), _jsx("input", { id: "proteome-positive-file", className: "sr-only-file-input", type: "file", accept: ".xlsx,.xls,.csv,.txt", onChange: (event) => {
                                                            const files = event.target.files;
                                                            const selectedFile = files && files.length > 0 ? files[0] : null;
                                                            setPositiveFile(selectedFile);
                                                            setPositiveFileName(selectedFile?.name ?? null);
                                                        } }), _jsx("div", { className: `upload-file-display ${positiveFileName ? "upload-file-display--selected" : ""}`, children: positiveFileName ?? "No file selected" })] }), _jsxs("div", { className: "upload-card", children: [_jsxs("div", { className: "upload-card__top", children: [_jsxs("div", { children: [_jsx("span", { className: "upload-card__label", children: "Negative significance file" }), _jsx("p", { className: "upload-card__helper", children: "Upload the matching negative significance file with the same k-mer length." })] }), _jsx("label", { className: "upload-trigger", htmlFor: "proteome-negative-file", children: "Choose file" })] }), _jsx("input", { id: "proteome-negative-file", className: "sr-only-file-input", type: "file", accept: ".xlsx,.xls,.csv,.txt", onChange: (event) => {
                                                            const files = event.target.files;
                                                            const selectedFile = files && files.length > 0 ? files[0] : null;
                                                            setNegativeFile(selectedFile);
                                                            setNegativeFileName(selectedFile?.name ?? null);
                                                        } }), _jsx("div", { className: `upload-file-display ${negativeFileName ? "upload-file-display--selected" : ""}`, children: negativeFileName ?? "No file selected" })] }), _jsxs("div", { className: "upload-card", children: [_jsxs("div", { className: "upload-card__top", children: [_jsxs("div", { children: [_jsx("span", { className: "upload-card__label", children: "Proteome FASTA file" }), _jsx("p", { className: "upload-card__helper", children: "Local FASTA only. Proteome Mapping will build the required k-mer index from this file." })] }), _jsx("label", { className: "upload-trigger", htmlFor: "proteome-fasta-file", children: "Choose file" })] }), _jsx("input", { id: "proteome-fasta-file", className: "sr-only-file-input", type: "file", accept: ".fasta,.fa,.faa,.txt", onChange: (event) => {
                                                            const files = event.target.files;
                                                            const selectedFile = files && files.length > 0 ? files[0] : null;
                                                            setProteomeFile(selectedFile);
                                                            setProteomeFileName(selectedFile?.name ?? null);
                                                        } }), _jsx("div", { className: `upload-file-display ${proteomeFileName ? "upload-file-display--selected" : ""}`, children: proteomeFileName ?? "No file selected" })] })] })] }), _jsxs("section", { className: "proteome-section", children: [_jsxs("div", { className: "proteome-section__header", children: [_jsx("h3", { children: "Run Settings" }), _jsx("p", { children: "Set the output naming and filtering options for this mapping job." })] }), _jsxs("label", { className: "form-field proteome-full-width-field", children: [_jsx("span", { children: "Output folder name" }), _jsx("input", { type: "text", value: outputFolderName, onChange: (event) => setOutputFolderName(event.target.value), placeholder: "Optional. Leave blank to auto-generate" }), _jsx("small", { children: "Optional. If left blank, a dataset name will be generated from the uploaded significance file and top_n." })] }), _jsxs("div", { className: "settings-grid", children: [_jsxs("label", { className: "form-field", children: [_jsx("span", { children: "top_n" }), _jsx("input", { type: "text", inputMode: "numeric", pattern: "[0-9]*", value: topNValue, onChange: (event) => setTopNValue(event.target.value), placeholder: "Leave empty for all" }), _jsx("small", { children: "Optional. Keep only the top N lowest-q rows from each file." })] }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "q_cutoff" }), _jsx("input", { type: "text", value: qCutoffValue, onChange: (event) => setQCutoffValue(event.target.value), placeholder: "0.01" }), _jsx("small", { children: "Rows above this q-value are excluded from the proteome mapping step." })] })] }), _jsxs("label", { className: "form-field form-checkbox settings-toggle", children: [_jsx("input", { type: "checkbox", checked: wildcards, onChange: (event) => setWildcards(event.target.checked) }), _jsx("span", { children: "Enable wildcard X matching" })] })] }), _jsxs("div", { className: "proteome-actions", children: [_jsx("button", { type: "submit", className: "primary-button", disabled: loading, children: loading ? "Running Proteome Mapping..." : "Start Proteome Mapping" }), result && _jsx(DownloadButton, { resultId: result.result_id, label: "Download Proteome Mapping bundle" })] })] }), error && _jsx(StatusBanner, { tone: "error", title: "Proteome Mapping failed", message: error }), result && !error && (_jsx(StatusBanner, { tone: "success", title: "Proteome Mapping complete", message: "The positive and negative proteome mapping outputs are ready to download." }))] }), result && outputFiles.length > 0 && (_jsx(SectionCard, { title: "Proteome Mapping Outputs", description: "These files are included in the downloadable ZIP bundle.", actions: _jsx(DownloadButton, { resultId: result.result_id, label: "Download combined outputs" }), children: _jsx("div", { className: "table-wrapper", children: _jsxs("table", { children: [_jsx("thead", { children: _jsxs("tr", { children: [_jsx("th", { children: "Setting" }), _jsx("th", { children: "Value" })] }) }), _jsxs("tbody", { children: [_jsxs("tr", { children: [_jsx("td", { children: "top_n" }), _jsx("td", { children: result.top_n ?? "All rows" })] }), _jsxs("tr", { children: [_jsx("td", { children: "Output folder name" }), _jsx("td", { children: result.output_folder_name })] }), _jsxs("tr", { children: [_jsx("td", { children: "Wildcard matching" }), _jsx("td", { children: result.wildcards ? "Enabled" : "Disabled" })] }), _jsxs("tr", { children: [_jsx("td", { children: "q_cutoff" }), _jsx("td", { children: result.q_cutoff })] }), _jsxs("tr", { children: [_jsx("td", { children: "Output files" }), _jsx("td", { children: outputFiles.join(", ") })] })] })] }) }) }))] }));
}
