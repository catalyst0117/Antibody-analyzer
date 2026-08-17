import { jsx as _jsx, jsxs as _jsxs, Fragment as _Fragment } from "react/jsx-runtime";
import { useEffect, useMemo, useState } from "react";
import { Link } from "react-router-dom";
import { apiClient } from "../api/client";
import { DownloadButton } from "../components/DownloadButton";
import { SectionCard } from "../components/SectionCard";
import { StatusBanner } from "../components/StatusBanner";
const KMER_SESSION_KEY = "kmer-analysis-page-state";
function loadKmerSessionState() {
    if (typeof window === "undefined") {
        return null;
    }
    const raw = window.sessionStorage.getItem(KMER_SESSION_KEY);
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
export function KmerAnalysisPage() {
    const savedState = loadKmerSessionState();
    const [inputMode, setInputMode] = useState(savedState?.inputMode ?? "merged");
    const [dataFile, setDataFile] = useState(null);
    const [positiveFile, setPositiveFile] = useState(null);
    const [negativeFile, setNegativeFile] = useState(null);
    const [dataFileName, setDataFileName] = useState(savedState?.dataFileName ?? null);
    const [positiveFileName, setPositiveFileName] = useState(savedState?.positiveFileName ?? null);
    const [negativeFileName, setNegativeFileName] = useState(savedState?.negativeFileName ?? null);
    const [positiveKeyword, setPositiveKeyword] = useState(savedState?.positiveKeyword ?? "AD");
    const [negativeKeyword, setNegativeKeyword] = useState(savedState?.negativeKeyword ?? "NC");
    const [kValue, setKValue] = useState(savedState?.kValue ?? "4");
    const [archiveName, setArchiveName] = useState(savedState?.archiveName ?? "");
    const [wildcards, setWildcards] = useState(savedState?.wildcards ?? "");
    const [normalize, setNormalize] = useState(savedState?.normalize ?? true);
    const [result, setResult] = useState(savedState?.result ?? null);
    const [taskId, setTaskId] = useState(savedState?.taskId ?? null);
    const [taskStatus, setTaskStatus] = useState(savedState?.taskStatus ?? null);
    const [loading, setLoading] = useState(savedState?.loading ?? false);
    const [error, setError] = useState(savedState?.error ?? null);
    useEffect(() => {
        if (typeof window === "undefined") {
            return;
        }
        const nextState = {
            inputMode,
            kValue,
            dataFileName,
            positiveFileName,
            negativeFileName,
            positiveKeyword,
            negativeKeyword,
            archiveName,
            wildcards,
            normalize,
            result,
            taskId,
            taskStatus,
            loading,
            error,
        };
        window.sessionStorage.setItem(KMER_SESSION_KEY, JSON.stringify(nextState));
    }, [
        archiveName,
        dataFileName,
        error,
        inputMode,
        kValue,
        loading,
        negativeFileName,
        negativeKeyword,
        normalize,
        positiveFileName,
        positiveKeyword,
        result,
        taskId,
        taskStatus,
        wildcards,
    ]);
    useEffect(() => {
        if (!taskId) {
            return;
        }
        let cancelled = false;
        const pollStatus = async () => {
            try {
                const response = await apiClient.get(`/analyze-kmers/${taskId}`);
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
            }
            catch (pollError) {
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
    const handleSubmit = async (event) => {
        event.preventDefault();
        if (inputMode === "merged" && !dataFile) {
            alert("Upload the merged cohort spreadsheet first.");
            return;
        }
        if (inputMode === "merged" && (!positiveKeyword.trim() || !negativeKeyword.trim())) {
            alert("Enter both positive and negative column keywords.");
            return;
        }
        if (inputMode === "separate" && (!positiveFile || !negativeFile)) {
            alert("Upload both positive and negative cohort files.");
            return;
        }
        const parsedK = Number.parseInt(kValue, 10);
        if (!Number.isInteger(parsedK) || parsedK < 4 || parsedK > 10) {
            setError("k-mer size must be an integer between 4 and 10.");
            return;
        }
        const wildcardTokens = wildcards
            .split(",")
            .map((token) => token.trim())
            .filter(Boolean);
        if (wildcardTokens.some((token) => !/^\d+$/.test(token))) {
            setError("Wildcard positions must be comma-separated integers.");
            return;
        }
        const wildcardValues = wildcardTokens.map(Number);
        if (new Set(wildcardValues).size !== wildcardValues.length) {
            setError("Wildcard positions must be unique.");
            return;
        }
        if (wildcardValues.some((position) => position < 0 || position >= parsedK)) {
            setError(`Wildcard positions must be between 0 and ${parsedK - 1}.`);
            return;
        }
        const combinationCount = 20 ** (parsedK - wildcardValues.length);
        if (combinationCount > 1000000) {
            setError(`This configuration requires ${combinationCount.toLocaleString()} prebuilt k-mers. ` +
                "Reduce k or add wildcard positions to stay at or below 1,000,000.");
            return;
        }
        const formData = new FormData();
        formData.append("input_mode", inputMode);
        if (inputMode === "merged" && dataFile) {
            formData.append("data_file", dataFile);
            formData.append("positive_keyword", positiveKeyword.trim());
            formData.append("negative_keyword", negativeKeyword.trim());
        }
        if (inputMode === "separate" && positiveFile && negativeFile) {
            formData.append("positive_file", positiveFile);
            formData.append("negative_file", negativeFile);
        }
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
            const response = await apiClient.post("/analyze-kmers", formData, {
                headers: { "Content-Type": "multipart/form-data" },
            });
            setTaskId(response.data.task_id);
        }
        catch (submitError) {
            const message = submitError instanceof Error ? submitError.message : "Unable to start analysis";
            setLoading(false);
            setError(message);
        }
    };
    const runs = useMemo(() => result?.runs ?? [], [result]);
    const proteomeMappingState = result && runs[0]
        ? {
            resultId: result.result_id,
            positiveFilename: runs[0].ad_filename,
            negativeFilename: runs[0].nc_filename,
            inputMode,
            positiveKeyword,
            negativeKeyword,
        }
        : undefined;
    const renderUploadCard = (inputId, label, helperText, selectedFileName, onFileChange) => (_jsxs("div", { className: "upload-card kmer-upload-card", children: [_jsxs("div", { className: "upload-card__top", children: [_jsxs("div", { children: [_jsx("span", { className: "upload-card__label", children: label }), _jsx("p", { className: "upload-card__helper", children: helperText })] }), _jsx("label", { className: "upload-trigger", htmlFor: inputId, children: "Choose file" })] }), _jsx("input", { id: inputId, className: "sr-only-file-input", type: "file", accept: ".xlsx,.csv", onChange: (event) => {
                    const files = event.target.files;
                    const selectedFile = files && files.length > 0 ? files[0] : null;
                    onFileChange(selectedFile);
                } }), _jsx("div", { className: `upload-file-display ${selectedFileName ? "upload-file-display--selected" : ""}`, children: selectedFileName ?? "No file selected" })] }));
    return (_jsxs("div", { className: "page-grid", children: [_jsxs(SectionCard, { title: "K-mer enrichment analysis", description: "Build the complete k-mer universe, apply product-based chi-square filtering to each sample, and run Mann\u2013Whitney U tests.", children: [_jsxs("form", { className: "form-grid", onSubmit: handleSubmit, children: [_jsxs("fieldset", { className: "form-field input-mode-group", children: [_jsx("legend", { children: "Input type" }), _jsxs("div", { className: "input-mode-options", children: [_jsxs("label", { className: `input-mode-option ${inputMode === "merged" ? "input-mode-option--active" : ""}`, children: [_jsx("input", { type: "radio", name: "kmer-input-mode", value: "merged", checked: inputMode === "merged", onChange: () => setInputMode("merged") }), _jsx("span", { children: "Merged cohort file" }), _jsx("small", { children: "Split columns using positive and negative keywords." })] }), _jsxs("label", { className: `input-mode-option ${inputMode === "separate" ? "input-mode-option--active" : ""}`, children: [_jsx("input", { type: "radio", name: "kmer-input-mode", value: "separate", checked: inputMode === "separate", onChange: () => setInputMode("separate") }), _jsx("span", { children: "Separate cohort files" }), _jsx("small", { children: "Upload positive and negative cohort files directly." })] })] })] }), inputMode === "merged" ? (_jsxs(_Fragment, { children: [renderUploadCard("kmer-merged-file", "Merged cohort file", "Original file should contain alternating cohort columns with sequence/count pairs.", dataFileName, (selectedFile) => {
                                        setDataFile(selectedFile);
                                        setDataFileName(selectedFile?.name ?? null);
                                    }), _jsxs("div", { className: "settings-grid", children: [_jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Positive column keyword" }), _jsx("input", { type: "text", value: positiveKeyword, onChange: (event) => setPositiveKeyword(event.target.value), placeholder: "AD" }), _jsx("small", { children: "Columns starting with this keyword are treated as the positive cohort." })] }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Negative column keyword" }), _jsx("input", { type: "text", value: negativeKeyword, onChange: (event) => setNegativeKeyword(event.target.value), placeholder: "NC" }), _jsx("small", { children: "Columns starting with this keyword are treated as the negative cohort." })] })] })] })) : (_jsxs("div", { className: "settings-grid", children: [renderUploadCard("kmer-positive-file", "Positive cohort file", "Upload a file containing positive cohort sequence/count columns.", positiveFileName, (selectedFile) => {
                                        setPositiveFile(selectedFile);
                                        setPositiveFileName(selectedFile?.name ?? null);
                                    }), renderUploadCard("kmer-negative-file", "Negative cohort file", "Upload a file containing negative cohort sequence/count columns.", negativeFileName, (selectedFile) => {
                                        setNegativeFile(selectedFile);
                                        setNegativeFileName(selectedFile?.name ?? null);
                                    })] })), _jsx("div", { className: "form-field form-field--inline", children: _jsxs("label", { children: [_jsx("span", { children: "k-mer size" }), _jsx("input", { type: "text", inputMode: "numeric", pattern: "[0-9]*", value: kValue, onChange: (event) => setKValue(event.target.value), placeholder: "4" }), _jsx("small", { children: "Enter one integer between 4 and 10." })] }) }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Download bundle name" }), _jsx("input", { type: "text", value: archiveName, onChange: (event) => setArchiveName(event.target.value), placeholder: "e.g. patient_cohort_k4_results" }), _jsx("small", { children: "Optional. Used as the ZIP filename when users download the result bundle." })] }), _jsxs("label", { className: "form-field", children: [_jsx("span", { children: "Wildcard positions" }), _jsx("input", { type: "text", placeholder: "e.g. 1,3", value: wildcards, onChange: (event) => setWildcards(event.target.value) }), _jsx("small", { children: "Leave empty for strict k-mers. Use unique, comma-separated, zero-indexed positions." })] }), _jsxs("label", { className: "form-field form-checkbox", children: [_jsx("input", { type: "checkbox", checked: normalize, onChange: (event) => setNormalize(event.target.checked) }), _jsx("span", { children: "Normalize the matrix and statistical values to proportions" }), _jsx("small", { children: "Uses each sample's pre-filter tiled total. Turn this off to keep raw filtered counts." })] }), _jsxs("div", { className: "form-actions", children: [_jsx("button", { type: "submit", className: "primary-button", disabled: loading, children: loading ? "Running analysis..." : "Start analysis" }), result && _jsx(DownloadButton, { resultId: result.result_id, label: "Download CSV bundle" })] })] }), taskStatus && loading && (_jsxs("div", { className: "progress-panel", "aria-live": "polite", children: [_jsxs("div", { className: "progress-panel__header", children: [_jsx("strong", { children: taskStatus.message }), _jsxs("span", { children: [taskStatus.progress, "%"] })] }), _jsx("div", { className: "progress-bar", role: "progressbar", "aria-valuenow": taskStatus.progress, "aria-valuemin": 0, "aria-valuemax": 100, children: _jsx("div", { className: "progress-bar__fill", style: { width: `${taskStatus.progress}%` } }) })] })), error && _jsx(StatusBanner, { tone: "error", title: "K-mer analysis failed", message: error }), result && !error && (_jsx(StatusBanner, { tone: "success", title: "K-mer analysis ready", message: "Download CSV exports for the completed k-mer run." }))] }), runs.length > 0 && result && (_jsx(SectionCard, { title: "Mann\u2013Whitney summaries", description: "Highlights the number of kmers enriched in each cohort for the selected k value.", actions: _jsxs("div", { className: "section-card__actions-group", children: [_jsx(Link, { to: "/module3", state: proteomeMappingState, className: "secondary-button", children: "Open Proteome Mapping" }), _jsx(DownloadButton, { resultId: result.result_id, label: "Download combined outputs" })] }), children: _jsx("div", { className: "table-wrapper", children: _jsxs("table", { children: [_jsx("thead", { children: _jsxs("tr", { children: [_jsx("th", { children: "k-mer size" }), _jsx("th", { children: "Total kmers" }), _jsx("th", { children: "Positive-elevated" }), _jsx("th", { children: "Negative-elevated" }), _jsx("th", { children: "Files" })] }) }), _jsx("tbody", { children: runs.map((run) => (_jsxs("tr", { children: [_jsx("td", { children: run.k }), _jsx("td", { children: run.total_kmers.toLocaleString() }), _jsx("td", { children: run.ad_elevated.toLocaleString() }), _jsx("td", { children: run.nc_elevated.toLocaleString() }), _jsx("td", { children: _jsxs("ul", { className: "file-list", children: [_jsx("li", { children: run.result_filename }), _jsx("li", { children: run.ad_filename }), _jsx("li", { children: run.nc_filename }), _jsx("li", { children: run.matrix_filename })] }) })] }, run.k))) })] }) }) }))] }));
}
