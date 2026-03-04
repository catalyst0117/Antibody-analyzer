import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
import { Link } from "react-router-dom";
const NAV_ITEMS = [
    { label: "FASTQ Processor", path: "/fastq" },
    { label: "K-mer Analysis", path: "/kmer" }
];
export function HeaderNav({ currentPath }) {
    return (_jsxs("header", { className: "app-header", children: [_jsxs("div", { className: "brand", children: [_jsx("span", { className: "brand-logo", "aria-hidden": true, children: "\uD83E\uDDEC" }), _jsxs("div", { children: [_jsx("h1", { children: "Antibody Analyzer" }), _jsx("p", { className: "brand-tagline", children: "Web interface for the sequencing modules" })] })] }), _jsx("nav", { className: "app-nav", children: NAV_ITEMS.map((item) => {
                    const active = currentPath.startsWith(item.path);
                    return (_jsx(Link, { to: item.path, className: active ? "nav-link active" : "nav-link", children: item.label }, item.path));
                }) })] }));
}
