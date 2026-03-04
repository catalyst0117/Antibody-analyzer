import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
export function StatusBanner({ tone, title, message }) {
    return (_jsxs("div", { className: `status-banner status-${tone}`, role: "status", children: [_jsx("strong", { children: title }), message && _jsx("span", { children: message })] }));
}
