import { jsx as _jsx, jsxs as _jsxs } from "react/jsx-runtime";
export function SectionCard({ title, description, actions, children }) {
    return (_jsxs("section", { className: "section-card", children: [_jsxs("header", { className: "section-card__header", children: [_jsxs("div", { children: [_jsx("h2", { children: title }), description && _jsx("p", { className: "section-card__description", children: description })] }), actions && _jsx("div", { className: "section-card__actions", children: actions })] }), _jsx("div", { className: "section-card__content", children: children })] }));
}
