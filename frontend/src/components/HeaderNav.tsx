import { Link } from "react-router-dom";

type HeaderNavProps = {
  currentPath: string;
};

const NAV_ITEMS = [
  { label: "FASTQ Processor", path: "/fastq" },
  { label: "K-mer Analysis", path: "/kmer" }
];

export function HeaderNav({ currentPath }: HeaderNavProps) {
  return (
    <header className="app-header">
      <div className="brand">
        <span className="brand-logo" aria-hidden>🧬</span>
        <div>
          <h1>Antibody Analyzer</h1>
          <p className="brand-tagline">Web interface for the sequencing modules</p>
        </div>
      </div>
      <nav className="app-nav">
        {NAV_ITEMS.map((item) => {
          const active = currentPath.startsWith(item.path);
          return (
            <Link key={item.path} to={item.path} className={active ? "nav-link active" : "nav-link"}>
              {item.label}
            </Link>
          );
        })}
      </nav>
    </header>
  );
}
