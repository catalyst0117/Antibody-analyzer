import { SectionCard } from "../components/SectionCard";

export function AboutPage() {
  return (
    <div className="page-grid about-page">
      <SectionCard
        title="About the Glabe Lab"
        description="Research context and project information for the Antibody Target Prediction Tool."
      >
        <div className="about-copy">
          <p>
            The Glabe Lab is a research group dedicated to studying disease-correlated antibodies.
            Led by Professor Charles Glabe, the K-mer Project began in 2016 to study disease
            populations and catalog potentially significant antibody sequences.
          </p>
          <p>
            The Antibody Target Prediction Tool analyzes patient sequencing data by measuring the
            prevalence of k-mers, which are amino-acid sequences of a defined length. Significant
            k-mers are then cross-referenced against a reference proteome to identify candidate
            proteins for further biological investigation.
          </p>
          <p>
            The tool is intended to support research into proteins and antibody responses that may
            help clarify disease mechanisms or inform future therapeutic studies.
          </p>
          <p>
            The project has been extended from tetramer-only analysis to support k-mers of any
            length, with the analysis modules integrated into the full-stack application and
            improvements to proteome mapping and result reporting.
          </p>
        </div>
      </SectionCard>

      <SectionCard title="Project Contact" description="Technical and project information.">
        <div className="about-contact">
          <p><strong>Developer:</strong> Chengyan Zhao</p>
          <p>
            <strong>Email:</strong>{" "}
            <a href="mailto:czhao27@uci.edu">czhao27@uci.edu</a>
          </p>
          <p>
            <strong>Repository:</strong>{" "}
            <a href="https://github.com/catalyst0117/Antibody-analyzer" target="_blank" rel="noreferrer">
              github.com/catalyst0117/Antibody-analyzer
            </a>
          </p>
          <p>
            For lab-related questions, contact the Glabe Lab at{" "}
            <a href="mailto:glabelab@uci.edu">glabelab@uci.edu</a>.
          </p>
        </div>
      </SectionCard>
    </div>
  );
}