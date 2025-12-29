import { PropsWithChildren } from "react";

type SectionCardProps = PropsWithChildren<{
  title: string;
  description?: string;
  actions?: React.ReactNode;
}>;

export function SectionCard({ title, description, actions, children }: SectionCardProps) {
  return (
    <section className="section-card">
      <header className="section-card__header">
        <div>
          <h2>{title}</h2>
          {description && <p className="section-card__description">{description}</p>}
        </div>
        {actions && <div className="section-card__actions">{actions}</div>}
      </header>
      <div className="section-card__content">{children}</div>
    </section>
  );
}
