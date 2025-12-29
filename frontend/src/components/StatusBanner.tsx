type StatusBannerProps = {
  tone: "info" | "success" | "error";
  title: string;
  message?: string;
};

export function StatusBanner({ tone, title, message }: StatusBannerProps) {
  return (
    <div className={`status-banner status-${tone}`} role="status">
      <strong>{title}</strong>
      {message && <span>{message}</span>}
    </div>
  );
}
