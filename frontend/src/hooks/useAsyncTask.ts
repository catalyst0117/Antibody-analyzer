import { useCallback, useState } from "react";

type AsyncState<T> = {
  execute: (...args: any[]) => Promise<T | undefined>;
  loading: boolean;
  error: string | null;
};

export function useAsyncTask<T>(task: (...args: any[]) => Promise<T>): AsyncState<T> {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const execute = useCallback(
    async (...args: any[]) => {
      setLoading(true);
      setError(null);
      try {
        const result = await task(...args);
        return result;
      } catch (err) {
        const message = err instanceof Error ? err.message : "Unexpected error";
        setError(message);
        return undefined;
      } finally {
        setLoading(false);
      }
    },
    [task]
  );

  return { execute, loading, error };
}
