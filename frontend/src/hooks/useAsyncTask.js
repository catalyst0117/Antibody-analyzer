import { useCallback, useState } from "react";
export function useAsyncTask(task) {
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const execute = useCallback(async (...args) => {
        setLoading(true);
        setError(null);
        try {
            const result = await task(...args);
            return result;
        }
        catch (err) {
            const message = err instanceof Error ? err.message : "Unexpected error";
            setError(message);
            return undefined;
        }
        finally {
            setLoading(false);
        }
    }, [task]);
    return { execute, loading, error };
}
