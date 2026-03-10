import axios from "axios";

const baseURL = import.meta.env.VITE_API_BASE_URL ?? "http://localhost:8000/api";

export const apiClient = axios.create({
  baseURL,
  headers: {
    "X-Requested-With": "XMLHttpRequest",
  },
});
