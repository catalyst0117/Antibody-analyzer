from __future__ import annotations

import shutil
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional


@dataclass
class StoredResult:
    summary: dict
    archive_path: Path
    download_name: str
    created_at: float


class ResultStore:
    """Stores generated analysis artifacts on disk and exposes download handles."""

    def __init__(self, base_dir: Optional[Path] = None, max_age_seconds: int = 60 * 60):
        self.base_dir = Path(base_dir or tempfile.mkdtemp(prefix="antibody_results_"))
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.max_age_seconds = max_age_seconds
        self._lock = threading.Lock()
        self._results: Dict[str, StoredResult] = {}

    def _cleanup_locked(self) -> None:
        expired: List[str] = []
        now = time.time()
        for result_id, result in list(self._results.items()):
            if now - result.created_at > self.max_age_seconds:
                expired.append(result_id)
        for result_id in expired:
            stored = self._results.pop(result_id, None)
            if stored:
                stored.archive_path.unlink(missing_ok=True)

    def _normalize_download_name(self, download_name: Optional[str]) -> str:
        name = (download_name or "results").strip()
        if name.lower().endswith(".zip"):
            name = name[:-4]
        name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
        if not name:
            name = "results"
        return f"{name}.zip"

    def create_result(
        self,
        summary: dict,
        files: Iterable[Path],
        download_name: Optional[str] = None,
        archive_root_name: Optional[str] = None,
    ) -> str:
        archive_dir = self.base_dir / uuid.uuid4().hex
        archive_dir.mkdir(parents=True, exist_ok=True)
        bundle_root_name = self._normalize_download_name(archive_root_name or download_name).removesuffix(".zip")
        bundle_root_dir = archive_dir / bundle_root_name
        bundle_root_dir.mkdir(parents=True, exist_ok=True)
        for file_path in files:
            shutil.copy(file_path, bundle_root_dir / file_path.name)
        archive_path = shutil.make_archive(
            str(archive_dir),
            "zip",
            root_dir=archive_dir,
            base_dir=bundle_root_dir.name,
        )
        shutil.rmtree(archive_dir, ignore_errors=True)

        with self._lock:
            self._cleanup_locked()
            result_id = uuid.uuid4().hex
            self._results[result_id] = StoredResult(
                summary=summary,
                archive_path=Path(archive_path),
                download_name=self._normalize_download_name(download_name),
                created_at=time.time(),
            )
        return result_id

    def get_summary(self, result_id: str) -> Optional[dict]:
        with self._lock:
            stored = self._results.get(result_id)
            if not stored:
                return None
            return stored.summary

    def get_archive_path(self, result_id: str) -> Optional[Path]:
        with self._lock:
            stored = self._results.get(result_id)
            if not stored:
                return None
            return stored.archive_path

    def get_download_name(self, result_id: str) -> Optional[str]:
        with self._lock:
            stored = self._results.get(result_id)
            if not stored:
                return None
            return stored.download_name


__all__ = ["ResultStore"]
