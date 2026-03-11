from __future__ import annotations

import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Dict, Optional

from app.models.responses import KmerResponse


@dataclass
class KmerTask:
    task_id: str
    status: str
    progress: int
    message: str
    created_at: float
    result: Optional[KmerResponse] = None
    error: Optional[str] = None


class KmerTaskStore:
    def __init__(self, max_age_seconds: int = 60 * 60):
        self.max_age_seconds = max_age_seconds
        self._lock = threading.Lock()
        self._tasks: Dict[str, KmerTask] = {}

    def _cleanup_locked(self) -> None:
        now = time.time()
        expired = [
            task_id
            for task_id, task in self._tasks.items()
            if now - task.created_at > self.max_age_seconds
        ]
        for task_id in expired:
            self._tasks.pop(task_id, None)

    def create_task(self, message: str = "Queued analysis") -> str:
        with self._lock:
            self._cleanup_locked()
            task_id = uuid.uuid4().hex
            self._tasks[task_id] = KmerTask(
                task_id=task_id,
                status="queued",
                progress=0,
                message=message,
                created_at=time.time(),
            )
        return task_id

    def update_task(self, task_id: str, *, status: Optional[str] = None, progress: Optional[int] = None, message: Optional[str] = None) -> Optional[KmerTask]:
        with self._lock:
            task = self._tasks.get(task_id)
            if not task:
                return None
            if status is not None:
                task.status = status
            if progress is not None:
                task.progress = max(0, min(100, progress))
            if message is not None:
                task.message = message
            return task

    def complete_task(self, task_id: str, result: KmerResponse) -> Optional[KmerTask]:
        with self._lock:
            task = self._tasks.get(task_id)
            if not task:
                return None
            task.status = "succeeded"
            task.progress = 100
            task.message = "Analysis complete"
            task.result = result
            task.error = None
            return task

    def fail_task(self, task_id: str, error: str) -> Optional[KmerTask]:
        with self._lock:
            task = self._tasks.get(task_id)
            if not task:
                return None
            task.status = "failed"
            task.message = "Analysis failed"
            task.error = error
            return task

    def get_task(self, task_id: str) -> Optional[KmerTask]:
        with self._lock:
            return self._tasks.get(task_id)


__all__ = ["KmerTask", "KmerTaskStore"]
