"""Retrying HTTP session shared by every public-data client.

The public biomedical APIs used here (TCIA/NBIA, GDC, NCBI E-utilities) all
throttle aggressively and drop connections under load, so every request goes
through the same bounded exponential-backoff wrapper.
"""

from __future__ import annotations

import time
from typing import Any

import requests

DEFAULT_TIMEOUT = 180
DEFAULT_RETRIES = 5
USER_AGENT = "mri-cd8-til/0.1 (public-data inventory; contact via repository)"

#: Statuses worth retrying: throttling plus transient gateway failures.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})


class PublicDataSession:
    """A :class:`requests.Session` with retry/backoff and a polite User-Agent."""

    def __init__(
        self,
        retries: int = DEFAULT_RETRIES,
        timeout: int = DEFAULT_TIMEOUT,
        backoff_base: float = 2.0,
    ) -> None:
        self.retries = retries
        self.timeout = timeout
        self.backoff_base = backoff_base
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": USER_AGENT})

    def get(self, url: str, **kwargs: Any) -> requests.Response:
        kwargs.setdefault("timeout", self.timeout)
        last_exc: Exception | None = None
        for attempt in range(self.retries):
            try:
                response = self._session.get(url, **kwargs)
                if response.status_code in RETRY_STATUS:
                    raise requests.HTTPError(
                        f"retryable status {response.status_code}", response=response
                    )
                response.raise_for_status()
                return response
            except Exception as exc:  # noqa: BLE001 - deliberately broad, we retry
                last_exc = exc
                if attempt == self.retries - 1:
                    break
                time.sleep(self.backoff_base**attempt)
        assert last_exc is not None
        raise last_exc

    def get_json(self, url: str, **kwargs: Any) -> Any:
        return self.get(url, **kwargs).json()

    def get_text(self, url: str, **kwargs: Any) -> str:
        return self.get(url, **kwargs).text
