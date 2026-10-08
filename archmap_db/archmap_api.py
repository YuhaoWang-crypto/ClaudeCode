"""Minimal client for the public ArchMap backend (https://www.archmap.bio).

The web app (React SPA) talks to a REST backend.  Anonymous visitors get a
short-lived JWT from ``/temp_auth``; with it, the same endpoints the
"References -> Atlases" page uses are readable:

    GET  /atlases                      reference-atlas metadata
    GET  /models                       mapping models (scVI / scANVI / scPoli)
    GET  /scvi-atlases                 scvi-hub (HuggingFace) atlases
    POST /file_download/atlas_files    {atlasId} -> [{fileName, presignedUrl}]

Presigned URLs are GCS signed links valid for 7 days; they are never
written to the database or the committed snapshots.
"""
from __future__ import annotations

import json
import re
import subprocess
import time
import urllib.request

API = "https://custom-helix-329116.ey.r.appspot.com/v1"


class ArchMap:
    def __init__(self, api: str = API, retries: int = 4):
        self.api = api
        self.retries = retries
        self._jwt = None

    def _request(self, path, data=None, auth=True):
        url = f"{self.api}{path}"
        headers = {"Content-Type": "application/json"}
        if auth:
            headers["Authorization"] = f"Bearer {self.jwt}"
        body = json.dumps(data).encode() if data is not None else None
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(url, data=body, headers=headers)
                with urllib.request.urlopen(req, timeout=120) as r:
                    return json.loads(r.read())
            except Exception:
                if attempt == self.retries - 1:
                    raise
                time.sleep(2 ** (attempt + 1))

    @property
    def jwt(self):
        if self._jwt is None:
            self._jwt = self._request("/temp_auth", auth=False)["jwt"]
        return self._jwt

    def atlases(self):
        return self._request("/atlases")

    def models(self):
        return self._request("/models")

    def scvi_atlases(self):
        return self._request("/scvi-atlases")

    def atlas_files(self, atlas_id):
        """[{fileName, presignedUrl}] for one atlas (data, counts, model)."""
        return self._request("/file_download/atlas_files", {"atlasId": atlas_id})


def remote_size(url: str) -> int | None:
    """Size of a presigned object via a 1-byte range GET (HEAD is not signed)."""
    out = subprocess.run(
        ["curl", "-sS", "-m", "60", "-r", "0-0", "-D", "-", "-o", "/dev/null", url],
        capture_output=True, text=True,
    ).stdout
    m = re.search(r"content-range: bytes 0-0/(\d+)", out, re.I)
    return int(m.group(1)) if m else None
