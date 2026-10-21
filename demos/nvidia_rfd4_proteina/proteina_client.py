"""Minimal client for NVIDIA's hosted RFD4-Proteina NIM endpoint.

Standard library only -- no pip install needed. The endpoint is synchronous:
it holds the connection open (see NVCF-POLL-SECONDS) and returns the finished
designs in the response body, so there is no job id to poll.

API key is read from the NVIDIA_API_KEY environment variable, or from a
``.env`` file next to this module. The key is never written to stdout.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

ENDPOINT = "https://health.api.nvidia.com/v1/biology/ipd/rfd4-proteina/generative"

# Statuses worth another attempt. 429 is the common one: the free hosted
# endpoint is rate limited per key, and parallel calls trip it immediately.
RETRY_STATUSES = frozenset({408, 409, 425, 429, 500, 502, 503, 504})

TASKS = (
    "generation",
    "motif_scaffolding",
    "enzyme_design",
    "small_molecule_binder",
    "antibody_design",
    "partial_diffusion",
    "inverse_folding",
)

RECYCLE_STRATEGIES = ("thinking_low", "thinking_medium", "thinking_high")


class ProteinaError(RuntimeError):
    """An error reported by the endpoint, or a transport failure."""

    def __init__(self, message: str, status: int | None = None, body: str = ""):
        super().__init__(message)
        self.status = status
        self.body = body


def load_api_key(env_file: str | os.PathLike[str] | None = None) -> str:
    """Return the API key from the environment, falling back to a .env file."""
    key = os.environ.get("NVIDIA_API_KEY", "").strip()
    if key:
        return key

    candidates = [Path(env_file)] if env_file else [Path(__file__).with_name(".env")]
    for path in candidates:
        if not path.is_file():
            continue
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.strip().removeprefix("export ").strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            name, _, value = line.partition("=")
            if name.strip() == "NVIDIA_API_KEY":
                return value.strip().strip("'\"")

    raise ProteinaError(
        "No API key found. Set NVIDIA_API_KEY, or put it in a .env file "
        "next to proteina_client.py (see .env.example). Get a key at "
        "https://build.nvidia.com/ipd/rfd4-proteina"
    )


def inline_file(path: str | os.PathLike[str]) -> dict[str, str]:
    """Read an mmCIF file into the endpoint's InlineFile shape.

    The endpoint takes file *contents*, not paths, and rejects PDB format --
    it wants a complete mmCIF with an ``_atom_site`` loop.
    """
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if "_atom_site" not in text:
        raise ProteinaError(f"{p} has no _atom_site loop; is it really mmCIF?")
    return {"filename": p.name, "content": text}


def build_generation_request(
    template: dict[str, str],
    *,
    example_id: str = "design",
    replicates: int = 1,
    seed: int | None = None,
    recycle_strategy: str = "thinking_low",
    num_steps: int | None = None,
    metric_tags: tuple[str, ...] = ("design",),
) -> dict[str, Any]:
    """Unconditional backbone + sequence generation from an annotated mmCIF.

    The ``annotated_cif`` dialect carries the design conditions inside the
    file itself, and the design length follows the template's token count.
    """
    sampling: dict[str, Any] = {
        "num_replicates": replicates,
        "recycle_strategy": recycle_strategy,
    }
    if seed is not None:
        sampling["seed"] = seed
    if num_steps is not None:
        sampling["num_steps"] = num_steps
        # num_steps and recycle_strategy coexist; recycle_schedule is the
        # mutually exclusive one.

    return {
        "task": "generation",
        "input": {
            "format": "annotated_cif",
            "example_id": example_id,
            "file": template,
        },
        "sampling": sampling,
        "metrics": {"tags": list(metric_tags)},
        "output": {"structure_format": "mmcif", "include_sequence": True},
    }


def build_scaffold_request(
    template: dict[str, str],
    *,
    contig: str,
    motif_select: str,
    example_id: str = "scaffold",
    replicates: int = 1,
    seed: int | None = None,
    recycle_strategy: str = "thinking_low",
    metric_tags: tuple[str, ...] = ("design",),
) -> dict[str, Any]:
    """Motif scaffolding via the ``condition_spec`` dialect.

    ``contig`` is a comma-separated string of segments: a bare number is a
    designed stretch of that length, ``A17-131`` keeps those residues of
    chain A, and ``/0`` is a chain break. The condition spec must name an
    input structure, and that structure only exists to be referenced by the
    contig's preserved segments -- a contig that designs everything from
    scratch belongs in ``build_generation_request`` instead.
    """
    return {
        "task": "motif_scaffolding",
        "input": {
            "format": "condition_spec",
            "example_id": example_id,
            "spec": {
                "input": template["filename"],
                "contig": contig,
                # C_CRDGRP group 1: hold the motif's backbone coordinates.
                "conditions": {"C_CRDGRP": {"1": [{"select": motif_select}]}},
            },
            "files": [template],
        },
        "sampling": _sampling(replicates, seed, recycle_strategy),
        "metrics": {"tags": list(metric_tags)},
        "output": {"structure_format": "mmcif", "include_sequence": True},
    }


def _sampling(replicates: int, seed: int | None, recycle_strategy: str) -> dict[str, Any]:
    sampling: dict[str, Any] = {
        "num_replicates": replicates,
        "recycle_strategy": recycle_strategy,
    }
    if seed is not None:
        sampling["seed"] = seed
    return sampling


def motif_selection(chain: str, first: int, last: int) -> str:
    """Backbone-atom selection for a residue range.

    Residue ranges must be bracketed: ``sel('A/*/[1-15]/[N,CA,C,O]')``.
    """
    return f"sel('{chain}/*/[{first}-{last}]/[N,CA,C,O]')"


def invoke(
    payload: dict[str, Any],
    *,
    api_key: str,
    timeout: float = 900.0,
    poll_seconds: int = 300,
    max_attempts: int = 6,
    base_backoff: float = 20.0,
    include_runtime_metrics: bool = False,
    on_retry=None,
) -> dict[str, Any]:
    """POST one design request and return the parsed response.

    Retries throttling and transient server errors with exponential backoff.
    """
    body = json.dumps(payload).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        # Ask NVCF to hold the connection rather than hand back a polling job.
        "NVCF-POLL-SECONDS": str(poll_seconds),
    }
    if include_runtime_metrics:
        headers["X-Include-Runtime-Metrics"] = "true"

    last_error: ProteinaError | None = None
    for attempt in range(1, max_attempts + 1):
        request = urllib.request.Request(ENDPOINT, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            raw = exc.read().decode("utf-8", "replace")
            last_error = ProteinaError(_explain(exc.code, raw), status=exc.code, body=raw)
            retriable = exc.code in RETRY_STATUSES
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = ProteinaError(f"transport error: {exc}")
            retriable = True

        if not retriable or attempt == max_attempts:
            raise last_error

        delay = base_backoff * (2 ** (attempt - 1))
        if on_retry:
            on_retry(attempt, max_attempts, delay, last_error)
        time.sleep(delay)

    raise last_error if last_error else ProteinaError("request failed")


def _explain(status: int, raw: str) -> str:
    """Turn an error body into one readable line."""
    detail = raw.strip()
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        parsed = None
    if isinstance(parsed, dict):
        err = parsed.get("error")
        if isinstance(err, dict):
            detail = f"{err.get('code', 'error')}: {err.get('message', '')}".strip()
        elif isinstance(err, str):
            detail = err
        elif parsed.get("title"):
            detail = str(parsed["title"])

    hint = {
        401: " (key rejected -- check NVIDIA_API_KEY)",
        403: " (key is malformed, expired, or lacks access to this model)",
        404: " (wrong endpoint path)",
        422: " (request body did not match the schema)",
        429: " (rate limited -- the free endpoint allows very little concurrency)",
    }.get(status, "")
    return f"HTTP {status}{hint}: {detail[:600]}"


def check_credentials(api_key: str, timeout: float = 30.0) -> str:
    """Probe the endpoint with a deliberately empty body.

    A 422 is the success case: it proves the key authenticated and the route
    exists, without spending GPU time on a design.
    """
    try:
        invoke({}, api_key=api_key, timeout=timeout, max_attempts=1)
    except ProteinaError as exc:
        if exc.status == 422:
            return "ok: endpoint reachable and key accepted (empty body correctly rejected)"
        raise
    return "ok: endpoint accepted an empty body (unexpected, but reachable)"
