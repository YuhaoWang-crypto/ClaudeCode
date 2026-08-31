"""The campaign submit gate — the single mechanical dispatch authority.

This is the module the prompt's SETUP sub-agent writes to
``/state/lib/submit_gate.py`` and every submitter imports.  One ``True``
authorises exactly one container start; never per-batch, never per-loop.

The gate is **FAIL-CLOSED**: the entire body of :func:`submit_gate` is wrapped in
``try/except`` and any exception returns a falsey decision with
``reason="GATE_ERROR"``.  A delayed dispatch costs seconds; an ungated H100
dispatch costs unrecoverable dollars.

Steps, in the prompt's order:

======  =====================================================================
(a)     read ``governor.json`` fresh (never a once-at-start copy; a copy cached
        < 60 s by file mtime is fresh)
(b)     missing, or ``utc_now - set_at > 25 min`` -> ceiling 0,
        ``STALE_OR_MISSING_GOVERNOR`` -- the submitter-side dead-man switch
(c)     ``ceiling == 0`` -> False immediately, no retry
(d)     ``timeout_s > 1800`` and ``max(pct_spent_metered, pct_spent_projected)
        > pct_elapsed + 5`` -> ``LONG_TIMEOUT_OVER_PACE``
(d2)    ``stage in {generate, gen_screen}`` with a target -> read
        ``ledger_agg.json`` fresh; ``gen >= 200 and gen > 2*max(scr, 50)``
        -> ``GEN_BACKLOG_CAP``, no retry
(e)     count live GPU sandboxes = all tagged sandboxes minus tagged-CPU minus
        the ``cpu_sandboxes.jsonl`` allowlist, via the Modal SDK (never the
        harness usage table), rate-limited to one SDK call per 8 s per process
(f)     ``live >= ceiling - max(10, floor(0.1*ceiling))`` -> jittered
        exponential backoff, False after 8 attempts
(g)     True only when (e) and (f) pass
======  =====================================================================

The validation-gate check (``/state/gates/{target}.json`` status PASS) is
applied here too when a scoring job passes ``target=``, per the prompt: a verbal
claim, Slack post, or in-memory variable does not satisfy the gate.

The two side-effecting dependencies -- reading ``/state`` and listing the live
fleet -- are injected (:class:`StateReader`, :class:`FleetCounter`) so the gate
logic is testable without Modal.  :class:`ModalStateReader` /
:class:`ModalFleetCounter` are the production implementations and import
``modal`` lazily.
"""

from __future__ import annotations

import json
import math
import os
import random
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Protocol

__all__ = [
    "GateDecision",
    "StateReader",
    "FleetCounter",
    "LocalStateReader",
    "ModalStateReader",
    "ModalFleetCounter",
    "configure",
    "submit_gate",
]

GOVERNOR_STALE_SECONDS = 25 * 60
LONG_TIMEOUT_S = 1800
PACE_OVERSHOOT_PP = 5.0
GEN_BACKLOG_MIN_GEN = 200
GEN_BACKLOG_SCR_FLOOR = 50
SDK_CALL_MIN_INTERVAL_S = 8.0
MAX_BACKOFF_ATTEMPTS = 8
BACKOFF_LADDER_S = (5.0, 10.0, 20.0, 40.0, 60.0, 60.0, 60.0)


@dataclass(frozen=True)
class GateDecision:
    """Truthy/falsey gate result that also carries the reason for the log.

    ``bool(decision)`` is the ``True``/``False`` the prompt specifies, so
    ``if not submit_gate(...): return`` behaves exactly as written, while
    ``decision.reason`` gives the submitter something to put in the ledger.
    """

    allowed: bool
    reason: str = "OK"
    live: int | None = None
    ceiling: int | None = None
    attempts: int = 0

    def __bool__(self) -> bool:  # pragma: no cover - trivial
        return self.allowed


# --------------------------------------------------------------------------- #
# injected dependencies
# --------------------------------------------------------------------------- #


class StateReader(Protocol):
    """Reads small JSON files from the campaign state volume."""

    def read_json(self, path: str) -> dict | None:
        """Return the parsed file, or ``None`` if it does not exist."""

    def read_lines(self, path: str) -> list[str]:
        """Return the file's lines, or ``[]`` if it does not exist."""


class FleetCounter(Protocol):
    """Counts live sandboxes for the campaign app."""

    def count_all(self, app_id: str, project_tag: str) -> int: ...

    def count_cpu_tagged(self, app_id: str, project_tag: str) -> int: ...


class LocalStateReader:
    """State reader backed by a local directory (tests, and the in-sandbox
    mount where ``/state`` *is* a local path)."""

    def __init__(self, root: str = "/state", cache_ttl_s: float = 60.0):
        self.root = root
        self.cache_ttl_s = cache_ttl_s
        self._cache: dict[str, tuple[float, float, object]] = {}

    def _path(self, path: str) -> str:
        return os.path.join(self.root, path.lstrip("/"))

    def _read_text(self, path: str) -> str | None:
        p = self._path(path)
        try:
            mtime = os.path.getmtime(p)
        except OSError:
            self._cache.pop(path, None)
            return None
        now = time.time()
        hit = self._cache.get(path)
        # "a copy cached <60s by file mtime is fresh, never a once-at-start copy"
        if hit and hit[0] == mtime and (now - hit[1]) < self.cache_ttl_s:
            return hit[2]  # type: ignore[return-value]
        try:
            with open(p, "r", encoding="utf-8") as fh:
                text = fh.read()
        except OSError:
            return None
        self._cache[path] = (mtime, now, text)
        return text

    def read_json(self, path: str) -> dict | None:
        text = self._read_text(path)
        if text is None:
            return None
        return json.loads(text)

    def read_lines(self, path: str) -> list[str]:
        text = self._read_text(path)
        if text is None:
            return []
        return [ln for ln in text.splitlines() if ln.strip()]


class ModalStateReader:
    """Production state reader.

    Inside a Modal container with ``/state`` mounted this is just a local read.
    From a Claude Science kernel it goes through ONE long-lived CPU sandbox
    (``sb.exec('cat', path)``) — never ``vol.read_file`` from a compute_provider
    kernel, which the proxy answers with a 403 that no network-access grant
    fixes.
    """

    def __init__(self, sandbox=None, root: str = "/state", cache_ttl_s: float = 60.0):
        self.sandbox = sandbox
        self.root = root
        self.cache_ttl_s = cache_ttl_s
        self._local = LocalStateReader(root, cache_ttl_s)
        self._cache: dict[str, tuple[float, str | None]] = {}

    def _cat(self, path: str) -> str | None:
        full = os.path.join(self.root, path.lstrip("/"))
        if self.sandbox is None:
            return self._local._read_text(path)
        now = time.time()
        hit = self._cache.get(path)
        if hit and (now - hit[0]) < self.cache_ttl_s:
            return hit[1]
        proc = self.sandbox.exec("cat", full)
        out = proc.stdout.read()
        text = out if isinstance(out, str) else out.decode("utf-8")
        if proc.wait() != 0:
            text = None
        self._cache[path] = (now, text)
        return text

    def read_json(self, path: str) -> dict | None:
        text = self._cat(path)
        if not text:
            return None
        return json.loads(text)

    def read_lines(self, path: str) -> list[str]:
        text = self._cat(path)
        if not text:
            return []
        return [ln for ln in text.splitlines() if ln.strip()]


class ModalFleetCounter:
    """Live-sandbox count via the Modal SDK, rate-limited to one call per 8 s.

    Counts ALL tagged sandboxes minus the ``kind=cpu`` tagged ones, exactly as
    the prompt's step (e) specifies — never via the harness ``compute_usage``
    table.
    """

    def __init__(self, min_interval_s: float = SDK_CALL_MIN_INTERVAL_S):
        self.min_interval_s = min_interval_s
        self._cache: dict[tuple[str, str, bool], tuple[float, int]] = {}

    def _list_len(self, app_id: str, project_tag: str, cpu_only: bool) -> int:
        import modal  # lazy: the gate is importable without Modal installed

        key = (app_id, project_tag, cpu_only)
        now = time.time()
        hit = self._cache.get(key)
        if hit and (now - hit[0]) < self.min_interval_s:
            return hit[1]
        tags = {"claude-science-project": project_tag}
        if cpu_only:
            tags["kind"] = "cpu"
        n = len(list(modal.Sandbox.list(app_id=app_id, tags=tags)))
        self._cache[key] = (now, n)
        return n

    def count_all(self, app_id: str, project_tag: str) -> int:
        return self._list_len(app_id, project_tag, cpu_only=False)

    def count_cpu_tagged(self, app_id: str, project_tag: str) -> int:
        return self._list_len(app_id, project_tag, cpu_only=True)


# --------------------------------------------------------------------------- #
# module configuration (SETUP wires this once per process)
# --------------------------------------------------------------------------- #

_STATE: StateReader = LocalStateReader()
_FLEET: FleetCounter | None = None
_SLEEP = time.sleep
_NOW = lambda: datetime.now(timezone.utc)  # noqa: E731


def configure(
    state: StateReader | None = None,
    fleet: FleetCounter | None = None,
    sleep=None,
    now=None,
) -> None:
    """Wire the gate's dependencies.  Called once per submitter process."""
    global _STATE, _FLEET, _SLEEP, _NOW
    if state is not None:
        _STATE = state
    if fleet is not None:
        _FLEET = fleet
    if sleep is not None:
        _SLEEP = sleep
    if now is not None:
        _NOW = now


def _parse_utc(value) -> datetime | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value if value.tzinfo else value.replace(tzinfo=timezone.utc)
    text = str(value).strip().replace("Z", "+00:00")
    try:
        dt = datetime.fromisoformat(text)
    except ValueError:
        return None
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def _backoff_delay(attempt: int) -> float:
    """random.uniform(0,5)s then 5, 10, 20 ... capped 60, each x U(0.8, 1.2)."""
    if attempt == 0:
        return random.uniform(0.0, 5.0)
    base = BACKOFF_LADDER_S[min(attempt - 1, len(BACKOFF_LADDER_S) - 1)]
    return base * random.uniform(0.8, 1.2)


# --------------------------------------------------------------------------- #
# the gate
# --------------------------------------------------------------------------- #


def submit_gate(
    app_id: str,
    project_tag: str,
    timeout_s: float,
    gpu_type: str | None = None,
    target: str | None = None,
    stage: str | None = None,
) -> GateDecision:
    """Authorise exactly one container start.  Fail-closed on any exception."""
    try:
        return _submit_gate_inner(
            app_id, project_tag, timeout_s, gpu_type, target, stage
        )
    except Exception:  # noqa: BLE001 - the gate is never fail-open
        return GateDecision(False, "GATE_ERROR")


def _submit_gate_inner(
    app_id: str,
    project_tag: str,
    timeout_s: float,
    gpu_type: str | None,
    target: str | None,
    stage: str | None,
) -> GateDecision:
    now = _NOW()

    # (a) read governor.json fresh
    gov = _STATE.read_json("governor.json")

    # (b) missing or stale governor is a dead governor
    if gov is None:
        return GateDecision(False, "STALE_OR_MISSING_GOVERNOR", ceiling=0)
    set_at = _parse_utc(gov.get("set_at"))
    if set_at is None or (now - set_at).total_seconds() > GOVERNOR_STALE_SECONDS:
        return GateDecision(False, "STALE_OR_MISSING_GOVERNOR", ceiling=0)

    ceiling = int(gov.get("ceiling", 0) or 0)

    # (c) ceiling == 0: halted, no retry
    if ceiling == 0:
        return GateDecision(False, "CEILING_ZERO", ceiling=0)

    # (d) long timeouts are refused while over pace
    if timeout_s > LONG_TIMEOUT_S:
        pct_spent = max(
            float(gov.get("pct_spent_metered", 0.0) or 0.0),
            float(gov.get("pct_spent_projected", 0.0) or 0.0),
        )
        pct_elapsed = float(gov.get("pct_elapsed", 0.0) or 0.0)
        if pct_spent > pct_elapsed + PACE_OVERSHOOT_PP:
            return GateDecision(False, "LONG_TIMEOUT_OVER_PACE", ceiling=ceiling)

    # (d2) generation may not outrun scoring
    if stage in {"generate", "gen_screen"} and target:
        agg = _STATE.read_json("ledger_agg.json") or {}
        cell = agg.get(target) or {}
        gen = int(cell.get("gen", 0) or 0)
        scr = int(cell.get("scr", 0) or 0)
        if gen >= GEN_BACKLOG_MIN_GEN and gen > 2 * max(scr, GEN_BACKLOG_SCR_FLOOR):
            return GateDecision(False, "GEN_BACKLOG_CAP", ceiling=ceiling)

    # validation gate: production scoring on a target is blocked until PASS
    if target and stage in {"screen", "intermediate", "final", "gen_screen"}:
        gate = _STATE.read_json(f"gates/{target}.json")
        if not gate or str(gate.get("status", "")).upper() != "PASS":
            return GateDecision(False, "VALIDATION_GATE_NOT_PASS", ceiling=ceiling)

    if _FLEET is None:
        return GateDecision(False, "GATE_ERROR", ceiling=ceiling)

    # (e)/(f) live count with headroom margin and jittered backoff
    margin = max(10, math.floor(0.1 * ceiling))
    allowlist = _cpu_allowlist_size()

    for attempt in range(MAX_BACKOFF_ATTEMPTS):
        live = (
            _FLEET.count_all(app_id, project_tag)
            - _FLEET.count_cpu_tagged(app_id, project_tag)
            - allowlist
        )
        live = max(live, 0)
        if live < ceiling - margin:
            # (g)
            return GateDecision(True, "OK", live=live, ceiling=ceiling, attempts=attempt)
        _SLEEP(_backoff_delay(attempt))

    return GateDecision(
        False, "AT_CEILING", live=live, ceiling=ceiling, attempts=MAX_BACKOFF_ATTEMPTS
    )


def _cpu_allowlist_size() -> int:
    """CPU singletons that recorded themselves to ``/state/cpu_sandboxes.jsonl``.

    Step (e) subtracts this allowlist as well, because ``host.compute`` dispatch
    paths cannot carry the ``kind=cpu`` tag.
    """
    ids: set[str] = set()
    for line in _STATE.read_lines("cpu_sandboxes.jsonl"):
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        sb_id = row.get("sandbox_id")
        if sb_id and not row.get("terminated"):
            ids.add(sb_id)
    return len(ids)
