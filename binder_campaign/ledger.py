"""The campaign ledgers.

Three append-only JSONL trees on the shared Modal volume, all written as a
**per-writer subfile tree** (``/ledger/{writer_frame_id}.jsonl``, one file per
writing frame, never a single shared file — Modal volume appends across
concurrent sandboxes are not atomic and silently clobber):

* the **design-count ledger** — one row per job completion,
  ``{ts_utc, job_id, target, structure_method, stage, n_generated, n_scored,
  gpu_seconds, writer_frame_id}``, idempotent on ``(job_id, stage)``;
* the **job-metadata ledger** (``/ledger/job_metadata/{launching_frame_id}.jsonl``)
  — one append per dispatch, written immediately after ``submit_gate()`` returns
  True on EVERY dispatch path.  BUDGET computes ``ratecard_usd`` from this
  ledger, never from provider-side usage summaries;
* the **deviations ledger** (``/ledger/deviations.jsonl``) —
  ``{at_utc, kind, target, what, action}``, one JSON object per line,
  append-only, never overwritten in place.

Campaign-level totals are *exactly* aggregations over these files, never counts
of output files:

``designs_generated`` = sum of ``n_generated``;
``designs_screened`` = sum of ``n_scored`` where ``stage == 'screen'``;
``designs_ranked``  = sum of ``n_scored`` where ``stage == 'final'``;
``intermediate`` rows count toward neither.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from glob import glob
from typing import Iterable, Iterator, Literal

__all__ = [
    "STAGES",
    "METHOD_FLOOR",
    "LedgerRow",
    "JobMetadataRow",
    "Deviation",
    "DesignCountLedger",
    "LedgerTotals",
]

Stage = Literal["generate", "gen_screen", "screen", "intermediate", "final"]
STAGES: tuple[str, ...] = ("generate", "gen_screen", "screen", "intermediate", "final")

#: Every asterisked structure-design method must contribute at least this many
#: backbones to every target's scored pool.
METHOD_FLOOR = 50


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class LedgerRow:
    """One design-count ledger row (the minimal schema, fixed at kickoff)."""

    job_id: str
    target: str
    structure_method: str
    stage: str
    n_generated: int
    n_scored: int
    gpu_seconds: float
    writer_frame_id: str
    ts_utc: str = field(default_factory=_utcnow)

    def __post_init__(self) -> None:
        if self.stage not in STAGES:
            raise ValueError(f"unknown stage {self.stage!r}; expected one of {STAGES}")
        if self.n_generated < 0 or self.n_scored < 0:
            raise ValueError("counts must be non-negative")

    @property
    def key(self) -> tuple[str, str]:
        """Writes are idempotent on ``(job_id, stage)``."""
        return (self.job_id, self.stage)


@dataclass(frozen=True)
class JobMetadataRow:
    job_id: str
    campaign: str
    target: str
    launching_frame_id: str
    gpu_type: str
    est_hourly_usd: float
    submitted_at_utc: str = field(default_factory=_utcnow)
    terminal_at_utc: str | None = None
    status: str = "running"


@dataclass(frozen=True)
class Deviation:
    kind: str
    target: str
    what: str
    action: str
    at_utc: str = field(default_factory=_utcnow)


@dataclass(frozen=True)
class LedgerTotals:
    designs_generated: int
    designs_screened: int
    designs_ranked: int
    gpu_seconds: float
    per_target: dict[str, dict[str, int]]
    per_target_method: dict[tuple[str, str], int]

    def floor_matrix(self, targets: Iterable[str], methods: Iterable[str]) -> dict:
        """The ``(target, method)`` structure-count matrix, recomputed every cycle.

        Every cell below :data:`METHOD_FLOOR` is an open obligation the DISPATCH
        sub-agent must fill or log as a deviation; a missing matrix is itself a
        deviation.
        """
        matrix: dict[str, dict[str, int]] = {}
        open_cells: list[dict] = []
        for t in targets:
            row = {}
            for m in methods:
                n = self.per_target_method.get((t, m), 0)
                row[m] = n
                if n < METHOD_FLOOR:
                    open_cells.append(
                        {"target": t, "structure_method": m, "n": n,
                         "shortfall": METHOD_FLOOR - n}
                    )
            matrix[t] = row
        return {
            "set_at": _utcnow(),
            "floor": METHOD_FLOOR,
            "matrix": matrix,
            "open_obligations": open_cells,
        }

    def ledger_agg(self) -> dict[str, dict[str, int]]:
        """``/state/ledger_agg.json``: ``{target: {gen, scr}}`` for submit_gate (d2)."""
        return {
            t: {"gen": v["designs_generated"], "scr": v["designs_screened"]}
            for t, v in self.per_target.items()
        }


class DesignCountLedger:
    """Reader/writer for the per-writer design-count subfile tree."""

    def __init__(self, root: str, writer_frame_id: str | None = None):
        self.root = root
        self.writer_frame_id = writer_frame_id
        os.makedirs(root, exist_ok=True)
        os.makedirs(os.path.join(root, "job_metadata"), exist_ok=True)
        self._seen: set[tuple[str, str]] = set()

    # ---- writing ----------------------------------------------------------- #
    @property
    def _own_path(self) -> str:
        if not self.writer_frame_id:
            raise RuntimeError("this ledger handle is read-only (no writer_frame_id)")
        return os.path.join(self.root, f"{self.writer_frame_id}.jsonl")

    def append(self, row: LedgerRow) -> bool:
        """Append one row.  Returns False if ``(job_id, stage)`` was already written.

        Exactly one writer per subfile, so idempotency only has to be enforced
        within this handle plus a re-scan of our own file on first use.
        """
        if not self._seen:
            self._seen = {
                (r["job_id"], r["stage"])
                for r in _iter_jsonl(self._own_path)
                if "job_id" in r and "stage" in r
            }
        if row.key in self._seen:
            return False
        with open(self._own_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(asdict(row), sort_keys=True) + "\n")
        self._seen.add(row.key)
        return True

    def append_job_metadata(self, row: JobMetadataRow) -> None:
        path = os.path.join(
            self.root, "job_metadata", f"{row.launching_frame_id}.jsonl"
        )
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    def append_deviation(self, dev: Deviation) -> None:
        path = os.path.join(self.root, "deviations.jsonl")
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(asdict(dev), sort_keys=True) + "\n")

    # ---- reading ----------------------------------------------------------- #
    def iter_rows(self) -> Iterator[dict]:
        """Every design-count row across the whole subfile tree, deduplicated.

        Rows are idempotent on ``(job_id, stage)``; a duplicate that slipped in
        (e.g. a retried writer frame) is counted once, so re-runs cannot
        double-count.
        """
        seen: set[tuple[str, str]] = set()
        for path in sorted(glob(os.path.join(self.root, "*.jsonl"))):
            if os.path.basename(path) == "deviations.jsonl":
                continue
            for row in _iter_jsonl(path):
                key = (row.get("job_id"), row.get("stage"))
                if key in seen:
                    continue
                seen.add(key)
                yield row

    def totals(self) -> LedgerTotals:
        """The single aggregation that every reported count comes from."""
        gen = scr = ranked = 0
        gpu_s = 0.0
        per_target: dict[str, dict[str, int]] = defaultdict(
            lambda: {"designs_generated": 0, "designs_screened": 0,
                     "designs_ranked": 0}
        )
        per_tm: dict[tuple[str, str], int] = defaultdict(int)

        for row in self.iter_rows():
            t = row.get("target", "")
            stage = row.get("stage", "")
            n_gen = int(row.get("n_generated", 0) or 0)
            n_scr = int(row.get("n_scored", 0) or 0)
            gpu_s += float(row.get("gpu_seconds", 0.0) or 0.0)

            gen += n_gen
            per_target[t]["designs_generated"] += n_gen
            if n_gen:
                per_tm[(t, row.get("structure_method", ""))] += n_gen

            if stage == "screen":
                scr += n_scr
                per_target[t]["designs_screened"] += n_scr
            elif stage == "final":
                ranked += n_scr
                per_target[t]["designs_ranked"] += n_scr
            elif stage == "gen_screen":
                # a coupled gen_screen job writes its n_generated row and its
                # n_scored row; the screen half counts as screened.
                scr += n_scr
                per_target[t]["designs_screened"] += n_scr

        return LedgerTotals(
            designs_generated=gen,
            designs_screened=scr,
            designs_ranked=ranked,
            gpu_seconds=gpu_s,
            per_target=dict(per_target),
            per_target_method=dict(per_tm),
        )

    def deviations(self) -> list[dict]:
        """The deduplicated union of deviation rows, one object per line."""
        rows: list[dict] = []
        seen: set[str] = set()
        for path in sorted(glob(os.path.join(self.root, "**", "deviations*.jsonl"),
                                recursive=True)):
            for row in _iter_jsonl(path):
                key = json.dumps(row, sort_keys=True)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(row)
        return rows

    def ratecard_usd(self, now: datetime | None = None) -> float:
        """Rate card x elapsed lifetime per job, from the job-metadata ledger only.

        Terminal jobs are billed at their final lifetime.  This is BUDGET's
        cross-check figure; it is never posted as canonical.
        """
        now = now or datetime.now(timezone.utc)
        total = 0.0
        for path in sorted(
            glob(os.path.join(self.root, "job_metadata", "*.jsonl"))
        ):
            for row in _iter_jsonl(path):
                start = _parse(row.get("submitted_at_utc"))
                if start is None:
                    continue
                end = _parse(row.get("terminal_at_utc")) or now
                hours = max((end - start).total_seconds() / 3600.0, 0.0)
                total += hours * float(row.get("est_hourly_usd", 0.0) or 0.0)
        return total


def _iter_jsonl(path: str) -> Iterator[dict]:
    if not os.path.exists(path):
        return
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def _parse(value) -> datetime | None:
    if not value:
        return None
    try:
        dt = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return None
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
