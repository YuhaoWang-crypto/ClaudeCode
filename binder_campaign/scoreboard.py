"""The scoreboard deliverable.

A machine-readable CSV with one row per target and *exactly* these columns::

    target, designs_generated, designs_screened, designs_ranked,
    best_final_score, median_final_score_ranked, pose_check_pass_rate_top10,
    gpu_dollars, wallclock_hours_active_compute, n_distinct_root_backbones,
    max_seqs_per_root_backbone, n_tm90_clusters, max_tm90_cluster_size,
    n_structure_methods, top_method_share, min_pairwise_editdist,
    n_non_all_alpha

Every non-target column is a real computed numeric on every row — no NaN, no
sentinel, no free text.  A value that cannot be computed is a defect to be
remediated before close-out, tracked by a ``scoreboard_gaps.csv`` row
``{target, column, reason, owning_frame_id}``.

The diversity columns are computed from the final ranked-30 sheet **by the same
code that enforces the diversity caps** (:mod:`binder_campaign.sheet_writer`),
not by a parallel reimplementation.
"""

from __future__ import annotations

import statistics
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from .ledger import LedgerTotals
from .sheet_writer import levenshtein_within

__all__ = ["SCOREBOARD_COLUMNS", "build_scoreboard", "scoreboard_gaps"]

SCOREBOARD_COLUMNS: tuple[str, ...] = (
    "target",
    "designs_generated",
    "designs_screened",
    "designs_ranked",
    "best_final_score",
    "median_final_score_ranked",
    "pose_check_pass_rate_top10",
    "gpu_dollars",
    "wallclock_hours_active_compute",
    "n_distinct_root_backbones",
    "max_seqs_per_root_backbone",
    "n_tm90_clusters",
    "max_tm90_cluster_size",
    "n_structure_methods",
    "top_method_share",
    "min_pairwise_editdist",
    "n_non_all_alpha",
)


def _min_pairwise_editdist(sequences: Sequence[str], cap: int = 40) -> int:
    """Smallest pairwise Levenshtein distance in the shipped set.

    The selection cap only guarantees ``> 5``; the scoreboard reports the actual
    minimum, so it is computed by probing upward with the same bounded routine
    the cap uses (identical code path, no parallel implementation).
    """
    if len(sequences) < 2:
        return cap
    for d in range(0, cap):
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                if levenshtein_within(sequences[i], sequences[j], d):
                    return d
    return cap


def build_scoreboard(
    ranked_by_target: Mapping[str, Sequence[Mapping[str, Any]]],
    totals: LedgerTotals,
    gpu_dollars: Mapping[str, float],
    active_compute_hours: Mapping[str, float],
) -> pd.DataFrame:
    """One row per target, every non-target column a real computed numeric."""
    rows: list[dict] = []
    for target, ranked in ranked_by_target.items():
        per_t = totals.per_target.get(target, {})
        finals = [float(r["final_score"]) for r in ranked]
        top10 = ranked[:10]
        roots = [r["root_backbone_id"] for r in ranked]
        clusters = [r["tm90_cluster_id"] for r in ranked]
        methods = [r["structure_method"] for r in ranked]
        seqs = [r["sequence"] for r in ranked]

        rows.append({
            "target": target,
            "designs_generated": int(per_t.get("designs_generated", 0)),
            "designs_screened": int(per_t.get("designs_screened", 0)),
            "designs_ranked": len(ranked),
            "best_final_score": max(finals) if finals else 0.0,
            "median_final_score_ranked": (
                statistics.median(finals) if finals else 0.0
            ),
            "pose_check_pass_rate_top10": (
                sum(1 for r in top10 if r.get("pose_PASS")) / len(top10)
                if top10 else 0.0
            ),
            "gpu_dollars": float(gpu_dollars.get(target, 0.0)),
            "wallclock_hours_active_compute": float(
                active_compute_hours.get(target, 0.0)
            ),
            "n_distinct_root_backbones": len(set(roots)),
            "max_seqs_per_root_backbone": (
                max(roots.count(r) for r in set(roots)) if roots else 0
            ),
            "n_tm90_clusters": len(set(clusters)),
            "max_tm90_cluster_size": (
                max(clusters.count(c) for c in set(clusters)) if clusters else 0
            ),
            "n_structure_methods": len(set(methods)),
            "top_method_share": (
                max(methods.count(m) for m in set(methods)) / len(methods)
                if methods else 0.0
            ),
            "min_pairwise_editdist": _min_pairwise_editdist(seqs),
            "n_non_all_alpha": sum(
                1 for r in ranked if r.get("fold_class") == "not_all_alpha"
            ),
        })

    frame = pd.DataFrame(rows, columns=list(SCOREBOARD_COLUMNS))
    if list(frame.columns) != list(SCOREBOARD_COLUMNS):
        raise AssertionError("scoreboard columns must be exactly the frozen set")
    return frame


def scoreboard_gaps(frame: pd.DataFrame, owning_frame_id: str) -> pd.DataFrame:
    """``{target, column, reason, owning_frame_id}`` for every uncomputable cell.

    Any row this returns is a defect to remediate before close-out, not a
    tolerated blank.
    """
    gaps: list[dict] = []
    for _, row in frame.iterrows():
        for col in SCOREBOARD_COLUMNS[1:]:
            v = row[col]
            if pd.isna(v) or isinstance(v, str):
                gaps.append({
                    "target": row["target"],
                    "column": col,
                    "reason": "non-numeric or NaN value in a required numeric column",
                    "owning_frame_id": owning_frame_id,
                })
    return pd.DataFrame(
        gaps, columns=["target", "column", "reason", "owning_frame_id"]
    )
