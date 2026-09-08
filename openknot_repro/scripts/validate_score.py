"""Recompute the OpenKnot score for every released design and account for every difference.

This is the correctness gate for ``openknot/score.py``.  A faithful
re-implementation should reproduce the released ``target_openknot_score``
column exactly; where it does not, the difference has to be explained rather
than tolerated.  Two confounds are known and are checked for explicitly:

1. Two puzzles (W09 in rounds 1-2, P06 in round 3) were released with two
   slightly different target structures at different times, and the v4.5.0
   release notes say those designs were rescored against the better of the
   two.  So a design is also scored against every same-length target structure
   released for its puzzle, and the best is taken.

2. Reactivities in the release are rounded to three decimals while scoring ran
   on unrounded values.  A residue whose rounded value lands exactly on its
   threshold (0.125 unpaired, 0.5 paired, 0.25 crossed) can therefore fall
   either way.  Each such residue is counted, and the achievable score
   interval is compared against the published value.

Usage:  python scripts/validate_score.py [--nrows N]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.data import load_benchmark  # noqa: E402
from openknot.score import (  # noqa: E402
    CPQ_THRESHOLD,
    ECS_PAIRED_THRESHOLD,
    ECS_UNPAIRED_THRESHOLD,
    crossed_pair_scores,
    crossing_residues,
    embed_target,
    eterna_classic_score,
    filter_singlet_pairs,
)

RESULTS = Path(__file__).resolve().parent.parent / "results"
EXACT = 1e-6


def target_variants(meta: pd.DataFrame) -> dict[str, list[str]]:
    """Distinct target structures released per puzzle."""
    variants: dict[str, list[str]] = defaultdict(list)
    for puzzle, group in meta.groupby("puzzle"):
        for structure in group["target_structure"].dropna().unique():
            variants[puzzle].append(structure)
    return variants


def tie_envelope(structure: str, data, start: int, end: int) -> tuple[int, float]:
    """Count residues sitting exactly on a scoring threshold, and the score span they cover.

    Returns (n_ties, max_score_swing): flipping every tie residue from miss to
    hit moves the OpenKnot score by at most ``max_score_swing``.
    """
    n_region = end - start + 1
    ecs_ties = 0
    for i in range(start, end + 1):
        value = data[i]
        if value is None or not np.isfinite(value):
            continue
        if structure[i] == ".":
            if value == ECS_UNPAIRED_THRESHOLD:
                ecs_ties += 1
        elif value == ECS_PAIRED_THRESHOLD:
            ecs_ties += 1

    crossed = crossing_residues(filter_singlet_pairs(structure))
    crossed_in_region = [i for i in crossed if start <= i <= end]
    cpq_ties = sum(
        1 for i in crossed_in_region if np.isfinite(data[i]) and data[i] == CPQ_THRESHOLD
    )
    swing = 0.5 * 100.0 * ecs_ties / n_region
    if crossed_in_region:
        swing += 0.5 * 100.0 * cpq_ties / len(crossed_in_region)
    return ecs_ties + cpq_ties, swing


def score_all(meta: pd.DataFrame, reactivity: np.ndarray) -> pd.DataFrame:
    variants = target_variants(meta)
    n = len(meta)
    own = np.full(n, np.nan)
    best = np.full(n, np.nan)
    ecs_col = np.full(n, np.nan)
    cpq_col = np.full(n, np.nan)
    ties = np.zeros(n, dtype=int)
    swing = np.zeros(n)
    embed_cache: dict[tuple, str | None] = {}

    for i, row in enumerate(meta.itertuples(index=False)):
        target = row.target_structure
        if not isinstance(target, str):
            continue
        start, end = int(row.sub_start) - 1, int(row.sub_end) - 1
        full_length = len(row.sequence)
        data = reactivity[i][:full_length]  # matrix is padded to the longest construct

        scores = []
        for candidate in variants[row.puzzle]:
            if len(candidate) != end - start + 1:
                continue  # a different design length, not an alternate target
            key = (candidate, start, full_length)
            if key not in embed_cache:
                try:
                    embed_cache[key] = embed_target(candidate, start + 1, end + 1, full_length)
                except ValueError:
                    embed_cache[key] = None
            full = embed_cache[key]
            if full is None:
                continue
            ecs = eterna_classic_score(full, data, start, end)
            cpq = crossed_pair_scores(full, data, start, end)[1]
            scores.append((0.5 * ecs + 0.5 * cpq, ecs, cpq, full))
            if candidate == target:
                own[i] = scores[-1][0]
                ecs_col[i] = ecs
                cpq_col[i] = cpq
        if not scores:
            continue
        top = max(scores, key=lambda s: s[0])
        best[i] = top[0]
        ties[i], swing[i] = tie_envelope(top[3], data, start, end)

    out = meta.copy()
    out["ECS_recomputed"] = ecs_col
    out["CPQ_recomputed"] = cpq_col
    out["OKS_recomputed"] = own
    out["OKS_best_variant"] = best
    out["n_threshold_ties"] = ties
    out["tie_swing"] = swing
    out["delta"] = out["OKS_recomputed"] - out["target_openknot_score"]
    out["delta_best"] = out["OKS_best_variant"] - out["target_openknot_score"]
    return out


def report(scored: pd.DataFrame) -> pd.DataFrame:
    valid = scored.dropna(subset=["OKS_recomputed", "target_openknot_score"]).copy()
    total = len(valid)
    exact = valid["delta"].abs() < EXACT
    exact_best = valid["delta_best"].abs() < EXACT
    # Strict comparison counts a residue sitting exactly on its threshold as a
    # miss, so the recomputed score is a lower bound: the published value may be
    # up to `tie_swing` higher, and no lower.
    explained = exact_best | (
        (valid["delta_best"] <= EXACT) & (valid["delta_best"] >= -valid["tie_swing"] - EXACT)
    )
    unexplained = valid[~explained]

    rows = [
        ("designs scored", total, ""),
        ("exact vs own target structure", int(exact.sum()), f"{100 * exact.mean():.3f}%"),
        (
            "exact vs best released target variant",
            int(exact_best.sum()),
            f"{100 * exact_best.mean():.3f}%",
        ),
        (
            "explained (exact, or within rounding ties)",
            int(explained.sum()),
            f"{100 * explained.mean():.3f}%",
        ),
        ("unexplained", int(len(unexplained)), f"{100 * len(unexplained) / total:.3f}%"),
    ]
    summary = pd.DataFrame(rows, columns=["quantity", "count", "fraction"])
    print(summary.to_string(index=False))
    if len(unexplained):
        print("\nunexplained residuals by puzzle:")
        print(unexplained.groupby(["round", "puzzle"])["delta_best"].agg(["size", "mean"]).head(20))
    return summary, unexplained


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nrows", type=int, default=None)
    args = parser.parse_args()

    meta, reactivity = load_benchmark(nrows=args.nrows)
    print(f"loaded {len(meta)} designs, reactivity matrix {reactivity.shape}\n")

    scored = score_all(meta, reactivity)
    summary, unexplained = report(scored)

    RESULTS.mkdir(exist_ok=True)
    keep = [
        "id", "round", "puzzle", "method", "SN_filter", "design_length",
        "target_openknot_score", "ECS_recomputed", "CPQ_recomputed", "OKS_recomputed",
        "OKS_best_variant", "delta", "delta_best", "n_threshold_ties", "tie_swing",
        "RNet_F1", "RNet_F1_crossed_pair",
    ]
    scored[keep].to_csv(RESULTS / "recomputed_scores.csv", index=False)
    summary.to_csv(RESULTS / "score_validation.csv", index=False)
    unexplained[keep].to_csv(RESULTS / "score_unexplained.csv", index=False)
    print(f"\nwrote {RESULTS}/recomputed_scores.csv, score_validation.csv, score_unexplained.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
