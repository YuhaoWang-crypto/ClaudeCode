"""Per-method design success rates, the quantity the paper's bar charts report.

A method succeeds on a target when at least one of its designs for that target
reaches an OpenKnot score above 90 (the paper's cutoff, set from the ~10% error
rate of SHAPE-directed secondary structure assignment).

Two filters matter and are reported separately:

* ``SN_filter`` -- the release's own signal-to-noise gate (reads > 100 and
  signal/noise > 1).  The dataset README says to ignore designs that fail it.
* ``RNet_F1 >= 0.8`` -- the secondary-structure accuracy filter the paper applies
  before the solid bars of Fig. 2C/2F, to drop designs whose SHAPE profile fits
  the target while their actual base pairing does not.

Rates are computed twice, once from the released ``target_openknot_score`` and
once from the score recomputed by ``openknot/score.py``, so the two can be
compared directly.

Usage:  python scripts/success_rates.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.data import AI_METHODS, BASELINE_METHOD, HUMAN_METHOD  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"
SUCCESS_CUTOFF = 90.0
RNET_F1_CUTOFF = 0.8
# 'M2' rows are mutant libraries built on top of designs, not a design method.
NON_METHODS = {"M2"}


def load_scores() -> pd.DataFrame:
    path = RESULTS / "recomputed_scores.csv"
    if not path.exists():
        raise FileNotFoundError("run scripts/validate_score.py first")
    scores = pd.read_csv(path)
    return scores[~scores["method"].isin(NON_METHODS)].copy()


def per_target_best(scores: pd.DataFrame, score_column: str, rnet_filter: bool) -> pd.DataFrame:
    """Best score each method reached on each target, under the chosen filters."""
    subset = scores[scores["SN_filter"] == 1]
    if rnet_filter:
        subset = subset[subset["RNet_F1"] >= RNET_F1_CUTOFF]
    return (
        subset.groupby(["round", "method", "puzzle"])[score_column]
        .max()
        .rename("best_score")
        .reset_index()
    )


def rates(scores: pd.DataFrame, score_column: str, rnet_filter: bool) -> pd.DataFrame:
    """Success rate per (round, method), over the targets that method submitted to."""
    best = per_target_best(scores, score_column, rnet_filter)
    best["success"] = best["best_score"] > SUCCESS_CUTOFF
    targets_per_round = scores.groupby("round")["puzzle"].nunique().to_dict()

    out = []
    for (rnd, method), group in best.groupby(["round", "method"]):
        n_submitted = len(group)
        n_success = int(group["success"].sum())
        rate = n_success / n_submitted if n_submitted else np.nan
        se = np.sqrt(rate * (1 - rate) / n_submitted) if n_submitted else np.nan
        out.append(
            {
                "round": rnd,
                "method": method,
                "score_column": score_column,
                "rnet_filter": rnet_filter,
                "n_targets_in_round": targets_per_round[rnd],
                "n_targets_submitted": n_submitted,
                "n_success": n_success,
                "success_rate": 100 * rate,
                "standard_error": 100 * se,
            }
        )
    return pd.DataFrame(out)


def grouped_rates(scores: pd.DataFrame, score_column: str, rnet_filter: bool) -> pd.DataFrame:
    """Targets solved by *any* AI method vs by Eterna participants.

    This is the paper's headline comparison ("19/20 targets in Round 3 and in
    Round 4 for both AI methods and Eterna participants").
    """
    best = per_target_best(scores, score_column, rnet_filter)
    groups = {
        "AI (any method)": AI_METHODS,
        "Eterna (human)": [HUMAN_METHOD],
        BASELINE_METHOD: [BASELINE_METHOD],
    }
    out = []
    for label, members in groups.items():
        subset = best[best["method"].isin(members)]
        for rnd, group in subset.groupby("round"):
            solved = group.groupby("puzzle")["best_score"].max() > SUCCESS_CUTOFF
            n = len(solved)
            rate = solved.mean() if n else np.nan
            out.append(
                {
                    "round": rnd,
                    "group": label,
                    "score_column": score_column,
                    "rnet_filter": rnet_filter,
                    "n_targets": n,
                    "n_solved": int(solved.sum()),
                    "success_rate": 100 * rate,
                    "standard_error": 100 * np.sqrt(rate * (1 - rate) / n) if n else np.nan,
                }
            )
    return pd.DataFrame(out)


def main() -> int:
    scores = load_scores()
    per_method = []
    per_group = []
    for score_column in ("target_openknot_score", "OKS_recomputed"):
        for rnet_filter in (False, True):
            per_method.append(rates(scores, score_column, rnet_filter))
            per_group.append(grouped_rates(scores, score_column, rnet_filter))
    per_method = pd.concat(per_method, ignore_index=True)
    per_group = pd.concat(per_group, ignore_index=True)

    per_method.to_csv(RESULTS / "success_rates_by_method.csv", index=False)
    per_group.to_csv(RESULTS / "success_rates_by_group.csv", index=False)

    published = per_method[per_method["score_column"] == "target_openknot_score"]
    recomputed = per_method[per_method["score_column"] == "OKS_recomputed"]
    merged = published.merge(
        recomputed,
        on=["round", "method", "rnet_filter"],
        suffixes=("_published", "_recomputed"),
    )
    gap = (merged["success_rate_published"] - merged["success_rate_recomputed"]).abs()
    print(
        f"per-method success rates: {len(merged)} (round, method, filter) cells; "
        f"max |published - recomputed| = {gap.max():.3f} percentage points"
    )

    headline = per_group[
        (per_group["score_column"] == "target_openknot_score") & per_group["rnet_filter"]
    ]
    print("\ntargets solved (OpenKnot score > 90, RNet F1 >= 0.8, published scores):")
    for _, row in headline.sort_values(["round", "group"]).iterrows():
        print(
            f"  round {row['round']}  {row['group']:<18s} "
            f"{row['n_solved']:2d}/{row['n_targets']:2d}"
        )
    print(f"\nwrote {RESULTS}/success_rates_by_method.csv and success_rates_by_group.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
