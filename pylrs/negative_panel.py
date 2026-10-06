"""Design the measured-negative panel: what to actually order, and what it buys.

The binding constraint on this whole dataset is that **no published campaign
reports a variant shown not to work**. Every label is a survivor. That is why
`baseline.py` Task B has to invent its negatives, and why its AUC of 0.98 means
nothing: the simplest possible model wins it, because the task is really "does
this look like a selected clone" rather than "does this recognise the ncAA".

Measuring negatives is the cheapest experiment that changes what the model can
claim. This module turns that into a concrete plate.

## The design decision that matters

The instinct is to order the top-ranked variants. Don't — not only those. A
panel of predicted positives can measure a hit rate but **cannot calibrate the
model**: with no low-scoring wells you never learn whether a low score means
anything, so you cannot later use the score to *exclude* candidates, which is
most of its value.

This emits a **score-stratified** panel instead: wells spread across the score
range, so the readout is a calibration curve (score vs hit rate) rather than a
single number. It costs the same number of wells.

## Why the strata are what they are

* **top** — doubles as the pick-list, so the panel is not a detour.
* **upper-mid / lower-mid** — where the decision boundary actually sits; these
  wells carry the most information per well about where to set a threshold.
* **bottom** — the model's confident negatives. If these come back positive, the
  score is not usable for exclusion, which is a result worth knowing early.
* **random library** — unselected members, i.e. exactly the decoys `baseline.py`
  currently assumes. Measuring these tests the assumption itself.

Run:  python3 pylrs/negative_panel.py --ncaa pAzF --wells 24
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))

AA20 = "ACDEFGHIKLMNPQRSTVWY"

#: (stratum, fraction of wells, what the stratum is for)
STRATA = [
    ("top", 0.25, "predicted hits; doubles as the pick-list"),
    ("upper-mid", 0.20, "near the decision boundary"),
    ("lower-mid", 0.20, "near the decision boundary, below it"),
    ("bottom", 0.20, "the model's confident negatives"),
    ("random-library", 0.15, "unselected library members; tests the decoy assumption"),
]


def auc_standard_error(auc: float, n_pos: int, n_neg: int) -> float:
    """Hanley-McNeil standard error of an AUC estimate."""
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
           + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)
    return math.sqrt(max(var, 0.0))


def random_library(wt: str, positions: list[int], n: int, rng, exclude: set) -> list[str]:
    out = []
    while len(out) < n:
        muts = [f"{wt[p - 1]}{p}{AA20[rng.integers(20)]}" for p in positions]
        muts = [m for m in muts if m[0] != m[-1]]
        key = "/".join(muts)
        if muts and key not in exclude:
            exclude.add(key)
            out.append(key)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/pazf_recombination_ranked.csv"))
    ap.add_argument("--campaigns", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/literature_campaigns.csv"))
    ap.add_argument("--variants", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/tyrrs_variants.csv"))
    ap.add_argument("--ncaa", default="pAzF")
    ap.add_argument("--wells", type=int, default=24)
    ap.add_argument("--positions", default="32,107,158,159,162")
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("pylrs/data"))
    args = ap.parse_args(argv)

    from tyrrs_dataset import wild_type

    rng = np.random.default_rng(0)
    wt = wild_type(args.campaigns)
    positions = [int(p) for p in args.positions.split(",")]
    ranked = pd.read_csv(args.ranked).sort_values("P_pAzF", ascending=False,
                                                  ignore_index=True)
    known = set(pd.read_csv(args.variants).query("label == 1")["mutations"])

    rows = []
    n_ranked = len(ranked)
    for name, fraction, purpose in STRATA:
        k = max(1, round(args.wells * fraction))
        if name == "random-library":
            picks = random_library(wt, positions, k, rng, set(known))
            for m in picks:
                rows.append({"stratum": name, "mutations": m, "score": np.nan,
                             "purpose": purpose})
            continue
        lo, hi = {"top": (0.0, 0.02), "upper-mid": (0.25, 0.35),
                  "lower-mid": (0.55, 0.65), "bottom": (0.95, 1.0)}[name]
        band = ranked.iloc[int(lo * n_ranked):max(int(hi * n_ranked), int(lo * n_ranked) + k)]
        band = band[~band["mutations"].isin(known)]
        take = band.head(k) if name == "top" else band.sample(
            min(k, len(band)), random_state=0)
        for row in take.itertuples():
            rows.append({"stratum": name, "mutations": row.mutations,
                         "score": row.P_pAzF, "purpose": purpose})

    panel = pd.DataFrame(rows)
    panel.insert(0, "well", [f"{chr(65 + i // 12)}{i % 12 + 1}" for i in range(len(panel))])
    panel.insert(1, "ncAA", args.ncaa)
    panel["n_mut"] = panel["mutations"].str.count("/") + 1

    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / f"negative_panel_{args.ncaa}.csv"
    panel.to_csv(target, index=False)

    print(f"Panel for {args.ncaa}: {len(panel)} wells -> {target}\n")
    print(panel.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    print("\n" + "=" * 72)
    print("What this buys, in numbers")
    print("=" * 72)
    n_known_pos = int((pd.read_csv(args.variants)
                       .query("label == 1 and ncAA == @args.ncaa")).shape[0])
    print(f"Known {args.ncaa} positives already in hand: {n_known_pos}")
    print("\nAUC precision if the panel returns the stated split "
          "(Hanley-McNeil 95% CI):")
    print(f"  {'measured negatives':>20}  {'AUC 0.70':>16}  {'AUC 0.85':>16}")
    for n_neg in (10, 15, 20, 30, 50):
        cis = []
        for auc in (0.70, 0.85):
            se = auc_standard_error(auc, n_known_pos, n_neg)
            cis.append(f"±{1.96 * se:.3f}")
        print(f"  {n_neg:>20}  {cis[0]:>16}  {cis[1]:>16}")
    print(f"\n  With {n_known_pos} positives, ~20 measured negatives distinguishes a")
    print("  working model (0.85) from a useless one (0.5); it does NOT resolve")
    print("  0.70 from 0.85. Ranking several ncAAs at ~20 each and pooling is the")
    print("  cheaper route to a tight number than deepening one ncAA.")

    print("\n" + "=" * 72)
    print("Order of operations")
    print("=" * 72)
    print("1. Run this panel for pAzF first: it is the ncAA with the most known")
    print("   positives (7), so the AUC it yields is the best-powered one available.")
    print("2. Report BOTH outcomes -- hits in `top` and non-hits in `bottom` are")
    print("   each a result. A `bottom` well that comes back positive is the most")
    print("   informative single well on the plate.")
    print("3. Feed every measured negative back in and re-run baseline.py Task B.")
    print("   That is the point at which its AUC stops being an artifact.")
    print("\n⚠️ The `top` stratum is drawn from a model whose own calibration misses")
    print("   3 of 7 known pAzF clones (see pazf_calibration.csv). Expect the panel")
    print("   to confirm the dominant motif and to say little about the rest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
