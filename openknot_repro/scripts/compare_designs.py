"""Compare the designs generated here with the ones the paper submitted.

The designs made by `scripts/grnade_design.py` can only be scored in silico, so
comparing their simulated score against the paper's *experimental* score would
confound two different things: a much smaller sampling budget, and the gap
between a predicted and a measured SHAPE profile.

This script separates them by scoring the released designs the same way. For
each Round 3 / Round 4 target it takes the best released design of each method
(by experimental OpenKnot score), computes its RNet-simulated score, and puts
that beside the best simulated score reached here. The simulated-vs-simulated
column is the like-for-like one; the experimental column says what the measured
outcome was.

Usage:  python scripts/compare_designs.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.data import BENCH_CSV, META_COLUMNS  # noqa: E402
from openknot.rnet import RNet  # noqa: E402
from openknot.score import crossed_pair_scores, eterna_classic_score  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"
DESIGNS = RESULTS / "grnade_designs.csv"
CACHE = RESULTS / "rnet_predictions.csv"
OUT = RESULTS / "design_comparison.csv"

COMPARED_METHODS = ["gRNAde", "gRNAde-no3d", "Eterna", "Starting sequence"]


def simulated_score(target: str, reactivity: np.ndarray) -> float:
    data = list(reactivity)
    ecs = eterna_classic_score(target, data, 0, len(target) - 1)
    cpq = crossed_pair_scores(target, data, 0, len(target) - 1)[1]
    return 0.5 * ecs + 0.5 * cpq


def best_released_designs(meta: pd.DataFrame, puzzles: set[str]) -> pd.DataFrame:
    """Best released design per (puzzle, method), by experimental OpenKnot score."""
    usable = meta[
        meta["puzzle"].isin(puzzles)
        & meta["method"].isin(COMPARED_METHODS)
        & (meta["SN_filter"] == 1)
    ].dropna(subset=["design_sequence", "target_structure"])
    usable = usable[usable["design_sequence"].str.fullmatch(r"[ACGU]+")]
    usable = usable[usable["design_sequence"].str.len() == usable["target_structure"].str.len()]
    index = usable.groupby(["puzzle", "method"])["target_openknot_score"].idxmax()
    return usable.loc[index].reset_index(drop=True)


def main() -> int:
    if not DESIGNS.exists():
        print("run scripts/grnade_design.py first", file=sys.stderr)
        return 1
    mine = pd.read_csv(DESIGNS)
    puzzles = set(mine["puzzle"])
    print(f"{len(mine)} designs generated here across {len(puzzles)} targets")

    meta = pd.read_csv(BENCH_CSV, usecols=META_COLUMNS, low_memory=False)
    released = best_released_designs(meta, puzzles)
    print(f"{len(released)} released designs to score "
          f"({released['method'].value_counts().to_dict()})")

    cache = {}
    if CACHE.exists():
        previous = pd.read_csv(CACHE)
        cache = dict(zip(previous["design_sequence"], previous["simulated_openknot_score"]))
        print(f"{len(cache)} simulated scores available from the RNet run")

    todo = released[~released["design_sequence"].isin(cache)]
    if len(todo):
        print(f"scoring {len(todo)} released designs with RNet...")
        rnet = RNet(load_ss=False)
        for begin in range(0, len(todo), 4):
            batch = todo.iloc[begin : begin + 4]
            reactivities = rnet.reactivity(list(batch["design_sequence"]))
            for (_, row), reactivity in zip(batch.iterrows(), reactivities):
                cache[row["design_sequence"]] = simulated_score(
                    row["target_structure"], reactivity[:, 0]
                )
            print(f"  {min(begin + 4, len(todo))}/{len(todo)}", flush=True)

    released["simulated_openknot_score"] = released["design_sequence"].map(cache)

    rows = []
    for puzzle in sorted(puzzles):
        ours = mine[mine["puzzle"] == puzzle]
        theirs = released[released["puzzle"] == puzzle]
        row = {
            "puzzle": puzzle,
            "round": int(ours["round"].iloc[0]),
            "length": int(ours["length"].iloc[0]),
            "n_designs_here": len(ours),
            "our_gRNAde_best_simulated": ours["simulated_openknot_score"].max(),
            "our_gRNAde_best_rnet_F1": ours["rnet_F1"].max(),
        }
        for method in COMPARED_METHODS:
            block = theirs[theirs["method"] == method]
            key = method.replace(" ", "_").replace("-", "_")
            row[f"released_{key}_experimental"] = (
                block["target_openknot_score"].iloc[0] if len(block) else np.nan
            )
            row[f"released_{key}_simulated"] = (
                block["simulated_openknot_score"].iloc[0] if len(block) else np.nan
            )
        rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)

    print("\nper-target comparison (simulated = OpenKnot score from RNet reactivity):")
    columns = [
        "puzzle", "length", "our_gRNAde_best_simulated",
        "released_gRNAde_simulated", "released_gRNAde_experimental",
        "released_Eterna_experimental",
    ]
    print(table[columns].round(1).to_string(index=False))

    ours = table["our_gRNAde_best_simulated"]
    theirs = table["released_gRNAde_simulated"]
    both = table.dropna(subset=["our_gRNAde_best_simulated", "released_gRNAde_simulated"])
    print("\nlike-for-like (simulated score, our small sample vs their submitted design):")
    print(f"  ours  mean {ours.mean():.1f}   >90 on {int((ours > 90).sum())}/{len(table)} targets")
    print(f"  theirs mean {theirs.mean():.1f}   >90 on {int((theirs > 90).sum())}/{len(table)} targets")
    print(f"  ours >= theirs on {int((both['our_gRNAde_best_simulated'] >= both['released_gRNAde_simulated']).sum())}"
          f"/{len(both)} targets")
    experimental = table["released_gRNAde_experimental"]
    print(f"\n  for reference, their submitted designs measured >90 on "
          f"{int((experimental > 90).sum())}/{experimental.notna().sum()} targets")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
