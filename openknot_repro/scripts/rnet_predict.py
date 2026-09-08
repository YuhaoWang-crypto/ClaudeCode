"""Run RibonanzaNet over released designs: secondary structure and SHAPE reactivity.

For each design this records

  * the RNet secondary structure and its F1 against the target, to compare with
    the released ``RNet_structure`` / ``RNet_F1`` columns;
  * the RNet-predicted 2A3 reactivity profile, and the OpenKnot score computed
    from it instead of from the measurement -- the paper's Fig. S4 claim that a
    simulated score tracks the experimental one, especially for poor designs.

RNet is run on the design sequence alone, without the experimental pads. The
released ``RNet_structure`` is design-length, and on a 24-design probe,
predicting on the padded construct and then slicing agreed with the release
distinctly worse (0.90 mean pair-F1) than predicting on the design alone
(0.93 in that probe, 0.999 across the full run).

Writes results incrementally so a long CPU run can be resumed.

Usage:
  python scripts/rnet_predict.py --per-round 100          # stratified sample
  python scripts/rnet_predict.py --ids-from FILE          # specific designs
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.data import BENCH_CSV, META_COLUMNS  # noqa: E402
from openknot.rnet import RNet  # noqa: E402
from openknot.score import (  # noqa: E402
    crossed_pair_f1,
    crossed_pair_scores,
    eterna_classic_score,
    pair_f1,
)

RESULTS = Path(__file__).resolve().parent.parent / "results"
OUT = RESULTS / "rnet_predictions.csv"



def simulated_openknot_score(target: str, reactivity: np.ndarray) -> float:
    """OpenKnot score computed from predicted rather than measured reactivity."""
    data = list(reactivity)
    ecs = eterna_classic_score(target, data, 0, len(target) - 1)
    cpq = crossed_pair_scores(target, data, 0, len(target) - 1)[1]
    return 0.5 * ecs + 0.5 * cpq


def select(meta: pd.DataFrame, per_round: int, seed: int) -> pd.DataFrame:
    usable = meta.dropna(subset=["design_sequence", "target_structure", "RNet_structure"])
    usable = usable[usable["design_sequence"].str.fullmatch(r"[ACGU]+")]
    usable = usable[usable["design_sequence"].str.len() == usable["target_structure"].str.len()]
    parts = [
        group.sample(min(per_round, len(group)), random_state=seed)
        for _, group in usable.groupby("round")
    ]
    return pd.concat(parts).sort_values(["round", "puzzle"]).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-round", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    meta = pd.read_csv(BENCH_CSV, usecols=META_COLUMNS, low_memory=False)
    chosen = select(meta, args.per_round, args.seed)

    done: set[str] = set()
    if args.out.exists():
        previous = pd.read_csv(args.out)
        done = set(previous["design_sequence"])
        print(f"resuming: {len(done)} designs already predicted")
    todo = chosen[~chosen["design_sequence"].isin(done)]
    print(f"{len(todo)} designs to predict (of {len(chosen)} selected)")
    if todo.empty:
        return 0

    rnet = RNet()
    args.out.parent.mkdir(exist_ok=True)
    header = not args.out.exists()
    # group by length so padding never dominates a batch
    todo = todo.assign(_len=todo["design_sequence"].str.len()).sort_values("_len")
    start = time.time()
    processed = 0

    for begin in range(0, len(todo), args.batch_size):
        batch = todo.iloc[begin : begin + args.batch_size]
        sequences = list(batch["design_sequence"])
        structures = rnet.structures(sequences)
        reactivities = rnet.reactivity(sequences)

        rows = []
        for (_, row), structure, reactivity in zip(batch.iterrows(), structures, reactivities):
            target = row["target_structure"]
            rows.append(
                {
                    "id": row["id"],
                    "round": row["round"],
                    "puzzle": row["puzzle"],
                    "method": row["method"],
                    "design_sequence": row["design_sequence"],
                    "design_length": len(row["design_sequence"]),
                    "SN_filter": row["SN_filter"],
                    "target_openknot_score": row["target_openknot_score"],
                    "released_RNet_structure": row["RNet_structure"],
                    "released_RNet_F1": row["RNet_F1"],
                    "released_RNet_F1_crossed_pair": row["RNet_F1_crossed_pair"],
                    "rnet_structure": structure,
                    "rnet_F1": pair_f1(structure, target),
                    "rnet_F1_crossed_pair": crossed_pair_f1(structure, target),
                    "agreement_with_release": pair_f1(structure, row["RNet_structure"]),
                    "simulated_openknot_score": simulated_openknot_score(
                        target, reactivity[:, 0]
                    ),
                    "mean_predicted_2A3": float(reactivity[:, 0].mean()),
                }
            )
        pd.DataFrame(rows).to_csv(args.out, mode="a", header=header, index=False)
        header = False
        processed += len(rows)
        rate = (time.time() - start) / processed
        print(
            f"  {processed}/{len(todo)} designs  ({rate:.1f} s/design, "
            f"~{rate * (len(todo) - processed) / 60:.0f} min left)",
            flush=True,
        )

    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
