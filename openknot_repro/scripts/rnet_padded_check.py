"""Does the flanking pad explain why the simulated score is too generous?

The experimental OpenKnot score is measured on the padded construct and scored
over the design region. The simulated score in `scripts/rnet_predict.py` is
computed on the isolated design sequence, because that is what the release's
own RNet columns were computed on. Those are not the same molecule: 40-90 extra
nucleotides can change how a design folds.

This re-predicts a subsample on the full padded construct, scores the same
design region, and reports whether agreement with the measurement improves.

Usage:  python scripts/rnet_padded_check.py --limit 150
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.data import BENCH_CSV, META_COLUMNS  # noqa: E402
from openknot.rnet import RNet  # noqa: E402
from openknot.score import crossed_pair_scores, embed_target, eterna_classic_score  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"
OUT = RESULTS / "rnet_padded_check.csv"


def spearman(a, b) -> float:
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=150)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    predictions = pd.read_csv(RESULTS / "rnet_predictions.csv")
    meta = pd.read_csv(BENCH_CSV, usecols=META_COLUMNS, low_memory=False)
    merged = predictions.merge(
        meta[["design_sequence", "sequence", "sub_start", "sub_end", "target_structure"]]
        .drop_duplicates(subset=["design_sequence"]),
        on="design_sequence",
        how="inner",
    )
    merged = merged[merged["sequence"].str.fullmatch(r"[ACGU]+")]

    done: set[str] = set()
    if args.out.exists():
        done = set(pd.read_csv(args.out)["design_sequence"])
    todo = merged[~merged["design_sequence"].isin(done)].head(args.limit - len(done))
    print(f"{len(todo)} constructs to predict (padded), {len(done)} already done")

    if len(todo):
        rnet = RNet(load_ss=False)
        header = not args.out.exists()
        start = time.time()
        for begin in range(0, len(todo), 2):
            batch = todo.iloc[begin : begin + 2]
            reactivities = rnet.reactivity(list(batch["sequence"]))
            rows = []
            for (_, row), reactivity in zip(batch.iterrows(), reactivities):
                start_idx, end_idx = int(row["sub_start"]) - 1, int(row["sub_end"]) - 1
                full = embed_target(
                    row["target_structure"], start_idx + 1, end_idx + 1, len(row["sequence"])
                )
                data = list(reactivity[:, 0])
                ecs = eterna_classic_score(full, data, start_idx, end_idx)
                cpq = crossed_pair_scores(full, data, start_idx, end_idx)[1]
                rows.append({
                    "design_sequence": row["design_sequence"],
                    "puzzle": row["puzzle"],
                    "round": row["round"],
                    "target_openknot_score": row["target_openknot_score"],
                    "simulated_unpadded": row["simulated_openknot_score"],
                    "simulated_padded": 0.5 * ecs + 0.5 * cpq,
                })
            pd.DataFrame(rows).to_csv(args.out, mode="a", header=header, index=False)
            header = False
            elapsed = time.time() - start
            print(f"  {begin + len(batch)}/{len(todo)} ({elapsed / (begin + len(batch)):.1f} s each)",
                  flush=True)

    table = pd.read_csv(args.out)
    measured = table["target_openknot_score"]
    print(f"\n{len(table)} designs scored both ways")
    for label in ("simulated_unpadded", "simulated_padded"):
        values = table[label]
        error = values - measured
        print(f"  {label:20s} Spearman {spearman(values, measured):+.3f}  "
              f"Pearson {np.corrcoef(values, measured)[0, 1]:+.3f}  "
              f"mean signed error {error.mean():+6.2f}  mean |error| {error.abs().mean():5.2f}")
    poor = table[measured < 70]
    if len(poor) > 3:
        print(f"\n  on the {len(poor)} designs that measured below 70:")
        for label in ("simulated_unpadded", "simulated_padded"):
            error = poor[label] - poor["target_openknot_score"]
            print(f"    {label:20s} mean signed error {error.mean():+6.2f}")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
