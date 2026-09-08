"""Design sequences for the OpenKnot targets with gRNAde, and score them with RNet.

This re-runs the design half of the paper for one of its three AI methods, on
the same Round 3 and Round 4 targets, and scores the results the only way that
is possible without a wet lab: an OpenKnot score computed from RNet-predicted
reactivity instead of a measurement.

Sampling budget is the honest caveat. The paper's gRNAde submissions were
selected from up to a million samples per target; this runs a few per target on
CPU. To keep the comparison interpretable, `scripts/compare_designs.py` scores
the *released* gRNAde designs the same way, which separates "smaller search"
from "different model or scorer".

Usage:
  python scripts/grnade_design.py --rounds 3 --samples 8
  python scripts/grnade_design.py --rounds 3 4 --samples 8 --mode 3d
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

from openknot.grnade import GRNAde, load_targets  # noqa: E402
from openknot.rnet import RNet  # noqa: E402
from openknot.score import (  # noqa: E402
    crossed_pair_f1,
    crossed_pair_scores,
    eterna_classic_score,
    pair_f1,
)

RESULTS = Path(__file__).resolve().parent.parent / "results"
OUT = RESULTS / "grnade_designs.csv"
TEMPERATURES = (0.1, 0.5)


def simulated_score(target: str, reactivity: np.ndarray) -> float:
    data = list(reactivity)
    ecs = eterna_classic_score(target, data, 0, len(target) - 1)
    cpq = crossed_pair_scores(target, data, 0, len(target) - 1)[1]
    return 0.5 * ecs + 0.5 * cpq


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, nargs="+", default=[3])
    parser.add_argument("--samples", type=int, default=8, help="samples per temperature")
    parser.add_argument("--mode", choices=["2d", "3d"], default="2d")
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    targets = [t for t in load_targets() if t["round"] in args.rounds]
    print(f"{len(targets)} targets in round(s) {args.rounds}, mode {args.mode}")

    done: set[str] = set()
    if args.out.exists():
        done = set(pd.read_csv(args.out)["puzzle"])
        print(f"resuming: {len(done)} targets already designed")
    targets = [t for t in targets if t["Puzzle"] not in done]
    if not targets:
        print("nothing to do")
        return 0

    grnade = GRNAde(mode=args.mode)
    rnet = RNet()
    header = not args.out.exists()
    start = time.time()

    for index, target in enumerate(targets, 1):
        sec_struct = target["Dot-bracket"]
        native = target["Sequence"]
        if args.mode == "3d":
            pdb = next(Path(target["structure_dir"]).glob("*.pdb"), None)
            if pdb is None:
                print(f"  {target['Puzzle']}: no PDB, skipping")
                continue
            data, _ = grnade.featurize_pdb(pdb, sec_struct)
        else:
            data = grnade.featurize_2d(native, sec_struct)

        designs = []
        for temperature in TEMPERATURES:
            for sequence in grnade.sample(data, args.samples, temperature):
                designs.append((sequence, temperature))

        sequences = [d[0] for d in designs]
        structures = rnet.structures(sequences)
        reactivities = rnet.reactivity(sequences)

        rows = []
        for (sequence, temperature), structure, reactivity in zip(
            designs, structures, reactivities
        ):
            rows.append({
                "puzzle": target["Puzzle"],
                "round": target["round"],
                "title": target["Title"],
                "mode": args.mode,
                "temperature": temperature,
                "sequence": sequence,
                "length": len(sequence),
                "identity_to_native": float(
                    np.mean([a == b for a, b in zip(sequence, native)])
                ),
                "target_structure": sec_struct,
                "rnet_structure": structure,
                "rnet_F1": pair_f1(structure, sec_struct),
                "rnet_F1_crossed_pair": crossed_pair_f1(structure, sec_struct),
                "simulated_openknot_score": simulated_score(sec_struct, reactivity[:, 0]),
            })
        pd.DataFrame(rows).to_csv(args.out, mode="a", header=header, index=False)
        header = False
        best = max(r["simulated_openknot_score"] for r in rows)
        elapsed = time.time() - start
        print(
            f"  [{index}/{len(targets)}] {target['Puzzle']} ({len(sec_struct)} nt): "
            f"{len(rows)} designs, best simulated score {best:.1f}, "
            f"{elapsed / index / 60:.1f} min/target",
            flush=True,
        )

    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
