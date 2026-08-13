"""Demo 2 — validate the predictions against measured data you held out.

    python demos/demo2_validate.py

Demo 1 produces numbers.  This one asks whether to believe them: it predicts a
held-out cell line from its controls alone, then scores those predictions
against what was actually measured, with the Virtual Cell Challenge's own
metrics and against the baselines that matter.

The comparison that counts is not "beats predicting nothing" — it is:

  * ``control (Δ=0)``            predict no change.  Very hard to beat on error.
  * ``global mean``             the challenge's own reference model.
  * ``naive transfer``          copy the effect measured in other cell lines.
                                This already uses cross-context information and
                                is the honest thing to beat.

Run it on a different ``--target`` to check the conclusion is not one lucky
cell line.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root
import argparse

import numpy as np

from virtualcell import metrics as M
from virtualcell.benchmark import de_truth, evaluate
from virtualcell.data import load_all, shared_perturbations
from virtualcell.model import (ContextTransferModel, ControlBaseline,
                               GlobalMeanBaseline, NaiveTransferBaseline)
from virtualcell.predict import DEFAULT_HYPER


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="K562",
                    choices=["K562", "RPE1", "HepG2", "Jurkat"])
    ap.add_argument("--n-knockdowns", type=int, default=300,
                    help="how many held-out knockdowns to score")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    lines, genes, symbols = load_all()
    target = next(cl for cl in lines if cl.name == args.target)
    sources = [cl for cl in lines if cl.name != args.target]
    common = shared_perturbations(lines)

    rng = np.random.default_rng(args.seed)
    perts = np.sort(rng.choice(common, min(args.n_knockdowns, common.size),
                               replace=False))

    print(f"held out : {target.name}  (model sees only its control profile)")
    print(f"trained on: {[s.name for s in sources]}")
    print(f"scoring  : {perts.size} knockdowns x {genes.size} genes\n")

    truth = de_truth(target, perts)
    n_de = truth.sig.sum(axis=1)
    print(f"measured differential expression: median {np.median(n_de):.0f} genes "
          f"per knockdown, {np.mean(n_de >= 10) * 100:.0f}% have >=10\n")

    models = {
        "control (delta=0)": ControlBaseline(),
        "global mean [challenge baseline]": GlobalMeanBaseline(),
        "naive transfer": NaiveTransferBaseline(),
        "ContextTransfer (ours)": ContextTransferModel(DEFAULT_HYPER),
    }
    scored = {}
    for name, model in models.items():
        model.fit(sources, target.mu, symbols)
        scored[name] = evaluate(model.predict(perts), target, perts, truth,
                                symbols)

    base = scored["global mean [challenge baseline]"]
    head = (f"{'model':<34}{'discrim':>9}{'DE@100':>9}{'MAE':>9}"
            f"{'direction':>11}{'score':>8}{'balanced':>10}")
    print(head); print("-" * len(head))
    for name, m in scored.items():
        print(f"{name:<34}{m['discrimination_score_l1']:>9.3f}"
              f"{m['overlap_at_100']:>9.3f}{m['mae']:>9.4f}"
              f"{m['de_direction_match']:>11.3f}"
              f"{M.vcc_score(m, base)['avg_score']:>8.3f}"
              f"{M.vcc_score(m, base, clip=False)['avg_score']:>+10.3f}")

    ours, naive = scored["ContextTransfer (ours)"], scored["naive transfer"]
    print("\nreading this:")
    print(f"  discrimination 0.500 is chance; both transfer models are well above it")
    thresholds = [("discrimination", "discrimination_score_l1", 1),
                  ("DE overlap", "overlap_at_100", 1), ("MAE", "mae", -1)]
    for label, key, sign in thresholds:
        ok_o = sign * (ours[key] - base[key]) > 0
        ok_n = sign * (naive[key] - base[key]) > 0
        print(f"  {label:<15} vs baseline -- ours: {'better' if ok_o else 'WORSE':<6} "
              f" naive: {'better' if ok_n else 'WORSE'}")
    print("  the challenge enforces minimum thresholds on every metric, so a model")
    print("  failing one is in a different position from one failing none.")


if __name__ == "__main__":
    main()
