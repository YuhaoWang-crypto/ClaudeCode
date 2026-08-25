"""Does a predicted profile's nearest neighbour share the perturbation, or the cell line?

Borrowed from PathoROB (GenBio-PathFM, bioRxiv 2026.03.17.712534), where it
exposed a failure that accuracy was blind to: histopathology models sitting at
~98% balanced accuracy had robustness indices from 0.04 to 0.87, because their
embeddings were organised by *scanner site* rather than by biology.

The same question is the one that matters here, with site replaced by cell line.
Take every (perturbation, cell line) profile, find its nearest neighbour among
the others, and ask what the neighbour shares:

* neighbour shares the **perturbation** -> the representation is organised by
  biology, and cross-context transfer is a matter of decoding it;
* neighbour shares the **cell line** -> the representation is organised by
  context, and no amount of decoder tuning fixes transfer, because the signal
  being transferred is not in there.

Run it on the measured data first.  That number is a property of perturbation
biology, not of any model, and it bounds what a model built on these profiles
can do -- exactly as the replicate ceiling bounds the correlation metrics.

    python -m virtualcell.robustness_index
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .data import load_all, shared_perturbations


def robustness_index(profiles: np.ndarray, perturbation: np.ndarray,
                     context: np.ndarray, metric: str = "correlation"
                     ) -> dict:
    """Fraction of profiles whose nearest neighbour shares the perturbation.

    ``profiles`` is (n, genes); ``perturbation`` and ``context`` are (n,) labels.
    Neighbours within the same context are excluded, because a profile's nearest
    same-line neighbour is trivially same-line and would answer a different
    question.  So this asks: **given that the neighbour is in another cell line,
    is it the same perturbation?**
    """
    from scipy.spatial.distance import cdist

    D = cdist(profiles, profiles, metric=metric)
    np.fill_diagonal(D, np.inf)
    same_context = context[:, None] == context[None, :]
    D = np.where(same_context, np.inf, D)          # cross-context neighbours only

    nn = np.argmin(D, axis=1)
    hit = perturbation[nn] == perturbation
    # chance: one correct partner among the candidates in other contexts
    n_other = (~same_context).sum(axis=1)
    chance = float(np.mean(1.0 / np.maximum(n_other, 1)))
    return {"robustness_index": float(np.mean(hit)), "chance": chance,
            "n_profiles": int(profiles.shape[0]),
            "lift_over_chance": float(np.mean(hit) / max(chance, 1e-12))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-perturbations", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path,
                    default=Path("results/robustness_index.json"))
    args = ap.parse_args()

    lines, _, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(args.seed)
    perts = np.sort(rng.choice(common, min(args.n_perturbations, len(common)),
                               replace=False))

    rows, plabel, clabel = [], [], []
    for cl in lines:
        idx = {n: i for i, n in enumerate(cl.names)}
        for p in perts:
            rows.append(cl.delta[idx[p]])
            plabel.append(p)
            clabel.append(cl.name)
    X = np.stack(rows)
    plabel = np.array(plabel)
    clabel = np.array(clabel)

    out = {"measured effects": robustness_index(X, plabel, clabel)}

    # The same question on raw expression rather than effects, which is what a
    # context encoder sees.  Expected to be far worse: baseline expression is
    # dominated by lineage, which is the point of the comparison.
    basal = []
    for cl in lines:
        for _ in perts:
            basal.append(cl.mu)
    out["control profiles (no perturbation signal)"] = robustness_index(
        np.stack(basal) + 1e-9 * np.random.default_rng(0).normal(
            size=(len(basal), symbols.size)), plabel, clabel)

    print(f"{perts.size} perturbations x {len(lines)} cell lines = "
          f"{X.shape[0]} profiles, {symbols.size:,} genes")
    print(f"\n{'representation':<44}{'RI':>8}{'chance':>9}{'lift':>8}")
    print("-" * 69)
    for name, r in out.items():
        print(f"{name:<44}{r['robustness_index']:>8.3f}{r['chance']:>9.4f}"
              f"{r['lift_over_chance']:>8.1f}x")

    ri = out["measured effects"]["robustness_index"]
    print(f"\nOf {X.shape[0]} measured effect profiles, {ri:.1%} have their "
          f"nearest neighbour\nin another cell line be the same knockdown. "
          f"That is {out['measured effects']['lift_over_chance']:.0f}x chance, "
          f"and it is\nmeasured data -- so it bounds any model built on these "
          f"profiles.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
