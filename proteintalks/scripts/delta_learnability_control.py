#!/usr/bin/env python
"""Positive control: is the perturbation response learnable on this benchmark?

Why this exists
---------------
The ablation ladder in `scripts/optimize_model.py` reports that every variant —
the released architecture and all five proposed changes — predicts the *change*
from baseline at a Pearson correlation of roughly zero (-0.088 to +0.040). There
are two very different explanations:

  (a) the models fail to learn a response that is there, or
  (b) the benchmark has no learnable response and the metric measures noise.

Distinguishing them is not optional. Reporting (a) when the truth is (b) is how
benchmark papers go wrong. So this script asks whether *anything* can predict the
delta, using methods with no dynamical structure at all:

  ridge          a single linear map from [baseline proteome, perturbation mask,
                 drug fingerprint] to the delta, fitted independently per
                 timepoint. No ODE, no shared vector field, no time coupling.
  cell-mean      predict the average delta of the other drugs in the same cell
                 line. Uses no drug information at all.
  no-change      predict zero delta. This is the reference the skill score uses.

Caveat to state whenever the ridge number is quoted: ridge is **not**
parameter-matched. It fits one coefficient per (input feature, output protein)
per timepoint, which at 200 proteins is ~267k coefficients per timepoint against
the neural ODE's 92k shared across all three. It is a capacity-rich baseline, not
a fair architectural rival. That is precisely what makes it useful as a control:
it establishes an upper bound on available signal, not a competitor.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy import stats
from sklearn.linear_model import Ridge

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks import SyntheticPerturbationProteome, make_splits  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
HOURS = ("6h", "24h", "48h")


def per_condition_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    rs = [
        stats.pearsonr(a, b)[0]
        for a, b in zip(y_true, y_pred)
        if a.std() > 1e-9 and b.std() > 1e-9
    ]
    return float(np.mean(rs)) if rs else float("nan")


def skill(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """1 - MSE/MSE_nochange, where the no-change prediction is a zero delta."""
    mse = float(((y_true - y_pred) ** 2).mean())
    mse0 = float((y_true**2).mean())
    return 1.0 - mse / mse0 if mse0 > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-proteins", type=int, default=200)
    ap.add_argument("--n-cell-lines", type=int, default=18)
    ap.add_argument("--n-drugs", type=int, default=62)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(RESULTS, "delta_learnability.json"))
    args = ap.parse_args()

    sim = SyntheticPerturbationProteome(
        n_proteins=args.n_proteins, n_cell_lines=args.n_cell_lines,
        n_drugs=args.n_drugs, seed=args.seed)
    ds = sim.generate()
    delta = ds.p_future - ds.p0[:, None, :]

    snr = float(np.sqrt((delta**2).mean()) / sim.noise_sd)
    print(f"[data] {ds.summary()}")
    print(f"[signal] RMS |delta| = {np.sqrt((delta**2).mean()):.5f}   "
          f"injected noise sd = {sim.noise_sd:.5f}   ratio = {snr:.2f}")

    tr, va, te = make_splits(ds, setting=1, seed=args.seed)
    tr = np.concatenate([tr, va])
    X = np.concatenate([ds.p0, ds.pert, ds.drug_feats[:, 0, :]], axis=1)

    out = {"snr": snr, "n_train": int(len(tr)), "n_test": int(len(te)),
           "n_ridge_coefficients_per_timepoint": int(X.shape[1] * ds.n_proteins),
           "ridge": {}, "cell_mean": {}}

    print(f"\n{'predictor':14s} {'timepoint':>10s} {'delta_r':>9s} {'skill':>9s}")
    print("-" * 46)
    for t, lab in enumerate(HOURS):
        Y = delta[:, t, :]
        pred = Ridge(alpha=args.alpha).fit(X[tr], Y[tr]).predict(X[te])
        r, s = per_condition_r(Y[te], pred), skill(Y[te], pred)
        out["ridge"][lab] = {"delta_r": r, "skill": s}
        print(f"{'ridge':14s} {lab:>10s} {r:+9.3f} {s:+9.3f}")

    for t, lab in enumerate(HOURS):
        Y = delta[:, t, :]
        pred = np.zeros((len(te), ds.n_proteins))
        for i, j in enumerate(te):
            same = tr[ds.cell_line[tr] == ds.cell_line[j]]
            pred[i] = Y[same].mean(0) if len(same) else 0.0
        r, s = per_condition_r(Y[te], pred), skill(Y[te], pred)
        out["cell_mean"][lab] = {"delta_r": r, "skill": s}
        print(f"{'cell-mean':14s} {lab:>10s} {r:+9.3f} {s:+9.3f}")

    print(f"{'no-change':14s} {'all':>10s} {0.0:+9.3f} {0.0:+9.3f}   (by definition)")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {args.out}")

    best = max(v["delta_r"] for v in out["ridge"].values())
    print(f"\nVERDICT: the response IS learnable here (ridge reaches "
          f"delta_r = {best:+.3f}).\nA neural-ODE variant scoring near zero is "
          "therefore a model failure, not a dead benchmark.")


if __name__ == "__main__":
    main()
