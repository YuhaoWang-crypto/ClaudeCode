"""Committor validation by direct shooting -- the standard, assumption-free test.

The committor learned by Gen-COMPAS predicts, for a configuration x, the
probability of reaching state B before state A.  That prediction is tested
here the only way it can be tested: by launching many fresh trajectories from
the same configuration with independent Maxwell-Boltzmann velocities and
counting how many reach B first.

A correct model gives (a) predicted q on the diagonal against the measured
frequency, and (b) for the predicted transition-state ensemble, a histogram of
measured committors peaked at 1/2.

None of these trajectories is part of the training data.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import mdops  # noqa: E402
from analysis import load_store, fit_committor  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("--n-structures", type=int, default=30)
    ap.add_argument("--n-shots", type=int, default=15)
    ap.add_argument("--ps", type=float, default=20.0)
    ap.add_argument("--rng", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.rng)
    store = load_store(args.tag)
    qnet, _ = fit_committor(store)
    feats = common.featurize(store["coords"])
    q = qnet.predict(feats)

    # stratified selection across the whole committor range, with extra
    # weight on the predicted transition-state ensemble
    bins = [(0.0, 0.1), (0.1, 0.3), (0.3, 0.45), (0.45, 0.55),
            (0.55, 0.7), (0.7, 0.9), (0.9, 1.0)]
    per = max(1, args.n_structures // len(bins))
    picks = []
    for lo, hi in bins:
        cand = np.where((q >= lo) & (q < hi))[0]
        if len(cand) == 0:
            continue
        take = min(per, len(cand))
        picks.extend(rng.choice(cand, take, replace=False).tolist())
    picks = np.array(sorted(set(picks)))

    engine = mdops.Engine(threads=1, seed=args.rng)
    rows = []
    ns = 0.0
    for k, i in enumerate(picks):
        x0 = store["coords"][i]
        outs = []
        for s in range(args.n_shots):
            r = mdops.shoot(engine, x0, total_ps=args.ps, save_ps=1.0,
                            seed=int(rng.integers(1, 2 ** 30)))
            if r is not None and r["outcome"] >= 0:
                outs.append(r["outcome"])
        ns = engine.ns_used()
        if not outs:
            continue
        emp = float(np.mean(outs))
        rows.append(dict(index=int(i), q_pred=float(q[i]), q_emp=emp,
                         n=len(outs),
                         phi=float(common.phi_psi(x0)[0]),
                         psi=float(common.phi_psi(x0)[1])))
        print(f"  [{k + 1}/{len(picks)}] q_pred={q[i]:.3f}  "
              f"q_emp={emp:.3f} ({len(outs)} committed shots)  "
              f"phi={rows[-1]['phi']:6.1f}", flush=True)

    qp = np.array([r["q_pred"] for r in rows])
    qe = np.array([r["q_emp"] for r in rows])
    n = np.array([r["n"] for r in rows])

    mae = float(np.mean(np.abs(qp - qe)))
    corr = float(np.corrcoef(qp, qe)[0, 1]) if len(qp) > 2 else float("nan")
    # expected MAE if the model were perfect, given binomial noise
    exp_mae = float(np.mean(np.sqrt(qp * (1 - qp) / n) * np.sqrt(2 / np.pi)))

    tse = np.abs(qp - 0.5) < 0.05
    print(f"\n[validate] {len(rows)} structures, {args.n_shots} shots each, "
          f"{ns:.2f} ns of fresh MD")
    print(f"[validate] mean |q_pred - q_emp| = {mae:.3f} "
          f"(binomial noise floor {exp_mae:.3f})")
    print(f"[validate] correlation = {corr:.3f}")
    if tse.sum():
        print(f"[validate] predicted TSE (|q-0.5|<0.05): {tse.sum()} structures, "
              f"measured committor {qe[tse].mean():.3f} "
              f"+/- {qe[tse].std():.3f}")

    out = dict(rows=rows, mae=mae, expected_mae=exp_mae, corr=corr,
               ns=float(ns), n_shots=args.n_shots, ps=args.ps,
               tse_mean=float(qe[tse].mean()) if tse.sum() else None,
               tse_std=float(qe[tse].std()) if tse.sum() else None,
               tse_n=int(tse.sum()))
    with open(os.path.join(common.RESULTS,
                           f"validation_{args.tag}.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[validate] wrote results/validation_{args.tag}.json")


if __name__ == "__main__":
    main()
