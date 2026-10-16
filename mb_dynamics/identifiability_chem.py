"""Post-hoc diagnostics (not pre-registered; reported as such).

1. Two-population control: re-evaluate the same route at delta = 0 (homogeneous rate k, same chemistry point).
2. Non-identifiability of held-out predictions: for each chemistry and fold, collect every chemistry
   parameter point reached by any start of any route (R0, R1, R2, 2pop; extra route parameters dropped),
   re-evaluate it as R0 on the training curves (Q re-profiled), keep the points within 5 % of the best
   training cost, and report the range of their held-out NRMSE.
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import fit_cv as fc
import run_chem as rc
from decompose_gain import mean_test, train_cost

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def eval_point(args, route="R0"):
    chem, x, train, test = args
    if not rc.KIN:
        rc.init()
    cr = rc.ChemRoute(chem, route)
    c, Q = train_cost(cr, x, train)
    t = mean_test(fc.metrics(test, rc.predict(cr, x, Q, test))) if test else float("nan")
    return c, t


if __name__ == "__main__":
    res = json.load(open(os.path.join(RES, "chem_results_dt0.5.json")))
    rows = [r for r in res["fitB"] if r["K_p"] == 10.0]
    jobs, keys = [], []
    for r in rows:
        nb = len(rc.BASE[r["chem"]][0])
        for x in r["fit"]["all_x"]:
            jobs.append((r["chem"], list(x[:nb]), [tuple(k) for k in r["train"]], [tuple(k) for k in r["test"]]))
            keys.append((r["chem"], r["fold"], r["route"]))
    rc.init()
    with Pool(4, initializer=rc.init) as p:
        vals = p.map(eval_point, jobs, chunksize=4)
    out = {}
    for (chem, fold, route), (c, t) in zip(keys, vals):
        out.setdefault(f"{chem}|{fold}", []).append(dict(route=route, cost=c, test=t))
    summary = {}
    print("chem fold       n_pts best_train  n_within5%  held-out NRMSE range among near-optimal (R0 at those points)")
    for k, pts in sorted(out.items()):
        best = min(p["cost"] for p in pts)
        near = [p for p in pts if p["cost"] <= 1.05 * best]
        tt = [p["test"] for p in near if np.isfinite(p["test"])]
        summary[k] = dict(best_cost=best, n_points=len(pts), n_near=len(near),
                          test_min=min(tt) if tt else None, test_max=max(tt) if tt else None)
        rng = f"{100*min(tt):.1f}% – {100*max(tt):.1f}%" if tt else "—"
        print(f"{k:16s} {len(pts):3d}   {best:.4f}     {len(near):3d}       {rng}")
    # two-population control at delta = 0
    dec = []
    for r in rows:
        if r["route"] != "2pop":
            continue
        nb = len(rc.BASE[r["chem"]][0])
        # same route, same chemistry point, heterogeneity switched off (delta = 0)
        c0, t0 = eval_point((r["chem"], list(r["fit"]["x"][:nb]) + [0.0], [tuple(k) for k in r["train"]],
                             [tuple(k) for k in r["test"]]), route="2pop")
        dec.append(dict(chem=r["chem"], fold=r["fold"], delta=r["fit"]["x"][nb], fit_cost=r["fit"]["cost"],
                        fit_test=mean_test(r["test_metrics"]), delta0_cost=c0, delta0_test=t0))
    print("\n2pop: chem fold  delta | train cost 2pop / at delta=0 | held-out 2pop / at delta=0")
    for d in sorted(dec, key=lambda d: (d["chem"], d["fold"])):
        print(f"{d['chem']} {d['fold']:10s} {d['delta']:.2f} | {d['fit_cost']:.4f} / {d['delta0_cost']:.4f} | "
              f"{100*d['fit_test']:.1f} / {100*d['delta0_test']:.1f}")
    json.dump(dict(near_optimal=summary, points=out, two_pop_delta0=dec),
              open(os.path.join(RES, "chem_identifiability.json"), "w"), indent=1, default=float)
