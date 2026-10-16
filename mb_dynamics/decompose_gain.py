"""Decompose any R1/R2-vs-R0 held-out gain into (i) a different chemistry optimum and (ii) conformational dynamics.

For every R1/R2 fit-B result: re-evaluate (a) R0 at the R1/R2 chemistry parameters (dynamics removed, same
chemistry point; Q re-profiled on the same training curves) and (b) the same route at the same parameters with
s = 1.  If (a) reproduces the R1/R2 held-out error, the gain over the R0 fit is an optimizer/basin effect, not dynamics.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

import fit_cv as fc
import run_chem as rc

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def mean_test(m):
    return float(np.mean([v["nrmse_net"] for v in m.values()])) if m else float("nan")


def train_cost(cr, x, train):
    obs = fc.data()
    mod = rc.model_curves(cr, x, train)
    Q = fc.profile_Q(obs, mod, train)
    r = fc.residuals(obs, mod, train, Q)
    return float(np.sum(r ** 2)), Q


def one(r):
    rc.init() if not rc.KIN else None
    train, test = [tuple(k) for k in r["train"]], [tuple(k) for k in r["test"]]
    nb = len(rc.BASE[r["chem"]][0])
    x = r["fit"]["x"]
    out = dict(chem=r["chem"], route=r["route"], fold=r["fold"], log10_s=x[nb], fit_cost=r["fit"]["cost"],
               fit_test=mean_test(r["test_metrics"]))
    cr0 = rc.ChemRoute(r["chem"], "R0")
    c0, Q0 = train_cost(cr0, x[:nb], train)
    out["R0_at_params_cost"] = c0
    out["R0_at_params_test"] = mean_test(fc.metrics(test, rc.predict(cr0, x[:nb], Q0, test))) if test else float("nan")
    crr = rc.ChemRoute(r["chem"], r["route"])
    c1, Q1 = train_cost(crr, list(x[:nb]) + [0.0], train)
    out["route_s1_cost"] = c1
    out["route_s1_test"] = mean_test(fc.metrics(test, rc.predict(crr, list(x[:nb]) + [0.0], Q1, test))) if test else float("nan")
    out["chem_params"] = x[:nb]
    return out


if __name__ == "__main__":
    res = json.load(open(os.path.join(RES, "chem_results_dt0.5.json")))
    rows = [r for r in res["fitB"] if r["route"] in ("R1", "R2") and r["K_p"] == 10.0]
    r0 = {(r["chem"], r["fold"]): r for r in res["fitB"] if r["route"] == "R0" and r["K_p"] == 10.0}
    rc.init()
    with Pool(4, initializer=rc.init) as p:
        dec = p.map(one, rows, chunksize=1)
    for d in dec:
        b = r0[(d["chem"], d["fold"])]
        d["R0_fit_cost"] = b["fit"]["cost"]
        d["R0_fit_test"] = mean_test(b["test_metrics"])
    json.dump(dec, open(os.path.join(RES, "chem_gain_decomposition.json"), "w"), indent=1, default=float)
    print("chem route fold        log10_s | train cost: R0fit  R0@params route@s=1 route | held-out NRMSE: R0fit R0@params route@s=1 route")
    for d in sorted(dec, key=lambda d: (d["chem"], d["route"], d["fold"])):
        print(f"{d['chem']} {d['route']} {d['fold']:10s} {d['log10_s']:+.2f} | {d['R0_fit_cost']:.4f} {d['R0_at_params_cost']:.4f} "
              f"{d['route_s1_cost']:.4f} {d['fit_cost']:.4f} | {100*d['R0_fit_test']:.1f} {100*d['R0_at_params_test']:.1f} "
              f"{100*d['route_s1_test']:.1f} {100*d['fit_test']:.1f}")
