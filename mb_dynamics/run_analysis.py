"""Three-route comparison on the 2017 DNA20-MB SWV benchmark.

    python run_analysis.py --dt 0.5            # needs md_runs/dt{dt}_s10[1-4].npz and results/r2_model_dt{dt}.npz
Writes results/cv_results_dt{dt}.json (all fits, predictions, metrics, checks).
"""
import argparse
import json
import os
import time
import zlib
from multiprocessing import Pool

import numpy as np

import fit_cv as fc
import md_analysis as ma
import swv_backend as sb

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
KT = 1.380649e-23 * 298.15
ETA = 0.890e-3
A_MB = 0.35e-9
D_SE_NM2_S = KT / (6 * np.pi * ETA * A_MB) * 1e18      # Stokes-Einstein, MB sphere in water (nm^2/s)

ROUTES = {}


def build_routes(dt, beta):
    names = [f"dt{dt}_s{s}" for s in (101, 102, 103, 104)]
    z = np.concatenate([ma.load(n)["z"].ravel() for n in names])
    w = np.exp(-beta * (z - ma.Z_MIN))
    r2 = np.load(os.path.join(RES, f"r2_model_dt{dt}.npz"))
    kin2 = sb.Kinetics(r2["pooled_pi"], r2[f"pooled_w{beta:g}"], r2["pooled_L"], name="R2")
    out = {
        "R0": fc.Route("R0", sb.fast_eq_kinetics([1.0], [w.mean()])),
        "R1": fc.Route("R1", sb.smoluchowski_kinetics(z, 40, D_SE_NM2_S, beta, ma.Z_MIN), has_s=True),
        "R1_80": fc.Route("R1_80", sb.smoluchowski_kinetics(z, 80, D_SE_NM2_S, beta, ma.Z_MIN), has_s=True),
        "R2": fc.Route("R2", kin2, has_s=True),
        "2pop": fc.Route("2pop", two_pop=True),
    }
    info = dict(w_mean=float(w.mean()), n_frames=int(z.size), D_SE_nm2_per_s=D_SE_NM2_S,
                r2_k=int(r2["k"]), r2_lag_ps=float(r2["lag_ps"]), r2_eligible=bool(r2["eligible"]))
    return out, info


def _init(dt):
    for beta in (2.0, 5.0, 10.0):
        ROUTES[beta] = build_routes(dt, beta)


def task_fitB(args):
    beta, rname, fkind, fname, train, test = args
    t = time.time()
    route = ROUTES[beta][0][rname]
    best = fc.fit_route(route, train, fc.starts_for(route, n=6, seed=zlib.crc32(f'{rname}|{fname}|{beta}'.encode()) % 1000))
    allk = train + test
    pred = fc.predict(route, best["x"], best["Q"], allk)
    res = dict(kind="B", beta=beta, route=rname, fold_kind=fkind, fold=fname, train=train, test=test,
               fit=best, param_names=route.names(), train_metrics=fc.metrics(train, pred),
               test_metrics=fc.metrics(test, pred) if test else {},
               curves={f"f{int(k[0])}_A{int(round(k[1]*1000))}": {c: pred[k][c].tolist() for c in fc.COMPS} for k in allk},
               wall_s=time.time() - t)
    if rname == "R1":       # grid check: same parameters on 80 cells
        r80 = ROUTES[beta][0]["R1_80"]
        p80 = fc.predict(r80, best["x"], best["Q"], allk)
        res["grid80_max_rel_diff"] = float(max(np.abs(p80[k][c] - pred[k][c]).max() / np.ptp(pred[k][c])
                                              for k in allk for c in fc.COMPS))
        res["grid80_test_metrics"] = fc.metrics(test, p80) if test else {}
    return res


def task_fitA(args):
    beta, fkind, fname, train, test = args
    t = time.time()
    routes = [ROUTES[beta][0][n] for n in ("R0", "R1", "R2")]
    starts = fc.starts_for(routes[0], n=6, seed=zlib.crc32(f"A|{fname}".encode()) % 1000)
    best = fc.fit_shared(routes, train, starts)
    allk = train + test
    out = dict(kind="A", beta=beta, fold_kind=fkind, fold=fname, train=train, test=test, fit=best, routes={})
    for r in routes:
        x = best["x"] + ([0.0] if r.has_s else [])
        pred = fc.predict(r, x, best["Q"], allk)
        out["routes"][r.name] = dict(train_metrics=fc.metrics(train, pred), test_metrics=fc.metrics(test, pred) if test else {},
                                     curves={f"f{int(k[0])}_A{int(round(k[1]*1000))}": {c: pred[k][c].tolist() for c in fc.COMPS} for k in allk})
    out["wall_s"] = time.time() - t
    return out


def task_profile(args):
    beta, rname, pname, value = args
    route = ROUTES[beta][0][rname]
    allk = [(f, a) for f in fc.FREQS for a in fc.AMPS]
    names = route.names()
    i = names.index(pname)
    lo, hi = np.array(route.bounds()).T
    sub = fc.Route(rname, route.kin, has_s=route.has_s, two_pop=route.two_pop)
    # fix parameter i by narrowing its bounds
    b = route.bounds()
    b[i] = (value - 1e-6, value + 1e-6)
    sub.bounds = lambda b=b: b
    starts = fc.starts_for(route, n=4, seed=int(1000 * value) % 997)
    for s in starts:
        s[i] = value
    best = fc.fit_route(sub, allk, starts)
    return dict(kind="profile", beta=beta, route=rname, param=pname, value=value, cost=best["cost"], x=best["x"], Q=best["Q"])


def timescale_diagnostics(beta):
    routes, info = ROUTES[beta]
    out = {}
    for rn in ("R1", "R1_80", "R2"):
        kin = routes[rn].kin
        lam = np.sort(np.linalg.eigvalsh(kin.S))[::-1]
        out[rn] = dict(slowest_relaxation_s=float(-1 / lam[1]), fastest_s=float(-1 / lam[-1]), n_states=kin.M,
                       w_min=float(kin.w.min()), w_max=float(kin.w.max()))
    out["half_pulse_s"] = {f"{int(f)}Hz": 0.5 / f for f in fc.FREQS}
    return out


def main(dt, procs=4, quick=False):
    t0 = time.time()
    _init(dt)
    folds = []
    for kind in ("full", "freq", "amp"):
        for name, train, test in fc.folds(kind):
            folds.append((kind, name, train, test))
    tasks_B = [(5.0, r, k, n, tr, te) for r in ("R0", "R1", "R2", "2pop") for (k, n, tr, te) in folds]
    tasks_B += [(b, r, k, n, tr, te) for b in (2.0, 10.0) for r in ("R0", "R1", "R2") for (k, n, tr, te) in folds
                if k in ("full", "freq")]
    tasks_A = [(5.0, k, n, tr, te) for (k, n, tr, te) in folds]
    if quick:
        tasks_B, tasks_A = tasks_B[:2], tasks_A[:1]
    results = dict(dt=dt, info={b: ROUTES[b][1] for b in ROUTES}, timescales={b: timescale_diagnostics(b) for b in ROUTES})
    with Pool(procs, initializer=_init, initargs=(dt,)) as pool:
        fb = pool.map_async(task_fitB, tasks_B, chunksize=1)
        fa = pool.map_async(task_fitA, tasks_A, chunksize=1)
        results["fitB"] = fb.get()
        fc.log(f"fit B done ({len(tasks_B)} fits)")
        results["fitA"] = fa.get()
        fc.log("fit A done")
        # profiles on full data around the best full-data optimum
        prof = []
        if not quick:
            fullB = {r["route"]: r for r in results["fitB"] if r["fold"] == "full" and r["beta"] == 5.0}
            k0 = fullB["R0"]["fit"]["x"][0]
            for v in np.round(np.arange(k0 - 3, k0 + 3.01, 0.5), 3):
                if fc.LOGK_BOUNDS[0] <= v <= fc.LOGK_BOUNDS[1]:
                    prof.append((5.0, "R0", "log10_kc", float(v)))
            for rn in ("R1", "R2"):
                for v in np.arange(-3, 3.01, 1.0):
                    prof.append((5.0, rn, "log10_s", float(v)))
        results["profiles"] = pool.map(task_profile, prof, chunksize=1)
    results["wall_s"] = time.time() - t0
    os.makedirs(RES, exist_ok=True)
    out = os.path.join(RES, f"cv_results_dt{dt}{'_quick' if quick else ''}.json")
    json.dump(results, open(out, "w"), default=float)
    fc.log(f"wrote {out} in {results['wall_s']:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    main(a.dt, a.procs, a.quick)
