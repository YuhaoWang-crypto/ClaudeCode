"""Round 2: three-route comparison under stepwise MB chemistry (C1 EE, C2 E-C-E).

    python run_chem.py              # writes results/chem_results_dt0.5.json
Data, preprocessing, folds, objective, Q profiling and metrics are those of round 1 (fit_cv.py).
Stage 1: R0 for every chemistry/fold (4 LHS starts) and the K_p sensitivity of C2-R0.
Stage 2: R1 (20 cells; 40/80-cell recheck at the optimum), R2 (conditional), two-population control
(3 LHS + 1 warm start from the same-fold R0 solution) and fit A (warm start + 1 LHS).
"""
import json
import os
import time
import zlib
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

import chem_backend as cb
import fit_cv as fc
import md_analysis as ma
import swv_backend as sb
from run_analysis import D_SE_NM2_S

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DT = 0.5
BETA = 5.0
DE_B = (-0.3, 0.3)
LOGKP_B = (0.0, 8.0)
BASE = {
    "C1": (["log10_k1", "log10_k2", "alpha", "Em", "dE"], [fc.LOGK_BOUNDS, fc.LOGK_BOUNDS, fc.ALPHA_BOUNDS, fc.E0_BOUNDS, DE_B]),
    "C2": (["log10_k1", "log10_k2", "alpha", "Em", "dE", "log10_kp"],
           [fc.LOGK_BOUNDS, fc.LOGK_BOUNDS, fc.ALPHA_BOUNDS, fc.E0_BOUNDS, DE_B, LOGKP_B]),
}
KIN = {}


def init():
    names = [f"dt{DT}_s{s}" for s in (101, 102, 103, 104)]
    z = np.concatenate([ma.load(n)["z"].ravel() for n in names])
    w = np.exp(-BETA * (z - ma.Z_MIN))
    r2 = np.load(os.path.join(RES, f"r2_model_dt{DT}.npz"))
    KIN["R0"] = sb.fast_eq_kinetics([1.0], [w.mean()])
    for n in (20, 40, 80):
        KIN[f"R1_{n}"] = sb.smoluchowski_kinetics(z, n, D_SE_NM2_S, BETA, ma.Z_MIN)
    KIN["R2"] = sb.Kinetics(r2["pooled_pi"], r2[f"pooled_w{BETA:g}"], r2["pooled_L"], name="R2")


class ChemRoute:
    def __init__(self, chem, route, K_p=10.0, r1_cells=20):
        self.chem_name, self.route, self.K_p, self.r1_cells = chem, route, K_p, r1_cells
        self.chem = cb.Chemistry(chem, K_p=K_p)

    def names(self):
        n = list(BASE[self.chem_name][0])
        if self.route in ("R1", "R2"):
            n.append("log10_s")
        if self.route == "2pop":
            n.append("delta")
        return n

    def bounds(self):
        b = list(BASE[self.chem_name][1])
        if self.route in ("R1", "R2"):
            b.append(fc.LOGS_BOUNDS)
        if self.route == "2pop":
            b.append(fc.DELTA_BOUNDS)
        return b

    def unpack(self, x, cells=None):
        x = list(x)
        p = dict(k1=10 ** x[0], k2=10 ** x[1], alpha=x[2], Em=x[3], dE=x[4])
        if self.chem_name == "C2":
            p["kp"] = 10 ** x[5]
        nb = len(BASE[self.chem_name][0])
        if self.route == "R0":
            kin = KIN["R0"]
        elif self.route == "R1":
            kin = KIN[f"R1_{cells or self.r1_cells}"].scaled(10 ** x[nb])
        elif self.route == "R2":
            kin = KIN["R2"].scaled(10 ** x[nb])
        else:
            kin = sb.two_population_kinetics(x[nb])
        return kin, p


def model_curves(cr, x, keys, cells=None):
    kin, p = cr.unpack(x, cells)
    sims = cb.simulate(kin, cr.chem, p, [sb.Waveform(f, a) for f, a in keys])
    out = {}
    for key, s in zip(keys, sims):
        E = fc.data()[key]["E"]
        idx = np.rint(-E / 0.001).astype(int)
        out[key] = {c: fc.baseline_subtract(E, s[c][idx] * 1e6) for c in fc.COMPS}
    return out


def lhs_starts(cr, n, seed):
    rng = np.random.default_rng(seed)
    lo, hi = np.array(cr.bounds()).T
    pts = []
    perms = [rng.permutation(n) for _ in lo]          # Latin hypercube: one stratum per start and dimension
    for i in range(n):
        u = np.array([(perms[d][i] + rng.random()) / n for d in range(len(lo))])
        p = lo + (hi - lo) * u
        p[0] = lo[0] + 1 + (hi[0] - lo[0] - 2) * u[0]
        p[1] = lo[1] + 1 + (hi[1] - lo[1] - 2) * u[1]
        p[2] = 0.3 + 0.4 * u[2]
        p[3] = -0.27 + 0.02 * (rng.random() - 0.5)
        pts.append(p)
    return pts


def fit(cr, train, starts):
    obs = fc.data()

    def fun(x):
        mod = model_curves(cr, x, train)
        return fc.residuals(obs, mod, train, fc.profile_Q(obs, mod, train))

    lo, hi = np.array(cr.bounds()).T
    sols = []
    for x0 in starts:
        x0 = np.clip(np.asarray(x0, float), lo + 1e-9, hi - 1e-9)
        try:
            r = least_squares(fun, x0, bounds=(lo, hi), diff_step=1e-3, x_scale=(hi - lo) / 10,
                              max_nfev=300, ftol=1e-10, xtol=1e-10)
            sols.append(dict(x=r.x.tolist(), cost=float(2 * r.cost), nfev=int(r.nfev), x0=x0.tolist()))
        except Exception as e:      # logged, not hidden
            sols.append(dict(x=x0.tolist(), cost=float("inf"), error=repr(e), x0=x0.tolist()))
    sols.sort(key=lambda s: s["cost"])
    best = dict(sols[0])
    mod = model_curves(cr, best["x"], train)
    best["Q"] = fc.profile_Q(obs, mod, train)
    span = (hi - lo).max()
    best["at_bound"] = [n for n, v, l, h in zip(cr.names(), best["x"], lo, hi) if min(abs(v - l), abs(v - h)) < 1e-3 * span]
    best["n_starts"] = len(sols)
    best["n_near_optimal"] = sum(s["cost"] <= best["cost"] * 1.01 for s in sols)
    best["all_costs"] = [s["cost"] for s in sols]
    best["all_x"] = [s["x"] for s in sols]
    return best


def predict(cr, x, Q, keys, cells=None):
    mod = model_curves(cr, x, keys, cells)
    return {k: {c: Q * mod[k][c] for c in fc.COMPS} for k in keys}


def _key(k):
    return f"f{int(k[0])}_A{int(round(k[1] * 1000))}"


def task_fit(args):
    chem, route, fname, train, test, warm, K_p = args
    t = time.time()
    cr = ChemRoute(chem, route, K_p=K_p)
    seed = zlib.crc32(f"{chem}|{route}|{fname}|{K_p}".encode()) % 10000
    starts = lhs_starts(cr, 4 if warm is None else 3, seed)
    if warm is not None:
        w = list(warm)
        if route in ("R1", "R2"):
            w = w + [0.0]
        elif route == "2pop":
            w = w + [0.5]
        starts.append(np.array(w))
    best = fit(cr, train, starts)
    allk = train + test
    pred = predict(cr, best["x"], best["Q"], allk)
    res = dict(kind="B", chem=chem, route=route, K_p=K_p, fold=fname, train=train, test=test, fit=best,
               param_names=cr.names(), n_params=len(cr.names()) + 1,
               train_metrics=fc.metrics(train, pred), test_metrics=fc.metrics(test, pred) if test else {},
               curves={_key(k): {c: pred[k][c].tolist() for c in fc.COMPS} for k in allk})
    if route == "R1":
        for n in (40, 80):
            pn = predict(cr, best["x"], best["Q"], allk, cells=n)
            res[f"grid{n}_max_rel_diff"] = float(max(np.abs(pn[k][c] - pred[k][c]).max() / np.ptp(pred[k][c])
                                                     for k in allk for c in fc.COMPS))
            res[f"grid{n}_test_metrics"] = fc.metrics(test, pn) if test else {}
    res["wall_s"] = time.time() - t
    return res


def task_fitA(args):
    chem, fname, train, test, warm = args
    t = time.time()
    obs = fc.data()
    crs = [ChemRoute(chem, r) for r in ("R0", "R1", "R2")]

    def ext(cr, x):
        return list(x) + ([0.0] if cr.route in ("R1", "R2") else [])

    def fun(x):
        ms = [model_curves(cr, ext(cr, x), train) for cr in crs]
        num = den = 0.0
        for m in ms:
            for k in train:
                for c in fc.COMPS:
                    amp = np.ptp(obs[k][c])
                    num += np.dot(obs[k][c], m[k][c]) / amp ** 2 / len(obs[k][c])
                    den += np.dot(m[k][c], m[k][c]) / amp ** 2 / len(obs[k][c])
        Q = float(np.clip(num / den, *fc.Q_BOUNDS)) if den > 0 else fc.Q_BOUNDS[0]
        return np.concatenate([fc.residuals(obs, m, train, Q) / np.sqrt(3) for m in ms]), Q

    lo, hi = np.array(BASE[chem][1]).T
    starts = [np.array(warm)] + lhs_starts(crs[0], 1, zlib.crc32(f"A|{chem}|{fname}".encode()) % 10000)
    sols = []
    for x0 in starts:
        x0 = np.clip(np.asarray(x0, float)[: len(lo)], lo + 1e-9, hi - 1e-9)
        r = least_squares(lambda x: fun(x)[0], x0, bounds=(lo, hi), diff_step=1e-3, x_scale=(hi - lo) / 10,
                          max_nfev=300, ftol=1e-10, xtol=1e-10)
        sols.append(dict(x=r.x.tolist(), cost=float(2 * r.cost)))
    sols.sort(key=lambda s: s["cost"])
    best = dict(sols[0], all_costs=[s["cost"] for s in sols])
    best["Q"] = fun(np.array(best["x"]))[1]
    allk = train + test
    out = dict(kind="A", chem=chem, fold=fname, train=train, test=test, fit=best, routes={})
    for cr in crs:
        pred = predict(cr, ext(cr, best["x"]), best["Q"], allk)
        out["routes"][cr.route] = dict(train_metrics=fc.metrics(train, pred),
                                       test_metrics=fc.metrics(test, pred) if test else {})
    out["wall_s"] = time.time() - t
    return out


def main(procs=4):
    t0 = time.time()
    init()
    folds = [f for kind in ("full", "freq", "amp") for f in fc.folds(kind)]
    with Pool(procs, initializer=init) as pool:
        s1 = [(ch, "R0", n, tr, te, None, 10.0) for ch in ("C1", "C2") for (n, tr, te) in folds]
        s1 += [("C2", "R0", n, tr, te, None, kp) for kp in (1.0, 100.0) for (n, tr, te) in folds
               if n == "full" or n.startswith("test_f")]
        r1 = pool.map(task_fit, s1, chunksize=1)
        fc.log(f"stage 1 done ({len(r1)} R0 fits)")
        warm = {(r["chem"], r["fold"]): r["fit"]["x"] for r in r1 if r["K_p"] == 10.0}
        s2 = [(ch, rt, n, tr, te, warm[(ch, n)], 10.0) for rt in ("R1", "R2", "2pop") for ch in ("C2", "C1")
              for (n, tr, te) in folds]
        sA = [(ch, n, tr, te, warm[(ch, n)]) for ch in ("C2", "C1") for (n, tr, te) in folds
              if n == "full" or n.startswith("test_f")]
        b2 = pool.map_async(task_fit, s2, chunksize=1)
        a2 = pool.map_async(task_fitA, sA, chunksize=1)
        r2 = b2.get()
        fc.log(f"stage 2 fit B done ({len(r2)} fits)")
        ra = a2.get()
        fc.log("stage 2 fit A done")
    res = dict(dt=DT, beta=BETA, fitB=r1 + r2, fitA=ra, wall_s=time.time() - t0)
    os.makedirs(RES, exist_ok=True)
    json.dump(res, open(os.path.join(RES, f"chem_results_dt{DT}.json"), "w"), default=float)
    fc.log(f"wrote chem_results in {res['wall_s']:.0f} s")


if __name__ == "__main__":
    main()
