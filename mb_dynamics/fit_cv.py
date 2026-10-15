"""Fits A (shared theta) and B (per-route calibration), grouped cross-validation, profiles.

Data: 9 digitized DNA20-MB SWV conditions (forward, backward, net), frozen baseline step applied
identically to data and model curves.  Q (= nFA Gamma) enters linearly and is profiled out
analytically (bounded 1e-12 ... 1e-5 C) on the training curves; held-out curves are predicted with
the training Q.
"""
import csv
import json
import os
import time
from collections import defaultdict

import numpy as np
from scipy.optimize import least_squares

import swv_backend as sb

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "dauphin2017", "swv_digitized.csv")
COMPS = ("fwd", "bwd", "net")
CSV_COMP = {"forward": "fwd", "backward": "bwd", "net": "net"}
FREQS = (20.0, 100.0, 500.0)
AMPS = (0.025, 0.050, 0.075)
Q_BOUNDS = (1e-12, 1e-5)
LOGK_BOUNDS = (-3.0, 8.0)
ALPHA_BOUNDS = (0.1, 0.9)
E0_BOUNDS = (-0.45, 0.0)
LOGS_BOUNDS = (-3.0, 3.0)
DELTA_BOUNDS = (0.0, float(np.log(1e4)))


def baseline_subtract(E, I):
    m = (E > -0.12) | (E < -0.41)
    p = np.polyfit(E[m], I[m], 1)
    return I - np.polyval(p, E)


def load_data():
    raw = defaultdict(dict)
    for r in csv.DictReader(open(DATA)):
        key = (float(r["f_Hz"]), round(float(r["Esw_V"]), 3))
        raw[key].setdefault(CSV_COMP[r["curve"]], []).append((float(r["E_V"]), float(r["I_uA"])))
    data = {}
    for key, d in raw.items():
        # common grid of the three components (a few points have no visible or imputable trace)
        common = set.intersection(*[set(e for e, _ in d[c]) for c in COMPS])
        E = np.array(sorted(common, reverse=True))
        cur = {}
        for c in COMPS:
            dd = dict(d[c])
            I = np.array([dd[e] for e in E])
            cur[c] = baseline_subtract(E, I)
        data[key] = dict(E=E, **cur)
    return data


DATA_CACHE = None


def data():
    global DATA_CACHE
    if DATA_CACHE is None:
        DATA_CACHE = load_data()
    return DATA_CACHE


def model_curves(kin, keys, k_c, alpha, e0):
    """Unit-Q model currents (A per C) at the data grid, with the same baseline step."""
    wfs = [sb.Waveform(f, a) for f, a in keys]
    sims = sb.simulate(kin, wfs, k_c, alpha, e0)
    out = {}
    for key, s in zip(keys, sims):
        E = data()[key]["E"]
        idx = np.rint(-E / 0.001).astype(int)
        out[key] = {c: baseline_subtract(E, s[c][idx] * 1e6) for c in COMPS}   # uA per C
    return out


def profile_Q(obs, mod, keys):
    num = den = 0.0
    for k in keys:
        for c in COMPS:
            amp = np.ptp(obs[k][c])
            num += np.dot(obs[k][c], mod[k][c]) / amp ** 2 / len(obs[k][c])
            den += np.dot(mod[k][c], mod[k][c]) / amp ** 2 / len(obs[k][c])
    if den <= 0:
        return Q_BOUNDS[0]
    return float(np.clip(num / den, *Q_BOUNDS))


def residuals(obs, mod, keys, Q):
    r = []
    for k in keys:
        for c in COMPS:
            amp = np.ptp(obs[k][c])
            r.append((obs[k][c] - Q * mod[k][c]) / amp / np.sqrt(3 * len(obs[k][c]) * len(keys)))
    return np.concatenate(r)


# ---------------------------------------------------------------------------------------------
# Route parameterizations.  Each returns (kin, k_c) for a parameter vector.
# ---------------------------------------------------------------------------------------------

class Route:
    def __init__(self, name, kin=None, has_s=False, two_pop=False):
        self.name, self.kin, self.has_s, self.two_pop = name, kin, has_s, two_pop

    def bounds(self):
        if self.two_pop:
            return [LOGK_BOUNDS, DELTA_BOUNDS, ALPHA_BOUNDS, E0_BOUNDS]
        b = [LOGK_BOUNDS, ALPHA_BOUNDS, E0_BOUNDS]
        return b + [LOGS_BOUNDS] if self.has_s else b

    def unpack(self, p):
        if self.two_pop:
            return sb.two_population_kinetics(p[1]), 10 ** p[0], p[2], p[3]
        kin = self.kin.scaled(10 ** p[3]) if self.has_s else self.kin
        return kin, 10 ** p[0], p[1], p[2]

    def names(self):
        if self.two_pop:
            return ["log10_k", "delta", "alpha", "E0"]
        return ["log10_kc", "alpha", "E0"] + (["log10_s"] if self.has_s else [])


def fit_route(route, train_keys, starts, fixed_shared=None):
    obs = data()

    def fun(p):
        kin, kc, a, e0 = route.unpack(p)
        mod = model_curves(kin, train_keys, kc, a, e0)
        Q = profile_Q(obs, mod, train_keys)
        return residuals(obs, mod, train_keys, Q)

    lo, hi = np.array(route.bounds()).T
    sols = []
    for x0 in starts:
        x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
        try:
            r = least_squares(fun, x0, bounds=(lo, hi), diff_step=1e-3, x_scale=(hi - lo) / 10,
                              max_nfev=400, ftol=1e-10, xtol=1e-10)
            sols.append(dict(x=r.x.tolist(), cost=float(2 * r.cost), nfev=int(r.nfev), x0=list(map(float, x0))))
        except Exception as e:      # pragma: no cover - logged, not hidden
            sols.append(dict(x=list(map(float, x0)), cost=float("inf"), error=repr(e), x0=list(map(float, x0))))
    sols.sort(key=lambda s: s["cost"])
    best = sols[0]
    kin, kc, a, e0 = route.unpack(np.array(best["x"]))
    mod = model_curves(kin, train_keys, kc, a, e0)
    best["Q"] = profile_Q(obs, mod, train_keys)
    lo_hit = [n for n, v, l in zip(route.names(), best["x"], lo) if abs(v - l) < 1e-3 * (hi - lo).max()]
    hi_hit = [n for n, v, h in zip(route.names(), best["x"], hi) if abs(v - h) < 1e-3 * (hi - lo).max()]
    best["at_bound"] = lo_hit + hi_hit
    best["Q_at_bound"] = best["Q"] in Q_BOUNDS
    # multi-start spread: distinct optima within 1 % of the best cost
    near = [s for s in sols if s["cost"] <= best["cost"] * 1.01]
    best["n_starts"] = len(sols)
    best["n_near_optimal"] = len(near)
    best["near_optimal_x"] = [s["x"] for s in near]
    best["all_costs"] = [s["cost"] for s in sols]
    return best


def predict(route, x, Q, keys):
    kin, kc, a, e0 = route.unpack(np.array(x))
    mod = model_curves(kin, keys, kc, a, e0)
    return {k: {c: Q * mod[k][c] for c in COMPS} for k in keys}


def metrics(keys, pred):
    obs = data()
    out = {}
    for k in keys:
        m = {}
        for c in COMPS:
            m[f"nrmse_{c}"] = float(np.sqrt(np.mean((obs[k][c] - pred[k][c]) ** 2)) / np.ptp(obs[k][c]))
        op, pp = obs[k]["net"].max(), pred[k]["net"].max()
        m["peak_rel_err_net"] = float((pp - op) / op)
        m["peak_E_obs"] = float(obs[k]["E"][obs[k]["net"].argmax()])
        m["peak_E_pred"] = float(obs[k]["E"][pred[k]["net"].argmax()])
        out[f"f{int(k[0])}_A{int(round(k[1] * 1000))}"] = m
    return out


def folds(kind):
    allk = [(f, a) for f in FREQS for a in AMPS]
    if kind == "freq":
        return [(f"test_f{int(f)}", [k for k in allk if k[0] != f], [k for k in allk if k[0] == f]) for f in FREQS]
    if kind == "amp":
        return [(f"test_A{int(a * 1000)}", [k for k in allk if k[1] != a], [k for k in allk if k[1] == a]) for a in AMPS]
    return [("full", allk, [])]


def starts_for(route, n=6, seed=0):
    rng = np.random.default_rng(seed)
    lo, hi = np.array(route.bounds()).T
    pts = []
    # informed start: contact rate giving an effective rate ~1e2-1e3 s^-1 is not known per route, so
    # starts span the k range on a grid and the remaining coordinates are Latin-hypercube samples
    kgrid = np.linspace(lo[0] + 1, hi[0] - 1, n)
    rng.shuffle(kgrid)
    for i in range(n):
        u = (rng.permutation(n)[i] + rng.random()) / n
        p = lo + (hi - lo) * rng.random(len(lo))
        p[0] = kgrid[i]
        if route.two_pop:
            p[2], p[3] = 0.3 + 0.4 * u, -0.27 + 0.02 * (rng.random() - 0.5)
        else:
            p[1], p[2] = 0.3 + 0.4 * u, -0.27 + 0.02 * (rng.random() - 0.5)
        pts.append(p)
    return pts


def fit_shared(routes, train_keys, starts):
    """Fit A: one (k_c, alpha, E0, Q) shared by all routes (s = 1); routes and curves equally weighted."""
    obs = data()

    def mods(p):
        return [model_curves(r.kin, train_keys, 10 ** p[0], p[1], p[2]) for r in routes]

    def fun(p):
        ms = mods(p)
        num = den = 0.0
        for m in ms:     # joint Q over routes
            for k in train_keys:
                for c in COMPS:
                    amp = np.ptp(obs[k][c])
                    num += np.dot(obs[k][c], m[k][c]) / amp ** 2 / len(obs[k][c])
                    den += np.dot(m[k][c], m[k][c]) / amp ** 2 / len(obs[k][c])
        Q = float(np.clip(num / den, *Q_BOUNDS)) if den > 0 else Q_BOUNDS[0]
        return np.concatenate([residuals(obs, m, train_keys, Q) / np.sqrt(len(routes)) for m in ms]), Q

    lo, hi = np.array([LOGK_BOUNDS, ALPHA_BOUNDS, E0_BOUNDS]).T
    sols = []
    for x0 in starts:
        x0 = np.clip(np.asarray(x0)[:3], lo + 1e-9, hi - 1e-9)
        r = least_squares(lambda p: fun(p)[0], x0, bounds=(lo, hi), diff_step=1e-3, x_scale=(hi - lo) / 10,
                          max_nfev=400, ftol=1e-10, xtol=1e-10)
        sols.append(dict(x=r.x.tolist(), cost=float(2 * r.cost), x0=x0.tolist()))
    sols.sort(key=lambda s: s["cost"])
    best = sols[0]
    best["Q"] = fun(np.array(best["x"]))[1]
    best["all_costs"] = [s["cost"] for s in sols]
    return best


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)
