"""Dataset-agnostic fitting / cross-validation utilities for pulsed-voltammetry route comparisons.

Data model: ``data[(f_Hz, Esw_V)] = {"E": array (V), "fwd": uA, "bwd": uA, "net": uA}`` after the
frozen preprocessing.  A *model function* maps (params, keys) -> same structure in uA per coulomb of Q.
Q (= n F A Gamma) enters linearly and is profiled out analytically on the training curves; held-out
curves are predicted with the training Q.

Typical use (see SKILL.md):
    data = load_csv("swv_digitized.csv")                      # or build the dict yourself
    mf   = make_model_fn(lambda p: (kin, chem, chem_params(p)), data)
    best = fit(mf, data, train_keys, bounds, lhs_starts(bounds, 4, seed=1))
    pred = predict(mf, best["x"], best["Q"], test_keys)
    m    = metrics(data, test_keys, pred)
"""
import csv
from collections import defaultdict

import numpy as np
from scipy.optimize import least_squares

COMPS = ("fwd", "bwd", "net")
Q_BOUNDS = (1e-12, 1e-5)


def tau_int(x, dt, c=5.0):
    """Integrated autocorrelation time (Sokal automatic window), in units of dt."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    if acf[0] <= 0:
        return dt * 0.5
    acf /= acf[0]
    tau = 0.5
    for m in range(1, n):
        tau += acf[m]
        if m >= c * tau:
            break
    return max(tau, 0.5) * dt


def baseline_subtract(E, I, hi_edge=-0.12, lo_edge=-0.41):
    """Straight line fitted to the two edge windows (E > hi_edge, E < lo_edge) and removed.
    Freeze the windows before fitting and apply the SAME operator to data and model curves."""
    m = (E > hi_edge) | (E < lo_edge)
    p = np.polyfit(E[m], I[m], 1)
    return I - np.polyval(p, E)


def load_csv(path, comp_names=None, **baseline_kw):
    """CSV columns: f_Hz, Esw_V, curve (forward|backward|net), E_V, I_uA. Uses the E grid common to all
    three components of a condition and applies the baseline step."""
    comp_names = comp_names or {"forward": "fwd", "backward": "bwd", "net": "net"}
    raw = defaultdict(dict)
    for r in csv.DictReader(open(path)):
        key = (float(r["f_Hz"]), round(float(r["Esw_V"]), 4))
        raw[key].setdefault(comp_names[r["curve"]], []).append((float(r["E_V"]), float(r["I_uA"])))
    out = {}
    for key, d in raw.items():
        common = set.intersection(*[set(e for e, _ in d[c]) for c in COMPS])
        E = np.array(sorted(common, reverse=True))
        out[key] = dict(E=E, **{c: baseline_subtract(E, np.array([dict(d[c])[e] for e in E]), **baseline_kw)
                                for c in COMPS})
    return out


def make_model_fn(build, data, de=0.001, e_start=0.0, e_end=None, **baseline_kw):
    """build(params) -> (kin, chem, chem_params) for chem_backend, or (kin, None, (k_c, alpha, E0)) for
    swv_backend.  Returns mf(params, keys) -> {key: {comp: uA per C}} on the data grid."""
    import chem_backend as cb
    import swv_backend as sb

    def mf(params, keys):
        kin, chem, cp = build(params)
        lo = min(data[k]["E"].min() for k in keys) if e_end is None else e_end
        wfs = [sb.Waveform(f, a, e_start=e_start, e_end=lo - 1e-9, de=de) for f, a in keys]
        sims = cb.simulate(kin, chem, cp, wfs) if chem is not None else sb.simulate(kin, wfs, *cp)
        out = {}
        for key, s in zip(keys, sims):
            E = data[key]["E"]
            idx = np.rint((e_start - E) / de).astype(int)
            out[key] = {c: baseline_subtract(E, s[c][idx] * 1e6, **baseline_kw) for c in COMPS}
        return out
    return mf


def profile_Q(data, mod, keys):
    num = den = 0.0
    for k in keys:
        for c in COMPS:
            amp = np.ptp(data[k][c])
            num += np.dot(data[k][c], mod[k][c]) / amp ** 2 / len(data[k][c])
            den += np.dot(mod[k][c], mod[k][c]) / amp ** 2 / len(data[k][c])
    return float(np.clip(num / den, *Q_BOUNDS)) if den > 0 else Q_BOUNDS[0]


def residuals(data, mod, keys, Q):
    """Each curve weighted equally; each component normalized by its observed peak-to-valley."""
    r = []
    for k in keys:
        for c in COMPS:
            r.append((data[k][c] - Q * mod[k][c]) / np.ptp(data[k][c]) / np.sqrt(3 * len(data[k][c]) * len(keys)))
    return np.concatenate(r)


def lhs_starts(bounds, n, seed, overrides=None):
    """Latin-hypercube starts inside bounds; overrides = {index: (lo, hi)} narrows a coordinate's start range."""
    rng = np.random.default_rng(seed)
    lo, hi = np.array(bounds, float).T
    perms = [rng.permutation(n) for _ in lo]
    pts = []
    for i in range(n):
        u = np.array([(perms[d][i] + rng.random()) / n for d in range(len(lo))])
        p = lo + (hi - lo) * u
        for d, (a, b) in (overrides or {}).items():
            p[d] = a + (b - a) * u[d]
        pts.append(p)
    return pts


def fit(mf, data, train, bounds, starts, max_nfev=300):
    """Multi-start bounded least squares with Q profiled; returns the best solution plus every start's
    optimum (needed for the best-basin reference and identifiability checks)."""
    def fun(x):
        mod = mf(x, train)
        return residuals(data, mod, train, profile_Q(data, mod, train))

    lo, hi = np.array(bounds, float).T
    sols = []
    for x0 in starts:
        x0 = np.clip(np.asarray(x0, float), lo + 1e-9, hi - 1e-9)
        r = least_squares(fun, x0, bounds=(lo, hi), diff_step=1e-3, x_scale=(hi - lo) / 10,
                          max_nfev=max_nfev, ftol=1e-10, xtol=1e-10)
        sols.append(dict(x=r.x.tolist(), cost=float(2 * r.cost)))
    sols.sort(key=lambda s: s["cost"])
    best = dict(sols[0])
    best["Q"] = profile_Q(data, mf(best["x"], train), train)
    span = (hi - lo).max()
    best["at_bound"] = [i for i, (v, a, b) in enumerate(zip(best["x"], lo, hi)) if min(abs(v - a), abs(v - b)) < 1e-3 * span]
    best["n_near_optimal"] = sum(s["cost"] <= 1.01 * best["cost"] for s in sols)
    best["all"] = sols
    return best


def predict(mf, x, Q, keys):
    mod = mf(x, keys)
    return {k: {c: Q * mod[k][c] for c in COMPS} for k in keys}


def metrics(data, keys, pred):
    """Per held-out curve: NRMSE of each component (RMSE / observed peak-to-valley) and net peak error."""
    out = {}
    for k in keys:
        m = {f"nrmse_{c}": float(np.sqrt(np.mean((data[k][c] - pred[k][c]) ** 2)) / np.ptp(data[k][c])) for c in COMPS}
        op, pp = data[k]["net"].max(), pred[k]["net"].max()
        m["peak_rel_err_net"] = float((pp - op) / op)
        m["meets_30"] = bool(m["nrmse_net"] <= 0.30 and abs(m["peak_rel_err_net"]) <= 0.30)
        out[k] = m
    return out


def grouped_folds(keys, by):
    """Leave-one-group-out folds. by = 0 (frequency) or 1 (amplitude). Returns [(name, train, test)]."""
    groups = sorted(set(k[by] for k in keys))
    tag = "f" if by == 0 else "Esw"
    return [(f"test_{tag}{g:g}", [k for k in keys if k[by] != g], [k for k in keys if k[by] == g]) for g in groups]


def best_basin_reference(r0_mf, data, train, test, candidate_points):
    """R0*: evaluate every candidate chemistry point (from any route's starts, extra route parameters
    dropped) as R0 on the training curves, take the lowest training cost, and report its held-out error.
    Also returns the held-out range among points within 5 % of the best training cost (identifiability)."""
    rows = []
    for x in candidate_points:
        mod = r0_mf(x, train)
        Q = profile_Q(data, mod, train)
        c = float(np.sum(residuals(data, mod, train, Q) ** 2))
        m = metrics(data, test, predict(r0_mf, x, Q, test)) if test else {}
        t = float(np.mean([v["nrmse_net"] for v in m.values()])) if m else float("nan")
        rows.append(dict(x=list(x), cost=c, test=t))
    best = min(rows, key=lambda r: r["cost"])
    near = [r["test"] for r in rows if r["cost"] <= 1.05 * best["cost"]]
    return dict(best=best, near_optimal_test_range=(min(near), max(near)), n_near=len(near), all=rows)
