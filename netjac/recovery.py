"""
Perturbation-recovery state definition: the evaluation and control toolkit.

The question this module answers, in the order a project should ask it:

  1. Does a perturbation / recovery readout define the state BETTER than the
     static readout? -> nested models M0 (static) < M1 (+ one stimulated
     endpoint) < M2 (+ dose or second point) < M3 (+ recovery kinetics),
     compared on the SAME out-of-fold predictions, with a paired bootstrap of
     the increment and the false-positive count at a fixed sensitivity
     (`incremental_classification`, `incremental_survival`).
  2. Which components are the core? -> `core_ranking` (critical-mode
     participation / total fragility come from `spectrum`; here the
     generic ranking of a response or Jacobian matrix).
  3. What is the smallest set of knobs that controls the network? ->
     `driver_nodes` (structural controllability, maximum matching: Liu,
     Slotine & Barabasi 2011) and `greedy_restoration_set` (functional: the
     fewest knobs whose optimal adjustment returns a defective system's
     recovery signature to the reference).

Everything is model-agnostic: features are plain matrices, the network is an
adjacency matrix, the restoration objective is any callable.
"""
import numpy as np


# ----------------------------------------------------------- classification
def fp_at_sensitivity(y, score, sens=0.8):
    """False positives (and FP rate) at the lowest threshold reaching `sens`."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(score)
    thr = np.quantile(s[y], 1 - sens)
    pred = s >= thr
    fp = int(np.sum(pred & ~y))
    return fp, fp / max(int((~y).sum()), 1), float(thr)


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def oof_predictions(X, y, groups=None, n_splits=5, seed=0, C=1.0, model=None):
    """Out-of-fold probabilities of a standardised logistic model; folds
    are split by `groups` (donor / subject) when given."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    X = np.asarray(X, float)
    y = np.asarray(y).astype(int)
    p = np.zeros(len(y))
    if groups is not None:
        splits = GroupKFold(n_splits=n_splits).split(X, y, groups)
    else:
        splits = StratifiedKFold(n_splits=n_splits, shuffle=True,
                                 random_state=seed).split(X, y)
    for tr, te in splits:
        m = model() if model else make_pipeline(
            StandardScaler(), LogisticRegression(C=C, max_iter=5000))
        m.fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p


def incremental_classification(models, y, groups=None, n_splits=5, n_boot=2000,
                               sens=0.8, seed=0, reps=1, C=1.0):
    """models: ordered dict name -> feature matrix (nested M0, M1, ...).
    Returns per-model AUC, FP at fixed sensitivity, and paired bootstrap of
    the increment over the previous model and over the first model."""
    y = np.asarray(y).astype(int)
    names = list(models)
    P = {}
    for k in names:
        P[k] = np.mean([oof_predictions(models[k], y, groups, n_splits, seed + r, C)
                        for r in range(reps)], axis=0)
    rng = np.random.default_rng(seed)
    n = len(y)
    if groups is not None:
        ug = np.unique(groups)
        gidx = {g: np.where(groups == g)[0] for g in ug}
    out = {}
    base = names[0]
    for i, k in enumerate(names):
        fp, fpr, _ = fp_at_sensitivity(y, P[k], sens)
        res = {"auc": auc(y, P[k]), "fp": fp, "fpr": fpr, "pred": P[k]}
        if i > 0:
            prev = names[i - 1]
            d_prev, d_base, d_fp = [], [], []
            for _ in range(n_boot):
                if groups is not None:
                    gs = rng.choice(ug, len(ug), replace=True)
                    b = np.concatenate([gidx[g] for g in gs])
                else:
                    b = rng.integers(0, n, n)
                if y[b].min() == y[b].max():
                    continue
                a_k, a_p, a_b = auc(y[b], P[k][b]), auc(y[b], P[prev][b]), auc(y[b], P[base][b])
                d_prev.append(a_k - a_p)
                d_base.append(a_k - a_b)
                d_fp.append(fp_at_sensitivity(y[b], P[k][b], sens)[1]
                            - fp_at_sensitivity(y[b], P[prev][b], sens)[1])
            res["d_prev"] = (np.mean(d_prev), *np.percentile(d_prev, [2.5, 97.5]),
                             float(np.mean(np.array(d_prev) <= 0)))
            res["d_base"] = (np.mean(d_base), *np.percentile(d_base, [2.5, 97.5]))
            res["d_fpr"] = (np.mean(d_fp), *np.percentile(d_fp, [2.5, 97.5]))
        out[k] = res
    return out


# ----------------------------------------------------------------- survival
def harrell_c(time, event, risk):
    """Harrell's concordance (ties in risk count 1/2)."""
    time, event, risk = map(np.asarray, (time, event, risk))
    o = np.argsort(time)
    time, event, risk = time[o], event[o], risk[o]
    num = den = 0.0
    for i in np.where(event == 1)[0]:
        later = time > time[i]
        if not later.any():
            continue
        r = risk[later]
        num += np.sum(r < risk[i]) + 0.5 * np.sum(r == risk[i])
        den += later.sum()
    return num / den


def oof_cox(X, time, event, n_splits=5, seed=0):
    from sklearn.model_selection import KFold
    import statsmodels.api as smapi
    X = np.asarray(X, float)
    lp = np.zeros(len(time))
    for tr, te in KFold(n_splits, shuffle=True, random_state=seed).split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1
        mod = smapi.PHReg(time[tr], (X[tr] - mu) / sd, status=event[tr])
        b = mod.fit().params
        if not np.all(np.isfinite(b)):              # rare-event folds
            b = mod.fit(method="bfgs", maxiter=500).params
        if not np.all(np.isfinite(b)):
            b = np.zeros(X.shape[1])
        lp[te] = ((X[te] - mu) / sd) @ b
    return lp


def incremental_survival(models, time, event, n_splits=5, n_boot=500, seed=0,
                         horizon=None, sens=0.8):
    """Nested Cox models compared by out-of-fold Harrell C (paired bootstrap)
    and, if `horizon` is given, by FP at fixed sensitivity for death before
    the horizon among people followed to it or dead before it."""
    time = np.asarray(time, float)
    event = np.asarray(event).astype(int)
    names = list(models)
    LP = {k: oof_cox(models[k], time, event, n_splits, seed) for k in names}
    rng = np.random.default_rng(seed)
    n = len(time)
    out = {}
    if horizon is not None:
        known = (event == 1) & (time <= horizon) | (time > horizon)
        yh = ((event == 1) & (time <= horizon))[known]
    for i, k in enumerate(names):
        res = {"c": harrell_c(time, event, LP[k]), "lp": LP[k]}
        if horizon is not None:
            res["auc_h"] = auc(yh, LP[k][known])
            res["fp_h"] = fp_at_sensitivity(yh, LP[k][known], sens)
        if i > 0:
            prev = names[i - 1]
            ds = []
            sub = rng.integers(0, n, (n_boot, min(n, 3000)))
            for b in sub:
                ds.append(harrell_c(time[b], event[b], LP[k][b])
                          - harrell_c(time[b], event[b], LP[prev][b]))
            res["d_prev"] = (np.mean(ds), *np.percentile(ds, [2.5, 97.5]),
                             float(np.mean(np.array(ds) <= 0)))
            if horizon is not None:
                dfp, dauc = [], []
                yk, Lk, Lp = yh, LP[k][known], LP[prev][known]
                for _ in range(n_boot):
                    b = rng.integers(0, len(yk), len(yk))
                    dauc.append(auc(yk[b], Lk[b]) - auc(yk[b], Lp[b]))
                    dfp.append(fp_at_sensitivity(yk[b], Lk[b], sens)[1]
                               - fp_at_sensitivity(yk[b], Lp[b], sens)[1])
                res["d_auc_h"] = (np.mean(dauc), *np.percentile(dauc, [2.5, 97.5]))
                res["d_fpr_h"] = (np.mean(dfp), *np.percentile(dfp, [2.5, 97.5]))
        out[k] = res
    return out


# --------------------------------------------------------- network control
def driver_nodes(adj):
    """Minimum driver-node set of a directed network by structural
    controllability: drivers = nodes left unmatched by a maximum matching of
    the bipartite out-copy -> in-copy graph (Liu, Slotine & Barabasi 2011).
    adj[i, j] != 0 means an edge j -> i (column = source)."""
    import networkx as nx
    from networkx.algorithms import bipartite
    A = np.asarray(adj) != 0
    n = A.shape[0]
    G = nx.Graph()
    G.add_nodes_from((("o", j) for j in range(n)), bipartite=0)
    G.add_nodes_from((("i", i) for i in range(n)), bipartite=1)
    tgt, src = np.nonzero(A)
    G.add_edges_from((("o", s), ("i", t)) for t, s in zip(tgt, src) if s != t)
    M = bipartite.hopcroft_karp_matching(G, top_nodes=[("o", j) for j in range(n)])
    matched_in = {v[1] for k, v in M.items() if k[0] == "o"}
    drivers = [i for i in range(n) if i not in matched_in]
    return drivers


def core_ranking(R, names, k=20):
    """Rank components by total directed influence: out-strength (how much
    their perturbation moves the rest) and in-strength (how much they are
    moved); returns the two top-k lists."""
    R = np.abs(np.asarray(R, float))
    np.fill_diagonal(R, 0)
    out_s, in_s = R.sum(0), R.sum(1)
    return ([(names[i], float(out_s[i])) for i in np.argsort(-out_s)[:k]],
            [(names[i], float(in_s[i])) for i in np.argsort(-in_s)[:k]])


def greedy_restoration_set(knobs, objective, max_size=4, x0=None, bounds=(-2, 2),
                           starts=(0.0,)):
    """Fewest knobs whose optimal (log-scale) adjustment minimises
    objective(dict knob -> log multiplier). Greedy forward selection; at
    each size the chosen knobs are re-optimised jointly.
    `starts`: initial values tried for the newly added knob (multi-start;
    objectives with a flat region near 0 need it).
    Returns a list of (set, adjustments, objective value)."""
    from scipy.optimize import minimize
    chosen, path = [], [((), {}, float(objective({})))]
    for _ in range(max_size):
        best = None
        for k in knobs:
            if k in chosen:
                continue
            ks = chosen + [k]
            f = lambda v: objective(dict(zip(ks, v)))
            prev = dict(path[-1][1])
            for s0 in starts:
                v0 = np.array([prev.get(q, s0 if q == k else 0.0) for q in ks])
                r = minimize(f, v0, method="L-BFGS-B", bounds=[bounds] * len(ks))
                if best is None or r.fun < best[2]:
                    best = (tuple(ks), dict(zip(ks, r.x)), float(r.fun))
        chosen = list(best[0])
        path.append(best)
    return path
