"""
State-restoration pipeline: from samples of a reference (balance) state and
a perturbed / disease state, possibly followed over time under an
intervention, to

  1  the STEADY STATE (balance reference, its noise floor, the disease
     displacement and whether the disease state is itself persistent),
  2  ORDER PARAMETERS (the disease coordinate phi, module coordinates, and
     the label-free collective modes that should agree with them),
  3  the MINIMAL CORE NETWORK (disrupted modules, their couplings, core genes
     replicated across groups),
  4  a module-level DYNAMICAL MODEL  dy/dt = J (y - 1) + b u  in coordinates
     where y = 1 is the disease steady state and y = 0 is balance
     (eigenvalues = recovery rates, slow mode = what lags; b = where the
     intervention acts directly),
  5  the MINIMAL ADJUSTMENT: fewest modules whose sustained adjustment moves
     the steady state back to balance (functional) and the structural driver
     set of J.

Validated on psoriasis (GSE85034), atopic dermatitis vs placebo (GSE130588)
and IBD (GSE16879) in demos/d15_state_pipeline.py; the lessons of D11-D14
are built in (projection not distance, leave-one-out cohort axis, noise
floor from repeated references, pre-specified exclusion, module = disrupted
genes only, placebo / natural course as null).

Input: `StateData` (numpy arrays; one row per sample).
"""
from dataclasses import dataclass, field

import numpy as np

from netjac import recovery as rc


# =================================================================== data
@dataclass
class StateData:
    X: np.ndarray                     # samples x features (log scale)
    features: np.ndarray              # feature names (genes)
    subject: np.ndarray               # subject id per sample
    kind: np.ndarray                  # "ref" (balance) or "state" per sample
    time: np.ndarray                  # time since perturbation / treatment start (state samples)
    arm: dict = field(default_factory=dict)       # subject -> arm ("drug", "placebo", ...)
    outcome: dict = field(default_factory=dict)   # subject -> bool (known positive)
    group: dict = field(default_factory=dict)     # subject -> replication group
    control_arms: tuple = ("placebo", "Placebo", "control")

    def subjects(self):
        return sorted(set(self.subject[self.kind == "state"]), key=str)


def _bh(p):
    p = np.nan_to_num(np.asarray(p, float), nan=1.0)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return q


# ========================================================== 1 steady state
def steady_state(D, thr=(1.0, 0.585, 0.3), q=0.05, min_genes=50, exclude_frac=0.25):
    """Balance reference per subject (own reference samples if any, else the
    cohort reference mean), disease genes (x0 - ref; BH q, |log2FC| > thr
    with fall-back thresholds), leave-one-out cohort axes, displacement in
    noise-floor units, persistence of the disease state in control arms."""
    from scipy.stats import ttest_1samp
    ref_all = D.X[D.kind == "ref"]
    cohort_ref = ref_all.mean(0)
    subs = D.subjects()
    ref, x0, personal, n_ref = {}, {}, {}, {}
    for s in subs:
        r = D.X[(D.subject == s) & (D.kind == "ref")]
        ref[s] = r.mean(0) if len(r) else cohort_ref
        personal[s] = len(r) > 0
        n_ref[s] = len(r)
        st = (D.subject == s) & (D.kind == "state")
        x0[s] = D.X[st][np.argmin(D.time[st])]
    dev = np.array([x0[s] - ref[s] for s in subs])
    t, p = ttest_1samp(dev, 0, axis=0)
    qv = _bh(p)
    for th in thr:
        dis = np.where((abs(dev.mean(0)) > th) & (qv < q))[0]
        if len(dis) >= min_genes:
            break
    tot = dev.sum(0)
    axes = {s: (tot - dev[k]) / (len(subs) - 1) for k, s in enumerate(subs)}
    a0 = {s: _proj(x0[s] - ref[s], axes[s], dis) for s in subs}
    med = np.median(list(a0.values()))
    excluded = [s for s in subs if a0[s] < exclude_frac * med]
    # noise floor along the axis: within-subject reference repeats if any,
    # else between-subject spread of reference samples (an upper bound)
    u = dev.mean(0)
    reps = [D.X[(D.subject == s) & (D.kind == "ref")] for s in subs if n_ref[s] >= 2]
    if reps:
        diffs = [_proj(r[1] - r[0], u, dis) / np.sqrt(2) for r in reps]
        noise, noise_kind = float(np.std(diffs)), "within-subject repeated reference"
    else:
        pr = np.array([_proj(r - cohort_ref, u, dis) for r in ref_all])
        noise, noise_kind = float(np.std(pr)), "between-subject reference spread (upper bound)"
    out = {"subjects": subs, "ref": ref, "x0": x0, "personal": personal, "dis": dis,
           "thr": th, "axes": axes, "a0": a0, "excluded": excluded,
           "noise": noise, "noise_kind": noise_kind,
           "displacement_z": float(np.median([a0[s] for s in subs if s not in excluded]) / max(noise, 1e-9)),
           "n_up": int((dev.mean(0)[dis] > 0).sum()), "cohort_dev": dev.mean(0)}
    return out


def _proj(v, u, idx):
    uu = u[idx] / np.linalg.norm(u[idx])
    return float(v[idx] @ uu)


# ===================================================== 2 order parameters
def phi(x, s, S, idx=None):
    """Disease coordinate: 1 = this subject's disease state, 0 = balance."""
    idx = S["dis"] if idx is None else idx
    return _proj(x - S["ref"][s], S["axes"][s], idx) / _proj(S["x0"][s] - S["ref"][s], S["axes"][s], idx)


def module_index(features, modules, S, thr=0.585, kmin=3):
    """Module genes that are actually disrupted (|cohort deviation| > thr);
    modules with < kmin such genes are not part of the state (None)."""
    gi = {g: i for i, g in enumerate(features)}
    out = {}
    for name, genes in modules.items():
        idx = [gi[g] for g in genes if g in gi and abs(S["cohort_dev"][gi[g]]) > thr]
        out[name] = np.array(idx) if len(idx) >= kmin else None
    return out


def auto_modules(D, S, k=6, seed=0):
    """Data-driven modules when no curated list exists: cluster disease
    genes by their deviation profiles across all state samples
    (correlation distance, average linkage)."""
    from scipy.cluster.hierarchy import fcluster, linkage
    rows = []
    for s in S["subjects"]:
        for x in D.X[(D.subject == s) & (D.kind == "state")]:
            rows.append((x - S["ref"][s])[S["dis"]])
    Z = np.array(rows)
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    L = linkage(Z.T, "average", metric="correlation")
    lab = fcluster(L, k, "maxclust")
    feats = D.features[S["dis"]]
    return {f"auto-{c} ({'up' if S['cohort_dev'][S['dis'][lab == c]].mean() > 0 else 'down'}, n={int((lab == c).sum())})":
            list(feats[lab == c]) for c in np.unique(lab)}


def order_parameters(D, S, modules, n_top=2000):
    """phi for every state sample (per-subject scale: 1 = own disease state),
    module coordinates (cohort scale: 1 = median patient's baseline, so a
    mild patient starts below 1 and noise is not divided by a small own
    baseline), and the label-free collective modes of the state deviations
    (do they coincide with phi?)."""
    from netjac import statmech
    midx = module_index(D.features, modules, S)
    keep = [s for s in S["subjects"] if s not in S["excluded"]]
    scale = {}
    for name, idx in midx.items():
        if idx is not None:
            u = S["cohort_dev"]
            scale[name] = float(np.median([_proj(S["x0"][s] - S["ref"][s], u, idx) for s in keep]))
    rows = []
    for s in S["subjects"]:
        if s in S["excluded"]:
            continue
        m = (D.subject == s) & (D.kind == "state")
        for x, t in zip(D.X[m], D.time[m]):
            r = {"subject": s, "time": float(t), "phi": phi(x, s, S)}
            for name, idx in midx.items():
                r[name] = (_proj(x - S["ref"][s], S["cohort_dev"], idx) / scale[name]
                           if idx is not None else np.nan)
            r["_dev"] = x - S["ref"][s]
            rows.append(r)
    dev = np.array([r["_dev"] for r in rows])
    top = np.argsort(-dev.var(0))[:n_top]
    cm = statmech.collective_modes(dev[:, top], n_null=10)
    md = statmech.modes(dev[:, top], 3)
    ph = np.array([r["phi"] for r in rows])
    r_pc = [float(abs(np.corrcoef(md["scores"][:, k], ph)[0, 1])) for k in range(3)]
    for r in rows:
        del r["_dev"]
    return {"rows": rows, "midx": midx, "k_null": cm["k_null"], "PR": cm["PR"], "top1": cm["top1"],
            "pc_r_phi": r_pc, "module_names": [n for n, i in midx.items() if i is not None]}


# ================================================= 3 minimal core network
def core_network(D, S, OP, rep_q=0.01):
    """Module couplings (partial correlations of module phis across state
    samples, Ledoit-Wolf) and core genes that are disrupted in EVERY
    replication group (if groups are given)."""
    from scipy.stats import ttest_1samp
    from sklearn.covariance import LedoitWolf
    names = OP["module_names"]
    Y = np.array([[r[n] for n in names] for r in OP["rows"]])
    Y = Y[~np.isnan(Y).any(1)]
    P = LedoitWolf().fit(Y).precision_
    pc = -P / np.sqrt(np.outer(np.diag(P), np.diag(P)))
    np.fill_diagonal(pc, 1)
    core = None
    groups = sorted(set(D.group[s] for s in S["subjects"] if s in D.group))
    if len(groups) >= 2:
        ok = np.ones(len(S["dis"]), bool)
        sign = None
        for g in groups:
            dev = np.array([(S["x0"][s] - S["ref"][s])[S["dis"]] for s in S["subjects"] if D.group.get(s) == g])
            t, p = ttest_1samp(dev, 0, axis=0)
            ok &= p < rep_q
            sign = np.sign(t) if sign is None else sign
            ok &= np.sign(t) == sign
        core = list(D.features[S["dis"][ok]])
    load = S["cohort_dev"][S["dis"]]
    top = list(D.features[S["dis"][np.argsort(-abs(load))[:30]]])
    return {"names": names, "partial": pc, "core_replicated": core, "top_axis_genes": top}


# ===================================================== 4 dynamical model
def _trajectories(OP, names, D, treated_fn):
    by = {}
    for r in OP["rows"]:
        by.setdefault(r["subject"], []).append(r)
    out = []
    for s, rs in by.items():
        rs = sorted(rs, key=lambda r: r["time"])
        Y = np.array([[r[n] for n in names] for r in rs])
        t = np.array([r["time"] for r in rs])
        if len(rs) < 2 or np.isnan(Y).any():
            continue
        out.append({"subject": s, "t": t - t[0], "Y": Y, "u": float(treated_fn(s)),
                    "control": D.arm.get(s) in D.control_arms})
    return out


def _J_from(theta, m):
    """Stable parametrisation: J = -(P P^T + eps I) + (Q - Q^T). The
    symmetric part is negative definite, so every eigenvalue has negative
    real part (a contractive prior: no transient amplification)."""
    k = m * (m + 1) // 2
    P = np.zeros((m, m))
    P[np.tril_indices(m)] = theta[:k]
    Q = np.zeros((m, m))
    Q[np.triu_indices(m, 1)] = theta[k:k + m * (m - 1) // 2]
    return -(P @ P.T + 1e-3 * np.eye(m)) + (Q - Q.T), k + m * (m - 1) // 2


def _simulate(J, b0, b, tr):
    """Exact solution of dy/dt = J (y - yD) + b0 + b u with yD = y(0) (the
    patient's untreated steady state): y(t) = y* + e^{Jt} (y0 - y*),
    y* = y0 - J^-1 (b0 + b u)."""
    from scipy.linalg import expm
    y0 = tr["Y"][0]
    ys = y0 - np.linalg.solve(J, b0 + b * tr["u"])
    return np.array([ys + expm(J * t) @ (y0 - ys) for t in tr["t"]])


def fit_dynamics(D, OP, names=None, treated_fn=None, lams=(0.03, 0.3, 3.0), drift=None):
    """Fit dy/dt = J (y - y_D) + b0 + b u to whole trajectories (exact
    exponential solution from each subject's first sample; y_D = that
    sample, i.e. the untreated patient is assumed to be at steady state).
    b0 = natural-course drift (placebo / standard care, identifiable only
    when a control arm exists), b = intervention-specific direct input.
    Ridge on the coupling parameters (lambda by leave-subject-out).
    Baselines (leave-subject-out): no change from y0, and the mean
    trajectory of the same arm at the same time (a strong baseline)."""
    from scipy.optimize import least_squares
    names = OP["module_names"] if names is None else names
    treated_fn = treated_fn or (lambda s: D.arm.get(s, "drug") not in D.control_arms)
    TR = _trajectories(OP, names, D, treated_fn)
    m = len(names)
    drift = any(t["control"] for t in TR) if drift is None else drift
    _, nJ = _J_from(np.zeros(10 ** 4), m)

    def unpack(th):
        J, _ = _J_from(th, m)
        b0 = th[nJ:nJ + m] if drift else np.zeros(m)
        b = th[nJ + m * drift:nJ + m * drift + m]
        return J, b0, b

    def theta0():
        P = 0.3 * np.eye(m)
        return np.r_[P[np.tril_indices(m)], np.zeros(m * (m - 1) // 2), np.zeros(m * drift + m)]

    def fit(trs, lam):
        k = m * (m + 1) // 2
        diag_idx = np.cumsum(np.r_[0, np.arange(2, m + 1)])[:m]     # diagonal entries of tril
        off = np.setdiff1d(np.arange(nJ), diag_idx)

        def res(th):
            J, b0, b = unpack(th)
            e = _eig(J)
            r = [(_simulate(J, b0, b, t, e)[1:] - t["Y"][1:]).ravel() for t in trs]
            return np.r_[np.concatenate(r), np.sqrt(lam) * th[off]]
        return least_squares(res, theta0(), method="trf", max_nfev=400).x

    subs = [t["subject"] for t in TR]
    cv = {}
    for lam in lams:
        err = []
        for i, t in enumerate(TR):
            th = fit(TR[:i] + TR[i + 1:], lam)
            J, b0, b = unpack(th)
            err.append(((_simulate(J, b0, b, t)[1:] - t["Y"][1:]) ** 2).ravel())
        cv[lam] = float(np.sqrt(np.mean(np.concatenate(err))))
    lam = min(cv, key=cv.get)
    J, b0, b = unpack(fit(TR, lam))
    e0, e1 = [], []
    for i, t in enumerate(TR):
        e0.append(((t["Y"][1:] - t["Y"][0]) ** 2).ravel())
        same = [o for j, o in enumerate(TR) if j != i and o["u"] == t["u"]]
        for k, tt in enumerate(t["t"][1:], 1):
            pool = [o["Y"][list(o["t"]).index(tt)] - o["Y"][0] for o in same if tt in o["t"]]
            pred = t["Y"][0] + (np.mean(pool, 0) if pool else 0)
            e1.append((pred - t["Y"][k]) ** 2)
    ev, V = np.linalg.eig(J)
    o = np.argsort(-ev.real)
    ev, V = ev[o], V[:, o]
    yD = np.mean([t["Y"][0] for t in TR], 0)
    return {"names": names, "J": J, "b0": b0, "b": b, "drift": drift, "lam": lam, "cv_rmse": cv[lam],
            "cv_all": cv, "rmse_no_change": float(np.sqrt(np.mean(np.concatenate(e0)))),
            "rmse_arm_mean": float(np.sqrt(np.mean(np.concatenate(e1)))),
            "eig": ev, "tau": -1 / ev.real, "slow_mode": np.abs(V[:, 0].real) / np.abs(V[:, 0].real).sum(),
            "yD": yD, "y_star": yD - np.linalg.solve(J, b0 + b),
            "n_subjects": len(TR), "n_points": int(sum(len(t["t"]) - 1 for t in TR))}


# ================================================== 5 minimal adjustment
def minimal_adjustment(dyn, target=0.9, max_size=None, J=None):
    """Fewest modules whose sustained direct adjustment u_S (B_S = unit
    vectors) moves the steady state from the disease state y = 1 to
    balance: y(u) = 1 - J^-1 B_S u, u by least squares.
    Restoration = 1 - |y| / |1|. Greedy forward selection; also every
    single-module lever with the pattern it produces, and the structural
    driver set of J (maximum matching)."""
    J = dyn["J"] if J is None else J
    names = dyn["names"]
    m = len(names)
    Ji = np.linalg.inv(J)
    one = dyn["yD"]                                 # cohort disease steady state

    def restore(S):
        if not S:
            return 0.0, one.copy(), np.zeros(0)
        A = Ji[:, S]
        u, *_ = np.linalg.lstsq(A, one, rcond=None)
        y = one - A @ u
        return float(1 - np.linalg.norm(y) / np.linalg.norm(one)), y, u
    singles = {}
    for k in range(m):
        r, y, u = restore([k])
        # propagation: restoration of the OTHER modules relative to their displacement
        oth = [i for i in range(m) if i != k]
        prop = float(1 - np.linalg.norm(y[oth]) / np.linalg.norm(one[oth])) if oth else 0.0
        singles[names[k]] = {"restoration": r, "pattern": 1 - y / one, "u": float(u[0]), "propagation": prop}
    chosen, path = [], []
    for _ in range(max_size or m):
        best = max((k for k in range(m) if k not in chosen), key=lambda k: restore(chosen + [k])[0])
        chosen.append(best)
        r, y, u = restore(chosen)
        path.append({"set": [names[k] for k in chosen], "restoration": r, "residual": dict(zip(names, y))})
        if r >= target:
            break
    thr = 0.1 * np.abs(J).max()
    A = (np.abs(J) > thr).astype(int)
    np.fill_diagonal(A, 0)
    drivers = [names[i] for i in rc.driver_nodes(A)]
    # the intervention actually applied, seen as a lever pattern
    applied = None
    if "b" in dyn and np.any(dyn["b"]):
        y = one - Ji @ dyn["b"]
        applied = {"restoration": float(1 - np.linalg.norm(y) / np.linalg.norm(one)),
                   "direct_target": names[int(np.argmin(dyn["b"]))],
                   "b": dict(zip(names, dyn["b"]))}
    return {"singles": singles, "path": path, "structural_drivers": drivers, "applied": applied}


# =================================================================== card
def run(D, modules, treated_fn=None, outcome_time=None):
    """All five steps; returns one dict (the 'state card')."""
    S = steady_state(D)
    OP = order_parameters(D, S, modules)
    CN = core_network(D, S, OP)
    card = {"S": {k: S[k] for k in ("thr", "noise", "noise_kind", "displacement_z", "n_up", "excluded")},
            "n_dis": len(S["dis"]), "n_subjects": len(S["subjects"]), "OP": OP, "CN": CN}
    times = sorted(set(r["time"] for r in OP["rows"]))
    if len(times) >= 2:
        dyn = fit_dynamics(D, OP, treated_fn=treated_fn)
        card["dyn"] = dyn
        card["adjust"] = minimal_adjustment(dyn)
        # persistence of the disease state without intervention
        ctrl = [r for r in OP["rows"] if D.arm.get(r["subject"]) in D.control_arms and r["time"] > 0]
        if ctrl:
            card["persistence"] = {t: float(np.median([1 - r["phi"] for r in ctrl if r["time"] == t]))
                                   for t in sorted(set(r["time"] for r in ctrl))}
    if D.outcome:
        card["outcome_auc"] = {}
        for t in times:
            rr = [r for r in OP["rows"] if r["time"] == t and r["subject"] in D.outcome]
            y = np.array([D.outcome[r["subject"]] for r in rr])
            if len(rr) > 5 and y.min() != y.max():
                card["outcome_auc"][t] = rc.auc(y, -np.array([r["phi"] for r in rr]))
    return card


def print_card(card, title=""):
    S, OP, CN = card["S"], card["OP"], card["CN"]
    print(f"\n=== state card: {title} ===")
    print(f"1 steady state: {card['n_subjects']} subjects; {card['n_dis']} disease genes (|log2FC|>{S['thr']}, "
          f"{S['n_up']} up); displacement {S['displacement_z']:.1f} x noise floor ({S['noise_kind']}); "
          f"excluded (no disease signature) {S['excluded']}")
    if "persistence" in card:
        print(f"  disease state without intervention (control arm), restoration by time: "
              + ", ".join(f"t={t:g}: {v:+.2f}" for t, v in card["persistence"].items()))
    print(f"2 order parameters: phi (1 = disease, 0 = balance); label-free modes above shuffle null {OP['k_null']}, "
          f"PR {OP['PR']:.1f}, PC1 {OP['top1']:.0%} of variance; |r(PC1..3, phi)| {np.round(OP['pc_r_phi'], 2)}")
    if card.get("outcome_auc"):
        print("  phi as outcome classifier (lower phi = better) by time: "
              + ", ".join(f"t={t:g}: {a:.2f}" for t, a in card["outcome_auc"].items()))
    print(f"3 core network: modules in the state {CN['names']}")
    if CN["core_replicated"] is not None:
        print(f"  core genes replicated in every group: {len(CN['core_replicated'])}")
    n = CN["names"]
    pc = CN["partial"]
    strong = sorted(((abs(pc[i, j]), n[i], n[j], pc[i, j]) for i in range(len(n)) for j in range(i + 1, len(n))),
                    reverse=True)[:5]
    print("  strongest module couplings (partial r): " + "; ".join(f"{a} - {b} {v:+.2f}" for _, a, b, v in strong))
    if "dyn" in card:
        d = card["dyn"]
        print(f"4 dynamics dy/dt = J(y - yD) + b0 + b u: {d['n_points']} follow-up points, {d['n_subjects']} subjects, "
              f"drift term {'on' if d['drift'] else 'off (no control arm)'}, ridge {d['lam']}; LOSO RMSE {d['cv_rmse']:.3f} "
              f"vs no-change {d['rmse_no_change']:.3f} vs same-arm mean trajectory {d['rmse_arm_mean']:.3f}")
        print(f"  time constants {np.round(d['tau'], 1)}; slow mode " +
              ", ".join(f"{a} {w:.2f}" for a, w in sorted(zip(d['names'], d['slow_mode']), key=lambda z: -z[1])[:3]))
        if d["drift"]:
            print("  natural-course drift b0: " + ", ".join(f"{a} {v:+.3f}" for a, v in zip(d['names'], d['b0'])))
        print("  direct input b (most negative = where the intervention acts): " +
              ", ".join(f"{a} {v:+.3f}" for a, v in sorted(zip(d['names'], d['b']), key=lambda z: z[1])))
        print("  cohort disease state yD -> predicted on-intervention steady state y*: " +
              ", ".join(f"{a} {u:+.2f}->{v:+.2f}" for a, u, v in zip(d['names'], d['yD'], d['y_star'])))
        a = card["adjust"]
        print("5 minimal adjustment (sustained direct adjustment of a module set):")
        for p in a["path"]:
            print(f"  {len(p['set'])} module(s) {p['set']} -> restoration {p['restoration']:.2f}")
        print("  single-module levers (restoration of whole state | of the other modules): " +
              ", ".join(f"{k} {v['restoration']:.2f}|{v['propagation']:.2f}" for k, v in
                                                    sorted(a["singles"].items(), key=lambda kv: -kv[1]["restoration"])))
        print(f"  structural driver modules (max matching on |J|): {a['structural_drivers']}")
        if a["applied"]:
            print(f"  applied intervention: model restoration {a['applied']['restoration']:.2f}, "
                  f"largest direct effect on {a['applied']['direct_target']}")
