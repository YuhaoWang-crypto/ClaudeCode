"""
D9  Micro -> macro: statistical physics of omics ensembles, and where the
    network Jacobian shows up in them.

The question (Boltzmann-style): from the DISTRIBUTION of micro-states (cells,
tissue samples, single-cell transcriptomes), can we derive macro state
variables (like magnetisation from spins, pressure from velocities), find the
couplings, tell positive from negative feedback, and find fragile points and
control levers?

Parts
  A  Theory, exact (all checks vs Lyapunov / closed forms):
     - 2-node loop: C12 = (k12 C22 + k21 C11)/(a+b). Positive feedback is
       visible in static correlations; negative feedback cancels itself.
     - Fluctuation-response: with detailed balance the response to a field
       equals the covariance column. A negative-feedback homeostat, acting
       through a hidden sensed variable, REVERSES the sign: partners that
       co-fluctuate with gene j go UP when j is knocked down ("anti-FDT").
  B  Insulin-glucose from micro to macro: a beta-cell POPULATION with
     heterogeneous glucose thresholds -> the macro secretion curve -> the
     glucose-insulin-beta-mass Jacobian (Topp 2000) -> when the system
     converges (compensation) vs diverges (beta-cell collapse), and which
     lever moves it back.
  C  8 indications, CPTAC tumour vs normal proteomes (6-10k proteins each):
     - macro descriptors (effective dimension, coupling strength, order
       parameters);
     - mode identity and universality;
     - prior-network couplings (complexes, signed edges, 2-loops);
     - all-pairs rewiring (3.5-5 x 10^7 pairs per indication) with a
       permutation null, before and after subtracting the macro fields;
     - Boltzmann "energy" vs grade.
  D  Directed Jacobian from perturbations: Replogle 2022 K562 Perturb-seq.
     - FDT vs anti-FDT classification of 1,970 knockdowns;
     - sparse modular response analysis with split-half replication;
     - feedback loops;
     - macro response modes and their universality (K562 vs RPE1).
  E  Islet single-cell RNA (Xin 2016, Segerstolpe 2016): can the donor-level
     beta-cell population descriptors be measured today? (negative)

Run:  python3 -m netjac.demos.d9_micro_to_macro [--part A,B,C,D,E]
"""
import argparse
import os
import pickle
import re
import time

import numpy as np
from scipy import linalg
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

from netjac import DATA_CACHE, FIG
from netjac import statmech as sm
from netjac import statespace as ss

RES = os.path.join(DATA_CACHE, "d9_results.pkl")


def _load_res():
    if os.path.exists(RES):
        return pickle.load(open(RES, "rb"))
    return {}


def _save_res(R):
    os.makedirs(DATA_CACHE, exist_ok=True)
    pickle.dump(R, open(RES, "wb"))


# ======================================================================= A
def part_A():
    out = {}
    # A1 closed form vs Lyapunov
    errs = []
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b = rng.uniform(0.3, 2, 2)
        k12, k21 = rng.uniform(-1.5, 1.5, 2)
        o = sm.ou_two_node(a, b, k12, k21, *rng.uniform(0.3, 2, 2))
        if o["stable"]:
            errs.append(abs(o["C"][0, 1] - o["C12_closed"]))
    out["closed_form_err"] = max(errs)
    # A2 loop visibility
    t = sm.loop_correlation_theory(20000)
    out["loop_theory"] = {k: np.abs(v) for k, v in t.items()}
    # A3 hidden-sensor homeostat: E (enzyme mRNA), M (sensed metabolite,
    # unmeasured), T (another target of the same TF). The TF drive falls
    # with M (SREBP-sterol / UPR-load / NRF1-proteasome logic) and TF noise
    # is shared by E and T (correlated noise, D_ET > 0).
    g, kf = 1.0, 1.0
    J = np.array([[-1.0, -g, 0.0],
                  [kf, -1.0, 0.0],
                  [0.0, -g, -1.0]])
    D = np.array([[1.0, 0.0, 0.8], [0.0, 0.2, 0.0], [0.8, 0.0, 1.0]])
    C = linalg.solve_continuous_lyapunov(J, -D)
    resp = -np.linalg.solve(J, np.array([-1.0, 0.0, 0.0]))   # knock E down
    out["homeostat"] = {"corr_ET": C[0, 2] / np.sqrt(C[0, 0] * C[2, 2]),
                        "dE": resp[0], "dT": resp[2]}
    # detailed-balance control: symmetric J, isotropic noise -> response = C
    Js = np.array([[-1.0, 0.4, 0.0], [0.4, -1.0, 0.3], [0.0, 0.3, -1.0]])
    Cs = linalg.solve_continuous_lyapunov(Js, -2 * np.eye(3))
    rs = -np.linalg.solve(Js, np.eye(3))
    out["fdt_symmetric_err"] = float(np.abs(rs - Cs).max())
    # co-regulated module (shared driver X, positive coupling E -> X):
    Jp = np.array([[-1.0, 0.0, 0.6], [0.0, -1.0, 0.6], [0.5, 0.0, -1.0]])
    Cp = linalg.solve_continuous_lyapunov(Jp, -np.eye(3))
    rp = -np.linalg.solve(Jp, np.array([-1.0, 0.0, 0.0]))
    out["coreg"] = {"corr_ET": Cp[0, 1] / np.sqrt(Cp[0, 0] * Cp[1, 1]),
                    "dE": rp[0], "dT": rp[1]}
    return out


# ======================================================================= B
TOPP = {"si": 0.72, "Eg0": 1.44, "R0": 864.0, "sigma": 43.2, "alpha": 20000.0,
        "k": 432.0, "d0": 0.06, "r1": 0.00084, "r2": 2.4e-06}


class BetaPopulation:
    """N beta cells, cell i secretes s(G; theta_i) = G^n / (theta_i^n + G^n)
    (a steep single-cell switch, n = 8), thresholds theta ~ LogNormal(m, w).
    Macro secretion per unit beta mass: sigma_pop(G) = sigma_max * mean_i s.
    basal: glucose-independent secretion fraction (sulfonylurea-like)."""

    def __init__(self, m, w, n=8.0, sigma_max=TOPP["sigma"], N=600,
                 drop_low=0.0, basal=0.0):
        q = (np.arange(N) + 0.5) / N
        from scipy.stats import norm
        th = m * np.exp(w * norm.ppf(q))
        if drop_low > 0:
            th = th[int(drop_low * N):]
        self.th, self.n, self.smax, self.basal = th, n, sigma_max, basal

    def sigma(self, G):
        G = np.atleast_1d(np.asarray(G, float))[:, None]
        s = G ** self.n / (self.th ** self.n + G ** self.n)
        v = self.smax * (self.basal + (1 - self.basal) * s.mean(1))
        return v if v.size > 1 else float(v[0])

    def dsigma(self, G, h=1e-3):
        return (self.sigma(G + h) - self.sigma(G - h)) / (2 * h)

    def density(self, G):
        """threshold density p_theta(G) (per mg/dl)."""
        from scipy.stats import gaussian_kde
        return float(gaussian_kde(self.th)(G)[0])


def fit_population(n=8.0):
    from scipy.optimize import minimize
    Gs = np.linspace(40, 600, 200)
    target = TOPP["sigma"] * Gs ** 2 / (TOPP["alpha"] + Gs ** 2)

    def loss(p):
        P = BetaPopulation(np.exp(p[0]), np.exp(p[1]), n)
        return np.mean((P.sigma(Gs) - target) ** 2)
    r = minimize(loss, [np.log(140), np.log(0.8)], method="Nelder-Mead",
                 options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 2000})
    m, w = np.exp(r.x)
    P = BetaPopulation(m, w, n)
    rel = np.sqrt(r.fun) / TOPP["sigma"]
    return m, w, rel, Gs, target, P.sigma(Gs)


def topp_state(P, p=TOPP, beta=None):
    """Equilibria of the macro model with population secretion P.
    Healthy state G=100 and saddle G=250 are fixed by beta turnover
    (r1 G - r2 G^2 = d0); beta* follows from compensation."""
    disc = np.sqrt(p["r1"] ** 2 - 4 * p["r2"] * p["d0"])
    Gh, Gs = (p["r1"] - disc) / (2 * p["r2"]), (p["r1"] + disc) / (2 * p["r2"])
    Ih = (p["R0"] / Gh - p["Eg0"]) / p["si"]
    Is = (p["R0"] / Gs - p["Eg0"]) / p["si"]
    bh = p["k"] * Ih / P.sigma(Gh)
    bs = p["k"] * Is / P.sigma(Gs)
    # fast G-I Jacobian at the healthy state
    J = np.array([[-(p["Eg0"] + p["si"] * Ih), -p["si"] * Gh],
                  [bh * P.dsigma(Gh), -p["k"]]])
    ev = np.linalg.eigvals(J)
    zeta = -np.trace(J) / (2 * np.sqrt(np.linalg.det(J)))
    tau = ss.critical_delay(J, 0, 1)
    return {"G_h": Gh, "G_s": Gs, "beta_h": bh, "beta_s": bs,
            "margin": 1 - bs / bh, "J": J, "eig": ev, "zeta": zeta,
            "tau_crit_min": tau[0] * 24 * 60 if tau else np.inf,
            "gain": bh * P.dsigma(Gh)}


def fast_glucose(P, beta, p=TOPP):
    """Quasi-steady glucose of the fast G-I loop at a given beta mass:
    R0 = (Eg0 + si * beta * sigma(G)/k) * G (monotone in G)."""
    from scipy.optimize import brentq
    f = lambda G: (p["Eg0"] + p["si"] * beta * P.sigma(G) / p["k"]) * G - p["R0"]
    return brentq(f, 1e-3, 5000)


def simulate_topp(P, beta0, days=1500, p=TOPP, cap=None):
    """Slow beta-mass dynamics on the fast G-I manifold (time in days).
    cap: optional ceiling on beta mass (absolute units); net growth is
    throttled by (1 - beta/cap). Topp's original model has no ceiling, so
    compensation can demand any mass."""
    from scipy.integrate import solve_ivp

    def f(t, y):
        G = fast_glucose(P, max(y[0], 1e-6), p)
        rate = -p["d0"] + p["r1"] * G - p["r2"] * G ** 2
        if cap is not None and rate > 0:
            rate *= max(0.0, 1 - y[0] / cap)
        return [rate * y[0]]
    s = solve_ivp(f, (0, days), [beta0], max_step=5, rtol=1e-6)
    G = np.array([fast_glucose(P, max(b, 1e-6), p) for b in s.y[0]])
    return s.t, s.y[0], G


def fate(P, beta0, p=TOPP, cap_mult=3.0, days=4000):
    """Long-run attractor with a beta-mass ceiling of cap_mult x baseline:
    'healthy' (G back to the set point), 'compensated hyperglycaemia'
    (beta pinned at the ceiling, G between set point and saddle) or
    'collapse' (beta -> 0)."""
    t, b, G = simulate_topp(P, beta0, days, p, cap=cap_mult * beta0)
    st = topp_state(P, p)
    if b[-1] < 0.05 * beta0:
        kind = "collapse"
    elif G[-1] < st["G_h"] + 3:
        kind = "healthy"
    else:
        kind = "compensated hyperglycaemia"
    return {"kind": kind, "G": float(G[-1]), "beta": float(b[-1] / beta0)}


def part_B(cap_mult=3.0):
    out = {"cap_mult": cap_mult}
    m0, w0, rel, Gs, target, fit = fit_population()
    out["fit"] = {"m": m0, "w": w0, "rel_rmse": rel, "G": Gs, "topp": target,
                  "pop": fit}
    P0 = BetaPopulation(m0, w0)
    base = topp_state(P0)
    out["base"] = base
    beta0 = base["beta_h"]
    # susceptibility = density of cells at threshold (steep-switch limit)
    Gg = np.linspace(40, 400, 181)
    out["chi"] = {"G": Gg, "dsigma": P0.dsigma(Gg),
                  "density": np.array([P0.density(g) for g in Gg]) * TOPP["sigma"]}
    # B1 phase diagram over the micro distribution (m, w), with beta frozen
    # at the healthy individual's mass beta0 (the acute step), the slow beta
    # loop then decides: G_fast < G_s -> beta grows/settles back (converges);
    # G_fast > G_s -> beta shrinks, glucose rises further (diverges).
    ms = np.geomspace(0.7, 9.0, 61) * m0
    ws = np.linspace(0.15, 1.6, 59)
    Gfast = np.array([[fast_glucose(BetaPopulation(mm, ww), beta0) for mm in ms]
                      for ww in ws])
    out["phase"] = {"m": ms, "w": ws, "Gfast": Gfast, "G_s": base["G_s"]}
    from scipy.optimize import brentq
    crit = {}
    for lab, ww in (("synchronised", 0.3 * w0), ("normal", w0), ("broad", 2 * w0)):
        crit[lab] = brentq(lambda x: fast_glucose(BetaPopulation(x * m0, ww),
                                                  beta0) - base["G_s"], 1, 60)
    out["critical_shift"] = crit
    crit_drop = brentq(lambda f: fast_glucose(BetaPopulation(m0, w0, drop_low=f),
                                              beta0) - base["G_s"], 0.05, 0.95)
    out["critical_drop"] = crit_drop
    # B2 named micro scenarios
    sc = {
        "healthy": dict(m=m0, w=w0),
        "thresholds +40%": dict(m=1.4 * m0, w=w0),
        "thresholds +80%": dict(m=1.8 * m0, w=w0),
        "broadened (w x2)": dict(m=m0, w=2 * w0),
        "synchronised (w x0.3)": dict(m=m0, w=0.3 * w0),
        "lose 20% most sensitive": dict(m=m0, w=w0, drop_low=0.2),
        "lose 20% + thresholds +40%": dict(m=1.4 * m0, w=w0, drop_low=0.2),
    }
    rows = {}
    for k, v in sc.items():
        P = BetaPopulation(**v)
        st = topp_state(P)
        Gf = fast_glucose(P, beta0)
        rows[k] = {"G_fast": Gf, "converges": Gf < base["G_s"],
                   "fate": fate(P, beta0, cap_mult=cap_mult),
                   "beta_needed": st["beta_h"] / beta0, "margin": st["margin"],
                   "zeta": st["zeta"], "tau_crit_min": st["tau_crit_min"],
                   "gain": st["gain"] / base["gain"]}
    out["scenarios"] = rows
    # B3 levers on a diverging population: thresholds x2 and the most
    # sensitive cells lost, with the loss set so that the acute glucose sits
    # at 280 mg/dl (just past the saddle at 250)
    drop = brentq(lambda f: fast_glucose(BetaPopulation(2 * m0, w0, drop_low=f),
                                         beta0) - 280.0, 0.01, 0.95)
    sick = dict(m=2 * m0, w=w0, drop_low=drop)
    out["sick"] = sick
    levers = {
        "none": (dict(sick), {}),
        "GLP-1RA (thresholds x0.75)": (dict(sick, m=sick["m"] * 0.75), {}),
        "sulfonylurea (basal 15%)": (dict(sick, basal=0.15), {}),
        "insulin sensitiser (si x1.5)": (dict(sick), {"si": 1.5}),
        "metformin (R0 x0.8)": (dict(sick), {"R0": 0.8}),
        "beta protection (d0 x0.7)": (dict(sick), {"d0": 0.7}),
    }
    lev = {}
    for k, (pp, mult) in levers.items():
        p = dict(TOPP)
        for q, f in mult.items():
            p[q] *= f
        P = BetaPopulation(**pp)
        st = topp_state(P, p)
        Gf = fast_glucose(P, beta0, p)
        # hypoglycaemia risk: the same drug given to the healthy population
        Ph = BetaPopulation(**dict(pp, m=m0 * pp["m"] / sick["m"], drop_low=0.0))
        Gh_drug = fast_glucose(Ph, beta0, p)
        lev[k] = {"G_fast": Gf, "G_s": st["G_s"], "converges": Gf < st["G_s"],
                  "fate": fate(P, beta0, p, cap_mult=cap_mult),
                  "beta_needed": st["beta_h"] / beta0,
                  "margin_to_saddle": st["G_s"] - Gf, "G_h": st["G_h"],
                  "healthy_on_drug": Gh_drug}
    out["levers"] = lev
    # B4 time courses
    tc = {}
    for k, pp in {"thresholds +80% (converges)": dict(m=1.8 * m0, w=w0),
                  "x2 + cell loss (diverges)": dict(sick),
                  "same + GLP-1RA": dict(sick, m=sick["m"] * 0.75),
                  "same + sensitiser": dict(sick)}.items():
        p = dict(TOPP)
        if "sensitiser" in k:
            p["si"] *= 1.5
        t, b, G = simulate_topp(BetaPopulation(**pp), beta0, days=2500, p=p,
                                cap=cap_mult * beta0)
        tc[k] = (t, b / beta0, G)
    out["timecourse"] = tc
    # B5 heterogeneity trade-off: gain and critical delay vs width
    wsc = np.linspace(0.1, 1.5, 40)
    out["width_scan"] = {"w": wsc, "zeta": [], "tau": [], "gain": []}
    for ww in wsc:
        st = topp_state(BetaPopulation(m0, ww))
        out["width_scan"]["zeta"].append(st["zeta"])
        out["width_scan"]["tau"].append(st["tau_crit_min"])
        out["width_scan"]["gain"].append(st["gain"])
    return out


# ======================================================================= C
def _equal_n(rng, mask, n):
    return rng.choice(np.where(mask)[0], n, replace=False)


def part_C(n_rep=5, n_perm=4):
    from netjac.data import cptac, gprofiler as gp, omnipath as op, tcga
    rng = np.random.default_rng(0)
    E = op.signed_edges()
    FB = op.feedback_pairs(E)
    CX = op.complexes()
    fbset = {(a, b) for a, b, *_ in FB} | {(b, a) for a, b, *_ in FB}
    oneway = [(a, b, s) for a, b, s in E if (a, b) not in fbset]
    out = {"desc": {}, "modes": {}, "coupling": {}, "rewire": {}, "loadings": {},
           "energy": {}}
    data = {}
    for name in cptac.STUDIES:
        t0 = time.time()
        d = cptac.load(name)
        X = cptac.impute(d["X"])
        t = d["tumour"]
        g = d["genes"]
        data[name] = (X, t, g, d["cases"])
        n = int(min(t.sum(), (~t).sum()))
        # C1 macro descriptors at equal n
        desc = {"T": [], "N": []}
        for rep in range(n_rep):
            for lab, m in (("T", t), ("N", ~t)):
                idx = _equal_n(rng, m, n)
                c = sm.collective_modes(X[idx], n_null=3, rng=rng)
                desc[lab].append([c["PR"], c["k_null"], c["top1"], c["mean_abs_r"],
                                  X[idx].var(0, ddof=1).mean(), c["entropy"]])
        out["desc"][name] = {k: np.array(v) for k, v in desc.items()}
        out["desc"][name]["n"] = n
        out["desc"][name]["p"] = X.shape[1]
        # C2 mode identity (top 3 per ensemble, full ensemble)
        out["modes"][name] = {}
        for lab, m in (("T", t), ("N", ~t)):
            M = sm.modes(X[m], 10)
            out["loadings"][(name, lab)] = (g, M["loadings"])
            ann = []
            for k in range(3):
                L = M["loadings"][:, k]
                o = np.argsort(L)
                ann.append((M["lam"][k] / X.shape[1],
                            gp.top_term(g[o[-100:]], g)[0],
                            gp.top_term(g[o[:100]], g)[0]))
            out["modes"][name][lab] = ann
        # C4 prior-network couplings at equal n
        gs = set(g)
        cp = op.index_pairs(op.complex_pairs(CX, g), g)
        ow_l = [(a, b, s) for a, b, s in oneway if a in gs and b in gs]
        ow = op.index_pairs(ow_l, g)
        ows = np.array([s for *_, s in ow_l])
        fb_l = [x for x in FB if x[0] in gs and x[1] in gs]
        fbi = op.index_pairs(fb_l, g)
        fbk = np.array([x[4] for x in fb_l])
        rp = rng.integers(0, len(g), (100000, 2))
        rp = rp[rp[:, 0] != rp[:, 1]]
        cpl = {}
        for lab, m in (("T", t), ("N", ~t)):
            idx = _equal_n(rng, m, n)
            for mode, Z in (("raw", sm.zscore(X[idx])), ("res5", sm.remove_modes(X[idx], 5))):
                f = lambda P: sm.corr_pairs(Z, P)
                cpl[(lab, mode)] = {"complex": f(cp), "random": f(rp),
                                    "stim": f(ow[ows > 0]), "inh": f(ow[ows < 0]),
                                    "posFB": f(fbi[fbk == "positive"]),
                                    "negFB": f(fbi[fbk == "negative"])}
        out["coupling"][name] = cpl
        # C5 all-pairs rewiring, before / after subtracting macro fields
        rw = {}
        for k in (0, 1, 3, 5):
            r = sm.rewiring_scan(X[t], X[~t], thr=5.0, n_perm=n_perm,
                                 k_remove=k, keep=400000 if k == 1 else 50,
                                 rng=rng)
            rw[k] = {"n_obs": r["n_obs"], "null": r["null"], "n_pairs": r["n_pairs"],
                     "fdr": r["fdr"]}
            if k == 1:
                top = r["top"]
                rw["top_k1"] = [(g[int(a)], g[int(b)], dz, ra, rb)
                                for a, b, dz, ra, rb in top]
        out["rewire"][name] = rw
        print(f"  C {name}: {time.time() - t0:.0f}s", flush=True)
    # C3 universality of order parameters across indications
    common = None
    for (name, lab), (g, L) in out["loadings"].items():
        common = set(g) if common is None else common & set(g)
    common = sorted(common)
    keys = [(nm, lab) for lab in ("T", "N") for nm in cptac.STUDIES]
    V = {}
    for key in keys:
        g, _ = out["loadings"][key]
        X, t, gg, _ = data[key[0]]
        ix = {x: i for i, x in enumerate(gg)}
        cols = [ix[x] for x in common]
        m = t if key[1] == "T" else ~t
        V[key] = sm.modes(X[m][:, cols], 10)["loadings"]
    Ov = np.array([[sm.subspace_overlap(V[a], V[b]) for b in keys] for a in keys])
    out["universality"] = {"keys": keys, "overlap": Ov, "n_common": len(common),
                           "random": 10 / len(common)}
    # C5b recurrence of residual (k=1) rewired pairs across indications
    cnt = {}
    for name in cptac.STUDIES:
        for a, b, dz, ra, rb in out["rewire"][name]["top_k1"]:
            key = (a, b) if a < b else (b, a)
            cnt.setdefault(key, []).append((name, np.sign(dz), ra, rb))
    rec = [(k, v) for k, v in cnt.items()
           if len(v) >= 3 and len({s for _, s, *_ in v}) == 1]
    rec.sort(key=lambda kv: -len(kv[1]))
    # null recurrence: shuffle gene labels within each indication
    out["recurrent"] = rec[:200]
    out["n_recurrent"] = len(rec)
    null_rec = []
    for rep in range(3):
        cntn = {}
        for name in cptac.STUDIES:
            g = data[name][2]
            perm = dict(zip(g, rng.permutation(g)))
            for a, b, dz, *_ in out["rewire"][name]["top_k1"]:
                a, b = perm[a], perm[b]
                key = (a, b) if a < b else (b, a)
                cntn.setdefault(key, []).append(np.sign(dz))
        null_rec.append(sum(1 for v in cntn.values() if len(v) >= 3 and len(set(v)) == 1))
    out["n_recurrent_null"] = null_rec
    # C6 Boltzmann energy relative to the normal ensemble vs grade / stage
    cl = tcga.clinical("CPTAC-3").set_index("patient")

    def num(x, table):
        if not isinstance(x, str):
            return np.nan
        x = x.replace("Stage ", "")
        for k, v in table:
            if x.startswith(k):
                return v
        return np.nan
    GT = [("G1", 1), ("G2", 2), ("G3", 3), ("G4", 4)]
    ST = [("IV", 4), ("III", 3), ("II", 2), ("I", 1)]
    for name in ("ccRCC", "LUAD", "LSCC", "HNSCC", "PDAC", "UCEC"):
        X, t, g, cases = data[name]
        N = X[~t]
        mu, sd = N.mean(0), N.std(0, ddof=1)
        sd[sd == 0] = 1
        Zn, Zt = (N - mu) / sd, (X[t] - mu) / sd
        sp = sm.ShrinkPrecision(Zn)
        energy = ((Zt ** 2).sum(1) - ((Zt @ sp.W) ** 2).sum(1)) / sp.c / 2 / X.shape[1]
        euclid = np.sqrt((Zt ** 2).mean(1))
        Vt = np.linalg.svd(Zn, full_matrices=False)[2][:10]
        soft = ((Zt @ Vt.T) ** 2).sum(1) / (Zt ** 2).sum(1)
        ct = cases[t]
        res = {}
        for lab, col, tab in (("grade", "grade", GT), ("stage", "stage", ST)):
            v = np.array([num(cl[col].get(c), tab) if c in cl.index else np.nan
                          for c in ct])
            m = np.isfinite(v)
            if m.sum() > 20 and len(np.unique(v[m])) > 1:
                res[lab] = {k: spearmanr(x[m], v[m]) for k, x in
                            (("energy", energy), ("euclid", euclid),
                             ("soft", soft))}
                res[lab]["n"] = int(m.sum())
        res["r_energy_euclid"] = float(np.corrcoef(energy, euclid)[0, 1])
        out["energy"][name] = res
    del out["loadings"]
    return out


# ======================================================================= D
FAMILIES = {
    "proteasome": r"^PSM[ABCD]\d+$",
    "tRNA synthetase": r"^([A-Z]ARS|EPRS|FARS[AB])$",
    "sterol synthesis": r"^(HMGCS1|HMGCR|MVK|PMVK|MVD|FDPS|FDFT1|SQLE|LSS|CYP51A1|"
                        r"IDI1|DHCR7|SC5D|MSMO1|NSDHL|HSD17B7|EBP|TM7SF2|DHCR24|SCAP|SREBF2)$",
    "ER chaperone / UPR": r"^(HSPA5|HSP90B1|PDIA\d|DNAJB11|CALR|CANX|HYOU1|SEL1L|"
                          r"DERL\d|SYVN1|DDOST|SEC61A1|SEC63)$",
    "actin (MRTF-SRF)": r"^(ACTB|ACTG1|CFL1|PFN1|CAPZB|CAPZA1|DSTN)$",
    "cytosolic ribosome": r"^RP[SL]\d+[A-Z]?\d*$",
}


def part_D(K=50, lasso_alpha=0.01):
    from sklearn.linear_model import Lasso
    from netjac.data import perturbseq as ps, omnipath as op, gprofiler as gp
    out = {}
    r = ps.reduce_singlecell()
    H, Nh, pert, gid, gene = (r["half_means"], r["half_n"], r["pert"],
                              r["gene_id"], r["gene"])
    pg = np.array([p.split("_")[1] for p in pert])
    pid = np.array([p.split("_")[-1] for p in pert])
    ix = {e: i for i, e in enumerate(gid)}
    ctrl = pg == "non-targeting"
    meas = np.array([i for i in range(len(pert)) if pid[i] in ix and not ctrl[i]])
    own = np.array([ix[pid[i]] for i in meas])
    n = len(meas)
    out["n_perturbed_measured"] = n
    out["n_ctrl_cells"] = len(r["ctrl_cells"])
    # D1 equilibrium ensemble: control cells -> covariance columns
    C = r["ctrl_cells"].astype(np.float64)
    C = (C - C.mean(0)) / C.std(0)
    Ccol = (C.T @ C[:, own]) / len(C)
    lam = np.linalg.eigvalsh(np.cov(C[:, np.argsort(-C.var(0))[:2000]].T))[::-1]
    out["ctrl_PR"] = sm.participation_ratio(lam)
    del C
    # D2 FDT ranking: response of gene j's own top-K co-fluctuating partners,
    # ranked against the partner sets of every other perturbed gene
    P = np.zeros((n, len(gid)), np.float32)
    for k in range(n):
        c = Ccol[:, k].copy()
        c[own[k]] = -np.inf
        P[k, np.argsort(-c)[:K]] = 1.0 / K

    def ranks(h, assign):
        Dm = H[meas, h].copy()
        Dm[np.arange(n), own] = 0
        M = Dm @ P.T
        return np.array([(M[k] < M[k, assign[k]]).mean() for k in range(n)])
    idn = np.arange(n)
    qA, qB = ranks(0, idn), ranks(1, idn)
    rng = np.random.default_rng(1)
    sh = rng.permutation(n)
    nA, nB = ranks(0, sh), ranks(1, sh)
    g = pg[meas]
    cls = {}
    for tt in (0.05, 0.02):
        cls[tt] = {"fdt": int(np.sum((qA < tt) & (qB < tt))),
                   "fdt_null": int(np.sum((nA < tt) & (nB < tt))),
                   "anti": int(np.sum((qA > 1 - tt) & (qB > 1 - tt))),
                   "anti_null": int(np.sum((nA > 1 - tt) & (nB > 1 - tt)))}
    out["fdt_counts"] = cls
    anti = (qA > 0.95) & (qB > 0.95)
    fdt = (qA < 0.05) & (qB < 0.05)
    out["anti_genes"] = sorted(set(g[anti]))
    out["fdt_genes"] = sorted(set(g[fdt]))
    fam = {}
    for k, pat in FAMILIES.items():
        mk = np.array([bool(re.match(pat, x)) for x in g])
        a, f, b = int(np.sum(anti & mk)), int(np.sum(fdt & mk)), int(mk.sum())
        pa = fisher_exact([[a, b - a], [anti.sum() - a, (~mk).sum() - anti.sum() + a]],
                          alternative="greater")[1]
        pf = fisher_exact([[f, b - f], [fdt.sum() - f, (~mk).sum() - fdt.sum() + f]],
                          alternative="greater")[1]
        fam[k] = {"tested": b, "anti": a, "p_anti": pa, "fdt": f, "p_fdt": pf}
    out["families"] = fam
    out["fdt_q"] = (qA, qB, g)
    # worked examples
    ex = {}
    for name in ("HMGCS1", "HSPA5", "PSMB6", "RPS3", "PHB"):
        w = np.where(g == name)[0]
        if not len(w):
            continue
        k = w[0]
        d = H[meas[k]].mean(0).copy()
        d[own[k]] = np.nan
        c = Ccol[:, k].copy()
        c[own[k]] = -np.inf
        top = np.argsort(-c)[:K]
        ex[name] = {"self": float(H[meas[k]].mean(0)[own[k]]),
                    "partners": [(gene[i], float(Ccol[i, k]), float(d[i])) for i in top[:12]],
                    "partner_mean": float(np.nanmean(d[top])),
                    "all_mean": float(np.nanmean(d)), "qA": qA[k], "qB": qB[k],
                    "scatter": (Ccol[:, k].copy(), d)}
    out["examples"] = ex
    # D3 sparse modular response analysis on knocked-down genes
    D = H[meas].mean(1)
    selfz = D[np.arange(n), own]
    best = {}
    for k in range(n):
        if g[k] not in best or selfz[k] < selfz[best[g[k]]]:
            best[g[k]] = k
    ks = np.array(sorted(best.values()))
    rep = np.array([np.corrcoef(H[meas[k], 0], H[meas[k], 1])[0, 1] for k in ks])
    sel = ks[(selfz[ks] < -0.5) & (rep > 0.2)]
    oS = own[sel]
    gS = g[sel]
    nS = len(sel)

    def RSS(Hh):
        M = Hh[meas[sel]][:, oS].T
        return M / np.diag(M)[None, :]
    RA, RB, Rm = RSS(H[:, 0]), RSS(H[:, 1]), RSS(H.mean(1))
    off = ~np.eye(nS, dtype=bool)

    def sparse_mra(R):
        o = np.zeros((nS, nS))
        for i in range(nS):
            cols = np.arange(nS) != i
            Xk = R[:, cols].T.copy()
            Xk[:, i] = 0
            o[i] = Lasso(alpha=lasso_alpha, fit_intercept=False,
                         max_iter=2000).fit(Xk, R[i, cols]).coef_
        return o
    sA, sB = sparse_mra(RA), sparse_mra(RB)
    repE = (sA != 0) & (sB != 0) & (np.sign(sA) == np.sign(sB)) & off
    nzA, nzB = ((sA != 0) & off).sum(), ((sB != 0) & off).sum()
    rloc = np.where(repE, (sA + sB) / 2, 0)
    gi = {x: i for i, x in enumerate(gS)}
    known_dir = np.zeros((nS, nS), bool)
    for a, b, s in op.signed_edges():
        if a in gi and b in gi:
            known_dir[gi[b], gi[a]] = True
    cpx = np.zeros((nS, nS), bool)
    for a, b in op.complex_pairs(op.complexes(), gS):
        cpx[gi[a], gi[b]] = cpx[gi[b], gi[a]] = True
    known = known_dir | known_dir.T | cpx
    Kn = int(repE.sum())
    topR = np.argsort(-np.abs(Rm[off]))[:Kn]
    topS = np.argsort(-(np.abs(Rm) + np.abs(Rm.T))[off])[:Kn]
    # MRA and raw responses (both halves), off-diagonal agreement
    out["mra"] = {"n_nodes": nS, "R_split_half_r": float(np.corrcoef(RA[off], RB[off])[0, 1]),
                  "edges_A": int(nzA), "edges_B": int(nzB), "edges_both": Kn,
                  "edges_chance": float(nzA * nzB / off.sum()),
                  "sign_agree": float(np.mean(np.sign(sA[repE]) == np.sign(sB[repE]))),
                  "prec_mra": float(known[off][repE[off]].mean()),
                  "prec_R": float(known[off][topR].mean()),
                  "prec_Rsym": float(known[off][topS].mean()),
                  "base": float(known[off].mean()),
                  "frac_activating": float(np.mean(rloc[repE] > 0))}
    # naive inversion (ridge) for contrast
    rnA, rnB = sm.mra_local(RA), sm.mra_local(RB)
    out["mra"]["naive_inverse_split_half_r"] = float(np.corrcoef(rnA[off], rnB[off])[0, 1])
    # mutual loops in the local (MRA) network
    mut = np.argwhere(np.triu(repE & repE.T, 1))
    kinds = {"mutual activation": 0, "mutual inhibition": 0, "negative": 0}
    for i, j in mut:
        s = np.sign(rloc[i, j]) * np.sign(rloc[j, i])
        kinds["negative" if s < 0 else ("mutual activation" if rloc[i, j] > 0
                                        else "mutual inhibition")] += 1
    out["mra"]["loops_local"] = kinds
    out["mra"]["loops_local_chance"] = float((repE[off].mean()) ** 2 * off.sum() / 2)
    # local Jacobian J = -I + r (uniform turnover) -> spectrum
    Jl = -np.eye(nS) + rloc
    ev = np.linalg.eigvals(Jl)
    out["mra"]["eig"] = ev
    out["mra"]["out_degree"] = sorted(zip(repE.sum(0), gS), reverse=True)[:12]
    # D4 loops from TOTAL responses (asymmetric information), split-half
    sel2 = ks[selfz[ks] < -0.5]
    o2 = own[sel2]
    g2 = g[sel2]
    A_ = H[meas[sel2]][:, 0][:, o2].T * np.sqrt(Nh[meas[sel2], 0])[None, :]
    B_ = H[meas[sel2]][:, 1][:, o2].T * np.sqrt(Nh[meas[sel2], 1])[None, :]
    off2 = ~np.eye(len(sel2), dtype=bool)
    thr = 4.0
    sig = (np.abs(A_) > thr) & (np.abs(B_) > thr) & (np.sign(A_) == np.sign(B_)) & off2
    disc = (np.abs(A_) > thr) & (np.abs(B_) > thr) & (np.sign(A_) != np.sign(B_)) & off2
    eff = -np.sign(A_ + B_)          # + : knocking j lowers i -> j activates i
    pairs = np.argwhere(np.triu(sig & sig.T, 1))
    kinds2 = {"mutual activation": [], "mutual inhibition": [], "negative": []}
    for i, j in pairs:
        s = eff[i, j] * eff[j, i]
        if s < 0:
            src_act = (g2[j], g2[i]) if eff[i, j] > 0 else (g2[i], g2[j])
            kinds2["negative"].append(src_act + (min(abs(A_[i, j] + B_[i, j]),
                                                     abs(A_[j, i] + B_[j, i])) / 2,))
        else:
            kinds2["mutual activation" if eff[i, j] > 0 else "mutual inhibition"].append(
                (g2[i], g2[j]))
    neg = sorted(kinds2["negative"], key=lambda x: -x[2])
    out["total_loops"] = {"n_nodes": len(sel2), "sig_edges": int(sig.sum()),
                          "discordant": int(disc.sum()),
                          "counts": {k: len(v) for k, v in kinds2.items()},
                          "chance": float(sig[off2].mean() ** 2 * off2.sum() / 2),
                          "neg_top": neg[:25]}
    # module-level summary of negative loops: activator side / repressed side
    act = {}
    for a, b, _ in neg:
        act.setdefault(a, set()).add(b)
    out["total_loops"]["neg_hubs"] = sorted(((len(v), k) for k, v in act.items()),
                                            reverse=True)[:12]
    # D5 macro response modes + universality K562 vs RPE1
    Bk = ps.load_bulk("K562_essential_bulk")
    strong = (Bk["obs"]["energy_test_p_value"].values < 0.01) & ~Bk["control"]
    Rr = Bk["X"][strong]
    U, s, Vt = np.linalg.svd(Rr - Rr.mean(0), full_matrices=False)
    lam = s ** 2
    Rn = sm.shuffle_columns(Rr, np.random.default_rng(0))
    s2 = np.linalg.svd(Rn - Rn.mean(0), compute_uv=False) ** 2
    out["resp_modes"] = {"n_pert": int(strong.sum()), "PR": sm.participation_ratio(lam),
                         "top10": float(lam[:10].sum() / lam.sum()),
                         "k_null": int((lam > s2[0]).sum()), "lam": lam[:40] / lam.sum(),
                         "ann": [(lam[k] / lam.sum(),
                                  gp.top_term(Bk["gene"][np.argsort(Vt[k])[-100:]], Bk["gene"])[0],
                                  gp.top_term(Bk["gene"][np.argsort(Vt[k])[:100]], Bk["gene"])[0])
                                 for k in range(5)]}
    Rp = ps.load_bulk("rpe1_bulk")
    sp_ = (Rp["obs"]["energy_test_p_value"].values < 0.01) & ~Rp["control"]
    cg = np.intersect1d(Bk["gene_id"], Rp["gene_id"])
    cpn = np.intersect1d(Bk["pert_gene"][strong], Rp["pert_gene"][sp_])

    def mat(Dd):
        gi_ = {x: i for i, x in enumerate(Dd["gene_id"])}
        cols = [gi_[x] for x in cg]
        return np.array([Dd["X"][Dd["pert_gene"] == p][:, cols].mean(0) for p in cpn])
    A2, B2 = mat(Bk), mat(Rp)
    pc = np.array([np.corrcoef(A2[i], B2[i])[0, 1] for i in range(len(cpn))])
    shp = np.random.default_rng(0).permutation(len(cpn))
    pn = np.array([np.corrcoef(A2[i], B2[shp[i]])[0, 1] for i in range(len(cpn))])
    uni = {}
    for k in (5, 10, 20):
        Va = np.linalg.svd(A2 - A2.mean(0), full_matrices=False)[2][:k].T
        Vb = np.linalg.svd(B2 - B2.mean(0), full_matrices=False)[2][:k].T
        uni[k] = (sm.subspace_overlap(Va, Vb), k / len(cg))
    out["universality"] = {"n_genes": len(cg), "n_pert": len(cpn),
                           "per_pert_r": pc, "per_pert_null": pn, "overlap": uni}
    return out


# ======================================================================= E
def part_E(n_cells=10, reps=20):
    from netjac.data import islet
    ID = ["INS", "IAPP", "G6PC2", "ABCC8", "MAFA", "PDX1", "PCSK1", "SLC30A8",
          "HADH", "NPTX2", "ERO1B", "DLK1", "PCSK2", "CHGB", "SCG5", "GLP1R",
          "SLC2A2", "GCK", "KCNJ11", "NKX6-1", "PAX6", "NEUROD1", "ISL1", "MAFB"]
    out = {}
    rng = np.random.default_rng(0)
    for nm, loader in (("Xin 2016", islet.load_xin), ("Segerstolpe 2016", islet.load)):
        d = loader()
        X = d["X"]
        g = list(d["genes"])
        expr = (X > 0).mean(0) > 0.5
        idg = [g.index(x) for x in ID if x in g and expr[g.index(x)]]
        rows = []
        for don in np.unique(d["donor"]):
            m = d["donor"] == don
            if m.sum() < n_cells:
                continue
            vals = []
            for _ in range(reps):
                idx = rng.choice(np.where(m)[0], n_cells, replace=False)
                lam = sm.corr_eigs(X[idx][:, expr])
                with np.errstate(invalid="ignore", divide="ignore"):
                    Cm = np.corrcoef(X[idx][:, idg].T)
                iu = np.triu_indices(len(idg), 1)
                vals.append((sm.participation_ratio(lam), np.nanmean(Cm[iu]),
                             X[idx][:, expr].var(0, ddof=1).mean()))
            rows.append((bool(d["t2d"][m][0]), *np.nanmean(vals, 0), int(m.sum())))
        rows = np.array(rows, dtype=float)
        t = rows[:, 0].astype(bool)
        res = {"n_ND": int((~t).sum()), "n_T2D": int(t.sum()),
               "cells": rows[:, 4], "gck_dropout": float((X[:, g.index("GCK")] == 0).mean())}
        for j, lab in enumerate(["PR", "identity-module r", "temperature"]):
            a, b = rows[~t, 1 + j], rows[t, 1 + j]
            res[lab] = (a.mean(), b.mean(), mannwhitneyu(a, b).pvalue)
        out[nm] = res
    return out


# ================================================================== report
def report(R):
    if "A" in R:
        A = R["A"]
        print("\nA  Theory")
        print(f"   2-node closed form vs Lyapunov: max err {A['closed_form_err']:.1e}")
        print("   median |r| for random stable loops: " + ", ".join(
            f"{k} {np.median(v):.2f}" for k, v in A["loop_theory"].items()))
        h, c = A["homeostat"], A["coreg"]
        print(f"   hidden-sensor homeostat: corr(E,T) = {h['corr_ET']:+.2f}, "
              f"knock E: dE {h['dE']:+.2f}, dT {h['dT']:+.2f}  (anti-FDT)")
        print(f"   co-regulated module: corr {c['corr_ET']:+.2f}, dE {c['dE']:+.2f}, "
              f"dT {c['dT']:+.2f}  (FDT-consistent)")
        print(f"   symmetric J + isotropic noise: |-J^-1 - C| = {A['fdt_symmetric_err']:.1e}")
    if "B" in R:
        B = R["B"]
        f, b = B["fit"], B["base"]
        print("\nB  Beta-cell population -> insulin-glucose Jacobian (Topp 2000)")
        print(f"   steep switches (n=8), thresholds LogNormal(median {f['m']:.0f} mg/dl, "
              f"width {f['w']:.2f}) reproduce Topp's G^2/(a+G^2): rel. RMSE {f['rel_rmse']:.3%}")
        print(f"   healthy: beta* {b['beta_h']:.0f}, saddle beta {b['beta_s']:.1f} -> "
              f"beta-loss margin {b['margin']:.1%}; fast-loop zeta {b['zeta']:.2f}, "
              f"critical insulin-action delay {b['tau_crit_min']:.0f} min")
        print("   acute threshold shift that tips the system (G_fast > 250): " + ", ".join(
            f"{k} x{v:.1f}" for k, v in B["critical_shift"].items()))
        print(f"   loss of most-sensitive cells that tips it: {B['critical_drop']:.0%}")
        for k, v in B["scenarios"].items():
            print(f"   {k:28s} G_fast {v['G_fast']:5.0f} {'converges' if v['converges'] else 'DIVERGES'}"
                  f" | beta needed x{v['beta_needed']:.2f} | gain x{v['gain']:.2f} | "
                  f"zeta {v['zeta']:.2f} | tau* {v['tau_crit_min']:.0f} min | with cap "
                  f"x{B['cap_mult']:.0f}: {v['fate']['kind']} (G {v['fate']['G']:.0f})")
        print(f"   diverging case: thresholds x2 + lose {B['sick']['drop_low']:.0%} most sensitive")
        for k, v in B["levers"].items():
            print(f"   {k:30s} G_fast {v['G_fast']:5.0f} vs saddle {v['G_s']:.0f}: "
                  f"{'CONVERGES' if v['converges'] else 'diverges'} (Topp, beta needed "
                  f"x{v['beta_needed']:.1f}); with cap: {v['fate']['kind']} (G {v['fate']['G']:.0f}); "
                  f"healthy person on it: {v['healthy_on_drug']:.0f} mg/dl")
        for k, (t, bb, G) in B["timecourse"].items():
            print(f"   time course {k:30s}: beta {bb[-1]:.2f} x, G {G[-1]:.0f} at {t[-1]:.0f} d")
    if "C" in R:
        C = R["C"]
        print("\nC  CPTAC tumour vs normal (equal n per ensemble, 5 subsamples)")
        lab = ["PR", "k>null", "top1", "mean|r|", "var", "entropy"]
        for nm, d in C["desc"].items():
            T, N = d["T"].mean(0), d["N"].mean(0)
            print(f"   {nm:6s} n={d['n']:3d} p={d['p']:5d} | " + " | ".join(
                f"{l} {a:.3g}/{b:.3g}" for l, a, b in zip(lab, T, N)))
        pr = np.array([[d["T"][:, 0].mean(), d["N"][:, 0].mean()] for d in C["desc"].values()])
        mr = np.array([[d["T"][:, 3].mean(), d["N"][:, 3].mean()] for d in C["desc"].values()])
        print(f"   PR higher in tumour in {np.sum(pr[:, 0] > pr[:, 1])}/8; mean|r| lower in "
              f"{np.sum(mr[:, 0] < mr[:, 1])}/8")
        for nm, m in C["modes"].items():
            for labx in ("N", "T"):
                v, hi, lo = m[labx][0]
                print(f"   {nm:6s} {labx} mode1 {v:.0%}: + {hi[:45]} | - {lo[:45]}")
        U = C["universality"]
        Ov, keys = U["overlap"], U["keys"]
        nT = len(keys) // 2
        tt = Ov[:nT, :nT][np.triu_indices(nT, 1)]
        nn = Ov[nT:, nT:][np.triu_indices(nT, 1)]
        tn = Ov[:nT, nT:][~np.eye(nT, dtype=bool)]
        print(f"   mode-subspace overlap (k=10, {U['n_common']} common proteins, random "
              f"{U['random']:.4f}): T-T {tt.mean():.3f}, N-N {nn.mean():.3f}, T-N(other) {tn.mean():.3f}, "
              f"same-indication T-N {np.mean(np.diag(Ov[:nT, nT:])):.3f}")
        print("   couplings (tumour, raw | res5): complex median r, stim vs inh mean r, "
              "|r| posFB vs negFB (MWU p)")
        pos_gt = 0
        for nm, cp in C["coupling"].items():
            s = []
            for lab2 in ("T", "N"):
                for mode in ("raw", "res5"):
                    c = cp[(lab2, mode)]
                    p = mannwhitneyu(np.abs(c["posFB"]), np.abs(c["negFB"]),
                                     alternative="greater").pvalue
                    s.append(f"{lab2}-{mode}: cx {np.median(c['complex']):.3f} "
                             f"st/inh {c['stim'].mean():.3f}/{c['inh'].mean():.3f} "
                             f"FB {np.median(np.abs(c['posFB'])):.3f}/{np.median(np.abs(c['negFB'])):.3f} p={p:.0e}")
                    if lab2 == "T" and mode == "res5":
                        pos_gt += np.median(np.abs(c["posFB"])) > np.median(np.abs(c["negFB"]))
            print(f"   {nm:6s} " + " ; ".join(s))
        print(f"   tumour res5: |r| posFB > negFB in {pos_gt}/8 indications")
        print("   all-pairs rewiring |dz|>5 (pairs, permutation null mean, FDR) by modes removed:")
        for nm, rw in C["rewire"].items():
            print(f"   {nm:6s} {rw[0]['n_pairs']:.2e} pairs | " + " | ".join(
                f"k={k}: {rw[k]['n_obs']} (null {np.mean(rw[k]['null']):.0f}, FDR {rw[k]['fdr']:.2f})"
                for k in (0, 1, 3, 5)))
        print(f"   pairs rewired (k=1) in >=3 indications with one sign: {C['n_recurrent']} "
              f"(gene-label null {C['n_recurrent_null']})")
        for (a, b), v in C["recurrent"][:15]:
            print(f"      {a}-{b}: {len(v)} indications ({', '.join(x[0] for x in v)}), "
                  f"{'gain' if v[0][1] > 0 else 'loss'} of coupling in tumour")
        print("   Boltzmann energy vs grade/stage (Spearman rho, p):")
        for nm, e in C["energy"].items():
            s = []
            for lab2 in ("grade", "stage"):
                if lab2 in e:
                    s.append(f"{lab2}(n={e[lab2]['n']}): " + " ".join(
                        f"{k} {e[lab2][k][0]:+.2f}({e[lab2][k][1]:.0e})"
                        for k in ("energy", "euclid", "soft")))
            print(f"   {nm:6s} r(E,euclid)={e['r_energy_euclid']:.2f} " + " ; ".join(s))
    if "D" in R:
        D = R["D"]
        print("\nD  Perturb-seq K562 (Replogle 2022): directed Jacobian from knockdowns")
        print(f"   {D['n_perturbed_measured']} knockdowns of measured genes, "
              f"{D['n_ctrl_cells']} control cells (equilibrium ensemble, PR {D['ctrl_PR']:.1f})")
        for tt, c in D["fdt_counts"].items():
            print(f"   both halves q<{tt}: FDT-consistent {c['fdt']} (null {c['fdt_null']}), "
                  f"anti-FDT {c['anti']} (null {c['anti_null']})")
        for k, v in D["families"].items():
            print(f"   {k:20s} tested {v['tested']:3d}: anti-FDT {v['anti']:2d} (p={v['p_anti']:.0e}), "
                  f"FDT {v['fdt']:2d} (p={v['p_fdt']:.0e})")
        for k, v in D["examples"].items():
            print(f"   {k}: self {v['self']:+.2f}; top-50 co-fluctuating partners mean change "
                  f"{v['partner_mean']:+.2f} (all genes {v['all_mean']:+.3f}); rank q {v['qA']:.2f}/{v['qB']:.2f}; "
                  + ", ".join(f"{a}({b:+.2f}->{c:+.2f})" for a, b, c in v["partners"][:6]))
        m = D["mra"]
        print(f"   sparse MRA on {m['n_nodes']} nodes: R split-half r {m['R_split_half_r']:.2f}; "
              f"naive inverse split-half r {m['naive_inverse_split_half_r']:.3f}; "
              f"lasso edges A {m['edges_A']} B {m['edges_B']} both {m['edges_both']} "
              f"(chance {m['edges_chance']:.0f}), sign agreement {m['sign_agree']:.3f}")
        print(f"   precision for known interactions: MRA {m['prec_mra']:.3f} | total |R| {m['prec_R']:.3f}"
              f" | |R|+|R^T| {m['prec_Rsym']:.3f} | base {m['base']:.3f}")
        print(f"   local loops {m['loops_local']} (chance {m['loops_local_chance']:.0f}); "
              f"Jacobian -I+r leading Re(lambda) {np.max(m['eig'].real):.3f}")
        tl = D["total_loops"]
        print(f"   total-response loops ({tl['n_nodes']} nodes, |z|>4 both halves, "
              f"discordant {tl['discordant']}): {tl['counts']} (chance {tl['chance']:.0f})")
        print("   negative-loop hubs (activator side, #partners): " + ", ".join(
            f"{k}({n})" for n, k in tl["neg_hubs"]))
        rm = D["resp_modes"]
        print(f"   response macro modes: {rm['n_pert']} strong knockdowns -> PR {rm['PR']:.1f}, "
              f"{rm['k_null']} modes above shuffle null, top-10 {rm['top10']:.0%}")
        for v, hi, lo in rm["ann"]:
            print(f"      {v:.1%}: + {hi[:50]} | - {lo[:50]}")
        u = D["universality"]
        print(f"   K562 vs RPE1 ({u['n_pert']} shared knockdowns, {u['n_genes']} genes): "
              f"per-knockdown response r median {np.median(u['per_pert_r']):.2f} "
              f"(mismatched {np.median(u['per_pert_null']):.2f}); subspace overlap " + ", ".join(
                  f"k={k}: {a:.2f} (random {b:.4f})" for k, (a, b) in u["overlap"].items()))
    if "E" in R:
        print("\nE  Islet beta-cell populations, donor level (10-cell subsamples)")
        for nm, e in R["E"].items():
            print(f"   {nm}: ND {e['n_ND']} vs T2D {e['n_T2D']} donors; GCK dropout "
                  f"{e['gck_dropout']:.0%}; " + "; ".join(
                      f"{k} {a:.3f}/{b:.3f} p={p:.2f}" for k, (a, b, p) in
                      ((k, e[k]) for k in ("PR", "identity-module r", "temperature"))))


# ================================================================= figures
def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    paths = []
    # ---------------------------------------------------- d9a: theory + beta
    if "A" in R and "B" in R:
        A, B = R["A"], R["B"]
        fig, ax = plt.subplots(2, 4, figsize=(19, 8.6))
        a = ax[0, 0]
        for k, c in (("positive", "C3"), ("one-way", "C7"), ("negative", "C0")):
            a.hist(A["loop_theory"][k], bins=50, range=(0, 1), histtype="step",
                   lw=2, color=c, density=True, label=f"{k} (median {np.median(A['loop_theory'][k]):.2f})")
        a.set_xlabel("|stationary correlation|")
        a.set_title("A  2-node loops: negative feedback hides\nitself in static correlations")
        a.legend(fontsize=8)
        a = ax[0, 1]
        h, c = A["homeostat"], A["coreg"]
        xs = np.arange(2)
        a.bar(xs - 0.18, [h["corr_ET"], c["corr_ET"]], 0.36, label="co-fluctuation corr(E,T)")
        a.bar(xs + 0.18, [h["dT"], c["dT"]], 0.36, label="response of T to E knockdown")
        a.axhline(0, color="k", lw=0.8)
        a.set_xticks(xs, ["hidden-sensor\nhomeostat", "co-regulated\nmodule"])
        a.set_title("A  fluctuation-response sign test:\nanti-FDT = negative feedback")
        a.legend(fontsize=8)
        f = B["fit"]
        a = ax[0, 2]
        P0 = BetaPopulation(f["m"], f["w"])
        for th in P0.th[::60]:
            a.plot(f["G"], TOPP["sigma"] * f["G"] ** 8 / (th ** 8 + f["G"] ** 8),
                   color="0.8", lw=0.8)
        a.plot(f["G"], f["topp"], "k", lw=3, label="Topp macro G²/(α+G²)")
        a.plot(f["G"], f["pop"], "C3--", lw=2,
               label=f"mean of 600 steep cells\n(RMSE {f['rel_rmse']:.1%})")
        a.set_xlabel("glucose (mg/dl)")
        a.set_ylabel("secretion per β mass")
        a.set_title("B  micro switches → macro dose-response")
        a.legend(fontsize=8)
        a = ax[0, 3]
        ch = B["chi"]
        a.plot(ch["G"], ch["dsigma"], "k", lw=2, label="macro gain dσ/dG")
        a.plot(ch["G"], ch["density"], "C3--", lw=2, label="σmax × threshold density")
        a.axvline(100, color="C2", ls=":", label="fasting set point")
        a.set_xlabel("glucose (mg/dl)")
        a.set_title("B  susceptibility = density of cells\nat threshold (loop gain J21)")
        a.legend(fontsize=8)
        a = ax[1, 0]
        ph = B["phase"]
        cs = a.contourf(ph["m"] / f["m"], ph["w"], np.minimum(ph["Gfast"], 600),
                        levels=np.linspace(80, 600, 27), cmap="magma_r")
        a.contour(ph["m"] / f["m"], ph["w"], ph["Gfast"], levels=[ph["G_s"]],
                  colors="cyan", linewidths=2.5)
        a.plot(1, f["w"], "w*", ms=14)
        a.set_xscale("log")
        a.set_xticks([1, 2, 3, 5, 8], ["1", "2", "3", "5", "8"])
        a.minorticks_off()
        a.set_xlabel("threshold shift (× healthy median)")
        a.set_ylabel("threshold heterogeneity w")
        a.set_title("B  converge | diverge boundary (cyan)\nacute glucose at fixed β mass")
        plt.colorbar(cs, ax=a, label="G_fast (mg/dl)")
        a = ax[1, 1]
        for k, (t, bb, G) in B["timecourse"].items():
            a.plot(t, G, lw=2, label=k)
        a.axhline(250, color="k", ls=":", lw=1)
        a.set_ylim(0, 620)
        a.set_xlabel("days")
        a.set_ylabel("glucose (mg/dl)")
        a.set_title("B  slow β-mass loop: compensation\nvs collapse, and rescue")
        a.legend(fontsize=7)
        a = ax[1, 2]
        ws = B["width_scan"]
        a.plot(ws["w"], ws["zeta"], "k", lw=2, label="damping ζ (fast G–I loop)")
        a2 = a.twinx()
        a2.plot(ws["w"], np.minimum(ws["tau"], 2000), "C3", lw=2)
        a2.set_ylabel("critical insulin-action delay (min)", color="C3")
        a.axvline(f["w"], color="C2", ls=":")
        a.set_xlabel("threshold heterogeneity w")
        a.set_ylabel("ζ")
        a.set_title("B  heterogeneity buys damping:\nsynchronised cells → oscillation risk")
        a = ax[1, 3]
        lev = B["levers"]
        names = list(lev)
        col = {"healthy": "C2", "compensated hyperglycaemia": "C1", "collapse": "C3"}
        a.barh(range(len(names)), [lev[k]["G_fast"] for k in names],
               color=[col[lev[k]["fate"]["kind"]] for k in names])
        for i, k in enumerate(names):
            a.text(5, i, f"{lev[k]['fate']['kind']} (G {lev[k]['fate']['G']:.0f})",
                   va="center", fontsize=7, color="w")
        for i, k in enumerate(names):
            a.plot(lev[k]["G_s"], i, "k|", ms=18, mew=2)
        a.set_yticks(range(len(names)), names, fontsize=8)
        a.set_xlabel("acute glucose (mg/dl); | = saddle")
        a.set_title("B  levers on a diverging population\n(bar = acute G; text = fate with β ≤ 3×)")
        fig.tight_layout()
        p = os.path.join(FIG, "d9a_theory_beta_population.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        paths.append(p)
    # ---------------------------------------------------- d9b: CPTAC
    if "C" in R:
        C = R["C"]
        names = list(C["desc"])
        fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
        x = np.arange(len(names))
        for a, j, lab in ((ax[0, 0], 0, "participation ratio (effective dimension)"),
                          (ax[0, 1], 3, "mean |r| (coupling strength)")):
            for off_, key, c in ((-0.2, "N", "C0"), (0.2, "T", "C3")):
                v = np.array([C["desc"][n][key][:, j] for n in names])
                a.bar(x + off_, v.mean(1), 0.4, yerr=v.std(1), color=c,
                      label="normal" if key == "N" else "tumour")
            a.set_xticks(x, names, rotation=30)
            a.set_title(f"C  {lab}\n(equal n per ensemble)")
            a.legend(fontsize=8)
        a = ax[0, 2]
        U = C["universality"]
        im = a.imshow(U["overlap"], cmap="viridis", vmin=0, vmax=0.6)
        lbl = [f"{n}-{l}" for n, l in U["keys"]]
        a.set_xticks(range(len(lbl)), lbl, rotation=90, fontsize=7)
        a.set_yticks(range(len(lbl)), lbl, fontsize=7)
        plt.colorbar(im, ax=a)
        a.set_title(f"C  order-parameter universality\n(top-10 mode overlap; random {U['random']:.3f})")
        a = ax[1, 0]
        for n in names:
            rw = C["rewire"][n]
            ks = [0, 1, 3, 5]
            a.semilogy(ks, [max(rw[k]["n_obs"], 1) for k in ks], "o-", label=n)
        nul = np.mean([np.mean(C["rewire"][n][k]["null"]) for n in names for k in (0, 1, 3, 5)])
        a.axhline(nul, color="k", ls="--", label="permutation null (mean)")
        a.set_xlabel("collective modes subtracted per ensemble")
        a.set_ylabel("rewired pairs |Δz| > 5")
        a.set_title("C  all-pairs rewiring is mostly\na change of a few macro fields")
        a.legend(fontsize=7, ncol=2)
        a = ax[1, 1]
        pos = [np.median(np.abs(C["coupling"][n][("T", "res5")]["posFB"])) for n in names]
        neg = [np.median(np.abs(C["coupling"][n][("T", "res5")]["negFB"])) for n in names]
        one = [np.median(np.abs(np.r_[C["coupling"][n][("T", "res5")]["stim"],
                                      C["coupling"][n][("T", "res5")]["inh"]])) for n in names]
        rnd = [np.median(np.abs(C["coupling"][n][("T", "res5")]["random"])) for n in names]
        a.bar(x - 0.3, pos, 0.2, color="C3", label="positive 2-loops")
        a.bar(x - 0.1, one, 0.2, color="C7", label="one-way edges")
        a.bar(x + 0.1, neg, 0.2, color="C0", label="negative 2-loops")
        a.bar(x + 0.3, rnd, 0.2, color="0.8", label="random pairs")
        a.set_xticks(x, names, rotation=30)
        a.set_ylabel("median |r| (tumour, 5 modes removed)")
        a.set_title("C  OmniPath loops in co-abundance")
        a.legend(fontsize=7)
        a = ax[1, 2]
        e = C["energy"]
        nm = [n for n in e if "grade" in e[n]]
        for i, n in enumerate(nm):
            for j, (k, c) in enumerate((("energy", "C1"), ("euclid", "C7"), ("soft", "C2"))):
                a.bar(i + (j - 1) * 0.27, e[n]["grade"][k][0], 0.27, color=c,
                      label=k if i == 0 else None)
        a.axhline(0, color="k", lw=0.8)
        a.set_xticks(range(len(nm)), nm, rotation=30)
        a.set_ylabel("Spearman ρ with grade")
        a.set_title("C  Boltzmann energy vs plain distance\nvs soft-mode fraction")
        a.legend(fontsize=8)
        fig.tight_layout()
        p = os.path.join(FIG, "d9b_cptac_indications.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        paths.append(p)
    # ---------------------------------------------------- d9c: Perturb-seq
    if "D" in R:
        D = R["D"]
        fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
        for a, nm in ((ax[0, 0], "HMGCS1"), (ax[0, 1], "PHB")):
            if nm not in D["examples"]:
                continue
            c, d = D["examples"][nm]["scatter"]
            m = np.isfinite(d)
            a.scatter(c[m], d[m], s=3, alpha=0.3, color="0.5")
            top = np.argsort(-np.where(m, c, -np.inf))[:50]
            a.scatter(c[top], d[top], s=12, color="C3", label="top-50 co-fluctuating")
            a.axhline(0, color="k", lw=0.8)
            a.axvline(0, color="k", lw=0.8)
            a.set_xlabel(f"correlation with {nm} in control cells (fluctuation)")
            a.set_ylabel(f"change after {nm} knockdown (response, SD)")
            kind = "anti-FDT: partners RISE (homeostat)" if D["examples"][nm]["partner_mean"] > 0 \
                else "FDT-consistent: partners fall"
            a.set_title(f"D  {nm}: {kind}")
            a.legend(fontsize=8)
        a = ax[0, 2]
        fam = D["families"]
        ks = list(fam)
        y = np.arange(len(ks))
        a.barh(y - 0.2, [fam[k]["anti"] / fam[k]["tested"] for k in ks], 0.4, color="C3",
               label="anti-FDT (homeostat)")
        a.barh(y + 0.2, [fam[k]["fdt"] / fam[k]["tested"] for k in ks], 0.4, color="C0",
               label="FDT-consistent")
        for i, k in enumerate(ks):
            a.text(0.02 + max(fam[k]["anti"], fam[k]["fdt"]) / fam[k]["tested"], i,
                   f"n={fam[k]['tested']}", va="center", fontsize=7)
        a.set_yticks(y, ks, fontsize=8)
        a.set_xlabel("fraction of family knockdowns")
        a.set_title("D  pre-specified gene families")
        a.legend(fontsize=8)
        a = ax[1, 0]
        m = D["mra"]
        a.bar(range(4), [m["base"], m["prec_R"], m["prec_Rsym"], m["prec_mra"]],
              color=["0.7", "C7", "C7", "C3"])
        a.set_xticks(range(4), ["all pairs", "top |R|\n(total)", "top |R|+|Rᵀ|", "sparse MRA\n(direct)"])
        a.set_ylabel("fraction known (OmniPath / complexes)")
        a.set_title(f"D  direct vs total effects ({m['n_nodes']} nodes,\n"
                    f"{m['edges_both']} split-half edges, chance {m['edges_chance']:.0f})")
        a = ax[1, 1]
        rm = D["resp_modes"]
        a.semilogy(np.arange(1, 41), rm["lam"], "ko-", ms=3)
        a.axvline(rm["k_null"] + 0.5, color="C3", ls="--", label=f"{rm['k_null']} modes > null")
        a.set_xlabel("mode")
        a.set_ylabel("variance fraction")
        a.set_title(f"D  macro response modes (PR {rm['PR']:.0f} of 8,563 genes)")
        a.legend(fontsize=8)
        a = ax[1, 2]
        u = D["universality"]
        a.hist(u["per_pert_r"], bins=40, alpha=0.7, label="same knockdown")
        a.hist(u["per_pert_null"], bins=40, alpha=0.7, label="mismatched")
        a.set_xlabel("response correlation K562 vs RPE1")
        ov = u["overlap"][5]
        a.set_title(f"D  universality across cell types\n(top-5 mode overlap {ov[0]:.2f} vs random {ov[1]:.4f})")
        a.legend(fontsize=8)
        fig.tight_layout()
        p = os.path.join(FIG, "d9c_perturbseq_jacobian.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        paths.append(p)
    return paths


def main(parts="A,B,C,D,E", fresh=False):
    R = {} if fresh else _load_res()
    for p in parts.split(","):
        p = p.strip()
        t0 = time.time()
        R[p] = {"A": part_A, "B": part_B, "C": part_C, "D": part_D, "E": part_E}[p]()
        print(f"[part {p} done in {time.time() - t0:.0f}s]", flush=True)
        _save_res(R)
    report(R)
    for p in figures(R):
        print("figure:", p)
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="A,B,C,D,E")
    ap.add_argument("--report", action="store_true", help="only print + plot cached results")
    a = ap.parse_args()
    if a.report:
        R = _load_res()
        report(R)
        print(figures(R))
    else:
        main(a.part)
