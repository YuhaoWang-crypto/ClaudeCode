"""
Reusable helpers for the de-novo multi-modal drug-discovery campaign skill.
Deterministic / computational steps only — MCP data pulls and GPU dispatch are
documented in SKILL.md and reference/*.md (they run in the repl / compute_provider
kernels, not here).

Rigor discipline (from the network-biomarker layer): every function returns
computed numbers, never asserted constants. Label rigorous vs hypothesis in prose.
"""
import numpy as np


# ── STAGE 1: indication composite scoring ───────────────────────────────────
def composite_score(row, weights=None):
    """Weighted composite of an indication's 0-5 dimension scores.
    row: dict with keys funding,unmet,multimodal,biomarker,genetic_support,
         clinical_activity,crowding (each 0-5; crowding is inverse-desirable).
    Returns float composite (higher = more attractive)."""
    if weights is None:
        weights = {"funding":0.22,"unmet":0.18,"multimodal":0.20,"biomarker":0.15,
                   "genetic_support":0.12,"clinical_activity":0.08,"crowding":0.05}
    return round(sum(row[k]*w for k,w in weights.items()), 3)


def normalize_log(values):
    """Min-max normalize log10(1+x) to 0-5 (for target/drug counts)."""
    v = np.log10(1.0 + np.asarray(values, float))
    lo, hi = v.min(), v.max()
    if hi <= lo: return np.full_like(v, 2.5)
    return 5.0 * (v - lo) / (hi - lo)


# ── network-biomarker layer: fibration → CRNT δ → critical slowing ──────────
def input_fibers(edges):
    """Minimal balanced colouring by signed in-edge signature (Morone-Leifer-Makse).
    edges: list of (source, target, sign) with sign in {'+','-'}.
    Returns dict {fiber_signature: [nodes]} — nodes that synchronize."""
    nodes = set()
    for s,t,_ in edges: nodes.add(s); nodes.add(t)
    preds = {n: [] for n in nodes}
    for s,t,sign in edges: preds[t].append((s,sign))
    sig = {n: tuple(sorted(preds[n])) for n in nodes}
    fibers = {}
    for n,sg in sig.items(): fibers.setdefault(sg, []).append(n)
    return {str(k): sorted(v) for k,v in fibers.items()}


def crnt_deficiency(n_complexes, linkage_classes, stoich_rank):
    """Feinberg deficiency δ = n − ℓ − s. δ=0+weakly-reversible ⇒ no bistability;
    δ≥1 ⇒ bistability structurally permitted. Compute on the EFFECTIVE network
    (strip chemostatted species first — see reference/methodology gotcha)."""
    return int(n_complexes) - int(linkage_classes) - int(stoich_rank)


def schlogl_states(k1, k2, k3, p, xmax=12.0, ngrid=4000):
    """Steady states + stability of the canonical 1-species bistable CRN
    dx/dt = k1 x^2 - k2 x^3 + k3 - p x. Returns (roots, is_stable) sorted."""
    from scipy.optimize import brentq
    f = lambda x: k1*x**2 - k2*x**3 + k3 - p*x
    xs = np.linspace(0, xmax, ngrid); fx = f(xs); roots = []
    for i in range(len(xs)-1):
        if fx[i]*fx[i+1] < 0: roots.append(brentq(f, xs[i], xs[i+1]))
    roots = sorted(roots)
    stab = [ (f(r+1e-4)-f(r-1e-4))/2e-4 < 0 for r in roots ]
    return roots, stab


def critical_slowing(k1, k2, k3, p, branch="lo", dt=1.0):
    """Early-warning stats on a chosen branch of the bistable switch.
    Returns dict(nrf2, lam, recovery_tau, variance, ar1) — as drive p → fold,
    lam→0, tau→∞, variance & ar1 rise. branch='lo' tracks the disappearing OFF
    (sensitive) branch. Returns None if the branch is absent at this p."""
    roots, stab = schlogl_states(k1, k2, k3, p)
    if len(roots) < 2: return None
    x = roots[0] if branch == "lo" else roots[-1]
    f = lambda z: k1*z**2 - k2*z**3 + k3 - p*z
    lam = (f(x+1e-4)-f(x-1e-4))/2e-4
    if lam >= 0: return None
    return dict(nrf2=round(x,3), lam=round(lam,4), recovery_tau=round(-1/lam,2),
                variance=round(1/(2*abs(lam)),3), ar1=round(float(np.exp(lam*dt)),4))


def langevin_earlywarning(k1, k2, k3, p, sigma=0.06, T=4000, dt=0.05, seed=0, burn=1000):
    """Simulate the switch as a Langevin SDE and MEASURE variance + lag-1 AR
    (the rigorous, model-free statistic). Returns dict(var, ar1, series)."""
    rng = np.random.default_rng(seed)
    roots, _ = schlogl_states(k1, k2, k3, p)
    x = roots[0] if roots else 0.3
    f = lambda z: k1*z**2 - k2*z**3 + k3 - p*z
    xs = np.empty(T)
    for i in range(T):
        x = max(x + f(x)*dt + sigma*np.sqrt(dt)*rng.standard_normal(), 1e-6); xs[i] = x
    r = xs[burn:] - xs[burn:].mean()
    ar1 = float(np.corrcoef(r[:-1], r[1:])[0,1])
    return dict(var=float(r.var()), ar1=ar1, series=xs[burn::5].tolist())


# ── STAGE 8: consolidated lead ranking ──────────────────────────────────────
def rank_leads(leads, weights=None):
    """Rank executed leads by weighted composite of 0-5 dims.
    leads: list of dicts with feasibility,potency_evidence,safety,developability.
    Returns the list sorted desc with 'composite' and 'rank' added."""
    if weights is None:
        weights = {"feasibility":0.30,"potency_evidence":0.30,"safety":0.22,"developability":0.18}
    for L in leads:
        L["composite"] = round(sum(L[k]*w for k,w in weights.items()), 3)
    out = sorted(leads, key=lambda L: -L["composite"])
    for i,L in enumerate(out, 1): L["rank"] = i
    return out


# ── STAGE 9: siRNA-LNP knockdown (Mihaila 5-ODE, delegated to lnp skill) ────
def lnp_knockdown_from_escape(k2_escape):
    """Wrapper reminder: use the lnp-delivery-kinetics skill's mihaila_simulate
    for the validated 5-ODE model. This returns the escape→knockdown ordering
    note only. See SKILL.md ## Stage 9."""
    return {"note": "load lnp-delivery-kinetics; kd = 1 - M[-1]/M[0] from mihaila_simulate(k2=k2_escape)",
            "k2_escape": k2_escape}
