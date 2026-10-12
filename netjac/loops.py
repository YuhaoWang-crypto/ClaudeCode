"""
loops — from pairwise interactions to whole-network stability.

For a system dx/dt = f(x) (any dimension, given as a Python function or an
SBML model through core.KineticModel) this module answers:

  1. Where are the equilibria, and what is each one?  (node / focus /
     saddle / centre / non-hyperbolic) -> System.equilibria, analyse_point
  2. What does each k-variable SUB-system do if the other variables are
     clamped?  -> principal_subsystems (every 2x2, 3x3, ... principal block)
  3. How do the sub-systems COMBINE into the whole?  Exactly through the
     characteristic polynomial
         det(lambda I - J) = lambda^n + a_1 lambda^(n-1) + ... + a_n
     whose coefficients are sums over sets of DISJOINT feedback loops
     (Levins 1974 loop analysis):
         a_k = sum_{disjoint loop sets L covering k nodes} (-1)^|L| prod gain
     (equivalently a_k = (-1)^k * sum of all k x k principal minors).
     Feedback at level k: F_k = -a_k.  Local stability <=> Routh-Hurwitz on
     (a_1..a_n); necessary: every F_k < 0.  -> feedback_levels,
     loop_decomposition, hurwitz
  4. Sufficient shortcuts: diagonal dominance (Gershgorin), sign-pattern
     stability (only the SIGNS of interactions known) -> gershgorin,
     sign_pattern_stability
  5. How does the picture change with a parameter? -> bifurcation_scan
  6. Global behaviour from many initial states (linearisation is LOCAL) ->
     global_portrait
  7. Is bigger better? -> may_random (May 1972 complexity-stability)

Rigour: (1)-(4) are exact linear algebra at an equilibrium; (5)-(6) are
numerical; what a model's equations mean physiologically is the modeller's
assumption.
"""
from itertools import combinations

import numpy as np


# ---------------------------------------------------------------- systems
class System:
    def __init__(self, name, f, names, params=None, jac=None, nonneg=True):
        self.name, self._f, self.names = name, f, list(names)
        self.p = dict(params or {})
        self._jac = jac
        self.n = len(names)
        self.nonneg = nonneg

    def rhs(self, x, p=None):
        return np.asarray(self._f(np.asarray(x, float), {**self.p, **(p or {})}),
                          float)

    def jacobian(self, x, p=None, h=1e-7):
        if self._jac is not None:
            return np.asarray(self._jac(np.asarray(x, float),
                                        {**self.p, **(p or {})}), float)
        x = np.asarray(x, float)
        J = np.zeros((self.n, self.n))
        for k in range(self.n):
            d = h * max(1.0, abs(x[k]))
            xp, xm = x.copy(), x.copy()
            xp[k] += d
            xm[k] -= d
            J[:, k] = (self.rhs(xp, p) - self.rhs(xm, p)) / (2 * d)
        return J

    def equilibria(self, box, n_starts=300, p=None, seed=0, tol=1e-9):
        """Multi-start root finding inside box = [(lo, hi), ...]."""
        from scipy.optimize import root
        rng = np.random.default_rng(seed)
        found = []
        lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
        starts = [lo + (hi - lo) * rng.random(self.n) for _ in range(n_starts)]
        starts += [lo.copy()]                       # include the corner
        for x0 in starts:
            scale = 1.0 + np.abs(self.rhs(x0, p)).max()   # residual scale
            s = root(lambda x: self.rhs(x, p), x0, method="hybr", tol=1e-12)
            x = s.x
            if not np.isfinite(x).all():
                continue
            if np.abs(self.rhs(x, p)).max() > 1e-8 * scale:
                continue
            if self.nonneg and (x < -1e-9).any():
                continue
            x = np.where(np.abs(x) < 1e-10, 0.0, x)
            if not any(np.allclose(x, y, rtol=1e-5, atol=1e-7) for y in found):
                found.append(x)
        return sorted(found, key=lambda v: tuple(v))

    def simulate(self, x0, t_end, n=2000, p=None):
        from scipy.integrate import solve_ivp
        t = np.linspace(0, t_end, n)
        s = solve_ivp(lambda _t, x: self.rhs(x, p), (0, t_end), x0, t_eval=t,
                      method="LSODA", rtol=1e-8, atol=1e-10)
        return s.t, s.y.T


# ------------------------------------------------------- point analysis
def classify_eigs(vals, tol=1e-9):
    re = vals.real
    scale = max(np.median(np.abs(vals)), 1e-12)
    n_pos = int((re > tol * scale).sum())
    n_zero = int((np.abs(re) <= tol * scale).sum())
    lead = vals[np.argmax(re)]
    cplx = abs(lead.imag) > tol * scale
    if n_zero and n_pos == 0:
        kind = ("centre (non-hyperbolic: linearisation cannot decide; "
                "nonlinear terms do)" if cplx else
                "non-hyperbolic (zero eigenvalue: bifurcation point)")
    elif n_pos == 0:
        kind = "stable focus (damped oscillation)" if cplx else \
            "stable node (monotone return)"
    elif n_pos == len(vals):
        kind = "unstable focus (spirals out)" if cplx else "unstable node"
    else:
        up = vals[re > tol * scale]
        if n_pos == 2 and abs(up[0].imag) > tol * scale:
            kind = ("unstable focus pair (oscillatory instability -> limit "
                    "cycle if bounded; confirm by simulation)")
        else:
            kind = f"saddle ({n_pos} unstable direction{'s' * (n_pos > 1)})"
    return {"kind": kind, "lead": lead, "n_unstable": n_pos,
            "stable": n_pos == 0 and n_zero == 0,
            "tau": 1 / abs(lead.real) if lead.real != 0 else np.inf,
            "period": 2 * np.pi / abs(lead.imag) if cplx else None}


def analyse_point(J):
    vals = np.linalg.eigvals(J)
    out = classify_eigs(vals)
    out.update(eigenvalues=vals, trace=np.trace(J), det=np.linalg.det(J))
    return out


# ------------------------------------------------ sub-systems and loops
def principal_subsystems(J, k):
    """Every k x k principal block = the k-variable sub-network with the
    others clamped at their equilibrium values."""
    n = len(J)
    rows = []
    for idx in combinations(range(n), k):
        sub = J[np.ix_(idx, idx)]
        a = analyse_point(sub)
        rows.append({"idx": idx, "eigs": a["eigenvalues"], "kind": a["kind"],
                     "stable": a["stable"]})
    return rows


def char_coeffs(J):
    """a_1..a_n of det(lambda I - J) = lambda^n + a_1 lambda^(n-1) + ..."""
    return np.real(np.poly(J))[1:]


def feedback_levels(J):
    """Levins feedback F_k = -a_k (k = 1..n)."""
    return -char_coeffs(J)


def elementary_loops(J, tol=1e-12):
    """All elementary feedback loops (self-loops included) of the
    interaction digraph j -> i (J[i, j] != 0), with loop gain."""
    import networkx as nx
    n = len(J)
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(n):
            if i != j and abs(J[i, j]) > tol:
                G.add_edge(j, i)
    loops = [((i,), J[i, i]) for i in range(n) if abs(J[i, i]) > tol]
    for cyc in nx.simple_cycles(G):
        g = 1.0
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            g *= J[b, a]
        loops.append((tuple(cyc), g))
    return loops


def loop_decomposition(J, names=None, tol=1e-12):
    """a_k as a sum over sets of disjoint loops (Levins). Returns
    {k: [(description, contribution), ...]} and verifies against np.poly."""
    loops = elementary_loops(J, tol)
    n = len(J)
    terms = {k: [] for k in range(1, n + 1)}

    def rec(start, used, chosen):
        size = sum(len(l[0]) for l in chosen)
        if chosen:
            contrib = (-1) ** len(chosen) * np.prod([l[1] for l in chosen])
            terms[size].append((chosen, contrib))
        for i in range(start, len(loops)):
            nodes = set(loops[i][0])
            if nodes & used:
                continue
            rec(i + 1, used | nodes, chosen + [loops[i]])
    rec(0, set(), [])
    nm = names or [f"x{i}" for i in range(n)]
    out = {}
    for k, lst in terms.items():
        out[k] = [(" * ".join("(" + "->".join(nm[v] for v in l[0]) +
                              (f"->{nm[l[0][0]]}" if len(l[0]) > 1 else "")
                              + f": {l[1]:+.3g})" for l in ch), c)
                  for ch, c in lst]
    a = char_coeffs(J)
    check = max(abs(sum(c for _, c in out[k]) - a[k - 1]) for k in out)
    return out, check


def hurwitz(a):
    """Hurwitz determinants for lambda^n + a1 lambda^(n-1) + ... + an.
    All > 0  <=>  every root has negative real part."""
    n = len(a)
    c = np.r_[1.0, a]
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            k = 2 * j - i + 1
            if 0 <= k <= n:
                H[i, j] = c[k]
    return [np.linalg.det(H[:m, :m]) for m in range(1, n + 1)]


def gershgorin(J):
    """Sufficient stability: every row's self-damping beats the summed
    magnitude of incoming interactions (J_ii + sum_j |J_ij| < 0)."""
    m = np.diag(J) + np.abs(J).sum(1) - np.abs(np.diag(J))
    return bool((m < 0).all()), m


def sign_pattern_stability(J, n=4000, seed=0, spread=1.0):
    """Fraction of random magnitudes (log-normal, sd=spread) with the SAME
    sign pattern as J that are stable. 1.0 ~ qualitatively (sign-) stable:
    the signs alone guarantee convergence."""
    rng = np.random.default_rng(seed)
    S = np.sign(J)
    ok = 0
    for _ in range(n):
        M = S * np.exp(rng.normal(0, spread, J.shape))
        ok += np.linalg.eigvals(M).real.max() < 0
    return ok / n


# ------------------------------------------------------- bifurcations
def bifurcation_scan(sys_, pname, values, x0):
    """Follow one equilibrium branch along a parameter; flag where the
    leading eigenvalue's real part changes sign (real -> saddle-node /
    transcritical; complex -> Hopf)."""
    from scipy.optimize import root
    rows, x = [], np.asarray(x0, float)
    last_step = None
    for v in values:
        s = root(lambda z: sys_.rhs(z, {pname: v}), x, method="hybr")
        scale = 1.0 + np.abs(sys_.rhs(x, {pname: v})).max()
        step = np.linalg.norm(s.x - x)
        # a hop to another branch = a step far larger than the recent ones
        meaningful = (last_step is not None and
                      last_step > 1e-6 * (1 + np.linalg.norm(x)))
        jump = (meaningful and step > 10 * last_step and
                step > 0.05 * (1 + np.linalg.norm(x)))
        # judge by the residual, not s.success: hybr reports "no progress"
        # when it starts ON the solution
        if (jump or not np.isfinite(s.x).all() or
                np.abs(sys_.rhs(s.x, {pname: v})).max() > 1e-8 * scale):
            rows.append((v, None, None))
            continue
        last_step = step
        x = s.x
        vals = np.linalg.eigvals(sys_.jacobian(x, {pname: v}))
        rows.append((v, x.copy(), vals))
    events = []
    valid = [i for i, r in enumerate(rows) if r[2] is not None]
    if valid and valid[-1] < len(rows) - 1:
        i = valid[-1]
        events.append((rows[i + 1][0], "equilibrium branch disappears after "
                       f"{rows[i][0]:.6g} (saddle-node: collides with a saddle)"))
    for (v0, x0_, e0), (v1, x1_, e1) in zip(rows, rows[1:]):
        if e0 is None or e1 is None:
            continue
        l0, l1 = e0[np.argmax(e0.real)], e1[np.argmax(e1.real)]
        if np.sign(l0.real) != np.sign(l1.real) and l0.real * l1.real < 0:
            kind = "Hopf" if abs(l1.imag) > 1e-9 and abs(l0.imag) > 1e-9 \
                else "real eigenvalue crosses 0 (saddle-node / transcritical)"
            events.append((0.5 * (v0 + v1), kind))
        elif np.abs(e1.real).min() < 1e-6 * max(np.abs(e1).max(), 1e-12):
            # an eigenvalue touches 0 and the branch continues stable:
            # typical transcritical EXCHANGE of stability (two equilibria
            # collide and swap) -- continuation may hop onto the other branch
            events.append((v1, "zero eigenvalue: transcritical / branch "
                               "exchange"))
    return rows, events


# ---------------------------------------------------------- global view
def global_portrait(sys_, box, n_ic=60, t_end=200.0, seed=1, p=None):
    """Simulate from random initial states; classify where each ends:
    an equilibrium (which), a sustained oscillation, or divergence."""
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    eq = sys_.equilibria(box, p=p)
    res = []
    for _ in range(n_ic):
        x0 = lo + (hi - lo) * rng.random(sys_.n)
        t, X = sys_.simulate(x0, t_end, 4000, p)
        if len(X) == 0 or not np.isfinite(X).all() or np.abs(X).max() > 1e6:
            res.append("diverges")
            continue
        tail = X[len(X) // 2:]
        amp = (tail.max(0) - tail.min(0)).max()
        if amp > 1e-3 * (np.abs(tail).max() + 1e-9):
            res.append("oscillation")
            continue
        d = [np.linalg.norm(X[-1] - e) for e in eq]
        res.append(f"eq{int(np.argmin(d))}" if d and min(d) < 1e-3 *
                   (1 + np.linalg.norm(X[-1])) else "slow / unresolved")
    from collections import Counter
    return eq, Counter(res)


# --------------------------------------------------------------- May
def may_random(ns, C=0.3, sigma=0.5, d=1.0, trials=300, seed=0,
               structured=False):
    """P(stable) of random community matrices: diagonal -d, off-diagonal
    entries present with prob C, N(0, sigma). May (1972): stable w.h.p. iff
    sigma*sqrt(n*C) < d.  structured=True: pair every interaction with an
    opposite-signed reciprocal (predator-prey / negative-feedback pairs),
    the arrangement that is known to favour stability."""
    rng = np.random.default_rng(seed)
    out = []
    for n in ns:
        ok = 0
        for _ in range(trials):
            A = rng.normal(0, sigma, (n, n)) * (rng.random((n, n)) < C)
            if structured:
                # every interaction gets an opposite-signed reciprocal of
                # independent magnitude (predator-prey / negative-feedback
                # pairs; cf. Allesina & Tang 2012)
                U = np.triu(A, 1)
                R = np.abs(rng.normal(0, sigma, (n, n))) * (U != 0)
                A = U + (-np.sign(U) * np.triu(R, 1)).T
            np.fill_diagonal(A, -d)
            ok += np.linalg.eigvals(A).real.max() < 0
        out.append(ok / trials)
    return np.array(out)


# ------------------------------------------------- SBML (any) -> System
def sbml_system(source, name=None):
    """Wrap ANY SBML model (reaction-based, rate-rule ODEs, or mixed) as a
    System. State = rate-rule variables + floating species not fixed by
    rules; dx/dt is evaluated by libRoadRunner, so the curated equations are
    used verbatim. Returns (System, roadrunner instance)."""
    import roadrunner
    from netjac.core import fetch_biomodel, _load_sbml_string, quiet
    roadrunner.Logger.disableLogging()
    path = fetch_biomodel(source) if source.startswith("BIOMD") else source
    sbml, doc = _load_sbml_string(path)
    rr = roadrunner.RoadRunner(sbml)
    m = doc.getModel()
    roc = rr.getRatesOfChangeNamedArray()
    # every ODE state (rate rules included); names come back as "x'"
    names = [c[:-1] if c.endswith("'") else c for c in roc.colnames]
    sel = [f"{nm}'" for nm in names]
    defaults = {p: rr[p] for p in rr.model.getGlobalParameterIds()}

    changed = set()

    def f(x, p):
        # parameters not passed in this call go back to their SBML defaults
        # (roadrunner is stateful: a scan must not leak into later calls)
        for k in list(changed):
            if k not in p:
                rr[k] = defaults[k]
                changed.discard(k)
        for k, v in p.items():
            if k in defaults and rr[k] != v:
                rr[k] = float(v)
                changed.add(k)
        for nm, v in zip(names, x):
            rr[nm] = float(v)
        try:
            return np.asarray(rr.getRatesOfChange(), float)
        except RuntimeError:            # species driven by rate rules
            return np.array([rr[q] for q in sel], float)

    S = System(name or m.getName() or source, f, names, params={},
               nonneg=False)
    S.rr, S.defaults, S.quiet = rr, defaults, quiet

    def x0():
        return np.array([rr[nm] for nm in names])
    S.x_initial = x0()

    def attractor(t_end=2000.0, n=20000):
        """Simulate from the SBML initial state; report the late behaviour."""
        rr.reset()
        with quiet():
            r = rr.simulate(0, t_end, n, selections=["time"] + names)
        X = np.asarray(r)[:, 1:]
        t = np.asarray(r)[:, 0]
        tail = X[len(X) // 2:]
        tt = t[len(X) // 2:]
        rel_amp = (tail.max(0) - tail.min(0)) / (np.abs(tail).max(0) + 1e-12)
        j = int(np.argmax(rel_amp))
        y = tail[:, j] - np.polyval(np.polyfit(tt, tail[:, j], 1), tt)
        zc = np.where((y[:-1] < 0) & (y[1:] >= 0))[0]
        # oscillation = repeated up-crossings of the detrended tail, not just
        # a slow drift that has not settled within the horizon
        osc = rel_amp.max() > 1e-3 and len(zc) >= 3
        period = np.median(np.diff(tt[zc])) if osc else None
        state = ("oscillation" if osc else
                 "slow transient (not settled)" if rel_amp.max() > 1e-3
                 else "steady")
        return {"t": t, "X": X, "oscillating": bool(osc), "period": period,
                "x_end": X[-1], "rel_amp": rel_amp, "state": state}
    S.attractor = attractor
    return S, rr
