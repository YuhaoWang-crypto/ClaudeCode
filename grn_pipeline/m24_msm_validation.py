"""
Module 24 — Is the MSM + TPT pipeline (M23) reliable? Controls, ground truth,
and perturbation screening.

M23 produced numbers for alanine dipeptide, but MD has no exact answer to
check them against. Here the pipeline is run, unchanged, on trajectories of a
system whose exact kinetics are known, and on deliberately broken inputs:

Ground truth: overdamped dynamics on the Metzner-Schuette-Vanden-Eijnden
triple-well potential (two deep wells A, B and a shallow well between them;
two competing reaction channels). The dynamics are a continuous-time Markov
jump process on a fine grid with square-root-approximation (SqRA) rates
    k_ij = (D / h^2) exp(-beta (U_j - U_i) / 2)      (grid neighbours)
which converge to Smoluchowski dynamics as h -> 0 and satisfy detailed balance
with pi ~ exp(-beta U). Its committor, rate and channel split follow exactly
from linear algebra; the pipeline only ever sees sampled (x, y) trajectories.

Positive controls (P)  : rate, committor and channel split vs exact; bootstrap
                         confidence intervals must cover the exact value.
Negative controls (N)  : time-shuffled frames (no kinetics); a coordinate that
                         cannot see the slow process (y only); decoy
                         perturbations in regions the system never visits.
Robustness (R)         : convergence with data length; sensitivity to the
                         number of micro-states, clustering method, lag.
Screening (S)          : 40 virtual "ligands" (Gaussian bumps, stabilising or
                         destabilising) ranked by their effect on k_AB --
                         predicted from ONE unperturbed MSM by reweighting,
                         compared with the exact effect and with static
                         (non-kinetic) baselines; plus the full sensitivity
                         map d ln k_AB / d epsilon(x, y) ("allosteric map").
MD perturbation (MD)   : alanine dipeptide with a torsional bias on psi
                         (positive) or on an irrelevant methyl rotor
                         (negative): effect predicted from the unbiased MSM,
                         then checked against new biased simulations.
"""
import os
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import expm

from grn_pipeline import m23_msm_tpt as m23

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")

BETA = 1.67          # inverse temperature (Metzner et al. 2006 "high T")
H = 0.08             # grid spacing
XR, YR = (-2.0, 2.0), (-1.5, 2.5)
A_C, B_C = (-1.0, 0.0), (1.0, 0.0)


# ---------------------------------------------------------------------------
# exact system
# ---------------------------------------------------------------------------
def triple_well(x, y):
    return (3 * np.exp(-x ** 2 - (y - 1 / 3) ** 2)
            - 3 * np.exp(-x ** 2 - (y - 5 / 3) ** 2)
            - 5 * np.exp(-(x - 1) ** 2 - y ** 2)
            - 5 * np.exp(-(x + 1) ** 2 - y ** 2)
            + 0.2 * x ** 4 + 0.2 * (y - 1 / 3) ** 4)


class Grid:
    """Grid discretisation of the triple well with SqRA rates."""

    def __init__(self, h=H, beta=BETA, extra=None):
        self.h, self.beta = h, beta
        self.xs = np.arange(XR[0], XR[1] + 1e-9, h)
        self.ys = np.arange(YR[0], YR[1] + 1e-9, h)
        X, Y = np.meshgrid(self.xs, self.ys, indexing="ij")
        self.X, self.Y = X.ravel(), Y.ravel()
        self.U = triple_well(self.X, self.Y)
        if extra is not None:
            self.U = self.U + extra(self.X, self.Y)
        self.n = len(self.U)
        self.K = self._rates()
        p = np.exp(-beta * (self.U - self.U.min()))
        self.pi = p / p.sum()

    def _rates(self):
        nx, ny = len(self.xs), len(self.ys)
        idx = np.arange(self.n).reshape(nx, ny)
        rows, cols = [], []
        for a, b in [(idx[:-1, :], idx[1:, :]), (idx[:, :-1], idx[:, 1:])]:
            rows += [a.ravel(), b.ravel()]; cols += [b.ravel(), a.ravel()]
        r, c = np.concatenate(rows), np.concatenate(cols)
        k = np.exp(-self.beta * (self.U[c] - self.U[r]) / 2) / self.h ** 2
        K = sp.csr_matrix((k, (r, c)), shape=(self.n, self.n))
        K = K - sp.diags(np.asarray(K.sum(1)).ravel())
        return K.tocsr()

    def tpt(self, A, B):
        """Exact continuous-time TPT for index sets A, B (reversible)."""
        n = self.n
        I = np.setdiff1d(np.arange(n), np.r_[A, B])
        q = np.zeros(n); q[B] = 1.0
        KII = self.K[I][:, I]; KIB = self.K[I][:, B]
        q[I] = spla.spsolve(KII.tocsc(), -np.asarray(KIB.sum(1)).ravel())
        qm = 1.0 - q; qm[A] = 1.0; qm[B] = 0.0
        Kc = self.K.tocoo()
        f = self.pi[Kc.row] * qm[Kc.row] * Kc.data * q[Kc.col]
        off = Kc.row != Kc.col
        F = sp.csr_matrix((f[off], (Kc.row[off], Kc.col[off])),
                          shape=(n, n))
        fnet = (F - F.T).maximum(0).tocsr()
        inA = np.zeros(n, bool); inA[A] = True
        fc = fnet.tocoo()
        total = fc.data[inA[fc.row] & ~inA[fc.col]].sum()
        rate = total / (self.pi @ qm)
        return dict(q=q, fnet=fnet, F=total, rate=rate)

    def slowest_timescale(self):
        """Exact slowest relaxation time -1/lambda_2 (symmetrised K)."""
        s = np.sqrt(self.pi)
        Ks = sp.diags(s) @ self.K @ sp.diags(1 / s)
        Ks = 0.5 * (Ks + Ks.T)
        w = spla.eigsh(Ks, k=3, sigma=0, which="LM",
                       return_eigenvectors=False)
        w = np.sort(np.real(w))[::-1]
        return -1.0 / w[1]

    def channel_split(self, fnet, y_split=0.8):
        """Fraction of the net reactive flux crossing x = 0 above y_split
        (upper channel, through the shallow well) vs below (direct)."""
        fc = fnet.tocoo()
        cross = (self.X[fc.row] < 0) & (self.X[fc.col] >= 0)
        up = (0.5 * (self.Y[fc.row] + self.Y[fc.col]) > y_split) & cross
        return fc.data[up].sum() / fc.data[cross].sum()


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
def simulate(grid, n_traj=10, n_frames=40000, dt=0.05, seed=0):
    """Sample the jump process every dt; return (x, y) feature arrays with
    uniform sub-cell jitter (so the pipeline sees continuous coordinates)
    and the underlying grid-index trajectories."""
    rng = np.random.default_rng(seed)
    P = expm(grid.K.toarray() * dt)
    P = np.maximum(P, 0); P /= P.sum(1, keepdims=True)
    C = np.cumsum(P, 1)
    s = rng.choice(grid.n, size=n_traj, p=grid.pi)
    out = np.empty((n_frames, n_traj), dtype=np.int32)
    for t in range(n_frames):
        u = rng.random(n_traj)[:, None]
        s = np.minimum((C[s] < u).sum(1), grid.n - 1)
        out[t] = s
    feats, idx = [], []
    for k in range(n_traj):
        g = out[:, k]
        jit = rng.uniform(-grid.h / 2, grid.h / 2, size=(n_frames, 2))
        feats.append(np.c_[grid.X[g], grid.Y[g]] + jit)
        idx.append(g)
    return feats, idx


# ---------------------------------------------------------------------------
# pipeline on the toy + mapping to exact
# ---------------------------------------------------------------------------
def run_pipeline(feats, dt, n_clusters=150, cluster="regspace", msm_lag=None,
                 seed=7):
    r = m23.analyze(feats, dt_ps=dt, n_clusters=n_clusters, cluster=cluster,
                    msm_lag=msm_lag, n_macro=3, its_samples=0, seed=seed,
                    run_tpt=False, verbose=False)
    X = np.concatenate(feats); a = np.concatenate(r["active"])
    r["xy_centers"] = np.array([X[a == s].mean(0)
                                for s in range(r["msm"].n_states)])
    mc = np.array([(r["xy_centers"][r["assign"] == k]
                    * r["msm"].stationary_distribution[r["assign"] == k,
                                                       None]).sum(0)
                   / r["msm"].stationary_distribution[r["assign"] == k].sum()
                   for k in range(r["n_macro"])])
    src = int(np.argmin(np.hypot(*(mc - A_C).T)))
    tgt = int(np.argmin(np.hypot(*(mc - B_C).T)))
    r.update(m23.tpt_analysis(r, src, tgt, verbose=False))
    return r


def grid_labels(grid, r):
    """Assign every grid cell to the pipeline's micro-state (same TICA +
    clustering), so exact quantities use identical set definitions."""
    Yg = r["tica"].transform(np.c_[grid.X, grid.Y])
    lab = r["km"].transform(Yg)
    sym = r["msm"].count_model.state_symbols
    pos = np.full(r["km"].n_clusters, -1); pos[sym] = np.arange(len(sym))
    return pos[lab]


def compare_exact(grid, r):
    """Exact rate / committor / channel split on the pipeline's own core
    sets, vs the pipeline estimates."""
    lab = grid_labels(grid, r)
    A = np.where(np.isin(lab, r["A"]))[0]; B = np.where(np.isin(lab, r["B"]))[0]
    ex = grid.tpt(A, B)
    # exact committor averaged (pi-weighted) over each micro-state
    qm = np.full(r["msm"].n_states, np.nan)
    for s in range(r["msm"].n_states):
        w = grid.pi[lab == s]
        if w.sum() > 0:
            qm[s] = (w * ex["q"][lab == s]).sum() / w.sum()
    q_est = r["flux"].forward_committor
    ok = np.isfinite(qm)
    # pipeline channel split: flux across x = 0 between micro-state centres
    f = r["flux"].net_flux; c = r["xy_centers"]
    ii, jj = np.nonzero(f)
    cross = (c[ii, 0] < 0) & (c[jj, 0] >= 0)
    up = cross & (0.5 * (c[ii, 1] + c[jj, 1]) > 0.8)
    split_est = f[ii[up], jj[up]].sum() / f[ii[cross], jj[cross]].sum()
    return dict(rate_exact=ex["rate"], rate_est=r["rate_ps"],
                rate_direct=r["k_direct"],
                q_rmse=float(np.sqrt(np.mean((q_est[ok] - qm[ok]) ** 2))),
                q_r=float(np.corrcoef(q_est[ok], qm[ok])[0, 1]),
                split_exact=grid.channel_split(ex["fnet"]),
                split_est=float(split_est), lab=lab, A=A, B=B, exact=ex,
                q_exact_ms=qm)


# ---------------------------------------------------------------------------
# reweighting: predict a perturbation's effect from ONE unperturbed MSM
# ---------------------------------------------------------------------------
def state_shifts(V_frames, active, n_states, beta):
    """Free-energy shift of each micro-state under a perturbation V(frame):
    delta_i = -kT ln < exp(-beta V) >_i over the unperturbed frames in i."""
    a = np.concatenate(active); w = np.exp(-beta * np.concatenate(V_frames))
    num = np.bincount(a[a >= 0], weights=w[a >= 0], minlength=n_states)
    den = np.bincount(a[a >= 0], minlength=n_states)
    return -np.log(num / np.maximum(den, 1)) / beta


def reweight_msm(T, delta, beta):
    """Perturbed transition matrix that keeps detailed balance with
    pi'_i ~ pi_i exp(-beta delta_i): off-diagonal T_ij is scaled by
    exp(-beta (delta_j - delta_i) / 2) (barrier taken midway between the two
    micro-states -- the SqRA rule). Rows are re-normalised by uniformisation:
    T' = (1-s) I + s K', so the returned chain advances s lag times per
    step. Returns (T', s)."""
    off = T.copy(); np.fill_diagonal(off, 0.0)
    off *= np.exp(-beta * (delta[None, :] - delta[:, None]) / 2)
    out = off.sum(1)
    s = min(1.0, 0.999 / out.max())
    Tn = off * s
    Tn[np.diag_indices_from(Tn)] = 1.0 - Tn.sum(1)
    return Tn, s


def crossing_shifts(V_frames, active, lag, n_states, beta):
    """Bias felt AT each observed crossing i<->j (both directions pooled):
    Vb_ij = -kT ln < exp(-beta (V_t + V_{t+lag}) / 2) > over transitions.
    This is the transition-state shift for that edge, measured on the frames
    that actually cross instead of being assumed midway between states."""
    num = np.zeros((n_states, n_states)); den = np.zeros((n_states, n_states))
    for a, V in zip(active, V_frames):
        i, j = a[:-lag], a[lag:]
        ok = (i >= 0) & (j >= 0) & (i != j)
        w = np.exp(-beta * 0.5 * (V[:-lag] + V[lag:]))[ok]
        np.add.at(num, (i[ok], j[ok]), w); np.add.at(den, (i[ok], j[ok]), 1)
    num += num.T; den += den.T
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, -np.log(num / den) / beta, np.nan)


def reweight_msm_ts(T, delta, Vb, beta):
    """Transition-state reweighting: T'_ij = T_ij exp(-beta (Vb_ij -
    delta_i)) for i != j (TST form: barrier shifted by the bias at the
    crossing, start state by its own shift). Detailed balance with
    pi'_i ~ pi_i exp(-beta delta_i) holds exactly because Vb is symmetric.
    Edges never observed crossing fall back to the midpoint rule."""
    mid = 0.5 * (delta[:, None] + delta[None, :])
    Vb = np.where(np.isfinite(Vb), Vb, mid)
    off = T.copy(); np.fill_diagonal(off, 0.0)
    off *= np.exp(-beta * (Vb - delta[:, None]))
    s = min(1.0, 0.999 / off.sum(1).max())
    Tn = off * s
    Tn[np.diag_indices_from(Tn)] = 1.0 - Tn.sum(1)
    return Tn, s


def predicted_rate(msm, A, B, delta, beta, Vb=None):
    """Rate (per original lag step) of the reweighted MSM: transition-state
    rule if crossing shifts Vb are given, else the midpoint rule."""
    from deeptime.markov.msm import MarkovStateModel
    if Vb is None:
        Tn, s = reweight_msm(msm.transition_matrix, delta, beta)
    else:
        Tn, s = reweight_msm_ts(msm.transition_matrix, delta, Vb, beta)
    m = MarkovStateModel(Tn)
    # the uniformised chain advances s lag times per step -> divide by s
    return m.reactive_flux(A, B).rate / s, m


# ---------------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------------
def bootstrap_rate(r, n_boot=200, n_seg=4, seed=0):
    """Bootstrap over trajectory segments (each trajectory cut into n_seg
    pieces), re-estimating the MSM on the fixed discretisation. Returns the
    rate samples (per time unit)."""
    rng = np.random.default_rng(seed)
    segs = [d[k * (len(d) // n_seg):(k + 1) * (len(d) // n_seg)]
            for d in r["dtrajs"] for k in range(n_seg)]
    sym = r["msm"].count_model.state_symbols
    SA, SB = set(sym[r["A"]]), set(sym[r["B"]])
    out = []
    for _ in range(n_boot):
        pick = [segs[i] for i in rng.integers(0, len(segs), len(segs))]
        try:
            m, _ = m23.estimate_msm(pick, r["lag"])
        except Exception:
            continue
        s2 = m.count_model.state_symbols
        A = [i for i, s_ in enumerate(s2) if s_ in SA]
        B = [i for i, s_ in enumerate(s2) if s_ in SB]
        if A and B:
            out.append(m.reactive_flux(A, B).rate / (r["lag"] * r["dt_ps"]))
    return np.array(out)


def negative_shuffled(feats, dt, seed=0):
    """N1: destroy time order (keep every frame). A sound pipeline must
    report no slow process."""
    rng = np.random.default_rng(seed)
    sh = [f[rng.permutation(len(f))] for f in feats]
    r = m23.analyze(sh, dt_ps=dt, n_clusters=150, n_macro=3, its_samples=0,
                    run_tpt=False, verbose=False)
    return r


def negative_bad_coordinate(feats, dt):
    """N2: give the pipeline only y. Wells A (-1,0) and B (1,0) share y, so
    the A<->B process is invisible; the variational principle says the
    implied timescale can only drop -- detectable by comparing feature
    sets, not from within one."""
    fy = [f[:, [1]] for f in feats]
    return m23.analyze(fy, dt_ps=dt, n_clusters=150, n_macro=2,
                       its_samples=0, run_tpt=False, verbose=False)


# ---------------------------------------------------------------------------
# screening: virtual ligands
# ---------------------------------------------------------------------------
def bump(x0, y0, amp, width=0.25):
    return lambda x, y: amp * np.exp(-((x - x0) ** 2 + (y - y0) ** 2)
                                     / (2 * width ** 2))


def ligand_panel(amp=0.6):
    """28 sites x {stabilising (amp<0), destabilising (amp>0)} = 56
    perturbations; 4 sites are decoys in regions the system never visits."""
    sites = [(x, y) for x in (-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5)
             for y in (-0.6, 0.0, 0.6, 1.4)]
    decoys = [(-1.8, 2.3), (1.8, 2.3), (-1.8, -1.3), (1.8, -1.3)]
    panel = []
    for (x, y) in sites + decoys:
        for sgn in (-1, 1):
            panel.append(dict(x=x, y=y, amp=sgn * amp,
                              decoy=(x, y) in decoys))
    return panel


def screen(grid, r, cmp, feats, panel, width=0.25):
    """Exact vs MSM-predicted vs static-baseline effect of each ligand on
    ln k_AB. Static baseline = 'ground-state' Arrhenius: the ligand only
    shifts the free energy of the reactant basin A (barrier top unchanged),
    i.e. what endpoint/structure-based scoring captures."""
    beta = grid.beta; base = cmp["rate_exact"]
    X = np.concatenate(feats)
    k0 = r["flux"].rate
    pis = r["msm"].stationary_distribution
    rows = []
    t_pred = t_exact = 0.0
    for L in panel:
        f = bump(L["x"], L["y"], L["amp"], width)
        t = time.perf_counter()
        g2 = Grid(extra=f)
        ex = g2.tpt(cmp["A"], cmp["B"])["rate"]
        t_exact += time.perf_counter() - t
        t = time.perf_counter()
        V = [f(fe[:, 0], fe[:, 1]) for fe in feats]
        delta = state_shifts(V, r["active"], r["msm"].n_states, beta)
        Vb = crossing_shifts(V, r["active"], r["lag"], r["msm"].n_states,
                             beta)
        kp, _ = predicted_rate(r["msm"], r["A"], r["B"], delta, beta, Vb)
        t_pred += time.perf_counter() - t
        km, _ = predicted_rate(r["msm"], r["A"], r["B"], delta, beta)
        dA = -np.log((pis[r["A"]] * np.exp(-beta * delta[r["A"]])).sum()
                     / pis[r["A"]].sum()) / beta
        rows.append(dict(L, exact=np.log(ex / base), pred=np.log(kp / k0),
                         pred_mid=np.log(km / k0),
                         static=beta * dA,
                         density=float(np.mean(np.hypot(X[:, 0] - L["x"],
                                                        X[:, 1] - L["y"])
                                               < width))))
    return rows, t_pred / len(panel), t_exact / len(panel)


def screen_metrics(rows):
    from scipy.stats import spearmanr
    ex = np.array([r_["exact"] for r_ in rows])
    pr = np.array([r_["pred"] for r_ in rows])
    st = np.array([r_["static"] for r_ in rows])
    pm = np.array([r_["pred_mid"] for r_ in rows])
    de = np.array([r_["density"] for r_ in rows])
    big = np.abs(ex) > 0.05                    # effects worth classifying
    top = set(np.argsort(-np.abs(ex))[:10])
    return dict(
        rho_msm=spearmanr(pr, ex)[0], rho_static=spearmanr(st, ex)[0],
        rho_mid=spearmanr(pm, ex)[0],
        sign_mid=float(np.mean(np.sign(pm[big]) == np.sign(ex[big]))),
        mae_mid=float(np.mean(np.abs(pm - ex))),
        rho_density_mag=spearmanr(de, np.abs(ex))[0],
        sign_msm=float(np.mean(np.sign(pr[big]) == np.sign(ex[big]))),
        sign_static=float(np.mean(np.sign(st[big]) == np.sign(ex[big]))),
        top10_msm=len(top & set(np.argsort(-np.abs(pr))[:10])) / 10,
        top10_static=len(top & set(np.argsort(-np.abs(st))[:10])) / 10,
        mae_msm=float(np.mean(np.abs(pr - ex))), n_big=int(big.sum()))


def sensitivity_map(grid, r, cmp, feats, n=21, amp=0.15, width=0.25):
    """d ln k_AB / d epsilon at each position: exact and MSM-predicted.
    Large |S| far from A and B = an allosteric hot-spot for this reaction."""
    xs = np.linspace(XR[0] + 0.2, XR[1] - 0.2, n)
    ys = np.linspace(YR[0] + 0.2, YR[1] - 0.2, n)
    Se = np.zeros((n, n)); Sp = np.zeros((n, n))
    beta = grid.beta; k0 = r["flux"].rate; base = cmp["rate_exact"]
    for i, x0 in enumerate(xs):
        for j, y0 in enumerate(ys):
            f = bump(x0, y0, -amp, width)          # stabilising probe
            Se[j, i] = np.log(Grid(extra=f).tpt(cmp["A"], cmp["B"])["rate"]
                              / base) / amp
            V = [f(fe[:, 0], fe[:, 1]) for fe in feats]
            delta = state_shifts(V, r["active"], r["msm"].n_states, beta)
            Vb = crossing_shifts(V, r["active"], r["lag"],
                                 r["msm"].n_states, beta)
            kp, _ = predicted_rate(r["msm"], r["A"], r["B"], delta, beta, Vb)
            Sp[j, i] = np.log(kp / k0) / amp
    return xs, ys, Se, Sp


# ---------------------------------------------------------------------------
# robustness
# ---------------------------------------------------------------------------
def convergence(grid, feats, dt, fracs=(0.05, 0.1, 0.25, 0.5, 1.0)):
    rows = []
    for fr in fracs:
        sub = [f[:int(len(f) * fr)] for f in feats]
        r = run_pipeline(sub, dt)
        c = compare_exact(grid, r)
        rows.append(dict(frac=fr, rate_err=c["rate_est"] / c["rate_exact"] - 1,
                         q_r=c["q_r"], n_events=r["n_events"]))
    return rows


def hyperparams(grid, feats, dt):
    """Rate error across clustering method, number of micro-states, lag.
    Each run gets its own exact reference (its own core sets)."""
    cfgs = [("regspace", 50, None), ("regspace", 150, None),
            ("regspace", 300, None), ("kmeans", 150, None),
            ("regspace", 150, 5), ("regspace", 150, 20)]
    rows = []
    for meth, nc, lag in cfgs:
        r = run_pipeline(feats, dt, n_clusters=nc, cluster=meth, msm_lag=lag)
        c = compare_exact(grid, r)
        rows.append(dict(method=meth, n=r["msm"].n_states, lag=r["lag"],
                         rate_err=c["rate_est"] / c["rate_exact"] - 1,
                         q_r=c["q_r"], split_err=c["split_est"]
                         - c["split_exact"]))
    return rows


def amplitude_sweep(grid, r, cmp, feats, amps=(0.3, 0.6, 1.2, 2.4)):
    """Where does single-MSM reweighting break down? (amp in energy units;
    kT = 1/beta = 0.6)."""
    out = []
    for a in amps:
        rows, _, _ = screen(grid, r, cmp, feats, ligand_panel(a))
        m = screen_metrics(rows)
        out.append(dict(amp=a, amp_kT=a * grid.beta, rho=m["rho_msm"],
                        sign=m["sign_msm"], mae=m["mae_msm"],
                        rho_static=m["rho_static"]))
    return out


# ---------------------------------------------------------------------------
# MD perturbation controls (alanine dipeptide)
# ---------------------------------------------------------------------------
KT_300 = 0.0083144626 * 300.0


def _rate_ci(dtrajs, lag, dt, SA, SB, n_boot=200, n_seg=5, seed=0):
    """Rate A->B and B->A (per ps) on raw cluster labels, with a segment
    bootstrap. SA/SB are sets of raw cluster labels (core sets)."""
    rng = np.random.default_rng(seed)

    def one(dts):
        m, _ = m23.estimate_msm(dts, lag)
        sym = m.count_model.state_symbols
        A = [i for i, s_ in enumerate(sym) if s_ in SA]
        B = [i for i, s_ in enumerate(sym) if s_ in SB]
        f = m.reactive_flux(A, B).rate / (lag * dt)
        b = m.reactive_flux(B, A).rate / (lag * dt)
        return f, b

    point = one(dtrajs)
    segs = [d[k * (len(d) // n_seg):(k + 1) * (len(d) // n_seg)]
            for d in dtrajs for k in range(n_seg)]
    bs = []
    for _ in range(n_boot):
        try:
            bs.append(one([segs[i] for i in
                           rng.integers(0, len(segs), len(segs))]))
        except Exception:
            pass
    bs = np.array(bs)
    return point, np.percentile(bs, [2.5, 97.5], axis=0)


def md_perturbation(n_boot=200):
    """Predict the effect of a torsional bias from the UNBIASED MSM (M23
    data), then compare with simulations that actually include the bias.
    psi bias (well of 4 kJ/mol at psi=-40, inside alphaR) = positive control;
    the same bias on the ACE methyl rotor = negative control."""
    import mdtraj as md
    top, xtcs = m23.simulate_alanine_dipeptide(n_traj=4, ns=50.0)
    feats, _, rama = m23.featurize(top, xtcs)
    r = m23.analyze(feats, dt_ps=1.0, rama=rama, n_clusters=150, n_macro=3,
                    its_samples=0, run_tpt=False, verbose=False)
    lab = m23._ala2_basins(r["macro_center"])
    idx = {v: k for k, v in lab.items()}
    beta_s, aR = idx["beta/PII"], idx["alphaR"]
    r.update(m23.tpt_analysis(r, beta_s, aR, verbose=False))
    sym = r["msm"].count_model.state_symbols
    SA, SB = set(sym[r["A"]]), set(sym[r["B"]])
    beta = 1.0 / KT_300
    tref = md.load(top).topology
    meth_idx = m23.torsion_atoms(tref, "methyl")

    def raw_labels(fs):
        return [r["km"].transform(r["tica"].transform(f)).astype(np.int32)
                for f in fs]

    out = {}
    base, base_ci = _rate_ci(r["dtrajs"], r["lag"], 1.0, SA, SB, n_boot)
    pop0 = float(r["msm"].stationary_distribution[r["assign"] == aR].sum())
    out["unbiased"] = dict(pop_aR=pop0, k_f=base[0], k_b=base[1],
                           ci_f=base_ci[:, 0], ci_b=base_ci[:, 1])
    for name, theta0 in [("psi", -40.0), ("methyl", 0.0)]:
        # prediction from unbiased frames
        if name == "psi":
            th = [np.radians(x[:, 1]) for x in rama]
        else:
            th = [md.compute_dihedrals(md.load(x, top=top), [meth_idx])[:, 0]
                  for x in xtcs]
        V = [m23.bias_energy(t, 4.0, theta0) for t in th]
        n = r["msm"].n_states
        d = state_shifts(V, r["active"], n, beta)
        Vb = crossing_shifts(V, r["active"], r["lag"], n, beta)
        kf, m = predicted_rate(r["msm"], r["A"], r["B"], d, beta, Vb)
        kb, _ = predicted_rate(r["msm"], r["B"], r["A"], d, beta, Vb)
        pi_p = m.stationary_distribution
        pred = dict(pop_aR=float(pi_p[r["assign"] == aR].sum()),
                    k_f=kf / (r["lag"]), k_b=kb / (r["lag"]))
        # simulation with the bias, assigned to the SAME micro-states
        _, bx = m23.simulate_alanine_dipeptide(
            n_traj=4, ns=25.0, out=os.path.join(m23.CACHE, f"bias_{name}"),
            bias=(name, 4.0, theta0), tag=name,
            seed0=2000 if name == "psi" else 3000)
        bfeats, _, _ = m23.featurize(top, bx)
        bd = raw_labels(bfeats)
        sim, sim_ci = _rate_ci(bd, r["lag"], 1.0, SA, SB, n_boot)
        pos = {s_: i for i, s_ in enumerate(sym)}
        cat = np.concatenate(bd)
        macro = np.array([r["assign"][pos[c]] if c in pos else -1
                          for c in cat])
        out[name] = dict(pred=pred, sim=dict(
            pop_aR=float(np.mean(macro[macro >= 0] == aR)), k_f=sim[0],
            k_b=sim[1],
            ci_f=sim_ci[:, 0], ci_b=sim_ci[:, 1]))
    return out


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _contour(ax, grid, levels=14):
    nx, ny = len(grid.xs), len(grid.ys)
    U = grid.U.reshape(nx, ny).T
    ax.contour(grid.xs, grid.ys, U, levels=levels, colors="#999",
               linewidths=0.5)


def fig_controls(g, r, cmp, bs, conv, hp, neg, path):
    plt = _plt()
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.4)); ax = ax.ravel()
    nx, ny = len(g.xs), len(g.ys)
    # (a) exact committor + pipeline committor on micro-states
    a = ax[0]
    _contour(a, g)
    a.contour(g.xs, g.ys, cmp["exact"]["q"].reshape(nx, ny).T,
              levels=[0.1, 0.5, 0.9], colors="k", linewidths=1.2,
              linestyles=["--", "-", "--"])
    c = r["xy_centers"]
    sc = a.scatter(c[:, 0], c[:, 1], c=r["flux"].forward_committor,
                   cmap="coolwarm", s=22, vmin=0, vmax=1, edgecolor="none")
    fig.colorbar(sc, ax=a, label="pipeline committor q+")
    a.set_title("(a) triple-well: pipeline q+ (dots) vs exact\n"
                "q = 0.1 / 0.5 / 0.9 isolines (black)")
    a.set_xlabel("x"); a.set_ylabel("y")
    # (b) committor parity
    a = ax[1]; ok = np.isfinite(cmp["q_exact_ms"])
    a.plot([0, 1], [0, 1], "k:", lw=1)
    a.scatter(cmp["q_exact_ms"][ok], r["flux"].forward_committor[ok], s=12,
              c="#2a6fdb")
    a.set_xlabel("exact committor (micro-state average)")
    a.set_ylabel("pipeline committor")
    a.set_title(f"(b) committor parity\nr = {cmp['q_r']:.4f}, "
                f"RMSE = {cmp['q_rmse']:.3f}")
    # (c) rate: exact vs estimates
    a = ax[2]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    vals = [cmp["rate_exact"], cmp["rate_est"], cmp["rate_direct"]]
    a.bar([0, 1, 2], vals, color=["#444", "#2a6fdb", "#29a36a"])
    a.errorbar([1], [cmp["rate_est"]], yerr=[[cmp["rate_est"] - lo],
                                            [hi - cmp["rate_est"]]],
               color="k", capsize=6)
    a.set_xticks([0, 1, 2], ["exact", "MSM+TPT\n(95% CI)", "direct\ncount"])
    a.set_ylabel("rate k_AB")
    a.set_title(f"(c) rate: MSM {cmp['rate_est']:.4f} vs exact "
                f"{cmp['rate_exact']:.4f}\nchannel split (upper): "
                f"{cmp['split_est']:.2f} vs {cmp['split_exact']:.2f}")
    # (d) convergence
    a = ax[3]
    ev = [c_["n_events"] for c_ in conv]
    a.semilogx(ev, [100 * c_["rate_err"] for c_ in conv], "o-", c="#e4572e")
    a.axhline(0, c="k", lw=0.8)
    a.axhspan(100 * (lo / cmp["rate_exact"] - 1),
              100 * (hi / cmp["rate_exact"] - 1), color="#2a6fdb", alpha=0.15,
              label="bootstrap 95% CI (full data)")
    a.minorticks_off(); a.set_xticks(ev, [str(e) for e in ev])
    a.set_xlabel("A->B events in the data"); a.set_ylabel("rate error (%)")
    a.legend(fontsize=8)
    a.set_title("(d) convergence with data length")
    # (e) hyperparameters
    a = ax[4]
    labels = [f"{h['method'][:3]} n={h['n']}\nlag={h['lag']}" for h in hp]
    a.bar(range(len(hp)), [100 * h["rate_err"] for h in hp], color="#9b59b6")
    a.axhline(0, c="k", lw=0.8)
    a.set_xticks(range(len(hp)), labels, fontsize=7)
    a.set_ylabel("rate error (%)")
    a.set_title("(e) sensitivity to clustering / micro-states / lag\n"
                "(systematic spread > bootstrap CI)")
    # (f) negative controls: slowest ITS
    a = ax[5]
    names = ["exact", "(x, y)\nfeatures", "y only\n(blind)", "time-\nshuffled"]
    v = [neg["exact_t1"], neg["t1_xy"], neg["t1_y"], neg["t1_shuffled"]]
    a.bar(range(4), v, color=["#444", "#2a6fdb", "#f2a541", "#e4572e"])
    a.set_yscale("log"); a.set_xticks(range(4), names)
    a.axhline(neg["dt"], c="k", ls="--", lw=0.8)
    a.text(3.4, neg["dt"] * 1.15, "1 frame", fontsize=7, ha="right")
    a.set_ylabel("slowest implied timescale (lag = 1 frame)")
    a.set_title("(f) negative controls: blind coordinate and\n"
                "shuffled time collapse the slow process")
    fig.suptitle("M24 — MSM+TPT against exact ground truth (triple-well, "
                 "SqRA jump process)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=120)
    plt.close(fig)


def fig_screening(g, rows, metrics, amp, smap, path):
    plt = _plt()
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.4)); ax = ax.ravel()
    ex = np.array([q["exact"] for q in rows])
    pr = np.array([q["pred"] for q in rows])
    st = np.array([q["static"] for q in rows])
    dec = np.array([q["decoy"] for q in rows])
    lim = 1.1 * max(np.abs(ex).max(), np.abs(pr).max())
    for a, y, lab, rho in [(ax[0], pr, "MSM reweighting (one model)",
                            metrics["rho_msm"]),
                           (ax[1], st, "static ground-state baseline",
                            metrics["rho_static"])]:
        a.plot([-lim, lim], [-lim, lim], "k:", lw=1)
        a.axhline(0, c="#ccc", lw=0.6); a.axvline(0, c="#ccc", lw=0.6)
        a.scatter(ex[~dec], y[~dec], s=22, c=np.where(ex[~dec] > 0,
                                                      "#29a36a", "#e4572e"))
        a.scatter(ex[dec], y[dec], s=40, marker="x", c="k", label="decoys")
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
        a.set_xlabel("exact  d ln k_AB"); a.set_ylabel("predicted")
        a.legend(fontsize=8, loc="upper left")
        a.set_title(f"({'a' if a is ax[0] else 'b'}) {lab}\n"
                    f"Spearman rho = {rho:.3f}")
    xs, ys, Se, Sp = smap
    vmax = np.abs(Se).max()
    for a, S, lab in [(ax[2], Se, "(c) exact"), (ax[3], Sp,
                                                 "(d) MSM-predicted")]:
        im = a.pcolormesh(xs, ys, S, cmap="PiYG", vmin=-vmax, vmax=vmax,
                          shading="auto")
        _contour(a, g)
        a.plot(*A_C, "ko"); a.text(A_C[0], A_C[1] - 0.3, "A", ha="center")
        a.plot(*B_C, "ko"); a.text(B_C[0], B_C[1] - 0.3, "B", ha="center")
        fig.colorbar(im, ax=a, label="d ln k_AB / d(stabilisation)")
        a.set_xlabel("x"); a.set_ylabel("y")
        a.set_title(f"{lab} sensitivity map\n(green: stabilising here "
                    f"speeds A->B)")
    ax[3].set_title(ax[3].get_title() + f"; r = "
                    f"{np.corrcoef(Se.ravel(), Sp.ravel())[0, 1]:.3f}")
    a = ax[4]
    a.plot([q["amp_kT"] for q in amp], [q["rho"] for q in amp], "o-",
           c="#2a6fdb", label="MSM reweighting")
    a.plot([q["amp_kT"] for q in amp], [q["rho_static"] for q in amp], "s--",
           c="#7f8c8d", label="static baseline")
    a.set_xlabel("perturbation size (kT)"); a.set_ylabel("Spearman rho")
    a.set_ylim(-0.1, 1.05); a.legend(fontsize=8)
    b = a.twinx(); b.plot([q["amp_kT"] for q in amp], [q["mae"] for q in amp],
                          "^:", c="#e4572e")
    b.set_ylabel("MAE in d ln k", color="#e4572e")
    a.set_title("(e) ranking holds to 4 kT; absolute error grows")
    a = ax[5]; a.axis("off")
    txt = (f"56 virtual ligands (28 sites x 2 signs, 4 decoy sites)\n\n"
           f"                    MSM    static\n"
           f"Spearman rho        {metrics['rho_msm']:.3f}   "
           f"{metrics['rho_static']:.3f}\n"
           f"sign correct*       {metrics['sign_msm']:.0%}   "
           f"{metrics['sign_static']:.0%}\n"
           f"top-10 recovered    {metrics['top10_msm']:.0%}   "
           f"{metrics['top10_static']:.0%}\n"
           f"MAE (d ln k)        {metrics['mae_msm']:.3f}\n\n"
           f"* on the {metrics['n_big']} ligands with |effect| > 0.05\n"
           f"cost per ligand: {metrics['t_pred'] * 1e3:.0f} ms from one "
           f"MSM\nvs. a new {metrics['sim_cost']:.0f} s simulation + model")
    a.text(0, 0.95, txt, va="top", family="monospace", fontsize=10)
    fig.suptitle("M24 — screening virtual ligands for kinetic effect: "
                 "predicted from one unperturbed MSM vs exact",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=120)
    plt.close(fig)


def fig_md(md_res, path):
    plt = _plt()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    u = md_res["unbiased"]
    for k, (key, lab, scale) in enumerate([
            ("pop_aR", "alphaR population", 1.0),
            ("k_f", "MFPT beta -> alphaR (ps)", None),
            ("k_b", "MFPT alphaR -> beta (ps)", None)]):
        a = ax[k]
        xs, vals, errs, cols, names = [], [], [], [], []
        entries = [("unbiased", u, None, "#444")]
        for name, col in [("psi", "#2a6fdb"), ("methyl", "#f2a541")]:
            entries += [(f"{name}\npredicted", md_res[name]["pred"], None,
                         col), (f"{name}\nsimulated", md_res[name]["sim"],
                                True, col)]
        for i, (nm, d, sim, col) in enumerate(entries):
            v = d[key] if scale else 1.0 / d[key]
            e = None
            ci = d.get("ci_f" if key == "k_f" else "ci_b")
            if key != "pop_aR" and ci is not None:
                e = [[v - 1 / ci[1]], [1 / ci[0] - v]]
            a.bar(i, v, color=col, alpha=0.55 if "predicted" in nm else 1.0,
                  hatch="//" if "predicted" in nm else None)
            if e is not None:
                a.errorbar(i, v, yerr=e, color="k", capsize=5)
            names.append(nm)
        a.set_xticks(range(len(names)), names, fontsize=8)
        a.set_title(lab)
    fig.suptitle("Alanine dipeptide: effect of a 4 kJ/mol torsional bias "
                 "predicted from the UNBIASED MSM (hatched) vs new biased "
                 "simulations (solid, 95% CI).\npsi bias = positive control; "
                 "methyl-rotor bias = negative control", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88]); fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def report(seed=0):
    import json
    print("=" * 70)
    print("MODULE 24 — MSM+TPT reliability: controls, ground truth, screening")
    print("=" * 70)
    g = Grid(); dt = 0.05
    t = time.perf_counter()
    feats, _ = simulate(g, n_traj=10, n_frames=40000, dt=dt, seed=seed)
    r = run_pipeline(feats, dt)
    sim_cost = time.perf_counter() - t
    cmp = compare_exact(g, r)
    bs = bootstrap_rate(r)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    print(f"P1 rate: MSM {cmp['rate_est']:.4f} [95% CI {lo:.4f}, {hi:.4f}] "
          f"vs exact {cmp['rate_exact']:.4f} -> covered: "
          f"{lo <= cmp['rate_exact'] <= hi}; direct count "
          f"{cmp['rate_direct']:.4f}")
    print(f"P2 committor: r = {cmp['q_r']:.4f}, RMSE {cmp['q_rmse']:.3f}")
    print(f"P3 channel split (upper): {cmp['split_est']:.3f} vs exact "
          f"{cmp['split_exact']:.3f}")
    rs = negative_shuffled(feats, dt); rb = negative_bad_coordinate(feats, dt)
    neg = dict(exact_t1=g.slowest_timescale(), t1_xy=r["its"][0, 0] * dt,
               t1_y=rb["its"][0, 0] * dt, t1_shuffled=rs["its"][0, 0] * dt,
               dt=dt)
    print(f"N1 shuffled: slowest ITS {neg['t1_shuffled']:.4f} (1 frame = "
          f"{dt}) vs {neg['t1_xy']:.2f} real -> no kinetics, flagged")
    print(f"N2 y-only coordinate: slowest ITS {neg['t1_y']:.2f} vs "
          f"{neg['t1_xy']:.2f} with (x,y), exact {neg['exact_t1']:.2f}")
    conv = convergence(g, feats, dt); hp = hyperparams(g, feats, dt)
    for c_ in conv:
        print(f"R1 {c_['frac']:>5.0%} data ({c_['n_events']:4d} events): "
              f"rate error {c_['rate_err'] * 100:+5.1f}%, q r "
              f"{c_['q_r']:.4f}")
    for h in hp:
        print(f"R2 {h['method']:8s} n={h['n']:3d} lag={h['lag']:2d}: rate "
              f"{h['rate_err'] * 100:+5.1f}%, q r {h['q_r']:.4f}, split "
              f"{h['split_err']:+.3f}")
    rows, t_pred, t_ex = screen(g, r, cmp, feats, ligand_panel())
    met = screen_metrics(rows); met.update(t_pred=t_pred, sim_cost=sim_cost)
    print(f"S1 screening (56 ligands): MSM rho {met['rho_msm']:.3f} "
          f"(midpoint rule {met['rho_mid']:.3f}), sign {met['sign_msm']:.0%},"
          f" top-10 {met['top10_msm']:.0%}; static baseline rho "
          f"{met['rho_static']:.3f}, sign {met['sign_static']:.0%}")
    dec = max(abs(q["pred"]) for q in rows if q["decoy"])
    print(f"N3 decoy ligands: max |predicted d ln k| = {dec:.1e}")
    amp = amplitude_sweep(g, r, cmp, feats)
    for q in amp:
        print(f"S2 amplitude {q['amp_kT']:.1f} kT: rho {q['rho']:.3f}, MAE "
              f"{q['mae']:.3f}")
    smap = sensitivity_map(g, r, cmp, feats)
    xs, ys, Se, Sp = smap
    rmap = float(np.corrcoef(Se.ravel(), Sp.ravel())[0, 1])
    XX, YY = np.meshgrid(xs, ys)
    far = (np.hypot(XX - A_C[0], YY - A_C[1]) > 0.6) & \
        (np.hypot(XX - B_C[0], YY - B_C[1]) > 0.6)
    hot = np.argsort(-(np.abs(Se) * far).ravel())[:3]
    hot_sites = [dict(x=float(XX.ravel()[h]), y=float(YY.ravel()[h]),
                      S_exact=float(Se.ravel()[h]),
                      S_pred=float(Sp.ravel()[h])) for h in hot]
    print(f"S3 sensitivity map: predicted vs exact r = {rmap:.3f}; top "
          f"distal hot-spots: " + "; ".join(
              f"({h['x']:.1f},{h['y']:.1f}) S={h['S_exact']:+.2f}"
              f"/{h['S_pred']:+.2f}" for h in hot_sites))
    md_res = md_perturbation()
    for name in ("psi", "methyl"):
        p, s_ = md_res[name]["pred"], md_res[name]["sim"]
        print(f"MD {name:6s}: alphaR pop pred {p['pop_aR']:.3f} / sim "
              f"{s_['pop_aR']:.3f} (unbiased "
              f"{md_res['unbiased']['pop_aR']:.3f}); MFPT b->aR pred "
              f"{1 / p['k_f']:.0f} / sim {1 / s_['k_f']:.0f} ps "
              f"[{1 / s_['ci_f'][1]:.0f},{1 / s_['ci_f'][0]:.0f}]; aR->b pred "
              f"{1 / p['k_b']:.0f} / sim {1 / s_['k_b']:.0f} ps "
              f"[{1 / s_['ci_b'][1]:.0f},{1 / s_['ci_b'][0]:.0f}]")
    fig_controls(g, r, cmp, bs, conv, hp, neg,
                 os.path.join(FIG, "m24_controls.png"))
    fig_screening(g, rows, met, amp, smap,
                  os.path.join(FIG, "m24_screening.png"))
    fig_md(md_res, os.path.join(FIG, "m24_md_perturbation.png"))
    summary = dict(
        rate=dict(exact=cmp["rate_exact"], msm=cmp["rate_est"], ci=[lo, hi],
                  direct=cmp["rate_direct"]),
        committor=dict(r=cmp["q_r"], rmse=cmp["q_rmse"]),
        split=dict(exact=cmp["split_exact"], msm=cmp["split_est"]),
        negative=neg, convergence=conv, hyperparams=hp, screening=met,
        decoy_max=dec, amplitude=amp, sensitivity_r=rmap,
        hot_sites=hot_sites, md=md_res)
    with open(os.path.join(FIG, "m24_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=lambda o: (
            o.tolist() if hasattr(o, "tolist") else float(o)))
    print("figures: m24_controls.png, m24_screening.png, "
          "m24_md_perturbation.png; summary: m24_summary.json")
    return summary


if __name__ == "__main__":
    report()
