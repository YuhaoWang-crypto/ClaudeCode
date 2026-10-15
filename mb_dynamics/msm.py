"""R2: conformational states and a reversible continuous-time generator from unbiased CG trajectories.

Generator (column convention, L[i, j] = rate j -> i) with pi fixed to trajectory weights and
non-negative symmetric equilibrium fluxes F_ij = L_ij pi_j = L_ji pi_i, F = exp(theta).
Estimated by maximizing the lag-tau transition log-likelihood sum C(j->i) log [exp(L tau)]_ij
(analytic gradient via the Daleckii-Krein formula).  No matrix logarithm is used, and no edge rate
is ever fitted to electrochemical current.
"""
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.optimize import minimize

FEATURES = ("z", "rg", "r_end", "nc")


def feature_matrix(run):
    return np.column_stack([run[k] for k in FEATURES])


class StateModel:
    def __init__(self, build_runs, k, seed=0, stride=5):
        X = np.vstack([feature_matrix(r)[::stride] for r in build_runs])
        self.mu, self.sd = X.mean(0), X.std(0)
        self.sd[self.sd == 0] = 1
        rng_state = np.random.default_rng(seed)
        best = None
        for rep in range(5):
            c, lab = kmeans2((X - self.mu) / self.sd, k, minit="++", iter=60,
                             seed=int(rng_state.integers(1 << 31)))
            inert = np.sum(((X - self.mu) / self.sd - c[lab]) ** 2)
            if np.bincount(lab, minlength=k).min() > 0 and (best is None or inert < best[0]):
                best = (inert, c)
        self.centres = best[1]
        self.k = k

    def assign(self, run):
        Xs = (feature_matrix(run) - self.mu) / self.sd
        d = ((Xs[:, None, :] - self.centres[None]) ** 2).sum(2)
        return d.argmin(1)


def count_matrix(dtrajs, lag, k):
    C = np.zeros((k, k))
    for s in dtrajs:
        if len(s) > lag:
            np.add.at(C, (s[:-lag], s[lag:]), 1)
    return C            # row convention: C[a, b] = transitions a -> b


def stationary_from(dtrajs, k):
    n = np.zeros(k)
    for s in dtrajs:
        n += np.bincount(s, minlength=k)
    return n / n.sum()


class Generator:
    def __init__(self, pi, F):
        self.pi = pi
        self.F = F
        k = len(pi)
        L = F / pi[None, :]
        np.fill_diagonal(L, 0)
        L -= np.diag(L.sum(0))
        self.L = L

    def S(self):
        sq = np.sqrt(self.pi)
        return (self.L * sq[None, :]) / sq[:, None]

    def propagator(self, t):
        sq = np.sqrt(self.pi)
        lam, V = np.linalg.eigh(0.5 * (self.S() + self.S().T))
        E = (V * np.exp(lam * t)) @ V.T
        return (E * sq[:, None]) / sq[None, :]      # column-stochastic: P[i, j] = prob j -> i

    def timescales(self):
        lam = np.sort(np.linalg.eigvalsh(0.5 * (self.S() + self.S().T)))[::-1]
        return -1.0 / lam[1:]


def fit_generator(C_row, pi, tau, theta0=None, fmin=1e-12):
    k = len(pi)
    iu = np.triu_indices(k, 1)
    sq = np.sqrt(pi)
    Ccol = C_row.T
    ntot = C_row.sum()
    if theta0 is None:
        Fs = 0.5 * (C_row + C_row.T) / max(ntot, 1) / tau
        f0 = np.maximum(Fs[iu], 1e-6 * max(Fs[iu].max(), 1e-30))
        theta0 = np.log(f0)
    lo = np.log(fmin)

    def unpack(th):
        F = np.zeros((k, k))
        F[iu] = np.exp(th)
        return F + F.T

    def nll(th):
        F = unpack(th)
        S = F / np.sqrt(np.outer(pi, pi))
        np.fill_diagonal(S, 0)
        S -= np.diag(F.sum(0) / pi)
        lam, V = np.linalg.eigh(S)
        E = (V * np.exp(lam * tau)) @ V.T
        P = (E * sq[:, None]) / sq[None, :]
        P = np.maximum(P, 1e-300)
        ll = np.sum(Ccol * np.log(P))
        G = Ccol / np.where(np.abs(E) > 1e-300, E, 1e-300)
        dl = lam[:, None] - lam[None, :]
        el = np.exp(lam * tau)
        with np.errstate(divide="ignore", invalid="ignore"):
            Phi = np.where(np.abs(dl) > 1e-12 * (np.abs(lam[:, None]) + 1e-300),
                           (el[:, None] - el[None, :]) / dl, tau * el[:, None])
        H = V @ (Phi * (V.T @ G @ V)) @ V.T
        H = 0.5 * (H + H.T)
        i, j = iu
        dth = np.exp(th) * (2 * H[i, j] / np.sqrt(pi[i] * pi[j]) - H[i, i] / pi[i] - H[j, j] / pi[j])
        return -ll / ntot, -dth / ntot

    res = minimize(nll, theta0, jac=True, method="L-BFGS-B",
                   bounds=[(lo, None)] * len(theta0), options=dict(maxiter=5000, ftol=1e-13, gtol=1e-9))
    gen = Generator(pi, unpack(res.x))
    gen.loglik = -res.fun * ntot
    gen.theta = res.x
    gen.converged = bool(res.success)
    return gen


def empirical_T(C_row):
    rs = C_row.sum(1, keepdims=True)
    return np.where(rs > 0, C_row / np.maximum(rs, 1), 0).T     # column-stochastic


def ck_test(gen, dtrajs_val, lag, mults=(2, 4, 8), main_pi=0.05, nboot=200, rng=None):
    """Compare exp(L k tau) with empirical k-tau transition matrices of independent trajectories.
    Returns per multiple: max |pred - emp| over main initial states (pi >= main_pi) and bootstrap SE."""
    rng = rng or np.random.default_rng(0)
    k = len(gen.pi)
    main = np.where(gen.pi >= main_pi)[0]
    out = []
    # blocks for bootstrap: 10 contiguous blocks per validation trajectory
    blocks = []
    for s in dtrajs_val:
        for b in np.array_split(s, 10):
            blocks.append(b)
    for m in mults:
        t = m * lag
        Cb = np.array([count_matrix([b], t, k) for b in blocks])
        Temp = empirical_T(Cb.sum(0))
        Pred = gen.propagator(t * gen.dt)
        dev = np.abs(Pred[:, main] - Temp[:, main])
        bs = []
        for _ in range(nboot):
            idx = rng.integers(0, len(blocks), len(blocks))
            bs.append(empirical_T(Cb[idx].sum(0))[:, main])
        se = np.std(bs, axis=0)
        imax = np.unravel_index(dev.argmax(), dev.shape)
        out.append(dict(mult=m, max_dev=float(dev.max()), se_at_max=float(se[imax]),
                        z_at_max=float(dev[imax] / max(se[imax], 1e-12))))
    return out


def ess_and_entries(dtrajs, k, lag):
    """Per state: decorrelated effective sample size of occupancy, number of committed entries
    (entry followed by >= lag frames of residence), and seeds that contain committed entries."""
    from md_analysis import tau_int
    ess = np.zeros(k)
    entries = np.zeros(k, int)
    seeds_with = np.zeros(k, int)
    for s in dtrajs:
        for i in range(k):
            ind = (s == i).astype(float)
            n_i = ind.sum()
            if n_i == 0:
                continue
            ti = tau_int(ind, 1.0)
            ess[i] += n_i / (2 * ti)
        # committed entries
        change = np.r_[True, s[1:] != s[:-1]]
        starts = np.where(change)[0]
        lens = np.diff(np.r_[starts, len(s)])
        got = np.zeros(k, bool)
        for st, ln in zip(starts[1:], lens[1:]):
            if ln >= lag:
                entries[s[st]] += 1
                got[s[st]] = True
        seeds_with += got
    return ess, entries, seeds_with


def committed_entries(dtrajs, k, lag):
    """Entries into each state followed by >= lag frames of residence (flicker not counted)."""
    entries = np.zeros(k, int)
    for s in dtrajs:
        change = np.r_[True, s[1:] != s[:-1]]
        starts = np.where(change)[0]
        lens = np.diff(np.r_[starts, len(s)])
        ok = lens[1:] >= lag
        entries += np.bincount(s[starts[1:]][ok], minlength=k)
    return entries
