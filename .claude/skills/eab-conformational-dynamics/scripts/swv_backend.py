"""Shared square-wave-voltammetry current backend for conformation x redox kinetics.

Column-generator convention (L[i, j] = rate j -> i, columns sum to zero):
    dp_R/dt = L p_R + K_red p_O - K_ox p_R
    dp_O/dt = L p_O - K_red p_O + K_ox p_R
K_red = diag(k_loc) e_r(E), K_ox = diag(k_loc) e_o(E),
e_r = exp(-alpha n f eta), e_o = exp((1-alpha) n f eta), eta = E - E0.

Because O and R share L, the total p_O + p_R obeys dp/dt = L p and stays at the equilibrium pi
when it starts there.  Only x = p_O is dynamic:
    dx/dt = (L - g(E) diag(k_loc)) x + e_o(E) diag(k_loc) pi,     g = e_r + e_o.
At fixed E the stationary point is Nernst, x_eq = pi * theta(E), theta = 1/(1 + exp(-n f eta)),
so the deviation u = x - x_eq obeys the homogeneous du/dt = (L - g diag(k_loc)) u, and the
current is I = Q * sum_i k_loc_i [e_r x_i - e_o (pi_i - x_i)] = Q * g * sum_i k_loc_i u_i
(cathodic positive).  Propagating u (not x) avoids catastrophic cancellation in the current.
For a reversible L, v = pi^{-1/2} u gives dv/dt = A v with symmetric A = S - g diag(k_loc),
S = pi^{-1/2} L pi^{1/2}.  -A = G^T G with G = [B ; diag(sqrt(g k_loc))], where B (one row per
edge, sqrt(F_ij) (e_i/sqrt(pi_i) - e_j/sqrt(pi_j)), F_ij = L_ij pi_j) is the square root of -S.
The spectrum is taken from the SVD of G: the slow (reactive) eigenvalue then has relative error
~ eps*sqrt(||A||/|lam|) instead of eps*||A||/|lam|, which matters when conformational exchange is
10^8-10^13 times faster than electron transfer.  A half-pulse is propagated exactly:
    v(h) = W exp(-sigma^2 h) W^T v(0).
On a potential change u jumps by pi (theta_old - theta_new).
The state is carried through every half-pulse and step; nothing is reset to equilibrium.
"""
import numpy as np

F_CONST = 96485.33212
R_GAS = 8.314462618
T_DEFAULT = 298.15
N_ELEC = 2


class Waveform:
    """Osteryoung SWV staircase from e_start towards e_end in steps de, each step = forward half-pulse
    (displaced by e_sw in the scan direction) then backward half-pulse, each 1/(2f) long, currents
    sampled at the end of each half.  Reduction scan (e_end < e_start, e.g. methylene blue): forward at
    E_k - e_sw.  Oxidation scan (e_end > e_start, e.g. ferrocene): forward at E_k + e_sw.
    Currents are always reported cathodic-positive; on an oxidation scan the forward current is negative.
    Potentials are handled in integer mV (exact for 1 mV steps and integer-mV amplitudes)."""

    def __init__(self, f_hz, esw_v, e_start=0.0, e_end=-0.45, de=0.001):
        self.f = float(f_hz)
        self.esw = float(esw_v)
        self.direction = -1 if e_end <= e_start else 1   # a single-point scan keeps the reduction convention
        nstep = int(round(abs(e_start - e_end) / de)) + 1
        self.e_stair = e_start + self.direction * de * np.arange(nstep)
        self.h = 0.5 / self.f
        self.mv_fwd = np.rint((self.e_stair + self.direction * self.esw) * 1000).astype(int)
        self.mv_bwd = np.rint((self.e_stair - self.direction * self.esw) * 1000).astype(int)
        self.e_start = e_start


class Kinetics:
    """Conformational model: equilibrium weights pi (M,), relative local rate w (M,) such that
    k_loc = k_c * w, and a reversible column generator L (M x M, s^-1) or None (no exchange)."""

    def __init__(self, pi, w, L=None, name=""):
        self.pi = np.asarray(pi, float) / np.sum(pi)
        self.w = np.asarray(w, float)
        self.M = len(self.pi)
        self.name = name
        if L is None:
            self.S = np.zeros((self.M, self.M))
            self.B = np.zeros((0, self.M))
        else:
            L = np.asarray(L, float)
            assert np.allclose(L.sum(0), 0, atol=1e-8 * max(1.0, np.abs(L).max())), "columns must sum to 0"
            off = L - np.diag(np.diag(L))
            assert (off >= -1e-12 * max(1.0, np.abs(L).max())).all(), "negative off-diagonal"
            sq = np.sqrt(self.pi)
            S = (L * sq[None, :]) / sq[:, None]
            asym = np.abs(S - S.T).max() / max(1e-300, np.abs(S).max())
            assert asym < 1e-6, f"generator not reversible w.r.t. pi (asym {asym:.2e})"
            self.S = 0.5 * (S + S.T)
            self._build_B(L)

    def _build_B(self, L):
        sq = np.sqrt(self.pi)
        rows = []
        if L is not None:
            Fl = np.asarray(L, float) * self.pi[None, :]
            Fl = 0.5 * (Fl + Fl.T)
            ii, jj = np.where(np.triu(Fl, 1) > 0)
            for i, j in zip(ii, jj):
                r = np.zeros(self.M)
                r[i] = 1.0 / sq[i]
                r[j] = -1.0 / sq[j]
                rows.append(np.sqrt(Fl[i, j]) * r)
        self.B = np.array(rows).reshape(-1, self.M)

    def scaled(self, s):
        k = Kinetics.__new__(Kinetics)
        k.__dict__.update(self.__dict__)
        k.S = self.S * s
        k.B = self.B * np.sqrt(s)
        return k


def _theta(mv, e0, T=T_DEFAULT):
    fn = N_ELEC * F_CONST / (R_GAS * T)
    return 1.0 / (1.0 + np.exp(np.clip(-fn * (mv / 1000.0 - e0), -700, 700)))


def _eig_bank(kin, k_c, alpha, e0, mv_values, T=T_DEFAULT):
    fn = N_ELEC * F_CONST / (R_GAS * T)
    E = mv_values / 1000.0
    eta = E - e0
    er = np.exp(np.clip(-alpha * fn * eta, -700, 700))
    eo = np.exp(np.clip((1 - alpha) * fn * eta, -700, 700))
    kloc = k_c * kin.w
    g = er + eo
    K = len(E)
    nB = kin.B.shape[0]
    G = np.zeros((K, nB + kin.M, kin.M))
    if nB:
        G[:, :nB, :] = kin.B[None]
    G[:, nB + np.arange(kin.M), np.arange(kin.M)] = np.sqrt(g[:, None] * kloc[None, :])
    _, sig, Wt = np.linalg.svd(G, full_matrices=False)
    return dict(g=g, rate=sig ** 2, W=np.transpose(Wt, (0, 2, 1)), kloc=kloc,
                theta=_theta(mv_values, e0, T))


def simulate(kin, wf_list, k_c, alpha, e0, Q=1.0, T=T_DEFAULT):
    """Return list of dicts (one per waveform) with forward, backward, net currents (A) vs E_stair."""
    mv_all = np.unique(np.concatenate([np.r_[w.mv_fwd, w.mv_bwd, int(round(w.e_start * 1000))] for w in wf_list]))
    bank = _eig_bank(kin, k_c, alpha, e0, mv_all, T)
    idx = {int(m): i for i, m in enumerate(mv_all)}
    sq = np.sqrt(kin.pi)
    gk = bank["g"][:, None] * bank["kloc"][None, :] * sq[None, :]   # current weights on v
    out = []
    for wf in wf_list:
        h = wf.h
        decay = np.exp(-bank["rate"] * h)
        th_prev = _theta(np.array([round(wf.e_start * 1000)]), e0, T)[0]
        v = np.zeros(kin.M)          # starts at Nernst equilibrium of the initial potential
        nk = len(wf.e_stair)
        i_f = np.empty(nk)
        i_b = np.empty(nk)
        for k in range(nk):
            for half, mv in ((0, wf.mv_fwd[k]), (1, wf.mv_bwd[k])):
                j = idx[int(mv)]
                v = v + sq * (th_prev - bank["theta"][j])
                th_prev = bank["theta"][j]
                W = bank["W"][j]
                v = W @ (decay[j] * (W.T @ v))
                cur = np.dot(gk[j], v)
                if half == 0:
                    i_f[k] = cur
                else:
                    i_b[k] = cur
        out.append(dict(E=wf.e_stair, fwd=Q * i_f, bwd=Q * i_b, net=Q * (i_f - i_b)))
    return out


# ---------------------------------------------------------------------------------------------
# Route builders
# ---------------------------------------------------------------------------------------------

def fast_eq_kinetics(pi, w):
    """R0: conformations at equilibrium at all times -> one 2-state O/R system, k_eff = <w>_pi."""
    pi = np.asarray(pi, float) / np.sum(pi)
    return Kinetics([1.0], [float(np.dot(pi, w))], None, name="R0")


def smoluchowski_kinetics(z_samples, n_cells, D_nm2_per_s, beta_nm, z_min, z_range=None, name="R1"):
    """R1: finite-volume (square-root approximation) generator on z with no-flux walls.
    pi_i from the histogram of the flexible-model trajectories (pseudo-count 0.5),
    w_i = <exp(-beta (z - z_min))> over the samples in cell i (cell-centre value if empty)."""
    z = np.asarray(z_samples, float)
    lo, hi = (z.min(), z.max()) if z_range is None else z_range
    edges = np.linspace(lo, hi + 1e-9, n_cells + 1)
    dz = edges[1] - edges[0]
    cnt, _ = np.histogram(z, edges)
    pi = (cnt + 0.5) / (cnt + 0.5).sum()
    ib = np.clip(np.digitize(z, edges) - 1, 0, n_cells - 1)
    ez = np.exp(-beta_nm * (z - z_min))
    sums = np.bincount(ib, weights=ez, minlength=n_cells)
    centres = 0.5 * (edges[1:] + edges[:-1])
    w = np.where(cnt > 0, sums / np.maximum(cnt, 1), np.exp(-beta_nm * (centres - z_min)))
    r = D_nm2_per_s / dz ** 2
    L = np.zeros((n_cells, n_cells))
    for i in range(n_cells - 1):
        L[i + 1, i] = r * np.sqrt(pi[i + 1] / pi[i])
        L[i, i + 1] = r * np.sqrt(pi[i] / pi[i + 1])
    L -= np.diag(L.sum(0))
    kin = Kinetics(pi, w, L, name=name)
    kin.centres = centres
    kin.dz = dz
    return kin


def two_population_kinetics(delta):
    """Non-MD control: two non-exchanging static populations (1/2, 1/2) with rates k e^{+-delta}."""
    return Kinetics([0.5, 0.5], [np.exp(delta), np.exp(-delta)], None, name="2pop")


def full_system_expm_check(kin, k_c, alpha, e0, E, h, x0, T=T_DEFAULT):
    """Independent check: propagate the full 2M system (p_O, p_R) with scipy expm."""
    from scipy.linalg import expm
    fn = N_ELEC * F_CONST / (R_GAS * T)
    eta = E - e0
    kl = k_c * kin.w
    Kr = np.diag(kl * np.exp(-alpha * fn * eta))
    Ko = np.diag(kl * np.exp((1 - alpha) * fn * eta))
    sq = np.sqrt(kin.pi)
    L = (kin.S * sq[:, None]) / sq[None, :]   # back to column generator
    G = np.block([[L - Kr, Ko], [Kr, L - Ko]])
    p = np.r_[x0, kin.pi - x0]
    return (expm(G * h) @ p)[: kin.M]
