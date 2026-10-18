"""
LF2 — Landscape–flux decomposition and entropy production on gene circuits.

Theory (J. Wang's landscape–flux framework):
    dx/dt = F(x) + sqrt(2D) ξ(t)
    stationary Fokker–Planck:  ∇·J = 0,   J = F P − D ∇P
    U = −ln P_ss  (the "potential landscape")
    F = −D ∇U + J/P_ss      (gradient part + rotational "flux velocity")
    entropy production rate  σ = ∫ J·D⁻¹·J / P dx  ≥ 0,  = 0 ⇔ detailed balance

Numerics: the Fokker–Planck operator is discretised as a continuous-time
Markov jump process on a grid with Scharfetter–Gummel rates, which obey local
detailed balance exactly (k+/k− = exp(F·h/D)).  The stationary π comes from a
sparse linear solve and σ from the exact Schnakenberg formula
    σ = ½ Σ_ij (π_i k_ij − π_j k_ji) ln(π_i k_ij / π_j k_ji),
which is ≥ 0 by construction and converges to the continuum value.

Cases
  0. gradient control F = −∇V            -> σ ≈ 0 (numerical floor)
  1. mutual activation A⇄B (the "looks like static equilibrium" example)
  2. activator–inhibitor (A ⟶ B ⊣ A)
  3. Brusselator swept through its Hopf point (healthy → critical → oscillating)
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from . import _data


def _B(a):
    """Bernoulli function x/(e^x − 1), stable near 0."""
    a = np.asarray(a, float)
    out = np.empty_like(a)
    small = np.abs(a) < 1e-6
    out[small] = 1 - a[small] / 2
    out[~small] = a[~small] / np.expm1(a[~small])
    return out


class Grid2D:
    def __init__(self, F, lo, hi, n=90, D=0.05):
        self.F, self.n, self.D = F, n, D
        self.x = np.linspace(lo[0], hi[0], n)
        self.y = np.linspace(lo[1], hi[1], n)
        self.hx, self.hy = self.x[1] - self.x[0], self.y[1] - self.y[0]
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing="ij")
        self._build()

    def idx(self, i, j):
        return i * self.n + j

    def _build(self):
        n, D = self.n, self.D
        rows, cols, vals = [], [], []
        self.edges = []  # (i_state, j_state, k_fwd, k_bwd, axis)
        for axis in (0, 1):
            h = self.hx if axis == 0 else self.hy
            if axis == 0:
                I, J = np.meshgrid(np.arange(n - 1), np.arange(n), indexing="ij")
                xm = (self.x[I] + self.x[I + 1]) / 2; ym = self.y[J]
                a_idx = I * n + J; b_idx = (I + 1) * n + J
            else:
                I, J = np.meshgrid(np.arange(n), np.arange(n - 1), indexing="ij")
                xm = self.x[I]; ym = (self.y[J] + self.y[J + 1]) / 2
                a_idx = I * n + J; b_idx = I * n + J + 1
            Fm = self.F(xm, ym)[axis]
            a = Fm * h / D
            kf = D / h ** 2 * _B(-a)
            kb = D / h ** 2 * _B(a)
            a_idx, b_idx, kf, kb = (v.ravel() for v in (a_idx, b_idx, kf, kb))
            self.edges.append((a_idx, b_idx, kf, kb, axis))
            rows += [a_idx, b_idx]; cols += [b_idx, a_idx]; vals += [kf, kb]
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        vals = np.concatenate(vals)
        N = n * n
        K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
        self.Q = K - sp.diags(np.asarray(K.sum(1)).ravel())

    def stationary(self):
        N = self.n ** 2
        A = self.Q.T.tolil()
        A[0, :] = 1.0
        b = np.zeros(N); b[0] = 1.0
        pi = spla.spsolve(A.tocsr(), b)
        pi = np.clip(pi, 1e-300, None); pi /= pi.sum()
        self.pi = pi
        self.P = pi.reshape(self.n, self.n) / (self.hx * self.hy)
        self.U = -np.log(self.P)
        return self

    def epr(self):
        s = 0.0
        for a, b, kf, kb, _ in self.edges:
            f = self.pi[a] * kf; r = self.pi[b] * kb
            m = (f > 0) & (r > 0)
            s += ((f[m] - r[m]) * np.log(f[m] / r[m])).sum()
        return s

    def flux_field(self):
        """Probability current density J at cell centres (average of edges)."""
        n = self.n
        Jx = np.zeros((n, n)); Jy = np.zeros((n, n))
        for a, b, kf, kb, axis in self.edges:
            net = (self.pi[a] * kf - self.pi[b] * kb)
            h_perp = self.hy if axis == 0 else self.hx
            dens = net / h_perp  # probability current per unit length
            ia, ja = np.divmod(a, n); ib, jb = np.divmod(b, n)
            T = Jx if axis == 0 else Jy
            np.add.at(T, (ia, ja), dens / 2); np.add.at(T, (ib, jb), dens / 2)
        self.Jx, self.Jy = Jx, Jy
        return Jx, Jy

    def flux_vs_gradient(self, mask_q=0.5):
        """Cosine between flux velocity V = J/P and ∇U over the high-probability
        region (top `mask_q` of probability mass)."""
        Jx, Jy = self.flux_field()
        Vx, Vy = Jx / self.P, Jy / self.P
        gUx, gUy = np.gradient(self.U, self.hx, self.hy)
        order = np.argsort(-self.P.ravel())
        cum = np.cumsum(self.P.ravel()[order]) * self.hx * self.hy
        m = np.zeros(self.P.size, bool); m[order[cum <= mask_q]] = True
        m = m.reshape(self.P.shape)
        # trim boundary
        m[:2] = m[-2:] = False; m[:, :2] = m[:, -2:] = False
        num = (Vx * gUx + Vy * gUy)[m]
        den = (np.hypot(Vx, Vy) * np.hypot(gUx, gUy))[m] + 1e-30
        w = self.P[m]
        cos_w = (np.abs(num / den) * w).sum() / w.sum()
        # magnitude ratio |V| / |D∇U|
        ratio = (np.hypot(Vx, Vy)[m] * w).sum() / (self.D * np.hypot(gUx, gUy)[m] * w).sum()
        return cos_w, ratio


# ---------------------------------------------------------------- circuits
def hill(z, K=1.0, n=4):
    z = np.maximum(z, 0)
    return z ** n / (K ** n + z ** n)


def F_gradient(x, y):
    # F = −∇V, V = double well in x + harmonic in y + coupling
    dVx = 4 * (x - 1.5) * ((x - 1.5) ** 2 - 0.5) + 0.3 * (y - 1.5)
    dVy = 2 * (y - 1.5) + 0.3 * (x - 1.5)
    return -dVx, -dVy


def F_mutual_activation(x, y, s=0.2, g=2.6):
    return s + g * hill(y) - x, s + g * hill(x) - y


def F_mutual_activation_asym(x, y, s=0.2, g=2.6, kb=0.33):
    # same wiring, but gene B is 3x slower (degradation rate kb)
    return s + g * hill(y) - x, kb * (s + g * hill(x) - y)


def F_activator_inhibitor(x, y, s=0.2, g=2.6):
    return s + g * (1 - hill(y)) - x, s + g * hill(x) - y


def make_brusselator(b, a=1.0):
    def F(x, y):
        return a - (b + 1) * x + x * x * y, b * x - x * x * y
    return F


# ---------------------------------------------------------------- run
def run(make_fig=True, n=90):
    print("=" * 72); print("LF2  Landscape–flux decomposition & entropy production"); print("=" * 72)
    D = 0.05
    cases = {
        "gradient control (F=-∇V)": Grid2D(F_gradient, (0, 0), (3, 3), n, D),
        "mutual activation A⇄B": Grid2D(F_mutual_activation, (0, 0), (4, 4), n, D),
        "mutual activation, B 3x slower": Grid2D(F_mutual_activation_asym, (0, 0), (4, 4), n, D),
        "activator–inhibitor A→B⊣A": Grid2D(F_activator_inhibitor, (0, 0), (3.2, 3.2), n, D),
    }
    res = {}
    for name, g in cases.items():
        g.stationary()
        res[name] = dict(epr=g.epr(), cos=g.flux_vs_gradient()[0],
                         ratio=g.flux_vs_gradient()[1], grid=g)
        print(f"  {name:32s} EPR σ = {res[name]['epr']:10.2e} k_B/time   "
              f"|V|/|D∇U| = {res[name]['ratio']:.3f}")

    # convergence check of σ with grid spacing (mutual activation)
    conv = [(m, Grid2D(F_mutual_activation, (0, 0), (4, 4), m, D).stationary().epr())
            for m in (50, 70, 90, 120)]
    print("  grid convergence (mutual activation): " +
          ", ".join(f"n={m}: σ={s:.4f}" for m, s in conv))

    # Brusselator sweep through Hopf b_c = 1 + a^2 = 2
    bs = np.array([1.2, 1.5, 1.7, 1.85, 1.95, 2.05, 2.2, 2.4])
    eprs = []
    for b in bs:
        g = Grid2D(make_brusselator(b), (0, 0), (4.5, 5.5), n, D).stationary()
        eprs.append(g.epr())
    eprs = np.array(eprs)
    print("  Brusselator (Hopf at b=2.0):")
    for b, e in zip(bs, eprs):
        tag = "  <- just before Hopf" if b == 1.95 else ("  <- just after" if b == 2.05 else "")
        print(f"     b={b:4.2f}  σ={e:8.3f}{tag}")
    # deterministic-limit check: for b<2 the fixed point is a stable focus,
    # rotation already present -> flux before the bifurcation
    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 5, figsize=(23, 4.5))
        Jref = None
        for k, name in enumerate(cases):
            g = res[name]["grid"]
            Jx, Jy = g.flux_field()
            Uc = np.clip(g.U - g.U.min(), 0, 12)
            cf = ax[k].contourf(g.X, g.Y, Uc, 30, cmap="magma")
            mag = np.hypot(Jx, Jy)
            if res[name]["epr"] < 1e-8:
                ax[k].text(0.5, 0.92, "J ≡ 0 (detailed balance)", color="k",
                           transform=ax[k].transAxes, ha="center")
            else:
                s_ = slice(3, None, 5)
                keep = (Uc < 6)[s_, s_]
                scale = mag[s_, s_][keep].max()
                ax[k].quiver(g.X[s_, s_][keep], g.Y[s_, s_][keep],
                             Jx[s_, s_][keep] / scale, Jy[s_, s_][keep] / scale,
                             color="cyan", scale=12, width=0.005)
                ax[k].text(0.02, 0.02, f"max|J| = {scale:.2g}", color="k",
                           transform=ax[k].transAxes, fontsize=8)
            ax[k].set_title(f"{name}\nσ = {res[name]['epr']:.2g}", fontsize=10)
            ax[k].set_xlabel("gene A"); ax[k].set_ylabel("gene B")
        fig.colorbar(cf, ax=ax[3], label="U − U_min (dark = deep valley)")
        ax[4].plot(bs, eprs, "o-", c="crimson"); ax[4].axvline(2.0, ls=":", c="k")
        ax[4].set_xlabel("Brusselator b  (Hopf at 2)"); ax[4].set_ylabel("EPR σ")
        ax[4].set_title("entropy production through a Hopf bifurcation")
        fig.suptitle("Potential U = −ln P_ss (colour) and probability flux J (cyan arrows, "
                     "drawn where U−U_min<6; each panel self-scaled)", y=1.0)
        fig.tight_layout(); fig.savefig(_data.fig("lf2_landscape_flux.png"), dpi=130)
        plt.close(fig)
    return dict(cases=res, conv=conv, brusselator=(bs, eprs))


if __name__ == "__main__":
    run()
