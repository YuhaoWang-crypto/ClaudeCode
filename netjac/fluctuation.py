"""
fluctuation — metabolite covariance <-> Jacobian (linear-noise approximation).

Forward (Steuer et al. 2003, Bioinformatics 19:1019):  around a stable steady
state with intrinsic noise of diffusion matrix D,
        Jz Cz + Cz Jz^T + Dz = 0                (continuous Lyapunov equation)
For the chemical Langevin equation D = N diag(v) N^T / Omega (Omega = system
size); the factor-2 convention of Sun et al. 2015 (JC + CJ^T = -2D) is the same
equation with D -> D/2.

Inverse (Nägele et al. 2014 PLoS One 9:e96188; COVRECON, Li et al. 2023
Bioinformatics 39:btad397): if C is measured (metabolomics replicates / cells)
and D is fixed by stoichiometry + fluxes, the Lyapunov equation is LINEAR in
the unknown elasticities eps[r,k] on a known sparsity pattern (which species
appear in which rate law) -> regularised least squares recovers the Jacobian,
i.e. the network's EDGES, from snapshot data. Known to be ill-conditioned
(Sun, Länger, Weckwerth 2015); we report recovery error honestly.
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov


def diffusion(model, v, omega=1.0):
    return model.N @ np.diag(np.abs(v)) @ model.N.T / omega


def covariance(model, x, omega=1.0, eps=None):
    """Exact LNA covariance in x-coordinates (rank n - n_cons)."""
    eps = model.elasticities(x) if eps is None else eps
    v = model.v(x)
    Jz = model.B.T @ model.N @ eps @ model.B
    Dz = model.B.T @ diffusion(model, v, omega) @ model.B
    Cz = solve_continuous_lyapunov(Jz, -Dz)
    return model.B @ Cz @ model.B.T, Cz


def sample_lna(model, x, n, omega=1.0, rng=None, eps=None):
    """Draw n independent stationary snapshots (e.g. single cells or
    replicate extracts) from the LNA Gaussian around x*."""
    rng = np.random.default_rng(rng)
    _, Cz = covariance(model, x, omega, eps)
    Lc = np.linalg.cholesky(Cz + 1e-18 * np.eye(len(Cz)))
    z = rng.standard_normal((n, len(Cz))) @ Lc.T
    return x[None, :] + z @ model.B.T


def sparsity(model, x, eps=None, tol=1e-9):
    eps = model.elasticities(x) if eps is None else eps
    scale = np.abs(eps).max()
    return np.abs(eps) > tol * scale


def inverse_elasticities(model, C, v, mask, omega=1.0, ridge=1e-8,
                         prior=None, alpha=0.0):
    """Solve  Jz Cz + Cz Jz^T + Dz = 0  for eps on `mask`, least squares.
    C: measured covariance (x-coords). Returns (eps_hat, condition number).
    prior/alpha: optional Tikhonov pull toward a reference elasticity matrix
    (e.g. the curated healthy model) -> estimates the DIFFERENTIAL Jacobian,
    which is far better conditioned than the absolute one."""
    B, N = model.B, model.N
    Cz = B.T @ C @ B
    Dz = B.T @ diffusion(model, v, omega) @ B
    m = Cz.shape[0]
    iu = np.triu_indices(m)
    A_cols = []
    idx = np.argwhere(mask)
    BN = B.T @ N                          # m x R
    for r, k in idx:
        # d Jz / d eps[r,k] = BN[:, r] outer B[k, :]
        dJ = np.outer(BN[:, r], B[k, :])
        G = dJ @ Cz + Cz @ dJ.T
        A_cols.append(G[iu])
    A = np.array(A_cols).T
    b = -Dz[iu]
    # column scaling + small ridge for conditioning
    s = np.linalg.norm(A, axis=0) + 1e-300
    As = A / s
    k = As.shape[1]
    s0 = np.zeros(k) if prior is None else prior[mask] * s
    lam = ridge + alpha
    sol = np.linalg.solve(As.T @ As + lam * np.eye(k), As.T @ b + lam * s0)
    eps_hat = np.zeros(mask.shape)
    eps_hat[mask] = sol / s
    cond = np.linalg.cond(As)
    return eps_hat, cond
