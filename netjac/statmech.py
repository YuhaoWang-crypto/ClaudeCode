"""
Micro -> macro: statistical-physics readouts of a population of micro-states.

The dictionary (what plays the role of what):

  micro-state          one sample's / one cell's vector x (proteins, mRNAs)
  ensemble             a cohort of tissue samples, or the control cells of a
                       Perturb-seq screen (unperturbed cells = equilibrium)
  coupling J_ij        the network Jacobian. For a linear(ised) network driven
                       by noise, the stationary covariance C solves
                       J C + C J^T + D = 0 (Lyapunov). For a Gaussian
                       "Boltzmann" ensemble p(x) ~ exp(-x^T C^-1 x / 2) the
                       couplings are the precision C^-1 -> partial
                       correlations.
  order parameter      a collective mode: an eigenvector of C whose eigenvalue
                       stands far above the random-matrix (Marchenko-Pastur /
                       column-shuffle) bulk. Like magnetisation, it is one
                       number per micro-state that summarises thousands of
                       coordinated components.
  effective dimension  participation ratio PR = (sum l)^2 / sum l^2: how many
                       independent macro variables the ensemble really has.
  external field h     a perturbation (drug, knockdown). Linear response
                       dx = -J^-1 B h.
  fluctuation-response If J is symmetric and the noise isotropic (detailed
                       balance), -J^-1 is proportional to C, so the RESPONSE to a
                       field on node j is proportional to column j of the
                       unperturbed COVARIANCE (fluctuation-dissipation).
                       Departures from it measure directed, non-equilibrium
                       regulation.
  frustration          signed couplings that cannot all be satisfied
                       (spin-glass); here, prior-knowledge signs vs data signs.

Two-node theory, used for the feedback-loop test:
  J = [[-a, k12], [k21, -b]] and noise D = diag(d1, d2) give
  C12 = (k12 C22 + k21 C11) / (a + b).
  The same weighted sum k12 C22 + k21 C11 appears when the variation comes from
  parameter differences between individuals (C = R S R^T, R = -J^-1).
  Consequences:
  - positive feedback (k12, k21 of equal sign) adds up: a visible, strong
    correlation of that sign;
  - negative feedback (opposite signs) cancels: a weak, sign-ambiguous
    correlation.
  So homeostatic loops hide themselves in static co-abundance data.
  `ou_two_node` verifies this against the Lyapunov solution.
"""
import numpy as np
from scipy import linalg


# ----------------------------------------------------------------- basics
def zscore(X):
    X = np.asarray(X, dtype=float)
    mu = X.mean(0)
    sd = X.std(0, ddof=1)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def corr_eigs(X):
    """Non-zero eigenvalues of the p x p correlation matrix of X (n x p),
    via the n x n Gram matrix (cheap when p >> n). Sorted descending."""
    Z = zscore(X)
    n = Z.shape[0]
    G = Z @ Z.T / (n - 1)
    lam = np.linalg.eigvalsh(G)[::-1]
    return lam[lam > 1e-10 * lam[0]]


def mp_bounds(n, p):
    """Marchenko-Pastur support of a correlation spectrum of independent
    variables (q = p/n)."""
    q = p / (n - 1)
    return (1 - np.sqrt(q)) ** 2, (1 + np.sqrt(q)) ** 2


def participation_ratio(lam):
    lam = np.asarray(lam, float)
    return lam.sum() ** 2 / (lam ** 2).sum()


def spectral_entropy(lam):
    p = np.asarray(lam, float) / np.sum(lam)
    p = p[p > 0]
    return float(-(p * np.log(p)).sum() / np.log(len(p)))


def shuffle_columns(X, rng):
    """Destroy every coupling, keep every marginal distribution."""
    Y = np.array(X, copy=True)
    for j in range(Y.shape[1]):
        Y[:, j] = Y[rng.permutation(Y.shape[0]), j]
    return Y


def collective_modes(X, n_null=10, rng=None, q=0.95):
    """Count modes above the column-shuffle null and report the macro
    descriptors of the ensemble."""
    rng = np.random.default_rng(0) if rng is None else rng
    lam = corr_eigs(X)
    n, p = X.shape
    null_top = np.array([corr_eigs(shuffle_columns(X, rng))[0]
                         for _ in range(n_null)])
    thr = np.quantile(null_top, q)
    Z = zscore(X)
    return {"lam": lam, "k_null": int((lam > thr).sum()),
            "k_mp": int((lam > mp_bounds(n, p)[1]).sum()),
            "null_top": null_top, "mp_plus": mp_bounds(n, p)[1],
            "PR": participation_ratio(lam), "top1": lam[0] / lam.sum(),
            "entropy": spectral_entropy(lam),
            "mean_abs_r": mean_abs_corr(Z),
            "n": n, "p": p}


def mean_abs_corr(Z, n_pairs=200000, rng=None):
    """Mean |r| over random gene pairs (a 'coupling strength' summary)."""
    rng = np.random.default_rng(1) if rng is None else rng
    n, p = Z.shape
    i = rng.integers(0, p, n_pairs)
    j = rng.integers(0, p, n_pairs)
    m = i != j
    r = (Z[:, i[m]] * Z[:, j[m]]).sum(0) / (n - 1)
    return float(np.abs(r).mean())


def modes(X, k):
    """Top-k collective modes: scores (n x k), loadings (p x k, orthonormal),
    eigenvalues (correlation units)."""
    Z = zscore(X)
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    lam = s ** 2 / (Z.shape[0] - 1)
    return {"scores": U[:, :k] * s[:k], "loadings": Vt[:k].T, "lam": lam[:k]}


def remove_modes(X, k):
    """Subtract the top-k collective modes ('mean-field' / composition
    fields) and return the re-standardised residual."""
    Z = zscore(X)
    if k == 0:
        return Z
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    R = Z - (U[:, :k] * s[:k]) @ Vt[:k]
    return zscore(R)


def subspace_overlap(Va, Vb):
    """Mean squared cosine between two k-dim subspaces (orthonormal columns,
    same row order): 1 = identical, k/p = random expectation."""
    k = min(Va.shape[1], Vb.shape[1])
    Qa, _ = np.linalg.qr(Va)
    Qb, _ = np.linalg.qr(Vb)
    return float(np.linalg.norm(Qa.T @ Qb) ** 2 / k)


# ------------------------------------------------------------- couplings
def corr_pairs(Z, pairs):
    n = Z.shape[0]
    return (Z[:, pairs[:, 0]] * Z[:, pairs[:, 1]]).sum(0) / (n - 1)


class ShrinkPrecision:
    """Ledoit-Wolf shrunk correlation  S* = (1-a) S + a I  and its inverse via
    Woodbury (never forms p x p). partial(pairs) gives the Gaussian-Boltzmann
    couplings  rho_ij = -P_ij / sqrt(P_ii P_jj)."""

    def __init__(self, Z, alpha=None):
        n, p = Z.shape
        if alpha is None:
            from sklearn.covariance import ledoit_wolf_shrinkage
            alpha = float(ledoit_wolf_shrinkage(Z, assume_centered=True,
                                                block_size=2000))
        self.alpha = alpha
        c = alpha                                    # target mu = 1 (corr)
        U = np.sqrt((1 - alpha) / (n - 1)) * Z.T     # p x n
        M = np.linalg.inv(c * np.eye(n) + U.T @ U)
        L = np.linalg.cholesky(M)
        self.W = U @ L                               # P = (I - W W^T) / c
        self.c = c
        self.diag = (1 - (self.W ** 2).sum(1)) / c

    def entries(self, pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        return -(self.W[i] * self.W[j]).sum(1) / self.c

    def partial(self, pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        return -self.entries(pairs) / np.sqrt(self.diag[i] * self.diag[j])


def fisher_dz(ra, na, rb, nb):
    ra = np.clip(ra, -0.999999, 0.999999)
    rb = np.clip(rb, -0.999999, 0.999999)
    return (np.arctanh(ra) - np.arctanh(rb)) / np.sqrt(1 / (na - 3) + 1 / (nb - 3))


def rewiring_scan(Xa, Xb, thr=5.0, n_perm=10, block=1500, keep=2000, rng=None,
                  k_remove=0):
    """All-pairs differential coupling between two ensembles.
    Counts pairs with |Fisher dz| > thr, against a label-permutation null
    (same sample sizes); returns the strongest `keep` pairs.
    k_remove > 0 first subtracts each ensemble's own top-k collective modes
    (inside every permutation too), so that what remains is rewiring of
    the couplings rather than a change of a few macro fields."""
    rng = np.random.default_rng(0) if rng is None else rng
    na, nb = len(Xa), len(Xb)
    Xall = np.vstack([Xa, Xb])

    def scan(A, B, collect):
        Za = remove_modes(A, k_remove).astype(np.float32)
        Zb = remove_modes(B, k_remove).astype(np.float32)
        p = Za.shape[1]
        cnt, best = 0, []
        for s in range(0, p, block):
            e = min(p, s + block)
            ra = Za[:, s:e].T @ Za / (na - 1)
            rb = Zb[:, s:e].T @ Zb / (nb - 1)
            dz = fisher_dz(ra, na, rb, nb)
            ii, jj = np.meshgrid(np.arange(s, e), np.arange(p), indexing="ij")
            up = jj > ii
            m = up & (np.abs(dz) > thr)
            cnt += int(m.sum())
            if collect and m.any():
                best.append(np.column_stack([ii[m], jj[m], dz[m], ra[m], rb[m]]))
        out = np.vstack(best) if best else np.zeros((0, 5))
        if len(out) > keep:
            out = out[np.argsort(-np.abs(out[:, 2]))[:keep]]
        return cnt, out

    n_obs, top = scan(Xa, Xb, True)
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(na + nb)
        c, _ = scan(Xall[perm[:na]], Xall[perm[na:]], False)
        null.append(c)
    p_tot = Xa.shape[1] * (Xa.shape[1] - 1) // 2
    return {"n_obs": n_obs, "null": np.array(null), "n_pairs": p_tot,
            "fdr": float(np.mean(null) / max(n_obs, 1)), "top": top}


# --------------------------------------------- two-node theory (analytic)
def ou_two_node(a, b, k12, k21, d1=1.0, d2=1.0):
    """Stationary covariance of dx = J x dt + noise: Lyapunov solution and the
    closed form C12 = (k12 C22 + k21 C11) / (a + b)."""
    J = np.array([[-a, k12], [k21, -b]])
    D = np.diag([d1, d2])
    stable = bool(np.all(np.linalg.eigvals(J).real < 0))
    if not stable:
        return {"C": None, "C12_closed": np.nan, "r": np.nan, "stable": False}
    C = linalg.solve_continuous_lyapunov(J, -D)
    closed = (k12 * C[1, 1] + k21 * C[0, 0]) / (a + b)
    r = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
    return {"C": C, "C12_closed": closed, "r": r, "stable": True}


def loop_correlation_theory(n=20000, rng=None, gain=(0.1, 2.0)):
    """Random stable 2-node loops: distribution of |r| for positive vs
    negative feedback (equal gains drawn from the same law)."""
    rng = np.random.default_rng(0) if rng is None else rng
    out = {"positive": [], "negative": [], "one-way": []}
    for _ in range(n):
        a, b = rng.uniform(0.5, 2, 2)
        g1, g2 = rng.uniform(*gain, 2)
        d1, d2 = rng.uniform(0.5, 2, 2)
        s = rng.choice([-1, 1])
        for kind, (k12, k21) in {"positive": (s * g1, s * g2),
                                  "negative": (g1, -g2),
                                  "one-way": (s * g1, 0.0)}.items():
            o = ou_two_node(a, b, k12, k21, d1, d2)
            if o["stable"]:
                out[kind].append(o["r"])
    return {k: np.array(v) for k, v in out.items()}


# ------------------------------------------------ perturbation inversion
def mra_local(R, ridge=0.0):
    """Modular response analysis (Kholodenko 2002). R[i, j] = response of
    node i to a perturbation of node j, normalised to node j's own change
    (R[j, j] = 1). Local (direct) response coefficients
        r = -diag(R^-1)^-1 R^-1,   r_ii = -1."""
    n = R.shape[0]
    Ri = np.linalg.solve(R.T @ R + ridge * np.eye(n), R.T) if ridge else np.linalg.inv(R)
    r = -Ri / np.diag(Ri)[:, None]
    return r


def mra_regression(Y, R_S, ridge=1e-2):
    """Direct effects of the perturbed set S on ANY readout gene.
    Y: (perturbations x readouts) responses, R_S: (perturbations x S)
    responses of the perturbed nodes themselves. Solves Y ~ R_S r^T, i.e.
    the readout's response is explained by the changes of the perturbed
    nodes, which removes paths that run through other perturbed nodes.
    Returns r (readouts x S)."""
    A = R_S.T @ R_S + ridge * np.trace(R_S.T @ R_S) / R_S.shape[1] * np.eye(R_S.shape[1])
    return np.linalg.solve(A, R_S.T @ Y).T


def fdt_score(response, cov_col, own):
    """Fluctuation-response agreement for one perturbed node j:
    Pearson r between the response vector and the j-th column of the
    unperturbed covariance (excluding node j itself)."""
    m = np.ones(len(response), bool)
    m[own] = False
    a, b = response[m], cov_col[m]
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])
