"""
Brain-hole 3 — Asymmetric partial-view query.

The library is full-transcriptome; the customer has 12 qPCR genes, a targeted
panel, or 4 secreted cytokines. No retraining: the same calibrated loadings W
are sliced to the observed rows and the linear-Gaussian posterior gives
    point estimate  z_hat = M_S r_S            (M_S precomputed = panel LUT)
    uncertainty     C_S = (W_S' D_S^-1 W_S + prior^-1)^-1
Programs whose 90% interval straddles a ternary threshold are MASKED, and
retrieval uses masked Hamming over the remaining bits — an under-determined
query returns "ambiguous between these N entries" instead of a confident
wrong answer.

Also here: a greedy A-optimal panel designer (pick the k genes that best
determine a set of target programs, accounting for nuisance factors).

⚠️ A cytokine-supernatant readout is a different modality; it maps into the
same coordinates only through a matrix fitted on PAIRED samples (here: 150
simulated paired samples). It constrains at most 4 directions.
"""
import numpy as np
from scipy.stats import chi2

from . import codec as C
from . import synth as S


def top_loaded(enc, programs, per=3):
    genes = []
    for p in programs:
        j = S.PROGRAMS.index(p)
        snr = np.abs(enc.W[:, j]) / np.sqrt(enc.D + enc.D_sys)
        genes += [g for g in np.argsort(-snr) if g not in genes][:per]
    return np.array(genes)


def design_panel(enc, programs, k, n_cells, n_cand=400):
    """Greedy A-optimal: minimise summed posterior variance of target programs."""
    tgt = [S.PROGRAMS.index(p) for p in programs]
    cand = np.argsort(-np.abs(enc.W[:, tgt]).max(1))[:n_cand]
    cov = np.eye(S.K_LATENT)
    chosen = []
    for _ in range(k):
        best, best_val, best_cov = None, np.inf, None
        for g in cand:
            if g in chosen:
                continue
            w = enc.W[g]
            cw = cov @ w
            new = cov - np.outer(cw, cw) / (enc.D[g] / n_cells + enc.D_sys[g] + w @ cw)
            val = np.trace(new[np.ix_(tgt, tgt)])
            if val < best_val:
                best, best_val, best_cov = g, val, new
        chosen.append(best)
        cov = best_cov
    return np.array(chosen)


class CytokineAssay:
    """Secreted IFN-g, TNF-a, IL-6, IL-10 (log concentration) as a 4-d readout."""
    PROGS = ["IFNg", "TNF_NFkB", "IL6_JAK_STAT3", "IL10_STAT3"]

    def __init__(self, rng, noise_sd=0.15):     # ~15% CV, typical ELISA / bead assay
        self.true_C = np.zeros((4, S.K_LATENT))
        for i, p in enumerate(self.PROGS):
            self.true_C[i, S.PROGRAMS.index(p)] = rng.uniform(0.8, 1.2)
        self.true_C[1, S.PROGRAMS.index("Inflammasome")] = 0.3
        self.true_C[2, S.PROGRAMS.index("TNF_NFkB")] = 0.3
        self.noise_sd = noise_sd

    def measure(self, z_true, rng):
        return z_true @ self.true_C.T + rng.normal(0, self.noise_sd, (len(z_true), 4))

    def fit_on_paired(self, z_full_hat, y, ridge=1e-2):
        A = z_full_hat
        self.C_hat = np.linalg.solve(A.T @ A + ridge * np.eye(A.shape[1]), A.T @ y).T
        self.D_hat = np.maximum((y - A @ self.C_hat.T).var(0), 1e-3)
        return self

    def view(self):
        return C.LinearGaussianView(self.C_hat, self.D_hat, n_cells=1, name="cytokine")


def sample_states(world, rng, n_samples, n_cells, leaf=None, a=None, mem_on=True):
    """Pseudobulk samples (all cells in a sample share the same state)."""
    A, M, E = S.random_factors(n_samples, rng)
    if not mem_on:
        M[:] = 0
    if a is not None:
        A = a
    leaves = rng.integers(0, len(S.LEAVES), n_samples) if leaf is None else np.full(n_samples, leaf)
    X = np.zeros((n_samples, world.G), np.float32)
    for i in range(n_samples):
        Xi = world.sample_norm(np.full(n_cells, leaves[i]), np.tile(A[i], (n_cells, 1)),
                               np.tile(M[i], (n_cells, 1)), np.tile(E[i], (n_cells, 1)),
                               np.full(n_cells, rng.integers(0, world.n_batches)), rng)
        X[i] = Xi.mean(0)
    Z = np.hstack([A, M, E])
    return X, leaves, Z


def eval_view(view, R_S, Z_true, targets):
    mean, sd = view.posterior(R_S)
    mod = C.SL_MOD
    err = mean[:, mod] - Z_true[:, mod]
    cover = (np.abs(err) <= 1.64 * sd[:, mod]).mean()
    t_hat, t_true = C.ternary(mean[:, mod]), C.ternary(Z_true[:, mod])
    det = C.determined_mask(mean[:, mod], sd[:, mod])
    shrink = view.identifiability()[mod]
    cov_t = view.kappa ** 2 * view.cov[np.ix_(targets, targets)]
    return dict(
        panel=view.name, n_obs=len(view.genes),
        rmse=float(np.sqrt((err ** 2).mean())),
        rmse_targets=float(np.sqrt(((mean[:, targets] - Z_true[:, targets]) ** 2).mean())),
        kappa=view.kappa,
        coverage90=float(cover),
        frac_determined=float(det.mean()),
        tern_acc_determined=float((t_hat == t_true)[det].mean()) if det.any() else float("nan"),
        tern_acc_all=float((t_hat == t_true).mean()),
        n_identifiable=int((shrink < 0.5).sum()),
        cond=view.condition_number(),
        # longest semi-axis of the joint 95% confidence ellipsoid over the target programs
        r95_targets=float(np.sqrt(np.linalg.eigvalsh(cov_t).max() * chi2.ppf(0.95, len(targets)))),
    ), mean, sd
