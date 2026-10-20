"""
Brain-hole 1 — Factorized State Barcode (64 bytes per cell or sample).

    byte  0..3   z_ID  tree path   (lineage, subtype, -, -)       uint8
    byte  4..15  z_ID  in-type embedding (12 PCs)                 int8
    byte 16..31  z_Mod "+1" bit-plane  (capacity 128 programs)
    byte 32..47  z_Mod "-1" bit-plane  (ternary: 2 bits / program)
    byte 48..55  z_Mem slow / memory variables                    int8
    byte 56..63  z_Env culture / dissociation artifacts           int8

How the four blocks are kept apart (the hard part):
  1. identity is removed first, by subtracting the reference centroid of the
     predicted leaf (hierarchical nearest-centroid on marker genes);
  2. program, memory and environment coordinates are estimated JOINTLY with
     one linear-Gaussian posterior over all 40 latent coordinates. Because
     the artifact columns are in the same design matrix, a gene shared by
     Dissociation_stress and TNF_NFkB is explained once, not twice. This is
     what a naive per-gene-set score cannot do.
  3. the loadings W are not learned unsupervised (unidentifiable — Locatello
     et al. 2019); each column is measured from a calibration contrast
     (stimulus vs control, warm vs cold dissociation, ...).
"""
import numpy as np

from . import synth as S

CODE_BYTES = 64
MOD_CAP = 128
TAU = 0.5            # ternary threshold, in "one calibration stimulus" units
Q_RANGE = 2.0        # int8 range for Mem / Env, in the same units
BLOCKS = {"ID": (0, 16), "MOD": (16, 48), "MEM": (48, 56), "ENV": (56, 64)}
SL_MOD = slice(S.P)
SL_MEM = slice(S.P, S.P + S.NM)
SL_ENV = slice(S.P + S.NM, S.K_LATENT)


# ---------------------------------------------------------------------------
class LinearGaussianView:
    """Posterior of z given an observed gene subset (the core of brain-hole 3).

        r_S = W_S z + eps,  eps ~ N(0, D_S / n_cells),  z ~ N(0, prior_sd^2)

    Everything that depends only on the panel (S) is precomputed once:
    projection matrix M_S (k x |S|) and posterior covariance C_S (k x k).
    That pair is the panel's lookup table.
    """

    def __init__(self, W, D, genes=None, n_cells=1, prior_sd=1.0, name="full", D_floor=None):
        genes = np.arange(W.shape[0]) if genes is None else np.asarray(genes)
        self.name, self.genes, self.kappa = name, genes, 1.0
        Ws = W[genes]
        noise = D[genes] / n_cells
        if D_floor is not None:          # batch / systematic variance: does not shrink with n
            noise = noise + D_floor[genes]
        Dinv = 1.0 / noise
        F = Ws.T @ (Ws * Dinv[:, None])          # Fisher information
        self.fisher = F
        prior_prec = np.eye(W.shape[1]) / np.square(prior_sd)
        self.cov = np.linalg.inv(F + prior_prec)
        self.M = self.cov @ (Ws.T * Dinv[None, :])
        self.prior_sd = np.broadcast_to(prior_sd, (W.shape[1],)).astype(float)

    def posterior(self, R_S):
        """R_S: (n, |S|) residuals on the panel genes. Returns mean, sd."""
        mean = R_S @ self.M.T
        return mean, self.kappa * np.sqrt(np.diag(self.cov))[None, :].repeat(len(R_S), 0)

    def calibrate(self, R_cal, Z_cal, coords, level=0.9):
        """Split-conformal inflation: scale sd so `level` of calibration errors
        fall inside the interval. Only ever inflates (kappa >= 1)."""
        self.kappa = 1.0
        mean, sd = self.posterior(R_cal)
        ratio = np.abs(mean[:, coords] - Z_cal[:, coords]) / sd[:, coords]
        from scipy.stats import norm
        self.kappa = max(1.0, float(np.quantile(ratio, level) / norm.ppf(0.5 + level / 2)))
        return self.kappa

    def identifiability(self):
        """Per-coordinate shrink ratio posterior_sd / prior_sd (1 = no info)."""
        return np.sqrt(np.diag(self.cov)) / self.prior_sd

    def condition_number(self):
        ev = np.linalg.eigvalsh(self.fisher)
        return np.inf if ev[0] <= 1e-9 * ev[-1] else ev[-1] / ev[0]


# ---------------------------------------------------------------------------
class DCSCEncoder:
    def fit(self, world, rng, n_ref=400, n_cal=200, cal_leaves=(0, 3, 6)):
        """Fit from simulated EXPERIMENTS only (never reads world.W_*)."""
        G, nleaf = world.G, len(S.LEAVES)
        zeros = lambda n, k: np.zeros((n, k))
        # reference atlas: unstimulated, clean-protocol cells per leaf
        cents, resid, batch_means = [], [], []
        for li in range(nleaf):
            leaf = np.full(n_ref, li)
            b = rng.integers(0, world.n_batches, n_ref)
            X = world.sample_norm(leaf, zeros(n_ref, S.P), zeros(n_ref, S.NM),
                                  zeros(n_ref, S.NE), b, rng)
            c = X.mean(0)
            cents.append(c)
            resid.append(X - c)
            batch_means += [X[b == k].mean(0) - c for k in np.unique(b)]
        self.centroid = np.array(cents)
        R = np.vstack(resid)
        self.D = R.var(0) + 1e-3                       # per-gene single-cell noise
        # between-batch variance of reference pseudobulks, minus its sampling part:
        # the floor that pseudobulk averaging cannot remove
        n_per = n_ref / world.n_batches
        self.D_sys = np.maximum(np.var(batch_means, 0) - self.D / n_per, 0) + 1e-5
        self.id_genes = np.arange(S.N_ID_GENES)
        self.lin_centroid = np.array([self.centroid[[i for i, (l, _) in enumerate(S.LEAVES) if l == lin]].mean(0)
                                      for lin in S.LINEAGES])
        # 12-PC in-type identity embedding
        Rid = R[:, self.id_genes]
        _, sv, vt = np.linalg.svd(Rid - Rid.mean(0), full_matrices=False)
        self.id_pcs = vt[:12]
        self.id_scale = sv[:12] / np.sqrt(len(Rid)) * 3 / 127

        # calibration contrasts: one column of W per latent coordinate
        W = np.zeros((G, S.K_LATENT))
        for j in range(S.K_LATENT):
            cols = []
            for li in cal_leaves:
                a, m, e = zeros(n_cal, S.P), zeros(n_cal, S.NM), zeros(n_cal, S.NE)
                if j < S.P:
                    a[:, j] = 1.0
                elif j < S.P + S.NM:
                    m[:, j - S.P] = 1.0
                else:
                    e[:, j - S.P - S.NM] = 1.0
                X = world.sample_norm(np.full(n_cal, li), a, m, e,
                                      rng.integers(0, world.n_batches, n_cal), rng)
                cols.append(X.mean(0) - self.centroid[li])
            W[:, j] = np.mean(cols, 0)
        # zero out pure-noise entries (|loading| below 3 s.e. of the contrast)
        se = np.sqrt(self.D * (1 / n_cal + 1 / n_ref) / len(cal_leaves))
        W[np.abs(W) < 3 * se[:, None]] = 0.0
        self.W = W
        self._views = {}
        return self

    # ---- identity --------------------------------------------------------
    def classify(self, X):
        Xi = X[:, self.id_genes]
        lin = _corr_argmax(Xi, self.lin_centroid[:, self.id_genes])
        leaf = np.empty(len(X), int)
        for l in range(len(S.LINEAGES)):
            idx = np.where(lin == l)[0]
            cand = [i for i, (ln, _) in enumerate(S.LEAVES) if ln == S.LINEAGES[l]]
            if len(idx):
                leaf[idx] = np.array(cand)[_corr_argmax(Xi[idx], self.centroid[cand][:, self.id_genes])]
        return leaf

    def view(self, n_cells=1, genes=None, name="full"):
        key = (n_cells, name)
        if key not in self._views:
            self._views[key] = LinearGaussianView(self.W, self.D, genes, n_cells, name=name,
                                                  D_floor=self.D_sys)
        return self._views[key]

    # ---- encode ------------------------------------------------------------
    def encode(self, X, leaf=None, n_cells=1):
        """X: (n, G) log-normalised rows (cells, or pseudobulk means of n_cells)."""
        leaf = self.classify(X) if leaf is None else np.asarray(leaf)
        R = X - self.centroid[leaf]
        z, sd = self.view(n_cells).posterior(R)
        id_emb = (R[:, self.id_genes] @ self.id_pcs.T) / self.id_scale
        return dict(leaf=leaf, z=z, sd=sd, id_emb=id_emb, code=pack(leaf, id_emb, z))

    def naive_scores(self, X, leaf):
        """Baseline: per-gene-set univariate projection (score_genes style)."""
        R = X - self.centroid[leaf]
        out = np.zeros((len(X), S.K_LATENT))
        for j in range(S.K_LATENT):
            g = np.argsort(-np.abs(self.W[:, j]))[:40]
            w = self.W[g, j]
            out[:, j] = R[:, g] @ w / (w @ w)
        return out


def _corr_argmax(X, C):
    Xc = X - X.mean(1, keepdims=True)
    Cc = C - C.mean(1, keepdims=True)
    Xc /= np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-9
    Cc /= np.linalg.norm(Cc, axis=1, keepdims=True) + 1e-9
    return np.argmax(Xc @ Cc.T, 1)


# ---------------------------------------------------------------------------
def ternary(z_mod, tau=TAU):
    return np.sign(z_mod) * (np.abs(z_mod) > tau)


def determined_mask(mean, sd, tau=TAU, zq=1.64):
    """True where the 90% interval sits inside one ternary bin."""
    lo, hi = mean - zq * sd, mean + zq * sd
    return (lo > tau) | (hi < -tau) | ((lo > -tau) & (hi < tau))


def _q8(v, scale=Q_RANGE):
    return np.clip(np.round(v / scale * 127), -127, 127).astype(np.int8)


def pack(leaf, id_emb, z):
    n = len(leaf)
    code = np.zeros((n, CODE_BYTES), np.uint8)
    code[:, 0] = np.array([S.LINEAGES.index(S.LEAVES[l][0]) for l in leaf]) + 1
    code[:, 1] = np.array([S.TREE[S.LEAVES[l][0]].index(S.LEAVES[l][1]) for l in leaf]) + 1
    code[:, 4:16] = np.clip(np.round(id_emb), -127, 127).astype(np.int8).view(np.uint8)
    t = ternary(z[:, SL_MOD])
    pos = np.zeros((n, MOD_CAP), bool)
    neg = np.zeros((n, MOD_CAP), bool)
    pos[:, :S.P] = t > 0
    neg[:, :S.P] = t < 0
    code[:, 16:32] = np.packbits(pos, axis=1, bitorder="little")
    code[:, 32:48] = np.packbits(neg, axis=1, bitorder="little")
    code[:, 48:48 + S.NM] = _q8(z[:, SL_MEM]).view(np.uint8)
    code[:, 56:56 + S.NE] = _q8(z[:, SL_ENV]).view(np.uint8)
    return code


def mask_bits(mask_mod):
    """Program mask (n, P) -> 32-byte bit mask aligned with the MOD block."""
    m = np.zeros((len(mask_mod), MOD_CAP), bool)
    m[:, :mask_mod.shape[1]] = mask_mod
    b = np.packbits(m, axis=1, bitorder="little")
    return np.hstack([b, b])


def unpack_mod(code):
    pos = np.unpackbits(code[:, 16:32], axis=1, bitorder="little")[:, :S.P].astype(int)
    neg = np.unpackbits(code[:, 32:48], axis=1, bitorder="little")[:, :S.P].astype(int)
    return pos - neg


def describe(code_row):
    """Human-readable barcode, e.g. Myeloid>CD14_Mono | IFNg+ TNF_NFkB+ Oxphos- | ..."""
    lin = S.LINEAGES[code_row[0] - 1]
    sub = S.TREE[lin][code_row[1] - 1]
    t = unpack_mod(code_row[None])[0]
    mods = " ".join(f"{S.PROGRAMS[p]}{'+' if t[p] > 0 else '-'}" for p in np.nonzero(t)[0]) or "quiescent"
    mem = code_row[48:48 + S.NM].view(np.int8) * Q_RANGE / 127
    env = code_row[56:56 + S.NE].view(np.int8) * Q_RANGE / 127
    memtxt = " ".join(f"{n}={v:.1f}" for n, v in zip(S.MEMORY, mem) if abs(v) > TAU) or "-"
    envtxt = " ".join(f"{n}={v:.1f}" for n, v in zip(S.ENV, env) if abs(v) > TAU) or "-"
    return f"{lin}>{sub} | {mods} | mem: {memtxt} | env: {envtxt}"
