"""
bioif.real.diffusion -- does a JOINT generative model over the 12 Tox21
endpoints buy anything over 12 independent logistic regressions?

Every predictive model in this repo so far predicts one endpoint at a time.
That is a modelling choice with a cost that is never priced: the 12 Tox21
endpoints are not independent, and the dependence is exactly where the
toxicology is. The sharpest case is the one `tox.py` already names as the
confound: SR-p53 is a genotoxic-stress reporter, SR-MMP is mitochondrial
membrane potential, i.e. a general-cytotoxicity readout. A p53-reporter
positive in a compound that also collapses MMP may mean nothing more than
"the cells were dying". An independent per-endpoint model *structurally
cannot* use a measured MMP value when predicting p53 -- there is no slot for
it. A joint model over the endpoint vector can, by conditioning.

So: fit a conditional DDPM over x0 in R^12 (labels mapped {0,1} -> {-1,+1}),
conditioned on a Morgan count fingerprint, and ask three questions that a
marginal AUROC table cannot answer on its own.

  (a) marginals  -- does per-endpoint predictive quality REGRESS? A joint
                    model that is worse at every single endpoint has bought
                    nothing, whatever its correlation structure looks like.
  (b) joint      -- does it reproduce the correlations BETWEEN endpoints?
  (c) conditional-- P(p53 | fingerprint, MMP=measured). The thing only a
                    joint model can express at all.

Three baselines, because the interesting comparison is not against the
obvious one:

  B1  12 independent logistic regressions. Cannot express (c) at all.
  B2  B1 + a Gaussian copula (multivariate probit) fitted on train-fold
      residuals. This is the baseline that matters: a copula is the cheap
      way to buy joint structure and a closed-form conditional, for ~70
      extra parameters on top of B1. A diffusion model has to beat B2, not
      B1.
  B3  B1's p53 model with the measured MMP label appended as a 2049th
      FEATURE. Also cheap, also able to use MMP -- but only discriminatively,
      for this one ordered pair, and only where MMP happens to be measured.
      A joint model gives every conditional at once; B3 gives one. Included
      because leaving it out would overstate what (c) proves.

VERDICT (measured, held-out scaffolds, seed 0): **the diffusion model does
NOT beat the baselines.** It regresses marginal quality badly (mean AUROC
0.646 vs 0.776 for logistic regression, losing on 11 of 12 endpoints), and
it is the worst of the three at reproducing the endpoint correlation matrix
(Frobenius error 1.65 vs 1.11 for the copula and 1.50 for independent LR).
On axis (c) -- the one axis where a joint model has a structural advantage --
conditioning on the measured MMP does move the diffusion model in the right
direction (+0.019 AUROC), but it starts so far behind that the conditioned
joint model (0.661) is still far below unconditioned logistic regression
(0.776), and the Gaussian copula gets a LARGER lift from the same
information (+0.052, to 0.828) at a tiny fraction of the cost. The
structural argument for joint modelling here is sound; this particular
instrument is the wrong way to cash it in. ~679k denoiser parameters
against ~37k observed training labels is roughly 18 parameters per label,
and the result looks like it.

Protocol, in brief (details at each call site):
  * Scaffold split train/calibration/test = 3726/1491/2236 compounds. Every
    choice -- epochs, diffusion steps T, hidden width, learning rate -- is
    made on the calibration fold. The test fold is scored once.
  * Missingness (7-26% per endpoint) is handled by MASKING the loss to
    observed coordinates and evaluating only observed labels. Missing
    coordinates are resampled from the train-fold prevalence each epoch so
    that the denoiser's INPUT is always a full 12-vector; they are never a
    loss target and never a scoring target.
  * Labels are ✅ when measured on the held-out test fold, ⚠️ when chosen on
    calibration or illustrative.
"""
from __future__ import annotations

import math
import os
import random
import time
from dataclasses import dataclass, field

import numpy as np
import torch
import torch.nn as nn
from scipy.special import ndtr, ndtri

from . import qsar, tox

SEED = 0
N_BITS = 2048
D = len(tox.TOX21_ENDPOINTS)                       # 12
FRACS = (0.5, 0.2, 0.3)
N_SAMPLES_TEST = 200                               # spec floor; see honesty note
N_SAMPLES_CAL = 48                                 # cheaper, calibration only
P53_J = tox.TOX21_ENDPOINTS.index(tox.P53)
MMP_J = tox.TOX21_ENDPOINTS.index(tox.CYTOTOX_CONTROL)

torch.set_num_threads(4)


def _seed_everything(seed: int = SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


# --------------------------------------------------------------------------
# Data: the label MATRIX and its observation mask, kept strictly separate
# --------------------------------------------------------------------------

@dataclass
class Data:
    X: np.ndarray                 # (n, 2048) float32 Morgan counts
    Y: np.ndarray                 # (n, 12)   int8, 0/1; meaningless where ~M
    M: np.ndarray                 # (n, 12)   bool, True = label observed
    smiles: list[str]
    scaffolds: np.ndarray         # (n,) object, Bemis-Murcko scaffold SMILES
    tr: list[int]
    ca: list[int]
    te: list[int]


def load_data(seed: int = SEED) -> Data:
    rows = tox.load_tox21()
    smiles = [r["smiles"] for r in rows]
    X, ok = qsar.featurize(smiles, n_bits=N_BITS)
    keep = np.where(ok)[0]
    rows = [rows[i] for i in keep]
    smiles = [smiles[i] for i in keep]
    X = X[keep]

    raw = np.array([[r[e] for e in tox.TOX21_ENDPOINTS] for r in rows], dtype=object)
    M = ~np.isin(raw, ["", "NA"])
    Y = np.zeros(raw.shape, dtype=np.int8)
    Y[M] = [int(float(v)) for v in raw[M]]

    # Scaffold split, never random: a congeneric dataset split at random
    # reports a number that does not survive new chemistry.
    tr, ca, te = qsar.scaffold_split(smiles, fracs=FRACS, seed=seed)
    scafs = np.array([qsar.scaffold(s) or f"__singleton_{i}"
                      for i, s in enumerate(smiles)], dtype=object)
    return Data(X, Y, M, smiles, scafs, tr, ca, te)


# --------------------------------------------------------------------------
# B1: independent per-endpoint logistic regression
# --------------------------------------------------------------------------

def fit_independent_lr(X, Y, M, idx, C: float = 1.0):
    """One LogisticRegression per endpoint, fitted on its observed rows only."""
    from sklearn.linear_model import LogisticRegression
    models = []
    for j in range(D):
        rows = [i for i in idx if M[i, j]]
        m = LogisticRegression(max_iter=5000, C=C, class_weight="balanced")
        m.fit(X[rows], Y[rows, j])
        models.append(m)
    return models


def lr_proba(models, X) -> np.ndarray:
    return np.column_stack([m.predict_proba(X)[:, 1] for m in models])


def crossfit_proba(X, Y, M, idx, scafs, n_folds: int = 2, seed: int = SEED):
    """
    Out-of-fold B1 probabilities on the TRAIN fold.

    The copula's correlation matrix is estimated from train-fold residuals.
    Using in-sample probabilities would be a self-own: an overfitted marginal
    explains each label too well on its own, the latent thresholds go too
    extreme, and the estimated residual correlation is attenuated towards
    zero -- which would hand the diffusion model an unearned win on axis (b).
    Cross-fitting within the train fold removes that bias. Scaffold-grouped,
    for the same reason the outer split is.
    """
    idx = list(idx)
    uniq = sorted({scafs[i] for i in idx})
    rng = np.random.default_rng(seed)
    assign = {s: int(k) for s, k in zip(uniq, rng.integers(0, n_folds, len(uniq)))}
    fold = np.array([assign[scafs[i]] for i in idx])
    P = np.full((len(idx), D), np.nan)
    for k in range(n_folds):
        fit_rows = [idx[a] for a in np.where(fold != k)[0]]
        pred_rows = np.where(fold == k)[0]
        if not len(pred_rows):
            continue
        models = fit_independent_lr(X, Y, M, fit_rows)
        P[pred_rows] = lr_proba(models, X[[idx[a] for a in pred_rows]])
    # a scaffold fold with no held-out rows cannot happen at this size, but a
    # NaN threshold would silently poison the copula fit, so refuse to guess
    assert np.isfinite(P).all(), "cross-fitting left unfilled rows"
    return P


# --------------------------------------------------------------------------
# B2: Gaussian copula (multivariate probit) over the B1 marginals
# --------------------------------------------------------------------------
# y_j = 1  <=>  z_j > t_j,  z ~ N(0, R),  t_ij = Phi^-1(1 - p_ij).
# The marginals are EXACTLY B1's, so axis (a) is identical by construction;
# R is the only new thing, ~66 free parameters.

_GL_X, _GL_W = np.polynomial.legendre.leggauss(48)
_GL_NODE = 0.5 * (_GL_X + 1.0)
_GL_WT = 0.5 * _GL_W


def bvn_sf(t1, t2, rho: float) -> np.ndarray:
    """
    P(Z1 > t1, Z2 > t2) for a standard bivariate normal with correlation rho.

    Integrated as P(Z1>t1) * E_w[ Phi(-(t2 - rho*Phi^-1(v(w)))/s) ] after the
    substitution v = Phi(u), which puts the quadrature on a FIXED [0,1] grid
    and so vectorises over compounds (each of which has its own thresholds).
    Checked against scipy's multivariate_normal.cdf to 1.3e-5 absolute.
    """
    t1 = np.asarray(t1, float)
    t2 = np.asarray(t2, float)
    s = math.sqrt(max(1e-12, 1.0 - rho * rho))
    lo = ndtr(t1)[..., None]
    v = np.clip(lo + (1.0 - lo) * _GL_NODE, 1e-12, 1 - 1e-12)
    integ = ndtr(-(t2[..., None] - rho * ndtri(v)) / s)
    return ndtr(-t1) * (integ * _GL_WT).sum(-1)


def _thresholds(P: np.ndarray) -> np.ndarray:
    return ndtri(1.0 - np.clip(P, 1e-6, 1 - 1e-6))


def fit_copula_corr(P: np.ndarray, Y: np.ndarray, M: np.ndarray,
                    grid: int = 31, cap: int = 2500, seed: int = SEED):
    """
    Pairwise composite-likelihood estimate of R, pair by pair, then projected
    to the nearest PSD correlation matrix.

    Pairwise (rather than full-likelihood) because with per-compound
    thresholds the full 12-dimensional probit likelihood needs a 12-d orthant
    probability per compound; the pairwise estimator is consistent and is all
    axes (b) and (c) actually use.
    """
    T = _thresholds(P)
    rng = np.random.default_rng(seed)
    rhos = np.linspace(-0.95, 0.95, grid)
    R = np.eye(D)
    for j in range(D):
        for k in range(j + 1, D):
            both = np.where(M[:, j] & M[:, k])[0]
            if len(both) > cap:
                both = rng.choice(both, cap, replace=False)
            if len(both) < 50:
                continue
            t1, t2 = T[both, j], T[both, k]
            y1, y2 = Y[both, j] == 1, Y[both, k] == 1
            q1, q2 = ndtr(-t1), ndtr(-t2)
            ll = np.empty(grid)
            for g, rho in enumerate(rhos):
                p11 = np.clip(bvn_sf(t1, t2, rho), 1e-12, 1.0)
                p10 = np.clip(q1 - p11, 1e-12, 1.0)
                p01 = np.clip(q2 - p11, 1e-12, 1.0)
                p00 = np.clip(1.0 - q1 - q2 + p11, 1e-12, 1.0)
                ll[g] = np.log(np.where(y1 & y2, p11,
                               np.where(y1 & ~y2, p10,
                               np.where(~y1 & y2, p01, p00)))).sum()
            g = int(np.argmax(ll))
            rho = rhos[g]
            if 0 < g < grid - 1:                      # parabolic refinement
                a, b, c = ll[g - 1], ll[g], ll[g + 1]
                den = a - 2 * b + c
                if den < 0:
                    rho += 0.5 * (a - c) / den * (rhos[1] - rhos[0])
            R[j, k] = R[k, j] = float(np.clip(rho, -0.97, 0.97))
    # nearest PSD correlation matrix: clip eigenvalues, renormalise diagonal
    w, V = np.linalg.eigh(R)
    R = V @ np.diag(np.clip(w, 1e-4, None)) @ V.T
    d = np.sqrt(np.diag(R))
    return R / np.outer(d, d)


def copula_samples(P: np.ndarray, R: np.ndarray, n: int, seed: int = SEED):
    """(n, n_compounds, 12) binary draws with B1 marginals and copula R."""
    T = _thresholds(P)
    L = np.linalg.cholesky(R)
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((n, P.shape[0], D)) @ L.T
    return (Z > T[None]).astype(np.int8)


def copula_conditional(P: np.ndarray, R: np.ndarray, j: int, k: int,
                       y_k: np.ndarray) -> np.ndarray:
    """
    P(y_j = 1 | y_k = observed value), closed form. This is the whole point
    of B2: joint structure AND an exact conditional for ~1 extra parameter.
    """
    T = _thresholds(P)
    tj, tk = T[:, j], T[:, k]
    rho = float(R[j, k])
    p11 = bvn_sf(tj, tk, rho)
    qj, qk = ndtr(-tj), ndtr(-tk)
    out = np.where(y_k == 1,
                   p11 / np.clip(qk, 1e-9, None),
                   (qj - p11) / np.clip(1.0 - qk, 1e-9, None))
    return np.clip(out, 0.0, 1.0)


# --------------------------------------------------------------------------
# The conditional DDPM
# --------------------------------------------------------------------------

def cosine_abar(T: int, s: float = 0.008) -> torch.Tensor:
    t = torch.arange(T + 1, dtype=torch.float64) / T
    f = torch.cos((t + s) / (1 + s) * math.pi / 2) ** 2
    return (f / f[0]).float()


@dataclass
class Schedule:
    T: int
    abar: torch.Tensor            # (T+1,) abar[0]=1
    beta: torch.Tensor            # (T,)  beta[i] for step i+1
    alpha: torch.Tensor

    @staticmethod
    def cosine(T: int) -> "Schedule":
        ab = cosine_abar(T)
        beta = (1.0 - ab[1:] / ab[:-1]).clamp(1e-5, 0.999)
        return Schedule(T, ab, beta, 1.0 - beta)


def t_embed(tn: torch.Tensor, dim: int = 64) -> torch.Tensor:
    """
    Sinusoidal embedding of the NORMALISED step tn = t/T in (0, 1], not of the
    integer index t.

    This matters and is easy to get wrong. The cosine schedule defines abar as
    a function of t/T alone, so a model conditioned on t/T sees a consistent
    noise level whatever T is, and the number of diffusion steps becomes a
    pure sampler-resolution knob that can be swept on the calibration fold
    without retraining. A model conditioned on the raw index t would read
    "t=25" as one noise level during T=50 training and a completely different
    one during T=100 sampling -- a silent train/inference mismatch that would
    make the T sweep meaningless.
    """
    half = dim // 2
    freqs = torch.exp(-math.log(10000.0) * torch.arange(half, dtype=torch.float32) / half)
    a = (tn.float() * 1000.0)[:, None] * freqs[None]
    return torch.cat([torch.sin(a), torch.cos(a)], dim=1)


class Denoiser(nn.Module):
    """
    Fingerprint encoder 2048 -> h, then [x_t (12) | t-emb (64) | cond (h)]
    -> h -> h -> 12 predicting epsilon. Deliberately the architecture the
    brief specifies, not a tuned one: the question is whether the FAMILY
    buys anything, and a hand-tuned special case would not answer it.
    """

    def __init__(self, h: int = 256, t_dim: int = 64, n_bits: int = N_BITS):
        super().__init__()
        self.t_dim = t_dim
        self.enc = nn.Sequential(nn.Linear(n_bits, h), nn.SiLU(), nn.LayerNorm(h))
        self.net = nn.Sequential(
            nn.Linear(D + t_dim + h, h), nn.SiLU(),
            nn.Linear(h, h), nn.SiLU(),
            nn.Linear(h, D),
        )

    def cond(self, X: torch.Tensor) -> torch.Tensor:
        return self.enc(X)

    def forward(self, x_t, tn, c):
        """tn is the normalised step t/T in (0, 1], not the integer index."""
        return self.net(torch.cat([x_t, t_embed(tn, self.t_dim), c], dim=1))

    @property
    def n_params(self) -> int:
        return sum(p.numel() for p in self.parameters())


def _features(X: np.ndarray) -> torch.Tensor:
    # log1p on Morgan COUNTS: a linear encoder fed raw counts is badly
    # conditioned. Monotone per-bit, so it cannot add information -- it just
    # stops the encoder from wasting capacity on scale. A choice made in the
    # diffusion model's favour.
    return torch.from_numpy(np.log1p(X))


def train_ddpm(data: Data, T: int, h: int, lr: float, epochs: int,
               batch: int = 256, seed: int = SEED, checkpoints=()):
    """
    Masked-MSE training. Returns {epoch: state_dict copy} at `checkpoints`
    plus the final model, so that "how many epochs" is one calibration-fold
    decision rather than one training run per candidate.
    """
    _seed_everything(seed)
    sch = Schedule.cosine(T)
    model = Denoiser(h=h)
    opt = torch.optim.Adam(model.parameters(), lr=lr)

    idx = np.asarray(data.tr)
    Xt = _features(data.X[idx])
    x0_obs = torch.from_numpy(np.where(data.Y[idx] == 1, 1.0, -1.0).astype(np.float32))
    mask = torch.from_numpy(data.M[idx])
    # train-fold prevalence, used ONLY to fill missing coordinates of the
    # denoiser's INPUT. Never a loss target (the loss is masked), never a
    # scoring target. Resampled every epoch, which is a one-draw stand-in for
    # multiple imputation.
    prev = torch.tensor([float(data.Y[[i for i in idx if data.M[i, j]], j].mean())
                         for j in range(D)], dtype=torch.float32)

    g = torch.Generator().manual_seed(seed)
    snaps: dict[int, dict] = {}
    n = len(idx)
    for ep in range(1, epochs + 1):
        perm = torch.randperm(n, generator=g)
        fill = (torch.rand((n, D), generator=g) < prev[None]).float() * 2 - 1
        x0 = torch.where(mask, x0_obs, fill)
        for b in range(0, n, batch):
            sl = perm[b:b + batch]
            xb, mb, cb = x0[sl], mask[sl], Xt[sl]
            t = torch.randint(1, T + 1, (len(sl),), generator=g)
            ab = sch.abar[t][:, None]
            eps = torch.randn((len(sl), D), generator=g)
            x_t = ab.sqrt() * xb + (1 - ab).sqrt() * eps
            pred = model(x_t, t.float() / T, model.cond(cb))
            loss = (((pred - eps) ** 2) * mb).sum() / mb.sum().clamp(min=1)
            opt.zero_grad()
            loss.backward()
            opt.step()
        if ep in checkpoints:
            snaps[ep] = {k: v.clone() for k, v in model.state_dict().items()}
    snaps[epochs] = {k: v.clone() for k, v in model.state_dict().items()}
    return model, snaps, sch


@torch.no_grad()
def sample_ddpm(model: Denoiser, X: np.ndarray, T: int, n_samples: int,
                clamp_j: int | None = None, clamp_x0: np.ndarray | None = None,
                seed: int = SEED, chunk: int = 200):
    """
    Ancestral sampling. Returns (n_samples, n_compounds, 12) float x0 draws.

    With `clamp_j` set this is inpainting-style conditional sampling: at every
    reverse step the clamped coordinate is overwritten with its KNOWN x0 value
    noised to that timestep, q(x_t | x0), so the denoiser sees a trajectory in
    which that endpoint is already decided and the other 11 are pulled towards
    it through the learned joint. That is the mechanism an independent model
    has no slot for.
    """
    sch = Schedule.cosine(T)
    model.eval()
    g = torch.Generator().manual_seed(seed + 1)
    out = np.empty((n_samples, len(X), D), dtype=np.float32)
    for a in range(0, len(X), chunk):
        b = min(a + chunk, len(X))
        c1 = model.cond(_features(X[a:b]))                 # once per compound
        c = c1.repeat_interleave(n_samples, dim=0)         # (m*S, h)
        m = b - a
        x = torch.randn((m * n_samples, D), generator=g)
        if clamp_j is not None:
            x0c = torch.from_numpy(
                clamp_x0[a:b].astype(np.float32)).repeat_interleave(n_samples)
        for i in range(T, 0, -1):
            tn = torch.full((m * n_samples,), i / T, dtype=torch.float32)
            eps = model(x, tn, c)
            ab, ab_prev = sch.abar[i], sch.abar[i - 1]
            x0h = ((x - (1 - ab).sqrt() * eps) / ab.sqrt()).clamp(-1, 1)
            beta, alpha = sch.beta[i - 1], sch.alpha[i - 1]
            mean = (ab_prev.sqrt() * beta / (1 - ab)) * x0h \
                + (alpha.sqrt() * (1 - ab_prev) / (1 - ab)) * x
            if i > 1:
                var = beta * (1 - ab_prev) / (1 - ab)
                x = mean + var.sqrt() * torch.randn(x.shape, generator=g)
            else:
                x = mean
            if clamp_j is not None:
                if i > 1:
                    x[:, clamp_j] = ab_prev.sqrt() * x0c \
                        + (1 - ab_prev).sqrt() * torch.randn(
                            (m * n_samples,), generator=g)
                else:
                    x[:, clamp_j] = x0c
        out[:, a:b, :] = x.view(m, n_samples, D).permute(1, 0, 2).numpy()
    return out


def samples_to_proba(S: np.ndarray) -> np.ndarray:
    """P(endpoint=1) = mean(sample > 0), as specified."""
    return (S > 0).mean(axis=0)


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------

def marginal_metrics(P: np.ndarray, Y: np.ndarray, M: np.ndarray, idx):
    """
    Per-endpoint (AUROC, AP, n, prevalence), scored ONLY on observed labels.

    `P` is (len(idx), 12) -- predictions for the rows of `idx`, in that order.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score
    rows = np.asarray(idx)
    assert P.shape == (len(rows), D), f"P is {P.shape}, expected {(len(rows), D)}"
    out = []
    for j in range(D):
        keep = M[rows, j]                       # positions within idx
        y, p = Y[rows[keep], j], P[keep, j]
        if len(set(y.tolist())) < 2:
            out.append((float("nan"), float("nan"), len(y), float("nan")))
            continue
        out.append((float(roc_auc_score(y, p)), float(average_precision_score(y, p)),
                    len(y), float(y.mean())))
    return out


def pairwise_corr_binary(B: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    Pairwise-complete correlation among 12 binary columns.

    B is (..., n, 12) -- either the single observed label matrix or a stack of
    S sampled matrices, in which case draws are pooled. The SAME observation
    mask is applied either way, so the observed and implied matrices are
    computed over identical compound sets per pair and the comparison is not
    confounded by which endpoints happen to be co-measured.
    """
    B = B.reshape(-1, B.shape[-2], D) if B.ndim == 3 else B[None]
    R = np.eye(D)
    for j in range(D):
        for k in range(j + 1, D):
            sel = M[:, j] & M[:, k]
            if sel.sum() < 20:
                R[j, k] = R[k, j] = np.nan
                continue
            a = B[:, sel, j].ravel().astype(float)
            b = B[:, sel, k].ravel().astype(float)
            if a.std() < 1e-9 or b.std() < 1e-9:
                R[j, k] = R[k, j] = 0.0
            else:
                R[j, k] = R[k, j] = float(np.corrcoef(a, b)[0, 1])
    return R


def frob(Ra: np.ndarray, Rb: np.ndarray) -> float:
    d = Ra - Rb
    return float(np.sqrt(np.nansum(d * d)))


# --------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------

@dataclass
class Marginals:
    name: str
    per_endpoint: list = field(default_factory=list)   # (auroc, ap, n, prev)

    @property
    def mean_auroc(self) -> float:
        return float(np.nanmean([a for a, *_ in self.per_endpoint]))

    @property
    def mean_ap(self) -> float:
        return float(np.nanmean([p for _, p, *_ in self.per_endpoint]))


@dataclass
class Joint:
    name: str
    frob: float
    rho_p53_mmp: float


@dataclass
class Conditional:
    name: str
    auroc: float
    ap: float
    note: str = ""


@dataclass
class Run:
    seed: int
    n: int
    n_tr: int
    n_cal: int
    n_te: int
    n_obs_train: int
    n_params: int
    wall: float
    choices: dict
    marginals: list = field(default_factory=list)
    joints: list = field(default_factory=list)
    conds: list = field(default_factory=list)
    extras: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# The experiment
# --------------------------------------------------------------------------

CAL_GRID = (
    # (hidden, lr) -- two points, not a search. Width and step size are the
    # two knobs most likely to be the difference between "the family cannot
    # do this" and "I undertrained it", so both get a calibration look.
    (256, 1e-3),
    (128, 3e-3),
)
CAL_EPOCHS = (60, 140, 240)
CAL_T = (20, 50, 100)


def run(seed: int = SEED, verbose: bool = True) -> Run:
    t_start = time.time()
    _seed_everything(seed)

    def say(msg):
        if verbose:
            print(f"  [{time.time() - t_start:6.1f}s] {msg}", flush=True)

    data = load_data(seed)
    X, Y, M = data.X, data.Y, data.M
    tr, ca, te = data.tr, data.ca, data.te
    n_obs_train = int(M[tr].sum())
    say(f"data n={len(X)} train={len(tr)} cal={len(ca)} test={len(te)} "
        f"observed train labels={n_obs_train}")

    # ---- B1 ---------------------------------------------------------------
    b1 = fit_independent_lr(X, Y, M, tr)
    P_te = lr_proba(b1, X[te])
    P_ca = lr_proba(b1, X[ca])
    say("B1 fitted (12 logistic regressions)")

    # ---- B2: copula on cross-fitted train residuals -----------------------
    P_tr_oof = crossfit_proba(X, Y, M, tr, data.scaffolds, n_folds=2, seed=seed)
    R_cop = fit_copula_corr(P_tr_oof, Y[tr], M[tr], seed=seed)
    say("B2 copula correlation fitted on out-of-fold train residuals")

    # ---- diffusion: all choices made on the calibration fold --------------
    best = None
    cal_log = []
    for (h, lr) in CAL_GRID:
        _, snaps, _ = train_ddpm(data, T=50, h=h, lr=lr,
                                 epochs=max(CAL_EPOCHS), seed=seed,
                                 checkpoints=CAL_EPOCHS)
        for ep in CAL_EPOCHS:
            m = Denoiser(h=h)
            m.load_state_dict(snaps[ep])
            S = sample_ddpm(m, X[ca], T=50, n_samples=N_SAMPLES_CAL, seed=seed)
            mm = Marginals("", marginal_metrics(samples_to_proba(S), Y, M, ca))
            cal_log.append({"h": h, "lr": lr, "epochs": ep, "T": 50,
                            "mean_auroc": mm.mean_auroc})
            if best is None or mm.mean_auroc > best["mean_auroc"]:
                best = dict(cal_log[-1], state=snaps[ep])
            say(f"cal h={h} lr={lr} ep={ep} T=50 -> mean AUROC "
                f"{mm.mean_auroc:.4f}")
    for T in CAL_T:
        if T == 50:
            continue
        m = Denoiser(h=best["h"])
        m.load_state_dict(best["state"])
        S = sample_ddpm(m, X[ca], T=T, n_samples=N_SAMPLES_CAL, seed=seed)
        mm = Marginals("", marginal_metrics(samples_to_proba(S), Y, M, ca))
        cal_log.append({"h": best["h"], "lr": best["lr"],
                        "epochs": best["epochs"], "T": T,
                        "mean_auroc": mm.mean_auroc})
        say(f"cal h={best['h']} ep={best['epochs']} T={T} -> mean AUROC "
            f"{mm.mean_auroc:.4f}")
        if mm.mean_auroc > best["mean_auroc"]:
            best = dict(cal_log[-1], state=best["state"])

    # NOTE: the T sweep reuses the T=50-trained weights. The schedule is
    # cosine in t/T, so T is a sampler-resolution knob here, not a retrain.
    model = Denoiser(h=best["h"])
    model.load_state_dict(best["state"])
    T_star = best["T"]
    say(f"chosen on calibration: h={best['h']} lr={best['lr']} "
        f"epochs={best['epochs']} T={T_star}")

    # ======================================================================
    # TEST FOLD -- touched once, from here down
    # ======================================================================
    S_un = sample_ddpm(model, X[te], T=T_star, n_samples=N_SAMPLES_TEST, seed=seed)
    P_dif = samples_to_proba(S_un)
    say(f"test: {N_SAMPLES_TEST} unconditional samples x {len(te)} compounds")

    # ---- axis (a) ---------------------------------------------------------
    marg = [Marginals("B1 independent LR", marginal_metrics(P_te, Y, M, te)),
            Marginals("B2 LR + copula", marginal_metrics(P_te, Y, M, te)),
            Marginals("DDPM joint", marginal_metrics(P_dif, Y, M, te))]

    # ---- axis (b) ---------------------------------------------------------
    Mte, Yte = M[te], Y[te]
    R_obs = pairwise_corr_binary(Yte.astype(float), Mte)
    rng = np.random.default_rng(seed)
    B_b1 = (rng.random((N_SAMPLES_TEST, len(te), D)) < P_te[None]).astype(np.int8)
    B_b2 = copula_samples(P_te, R_cop, N_SAMPLES_TEST, seed=seed)
    B_df = (S_un > 0).astype(np.int8)
    joints = []
    for name, B in (("B1 independent LR", B_b1), ("B2 LR + copula", B_b2),
                    ("DDPM joint", B_df)):
        R = pairwise_corr_binary(B, Mte)
        joints.append(Joint(name, frob(R, R_obs), float(R[P53_J, MMP_J])))
    say("axis (b) correlation matrices done")

    # ---- axis (c): p53 given a MEASURED MMP -------------------------------
    sel = np.where(Mte[:, P53_J] & Mte[:, MMP_J])[0]
    y_p53 = Yte[sel, P53_J]
    y_mmp = Yte[sel, MMP_J]
    te_sel = [te[i] for i in sel]

    def _m(score):
        from sklearn.metrics import average_precision_score, roc_auc_score
        return (float(roc_auc_score(y_p53, score)),
                float(average_precision_score(y_p53, score)))

    conds = [Conditional("B1 LR (cannot see MMP)", *_m(P_te[sel, P53_J]),
                         note="structurally unable to condition")]
    p_cop = copula_conditional(P_te, R_cop, P53_J, MMP_J, Yte[:, MMP_J])
    conds += [
        Conditional("B2 copula, unconditioned", *_m(P_te[sel, P53_J]),
                    note="identical to B1 by construction"),
        Conditional("B2 copula | MMP measured", *_m(p_cop[sel]),
                    note="closed-form Gaussian conditional"),
        Conditional("DDPM, unconditioned", *_m(P_dif[sel, P53_J])),
    ]
    x0_mmp = np.where(Yte[:, MMP_J] == 1, 1.0, -1.0)
    S_c = sample_ddpm(model, X[te_sel], T=T_star, n_samples=N_SAMPLES_TEST,
                      clamp_j=MMP_J, clamp_x0=x0_mmp[sel], seed=seed)
    p_dif_c = samples_to_proba(S_c)[:, P53_J]
    conds.append(Conditional("DDPM | MMP inpainted", *_m(p_dif_c),
                             note="MMP coord overwritten at every step"))
    say("axis (c) conditional sampling done")

    # B3: the cheap discriminative way to use MMP, for scale
    from sklearn.linear_model import LogisticRegression
    tr_sel = [i for i in tr if M[i, P53_J] and M[i, MMP_J]]
    Xa = np.hstack([X, np.where(Y[:, MMP_J] == 1, 1.0, 0.0)
                    .astype(np.float32)[:, None]])
    b3 = LogisticRegression(max_iter=5000, class_weight="balanced")
    b3.fit(Xa[tr_sel], Y[tr_sel, P53_J])
    conds.append(Conditional("B3 LR with MMP as a feature",
                             *_m(b3.predict_proba(Xa[te_sel])[:, 1]),
                             note="not joint; one conditional, trained for it"))

    # Diagnostic: is the DDPM's marginal loss the model or the 200-sample
    # grid? mean(x0) is continuous, so it separates the two.
    mean_x0 = Marginals("DDPM (mean x0 score)",
                        marginal_metrics(S_un.mean(axis=0), Y, M, te))

    return Run(
        seed=seed, n=len(X), n_tr=len(tr), n_cal=len(ca), n_te=len(te),
        n_obs_train=n_obs_train, n_params=model.n_params,
        wall=time.time() - t_start,
        choices={"h": best["h"], "lr": best["lr"], "epochs": best["epochs"],
                 "T": T_star, "cal_log": cal_log,
                 "cal_mean_auroc": best["mean_auroc"]},
        marginals=marg, joints=joints, conds=conds,
        extras={"R_obs": R_obs, "rho_obs_p53_mmp": float(R_obs[P53_J, MMP_J]),
                "n_axis_c": len(sel), "prev_p53_axis_c": float(y_p53.mean()),
                "n_mmp_pos": int((y_mmp == 1).sum()),
                "mean_x0": mean_x0, "copula_R": R_cop,
                "n_obs_test": int(Mte.sum())},
    )


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def report(res: Run | None = None) -> str:
    r = res if res is not None else run()
    ch = r.choices
    ex = r.extras
    by = {m.name: m for m in r.marginals}
    jb = {j.name: j for j in r.joints}
    cb = {c.name: c for c in r.conds}

    L = ["Joint diffusion vs independent logistic regression, 12 Tox21 endpoints",
         "=" * 74,
         f"  seed {r.seed}   wall-clock {r.wall / 60:.1f} min   torch "
         f"{torch.__version__}, 4 threads",
         f"  {r.n} compounds (unique InChIKey skeletons) x {D} endpoints",
         f"  scaffold split: train {r.n_tr} / calibration {r.n_cal} / test "
         f"{r.n_te}",
         f"  observed labels: {r.n_obs_train} train, {ex['n_obs_test']} test "
         f"(the rest are MISSING and are masked everywhere)",
         f"  DDPM denoiser: {r.n_params:,} parameters  ->  "
         f"{r.n_params / r.n_obs_train:.1f} parameters per observed training "
         f"label.",
         "  That ratio is the headline caveat: the model is heavily "
         "over-parameterised",
         "  for this dataset, and nothing below should be read as a statement "
         "about",
         "  diffusion models given 100x more data.",
         "",
         f"  ⚠️  chosen on CALIBRATION: hidden={ch['h']}, lr={ch['lr']}, "
         f"epochs={ch['epochs']}, T={ch['T']}",
         f"      (calibration mean AUROC {ch['cal_mean_auroc']:.4f}; "
         f"{len(ch['cal_log'])} configurations tried)",
         "  ✅  every number in the three tables below is the test fold, "
         "scored once.",
         ""]

    # ---- axis (a) ---------------------------------------------------------
    L += ["(a) MARGINALS -- per-endpoint AUROC on observed test labels   ✅",
          "-" * 74,
          f"  {'endpoint':<14}{'n':>6}{'prev':>7}{'B1/B2':>9}{'DDPM':>9}"
          f"{'delta':>9}"]
    b1m, dfm = by["B1 independent LR"], by["DDPM joint"]
    wins = 0
    for j, e in enumerate(tox.TOX21_ENDPOINTS):
        a1, _, n, pv = b1m.per_endpoint[j]
        a2, _, _, _ = dfm.per_endpoint[j]
        wins += int(a2 > a1)
        L.append(f"  {e:<14}{n:>6}{pv:>7.3f}{a1:>9.3f}{a2:>9.3f}"
                 f"{a2 - a1:>+9.3f}")
    L += [f"  {'MEAN':<14}{'':>6}{'':>7}{b1m.mean_auroc:>9.3f}"
          f"{dfm.mean_auroc:>9.3f}{dfm.mean_auroc - b1m.mean_auroc:>+9.3f}",
          f"  {'mean AP':<14}{'':>6}{'':>7}{b1m.mean_ap:>9.3f}"
          f"{dfm.mean_ap:>9.3f}{dfm.mean_ap - b1m.mean_ap:>+9.3f}",
          "",
          f"  B2's marginals ARE B1's: a copula re-couples the marginals "
          "without",
          "  changing them. So axis (a) is a two-way comparison, not three.",
          f"  DDPM beats B1 on {wins}/{D} endpoints.",
          f"  ⚠️  diagnostic: scoring the DDPM by mean(x0) instead of the "
          f"specified",
          f"      mean(sample>0) gives mean AUROC "
          f"{ex['mean_x0'].mean_auroc:.3f} -- so the "
          f"{'loss is not' if abs(ex['mean_x0'].mean_auroc - dfm.mean_auroc) < 0.02 else 'loss is partly'}"
          f" an artefact of",
          f"      the {N_SAMPLES_TEST}-sample resolution (1/{N_SAMPLES_TEST} "
          "probability grid).",
          ""]

    # ---- axis (b) ---------------------------------------------------------
    L += ["(b) JOINT STRUCTURE -- 12x12 endpoint correlation matrix   ✅",
          "-" * 74,
          "  Pairwise-complete phi correlations on test compounds vs the "
          "matrix each",
          f"  method implies ({N_SAMPLES_TEST} binary draws/compound, same "
          "observation mask).",
          "",
          f"  {'method':<24}{'||R-Robs||_F':>14}{'rho(p53,MMP)':>15}"
          f"{'error':>9}"]
    rho_obs = ex["rho_obs_p53_mmp"]
    for name in ("B1 independent LR", "B2 LR + copula", "DDPM joint"):
        j = jb[name]
        L.append(f"  {name:<24}{j.frob:>14.3f}{j.rho_p53_mmp:>15.3f}"
                 f"{j.rho_p53_mmp - rho_obs:>+9.3f}")
    L += [f"  {'OBSERVED':<24}{0.0:>14.3f}{rho_obs:>15.3f}{0.0:>+9.3f}",
          "",
          "  B1 is not at zero here: conditioning 12 models on the same "
          "fingerprint",
          "  already induces correlation between their outputs. That induced "
          "amount is",
          "  the real bar for axis (b), and it is a much higher bar than "
          "'independent",
          "  models predict independence', which is false.",
          ""]

    # ---- axis (c) ---------------------------------------------------------
    L += [f"(c) THE CONFOUND -- predicting SR-p53 when SR-MMP is known   ✅",
          "-" * 74,
          f"  {ex['n_axis_c']} test compounds with BOTH p53 and MMP observed; "
          f"p53 prevalence {ex['prev_p53_axis_c']:.3f}, "
          f"{ex['n_mmp_pos']} MMP-positive.",
          f"  Observed phi(p53, MMP) on the test fold = {rho_obs:.3f} -- the "
          "confound is real,",
          "  so there IS information in MMP to be had.",
          "",
          f"  {'predictor of p53':<30}{'AUROC':>8}{'AP':>8}  note"]
    for c in r.conds:
        L.append(f"  {c.name:<30}{c.auroc:>8.3f}{c.ap:>8.3f}  {c.note}")
    lift_cop = cb["B2 copula | MMP measured"].auroc - cb["B2 copula, unconditioned"].auroc
    lift_dif = cb["DDPM | MMP inpainted"].auroc - cb["DDPM, unconditioned"].auroc
    L += ["",
          f"  lift from conditioning on the measured MMP:",
          f"    B2 Gaussian copula  {lift_cop:+.3f} AUROC   (closed form, "
          f"~66 extra parameters)",
          f"    DDPM inpainting     {lift_dif:+.3f} AUROC   "
          f"({r.n_params:,} parameters)",
          ""]

    # ---- verdict ----------------------------------------------------------
    d_marg = dfm.mean_auroc - b1m.mean_auroc
    beats_b = jb["DDPM joint"].frob < min(jb["B1 independent LR"].frob,
                                          jb["B2 LR + copula"].frob)
    best_c = max(r.conds, key=lambda c: c.auroc)
    won = (d_marg >= -0.005) and beats_b and \
        best_c.name == "DDPM | MMP inpainted"
    L += ["=" * 74]
    if won:
        L.append("VERDICT: the joint diffusion model BEATS the baselines on all "
                 "three axes.")
    else:
        L += ["VERDICT: the joint diffusion model does NOT beat the baselines. "
              "It loses on",
              f"  (a) marginals  -- mean AUROC {dfm.mean_auroc:.3f} vs "
              f"{b1m.mean_auroc:.3f} for logistic regression "
              f"({d_marg:+.3f}), losing on {D - wins}/{D} endpoints;",
              f"  (b) joint      -- Frobenius error "
              f"{jb['DDPM joint'].frob:.3f} vs "
              f"{jb['B2 LR + copula'].frob:.3f} (copula) and "
              f"{jb['B1 independent LR'].frob:.3f} (independent LR);",
              f"  (c) conditional-- inpainting on the measured MMP does move "
              f"it the right way",
              f"      ({lift_dif:+.3f} AUROC), but from so far back that "
              f"{cb['DDPM | MMP inpainted'].auroc:.3f} is still below",
              f"      UNCONDITIONED logistic regression "
              f"({cb['B1 LR (cannot see MMP)'].auroc:.3f}), and the copula "
              f"extracts more from the",
              f"      same information ({lift_cop:+.3f}, to "
              f"{cb['B2 copula | MMP measured'].auroc:.3f}) for ~1/10,000 of "
              f"the parameters.",
              f"  Best predictor of p53 given a measured MMP on this fold: "
              f"{best_c.name} ({best_c.auroc:.3f}).",
              "",
              "  The structural argument for joint modelling survives this "
              "result: MMP",
              "  really does carry information about p53 "
              f"(phi = {rho_obs:.3f}), and both methods that",
              "  can condition on it gain from doing so. What fails is the "
              "instrument. A",
              "  denoiser with ~18 parameters per observed training label, "
              "learning a",
              "  12-dimensional distribution from ~3.7k compounds, spends its "
              "capacity on",
              "  a generative task nobody asked for and arrives at worse "
              "marginals; the",
              "  copula buys the same conditional structure by adding 66 "
              "numbers to a",
              "  model that already had good marginals.",
              "",
              "  This was not tuned until it won, and should not be. The "
              "calibration",
              "  sweep above is the whole search."]
    L += ["",
          "Scope limits (all of these are single points, not ranges):",
          "  one dataset (Tox21), one target class (nuclear-receptor and "
          "stress-response",
          "  reporter assays, all qHTS), one featurisation (2048-bit radius-2 "
          "Morgan",
          "  counts), one architecture (MLP denoiser, cosine schedule, "
          "ancestral",
          "  sampling), one seed, one scaffold split. Tox21 missingness is "
          "NOT missing",
          "  at random -- compounds are dropped from assays for assay-specific "
          "reasons --",
          "  and every method here treats it as if it were. A negative result "
          "on one",
          "  instrument at this data scale is not a negative result for joint "
          "modelling.",
          "=" * 74]
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
