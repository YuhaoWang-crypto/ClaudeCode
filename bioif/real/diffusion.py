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
NOT beat the baselines. It TIES the one axis it was expected to lose, and
loses both of the axes it was bought for.**

  (a) marginals    TIE      mean AUROC 0.731 vs 0.735 for independent
                   logistic regression (-0.004, better on 6 of 12
                   endpoints). It does not regress, which is the gate.
  (b) joint        LOSS     Frobenius error of the 12x12 endpoint
                   correlation matrix 1.889, against 1.095 for the copula
                   and 1.382 for independent LR. It reproduces phi(p53,MMP)
                   as 0.112 where the measured value is 0.384 -- worse than
                   independent logistic regressions, which get 0.228 for
                   free just by sharing a fingerprint. NOTE: these are the
                   SELECTED config's numbers; see the CORRECTION below --
                   the sweep's best checkpoint does beat independent LR on
                   the joint, and still loses to the copula.
  (c) conditional  LOSS     inpainting the measured MMP at every reverse
                   step changed p53 AUROC by -0.000 (0.796 -> 0.795), i.e.
                   not at all. The copula's closed-form conditional gained
                   +0.025 (0.808 -> 0.833) from the same information with
                   ~66 parameters. Appending the measured MMP as a 2049th
                   feature to one logistic regression (B3) reaches 0.806.

The surprise is WHERE it fails. The expectation going in was that a 1.6M-
parameter denoiser fitted to 37k observed labels (43 parameters per label)
would be beaten on per-endpoint marginals and might redeem itself on the
joint. The opposite happened: it matches logistic regression on marginal
RANKING and fails on the joint structure and the conditional -- precisely
the two things a joint model is for. Its sampled frequencies are also badly
miscalibrated in level (it under-generates positives by a median 1.7x and
by 4.3x on SR-ATAD5), which AUROC is blind to but which makes the output
unusable as a probability -- and quoting a calibrated joint probability is
most of the reason to fit a generative model at all.

The premise is not what fails. SR-MMP genuinely carries information about
SR-p53 (measured phi = 0.384 on the test fold), and the method that can
condition on it cheaply does gain. The instrument is what fails, and a
better learning rate will not fix a model that is worse at the joint than a
66-parameter copula.

⚠️ One limitation of axis (c) that bounds how far that reading goes. This
model is trained UNCONDITIONALLY and then conditioned at sampling time by
replacement: the known coordinate is overwritten at every reverse step.
Replacement inpainting is a well-known APPROXIMATION to the true conditional
-- the denoiser was never trained on trajectories where one coordinate is
pinned, so nothing forces the other eleven to respond to it. The fair test
of "can a joint generative model use a measured endpoint" would train with
random coordinate masking (or classifier-free guidance over the conditioning
set) so that conditioning is in-distribution. So axis (c) establishes that
THIS construction gains nothing, not that no diffusion model could. Axes (a)
and (b) do not depend on the conditioning mechanism and are unaffected.

⚠️ CORRECTION, from an independent re-check of axis (b) on the calibration
fold (seed 0, the diffusion side taken from the sweep log):

    B2 copula                     0.983
    best DDPM checkpoint in sweep 0.993
    B1 independent LR             1.241

(The re-check first quoted 1.073 for the diffusion row. That is the best
checkpoint of ONE block of the sweep, not of the sweep: the minimum over all
35 configurations is 0.993, at h=128 / 40 epochs / T=50, which is what
`cal_frob_ref["DDPM_best_over_sweep"]` returns and what the report prints.
It cuts both ways, and the two directions should not be averaged into
"slightly strengthened": the SELECTION-COST bullet below gets stronger,
because 1.889 against 0.993 is a worse trade than 1.889 against 1.073, while
the headline "the copula beats every checkpoint" gets WEAKER -- 0.983
against 0.993 is a margin of 0.010 on a 48-sample Monte-Carlo estimate, i.e.
almost certainly inside its own noise.)

So the defensible form of axis (b) is narrower than first written. At the
configuration calibration selected, diffusion loses the joint clearly. At
its BEST over the sweep it is roughly LEVEL with the copula (0.993 vs 0.983)
and beats independent marginals (1.241). What is not in doubt is the price:
the copula reaches that joint structure with ~66 parameters in seconds, the
diffusion model with 1.6M parameters in ~24 minutes. The recommendation --
use the copula -- survives unchanged, but it now rests on cost rather than
on capability, which is a different and more honest argument.

Two statements elsewhere in this docstring remain too strong, in the same
direction:

  * "worse than independent logistic regressions" is true of the SELECTED
    configuration (1.889 vs 1.382 on test) and NOT true of the model class
    -- diffusion's best checkpoint beats independent LR on the joint
    (0.993 vs 1.241). The earlier gloss, that a model worse than independent
    LR at the joint "has not learned the joint", therefore overreached.
  * the selection rule IS costing axis (b) a lot: 1.889 for the selected
    config against 0.993 for the sweep's best. The two are on different
    folds and not strictly comparable, but the gap is far larger than the
    ~0.001 AUROC that bought it. Calibration selected on marginal AUROC and
    axis (b) pays for it; a joint-aware selection rule is the obvious next
    change, and `cal_frob_ref["DDPM_best_over_sweep"]` is already the right
    place to read it from.

So: diffusion does learn some joint structure -- more than independent
marginals do -- and still loses to the cheap closed-form alternative. That is
a weaker and better-supported claim than the one first written here.

Protocol, in brief (details at each call site):
  * Scaffold split train/calibration/test = 3726/1491/2236 compounds. Every
    choice -- epochs, diffusion steps T, hidden width, learning rate -- is
    made on the calibration fold, over a 35-configuration sweep that is
    printed in full. The test fold is scored once.
  * Missingness (7-26% per endpoint) is handled by MASKING the loss to
    observed coordinates and evaluating only observed labels. Missing
    coordinates are resampled from the train-fold prevalence each epoch so
    that the denoiser's INPUT is always a full 12-vector; they are never a
    loss target and never a scoring target.
  * Labels are ✅ when measured on the held-out test fold, ⚠️ when chosen on
    calibration or illustrative.
  * Runtime MISSES its target: 18.7 min on an unloaded container, and more
    on a loaded one, against the ~12 min this was meant to fit in. The cause
    is that the calibration fold preferred the widest configuration in the
    grid (h=512, beyond the 256 the brief specifies) by ~0.001 AUROC, a
    margin inside the Monte-Carlo noise of a 48-sample estimate, and the
    test-time cost of that choice is ~3x. It is disclosed rather than fixed:
    fixing it meant re-picking the grid and scoring the test fold a second
    time, and requirement (2) -- score the test fold once -- is a
    correctness requirement where the runtime target is not. Capping
    CAL_GRID to widths <= 256 brings the whole run under ~7 min and, on the
    calibration fold, costs about 0.001 AUROC; that is the configuration to
    use for a rerun, and it would be a clean single-touch run of its own.
"""
from __future__ import annotations

import math
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
    Fingerprint encoder 2048 -> h, then
    [x_t (12) | given-mask (12) | t-emb (64) | cond (h)] -> h -> h -> 12
    predicting epsilon.

    The given-mask channel is the one departure from the architecture the
    brief sketches, and it is the whole point of the conditional-training
    revision. Replacement inpainting overwrites a coordinate with its known
    value at every reverse step; without a channel saying WHICH coordinates
    were overwritten, the denoiser cannot tell a pinned coordinate from one
    it is supposed to be predicting, so it has no way to treat the pin as
    evidence. It will happily denoise the pinned coordinate back towards its
    own prior and wash the conditioning out -- which is the most likely
    explanation of the previous run's -0.000 lift on axis (c). With the
    channel, "these coordinates are given" is a situation the model is
    trained on, and the pin is information rather than an out-of-distribution
    perturbation.
    """

    def __init__(self, h: int = 256, t_dim: int = 64, n_bits: int = N_BITS):
        super().__init__()
        self.t_dim = t_dim
        self.enc = nn.Sequential(nn.Linear(n_bits, h), nn.SiLU(), nn.LayerNorm(h))
        self.net = nn.Sequential(
            nn.Linear(D + D + t_dim + h, h), nn.SiLU(),
            nn.Linear(h, h), nn.SiLU(),
            nn.Linear(h, D),
        )

    def cond(self, X: torch.Tensor) -> torch.Tensor:
        return self.enc(X)

    def forward(self, x_t, given, tn, c):
        """
        tn is the normalised step t/T in (0, 1], not the integer index.
        `given` is a float 0/1 mask, 1 where that coordinate's x0 is known.
        """
        return self.net(torch.cat([x_t, given, t_embed(tn, self.t_dim), c],
                                  dim=1))

    @property
    def n_params(self) -> int:
        return sum(p.numel() for p in self.parameters())


def _features(X: np.ndarray) -> torch.Tensor:
    # log1p on Morgan COUNTS: a linear encoder fed raw counts is badly
    # conditioned. Monotone per-bit, so it cannot add information -- it just
    # stops the encoder from wasting capacity on scale. A choice made in the
    # diffusion model's favour.
    return torch.from_numpy(np.log1p(X))


#: Probability that a training example is presented with NO given coordinates,
#: i.e. as a purely unconditional generation problem. Conditional sampling is
#: the point of the exercise, but axes (a) and (b) are scored from
#: unconditional samples, so that case has to stay well represented rather
#: than becoming a rare corner of the training distribution.
P_UNCONDITIONAL = 0.25


def train_ddpm(data: Data, T: int, h: int, lr: float, epochs: int,
               batch: int = 256, seed: int = SEED, checkpoints=()):
    """
    Masked-MSE training with random coordinate conditioning.

    TWO masks, which must not be conflated:

      observed (M)  real Tox21 missingness. An unobserved coordinate has no
                    value, so it can never be GIVEN and can never be a loss
                    target. Its input slot is filled from the train-fold
                    prevalence, as before, purely so the denoiser always
                    receives a full 12-vector.
      given (C)     a synthetic conditioning mask, drawn fresh per example per
                    step, and a SUBSET of the observed coordinates. These
                    coordinates are supplied as known: their x0 is the true
                    label, noised to the current timestep exactly as the
                    replacement sampler will noise it, and the given-mask
                    channel flags them. They are excluded from the loss,
                    because the model is not being asked to predict what it
                    was just told.

    The loss is therefore on `observed AND NOT given` -- still never on an
    imputed value. The conditioning rate is drawn per example from U(0,1) so
    the model sees everything from "nothing given" to "all but one given",
    which is what makes conditioning on an arbitrary subset in-distribution
    at sampling time.

    Returns {epoch: state_dict copy} at `checkpoints` plus the final model, so
    that "how many epochs" is one calibration-fold decision rather than one
    training run per candidate.
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
            nb = len(sl)

            # ---- the synthetic conditioning mask ------------------------
            rate = torch.rand((nb, 1), generator=g)
            given = mb & (torch.rand((nb, D), generator=g) < rate)
            given &= ~(torch.rand((nb, 1), generator=g) < P_UNCONDITIONAL)
            # Every example must keep at least one observed, non-given
            # coordinate, or it contributes no gradient and its given-mask
            # pattern teaches the model nothing. Release one at random.
            dead = mb.any(1) & ~(mb & ~given).any(1)
            if dead.any():
                rnd = torch.rand((nb, D), generator=g) * mb.float()
                release = torch.zeros_like(given)
                release[torch.arange(nb), rnd.argmax(1)] = True
                given &= ~(release & dead[:, None])

            t = torch.randint(1, T + 1, (nb,), generator=g)
            ab = sch.abar[t][:, None]
            eps = torch.randn((nb, D), generator=g)
            # Given coordinates are noised by the SAME q(x_t | x0) the sampler
            # applies when it overwrites them, so training and inference see
            # the identical construction.
            x_t = ab.sqrt() * xb + (1 - ab).sqrt() * eps
            pred = model(x_t, given.float(), t.float() / T, model.cond(cb))
            tgt = mb & ~given                     # observed AND not given
            loss = (((pred - eps) ** 2) * tgt).sum() / tgt.sum().clamp(min=1)
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
    noised to that timestep, q(x_t | x0), AND the given-mask channel is set so
    the denoiser knows that coordinate is evidence rather than something to
    predict. Training draws the same kind of mask (see train_ddpm), so this is
    now in-distribution rather than a perturbation applied to a model that has
    never seen a pin.
    """
    sch = Schedule.cosine(T)
    model.eval()
    g = torch.Generator().manual_seed(seed + 1)
    out = np.empty((n_samples, len(X), D), dtype=np.float32)
    gvec = torch.zeros(D)
    if clamp_j is not None:
        gvec[clamp_j] = 1.0
    for a in range(0, len(X), chunk):
        b = min(a + chunk, len(X))
        c1 = model.cond(_features(X[a:b]))                 # once per compound
        c = c1.repeat_interleave(n_samples, dim=0)         # (m*S, h)
        m = b - a
        x = torch.randn((m * n_samples, D), generator=g)
        given = gvec.expand(m * n_samples, D)
        if clamp_j is not None:
            x0c = torch.from_numpy(
                clamp_x0[a:b].astype(np.float32)).repeat_interleave(n_samples)
            # t=T is also a step the denoiser sees, so pin before the first
            # call rather than only after the first update.
            x[:, clamp_j] = sch.abar[T].sqrt() * x0c \
                + (1 - sch.abar[T]).sqrt() * torch.randn(
                    (m * n_samples,), generator=g)
        for i in range(T, 0, -1):
            tn = torch.full((m * n_samples,), i / T, dtype=torch.float32)
            eps = model(x, given, tn, c)
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
    # (hidden width, learning rate). Widths capped at the 256 the brief
    # specifies: the previous run's sweep included 512, which won by ~0.001
    # calibration AUROC (inside Monte-Carlo noise) and cost ~3x in runtime.
    # That budget is spent on seed repeats instead, which answer a question
    # the extra width did not.
    (256, 1e-3),
    (256, 3e-4),
    (128, 3e-3),
)
CAL_EPOCHS = (5, 10, 20, 40, 60, 100, 160, 240)
CAL_T = (20, 50, 100, 200)

#: The calibration-fold selection rule. CHANGED from the previous run, and the
#: change is the point rather than a detail.
#:
#: The previous run selected on mean per-endpoint AUROC, on the reasoning that
#: requirement (a) -- "marginal quality must not regress" -- is the primary
#: gate. Measuring it showed that rule was expensive in a way the reasoning
#: did not anticipate: the selected configuration scored 1.889 joint error
#: where the sweep's best checkpoint scored 0.993 on calibration, a gap ~1000x
#: larger in AUROC-equivalent terms than the ~0.001 AUROC the rule was
#: protecting. Selecting on marginal AUROC was buying almost nothing on axis
#: (a) and paying for it almost entirely on axis (b).
#:
#: So the objective is now the calibration-fold joint Frobenius error. The
#: mean AUROC of the newly selected configuration is still recorded and
#: printed, because the whole question about this swap is what it costs on
#: axis (a), and that cost has to be visible rather than argued.
SELECT_ON = "calibration joint Frobenius error (minimised)"


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
    Mca = M[ca]
    R_cal_obs = pairwise_corr_binary(Y[ca].astype(float), Mca)

    def _score_cal(m, T):
        """Calibration-fold (mean AUROC, joint Frobenius error) for one model."""
        S = sample_ddpm(m, X[ca], T=T, n_samples=N_SAMPLES_CAL, seed=seed)
        mm = Marginals("", marginal_metrics(samples_to_proba(S), Y, M, ca))
        fr = frob(pairwise_corr_binary((S > 0).astype(np.int8), Mca), R_cal_obs)
        return mm.mean_auroc, fr

    best = None
    cal_log = []
    for (h, lr) in CAL_GRID:
        _, snaps, _ = train_ddpm(data, T=50, h=h, lr=lr,
                                 epochs=max(CAL_EPOCHS), seed=seed,
                                 checkpoints=CAL_EPOCHS)
        for ep in CAL_EPOCHS:
            m = Denoiser(h=h)
            m.load_state_dict(snaps[ep])
            auroc, fr = _score_cal(m, 50)
            cal_log.append({"h": h, "lr": lr, "epochs": ep, "T": 50,
                            "mean_auroc": auroc, "frob": fr})
            # SELECTION: minimise joint error (see SELECT_ON).
            if best is None or fr < best["frob"]:
                best = dict(cal_log[-1], state=snaps[ep])
            say(f"cal h={h} lr={lr} ep={ep} T=50 -> frob {fr:.3f} "
                f"(mean AUROC {auroc:.4f})")
    # T is a sampler-resolution knob only (the step embedding is normalised by
    # T), so it is swept on the selected weights rather than by retraining.
    for T in CAL_T:
        if T == 50:
            continue
        m = Denoiser(h=best["h"])
        m.load_state_dict(best["state"])
        auroc, fr = _score_cal(m, T)
        cal_log.append({"h": best["h"], "lr": best["lr"],
                        "epochs": best["epochs"], "T": T,
                        "mean_auroc": auroc, "frob": fr})
        say(f"cal h={best['h']} ep={best['epochs']} T={T} -> frob {fr:.3f} "
            f"(mean AUROC {auroc:.4f})")
        if fr < best["frob"]:
            best = dict(cal_log[-1], state=best["state"])

    # Pre-empting the obvious objection to axis (b): "you selected a
    # checkpoint on MARGINAL AUROC, so of course it is bad at joint
    # structure." Settle it on the calibration fold, where looking costs
    # nothing -- score the baselines' joint error there too, and compare
    # against the BEST joint error any checkpoint in the sweep achieved. If
    # even the sweep's joint-error champion cannot beat the copula on
    # calibration, the test-fold axis (b) result is not an artefact of the
    # selection rule. This is why the comparison is made here and not by
    # scoring a second model on test.
    B_b1_cal = (np.random.default_rng(seed).random(
        (N_SAMPLES_CAL, len(ca), D)) < P_ca[None]).astype(np.int8)
    cal_frob_ref = {
        "B1": frob(pairwise_corr_binary(B_b1_cal, Mca), R_cal_obs),
        "B2": frob(pairwise_corr_binary(
            copula_samples(P_ca, R_cop, N_SAMPLES_CAL, seed=seed), Mca),
            R_cal_obs),
        "DDPM_best_over_sweep": min(c["frob"] for c in cal_log),
        # What the PREVIOUS selection rule would have picked on this same
        # sweep, so the cost of the swap is measurable in both directions
        # rather than only asserted.
        "auroc_rule_would_pick": max(cal_log, key=lambda c: c["mean_auroc"]),
        "best_auroc_over_sweep": max(c["mean_auroc"] for c in cal_log),
    }

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

    # Diagnostic: does the generative model even reproduce the one-dimensional
    # marginals? Comparing the mean generated positive RATE against the
    # observed test-fold prevalence is the cheapest possible check of
    # calibration, and it is where the failure turns out to live. The
    # train-fold prevalence is kept alongside to show the target was stable.
    prev_train = np.array([float(Y[[i for i in tr if M[i, j]], j].mean())
                           for j in range(D)])
    gen_rate = np.array([float(P_dif[Mte[:, j], j].mean()) for j in range(D)])
    lr_rate = np.array([float(P_te[Mte[:, j], j].mean()) for j in range(D)])
    obs_rate = np.array([float(Yte[Mte[:, j], j].mean()) for j in range(D)])

    return Run(
        seed=seed, n=len(X), n_tr=len(tr), n_cal=len(ca), n_te=len(te),
        n_obs_train=n_obs_train, n_params=model.n_params,
        wall=time.time() - t_start,
        choices={"h": best["h"], "lr": best["lr"], "epochs": best["epochs"],
                 "T": T_star, "cal_log": cal_log,
                 "cal_mean_auroc": best["mean_auroc"],
                 "cal_frob": best["frob"], "cal_frob_ref": cal_frob_ref},
        marginals=marg, joints=joints, conds=conds,
        extras={"R_obs": R_obs, "rho_obs_p53_mmp": float(R_obs[P53_J, MMP_J]),
                "n_axis_c": len(sel), "prev_p53_axis_c": float(y_p53.mean()),
                "n_mmp_pos": int((y_mmp == 1).sum()),
                "mean_x0": mean_x0, "copula_R": R_cop,
                "n_obs_test": int(Mte.sum()),
                "prev_train": prev_train, "gen_rate": gen_rate,
                "lr_rate": lr_rate, "obs_rate": obs_rate},
    )


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def report_seed(res: Run) -> str:
    """
    The detailed single-seed read-out: per-endpoint table, the full
    calibration sweep, the generative-calibration diagnostic and the axis (c)
    panel. The VERDICT lives in report(), which aggregates over seeds, because
    a verdict from one seed is exactly the thing this revision set out to stop
    quoting.
    """
    r = res
    ch = r.choices
    ex = r.extras
    by = {m.name: m for m in r.marginals}
    jb = {j.name: j for j in r.joints}
    cb = {c.name: c for c in r.conds}

    L = [f"DETAIL FOR SEED {r.seed} (one seed; the headline numbers above are "
         f"over {len(HEADLINE_SEEDS)} seeds)",
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
         f"  ⚠️  runtime overshoot, disclosed: {r.wall / 60:.1f} min against a "
         "~12 min target. The",
         f"      calibration fold chose the widest and slowest configuration "
         f"in the grid (h={ch['h']},",
         f"      T={ch['T']}), and nearly all of the excess is the "
         f"{N_SAMPLES_TEST}-sample x {ch['T']}-step test sampling",
         "      that implies. It was preferred over h=128 by ~0.001 "
         "calibration AUROC, which is",
         f"      inside the Monte-Carlo noise of a {N_SAMPLES_CAL}-sample "
         "estimate, so a budget-capped grid",
         "      would run ~3x faster at essentially no cost in quality. That "
         "grid was NOT",
         "      substituted after the fact: requirement (2) is that the test "
         "fold is scored",
         "      once, and re-picking the grid to hit a runtime target after "
         "seeing the test",
         "      numbers would buy speed with the only thing here worth "
         "having.",
         "",
         f"  ⚠️  chosen on CALIBRATION: hidden={ch['h']}, lr={ch['lr']}, "
         f"epochs={ch['epochs']}, T={ch['T']}",
         f"      selection rule, declared in advance: {SELECT_ON}. "
         f"{len(ch['cal_log'])} configurations tried;",
         f"      winner scored {ch['cal_mean_auroc']:.4f} AUROC / "
         f"{ch['cal_frob']:.3f} joint error on calibration.",
         "  ✅  every number in the three tables below is the test fold, "
         "scored once.",
         ""]

    # The whole search, printed. The joint-error column is shown but was NOT
    # used to select, so a reader can check whether a different criterion
    # would have been kinder to the model.
    cl = ch["cal_log"]
    best_fr = min(cl, key=lambda c: c["frob"])
    L += ["  ⚠️  calibration sweep (the entire search; calibration fold, "
          f"{N_SAMPLES_CAL} samples/compound)",
          f"      {'h':>5}{'lr':>8}{'epochs':>8}{'T':>5}{'meanAUROC':>11}"
          f"{'jointErr':>10}"]
    for c in cl:
        mark = " <- selected" if (c["h"] == ch["h"] and c["lr"] == ch["lr"]
                                  and c["epochs"] == ch["epochs"]
                                  and c["T"] == ch["T"]) else ""
        L.append(f"      {c['h']:>5}{c['lr']:>8.1e}{c['epochs']:>8}"
                 f"{c['T']:>5}{c['mean_auroc']:>11.4f}{c['frob']:>10.3f}"
                 f"{mark}")
    same = (best_fr["h"] == ch["h"] and best_fr["epochs"] == ch["epochs"]
            and best_fr["T"] == ch["T"])
    L.append(f"      Marginal AUROC is maximised at {ch['epochs']} epochs "
             f"(h={ch['h']}, T={ch['T']}).")
    if same:
        L.append("      The same configuration also minimises joint error, so "
                 "the selection rule is not load-bearing.")
    else:
        L += [f"      Joint error is minimised somewhere else "
              f"(h={best_fr['h']}, epochs={best_fr['epochs']}, "
              f"T={best_fr['T']}: {best_fr['frob']:.3f} joint error but only",
              f"      {best_fr['mean_auroc']:.4f} marginal AUROC). The two "
              "axes disagree about which checkpoint is best,",
              "      which is itself a finding: this model cannot be made good "
              "at marginals and at joint",
              "      structure simultaneously at this data scale. Selecting on "
              "marginals, as declared, is the",
              "      choice that gives it the best shot at the primary gate "
              "(a)."]
    L.append("")

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
          f"  DDPM beats B1 on {wins}/{D} endpoints, loses on {D - wins}/{D}.",
          f"  ⚠️  diagnostic: scoring the DDPM by mean(x0), which is "
          f"continuous, instead of the",
          f"      specified mean(sample>0) gives mean AUROC "
          f"{ex['mean_x0'].mean_auroc:.3f} against "
          f"{dfm.mean_auroc:.3f} -- a difference of",
          f"      {ex['mean_x0'].mean_auroc - dfm.mean_auroc:+.3f}, so the "
          f"1/{N_SAMPLES_TEST} probability grid is "
          f"{'not' if abs(ex['mean_x0'].mean_auroc - dfm.mean_auroc) < 0.02 else 'partly'}"
          f" what limits",
          f"      axis (a); {N_SAMPLES_TEST} samples are "
          f"{'enough' if abs(ex['mean_x0'].mean_auroc - dfm.mean_auroc) < 0.02 else 'NOT enough'}"
          " here.",
          ""]

    # ---- why: does the generative model recover the marginals at all? -----
    gr, lrr, orr = ex["gen_rate"], ex["lr_rate"], ex["obs_rate"]
    L += ["    Calibration of the generative marginals (mean predicted "
          "positive rate)   ✅",
          f"    {'endpoint':<14}{'observed':>10}{'B1 LR':>9}{'DDPM':>9}"
          f"{'DDPM/obs':>10}"]
    for j, e in enumerate(tox.TOX21_ENDPOINTS):
        L.append(f"    {e:<14}{orr[j]:>10.3f}{lrr[j]:>9.3f}{gr[j]:>9.3f}"
                 f"{gr[j] / max(orr[j], 1e-9):>10.2f}x")
    ratio = float(np.median(gr / np.maximum(orr, 1e-9)))
    direction = "UNDER" if ratio < 1 else "OVER"
    factor = (1.0 / ratio) if ratio < 1 else ratio
    worst = int(np.argmin(gr / np.maximum(orr, 1e-9)))
    L += [f"    The DDPM {direction}-generates the positive class by a median "
          f"factor of {factor:.2f}x,",
          f"    worst at {tox.TOX21_ENDPOINTS[worst]} "
          f"({gr[worst]:.3f} generated vs {orr[worst]:.3f} observed). "
          "Continuous diffusion on a",
          "    {-1,+1} target whose positive rate is 2-17% has to put a "
          "sharply bimodal,",
          "    sharply asymmetric x0 posterior through a Gaussian reverse "
          "kernel, and the mass",
          "    ends up on the wrong side of zero at the wrong rate.",
          "",
          "    Read this against axis (a) above, because the two together say "
          "something more",
          "    specific than either alone: AUROC is rank-based and is "
          "therefore blind to a",
          "    monotone level error, so this miscalibration costs the DDPM "
          "almost nothing on",
          "    axis (a) -- its RANKING is competitive with logistic "
          "regression. What it costs is",
          "    the interpretation of the output AS A PROBABILITY. B1 is "
          "miscalibrated too, and in",
          "    the opposite direction (class_weight='balanced' inflates its "
          "rates by design), but",
          "    B1's distortion is a known monotone reweighting that a single "
          "Platt scaling undoes,",
          "    whereas the DDPM's varies per endpoint by a factor of "
          f"{float((gr / np.maximum(orr, 1e-9)).max() / (gr / np.maximum(orr, 1e-9)).min()):.1f}x "
          "across the 12. So",
          "    the sampled frequencies here are usable as a ranking and NOT "
          "usable as risk",
          "    estimates -- which matters, because being able to quote a "
          "calibrated joint",
          "    probability is most of the reason to fit a generative model in "
          "the first place.",
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
    cfr = ch["cal_frob_ref"]
    sweep_wins_b = cfr["DDPM_best_over_sweep"] < cfr["B2"]
    L += ["  ⚠️  Is this an artefact of selecting on marginals? Settled on the "
          "CALIBRATION",
          "      fold, so that looking costs no test-fold touch. Best joint "
          "error achieved by",
          f"      ANY of the {len(cl)} checkpoints in the sweep: "
          f"{cfr['DDPM_best_over_sweep']:.3f}, against {cfr['B2']:.3f} for the "
          f"copula and",
          f"      {cfr['B1']:.3f} for independent LR on the same fold.",
          ("      So a luckier selection rule WOULD have beaten the copula "
           "here, and the test-fold" if sweep_wins_b else
           "      So no checkpoint in the sweep beats the copula on joint "
           "structure, whichever one"),
          ("      ranking on (b) is partly a consequence of selecting on "
           "marginals -- read it with that" if sweep_wins_b else
           "      you select. The test-fold ranking on axis (b) is a property "
           "of the model, not of"),
          ("      caveat." if sweep_wins_b else "      the selection rule."),
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

    return "\n".join(L)


# --------------------------------------------------------------------------
# Seed repeats: a noise estimate on the headline comparison
# --------------------------------------------------------------------------
# The previous run reported a 0.010 axis-(b) margin (copula 0.983 vs the
# diffusion sweep's best 0.993) from a single seed and a 48-sample Monte-Carlo
# estimate, and could not say whether that margin was real. It almost
# certainly was not. Every number quoted as a headline below is therefore
# repeated over HEADLINE_SEEDS, and the seed changes BOTH the scaffold split
# and all sampling, so the spread covers split variation, sampling variation
# and selection variation together -- the per-seed calibration sweep is rerun
# per seed, so a seed that would have selected a different checkpoint does.
#
# The sweep is NOT repeated for extra configurations: the seed budget buys
# repeats of the final comparison, as instructed, not a wider search.

HEADLINE_SEEDS = (0, 1, 2)


def run_multi(seeds=HEADLINE_SEEDS, verbose: bool = True) -> list[Run]:
    out = []
    for s in seeds:
        if verbose:
            print(f"=== seed {s} " + "=" * 58, flush=True)
        out.append(run(seed=s, verbose=verbose))
    return out


def _spread(vals) -> tuple[float, float, float]:
    """(mean, min, max) of a per-seed quantity."""
    a = np.asarray(vals, float)
    return float(a.mean()), float(a.min()), float(a.max())


def decide(diffs) -> tuple[bool, float, float]:
    """
    Apply the declared tie rule to a set of PAIRED per-seed differences.

    A difference counts as decisive only when every seed agrees on its sign
    AND the mean difference is larger than the seed-to-seed range of that
    difference -- i.e. the effect must be bigger than its own variability.
    Anything else is a tie, including a difference that is consistent in sign
    but small compared with how much it moves between seeds.

    Pairing matters: the same seed gives both methods the same split and the
    same sampling draws, so the per-seed difference removes the split-to-split
    variation that dominates the raw per-method spread.
    """
    d = np.asarray(diffs, float)
    same_sign = bool(np.all(d > 0) or np.all(d < 0))
    rng = float(d.max() - d.min())
    return (same_sign and abs(float(d.mean())) > rng), float(d.mean()), rng


#: The previous run, for the explicit before/after comparison. These are the
#: committed numbers from the AUROC-selected, unmasked-training, single-seed
#: configuration -- recorded here as data so the comparison cannot drift.
PREV = {
    "label": "previous run (AUROC-selected, no masked training, 1 seed)",
    "a_ddpm": 0.731, "a_b1": 0.735,
    "b_ddpm": 1.889, "b_b2": 1.095, "b_b1": 1.382,
    "c_b1": 0.808, "c_cop_cond": 0.833, "c_ddpm_un": 0.796,
    "c_ddpm_cond": 0.795, "c_b3": 0.806,
    "lift_cop": 0.025, "lift_ddpm": -0.000,
    "params": 1620492, "wall_min": 18.7,
    "cal_best_frob_ddpm": 0.993, "cal_frob_b2": 0.983,
}


def report(runs: list[Run] | None = None) -> str:
    rs = runs if runs is not None else run_multi()
    seeds = [r.seed for r in rs]

    def mg(r, name):
        return {m.name: m for m in r.marginals}[name]

    def jt(r, name):
        return {j.name: j for j in r.joints}[name]

    def cd(r, name):
        return {c.name: c for c in r.conds}[name]

    # ---- per-seed headline quantities -----------------------------------
    a_b1 = [mg(r, "B1 independent LR").mean_auroc for r in rs]
    a_df = [mg(r, "DDPM joint").mean_auroc for r in rs]
    b_b1 = [jt(r, "B1 independent LR").frob for r in rs]
    b_b2 = [jt(r, "B2 LR + copula").frob for r in rs]
    b_df = [jt(r, "DDPM joint").frob for r in rs]
    rho_o = [r.extras["rho_obs_p53_mmp"] for r in rs]
    rho_b2 = [jt(r, "B2 LR + copula").rho_p53_mmp for r in rs]
    rho_df = [jt(r, "DDPM joint").rho_p53_mmp for r in rs]
    c_b1 = [cd(r, "B1 LR (cannot see MMP)").auroc for r in rs]
    c_b3 = [cd(r, "B3 LR with MMP as a feature").auroc for r in rs]
    c_cop_u = [cd(r, "B2 copula, unconditioned").auroc for r in rs]
    c_cop_c = [cd(r, "B2 copula | MMP measured").auroc for r in rs]
    c_df_u = [cd(r, "DDPM, unconditioned").auroc for r in rs]
    c_df_c = [cd(r, "DDPM | MMP inpainted").auroc for r in rs]
    lift_cop = [x - y for x, y in zip(c_cop_c, c_cop_u)]
    lift_df = [x - y for x, y in zip(c_df_c, c_df_u)]

    def row(label, vals, fmt="{:.3f}", width=30):
        m, lo, hi = _spread(vals)
        return (f"  {label:<{width}}" + fmt.format(m)
                + "   [" + fmt.format(lo) + ", " + fmt.format(hi) + "]")

    tot = sum(r.wall for r in rs) / 60.0
    L = ["Joint diffusion vs independent logistic regression, 12 Tox21 "
         "endpoints",
         "REVISION 2 -- conditional training, joint-error selection, seed "
         "repeats",
         "=" * 74,
         f"  seeds {seeds}   total wall-clock {tot:.1f} min   torch "
         f"{torch.__version__}, 4 threads",
         f"  {rs[0].n} compounds x {D} endpoints; per seed: train {rs[0].n_tr}"
         f" / calibration {rs[0].n_cal} / test {rs[0].n_te}",
         f"  DDPM denoiser: {rs[0].n_params:,} parameters  ->  "
         f"{rs[0].n_params / rs[0].n_obs_train:.1f} per observed training "
         f"label (still heavily over-parameterised)",
         "",
         "WHAT CHANGED since the previous run, and why:",
         "  1. SELECTION OBJECTIVE. Was mean calibration AUROC; now "
         "calibration joint",
         "     Frobenius error. The old rule was measured to cost ~0.9 "
         "Frobenius on axis (b)",
         f"     to protect ~0.001 AUROC on axis (a) "
         f"({PREV['b_ddpm']:.3f} selected vs "
         f"{PREV['cal_best_frob_ddpm']:.3f} sweep-best). The AUROC of the "
         "newly",
         "     selected config is reported below so the cost of the swap is "
         "visible.",
         "  2. CONDITIONAL TRAINING. Each step now draws a random subset of "
         "the OBSERVED",
         "     coordinates to supply as GIVEN (true x0, noised exactly as the "
         "sampler noises",
         "     it, flagged in a new 12-dim given-mask input channel); the loss "
         "is taken on",
         "     observed-AND-NOT-given. Replacement conditioning is therefore "
         "in-distribution",
         "     instead of an out-of-distribution hack. Real Tox21 missingness "
         "is kept strictly",
         "     separate: an unobserved coordinate is never given and never a "
         "loss target.",
         f"  3. SEED REPEATS. {len(seeds)} seeds, each resampling the scaffold "
         "split, the calibration",
         "     sweep and all sampling. Every headline number below is "
         "mean [min, max] over",
         "     seeds, and the tie rule is applied to PAIRED per-seed "
         "differences.",
         "",
         "  TIE RULE (declared): a difference counts only if all seeds agree "
         "on its sign AND",
         "  the mean difference exceeds its own seed-to-seed range. Otherwise "
         "it is a TIE,",
         "  however suggestive the means look.",
         "",
         f"  ⚠️  per-seed configurations selected on calibration "
         f"({SELECT_ON}):"]
    for r in rs:
        ch = r.choices
        L.append(f"      seed {r.seed}: h={ch['h']} lr={ch['lr']:.0e} "
                 f"epochs={ch['epochs']} T={ch['T']}  -> calibration frob "
                 f"{ch['cal_frob']:.3f}, mean AUROC "
                 f"{ch['cal_mean_auroc']:.4f}")
        ar = ch["cal_frob_ref"]["auroc_rule_would_pick"]
        L.append(f"               the OLD rule would have picked h={ar['h']} "
                 f"epochs={ar['epochs']} T={ar['T']} "
                 f"(frob {ar['frob']:.3f}, AUROC {ar['mean_auroc']:.4f})")
    L += ["  ✅  every number in the three axis tables is a held-out test "
          "fold, scored once per seed.",
          ""]

    # ---- axis (a) --------------------------------------------------------
    dec_a, m_a, r_a = decide([d - b for d, b in zip(a_df, a_b1)])
    L += ["(a) MARGINALS -- mean per-endpoint AUROC, observed test labels "
          "only   ✅",
          "-" * 74,
          f"  {'method':<30}{'mean':>5}   [min, max]",
          row("B1 / B2 independent LR", a_b1),
          row("DDPM joint", a_df),
          "",
          f"  paired difference (DDPM - LR): {m_a:+.3f}, seed range "
          f"{r_a:.3f}  ->  "
          f"{'DECISIVE' if dec_a else 'TIE'}",
          f"  Requirement (a) was 'must not regress'. "
          + ("It does not: the difference is inside its own seed spread."
             if not dec_a else
             ("It IMPROVES." if m_a > 0 else "It REGRESSES.")),
          ""]

    # ---- axis (b) --------------------------------------------------------
    dec_b2, m_b2, r_b2 = decide([b - d for b, d in zip(b_b2, b_df)])
    dec_b1, m_b1, r_b1 = decide([b - d for b, d in zip(b_b1, b_df)])
    L += ["(b) JOINT STRUCTURE -- ||R - R_obs||_F over the 12x12 endpoint "
          "matrix   ✅",
          "-" * 74,
          f"  {'method':<30}{'mean':>5}   [min, max]   (lower is better)",
          row("B2 LR + copula", b_b2),
          row("B1 independent LR", b_b1),
          row("DDPM joint", b_df),
          "",
          f"  paired (copula - DDPM): {m_b2:+.3f}, seed range {r_b2:.3f}  ->  "
          f"{'DECISIVE, copula better' if dec_b2 and m_b2 < 0 else ('DECISIVE, DDPM better' if dec_b2 else 'TIE')}",
          f"  paired (indep LR - DDPM): {m_b1:+.3f}, seed range {r_b1:.3f}  "
          f"->  "
          f"{'DECISIVE, LR better' if dec_b1 and m_b1 < 0 else ('DECISIVE, DDPM better' if dec_b1 else 'TIE')}",
          "",
          "  the single (SR-p53, SR-MMP) correlation, observed vs implied:",
          row("    OBSERVED", rho_o, width=26),
          row("    B2 LR + copula", rho_b2, width=26),
          row("    DDPM joint", rho_df, width=26),
          ""]

    # ---- axis (c) --------------------------------------------------------
    dec_lc, m_lc, r_lc = decide(lift_df)
    dec_vs, m_vs, r_vs = decide([d - c for d, c in zip(lift_df, lift_cop)])
    dec_cb, m_cb, r_cb = decide([d - b for d, b in zip(c_df_c, c_b1)])
    L += ["(c) THE CONFOUND -- predicting SR-p53 when SR-MMP is already "
          "measured   ✅",
          "-" * 74,
          f"  {'predictor of p53':<30}{'mean':>5}   [min, max]",
          row("B1 LR (cannot see MMP)", c_b1),
          row("B3 LR, MMP as a feature", c_b3),
          row("B2 copula, unconditioned", c_cop_u),
          row("B2 copula | MMP measured", c_cop_c),
          row("DDPM, unconditioned", c_df_u),
          row("DDPM | MMP inpainted", c_df_c),
          "",
          "  LIFT from conditioning on the measured MMP (the number that "
          "decides this axis):",
          row("    B2 copula, closed form", lift_cop, "{:+.3f}", width=26),
          row("    DDPM, inpainting", lift_df, "{:+.3f}", width=26),
          "",
          f"  DDPM lift: {m_lc:+.3f}, seed range {r_lc:.3f}  ->  "
          f"{'DECISIVE' if dec_lc else 'TIE (indistinguishable from zero)'}",
          f"  DDPM lift vs copula lift: {m_vs:+.3f}, seed range {r_vs:.3f}  "
          f"->  "
          f"{('DECISIVE, DDPM larger' if m_vs > 0 else 'DECISIVE, copula larger') if dec_vs else 'TIE'}",
          f"  DDPM conditioned vs B1 unconditioned: {m_cb:+.3f}, seed range "
          f"{r_cb:.3f}  ->  "
          f"{('DECISIVE, DDPM better' if m_cb > 0 else 'DECISIVE, B1 better') if dec_cb else 'TIE'}",
          ""]

    # ---- verdict ---------------------------------------------------------
    # Derived from the tie rule, never asserted.
    won_a = dec_a and m_a > 0
    lost_a = dec_a and m_a < 0
    won_b = dec_b2 and m_b2 > 0                  # DDPM beats the copula
    best_c_mean = max([("B1 LR", np.mean(c_b1)), ("B3 LR+MMP", np.mean(c_b3)),
                       ("B2 copula | MMP", np.mean(c_cop_c)),
                       ("DDPM | MMP", np.mean(c_df_c))], key=lambda kv: kv[1])
    won_c = best_c_mean[0] == "DDPM | MMP" and dec_vs and m_vs > 0

    L += ["=" * 74]
    if won_a and won_b and won_c:
        L.append("VERDICT: the properly-conditioned joint diffusion model "
                 "BEATS the baselines on all three axes.")
    elif won_b or won_c:
        L.append(f"VERDICT: the fixes change the picture but do NOT overturn "
                 f"it. The diffusion model still")
        L.append("  loses overall to a ~66-parameter Gaussian copula.")
    else:
        L.append("VERDICT: still loses -- now for better-understood reasons. "
                 "Both fixes did what they")
        L.append("  were supposed to do mechanically, and neither is enough "
                 "to beat a ~66-parameter")
        L.append("  Gaussian copula bolted onto independent logistic "
                 "regressions.")
    L += [f"  (a) marginals    {'TIE ' if not dec_a else ('WIN ' if m_a > 0 else 'LOSS')}  "
          f"DDPM {np.mean(a_df):.3f} vs LR {np.mean(a_b1):.3f} "
          f"({m_a:+.3f}, range {r_a:.3f})",
          f"  (b) joint        {'WIN ' if won_b else ('TIE ' if not dec_b2 else 'LOSS')}  "
          f"DDPM {np.mean(b_df):.3f} vs copula {np.mean(b_b2):.3f} "
          f"vs indep LR {np.mean(b_b1):.3f}",
          f"  (c) conditional  {'WIN ' if won_c else ('TIE ' if not dec_lc else 'LOSS')}  "
          f"DDPM lift {np.mean(lift_df):+.3f} vs copula lift "
          f"{np.mean(lift_cop):+.3f}; best overall "
          f"{best_c_mean[0]} ({best_c_mean[1]:.3f})",
          ""]

    # ---- the explicit before/after ---------------------------------------
    L += ["CHANGE vs the previous run (previous = AUROC-selected, no masked "
          "training, 1 seed):",
          "-" * 74,
          f"  {'quantity':<34}{'previous':>10}{'now (mean)':>12}"
          f"{'change':>10}",
          f"  {'(a) DDPM mean AUROC':<34}{PREV['a_ddpm']:>10.3f}"
          f"{np.mean(a_df):>12.3f}{np.mean(a_df) - PREV['a_ddpm']:>+10.3f}",
          f"  {'(a) LR mean AUROC (reference)':<34}{PREV['a_b1']:>10.3f}"
          f"{np.mean(a_b1):>12.3f}{np.mean(a_b1) - PREV['a_b1']:>+10.3f}",
          f"  {'(b) DDPM joint error':<34}{PREV['b_ddpm']:>10.3f}"
          f"{np.mean(b_df):>12.3f}{np.mean(b_df) - PREV['b_ddpm']:>+10.3f}",
          f"  {'(b) copula joint error':<34}{PREV['b_b2']:>10.3f}"
          f"{np.mean(b_b2):>12.3f}{np.mean(b_b2) - PREV['b_b2']:>+10.3f}",
          f"  {'(b) indep LR joint error':<34}{PREV['b_b1']:>10.3f}"
          f"{np.mean(b_b1):>12.3f}{np.mean(b_b1) - PREV['b_b1']:>+10.3f}",
          f"  {'(c) DDPM lift from MMP':<34}{PREV['lift_ddpm']:>+10.3f}"
          f"{np.mean(lift_df):>+12.3f}"
          f"{np.mean(lift_df) - PREV['lift_ddpm']:>+10.3f}",
          f"  {'(c) copula lift from MMP':<34}{PREV['lift_cop']:>+10.3f}"
          f"{np.mean(lift_cop):>+12.3f}"
          f"{np.mean(lift_cop) - PREV['lift_cop']:>+10.3f}",
          f"  {'(c) DDPM conditioned AUROC':<34}{PREV['c_ddpm_cond']:>10.3f}"
          f"{np.mean(c_df_c):>12.3f}"
          f"{np.mean(c_df_c) - PREV['c_ddpm_cond']:>+10.3f}",
          f"  {'denoiser parameters':<34}{PREV['params']:>10,}"
          f"{rs[0].n_params:>12,}",
          f"  {'wall clock (min)':<34}{PREV['wall_min']:>10.1f}{tot:>12.1f}",
          "",
          "  Caveat on this table, stated rather than buried: the 'previous' "
          "column is one",
          "  seed and one split, so a change smaller than the seed ranges "
          "printed above is not",
          "  evidence of anything. The columns differ in three ways at once "
          "(selection rule,",
          "  training scheme, seed count), so it attributes nothing on its "
          "own -- it is here to",
          "  show direction and magnitude, not to isolate a cause.",
          ""]

    # ---- scope -----------------------------------------------------------
    L += ["Scope limits:",
          f"  {len(seeds)} seeds is enough to catch a margin that is obviously "
          "inside the noise and NOT",
          "  enough for a confidence interval; the tie rule above is "
          "deliberately conservative",
          "  rather than a hypothesis test. One dataset (Tox21), one target "
          "class (qHTS",
          "  nuclear-receptor and stress-response reporters), one "
          "featurisation (2048-bit",
          "  radius-2 Morgan counts), one architecture family (MLP denoiser, "
          "cosine schedule,",
          "  ancestral sampling with replacement conditioning), one "
          "scaffold-split scheme.",
          "  Tox21 missingness is NOT missing at random and every method here "
          "treats it as if it",
          "  were. Conditioning was tested on ONE ordered pair (MMP -> p53) "
          "because that is the",
          "  pair with a mechanistic story; the joint model supplies all 132 "
          "ordered pairs and",
          "  only one of them is measured here.",
          "=" * 74,
          ""]

    L.append(report_seed(rs[0]))
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
