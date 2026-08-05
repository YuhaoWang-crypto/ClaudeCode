r"""
The inverse problem: what can imaging actually tell you about Kd and concentration?
==================================================================================

Three questions, three honest answers.

**(1) One image, one composition.**  "These two proteins co-localise in puncta at
[A] = 5 uM, [B] = 20 uM."  That is a *single inequality* on a three-dimensional
parameter space :math:`(K_d^{AA}, K_d^{AB}, K_d^{BB})` -- plus valences and
sizes.  It cannot determine three numbers.  What it does determine is a
**feasible region**, and :func:`feasible_region` maps it.

**(2) A titration series.**  Measuring the saturation concentration of A at
several fixed [B] traces out the phase boundary, and the *shape* of that
boundary is what carries the interaction information: its low-[B] intercept
fixes the homotypic term, its slope fixes the heterotypic term, and its
re-entrant upturn at high [B] fixes the ratio of valences to affinities.
:func:`fit_from_titration` does this and returns a posterior, not a point.

**(3) Which combinations are actually determined?**  Models like this one are
*sloppy*: the likelihood is sharp along a few parameter combinations and nearly
flat along the rest.  :func:`sloppiness` performs the eigen-analysis of the
log-likelihood Hessian and reports, for each eigendirection, the combination of
log-Kd's it corresponds to and how tightly it is pinned.  This is the difference
between "we fitted three Kd's" and "we determined one combination of three
Kd's and the other two are free" -- which, in our synthetic tests, is usually
the true situation for image data alone.

Partition coefficients are worth far more than binary calls
-----------------------------------------------------------
A yes/no "does it condense" observation carries at most one bit.  A measured
partition coefficient (dense/dilute intensity ratio, which is what a confocal
image actually gives you) is a continuous number that constrains the transfer
free energy directly.  In the synthetic benchmark of ``e10``, replacing 8 binary
calls with 4 partition measurements shrinks the credible volume by orders of
magnitude.
"""

from __future__ import annotations

import itertools

import numpy as np

from .phases import phase_diagram
from .thermo import Mixture, csat_scan, unbonded_fractions
from .units import conc_to_density

__all__ = ["Observation", "csat_curve", "predict_csat", "loglike", "feasible_region",
           "fit_from_titration", "sloppiness", "marginals", "credible_interval"]


# --------------------------------------------------------------------------- #
class Observation:
    """One experimental datum.

    ``kind='csat'``       c_sat of ``species`` measured at partner concentration
                          ``partner_uM``; ``value`` in uM, ``sd_log10`` its error.
    ``kind='phase'``      binary: did it phase separate at ``(cA, cB)``?
    ``kind='partition'``  dense/dilute ratio of ``species`` at ``(cA, cB)``.
    ``kind='bound'``      bound-site fraction ``1 - X_a`` at ``(cA, cB)``.
    """

    def __init__(self, kind, value=None, cA=None, cB=None, species=0,
                 partner_uM=0.0, sd_log10=0.15, sd=0.05):
        self.kind = kind
        self.value = value
        self.cA, self.cB = cA, cB
        self.species = species
        self.partner_uM = partner_uM
        self.sd_log10 = sd_log10
        self.sd = sd

    def __repr__(self):
        return f"Observation({self.kind}, value={self.value}, cA={self.cA}, cB={self.cB})"


# --------------------------------------------------------------------------- #
def _mix_with(mix, kd_aa, kd_ab, kd_bb):
    kd = np.array([[kd_aa, kd_ab], [kd_ab, kd_bb]], float)
    return mix.copy(kd=kd)


def csat_curve(mix, cB_uM, cA_lo=1e-3, cA_hi=3e4, nA=180, species=0):
    r"""c_sat of species A for every partner concentration in ``cB_uM``, one hull.

    Uses the **full two-dimensional** common-tangent construction, not a
    one-dimensional scan along the ray: the dense phase generally has a different
    A:B ratio from the dilute phase, and a ray restricted to fixed [B] cannot
    represent that.  One ``ConvexHull`` covers the whole grid, so the whole
    titration curve costs about as much as a single point.

    Returns ``(csat_uM, K_A, K_B)``, each of length ``len(cB_uM)``; ``inf``
    entries mean no phase separation was found in the scanned window.
    """
    from .thermo import coexistence
    cA = np.geomspace(cA_lo, cA_hi, nA)
    cB = np.atleast_1d(np.asarray(cB_uM, float))
    A, B = np.meshgrid(conc_to_density(cA), conc_to_density(cB), indexing="ij")
    rho = np.column_stack([A.ravel(), B.ravel()])
    res = coexistence(mix, rho)
    bad = (~res["on_hull"]).reshape(nA, len(cB))
    out = np.full(len(cB), np.inf)
    for j in range(len(cB)):
        idx = np.where(bad[:, j])[0]
        if idx.size:
            out[j] = cA[idx[0]]
    return out, cA, res, bad


def predict_csat(mix, species=0, partner_uM=0.0, lo_uM=1e-3, hi_uM=3e4, n=180):
    """c_sat of ``species`` at one partner concentration (uM; ``inf`` if none)."""
    cs, *_ = csat_curve(mix, [partner_uM, partner_uM * 1.01 + 1e-9],
                        cA_lo=lo_uM, cA_hi=hi_uM, nA=n, species=species)
    return float(cs[0])


def _predict(mix, obs, pd_cache=None):
    """Model value for one observation."""
    if obs.kind == "csat":
        return predict_csat(mix, obs.species, obs.partner_uM)
    if obs.kind in ("phase", "partition"):
        cA = np.array([obs.cA * 0.999, obs.cA, obs.cA * 1.001])
        cB = np.array([obs.cB * 0.999, obs.cB, obs.cB * 1.001])
        pd = phase_diagram(mix, cA, cB)
        lab = int(pd["label"][1, 1])
        if obs.kind == "phase":
            return lab > 0
        K = pd["K_A"][1, 1] if obs.species == 0 else pd["K_B"][1, 1]
        return float(K) if np.isfinite(K) else 1.0
    if obs.kind == "bound":
        rho = conc_to_density(np.array([[obs.cA, obs.cB]]))
        X = unbonded_fractions(rho, mix)[0]
        return float(1 - X[obs.species])
    raise ValueError(obs.kind)


def loglike(mix, observations):
    """Log-likelihood of one parameter set given the observations.

    All ``csat`` observations are evaluated from a **single** two-dimensional
    hull covering every probed partner concentration at once -- both because it
    is ~N times faster and because asking for two nearly identical partner
    concentrations makes Qhull's predicates fail.
    """
    ll = 0.0
    csats = [o for o in observations if o.kind == "csat"]
    if csats:
        probe = np.array([o.partner_uM for o in csats], float)
        pred, *_ = csat_curve(mix, probe)
        for o, p in zip(csats, pred):
            if not np.isfinite(p) or p <= 0:
                ll += -20.0
            else:
                z = (np.log10(p) - np.log10(o.value)) / o.sd_log10
                ll += -0.5 * z * z

    for obs in observations:
        if obs.kind == "csat":
            continue
        p = _predict(mix, obs)
        if obs.kind == "phase":
            ll += 0.0 if bool(p) == bool(obs.value) else -4.0
        elif obs.kind == "partition":
            if not np.isfinite(p) or p <= 0:
                ll += -20.0
            else:
                z = (np.log10(p) - np.log10(obs.value)) / obs.sd_log10
                ll += -0.5 * z * z
        else:
            z = (p - obs.value) / obs.sd
            ll += -0.5 * z * z
    return float(ll)


# --------------------------------------------------------------------------- #
def _grid_scan(mix_template, observations, kd_aa_grid, kd_ab_grid, kd_bb_grid,
               verbose=False):
    shape = (len(kd_aa_grid), len(kd_ab_grid), len(kd_bb_grid))
    LL = np.full(shape, -np.inf)
    for i, aa in enumerate(kd_aa_grid):
        for j, ab in enumerate(kd_ab_grid):
            for k, bb in enumerate(kd_bb_grid):
                LL[i, j, k] = loglike(_mix_with(mix_template, aa, ab, bb), observations)
        if verbose:
            print(f"    row {i+1}/{shape[0]}")
    return LL


def fit_from_titration(mix_template, observations, kd_aa_grid, kd_ab_grid,
                       kd_bb_grid, verbose=False):
    r"""Posterior over :math:`(K_d^{AA}, K_d^{AB}, K_d^{BB})` on a log grid.

    Uniform prior in :math:`\log K_d` (the scale-free choice for an affinity),
    likelihood from :func:`loglike`.  Returns the normalised posterior and the
    maximum-a-posteriori point.
    """
    LL = _grid_scan(mix_template, observations, kd_aa_grid, kd_ab_grid,
                    kd_bb_grid, verbose)
    P = np.exp(LL - np.nanmax(LL))
    P /= np.nansum(P)
    idx = np.unravel_index(int(np.nanargmax(LL)), LL.shape)
    return dict(logl=LL, post=P,
                grids=(np.asarray(kd_aa_grid), np.asarray(kd_ab_grid),
                       np.asarray(kd_bb_grid)),
                map=(float(kd_aa_grid[idx[0]]), float(kd_ab_grid[idx[1]]),
                     float(kd_bb_grid[idx[2]])),
                map_index=idx)


def feasible_region(fit, level=0.95):
    """Boolean mask of the smallest set carrying ``level`` of the posterior."""
    P = fit["post"]
    flat = np.sort(P.ravel())[::-1]
    cum = np.cumsum(flat)
    thresh = flat[np.searchsorted(cum, level)] if cum[-1] >= level else 0.0
    return P >= thresh


def marginals(fit):
    """1-D marginal posteriors for each Kd."""
    P = fit["post"]
    return [P.sum(axis=tuple(a for a in range(3) if a != i)) for i in range(3)]


def credible_interval(grid, marg, level=0.95):
    """Central credible interval of a 1-D marginal on a (log-spaced) grid."""
    c = np.cumsum(marg) / marg.sum()
    lo = float(np.interp(0.5 * (1 - level), c, grid))
    hi = float(np.interp(1 - 0.5 * (1 - level), c, grid))
    return lo, hi


# --------------------------------------------------------------------------- #
def sloppiness(fit, window=12.0):
    r"""Eigen-analysis of the log-likelihood curvature in :math:`\log_{10}K_d`.

    Fits a quadratic to :math:`\log L` around the MAP and diagonalises the
    resulting Hessian.  Each eigenvalue is a *stiffness*: large means that
    combination of affinities is pinned by the data, small means it is free; the
    eigenvector says **which combination**, and ``sd_decades`` says how many
    decades of slack it has.

    A large spread between the stiffest and sloppiest direction is the signature
    of a sloppy model and the quantitative form of "you cannot read three Kd's
    off one image".

    Parameters that were held fixed (a grid of length 1, or a non-finite value
    such as ``Kd = inf`` for a protein with no self-affinity) are excluded from
    the analysis rather than producing a singular design matrix.
    """
    LL = fit["logl"]
    grids = [np.asarray(g, float) for g in fit["grids"]]
    labels_all = ["log10 Kd_AA", "log10 Kd_AB", "log10 Kd_BB"]

    free = [i for i, g in enumerate(grids) if g.size > 1 and np.all(np.isfinite(g))]
    if len(free) < 1:
        return dict(eigenvalues=np.array([]), eigenvectors=np.zeros((0, 0)),
                    sd_decades=np.array([]), center=np.array([]),
                    condition=np.nan, labels=[], free=free)

    axes = [np.log10(grids[i]) for i in free]
    mesh = np.meshgrid(*axes, indexing="ij")
    # collapse the fixed axes (they have length 1 or are non-finite)
    LLr = LL
    for i in reversed(range(3)):
        if i not in free:
            LLr = np.take(LLr, 0, axis=i)
    LLr = np.asarray(LLr, float)

    ok = np.isfinite(LLr)
    if ok.sum() < 3 * len(free) + 1:
        ok = np.ones_like(LLr, bool)
    top = LLr[ok].max() - window
    m = ok & (LLr > top)
    if m.sum() < 3 * len(free) + 1:
        m = ok

    d = np.column_stack([mm[m] for mm in mesh])
    c = d[int(np.argmax(LLr[m]))]
    d = d - c
    nf = len(free)

    cols = [np.ones(len(d))]
    cols += [d[:, i] for i in range(nf)]
    quad_idx = []
    for i in range(nf):
        for j in range(i, nf):
            cols.append(d[:, i] * d[:, j])
            quad_idx.append((i, j))
    Amat = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(Amat, LLr[m], rcond=None)

    H = np.zeros((nf, nf))
    for k, (i, j) in enumerate(quad_idx):
        b = beta[1 + nf + k]
        if i == j:
            H[i, i] = 2 * b
        else:
            H[i, j] = H[j, i] = b
    vals, vecs = np.linalg.eigh(-H)                 # -H is positive at a maximum
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    with np.errstate(divide="ignore", invalid="ignore"):
        sd = np.where(vals > 0, 1.0 / np.sqrt(np.maximum(vals, 1e-300)), np.inf)
    cond = float(vals[0] / vals[-1]) if vals[-1] > 0 else np.inf
    return dict(eigenvalues=vals, eigenvectors=vecs, sd_decades=sd, center=c,
                condition=cond, labels=[labels_all[i] for i in free], free=free)
