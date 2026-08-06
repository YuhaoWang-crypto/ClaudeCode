r"""
From structure to function: what a pattern lets a cell *do*.
============================================================

:mod:`plife.bio.pathways` says which pattern a protein group forms.  This module
asks the next question -- **what capability does that pattern buy?** -- and
answers it with computed quantities rather than narration.  Each function below
is one link in a mechanism chain

    valence + Kd + concentration  ->  model quantity  ->  morphology  ->  capability

and each is checked in ``test_function.py``.

The five capabilities that fall out of percolation
--------------------------------------------------

**1. Ultrasensitivity without cooperativity.**  The gel fraction is exactly zero
below the percolation threshold and rises almost vertically above it, so a
graded input produces a switch-like output.  :func:`hill_coefficient` measures
it: a cross-linking pair of tetravalent proteins with ordinary,
*non-cooperative*, 10 uM bonds gives an effective Hill coefficient of 3-8.  No
allostery, no cooperative binding constants -- the sharpness is a property of
the connectivity transition itself.  This matters because ultrasensitivity is
usually *assumed* to need cooperative binding; percolation supplies it for free.

**2. Concentration buffering.**  Above threshold, added protein goes into the
network instead of into free monomer, so the free concentration is pinned.
:func:`free_monomer` computes it.  A cell whose expression level wanders by 3x
(which is what mRNA-to-protein scatter implies) still sees a nearly constant
*free* concentration -- the quantity that actually drives downstream binding.

**3. Co-concentration catalysis.**  Two reactants each enriched K-fold in a
condensate meet at a rate raised by :math:`K_E K_S` (a bimolecular rate goes as
the product).  :func:`rate_enhancement` does the bookkeeping including the
volume penalty, which is what stops the naive "10^4-fold!" claim from being
right.

**4. Exclusion-driven bistability.**  This is the interesting one.  If
phosphorylation creates the valence, and the resulting network sterically
excludes the phosphatase, then phosphorylation feeds back on itself and the
system becomes **bistable**: the same kinase/phosphatase ratio supports both an
off and an on state, with hysteresis.  :func:`phospho_bistability` finds the
window.  Nothing in the model was built to produce bistability -- it emerges
from percolation plus size exclusion, both of which were fixed by independent
measurements.

**5. Polymer-length sensing.**  When one partner's valence *is* its length
(dsDNA for cGAS, mRNA for stress granules, ubiquitin chains for p62), the
percolation threshold becomes a **length** threshold.
:func:`min_polymer_length` computes it.  Self/non-self discrimination by ligand
length is then not a special recognition mechanism but the valence gate.

What is a proxy and what is not
-------------------------------
``gel_fraction`` and ``free_monomer`` are exact for the tree (Bethe) topology
that Wertheim TPT1 already assumes, so they are as good as the free energy
itself.  ``rate_enhancement`` is exact bookkeeping given the partition
coefficients.  ``phospho_bistability`` is a *minimal* two-parameter feedback
model built on top of the computed gel fraction -- the bistability is a genuine
consequence of the computed curve, but the rate constants are illustrative and
the model has no spatial or kinetic detail.  Labelled accordingly at each call.
"""

from __future__ import annotations

import numpy as np

from .sites import SiteMixture, percolates, unbonded_fractions
from .units import conc_to_density, density_to_conc

__all__ = ["gel_fraction", "free_monomer", "dose_response", "hill_coefficient",
           "rate_enhancement", "phospho_bistability", "min_polymer_length",
           "buffering_quality", "competitive_occupancy", "conditional_valence"]


# --------------------------------------------------------------------------- #
#  the branching process, solved for its extinction probability
# --------------------------------------------------------------------------- #
def _bond_stats(conc_uM, sm: SiteMixture):
    rho = conc_to_density(np.atleast_2d(np.asarray(conc_uM, float)))
    X = unbonded_fractions(rho, sm)[0]
    D = sm.delta(rho)[0]
    pi = X[:, None] * D * (rho[0][sm.owner] * sm.n * X)[None, :]
    return rho, X, pi


def _sol_prob(conc_uM, sm: SiteMixture, itmax=4000, tol=1e-13):
    r"""Per-site probability :math:`t_k` that site *k* leads to a **finite** cluster.

    Fixed point of the multitype branching process:

    .. math::  t_k = X_k + \sum_j \pi_{kj}\prod_m t_m^{\,n_{s(j)m}-\delta_{jm}}

    -- a site is either unbonded (:math:`X_k`), or bonded to some *j*-site whose
    molecule must then terminate through all of *its* remaining sites.

    Below the percolation threshold :math:`t\equiv1` is the *only* root, and it
    is returned exactly rather than iterated towards.  That matters: the
    iteration converges at a rate set by the branching eigenvalue, so just short
    of the transition it stalls and leaves a residue of order 1e-6, which a
    downstream "is there a gel?" test would read as a real (tiny) network.  The
    spectral radius decides it exactly, so use it.
    """
    _, X, pi = _bond_stats(conc_uM, sm)
    A = sm.A
    same = (sm.owner[:, None] == sm.owner[None, :]).astype(float)
    mult = np.maximum(same * sm.n[None, :] - np.eye(A), 0.0)   # n_{s(j)m} - delta_jm
    if np.max(np.abs(np.linalg.eigvals(mult @ pi))) < 1.0:
        return np.ones(A)
    t = np.full(A, 0.5)
    for _ in range(itmax):
        with np.errstate(divide="ignore"):
            prod = np.exp(mult @ np.log(np.maximum(t, 1e-300)))
        tn = X + pi @ prod
        if np.max(np.abs(tn - t)) < tol:
            return np.clip(tn, 0.0, 1.0)
        t = 0.5 * (t + tn)
    return np.clip(t, 0.0, 1.0)


def gel_fraction(conc_uM, sm: SiteMixture):
    r"""Fraction of each species that belongs to the **infinite** network.

    A molecule is in the finite sol only if *every* one of its sites terminates,
    so :math:`g_a = 1-\prod_k t_k^{\,n_{ak}}`.  Exactly zero below the
    percolation threshold and rising almost vertically above it -- that shape is
    the origin of capability 1.
    """
    t = _sol_prob(conc_uM, sm)
    with np.errstate(divide="ignore"):
        logt = np.log(np.maximum(t, 1e-300))
    sol = np.array([np.exp(float((sm.n * logt)[sm.owner == a].sum()))
                    for a in range(sm.S)])
    return 1.0 - sol


def free_monomer(conc_uM, sm: SiteMixture):
    r"""Concentration (uM) of each species with **all** sites unbonded.

    :math:`c_a^{\rm free} = c_a\prod_k X_k^{\,n_{ak}}`.  This is the species that
    a downstream partner or a biosensor actually titrates, and above the gel
    point it stops tracking the total -- capability 2.
    """
    c = np.atleast_1d(np.asarray(conc_uM, float))
    _, X, _ = _bond_stats(c, sm)
    with np.errstate(divide="ignore"):
        logX = np.log(np.maximum(X, 1e-300))
    frac = np.array([np.exp(float((sm.n * logX)[sm.owner == a].sum()))
                     for a in range(sm.S)])
    return c * frac


# --------------------------------------------------------------------------- #
#  dose-response and its sharpness
# --------------------------------------------------------------------------- #
def dose_response(sm: SiteMixture, base_conc_uM, scale, observable="gel"):
    """Sweep a whole composition by ``scale`` and report gel fraction / branching.

    Returns ``dict(scale, conc, gel, branching, free)``; ``gel`` and ``free`` are
    ``(len(scale), S)``.
    """
    base = np.asarray(base_conc_uM, float)
    scale = np.atleast_1d(np.asarray(scale, float))
    gel = np.empty((len(scale), sm.S))
    free = np.empty((len(scale), sm.S))
    lam = np.empty(len(scale))
    for i, s in enumerate(scale):
        c = base * s
        gel[i] = gel_fraction(c, sm)
        free[i] = free_monomer(c, sm)
        lam[i] = percolates(c, sm)[1]
    return dict(scale=scale, conc=base[None, :] * scale[:, None],
                gel=gel, branching=lam, free=free)


def hill_coefficient(x, y, at=0.5):
    r"""Effective Hill coefficient of a dose-response curve at fractional output ``at``.

    .. math::  n_H = \frac{{\rm d}\ln[y/(1-y)]}{{\rm d}\ln x}

    evaluated where ``y`` first crosses ``at``.  This is the standard way to
    quote steepness, and it is what makes the percolation result comparable to
    the cooperative-binding literature: ``n_H = 1`` is a plain hyperbolic
    binding curve, haemoglobin is ~2.8, and the numbers this model produces from
    strictly non-cooperative bonds are larger than both.
    """
    x = np.asarray(x, float)
    y = np.clip(np.asarray(y, float), 1e-12, 1 - 1e-12)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[ok], y[ok]
    if len(x) < 3 or y.min() >= at or y.max() <= at:
        return np.nan
    i = int(np.argmax(y >= at))
    lo, hi = max(i - 1, 0), min(i + 1, len(x) - 1)
    if hi <= lo:
        return np.nan
    logit = np.log(y / (1 - y))
    return float((logit[hi] - logit[lo]) / (np.log(x[hi]) - np.log(x[lo])))


# --------------------------------------------------------------------------- #
#  capability 3: co-concentration
# --------------------------------------------------------------------------- #
def rate_enhancement(K_enzyme, K_substrate, phi_condensate=0.01,
                     viscosity_penalty=1.0):
    r"""Fold change in a bimolecular reaction flux caused by co-concentration.

    A bimolecular rate is :math:`k[E][S]`, so enriching both by :math:`K_E` and
    :math:`K_S` raises the *local* rate by :math:`K_EK_S` -- but only a fraction
    :math:`\phi` of the cell volume is condensate, so the **whole-cell** flux
    gain is

    .. math::  \frac{J_{\rm cond}}{J_{\rm bulk}}
               = \frac{\phi K_EK_S + (1-\phi)}{1}\Big/\eta

    with :math:`\eta` an optional slow-down for the crowded interior (diffusion
    in a dense condensate is typically 1-2 orders slower, which is why the naive
    :math:`K_EK_S` figure is an over-estimate).  Returns
    ``(whole_cell_fold, local_fold)``.
    """
    KE = np.asarray(K_enzyme, float)
    KS = np.asarray(K_substrate, float)
    phi = np.asarray(phi_condensate, float)
    local = KE * KS / max(viscosity_penalty, 1e-30)
    whole = (phi * local + (1.0 - phi))
    return whole, local


# --------------------------------------------------------------------------- #
#  capability 4: exclusion-driven bistability
# --------------------------------------------------------------------------- #
def phospho_bistability(gel_of_p, ratio, access_in_gel=0.05, n_p=201):
    r"""Steady states of a scaffold whose phosphorylation builds its own network.

    Minimal feedback loop, and every ingredient except the two rate constants is
    computed elsewhere in this package:

    * phosphorylation fraction :math:`p` sets the scaffold's valence, so it sets
      the gel fraction :math:`g(p)` -- supplied as the callable ``gel_of_p``;
    * the network **sterically excludes the phosphatase**, whose access is
      therefore :math:`\alpha(p)=1-g(p)\,(1-a)` with ``a = access_in_gel`` taken
      from the size-exclusion partition coefficient;
    * the kinase is not excluded (it is recruited, or simply smaller).

    Steady state of :math:`\dot p = k_{\rm kin}(1-p)-k_{\rm ppt}\,p\,\alpha(p)`:

    .. math::  R \;\equiv\; \frac{k_{\rm kin}}{k_{\rm ppt}}
               \;=\; \frac{p\,\alpha(p)}{1-p}

    The right-hand side is monotonic without feedback, so there is exactly one
    root.  With feedback :math:`\alpha` falls as :math:`p` rises, and the curve
    can turn over -- giving **three** roots, i.e. bistability with hysteresis.

    Returns ``dict(ratio, roots, n_states, stable, on, off, bistable,
    p_grid, rhs)``.  ``on``/``off`` are the outer (stable) branches.

    .. note::
       ⚠️ The rate constants are illustrative; the *bistability* is a genuine
       consequence of the computed ``gel_of_p`` curve and the measured exclusion
       factor, not of the numbers chosen for ``ratio``.
    """
    p = np.linspace(1e-4, 1 - 1e-4, n_p)
    g = np.array([float(gel_of_p(pi)) for pi in p])
    alpha = 1.0 - g * (1.0 - access_in_gel)
    rhs = p * alpha / (1.0 - p)                      # = R at steady state

    ratio = np.atleast_1d(np.asarray(ratio, float))
    out = []
    for R in ratio:
        d = rhs - R
        sign = np.sign(d)
        idx = np.where(np.diff(sign) != 0)[0]
        roots = []
        for i in idx:
            x0, x1 = p[i], p[i + 1]
            y0, y1 = d[i], d[i + 1]
            roots.append(float(x0 - y0 * (x1 - x0) / (y1 - y0)) if y1 != y0 else float(x0))
        roots = sorted(roots)
        # stability: dp/dt = k_ppt[(R)(1-p) - p alpha(p)], so a root is stable
        # where d/dp of that expression is negative, i.e. where rhs is INCREASING
        stable = []
        for r in roots:
            j = int(np.clip(np.searchsorted(p, r), 1, n_p - 1))
            stable.append(bool(rhs[j] > rhs[j - 1]))
        out.append(dict(R=float(R), roots=roots, stable=stable,
                        n_states=len(roots),
                        bistable=sum(stable) >= 2,
                        off=roots[0] if roots else np.nan,
                        on=roots[-1] if roots else np.nan))
    return dict(ratio=ratio, p_grid=p, gel=g, access=alpha, rhs=rhs, states=out,
                bistable_any=any(o["bistable"] for o in out),
                window=[float(rhs[np.argmax(rhs)]), float(np.min(rhs[np.argmax(rhs):]))]
                if np.argmax(rhs) < n_p - 1 else None)


# --------------------------------------------------------------------------- #
#  capability 5: polymer length as valence
# --------------------------------------------------------------------------- #
def min_polymer_length(build, lengths, conc_uM):
    r"""Shortest polymer that switches the network on.

    ``build(L)`` must return a :class:`SiteMixture` in which one partner's
    valence is set by the polymer length ``L``.  Returns
    ``dict(lengths, branching, gel, L_min)`` with ``L_min`` the first length
    whose branching eigenvalue reaches 1 (``inf`` if none does).

    This is what turns "how does cGAS ignore short self-DNA?" from a puzzle
    about recognition into arithmetic about valence.
    """
    lengths = np.atleast_1d(np.asarray(lengths, float))
    lam = np.empty(len(lengths))
    gel = np.empty(len(lengths))
    for i, L in enumerate(lengths):
        sm = build(float(L))
        lam[i] = percolates(conc_uM, sm)[1]
        gel[i] = gel_fraction(conc_uM, sm).max()
    on = lam >= 1.0
    return dict(lengths=lengths, branching=lam, gel=gel,
                L_min=float(lengths[on][0]) if on.any() else np.inf)


# --------------------------------------------------------------------------- #
#  allostery: a ligand that changes a *partner's* valence
# --------------------------------------------------------------------------- #
def competitive_occupancy(site_total_uM, ligand_totals_uM, kd_uM):
    r"""Exact fractional occupancy of one site class by several competing ligands.

    Solves site and ligand conservation simultaneously rather than assuming
    ligand excess:

    .. math::  S_T = s\Bigl(1+\sum_i \frac{L_{T,i}}{K_{d,i}+s}\Bigr),
               \qquad \theta_i = \frac{s\,L_{T,i}}{(K_{d,i}+s)\,S_T}

    (monotone in the free-site concentration :math:`s`, so a bisection is
    unconditionally safe).  Returns ``(theta, free_site_fraction)``.

    This exists because :class:`~plife.bio.sites.SiteMixture` cannot express
    **allostery** -- see :func:`conditional_valence`.
    """
    S_T = float(site_total_uM)
    L = np.atleast_1d(np.asarray(ligand_totals_uM, float))
    kd = np.atleast_1d(np.asarray(kd_uM, float))
    if S_T <= 0:
        return np.zeros_like(L), 1.0
    lo, hi = 0.0, S_T
    for _ in range(200):
        s = 0.5 * (lo + hi)
        f = s * (1.0 + np.sum(L / (kd + s))) - S_T
        lo, hi = (lo, s) if f > 0 else (s, hi)
    s = 0.5 * (lo + hi)
    theta = s * L / ((kd + s) * S_T)
    return theta, float(s / S_T)


def conditional_valence(n_closed, n_open, theta):
    r"""Mean-field allostery: valence interpolated by regulatory-site occupancy.

    .. math::  n_{\rm eff} = n_{\rm closed} + (n_{\rm open}-n_{\rm closed})\,\theta

    **Why this is needed at all.**  Wertheim TPT1 -- and therefore everything in
    :mod:`plife.bio.thermo` and :mod:`plife.bio.sites` -- treats every site as
    having a fixed multiplicity and a fixed Kd.  A protein whose binding valence
    *changes* when a partner docks somewhere else is doing allostery, and a flat
    model cannot represent it.  Get this wrong and the model does not merely
    lose accuracy, it predicts the opposite sign: with G3BP1's RGG valence held
    fixed, adding CAPRIN1 only competes for scarce RNA and *dissolves* the
    granule, whereas CAPRIN1 in fact promotes it by relieving G3BP1's IDR1
    autoinhibition.

    Usage is deliberately two-stage and explicit -- occupancy first, then
    network -- so the allosteric assumption stays visible instead of being
    buried in an effective Kd::

        theta, _ = competitive_occupancy(site_total, [c_activator, c_capper],
                                         [kd_act, kd_cap])
        n_eff = conditional_valence(n_closed, n_open, theta[0])
        sm = build_mixture(n_eff);  percolates(conc, sm)

    ⚠️ This is mean-field: it assumes the regulatory equilibrium is fast
    compared with network assembly and that valence responds linearly to
    occupancy.  Both are assumptions, not results.
    """
    return float(n_closed) + (float(n_open) - float(n_closed)) * np.asarray(theta, float)


# --------------------------------------------------------------------------- #
def buffering_quality(total_uM, free_uM):
    r"""Local slope :math:`{\rm d}\ln c_{\rm free}/{\rm d}\ln c_{\rm total}`.

    1 means no buffering (free tracks total); 0 means perfect buffering (free
    pinned however much total is added).  Report the value *above* threshold --
    that is where the capability lives.
    """
    t = np.asarray(total_uM, float)
    f = np.asarray(free_uM, float)
    ok = (t > 0) & (f > 0)
    return np.gradient(np.log(f[ok]), np.log(t[ok]))
