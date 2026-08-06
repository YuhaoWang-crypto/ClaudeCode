r"""
Physical units and the Kd -> model-parameter bridge.
====================================================

Everything downstream works in **nm, kT, and number-per-nm^3**.  This module is
the only place where laboratory quantities (Kd in M, concentration in uM,
molecular weight in kDa, mRNA in TPM) enter, so every assumption is visible in
one file.

The three conversions that matter
---------------------------------

**Concentration.**  ``1 M = 0.60221 nm^-3`` exactly (Avogadro / litre in nm^3).
A 10 uM protein is 6.0e-6 nm^-3 -- at a 2.5 nm radius that is a volume fraction
of 4e-4, i.e. *very* dilute.  This is why an ideal-gas reference works for the
dilute phase and why crowding corrections only matter inside condensates.

**Size.**  Compact globular proteins have a partial specific volume of
0.73 cm^3/g, so ``V = 1.212 A^3/Da`` and ``R = (3V/4pi)^(1/3)``.  A 50 kDa
protein comes out at 2.44 nm, which is right.  Disordered proteins are much
larger for the same mass -- use ``radius_idp`` or pass a measured R_h.

**Affinity.**  This is the one that people get wrong.  A dissociation constant
is an *equilibrium constant for forming one bond*, so in Wertheim's theory it is
already exactly the association volume:

.. math::  \Delta_{ab} \;=\; \frac{g^{HS}_{ab}(\sigma_{ab};\eta)}{K_{d,ab}^{(n)}}

with :math:`K_d^{(n)}` the Kd converted to nm^-3.  **No well-depth or
well-width convention is needed** -- that is the main reason ``thermo.py`` uses
Wertheim rather than an isotropic square well.  ``kd_to_epsilon`` below is
provided only for the isotropic/RPA route, and it *does* require a range
convention, which is why it takes ``delta`` explicitly and is labelled.

What is assumption, not measurement
-----------------------------------
* ``protein_from_mrna`` is a *scaling with a stated spread*, not a prediction.
  Across genes, mRNA explains roughly 40-60 % of protein-level variance; the
  function therefore returns a log-normal band and every downstream answer
  should be read as a probability, not a yes/no.
* Cytoplasmic viscosity is scale- and probe-dependent.  ``diffusion`` takes it
  as an argument and the default (5x water) is a placeholder, not a constant of
  nature.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "M_TO_NM3", "KT_KCAL", "radius_from_mw", "radius_idp", "conc_to_density",
    "density_to_conc", "kd_to_delta", "kd_to_epsilon", "epsilon_to_kd",
    "b2_from_kd", "volume_fraction", "diffusion", "protein_from_mrna",
    "copies_to_conc",
]

#: 1 mol/L expressed in molecules per nm^3
M_TO_NM3 = 6.02214076e23 / 1e24          # = 0.602214076
#: kT at 310 K in kcal/mol (body temperature)
KT_KCAL = 1.98720425e-3 * 310.0          # = 0.6160


# --------------------------------------------------------------------------- #
#  size
# --------------------------------------------------------------------------- #
def radius_from_mw(mw_kda, specific_volume_A3_per_Da=1.212):
    """Radius (nm) of a compact globular protein of mass ``mw_kda``.

    ``V = v_bar * MW`` with ``v_bar = 1.212 A^3/Da`` (0.73 cm^3/g), then
    ``R = (3V/4pi)^(1/3)``.  50 kDa -> 2.44 nm.
    """
    mw = np.asarray(mw_kda, float) * 1e3
    v_nm3 = mw * specific_volume_A3_per_Da * 1e-3      # A^3 -> nm^3
    return (3.0 * v_nm3 / (4.0 * np.pi)) ** (1.0 / 3.0)


def radius_idp(n_residues, nu=0.522, b_nm=0.254):
    r"""Radius of gyration (nm) of a disordered chain, ``R_g = b N^nu``.

    Defaults are the Marsh & Forman-Kay scaling fitted to SAXS/NMR on unfolded
    and intrinsically disordered proteins, :math:`R_g = 2.54\,N^{0.522}` A.  A
    300-residue IDP comes out at 5.0 nm against 4.1 nm for a compact globule of
    the same mass -- which is why disordered proteins reach the overlap
    concentration so much sooner, and why they dominate condensate scaffolds.

    ``nu`` is the solvent-quality exponent: 0.5 at theta, ~0.6 in good solvent.
    Note this returns :math:`R_g`, not :math:`R_h`; for a Gaussian coil
    :math:`R_h\approx0.64R_g`, so treating :math:`R_g` as the excluded-volume
    radius deliberately errs towards *over*-estimating the steric size.
    """
    return b_nm * np.asarray(n_residues, float) ** nu


def volume_fraction(density_nm3, radius_nm):
    """Volume fraction ``phi = rho * (4/3) pi R^3``."""
    r = np.asarray(radius_nm, float)
    return np.asarray(density_nm3, float) * (4.0 / 3.0) * np.pi * r**3


# --------------------------------------------------------------------------- #
#  concentration
# --------------------------------------------------------------------------- #
def conc_to_density(conc_uM):
    """uM -> molecules per nm^3."""
    return np.asarray(conc_uM, float) * 1e-6 * M_TO_NM3


def density_to_conc(density_nm3):
    """Molecules per nm^3 -> uM."""
    return np.asarray(density_nm3, float) / (1e-6 * M_TO_NM3)


def copies_to_conc(copies, cell_volume_um3=2000.0):
    """Copy number -> uM, given a cell (or compartment) volume in um^3.

    2000 um^3 is a typical mammalian cell; a HeLa nucleus is ~500 um^3, a yeast
    cell ~40 um^3.  1e6 copies in 2000 um^3 is 0.83 uM.
    """
    v_nm3 = np.asarray(cell_volume_um3, float) * 1e9
    return density_to_conc(np.asarray(copies, float) / v_nm3)


# --------------------------------------------------------------------------- #
#  affinity
# --------------------------------------------------------------------------- #
def kd_to_delta(kd_M, g_contact=1.0):
    r"""Wertheim association volume :math:`\Delta` (nm^3) from a measured Kd.

    .. math::  \Delta = g^{HS}(\sigma)\,/\,K_d^{(n)},\qquad
               K_d^{(n)} = K_d[\mathrm{M}] \times 0.6022\ \mathrm{nm^{-3}}

    ``g_contact`` is the hard-sphere contact value of the pair distribution
    function, which enhances association inside a crowded (dense) phase; pass 1
    for the dilute limit and ``thermo.g_contact_cs(eta)`` inside a condensate.

    This is an identity, not a model: Wertheim's :math:`\Delta` *is* the
    equilibrium constant for forming one bond.
    """
    kd_n = np.asarray(kd_M, float) * M_TO_NM3
    with np.errstate(divide="ignore"):
        return np.where(kd_n > 0, g_contact / np.maximum(kd_n, 1e-300), np.inf)


def kd_to_epsilon(kd_M, sigma_nm, delta_nm=0.5):
    r"""**Isotropic-model only.**  Square-well depth (in kT) reproducing ``kd``.

    .. math::  \frac1{K_d^{(n)}} = v_{\rm well}\left(e^{\beta\epsilon}-1\right),
               \qquad v_{\rm well}=\tfrac{4\pi}{3}\bigl[(\sigma+\delta)^3-\sigma^3\bigr]

    .. warning::
       The answer depends on the arbitrary range ``delta`` (through
       ``v_well``), and an isotropic well lets a particle make ~12 bonds where a
       single-site protein can make one.  Use this only for the RPA/length-scale
       route; use :func:`kd_to_delta` + Wertheim for anything thermodynamic.
    """
    kd_n = np.asarray(kd_M, float) * M_TO_NM3
    s = np.asarray(sigma_nm, float)
    v_well = (4 * np.pi / 3) * ((s + delta_nm) ** 3 - s**3)
    return np.log1p(1.0 / np.maximum(kd_n * v_well, 1e-300))


def epsilon_to_kd(eps_kT, sigma_nm, delta_nm=0.5):
    """Inverse of :func:`kd_to_epsilon` (isotropic square well)."""
    s = np.asarray(sigma_nm, float)
    v_well = (4 * np.pi / 3) * ((s + delta_nm) ** 3 - s**3)
    return 1.0 / (v_well * np.expm1(np.asarray(eps_kT, float)) * M_TO_NM3)


def b2_from_kd(kd_M, sigma_nm):
    r"""Osmotic second virial coefficient (nm^3) for hard spheres + one bond.

    .. math::  B_2 = \tfrac{2\pi}{3}\sigma^3 \;-\; \tfrac12 K_a,
               \qquad K_a = 1/K_d^{(n)}

    Exact at the two-body level and directly comparable to SLS / SAXS
    measurements, which makes it the cleanest experimental check of the mapping:
    if a measured B22 disagrees with this, the isotropic reduction is wrong for
    that pair (usually because the interaction is strongly anisotropic).
    """
    s = np.asarray(sigma_nm, float)
    ka = 1.0 / np.maximum(np.asarray(kd_M, float) * M_TO_NM3, 1e-300)
    return (2 * np.pi / 3) * s**3 - 0.5 * ka


# --------------------------------------------------------------------------- #
#  transport
# --------------------------------------------------------------------------- #
def diffusion(radius_nm, viscosity_Pa_s=5e-3, T_K=310.0):
    r"""Stokes-Einstein diffusion coefficient in um^2/s.

    .. math::  D = \frac{k_BT}{6\pi\eta R}

    The default viscosity (5 mPa.s, ~5x water) is a **placeholder** for
    cytoplasm; measured effective viscosities range from ~1x water for small
    solutes to >100x for large complexes, and cytoplasm is viscoelastic rather
    than Newtonian.  A 2.5 nm protein at 5 mPa.s gives D = 5.6 um^2/s.
    """
    kB = 1.380649e-23
    R_m = np.asarray(radius_nm, float) * 1e-9
    D_m2 = kB * T_K / (6 * np.pi * viscosity_Pa_s * R_m)
    return D_m2 * 1e12                                  # m^2/s -> um^2/s


# --------------------------------------------------------------------------- #
#  mRNA -> protein  (a scaling with a stated spread, NOT a prediction)
# --------------------------------------------------------------------------- #
def protein_from_mrna(tpm, copies_per_tpm=200.0, log10_sd=0.5,
                      cell_volume_um3=2000.0, n_samples=0, rng=None):
    r"""Turn a transcript level into a **concentration band**, not a number.

    ``copies_per_tpm`` is the median protein copies per TPM; 200 is a rough
    mammalian value (a 10 TPM transcript ~ 2000 copies ~ 1.7 nM in 2000 um^3).
    It varies by two orders of magnitude between genes with translation rate and
    protein half-life.

    ``log10_sd = 0.5`` encodes the residual scatter of protein about mRNA
    (~0.5 decades, i.e. a 3x band either way), consistent with the 40-60 % of
    protein variance that transcript level typically explains across genes.

    Returns ``(median_uM, lo_uM, hi_uM)`` for the central 95 % band, and the
    samples too if ``n_samples > 0``.  **Every downstream phase prediction made
    from mRNA should be reported as a probability over this band.**
    """
    tpm = np.asarray(tpm, float)
    med = copies_to_conc(tpm * copies_per_tpm, cell_volume_um3)
    lo = med * 10 ** (-1.96 * log10_sd)
    hi = med * 10 ** (+1.96 * log10_sd)
    if n_samples:
        rng = np.random.default_rng() if rng is None else rng
        s = med[..., None] * 10 ** (rng.normal(0, log10_sd, med.shape + (n_samples,)))
        return med, lo, hi, s
    return med, lo, hi
