"""
Exact reimplementation of the published LucCage coupled-equilibrium model.

Source: Quijano-Rubio, Yeh, Park et al., "De novo design of modular and tunable
protein biosensors", Nature 591:482-487 (2021), Supplementary Methods section 1
and the authors' released simulation code
    https://files.ipd.uw.edu/pub/de_novo_design_of_tunable_biosensors_2021/model_simulation.py

Why this file exists
--------------------
The quote's Step 1.3 asks for "free energy calculation for three thermodynamic
states, dG_open, dG_LT, dG_switch".  Two corrections fall out of reading the
actual literature:

1.  There is NO term called dG_switch in either the LucCage paper or the LOCKR
    paper that precedes it.  The published model has FOUR free energies:
        dG_open  cage -> open latch          (intramolecular, K_open dimensionless)
        dG_LT    latch/binder -> target      (K_LT is a dissociation constant, M)
        dG_CK    open cage -> key            (K_CK dissociation constant, M)
        dG_R     SmBiT -> LgBiT reconstitution (K_R = Kd_NanoBiT / C_eff = 0.19)

2.  dG_CK and dG_R exist only because LucCage is a TWO-COMPONENT, split-luciferase
    system.  An electrochemical SparkCage has no key and no luciferase, so both
    terms vanish -- and that is not a simplification, it is a loss.  See
    `key_compensation_penalty` below.

Species (names follow the authors' code)
    CLc    closed cage-latch                CLo    open cage-latch
    K      free key                         T      free target
    CLT    open cage-latch . target
    CLK    cage-latch . key, luciferase NOT reconstituted
    CLKT   cage-latch . key . target, NOT reconstituted
    CLKR   cage-latch . key, reconstituted        <- background signal
    CLKTR  cage-latch . key . target, reconstituted <- on signal

Observable (verbatim from the authors' code comment):
    fraction of reconstituted luciferase = (CLKTR + CLKR) / [lucCage]_total
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, fsolve

R_KCAL = 1.987204259e-3

# Baseline parameters, Quijano-Rubio 2021 Extended Data Fig. 1
DEFAULTS = dict(k_open=1.0e-3, k_ck=1.0e-8, k_lt=1.0e-9, k_r=0.19,
                cage_total=10e-9, key_total=100e-9)


def _species(clc, k_free, t_free, k_open, k_ck, k_lt, k_r):
    """All nine concentrations from the three independent ones."""
    clo = k_open * clc
    clt = clo * t_free / k_lt
    clk = k_free * clo / k_ck
    clkt = clt * k_free / k_ck
    clkr = clk / k_r
    clktr = clkt / k_r
    return clo, clt, clk, clkt, clkr, clktr


def solve(target_total, k_open=None, k_ck=None, k_lt=None, k_r=None,
          cage_total=None, key_total=None, guess=None):
    """Solve the 10-equation system at one total target concentration.

    Returns the fraction of reconstituted luciferase.  Solved in log space so
    the concentrations stay positive over the many decades the dose-response
    curve spans.
    """
    p = dict(DEFAULTS)
    for name, val in [("k_open", k_open), ("k_ck", k_ck), ("k_lt", k_lt),
                      ("k_r", k_r), ("cage_total", cage_total),
                      ("key_total", key_total)]:
        if val is not None:
            p[name] = val

    def residuals(log_v):
        clc, k_free, t_free = np.exp(log_v)
        clo, clt, clk, clkt, clkr, clktr = _species(
            clc, k_free, t_free, p["k_open"], p["k_ck"], p["k_lt"], p["k_r"])
        return [
            (clc + clo + clk + clkr + clkt + clktr + clt) / p["cage_total"] - 1.0,
            (k_free + clk + clkr + clkt + clktr) / p["key_total"] - 1.0,
            (t_free + clkt + clktr + clt) / max(target_total, 1e-30) - 1.0,
        ]

    if guess is None:
        guess = np.log([p["cage_total"] * 0.9, p["key_total"] * 0.9,
                        max(target_total, 1e-30) * 0.9])
    sol, _, ok, _ = fsolve(residuals, guess, full_output=True)
    clc, k_free, t_free = np.exp(sol)
    _, _, _, _, clkr, clktr = _species(clc, k_free, t_free, p["k_open"],
                                       p["k_ck"], p["k_lt"], p["k_r"])
    return {"fraction_reconstituted": float((clkr + clktr) / p["cage_total"]),
            "converged": bool(ok == 1), "log_guess": sol}


def dose_response(target_totals, **kw):
    """Fraction of reconstituted luciferase across a target titration."""
    out, guess = [], None
    for t in np.atleast_1d(target_totals):
        r = solve(float(t), guess=guess, **kw)
        guess = r["log_guess"]
        out.append(r["fraction_reconstituted"])
    return np.array(out)


def dynamic_range_percent(target_saturating=1e-5, **kw):
    """(R_max - R_min) / R_min * 100, the paper's own definition."""
    r_min = solve(1e-15, **kw)["fraction_reconstituted"]
    r_max = solve(target_saturating, **kw)["fraction_reconstituted"]
    return float((r_max - r_min) / r_min * 100.0)


# ---------------------------------------------------------------------------
# What the electrochemical architecture gives up
# ---------------------------------------------------------------------------

def key_compensation_factor(k_ck=None, key_total=None):
    """How much cage-opening the key pays for, expressed as a fold factor.

    LucCage opens when   dG_open - dG_CK - dG_LT << 0.
    The key's binding energy (negative dG_CK) pays part of the cost of opening
    the cage, so the target supplies only the remainder.  A single-component
    electrochemical sensor has no key: its condition is just
    dG_open - dG_LT << 0, and the target pays the full cost alone.

    NOTE ON A TEMPTING ERROR: the compensation is NOT exp(-dG_CK/RT), which for
    K_CK = 1e-8 M would suggest a factor of 1e8.  dG_CK is a standard-state free
    energy at 1 M key, and the key is present at 100 nM.  The compensation that
    is actually available is the mass-action term

        1 + [key] / K_CK

    which for the paper's 100 nM key and K_CK = 1e-8 M is only about 11-fold.
    Measured end-to-end against the full model (see validate.py) the penalty for
    dropping the key is roughly 40x in EC50 -- real and worth designing around,
    but four orders of magnitude smaller than the naive estimate.

    The design response is that SparkCage needs its own compensating interaction.
    The quote already implies one: item 1.2 budgets molecular dynamics for
    "binding sites on Latch to carbon electrodes".  That latch-electrode
    adsorption energy is the structural analogue of dG_CK and should be modelled
    as an explicit state, not left implicit.
    """
    p = dict(DEFAULTS)
    if k_ck is not None:
        p["k_ck"] = k_ck
    if key_total is not None:
        p["key_total"] = key_total
    return float(1.0 + p["key_total"] / p["k_ck"])


def dg_from_k(k, temperature_k=298.15):
    """Free energy from a dimensionless equilibrium constant or a Kd in M."""
    return float(-R_KCAL * temperature_k * np.log(k))
