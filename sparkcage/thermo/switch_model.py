"""
SparkCage three-state thermodynamic switch model.

Implements the coupled-equilibrium model behind quote item Phase I / Step 1 / (3):
"Free energy calculation for three thermodynamic states, dG_open, dG_LT, dG_switch".

This is the LOCKR / LucCage design principle (Langan 2019 Nature; Quijano-Rubio
2021 Nat Biotechnol) reduced to the single-component ELECTROCHEMICAL readout
used by SparkCage (no separate "key" component: the readout is the change in
methylene-blue / electrode electron-transfer distance when the latch is released).

States
------
    C          latch docked in the cage, methylene blue held close to the
               electrode  -> signal-ON for a "signal-off on binding" design,
               or signal-OFF for the inverse; the model is agnostic.
    O          latch thermally released, no analyte bound.
    OT         latch released and bound to the analyte (target).

Equilibria
----------
    C  <-> O          dG_open   ( > 0 ; cage holds the latch shut )
    O + T <-> OT      dG_LT     ( < 0 ; binder-analyte affinity )
    dG_switch = dG_open + dG_LT     net drive to switch at 1 M analyte

Population of the released (switched) sensor at analyte concentration T:

    K_open = exp(-dG_open / RT)
    K_LT   = exp(-dG_LT   / RT)  = 1 / Kd_LT

    f_switched(T) = ( K_open + K_open*K_LT*T ) / ( 1 + K_open + K_open*K_LT*T )

Analytic consequences (the whole design trade-off in two lines):

    EC50        = (1 + K_open) / (K_open * K_LT)  ~=  Kd_LT * exp(dG_open/RT)
    dyn. range  = f(inf)/f(0)   = (1 + K_open)/K_open  ~=  exp(dG_open/RT)

i.e. every kcal/mol added to dG_open buys ~5.5x dynamic range at 37 C and costs
exactly the same factor in limit of detection. That trade-off is the reason the
quote budgets a separate free-energy work package.

All energies in kcal/mol, concentrations in molar.
"""

from __future__ import annotations

import numpy as np

R_KCAL = 1.987204259e-3  # kcal / (mol K)


def rt(temperature_k: float = 310.15) -> float:
    """RT in kcal/mol. Default 37 C."""
    return R_KCAL * temperature_k


def k_from_dg(dg: float | np.ndarray, temperature_k: float = 310.15):
    """Equilibrium constant from a free energy in kcal/mol."""
    return np.exp(-np.asarray(dg, dtype=float) / rt(temperature_k))


def dg_from_kd(kd_molar: float | np.ndarray, temperature_k: float = 310.15):
    """Binding free energy (kcal/mol, negative) from a dissociation constant."""
    return rt(temperature_k) * np.log(np.asarray(kd_molar, dtype=float))


def fraction_switched(target_m, dg_open, dg_lt, temperature_k: float = 310.15):
    """Fraction of sensor in the released (O + OT) states.

    Parameters
    ----------
    target_m : free analyte concentration, molar (scalar or array)
    dg_open  : cage-latch opening free energy, kcal/mol (positive = caged)
    dg_lt    : latch/binder-analyte binding free energy, kcal/mol (negative)
    """
    t = np.asarray(target_m, dtype=float)
    k_open = k_from_dg(dg_open, temperature_k)
    k_lt = k_from_dg(dg_lt, temperature_k)
    num = k_open * (1.0 + k_lt * t)
    return num / (1.0 + num)


def dg_switch(dg_open: float, dg_lt: float) -> float:
    """Net switching free energy at 1 M analyte."""
    return dg_open + dg_lt


def ec50(dg_open, dg_lt, temperature_k: float = 310.15):
    """Analyte concentration giving half-maximal *change* in fraction switched."""
    k_open = k_from_dg(dg_open, temperature_k)
    k_lt = k_from_dg(dg_lt, temperature_k)
    return (1.0 + k_open) / (k_open * k_lt)


def dynamic_range(dg_open, temperature_k: float = 310.15):
    """f(saturating analyte) / f(no analyte) -- the fold signal change."""
    k_open = k_from_dg(dg_open, temperature_k)
    return (1.0 + k_open) / k_open


def optimal_dg_open(target_window_m, dg_lt, temperature_k: float = 310.15,
                    grid=None):
    """dG_open maximising signal change across a required analyte window.

    Picks dG_open that maximises  f(T_high) - f(T_low)  , i.e. the usable swing
    over the clinically required concentration range.
    """
    t_low, t_high = target_window_m
    if grid is None:
        grid = np.linspace(0.0, 12.0, 2401)
    swing = (fraction_switched(t_high, grid, dg_lt, temperature_k)
             - fraction_switched(t_low, grid, dg_lt, temperature_k))
    i = int(np.argmax(swing))
    return float(grid[i]), float(swing[i])


def sensor_table(dg_open, dg_lt, target_window_m=None, temperature_k=310.15):
    """One-row summary of a candidate design."""
    row = {
        "dG_open_kcal": float(dg_open),
        "dG_LT_kcal": float(dg_lt),
        "dG_switch_kcal": float(dg_switch(dg_open, dg_lt)),
        "Kd_binder_M": float(np.exp(np.asarray(dg_lt) / rt(temperature_k))),
        "EC50_M": float(ec50(dg_open, dg_lt, temperature_k)),
        "fold_dynamic_range": float(dynamic_range(dg_open, temperature_k)),
        "baseline_open_fraction": float(fraction_switched(0.0, dg_open, dg_lt,
                                                          temperature_k)),
    }
    if target_window_m is not None:
        lo, hi = target_window_m
        row["swing_over_window"] = float(
            fraction_switched(hi, dg_open, dg_lt, temperature_k)
            - fraction_switched(lo, dg_open, dg_lt, temperature_k))
    return row
