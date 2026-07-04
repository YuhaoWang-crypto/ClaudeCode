"""Mechanistic delivery-kinetics model (PLAN.md Phase 4).

A compartmental ODE model of the LNP → cell → protein cascade, following the
sequential first-order transfer picture of Müller et al. 2024 ("Kinetics of RNA-LNP
delivery and protein expression", Eur J Pharm Biopharm) and the LNP-siRNA ODE
models of Mihaila et al. 2017/2019. It turns biophysical rates (uptake, endosomal
escape, degradation, translation) into an interpretable expression time-course and
summary metrics — complementing the ML potency score from Track A.

States (amounts, arbitrary per-cell / per-dose units):
    L_ex    extracellular LNP (dose)
    L_endo  endosomal LNP (internalised)
    M       cytosolic free mRNA (escaped + released)
    P       protein

Fluxes (all first-order, /h):
    L_ex   --k_uptake-->  L_endo                      (corona-mediated adhesion + endocytosis)
    L_endo --k_escape-->  M   (x mRNA_per_lnp copies) (endosomal escape → cytosolic release)
    L_endo --k_lyso-->    ∅                            (lysosomal degradation)
    L_endo --k_recycle--> ∅                            (endosomal recycling / exocytosis, lost)
    M      --k_mdeg-->    ∅                            (cytosolic mRNA decay)
    M      --k_transl-->  P (production)               (translation)
    P      --k_pdeg-->    ∅                            (protein turnover)

Endosomal-escape EFFICIENCY = k_escape / (k_escape + k_lyso + k_recycle) — the
well-known ~1–2% bottleneck; this is the step the ionizable lipid (Track A) governs.

Because the cascade is linear, the total protein exposure has a closed form:
    AUC_P = dose · f_uptake · f_escape · mRNA_per_lnp · k_transl / (k_mdeg · k_pdeg)
where f_uptake = 1 (all L_ex eventually internalises here) and
      f_escape = k_escape/(k_escape+k_lyso+k_recycle).
`analytic_auc` returns this; `simulate` integrates the full time-course.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass
class DeliveryParams:
    # rates in per hour; literature-informed defaults (reporter = firefly luciferase)
    k_uptake: float = 0.3        # ~half internalised in ~2–3 h
    k_escape: float = 0.02       # slow escape …
    k_lyso: float = 1.0          # … vs fast lysosomal degradation
    k_recycle: float = 0.1       # endosomal recycling / exocytosis
    k_mdeg: float = 0.058        # cytosolic mRNA t½ ≈ 12 h
    k_transl: float = 10.0       # protein produced per mRNA per h (scale)
    k_pdeg: float = 0.23         # protein t½ ≈ 3 h (Fluc-like)
    mRNA_per_lnp: float = 100.0  # mRNA copies released per escaped LNP (yield)
    dose: float = 1.0            # extracellular LNP at t=0

    def escape_efficiency(self) -> float:
        return self.k_escape / (self.k_escape + self.k_lyso + self.k_recycle)


def _rhs(t, y, p: DeliveryParams):
    L_ex, L_endo, M, P = y
    dL_ex = -p.k_uptake * L_ex
    dL_endo = p.k_uptake * L_ex - (p.k_escape + p.k_lyso + p.k_recycle) * L_endo
    dM = p.mRNA_per_lnp * p.k_escape * L_endo - p.k_mdeg * M
    dP = p.k_transl * M - p.k_pdeg * P
    return [dL_ex, dL_endo, dM, dP]


def simulate(p: DeliveryParams, t_end: float = 96.0, n: int = 400) -> pd.DataFrame:
    """Integrate the cascade. Returns a DataFrame (t, L_ex, L_endo, M, P)."""
    t = np.linspace(0, t_end, n)
    sol = solve_ivp(_rhs, (0, t_end), [p.dose, 0.0, 0.0, 0.0], t_eval=t,
                    args=(p,), method="LSODA", rtol=1e-8, atol=1e-10)
    df = pd.DataFrame({"t": sol.t, "L_ex": sol.y[0], "L_endo": sol.y[1],
                       "M": sol.y[2], "P": sol.y[3]})
    return df


def analytic_auc(p: DeliveryParams) -> float:
    """Closed-form total protein exposure ∫P dt over all time."""
    return (p.dose * p.escape_efficiency() * p.mRNA_per_lnp
            * p.k_transl / (p.k_mdeg * p.k_pdeg))


def metrics(p: DeliveryParams, t_end: float = 96.0) -> dict:
    df = simulate(p, t_end=t_end)
    P = df["P"].values
    auc = float(np.trapezoid(P, df["t"].values))
    imax = int(np.argmax(P))
    return {
        "escape_efficiency": p.escape_efficiency(),
        "peak_protein": float(P[imax]),
        "t_peak_h": float(df["t"].values[imax]),
        "auc_protein_numeric": auc,
        "auc_protein_analytic": analytic_auc(p),
        # mRNA that ever reaches the cytosol (delivered dose fraction × copies)
        "cytosolic_mRNA_auc": float(np.trapezoid(df["M"].values, df["t"].values)),
    }
