"""
Endpoint #2 — mammalian DNA-damage response driving a GADD45a-GFP reporter.

Same seams as :class:`~genotox.core.SOSCore`, genuinely different biology.
Four differences that are not cosmetic:

1. **Delayed negative feedback.**  p53 induces Mdm2, which degrades p53, via
   transcription plus nuclear import.  That delay makes the damaged state
   *oscillate* rather than settle, so "p53 level" is not a well-defined
   quantity and the reporter integrates a pulse train.
2. **The reporter integrates.**  GFP matures slowly and is stable, so
   fluorescence tracks the time-integral of p53 rather than its amplitude.
   Two exposures with equal AUC but different pulse structure are
   indistinguishable at the plate reader.
3. **Aneugens enter through a different door.**  Spindle interference is not
   a DNA lesion, so it must not be projected onto the lesion variable — a
   comet readout calibrated against ``lesions`` would then report strand
   breaks that do not exist.  It gets its own state variable and reaches p53
   through the mitotic-surveillance route, which requires cells to be
   dividing.
4. **Cytostasis is endogenous.**  p53 -> p21 arrests the cycle, so the assay
   partly causes its own growth-gate failures.  In the bacterial core all
   growth inhibition came from outside.

Rigour label: topology is standard (p53-Mdm2 delayed feedback is textbook and
the pulsatile behaviour is a well-replicated experimental result).  The rate
constants here were tuned by parameter sweep to put the damaged-state pulse
period near the commonly reported few-hour range and to keep the undamaged
state stable; they are illustrative (H) and are NOT fitted to any published
time series.  Absolute EC values are not predictions.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .core import Exposure, SignalCore
from .damage import DamageFlux

_ZERO_FLUX = DamageFlux()


@dataclass
class P53Mdm2:
    """The p53-Mdm2 delayed-feedback loop, as a reusable 4-state sub-model.

    Extracted so the cytogenetic core (endpoint #3) can share one p53
    representation with the reporter core rather than carrying a second,
    silently different one.  It owns states (p53, Mdm2_mRNA, Mdm2_cyt,
    Mdm2_nuc) and knows nothing about what reads it.

    Parameters were chosen by sweep to put the damaged-state pulse period
    near the commonly reported few-hour range (~5.3 h) and keep the
    undamaged state stable.  Illustrative (H), not fitted to a time series.
    """

    n_states: int = 4
    s_p53: float = 0.06
    k_deg: float = 0.60
    K_deg: float = 0.05
    d_p53_basal: float = 0.0015
    K_atm: float = 0.30      # ATM shielding p53 from Mdm2-mediated degradation
    k_mdm2_txn: float = 0.030
    mdm2_leak: float = 0.04
    K_p53: float = 0.35
    h_p53: float = 6.0
    d_mdm2_mrna: float = 0.006
    k_mdm2_tsl: float = 0.060
    k_import: float = 0.015
    d_mdm2: float = 0.006
    psi_atm: float = 3.0     # ATM also destabilising Mdm2
    # p53 level giving 50% cycle arrest.  Set too low, a reporter core
    # arrests the culture before it can accumulate signal and every
    # genotoxicant self-gates into INCONCLUSIVE.
    K_arrest: float = 0.60
    h_arrest: float = 2.0

    seed: tuple = (0.05, 0.5, 1.0, 1.0)

    def rhs(self, y4, A: float):
        """dy/dt for the four states, given DDR kinase activity ``A``."""
        P, Mr, Mc, Mn = (max(v, 0.0) for v in y4)
        deg = self.k_deg * Mn * P / (P + self.K_deg) / (1.0 + A / self.K_atm)
        dP = self.s_p53 - deg - self.d_p53_basal * P
        prom = (self.mdm2_leak + (1.0 - self.mdm2_leak)
                * _hill(P, self.K_p53, self.h_p53))
        dMr = self.k_mdm2_txn * prom - self.d_mdm2_mrna * Mr
        dMc = (self.k_mdm2_tsl * Mr - self.k_import * Mc
               - self.d_mdm2 * Mc * (1.0 + self.psi_atm * A))
        dMn = self.k_import * Mc - self.d_mdm2 * Mn * (1.0 + self.psi_atm * A)
        return (dP, dMr, dMc, dMn)

    def arrest(self, P: float) -> float:
        """Fraction of maximal cycling rate left, via p53 -> p21."""
        return 1.0 / (1.0 + (max(P, 0.0) / self.K_arrest) ** self.h_arrest)


def _hill(x, K, h):
    x = max(x, 0.0)
    return x ** h / (K ** h + x ** h) if x > 0 else 0.0


@dataclass
class P53Core(SignalCore):
    """DNA damage / mitotic stress -> ATM-ATR -> p53-Mdm2 -> GADD45a-GFP."""

    name = "p53/GADD45a-GFP"
    state_names = ("D", "S", "p53", "Mdm2_mRNA", "Mdm2_cyt", "Mdm2_nuc",
                   "gfp_mRNA", "gfp_immature", "gfp_mature", "N")

    # -- channel weights ----------------------------------------------------
    # Contrast with SOSCore.sos_weights: the bacterial core is dominated by
    # replication-blocking bulky lesions, ATM by frank double-strand breaks.
    # Same upstream vector, different projection, different answer.
    dna_weights: dict = field(default_factory=lambda: {
        "dsb": 1.00, "icl": 0.90, "topo": 0.85, "bulky_adduct": 0.60,
        "ssb": 0.40, "alkylation": 0.35, "oxidative": 0.25,
        "aneugenic": 0.00,          # never a DNA lesion — see point 3 above
    })
    mitotic_weights: dict = field(default_factory=lambda: {
        "aneugenic": 1.00,
    })

    # -- lesion turnover ----------------------------------------------------
    v_repair: float = 1.2
    k_repair: float = 2.0
    k_repair_lin: float = 0.05
    # Spindle engagement is a fast-equilibrating occupancy, not an
    # accumulating lesion: drug on, drug off.  Modelled as a slow accumulator
    # it saturates the sensor at the lowest dose tested and every aneugen
    # returns the same induction ratio across three decades of dose.
    k_clear_mitotic: float = 0.50

    # -- damage sensing (ATM/ATR), treated as fast/algebraic ----------------
    K_dna: float = 0.30
    h_dna: float = 2.0
    K_mit: float = 2.0
    h_mit: float = 2.0
    w_dna: float = 1.0
    w_mit: float = 0.55   # mitotic surveillance is a real but weaker route

    # -- p53 / Mdm2 loop, shared with the cytogenetic core ------------------
    p53: P53Mdm2 = field(default_factory=P53Mdm2)

    # -- GADD45a-GFP reporter ------------------------------------------------
    gfp_leak: float = 0.05
    K_gadd: float = 0.25
    h_gadd: float = 3.0
    k_gfp_txn: float = 1.0
    d_gfp_mrna: float = 0.010
    k_gfp_tsl: float = 1.0
    k_maturation: float = 0.030   # GFP folding/oxidation; the reporter's lag
    d_gfp: float = 0.0            # stable; lost only by division

    # -- growth / viability ---------------------------------------------------
    mu_max: float = np.log(2) / 1440.0   # 24 h doubling
    K_tox: float = 25.0
    h_tox: float = 3.0
    N0: float = 1.0

    def __post_init__(self):
        self._basal = self._solve_basal()

    # -- helpers --------------------------------------------------------------
    def _hill(self, x, K, h):
        x = max(x, 0.0)
        return x ** h / (K ** h + x ** h) if x > 0 else 0.0

    def viability(self, D: float, extra_tox: float) -> float:
        v_dmg = 1.0 / (1.0 + (max(D, 0.0) / self.K_tox) ** self.h_tox)
        return v_dmg / (1.0 + max(extra_tox, 0.0))

    def gadd_promoter(self, P: float) -> float:
        """GADD45a promoter activity: a steep switch, not a linear sensor.

        This is why the reporter does not integrate p53 itself.  p53 excursions
        that stay below K_gadd contribute area to the p53 integral but almost
        nothing to the promoter integral, so the assay is structurally blind
        to low-level sustained p53 — see :func:`genotox.run_p53.integrator_check`.
        """
        return (self.gfp_leak + (1.0 - self.gfp_leak)
                * self._hill(P, self.K_gadd, self.h_gadd))

    def ddr_activity(self, D: float, S: float, cycling: float) -> float:
        """Combined ATM/ATR activity, capped at 1.

        ``cycling`` gates the mitotic route: a cell that is not dividing
        cannot trip mitotic surveillance.  It is deliberately computed from
        the damage-independent growth terms only, so this stays algebraic
        rather than forming a loop with p53-driven arrest.
        """
        a = (self.w_dna * self._hill(D, self.K_dna, self.h_dna)
             + self.w_mit * self._hill(S, self.K_mit, self.h_mit) * cycling)
        return min(a, 1.0)

    def _solve_basal(self) -> np.ndarray:
        """Undamaged steady state, obtained by relaxing the model."""
        exp = Exposure(flux=_ZERO_FLUX, duration_min=20000.0)
        y = np.array([0.0, 0.0, *self.p53.seed, 0.1, 0.1, 1.0, self.N0])
        from scipy.integrate import solve_ivp
        sol = solve_ivp(self.rhs, (0.0, 20000.0), y, args=(exp,),
                        method="LSODA", rtol=1e-9, atol=1e-11)
        y = sol.y[:, -1].copy()
        y[0] = y[1] = 0.0
        y[-1] = self.N0
        return y

    # -- model ----------------------------------------------------------------
    def basal_state(self) -> np.ndarray:
        return self._basal.copy()

    def rhs(self, t, y, exp: Exposure):
        D, S, P, Mr, Mc, Mn, Rg, Gi, Gm, N = y
        D, S, P = max(D, 0.0), max(S, 0.0), max(P, 0.0)

        v = self.viability(D, exp.extra_toxicity)
        cycling = v / (1.0 + max(exp.growth_inhibition, 0.0))
        # p53 -> p21 -> cycle arrest: the assay's own cytostasis
        mu = self.mu_max * cycling * self.p53.arrest(P)

        A = self.ddr_activity(D, S, cycling)

        repair = (self.v_repair * D / (self.k_repair + D)
                  + self.k_repair_lin * D)
        dD = exp.flux.project(self.dna_weights) - repair - mu * D
        dS = (exp.flux.project(self.mitotic_weights)
              - self.k_clear_mitotic * S - mu * S)

        dP, dMr, dMc, dMn = self.p53.rhs((P, Mr, Mc, Mn), A)

        prom_gadd = self.gadd_promoter(P)
        dRg = self.k_gfp_txn * prom_gadd * v - self.d_gfp_mrna * Rg - mu * Rg
        dGi = self.k_gfp_tsl * Rg * v - self.k_maturation * Gi - mu * Gi
        dGm = self.k_maturation * Gi - self.d_gfp * Gm - mu * Gm
        dN = mu * N
        return np.array([dD, dS, dP, dMr, dMc, dMn, dRg, dGi, dGm, dN])

    def observables(self, t, Y, exp):
        D, S, P, Mr, Mc, Mn, Rg, Gi, Gm, N = Y
        v = np.array([self.viability(d, exp.extra_toxicity) for d in D])
        cycling = v / (1.0 + max(exp.growth_inhibition, 0.0))
        A = np.array([self.ddr_activity(d, s, c)
                      for d, s, c in zip(D, S, cycling)])
        return {
            "t": t,
            "lesions": D,
            "mitotic_stress": S,
            "ddr_active": A,
            "p53": P,
            "mdm2": Mn,
            "gadd_promoter": np.array([self.gadd_promoter(x) for x in P]),
            "gfp_mRNA": Rg,
            "gfp_immature": Gi,
            # readout-layer contract
            "mature_gfp": Gm,
            "density": N,
            "viability": v,
        }
