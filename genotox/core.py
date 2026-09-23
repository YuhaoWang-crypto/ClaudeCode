"""
Signal-transduction cores.

A core owns the biology between a :class:`~genotox.damage.DamageFlux` and a
set of time-resolved internal variables.  It knows nothing about how the
damage was predicted and nothing about how the result will be measured.

Adding endpoint #2 (a p53/GADD45a mammalian reporter) means writing another
:class:`SignalCore`; adding endpoint #3 (comet) means writing another readout
against a core that already exposes a break-density observable.  Neither
touches the upstream layer.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from .damage import DamageFlux


@dataclass
class Exposure:
    """Everything a core needs to know about one well."""

    flux: DamageFlux
    extra_toxicity: float = 0.0      # lethal: suppresses biosynthesis + growth
    growth_inhibition: float = 0.0   # bacteriostatic: suppresses growth only
    s9: bool = False
    duration_min: float = 120.0


class SignalCore:
    """Base class for a damage -> signal model."""

    name = "abstract"
    state_names: tuple = ()

    def basal_state(self) -> np.ndarray:
        """Steady state under zero damage; also the t=0 condition."""
        raise NotImplementedError

    def rhs(self, t: float, y: np.ndarray, exp: Exposure) -> np.ndarray:
        raise NotImplementedError

    def observables(self, t: np.ndarray, Y: np.ndarray, exp: Exposure) -> dict:
        """Named trajectories for the readout layer to consume."""
        raise NotImplementedError

    def simulate(self, exp: Exposure, n_points: int = 241) -> dict:
        y0 = self.basal_state()
        t_eval = np.linspace(0.0, exp.duration_min, n_points)
        sol = solve_ivp(self.rhs, (0.0, exp.duration_min), y0, args=(exp,),
                        t_eval=t_eval, method="LSODA",
                        rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"{self.name} integration failed: {sol.message}")
        return self.observables(sol.t, sol.y, exp)


@dataclass
class SOSCore(SignalCore):
    """LexA/RecA SOS regulon driving a umuDC-lacZ reporter in S. typhimurium
    TA1535/pSK1002.

    State (5): lesion burden D, LexA repressor L, reporter mRNA M,
    beta-galactosidase per cell B, culture density N.

    Topology (matches the standard picture of the assay):
        lesion -> ssDNA -> RecA* -> LexA autocleavage -> umuDC promoter
        de-repression -> lacZ transcription/translation -> beta-gal.
    LexA also represses its own gene, giving the negative feedback that sets
    the basal set-point and the graded (not all-or-none) dose response.

    Rigour label: the topology is standard and mechanistic; the rate constants
    are plausible order-of-magnitude values (H), NOT fitted to a published
    induction dataset.  Absolute ECIR1.5 values from this core are therefore
    not predictions — relative behaviour and curve shape are what it is for.
    """

    name = "SOS/umuDC-lacZ"
    state_names = ("D", "LexA", "mRNA", "betagal", "N")

    # --- which lesion channels actually make ssDNA at a stalled fork -------
    # This dict is the whole reason the damage layer is a vector.  An
    # aneugen scores 0 here and the core is therefore structurally unable to
    # call it positive, which is the correct biology.
    sos_weights: dict = field(default_factory=lambda: {
        "bulky_adduct": 1.00,
        "icl": 1.00,
        "topo_bacterial": 0.90,
        "topo_mammalian": 0.00,
        "dsb": 0.80,
        "ssb": 0.60,
        "oxidative": 0.35,
        "alkylation": 0.25,
        "aneugenic": 0.00,
    })

    # --- lesion turnover ---------------------------------------------------
    v_repair: float = 1.2        # saturable repair Vmax (lesions/min)
    k_repair: float = 2.0        # its half-saturation (lesions)
    k_repair_lin: float = 0.05   # first-order repair/1-min^-1

    # --- RecA activation ---------------------------------------------------
    K_recA: float = 1.0
    h_recA: float = 2.0

    # --- LexA ---------------------------------------------------------------
    beta_lexA: float = 1.0       # max synthesis
    K_lexA: float = 1.0          # self-repression constant
    h_lexA: float = 2.0
    k_cleave: float = 3.0        # RecA*-stimulated autocleavage

    # --- reporter -----------------------------------------------------------
    # kumu_ratio = K_umu / L_basal.  Small ratio = tightly repressed "late"
    # SOS gene (umuDC).  Raising it toward 1 turns the same core into an
    # "early" gene (recA-like).  This one number is the gene-class knob.
    kumu_ratio: float = 0.25
    h_umu: float = 2.0
    leak: float = 0.03           # basal promoter activity fraction
    k_txn: float = 1.0
    d_mrna: float = 0.30
    k_tsl: float = 1.0
    d_prot: float = 0.0          # beta-gal is stable; lost only by dilution

    # --- growth / viability -------------------------------------------------
    mu_max: float = np.log(2) / 30.0   # 30 min doubling
    K_tox: float = 25.0                # lesion burden at 50% viability
    h_tox: float = 3.0
    N0: float = 0.05                   # OD600 at start of exposure

    def __post_init__(self):
        self._L_basal = self._solve_L_basal()
        self.K_umu = self.kumu_ratio * self._L_basal

    # -- helpers -------------------------------------------------------------
    def _solve_L_basal(self) -> float:
        """LexA steady state with no damage: synthesis = dilution."""
        from scipy.optimize import brentq
        f = lambda L: (self.beta_lexA / (1.0 + (L / self.K_lexA) ** self.h_lexA)
                       - self.mu_max * L)
        return brentq(f, 1e-6, 1e4)

    def promoter(self, L: float) -> float:
        return self.leak + (1.0 - self.leak) / (
            1.0 + (L / self.K_umu) ** self.h_umu)

    def viability(self, D: float, extra_tox: float) -> float:
        """Fraction of maximal biosynthetic/growth capacity.

        Two independent routes: lethal lesion burden, and damage-independent
        cytotoxicity carried through from the compound.  Without this term the
        model produces monotone induction curves and so cannot reproduce the
        roll-off that every real umu dose series shows at high dose.
        """
        v_dmg = 1.0 / (1.0 + (max(D, 0.0) / self.K_tox) ** self.h_tox)
        return v_dmg / (1.0 + max(extra_tox, 0.0))

    def influx(self, flux: DamageFlux) -> float:
        return flux.project(self.sos_weights)

    # -- model ---------------------------------------------------------------
    def basal_state(self) -> np.ndarray:
        L = self._L_basal
        M = self.k_txn * self.promoter(L) / (self.d_mrna + self.mu_max)
        B = self.k_tsl * M / (self.d_prot + self.mu_max)
        return np.array([0.0, L, M, B, self.N0])

    def rhs(self, t, y, exp: Exposure):
        D, L, M, B, N = y
        D = max(D, 0.0)

        v = self.viability(D, exp.extra_toxicity)
        # Growth is suppressed by everything that suppresses biosynthesis,
        # plus anything purely bacteriostatic.  The asymmetry matters: a
        # growth-arrested but translating cell keeps accumulating reporter
        # (no dilution), which is a documented source of inflated induction
        # ratios in the real assay — and the reason the growth gate exists.
        mu = self.mu_max * v / (1.0 + max(exp.growth_inhibition, 0.0))

        repair = (self.v_repair * D / (self.k_repair + D)
                  + self.k_repair_lin * D)
        dD = self.influx(exp.flux) - repair - mu * D

        recA = D ** self.h_recA / (self.K_recA ** self.h_recA + D ** self.h_recA)
        dL = (self.beta_lexA / (1.0 + (L / self.K_lexA) ** self.h_lexA)
              - self.k_cleave * recA * L - mu * L)

        dM = self.k_txn * self.promoter(L) * v - self.d_mrna * M - mu * M
        dB = self.k_tsl * M * v - self.d_prot * B - mu * B
        dN = mu * N
        return np.array([dD, dL, dM, dB, dN])

    def observables(self, t, Y, exp):
        D, L, M, B, N = Y
        v = np.array([self.viability(d, exp.extra_toxicity) for d in D])
        return {
            "t": t,
            "lesions": D,
            "LexA": L,
            "recA_active": D ** self.h_recA / (self.K_recA ** self.h_recA
                                               + D ** self.h_recA),
            "promoter": np.array([self.promoter(x) for x in L]),
            "mRNA": M,
            # readout layer contract: enzyme per cell + culture density
            "reporter_enzyme": B,
            "density": N,
            "viability": v,
        }
