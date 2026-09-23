"""
Endpoint #3 — comet and micronucleus, as two readouts of ONE core.

This is the case the layer split was built for.  The alkaline comet assay and
the cytokinesis-block micronucleus assay are not different biologies; they are
the same damaged cells measured with different instruments at different times:

  * comet reads break density in situ, after a few hours, with no division;
  * micronucleus reads what survived into the next mitosis, ~36 h later.

So this module adds one :class:`SignalCore` and two readouts, and reuses the
upstream layer, the p53 sub-model and the decision layer unchanged.

Why the lesion pool had to be split
-----------------------------------
The p53 core carried a single lumped ``lesions`` variable.  That is enough to
drive a transcriptional reporter and not enough for comet, because the
alkaline comet does **not** see bulky adducts directly.  What migrates is
strand breaks and alkali-labile sites.  Consequently:

  * an alkylating agent is comet-positive immediately -- its adducts *are*
    alkali-labile;
  * a bulky-adduct former is comet-positive only through the transient gaps
    that excision repair itself creates, so its signal is a repair artefact
    in the literal sense, and a cell that repairs faster looks *more* damaged;
  * an aneugen is comet-negative at every dose, because nothing is broken.

Why micronuclei need a mitosis
------------------------------
A micronucleus is scored after a division: an acentric fragment or a lagging
whole chromosome fails to reach a pole.  A compound that arrests the cycle
therefore produces *fewer* micronuclei, and since p53 -> p21 arrest is already
in the shared sub-model, the assay suppresses its own endpoint at high dose.
That is why the dose-response is expected to turn over, and why the validity
measure here is proliferation (CBPI), not cell number.

Centromere status is tracked separately throughout: fragments give
centromere-negative micronuclei, lagging chromosomes give centromere-positive
ones.  That single extra observable is what lets this endpoint do what the
step-2 reporter could not -- name the mechanism rather than just the call.

Rigour label: topology and the qualitative couplings are standard.  Every
rate constant is illustrative (H), not fitted.  Absolute micronucleus
frequencies and %tail values are not predictions.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .core import Exposure, SignalCore
from .damage import DamageFlux
from .p53 import P53Mdm2, _hill

_ZERO_FLUX = DamageFlux()


@dataclass
class CytogeneticCore(SignalCore):
    """Damage processing -> strand breaks -> chromosome damage -> micronuclei."""

    name = "cytogenetic (comet + CBMN)"
    state_names = ("adduct_bulky", "adduct_alkyl", "SSB", "DSB",
                   "acentric", "spindle",
                   "p53", "Mdm2_mRNA", "Mdm2_cyt", "Mdm2_nuc",
                   "cells_undivided", "cells_divided",
                   "MN_centromere_neg", "MN_centromere_pos")

    # -- channel routing ----------------------------------------------------
    # Each channel goes to the pool it physically belongs in.  Note that
    # `aneugenic` reaches neither break pool: an aneugen must be comet-
    # negative, and it is, structurally.
    to_bulky: dict = field(default_factory=lambda: {
        "bulky_adduct": 1.00, "icl": 0.85})
    to_alkyl: dict = field(default_factory=lambda: {
        "alkylation": 1.00})
    to_ssb: dict = field(default_factory=lambda: {
        "ssb": 1.00, "oxidative": 0.55})
    to_dsb: dict = field(default_factory=lambda: {
        "dsb": 1.00, "topo": 0.70})
    to_spindle: dict = field(default_factory=lambda: {
        "aneugenic": 1.00})

    # -- repair -------------------------------------------------------------
    k_ner: float = 0.020        # bulky adduct excision (slow)
    k_ber: float = 0.050        # alkyl adduct excision
    k_ligate: float = 0.100     # SSB rejoining (fast)
    k_rejoin: float = 0.010     # DSB rejoining
    f_misrepair: float = 0.06   # fraction of DSB rejoining that is wrong
    k_ssb_to_dsb: float = 0.004  # replication-dependent break conversion
    # Endogenous strand breaks: oxidative metabolism and ongoing repair keep
    # a small pool open at all times.  Without it the control comet reads
    # exactly zero and every induction ratio is a division by zero -- and a
    # real control slide reads a few percent tail, never none.
    ssb_spontaneous: float = 0.015
    k_clear_acentric: float = 0.0005
    k_clear_spindle: float = 0.50

    # -- DDR sensing --------------------------------------------------------
    K_dsb: float = 0.30
    K_ssb: float = 3.0
    K_adduct: float = 2.0
    K_spindle: float = 2.0
    w_dsb: float = 1.00
    w_ssb: float = 0.30
    w_adduct: float = 0.35
    w_spindle: float = 0.55

    # -- micronucleus formation, evaluated at each division -----------------
    # Background MN per division.  Part of the observed spontaneous rate now
    # arrives via the endogenous break pool above, so this is the remainder.
    p_spontaneous: float = 0.005
    p_max_clast: float = 0.45
    K_acentric: float = 3.0
    h_acentric: float = 2.0
    p_max_aneu: float = 0.55
    K_spindle_mn: float = 2.0
    h_spindle_mn: float = 2.0

    # -- shared p53 sub-model ------------------------------------------------
    p53: P53Mdm2 = field(default_factory=P53Mdm2)

    # -- growth / viability --------------------------------------------------
    mu_max: float = np.log(2) / 1440.0
    K_tox: float = 25.0
    h_tox: float = 3.0
    N0: float = 1.0

    def __post_init__(self):
        self._basal = self._solve_basal()

    # -- helpers -------------------------------------------------------------
    #: Per-lesion lethality weights (bulky, alkyl, SSB, DSB).  Counting an
    #: adduct as lethal as a double-strand break makes an adduct-former look
    #: acutely toxic at doses where it is only heavily adducted, which closes
    #: the comet's valid window before the assay can report anything.
    lethality: tuple = (0.30, 0.30, 0.50, 3.00)

    def burden(self, y) -> float:
        """Lethality-weighted lesion burden, for the shared viability term."""
        w = self.lethality
        return w[0] * y[0] + w[1] * y[1] + w[2] * y[2] + w[3] * y[3]

    def viability(self, load: float, extra_tox: float) -> float:
        v = 1.0 / (1.0 + (max(load, 0.0) / self.K_tox) ** self.h_tox)
        return v / (1.0 + max(extra_tox, 0.0))

    def ddr_activity(self, Adb, Ada, SSB, DSB, Sp, cycling) -> float:
        a = (self.w_dsb * _hill(DSB, self.K_dsb, 2.0)
             + self.w_ssb * _hill(SSB, self.K_ssb, 2.0)
             + self.w_adduct * _hill(Adb + Ada, self.K_adduct, 2.0)
             + self.w_spindle * _hill(Sp, self.K_spindle, 2.0) * cycling)
        return min(a, 1.0)

    def mn_probabilities(self, acentric, spindle, cycling) -> tuple:
        """Per-division probability of a centromere-neg / -pos micronucleus.

        The spontaneous background sits on the clastogenic side, which is
        where unstressed cells' micronuclei mostly come from, and it matters:
        without a non-zero control frequency there is no fold-induction to
        report.
        """
        p_c = (self.p_spontaneous
               + self.p_max_clast * _hill(acentric, self.K_acentric,
                                          self.h_acentric))
        p_a = self.p_max_aneu * _hill(spindle, self.K_spindle_mn,
                                      self.h_spindle_mn) * cycling
        return p_c, p_a

    def _solve_basal(self) -> np.ndarray:
        from scipy.integrate import solve_ivp
        exp = Exposure(flux=_ZERO_FLUX, duration_min=20000.0)
        y = np.array([0.0] * 6 + list(self.p53.seed) + [self.N0, 0.0, 0.0, 0.0])
        sol = solve_ivp(self.rhs, (0.0, 20000.0), y, args=(exp,),
                        method="LSODA", rtol=1e-9, atol=1e-11)
        y = sol.y[:, -1].copy()
        y[:6] = 0.0                      # damage pools start clean
        y[10:] = [self.N0, 0.0, 0.0, 0.0]  # and so does the cell bookkeeping
        return y

    # -- model ----------------------------------------------------------------
    def basal_state(self) -> np.ndarray:
        return self._basal.copy()

    def rhs(self, t, y, exp: Exposure):
        Adb, Ada, SSB, DSB, Fa, Sp = (max(v, 0.0) for v in y[:6])
        P, Mr, Mc, Mn = y[6:10]
        Nu, Nd, MNc, MNp = y[10:14]

        v = self.viability(self.burden(y), exp.extra_toxicity)
        cycling = v / (1.0 + max(exp.growth_inhibition, 0.0))
        mu = self.mu_max * cycling * self.p53.arrest(max(P, 0.0))
        # replication-associated break conversion needs cells in S phase; use
        # the achieved cycling rate as the proxy so an arrested culture stops
        # converting single-strand damage into double-strand breaks
        rep = mu / self.mu_max

        f = exp.flux
        # Excision repair removes adducts and, in doing so, *creates* the
        # transient single-strand gaps the alkaline comet detects.
        ner = self.k_ner * Adb
        ber = self.k_ber * Ada

        dAdb = f.project(self.to_bulky) - ner - mu * Adb
        dAda = f.project(self.to_alkyl) - ber - mu * Ada
        dSSB = (f.project(self.to_ssb) + self.ssb_spontaneous + ner + ber
                - self.k_ligate * SSB
                - self.k_ssb_to_dsb * SSB * rep - mu * SSB)
        dDSB = (f.project(self.to_dsb) + self.k_ssb_to_dsb * SSB * rep
                - self.k_rejoin * DSB - mu * DSB)
        dFa = (self.f_misrepair * self.k_rejoin * DSB
               - self.k_clear_acentric * Fa - mu * Fa)
        dSp = (f.project(self.to_spindle) - self.k_clear_spindle * Sp
               - mu * Sp)

        A = self.ddr_activity(Adb, Ada, SSB, DSB, Sp, cycling)
        dP, dMr, dMc, dMn = self.p53.rhs((P, Mr, Mc, Mn), A)

        # Cytochalasin-B block: cells divide their nucleus once, then stop.
        p_c, p_a = self.mn_probabilities(Fa, Sp, cycling)
        div = mu * Nu
        dNu, dNd = -div, div
        dMNc, dMNp = div * p_c, div * p_a

        return np.array([dAdb, dAda, dSSB, dDSB, dFa, dSp,
                         dP, dMr, dMc, dMn,
                         dNu, dNd, dMNc, dMNp])

    def observables(self, t, Y, exp):
        Adb, Ada, SSB, DSB, Fa, Sp = Y[:6]
        P = Y[6]
        Nu, Nd, MNc, MNp = Y[10:14]
        load = np.array([self.burden(Y[:, i]) for i in range(Y.shape[1])])
        v = np.array([self.viability(l, exp.extra_toxicity) for l in load])
        return {
            "t": t,
            "adduct_bulky": Adb,
            # alkyl adducts are themselves alkali-labile, so they migrate
            "alkali_labile": Ada,
            "ssb": SSB,
            "dsb": DSB,
            "acentric": Fa,
            "spindle": Sp,
            "p53": P,
            # lesion-burden name kept so damage_probe works on this core too
            "lesions": load,
            "cells_undivided": Nu,
            "cells_divided": Nd,
            "mn_centromere_neg": MNc,
            "mn_centromere_pos": MNp,
            "density": Nu + Nd,
            "viability": v,
        }
