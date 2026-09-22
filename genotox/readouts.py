"""
Readout layer — turns core trajectories into what an instrument reports.

Kept separate from the cores so that one simulated biology can be measured
several ways, which is the situation the real kits are in: the same DNA
damage is read as beta-gal absorbance (umu), tail moment (comet), or
micronucleus frequency, and the differences between those endpoints are
partly differences of *measurement*, not of biology.

A readout consumes the observables dict a core emits and returns a flat dict
of scalars.  It must declare which observables it needs, so a mismatch fails
loudly at wiring time instead of silently reading a zero.
"""
from __future__ import annotations

from dataclasses import dataclass


class Readout:
    name = "abstract"
    requires: tuple = ()

    def check(self, obs: dict) -> None:
        missing = [k for k in self.requires if k not in obs]
        if missing:
            raise KeyError(
                f"readout {self.name!r} needs observables {missing} which "
                f"core does not expose (has: {sorted(obs)})")

    def __call__(self, obs: dict) -> dict:
        raise NotImplementedError


@dataclass
class BetaGalONPG(Readout):
    """The umu test's actual measurement.

    beta-galactosidase hydrolyses ONPG to o-nitrophenol, quantified at
    410 nm.  Absorbance is proportional to *total* enzyme in the well, i.e.
    to (enzyme per cell) x (cell density), which is why the assay must be
    normalised by growth (A600) before any induction ratio is formed.  That
    normalisation is the step where a purely cytotoxic compound would
    otherwise masquerade as an inducer.

    Substrate depletion is not modelled (assumes ONPG saturating over the
    incubation); at very high enzyme this would overestimate A410.
    """

    name = "umu_betagal_ONPG"
    requires = ("reporter_enzyme", "density")
    onpg_min: float = 30.0       # ONPG incubation
    gain: float = 1.0            # lumps kcat, path length, extinction coeff.

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        B_end = float(obs["reporter_enzyme"][-1])
        N_end = float(obs["density"][-1])
        a410 = self.gain * B_end * N_end * self.onpg_min
        a600 = N_end
        # Miller-style specific activity: absorbance per unit biomass per min
        units = a410 / (a600 * self.onpg_min) if a600 > 0 else 0.0
        return {"A410": a410, "A600": a600, "units": units,
                "enzyme_per_cell": B_end}


@dataclass
class GrowthReadout(Readout):
    """Culture density, used as the assay's own validity gate.

    The umu protocol discards wells whose growth factor falls below ~0.5,
    because at that point the induction ratio is no longer interpretable.
    Reporting the gate alongside the signal is not optional bookkeeping — an
    unreported failed gate is how a cytotoxic compound gets called genotoxic.
    """

    name = "growth"
    requires = ("density", "viability")

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        return {"density_end": float(obs["density"][-1]),
                "viability_end": float(obs["viability"][-1])}


@dataclass
class DamageProbe(Readout):
    """Internal state, exposed as a pseudo-measurement.

    Not something the umu kit reports.  It exists because lesion burden is
    the quantity shared with comet (break density) and with gamma-H2AX foci,
    so it is the natural place to cross-calibrate endpoint #3 against this
    core rather than re-deriving damage from scratch.
    """

    name = "damage_probe"
    requires = ("lesions", "recA_active", "promoter")

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        return {"lesions_end": float(obs["lesions"][-1]),
                "lesions_max": float(obs["lesions"].max()),
                "recA_active_end": float(obs["recA_active"][-1]),
                "promoter_end": float(obs["promoter"][-1])}


#: Readouts available by name.  Endpoint #2/#3 register here.
REGISTRY = {r.name: r for r in (BetaGalONPG(), GrowthReadout(), DamageProbe())}


def apply_readouts(obs: dict, names=None) -> dict:
    names = list(REGISTRY) if names is None else list(names)
    out = {}
    for n in names:
        if n not in REGISTRY:
            raise KeyError(f"unknown readout {n!r}; have {sorted(REGISTRY)}")
        out.update(REGISTRY[n](obs))
    return out
