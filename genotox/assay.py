"""
Composition layer — bolts an upstream source, a signal core and a set of
readouts into one runnable virtual assay.

    VirtualAssay(source=TabulatedSource(), core=SOSCore(),
                 readouts=["umu_betagal_ONPG", "growth"])

Swapping any one of the three is a constructor argument, which is the whole
point of the split: endpoint #2 changes ``core``, a predictive chemistry layer
changes ``source``, and a new instrument changes ``readouts``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .core import Exposure, SignalCore, SOSCore
from .damage import Compound, DamageSource, TabulatedSource
from .readouts import apply_readouts


@dataclass
class VirtualAssay:
    source: DamageSource = field(default_factory=TabulatedSource)
    core: SignalCore = field(default_factory=SOSCore)
    readouts: tuple = ("umu_betagal_ONPG", "growth", "damage_probe")
    duration_min: float = 120.0
    #: trajectory sampling.  Matters for cores whose signal oscillates: too
    #: coarse a grid and PulseProbe silently miscounts peaks.
    n_points: int = 241

    @property
    def description(self) -> str:
        return (f"{self.core.name} | source: {self.source.provenance} | "
                f"readouts: {', '.join(self.readouts)}")

    def well(self, compound: Compound, dose_uM: float, s9: bool = False) -> dict:
        """Simulate a single well and return its measurements."""
        exp = Exposure(
            flux=self.source.flux(compound, dose_uM, s9=s9),
            extra_toxicity=self.source.toxicity(compound, dose_uM),
            growth_inhibition=self.source.growth_inhibition(compound, dose_uM),
            s9=s9,
            duration_min=self.duration_min,
        )
        obs = self.core.simulate(exp, n_points=self.n_points)
        res = apply_readouts(obs, self.readouts)
        res.update(dose_uM=dose_uM, s9=s9, compound=compound.name)
        return res

    def control(self, compound: Compound, s9: bool = False) -> dict:
        """Solvent control: same conditions, zero dose."""
        return self.well(compound, 0.0, s9=s9)

    def trajectory(self, compound: Compound, dose_uM: float,
                   s9: bool = False) -> dict:
        """Full time course, for plotting and for debugging the core."""
        exp = Exposure(
            flux=self.source.flux(compound, dose_uM, s9=s9),
            extra_toxicity=self.source.toxicity(compound, dose_uM),
            growth_inhibition=self.source.growth_inhibition(compound, dose_uM),
            s9=s9,
            duration_min=self.duration_min,
        )
        return self.core.simulate(exp, n_points=self.n_points)
