"""
Upstream layer — the DNA-damage state vector and the sources that produce it.

This is the seam between chemistry and biology.  Everything upstream (QSAR,
quantum descriptors, docking, metabolic activation) has exactly one job: emit
a :class:`DamageFlux`.  Everything downstream (SOS, p53/DDR, comet, MN) reads
that same object.  No assay-specific information crosses this boundary.

Why a vector and not a scalar
-----------------------------
A single "genotoxic potency" number cannot drive more than one assay
correctly, because different endpoints read different lesion channels:

  * umu / SOS          <- ssDNA at stalled forks (bulky, ICL, SSB, some oxid.)
  * comet (alkaline)   <- SSB + ALS + DSB, directly
  * micronucleus       <- DSB (clastogenic route) OR ``aneugenic`` (spindle),
                          two mechanisms that do not share a common lesion
  * chromosome aberr.  <- DSB misrepair

A spindle poison has ``aneugenic > 0`` and every DNA channel at zero.  With a
scalar input, a model is forced to call it SOS-positive; with this vector the
SOS core weights ``aneugenic`` at 0 and correctly returns negative.  That
behaviour is exercised in :mod:`genotox.run_umu`.

Rigour label
------------
The *structure* of this layer is mechanistic.  The *numbers* attached to the
demo compounds in :data:`DEMO_COMPOUNDS` are order-of-magnitude illustrative
values chosen to put the responses in the usual experimental range.  They are
NOT fitted to any measured dataset and must not be read as predictions for
those chemicals.  Marked (H) = hypothesis throughout.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields, replace

# Lesion channels, in a fixed order.  Extend here, not in the cores.
CHANNELS = (
    "bulky_adduct",   # NER substrate; replication-blocking
    "alkylation",     # N7-/O6-alkylG
    "oxidative",      # 8-oxoG and friends
    "ssb",            # single-strand breaks / alkali-labile sites
    "dsb",            # double-strand breaks
    "icl",            # interstrand crosslinks
    "topo",           # trapped topoisomerase cleavage complexes
    "aneugenic",      # spindle / kinetochore interference — NOT a DNA lesion
)


@dataclass(frozen=True)
class DamageFlux:
    """Lesion production rate per cell, per minute, by channel.

    Units are "lesion-equivalents per cell per minute".  They are internally
    consistent across channels only in the weak sense that the downstream
    cores apply their own per-channel efficiencies; do not compare channels
    to each other as if they were on one absolute scale.
    """

    bulky_adduct: float = 0.0
    alkylation: float = 0.0
    oxidative: float = 0.0
    ssb: float = 0.0
    dsb: float = 0.0
    icl: float = 0.0
    topo: float = 0.0
    aneugenic: float = 0.0

    def scaled(self, k: float) -> "DamageFlux":
        return replace(self, **{f.name: getattr(self, f.name) * k
                                for f in fields(self)})

    def project(self, weights: dict) -> float:
        """Collapse the vector onto one scalar using per-channel weights.

        Each downstream core owns its own weights, which is precisely how the
        same upstream chemistry yields different answers per endpoint.
        """
        return sum(getattr(self, c) * weights.get(c, 0.0) for c in CHANNELS)

    def as_dict(self) -> dict:
        return {c: getattr(self, c) for c in CHANNELS}

    def nonzero(self) -> dict:
        return {c: v for c, v in self.as_dict().items() if v}


@dataclass
class Compound:
    """A test article as the upstream layer sees it.

    ``per_uM`` is the lesion flux produced by 1 uM of the *proximate* (already
    activated, if activation is needed) form.  ``activation`` says how much of
    the parent is direct-acting versus S9-dependent.
    """

    name: str
    per_uM: DamageFlux
    direct_fraction: float = 1.0     # active without metabolic activation
    s9_fraction: float = 0.0         # additional activity unlocked by S9
    # Two distinct, separately measurable liabilities.  Collapsing them into
    # one multiplier makes every strong cytotoxicant look non-inducing,
    # because shutting down translation also shuts down the reporter.
    cytotoxic_per_uM: float = 0.0    # lethal: blocks biosynthesis AND growth
    cytostatic_per_uM: float = 0.0   # bacteriostatic: blocks growth only
    note: str = ""

    def flux(self, dose_uM: float, s9: bool = False) -> DamageFlux:
        eff = self.direct_fraction + (self.s9_fraction if s9 else 0.0)
        return self.per_uM.scaled(dose_uM * eff)

    def extra_toxicity(self, dose_uM: float) -> float:
        return self.cytotoxic_per_uM * dose_uM

    def growth_inhibition(self, dose_uM: float) -> float:
        return self.cytostatic_per_uM * dose_uM


class DamageSource:
    """Base class for anything that turns a compound + dose into a flux.

    Subclass this to plug a real upstream in.  The contract is one method.
    A QSAR/descriptor-based source would implement :meth:`flux` by predicting
    per-channel reactivity from structure instead of reading a table.
    """

    #: free-text provenance, printed in reports so a run is never ambiguous
    provenance = "abstract"

    def flux(self, compound: Compound, dose_uM: float, s9: bool = False) -> DamageFlux:
        raise NotImplementedError

    def toxicity(self, compound: Compound, dose_uM: float) -> float:
        return compound.extra_toxicity(dose_uM)

    def growth_inhibition(self, compound: Compound, dose_uM: float) -> float:
        return compound.growth_inhibition(dose_uM)


class TabulatedSource(DamageSource):
    """Reference source: reads the per-channel rates off the Compound itself.

    Deliberately trivial.  Its purpose is to let the downstream machinery be
    developed and tested against a *known* input, so that when a predictive
    source replaces it the change in behaviour is attributable.
    """

    provenance = "tabulated, illustrative (H) — not fitted to measured data"

    def flux(self, compound, dose_uM, s9=False):
        return compound.flux(dose_uM, s9=s9)


class QSARSource(DamageSource):
    """Placeholder for the predictive upstream.

    Left unimplemented on purpose.  The honest blocker is that the public
    labelled data for most channels is the *assay outcome*, not the lesion
    flux — training on Ames labels and then feeding an Ames model is circular.
    Channels that can be derived from first principles or from target-binding
    data (electrophilicity -> alkylation/bulky; tubulin binding -> aneugenic;
    Top2 pharmacophore -> topo) are the ones to implement first.
    """

    provenance = "not implemented"

    def flux(self, compound, dose_uM, s9=False):
        raise NotImplementedError(
            "QSARSource is a declared seam, not a working model. "
            "Implement per-channel prediction before using it; see docstring."
        )


# --------------------------------------------------------------------------
# Demo compounds.  ALL NUMBERS ILLUSTRATIVE (H).
# Chosen to span the four behaviours the architecture must distinguish.
# --------------------------------------------------------------------------
DEMO_COMPOUNDS = [
    Compound(
        name="direct-acting bulky (4NQO-like)",
        per_uM=DamageFlux(bulky_adduct=0.55, oxidative=0.20, ssb=0.05),
        direct_fraction=1.0, s9_fraction=0.0,
        cytotoxic_per_uM=0.0,
        note="strong SOS inducer without activation",
    ),
    Compound(
        name="promutagen (2AA-like)",
        per_uM=DamageFlux(bulky_adduct=0.50),
        direct_fraction=0.002, s9_fraction=0.85,
        note="near-silent without S9; the activation layer decides the call",
    ),
    Compound(
        name="alkylating agent (MMS-like)",
        per_uM=DamageFlux(alkylation=0.40, ssb=0.25),
        direct_fraction=1.0,
        note="SOS via BER intermediates -> gaps; weaker per lesion than bulky",
    ),
    Compound(
        name="aneugen (colchicine-like)",
        per_uM=DamageFlux(aneugenic=1.20),
        direct_fraction=1.0,
        # Spindle poisons are strongly antiproliferative -- they arrest cells
        # in mitosis.  Omitting that gives an aneugen with an unbounded valid
        # window, which no real one has.
        cytostatic_per_uM=0.08,
        note="NO DNA lesion; must be SOS-negative and MN-positive",
    ),
    Compound(
        name="masked genotoxicant",
        per_uM=DamageFlux(bulky_adduct=0.18),
        direct_fraction=1.0,
        cytostatic_per_uM=3.0,
        note="genotoxic, but kills the culture before the signal clears the "
             "threshold -> the assay cannot decide, and must say so",
    ),
    Compound(
        name="non-genotoxic cytotoxicant",
        per_uM=DamageFlux(),
        direct_fraction=1.0,
        cytotoxic_per_uM=0.030,
        note="kills growth without inducing; the assay's false-positive trap",
    ),
]
