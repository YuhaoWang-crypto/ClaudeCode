"""
bioif.adapters_demo -- worked adapters for the demo chain.

⚠️  EVERY NUMERIC CONSTANT IN THIS FILE IS ILLUSTRATIVE.

These are stubs with the *shape* of the real maps (direction, saturation,
which hops are noisy relative to which) so that the plumbing, the refusals
and the variance attribution can be exercised end to end. They are not
fitted, they are not validated, and no output of this module is a prediction
about any real gene, compound or cell line. Swapping any one of them for a
real model is a `Registry.register` call plus a domain predicate -- which is
the point being demonstrated.

The demo wires two arms that must meet:

    genetic arm    variant -> dPSI -> RNA lfc -> protein lfc -> residual activity
    chemical arm   compound -> pIC50 -> IC50 -> occupancy -> residual activity
                                                            \
                                                             -> fitness effect

Forcing both into `residual_target_activity` is the only reason a model
trained on genetic perturbation can score a compound at all. It is also
where the largest unpoliced assumption in this kind of pipeline lives, so
that bridge is flagged loudly rather than hidden in a join.
"""
from __future__ import annotations

import math
from typing import Sequence

from .adapter import Adapter, BRIDGE, EMPIRICAL
from .core import Context, Entity, Verdict, EXTRAPOLATE, IN_DOMAIN
from .quantities import (DELTA_PSI, FITNESS, OCCUPANCY, PROTEIN_LFC,
                         RESIDUAL_ACTIVITY, RNA_LFC)

# A stand-in for the identifier resolution service a real deployment needs.
# In production this is the single most load-bearing piece of infrastructure
# in the whole system and it is not a dict.
_VARIANT_TO_GENE = {"NC_000017.11:g.7676154G>A": ("ENSG00000141510", "P04637")}


def _clamp(x, lo, hi):
    return min(max(x, lo), hi)


class PsiToTranscriptDosage(Adapter):
    """dPSI -> functional transcript dosage, log2 vs isogenic WT."""
    name = "psi->transcript_dosage"
    version = "0.1-stub"
    kind = BRIDGE
    consumes = DELTA_PSI
    produces = RNA_LFC
    assumption = ("exon loss is assumed frame-disrupting and NMD-cleared, "
                  "and the variant heterozygous; neither is checked")

    def map_entity(self, e: Entity) -> Entity:
        gene, _ = _VARIANT_TO_GENE.get(e.id, ("ENSG00000000000", ""))
        return Entity("ensembl_gene", gene, build="GRCh38.p14/r112")

    def map_context(self, c: Context, e: Entity) -> Context:
        return c.with_cov("perturbagen", e.key)

    def _forward(self, xs: Sequence[float], claim, rng, noise):
        out = []
        for x in xs:
            # NMD efficiency is itself uncertain and transcript-specific.
            nmd = _clamp(rng.gauss(0.75, 0.20) if noise else 0.75, 0.0, 1.0)
            lost = 0.5 * max(0.0, -x) * nmd          # heterozygous
            out.append(math.log2(max(1.0 - lost, 1e-3)))
        return out


class RnaToProtein(Adapter):
    """
    Transcript lfc -> protein lfc.

    This is the hop that people wire with an `=` sign and that the variance
    attribution will show owning most of the final uncertainty. Post-
    transcriptional buffering, autoregulation and complex stoichiometry all
    live in this residual.
    """
    name = "rna->protein"
    version = "0.1-stub"
    kind = BRIDGE
    consumes = RNA_LFC
    produces = PROTEIN_LFC
    assumption = ("transcript->protein coupling taken as linear with a "
                  "single slope; buffering and complex stoichiometry unmodelled")

    def map_entity(self, e: Entity) -> Entity:
        uni = next((u for g, u in _VARIANT_TO_GENE.values()
                    if g == e.id), "P00000")
        return Entity("uniprot", uni, build="isoform-1")

    def _forward(self, xs, claim, rng, noise):
        slope = 0.60
        return [slope * x + (rng.gauss(0.0, 0.55) if noise else 0.0)
                for x in xs]


class ProteinToResidualActivity(Adapter):
    """Protein lfc -> residual pathway-level activity."""
    name = "protein->residual_activity"
    version = "0.1-stub"
    kind = BRIDGE
    consumes = PROTEIN_LFC
    produces = RESIDUAL_ACTIVITY
    assumption = ("activity taken proportional to abundance: no spare "
                  "capacity, no paralogue compensation, no scaffold role")

    def _forward(self, xs, claim, rng, noise):
        return [_clamp(2.0 ** x * (rng.gauss(1.0, 0.10) if noise else 1.0),
                       0.0, 1.0) for x in xs]


class OccupancyToResidualActivity(Adapter):
    """
    Fractional occupancy -> fractional functional inhibition.

    Same target quantity as the genetic arm, by construction. Whether the two
    are actually interchangeable is the scientific question the flag names;
    the interface's job is to make sure nobody reads the converged number
    without seeing it.
    """
    name = "occupancy->residual_activity"
    version = "0.1-stub"
    kind = BRIDGE
    consumes = OCCUPANCY
    produces = RESIDUAL_ACTIVITY
    assumption = ("equilibrium occupancy equated with functional inhibition, "
                  "and pharmacological inhibition equated with loss of gene "
                  "product; residence time and scaffolding ignored")
    #: The target this compound is being scored against. At a convergence
    #: node the two arms must agree on the ENTITY; what differs between them
    #: (a variant here, a compound at a dose there) moves into context.
    target = Entity("uniprot", "P04637", build="isoform-1")

    def map_entity(self, e: Entity) -> Entity:
        return self.target

    def map_context(self, c: Context, e: Entity) -> Context:
        return c.with_cov("perturbagen", e.key)

    def _forward(self, xs, claim, rng, noise):
        return [_clamp(1.0 - x * (_clamp(rng.gauss(0.85, 0.12), 0.0, 1.0)
                                  if noise else 0.85), 0.0, 1.0) for x in xs]


class ActivityToFitness(Adapter):
    """
    Residual activity -> fitness effect, in a named cell system.

    The only adapter here with a declared calibration set, so it is the one
    that demonstrates the context check: asked about a system it was not
    fitted on, it extrapolates loudly instead of answering quietly.
    """
    name = "residual_activity->fitness"
    version = "0.1-stub"
    kind = EMPIRICAL
    consumes = RESIDUAL_ACTIVITY
    produces = FITNESS
    requires_context = ("system",)
    calibrated_systems = ("depmap:ACH-000019", "depmap:ACH-000552")

    def _forward(self, xs, claim, rng, noise):
        amp = 1.20                       # dependency strength of this target
        return [-amp * (1.0 - x) ** 1.5 +
                (rng.gauss(0.0, 0.15) if noise else 0.0) for x in xs]


DEMO_ADAPTERS = [PsiToTranscriptDosage(), RnaToProtein(),
                 ProteinToResidualActivity(), OccupancyToResidualActivity(),
                 ActivityToFitness()]
