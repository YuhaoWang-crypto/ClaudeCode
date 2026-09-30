"""
bioif.adapter -- the edge type.

An Adapter is one hop from one Quantity to another. It must declare, in
machine-readable form:

  * what it consumes and produces (Quantity keys, so the reference state is
    part of the type),
  * how it maps the entity (a variant hop lands on a gene; a gene hop lands
    on a protein),
  * which context axes it *requires* (a Hill occupancy step is meaningless
    without a dose, so it refuses rather than defaulting to one),
  * where it is in domain,
  * and whether it is a measurement-like map or an assumption-bearing bridge.

The split in that last point is the one that matters for long chains. Some
hops are arithmetic (pIC50 -> IC50 in nM: lossless). Some hops are a
scientific assumption wearing a function's clothes (RNA fold change ->
protein fold change). A chain that does not distinguish them will report the
same confidence for both.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from .core import (Claim, Context, Entity, Estimate, Provenance, Quantity,
                   Verdict, IN_DOMAIN, EXTRAPOLATE, REFUSE,
                   CALIBRATED_PREDICTION, MEASURED, MECHANISTIC_HYPOTHESIS,
                   weakest)

#: Edge kinds.
#:   'coercion'  -- algebraic identity, no information added or lost
#:   'empirical' -- fitted on data, carries its own residual noise
#:   'bridge'    -- crosses a level of biology on an assumption; the assumption
#:                  is named and flagged onto every downstream claim
COERCION, EMPIRICAL, BRIDGE = "coercion", "empirical", "bridge"


class Adapter:
    """Base class. Subclasses implement `_forward` and optionally `domain`."""

    name: str = "adapter"
    version: str = "0"
    kind: str = EMPIRICAL
    consumes: Quantity
    produces: Quantity
    #: context axes that must be present on the incoming claim
    requires_context: tuple[str, ...] = ()
    #: context covariate keys that must be present (for conditioning facts
    #: that have no dedicated field, e.g. the [S]/Km a Cheng-Prusoff
    #: conversion needs and that an affinity record does not carry)
    requires_covariates: tuple[str, ...] = ()
    #: systems this edge was fitted on; empty tuple means "system-agnostic"
    calibrated_systems: tuple[str, ...] = ()
    #: named assumption, surfaced as a flag on all downstream claims
    assumption: str = ""

    @property
    def max_evidence(self) -> str:
        """
        The strongest evidence level this edge can emit, whatever it was fed.

        Coercions are algebra, so they preserve whatever came in. An empirical
        fit can at best yield a calibrated prediction. A bridge crosses a level
        of biology on an assumption, so its output is a mechanistic hypothesis
        even when every input to it was measured.
        """
        if self.kind == COERCION:
            return MEASURED                    # i.e. imposes no ceiling
        if self.kind == BRIDGE:
            return MECHANISTIC_HYPOTHESIS
        return CALIBRATED_PREDICTION

    # -- identity -----------------------------------------------------------
    def map_entity(self, e: Entity) -> Entity:
        """Default: the claim stays about the same entity."""
        return e

    def map_context(self, c: Context, e: Entity) -> Context:
        """Default: the conditioning is unchanged. Adapters that collapse two
        arms onto a shared entity override this to push the thing that
        differed -- the perturbagen -- into context."""
        return c

    # -- domain -------------------------------------------------------------
    def domain(self, claim: Claim) -> Verdict:
        """
        Default domain policy, which every subclass inherits for free:
        missing required context is a refusal, and a system outside the
        calibration set is an extrapolation (widened, not silently accepted).
        """
        for axis in self.requires_context:
            if getattr(claim.context, axis, None) in (None, ""):
                return Verdict(REFUSE,
                               f"{self.name} requires context.{axis}, which "
                               f"the incoming claim does not carry")
        have = dict(claim.context.covariates)
        for cov in self.requires_covariates:
            if cov not in have:
                return Verdict(REFUSE,
                               f"{self.name} requires the covariate {cov!r}, "
                               f"which the incoming claim does not carry")
        if self.calibrated_systems and claim.context.system:
            if claim.context.system not in self.calibrated_systems:
                return Verdict(
                    EXTRAPOLATE,
                    f"{self.name} was calibrated on "
                    f"{', '.join(self.calibrated_systems)}; asked about "
                    f"{claim.context.system}",
                    inflate=2.0)
        return Verdict(IN_DOMAIN)

    # -- the map ------------------------------------------------------------
    def _forward(self, xs: Sequence[float], claim: Claim, rng,
                 noise: bool) -> list[float]:
        raise NotImplementedError

    def apply(self, claim: Claim, rng, noise: bool = True,
              inflate: float = 1.0) -> Claim:
        if claim.quantity.key != self.consumes.key:
            raise TypeError(
                f"{self.name} consumes {self.consumes} but was handed "
                f"{claim.quantity}")
        ys = self._forward(claim.estimate.samples, claim, rng,
                           noise and inflate >= 0)
        if noise and inflate != 1.0:
            # widen this link's own contribution around its central tendency
            m = sum(ys) / len(ys)
            ys = [m + (y - m) * inflate for y in ys]
        out = Claim(
            entity=self.map_entity(claim.entity),
            quantity=self.produces,
            context=self.map_context(claim.context, claim.entity),
            estimate=Estimate(tuple(ys)),
            evidence=weakest(claim.evidence, self.max_evidence),
            provenance=claim.provenance + (
                Provenance(self.name, self.version,
                           ",".join(self.calibrated_systems) or "n/a"),),
            flags=claim.flags,
        ).clipped()
        if self.assumption:
            out = out.add_flag(f"[{self.kind}] {self.name}: {self.assumption}")
        return out


class Coercion(Adapter):
    """
    An algebraic re-expression. Lossless, noiseless, and therefore safe to
    insert automatically when routing a chain -- which is what makes model
    swapping cheap: a new model that reports Ki in nM instead of pIC50 needs
    no new plumbing.
    """
    kind = COERCION

    def __init__(self, name: str, consumes: Quantity, produces: Quantity,
                 fn: Callable[[float, Claim], float], version: str = "1",
                 requires_context: tuple[str, ...] = ()):
        self.name = name
        self.consumes = consumes
        self.produces = produces
        self.fn = fn
        self.version = version
        self.requires_context = requires_context

    def _forward(self, xs, claim, rng, noise):
        return [self.fn(x, claim) for x in xs]
