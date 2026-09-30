"""
bioif.core -- the interchange contract.

The premise of this package is that heterogeneous biological models fail to
compose for reasons that are NOT file formats. A tensor of floats is easy to
move. What is hard is that every model emits a *different kind of number*,
conditioned on a *different context*, with a *different reference state*, and
with no statement of where it stops being valid.

So the contract carries six things, and a link between two models only
type-checks if all six line up (or an explicit, assumption-bearing bridge is
inserted):

  1. Entity     -- WHAT the number is about, in a resolvable namespace+build.
  2. Quantity   -- WHAT KIND of number it is: name, unit, scale, reference.
  3. Context    -- the conditioning: system, dose, time, assay.
  4. Estimate   -- a distribution, never a bare point.
  5. Domain     -- an explicit in/out-of-domain verdict, with refusal allowed.
  6. Provenance -- model, version, what it was calibrated on.

Everything else in bioif is machinery around these six fields.
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field, replace
from typing import Sequence

# --------------------------------------------------------------------------
# 1. Entity -- identity with a namespace and a build
# --------------------------------------------------------------------------

#: Namespaces we know how to reason about. The `build` field is not optional
#: decoration: an Ensembl gene id without an assembly, or a UniProt accession
#: without an isoform, is the single most common silent join error when wiring
#: a DNA-level model to a protein-level model.
NAMESPACES = {
    "hgvs_g":       "genomic HGVS, build required (GRCh38)",
    "ensembl_gene": "ENSG, assembly/annotation release required",
    "ensembl_tx":   "ENST, annotation release required",
    "uniprot":      "UniProt accession, isoform suffix required",
    "inchikey":     "InChIKey (27 char), tautomer/salt-stripped",
    "depmap":       "DepMap model id (ACH-######)",
    "cl":           "Cell Ontology term",
    "efo":          "EFO / disease term",
}


@dataclass(frozen=True)
class Entity:
    """A thing a claim is about."""
    ns: str
    id: str
    build: str = ""

    def __post_init__(self):
        if self.ns not in NAMESPACES:
            raise ValueError(f"unknown namespace {self.ns!r}")

    @property
    def key(self) -> str:
        return f"{self.ns}:{self.id}" + (f"@{self.build}" if self.build else "")

    def __str__(self) -> str:
        return self.key


#: Namespaces where omitting `build` is a correctness bug rather than a style
#: nit, because the identifier is not stable without it.
BUILD_REQUIRED = {"hgvs_g", "ensembl_gene", "ensembl_tx", "uniprot"}


# --------------------------------------------------------------------------
# 2. Quantity -- the semantic type of the number
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Quantity:
    """
    What kind of number this is.

    `reference` is the field that catches the errors nobody catches by eye:
    a log2 fold change "vs DMSO at 6 h" and a log2 fold change "vs isogenic
    WT" have the same name, unit and scale, and are not the same quantity.
    Two models that both emit "log2FC" are not composable until this agrees.
    """
    name: str
    unit: str                       # 'dimensionless' | 'nM' | 'fraction' | ...
    scale: str                      # 'linear' | 'log2' | 'log10' | 'logit'
    reference: str = ""             # the contrast; "" means absolute
    support: tuple[float, float] = (-math.inf, math.inf)

    @property
    def key(self) -> str:
        return f"{self.name}|{self.unit}|{self.scale}|{self.reference}"

    def __str__(self) -> str:
        r = f" vs {self.reference}" if self.reference else ""
        return f"{self.name} [{self.unit}, {self.scale}]{r}"


# --------------------------------------------------------------------------
# 3. Context -- the conditioning that makes or breaks transfer
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Context:
    """
    The conditions the number is conditioned on.

    Most cross-model composition failures in practice are not wrong maths,
    they are a context mismatch nobody declared: model A predicts in K562,
    model B was fitted in HepG2, and the pipeline joins them on gene id.
    """
    system: str = ""                # 'depmap:ACH-000001' | 'cl:0000236' | ''
    dose_uM: float | None = None
    time_h: float | None = None
    assay: str = ""
    covariates: tuple[tuple[str, str], ...] = ()

    def with_(self, **kw) -> "Context":
        return replace(self, **kw)

    def with_cov(self, key: str, value: str) -> "Context":
        """
        Record a conditioning fact that has no dedicated field.

        This is where a perturbagen goes once a chain converges: at a shared
        node, two arms must agree on the *entity*, and what differs between
        them -- a variant on one side, a compound at a dose on the other --
        becomes context, not identity.
        """
        return replace(self, covariates=self.covariates + ((key, value),))

    def describe(self) -> str:
        bits = [f"system={self.system or '-'}"]
        if self.dose_uM is not None:
            bits.append(f"dose={self.dose_uM:g}uM")
        if self.time_h is not None:
            bits.append(f"t={self.time_h:g}h")
        if self.assay:
            bits.append(f"assay={self.assay}")
        bits += [f"{k}={v}" for k, v in self.covariates]
        return " ".join(bits)


# --------------------------------------------------------------------------
# 4. Estimate -- a distribution, never a bare point
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Estimate:
    """
    A Monte Carlo sample set.

    Samples rather than (mu, sigma) because every interesting bridge in a
    biological chain is nonlinear -- Hill functions, thresholds, log/linear
    hops -- and moment propagation through those is wrong in exactly the
    regime you care about (near a threshold).
    """
    samples: tuple[float, ...]

    @staticmethod
    def point(value: float, n: int = 1) -> "Estimate":
        return Estimate(tuple([float(value)] * n))

    @staticmethod
    def normal(mu: float, sd: float, n: int, rng) -> "Estimate":
        return Estimate(tuple(rng.gauss(mu, sd) for _ in range(n)))

    @property
    def mean(self) -> float:
        return statistics.fmean(self.samples)

    @property
    def sd(self) -> float:
        return statistics.pstdev(self.samples) if len(self.samples) > 1 else 0.0

    @property
    def var(self) -> float:
        return self.sd ** 2

    def quantile(self, q: float) -> float:
        s = sorted(self.samples)
        if len(s) == 1:
            return s[0]
        pos = q * (len(s) - 1)
        lo = int(math.floor(pos))
        hi = min(lo + 1, len(s) - 1)
        return s[lo] + (s[hi] - s[lo]) * (pos - lo)

    def ci(self, level: float = 0.90) -> tuple[float, float]:
        a = (1.0 - level) / 2.0
        return self.quantile(a), self.quantile(1.0 - a)

    def summary(self, level: float = 0.90) -> str:
        lo, hi = self.ci(level)
        return f"{self.mean:+.3f} [{lo:+.3f}, {hi:+.3f}]"


# --------------------------------------------------------------------------
# 5. Domain -- the right to refuse
# --------------------------------------------------------------------------

#: An adapter that cannot say "I don't know" will silently extrapolate, and a
#: long chain of silent extrapolations is exactly how a pipeline produces a
#: confident number with no support. Three verdicts, not two.
IN_DOMAIN, EXTRAPOLATE, REFUSE = "in_domain", "extrapolate", "refuse"


@dataclass(frozen=True)
class Verdict:
    status: str
    reason: str = ""
    #: multiplicative inflation applied to this link's own noise when the
    #: status is EXTRAPOLATE -- an honest widening, not a correction.
    inflate: float = 1.0


# --------------------------------------------------------------------------
# 6. Provenance
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Provenance:
    model: str
    version: str
    calibrated_on: str = ""
    note: str = ""

    def __str__(self) -> str:
        return f"{self.model}@{self.version}"


# --------------------------------------------------------------------------
# The claim: the unit that moves between models
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Claim:
    entity: Entity
    quantity: Quantity
    context: Context
    estimate: Estimate
    provenance: tuple[Provenance, ...] = ()
    #: Free-text markers that survive the whole chain. An assumption made at
    #: link 2 must still be visible on the final number, or the consumer will
    #: read the output as if it were measured.
    flags: tuple[str, ...] = ()

    def derive(self, **kw) -> "Claim":
        return replace(self, **kw)

    def add_flag(self, f: str) -> "Claim":
        return replace(self, flags=self.flags + (f,))

    def clipped(self) -> "Claim":
        lo, hi = self.quantity.support
        s = tuple(min(max(v, lo), hi) for v in self.estimate.samples)
        return replace(self, estimate=Estimate(s))

    def describe(self) -> str:
        return (f"{self.entity}  {self.quantity}\n"
                f"    context: {self.context.describe()}\n"
                f"    value:   {self.estimate.summary()}")


def check_identity(entity: Entity) -> list[str]:
    """Static lint on an identifier: the cheapest bug class to kill."""
    problems = []
    if entity.ns in BUILD_REQUIRED and not entity.build:
        problems.append(
            f"{entity.ns} identifier {entity.id!r} has no build/release; "
            "joins against another model are not reproducible")
    if entity.ns == "inchikey" and len(entity.id) != 27:
        problems.append(f"InChIKey {entity.id!r} is not 27 characters")
    return problems
