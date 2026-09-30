"""
bioif.real.sources -- real measured claims, and the bridge that will not fire.

`measured_affinity_claims` turns ChEMBL activity records into Claims. The
important design decision is that it returns a **list**, one claim per assay
record, rather than one pooled number per compound. Under the contract, two
measurements taken in different assays are two claims about different
contexts; pooling them is a modelling choice that has to be made explicitly
and flagged, which is what `naive_pooled_claim` does so the two can be
compared.

`KiToIC50` implements the conversion everyone reaches for. It is a bridge,
not a coercion, and it refuses: Cheng-Prusoff needs the substrate
concentration and Km, and an affinity record does not carry them. This is the
single clearest case of a hop that *looks* like arithmetic and is not.
"""
from __future__ import annotations

import math
from collections import defaultdict

from ..adapter import Adapter, BRIDGE
from ..core import (Claim, Context, Entity, Estimate, Provenance,
                    MEASURED, INFERRED_ASSOCIATION)
from ..quantities import BY_STANDARD_TYPE, PIC50, PKI
from . import chembl

#: The demo target. KRAS, resolved (not guessed) via chembl.resolve_target.
KRAS = Entity("uniprot", "P01116", build="isoform-1")


def _rows(target: str, molecule: str | None, std_type: str | None,
          snapshot: bool) -> list[dict]:
    rows = (chembl.load_snapshot() if snapshot
            else chembl.fetch_activities(target, molecule))
    out = [r for r in rows if r.get("target_chembl_id") == target]
    if molecule:
        out = [r for r in out if r.get("molecule_chembl_id") == molecule]
    if std_type:
        out = [r for r in out if r.get("standard_type") == std_type]
    return [r for r in out if r.get("pchembl_value") not in (None, "")]


def measured_affinity_claims(molecule: str = chembl.DEMO_MOLECULE,
                             target: str = chembl.DEMO_TARGET,
                             std_type: str | None = None,
                             entity: Entity = KRAS,
                             snapshot: bool = True,
                             n: int = 4000) -> list[Claim]:
    """
    One Claim per measured activity record.

    Each carries the assay it came from, the standard_type it actually is
    (pIC50 and pKi are different quantities), and evidence=measured.

    A ChEMBL record reports no measurement error, so the value is carried as
    an exact point replicated `n` times. Any interval seen downstream is
    therefore contributed entirely by the adapters, not by the measurement --
    which is the honest reading, and makes the per-assay spread in the demo
    attributable to assay identity alone.
    """
    claims = []
    for r in _rows(target, molecule, std_type, snapshot):
        q = BY_STANDARD_TYPE.get(r["standard_type"])
        if q is None:
            continue                      # unmodelled readout; do not coerce
        ctx = (Context(assay=r["assay_chembl_id"])
               .with_cov("perturbagen", f"chembl:{r['molecule_chembl_id']}")
               .with_cov("standard_type", r["standard_type"])
               .with_cov("relation", r.get("standard_relation") or "?")
               .with_cov("document", r.get("document_chembl_id") or "?"))
        claims.append(Claim(
            entity=entity, quantity=q, context=ctx,
            estimate=Estimate.point(float(r["pchembl_value"]), n),
            evidence=MEASURED,
            provenance=(Provenance("ChEMBL", "REST/34",
                                   r["assay_chembl_id"],
                                   (r.get("assay_description") or "")[:110]),),
        ))
    return claims


def naive_pooled_claim(claims: list[Claim], n: int = 6000) -> Claim:
    """
    What a pipeline does when it asks a database for "the affinity".

    Pools every record for the compound into one empirical distribution,
    across assay types and across assays. The contract cannot stop you doing
    this -- it can only make you do it in the open, so the flag and the
    demoted evidence level travel with the number.
    """
    if not claims:
        raise ValueError("no claims to pool")
    vals = [c.estimate.mean for c in claims]
    samples = tuple(vals[i % len(vals)] for i in range(n))
    types = sorted({dict(c.context.covariates)["standard_type"] for c in claims})
    assays = sorted({c.context.assay for c in claims})
    base = claims[0]
    return Claim(
        entity=base.entity, quantity=PIC50,
        context=Context(assay=f"POOLED({len(assays)} assays)")
        .with_cov("perturbagen", dict(base.context.covariates)["perturbagen"])
        .with_cov("pooled_types", "+".join(types)),
        estimate=Estimate(samples),
        # pooling across assays and readout types is no longer a measurement
        # of anything; it is an association over a heterogeneous set
        evidence=INFERRED_ASSOCIATION,
        provenance=tuple(p for c in claims for p in c.provenance)[:1],
        flags=(f"[pooling] {len(claims)} records from {len(assays)} assays and "
               f"readout types {'+'.join(types)} were pooled into one pIC50; "
               f"assay identity discarded",),
    )


class KiToIC50(Adapter):
    """
    pKi -> pIC50 by Cheng-Prusoff.

    Looks like unit algebra. Is not: IC50 = Ki * (1 + [S]/Km) for competitive
    inhibition, so the conversion needs the substrate concentration and Km of
    the assay being converted *into*, plus the assumption that the mechanism
    is competitive. An affinity record carries none of those, so this adapter
    refuses unless they are supplied explicitly.
    """
    name = "pKi->pIC50 (Cheng-Prusoff)"
    version = "0.1"
    kind = BRIDGE
    consumes = PKI
    produces = PIC50
    requires_covariates = ("substrate_over_km", "mechanism")
    assumption = ("Cheng-Prusoff assumes competitive inhibition and a known "
                  "[S]/Km for the target assay; neither is in an affinity record")

    def _forward(self, xs, claim, rng, noise):
        cov = dict(claim.context.covariates)
        if cov.get("mechanism") != "competitive":
            raise ValueError("Cheng-Prusoff applied to a non-competitive "
                             "mechanism")
        ratio = float(cov["substrate_over_km"])
        # pIC50 = pKi - log10(1 + [S]/Km)
        return [x - math.log10(1.0 + ratio) for x in xs]


def group_by_assay(claims: list[Claim]) -> dict[str, list[Claim]]:
    g: dict[str, list[Claim]] = defaultdict(list)
    for c in claims:
        g[c.context.assay].append(c)
    return dict(g)
