"""
bioif.selftest -- tests for the interface guarantees.

Run: python3 -m bioif.selftest

These assert the *properties the layer exists to provide*, not the values of
the stub models: a type mismatch is caught, a missing context refuses, an
out-of-domain question widens, routing finds and prefers the cheaper path,
unit algebra round-trips exactly, support is enforced, and assumption flags
survive to the endpoint.
"""
from __future__ import annotations

import random
import traceback

from .adapters_demo import DEMO_ADAPTERS, RnaToProtein
from .chain import run_chain
from .core import (Claim, Context, Entity, Estimate, Provenance,
                   check_identity, EXTRAPOLATE, IN_DOMAIN)
from .demo import CALIBRATED, UNCALIBRATED, chemical_source, genetic_source
from .quantities import (DELTA_PSI, FITNESS, IC50_NM, PIC50, PKD, PKI,
                         PROTEIN_LFC, RESIDUAL_ACTIVITY, RNA_LFC)
from .registry import Registry

_results: list[tuple[str, bool, str]] = []


def check(name):
    def deco(fn):
        try:
            fn()
            _results.append((name, True, ""))
        except AssertionError as e:
            _results.append((name, False, str(e) or "assertion failed"))
        except Exception:
            _results.append((name, False, traceback.format_exc(limit=2)))
        return fn
    return deco


def _reg() -> Registry:
    r = Registry()
    for a in DEMO_ADAPTERS:
        r.register(a)
    return r


@check("quantity mismatch is a type error, not a silent join")
def t_type():
    a = RnaToProtein()
    bad = Claim(Entity("uniprot", "P04637", "isoform-1"), PROTEIN_LFC,
                CALIBRATED, Estimate.point(0.0, 8))
    try:
        a.apply(bad, random.Random(0))
    except TypeError:
        return
    raise AssertionError("adapter accepted a claim of the wrong quantity")


@check("same name+unit but different reference state does not type-check")
def t_reference():
    from .core import Quantity
    other = Quantity(RNA_LFC.name, RNA_LFC.unit, RNA_LFC.scale,
                     reference="DMSO at 6h")
    assert other.key != RNA_LFC.key, "reference state is not part of the type"


@check("missing required context refuses the chain")
def t_refuse():
    reg = _reg()
    path = reg.route(PIC50, FITNESS)
    res = run_chain(path, chemical_source(CALIBRATED), seed=1)
    assert not res.ok, "chain answered without a dose"
    assert "dose_uM" in res.refusal_reason, res.refusal_reason


@check("out-of-calibration system widens instead of answering quietly")
def t_extrapolate():
    reg = _reg()
    path = reg.route(DELTA_PSI, FITNESS)
    a = run_chain(path, genetic_source(CALIBRATED), seed=1)
    b = run_chain(path, genetic_source(UNCALIBRATED), seed=1)
    assert a.ok and b.ok
    assert b.traces[-1].verdict.status == EXTRAPOLATE
    assert a.traces[-1].verdict.status == IN_DOMAIN
    wa = a.final.estimate.ci()[1] - a.final.estimate.ci()[0]
    wb = b.final.estimate.ci()[1] - b.final.estimate.ci()[0]
    assert wb > wa * 1.3, f"interval did not widen: {wa:.3f} -> {wb:.3f}"


@check("routing finds both arms to the shared endpoint")
def t_route():
    reg = _reg()
    g = reg.route(DELTA_PSI, FITNESS)
    c = reg.route(PIC50, FITNESS)
    assert g and c, "no route found"
    assert [x.name for x in g][-1] == "residual_activity->fitness"
    assert [x.name for x in c][-1] == "residual_activity->fitness"
    # the chemical arm's two unit hops were never written by hand
    assert any(x.kind == "coercion" for x in c)


@check("among equal-length routes, fewer assumption-bearing bridges wins")
def t_route_cost():
    from .adapter import Adapter, BRIDGE

    class BridgeToFitness(Adapter):
        """Same endpoints as ActivityToFitness, one hop, but a bridge."""
        name, version, kind = "residual->fitness (bridge)", "0", BRIDGE
        consumes, produces = RESIDUAL_ACTIVITY, FITNESS
        assumption = "a large leap offered as a shortcut"

        def _forward(self, xs, claim, rng, noise):
            return [-x for x in xs]

    reg = _reg()
    reg.register(BridgeToFitness())
    path = reg.route(RESIDUAL_ACTIVITY, FITNESS)
    assert len(path) == 1, path
    assert path[0].kind != BRIDGE, \
        f"router took the bridge {path[0].name!r} over an equal-length fit"
    # and the full chain still ends on the empirical edge
    assert reg.route(DELTA_PSI, FITNESS)[-1].kind != BRIDGE


@check("unit coercions round-trip exactly and add no variance")
def t_units():
    reg = Registry()
    rng = random.Random(0)
    src = Claim(Entity("inchikey", "RYYVLZVUVIJVGH-UHFFFAOYSA-N"), PIC50,
                CALIBRATED, Estimate.normal(7.4, 0.5, 500, rng))
    path = reg.route(PIC50, IC50_NM) + reg.route(IC50_NM, PIC50)
    out = src
    for a in path:
        out = a.apply(out, rng)
    for x, y in zip(src.estimate.samples, out.estimate.samples):
        assert abs(x - y) < 1e-9, f"round trip drifted: {x} -> {y}"


@check("support bounds are enforced on every claim")
def t_support():
    reg = _reg()
    path = reg.route(DELTA_PSI, FITNESS)
    res = run_chain(path, genetic_source(CALIBRATED), seed=3)
    for t in res.traces:
        lo, hi = t.out.quantity.support
        assert all(lo - 1e-9 <= v <= hi + 1e-9 for v in t.out.estimate.samples)


@check("assumption flags survive to the endpoint")
def t_flags():
    reg = _reg()
    path = reg.route(DELTA_PSI, FITNESS)
    res = run_chain(path, genetic_source(CALIBRATED), seed=3)
    assert len(res.final.flags) == 3, res.final.flags
    assert any("rna->protein" in f for f in res.final.flags)


@check("provenance records every hop")
def t_provenance():
    reg = _reg()
    path = reg.route(DELTA_PSI, FITNESS)
    res = run_chain(path, genetic_source(CALIBRATED), seed=3)
    assert len(res.final.provenance) == len(path) + 1


@check("variance attribution finds the weakest link")
def t_attribution():
    reg = _reg()
    path = reg.route(DELTA_PSI, FITNESS)
    res = run_chain(path, genetic_source(CALIBRATED), seed=3)
    shares = {t.adapter.name: t.var_share for t in res.traces}
    assert shares["rna->protein"] > shares["protein->residual_activity"], shares
    assert shares["rna->protein"] > 0.2, shares


@check("both arms converge on the same entity and quantity")
def t_convergence():
    reg = _reg()
    g = run_chain(reg.route(DELTA_PSI, FITNESS),
                  genetic_source(CALIBRATED), seed=3)
    c = run_chain(reg.route(PIC50, FITNESS),
                  chemical_source(CALIBRATED.with_(dose_uM=0.3)), seed=3)
    assert g.ok and c.ok
    gn = [t.out for t in g.traces
          if t.out.quantity.key == RESIDUAL_ACTIVITY.key][0]
    cn = [t.out for t in c.traces
          if t.out.quantity.key == RESIDUAL_ACTIVITY.key][0]
    assert gn.entity == cn.entity, (gn.entity, cn.entity)
    # and the perturbagen that differed is in context, not in the identifier
    assert dict(gn.context.covariates)["perturbagen"] != \
        dict(cn.context.covariates)["perturbagen"]


@check("identifier lint rejects a build-less accession")
def t_lint():
    assert check_identity(Entity("ensembl_gene", "ENSG00000141510"))
    assert not check_identity(
        Entity("ensembl_gene", "ENSG00000141510", "GRCh38.p14/r112"))
    assert check_identity(Entity("inchikey", "TOO-SHORT"))


# ==========================================================================
# Tests that only exist because the affinity source is now REAL data.
# These run off the committed ChEMBL snapshot, so they are offline and
# deterministic; set BIOIF_REFRESH=1 to re-fetch.
# ==========================================================================

@check("REAL: a gene symbol does not resolve to a target")
def t_resolve_refuses():
    from .real import chembl
    amb = chembl.resolve_target("KRAS", organism="", target_type="")
    assert not amb.ok, "resolver guessed instead of refusing"
    assert len(amb.candidates) > 1
    assert len({h.target_type for h in amb.candidates}) > 1, \
        "expected candidates of several target types"
    ok = chembl.resolve_target("KRAS")
    assert ok.ok and ok.resolved.target_type == "SINGLE PROTEIN"


@check("REAL: measured records become one claim per assay, not one number")
def t_real_source_is_a_set():
    from .real import sources
    cl = sources.measured_affinity_claims()
    assert len(cl) > 5, f"only {len(cl)} records"
    assert len(sources.group_by_assay(cl)) > 1
    assert all(c.evidence == "measured" for c in cl)
    assert all(c.context.assay for c in cl), "a claim lost its assay"


@check("REAL: pIC50 and pKd do not silently become the same quantity")
def t_real_types_are_distinct():
    from .real import sources
    cl = sources.measured_affinity_claims()
    keys = {c.quantity.key for c in cl}
    assert len(keys) > 1, "all records collapsed to one quantity"
    reg = _reg()
    assert reg.route(PKD, FITNESS) is None, \
        "a pKd found a route into a pIC50 chain"
    assert reg.route(PIC50, FITNESS) is not None


@check("REAL: Cheng-Prusoff refuses without mechanism and [S]/Km")
def t_cheng_prusoff():
    from .real import sources
    cl = [c for c in sources.measured_affinity_claims()
          if c.quantity.key == PKD.key]
    assert cl, "no pKd records in the snapshot"
    k = cl[0].derive(quantity=PKI)
    a = sources.KiToIC50()
    assert a.domain(k).status == "refuse"
    supplied = k.derive(context=k.context
                        .with_cov("substrate_over_km", "10")
                        .with_cov("mechanism", "competitive"))
    assert a.domain(supplied).status == "in_domain"
    out = a.apply(supplied, random.Random(0))
    # pIC50 = pKi - log10(1 + [S]/Km); [S]/Km = 10 -> shift of log10(11)
    import math
    assert abs((k.estimate.mean - out.estimate.mean) - math.log10(11)) < 1e-9


@check("REAL: evidence degrades from measured and never recovers")
def t_evidence_monotone():
    from .core import EVIDENCE_ORDER
    from .real import sources
    from .demo_real import build_registry
    reg = build_registry()
    path = reg.route(PIC50, FITNESS)
    src = [c for c in sources.measured_affinity_claims()
           if c.quantity.key == PIC50.key][0]
    src = src.derive(context=Context(
        assay=src.context.assay, covariates=src.context.covariates,
        dose_uM=0.3, system="depmap:ACH-000019", time_h=120.0))
    res = run_chain(path, src, seed=1)
    assert res.ok
    levels = [src.evidence] + [t.out.evidence for t in res.traces]
    idx = [EVIDENCE_ORDER.index(x) for x in levels]
    assert idx == sorted(idx), f"evidence improved along the chain: {levels}"
    assert levels[0] == "measured"
    assert levels[-1] == "mechanistic_hypothesis"


@check("REAL: pooling across assays is flagged and demotes the evidence")
def t_pooling_is_visible():
    from .real import sources
    cl = sources.measured_affinity_claims()
    pooled = sources.naive_pooled_claim(cl)
    assert pooled.flags and "pooling" in pooled.flags[0]
    assert pooled.evidence == "inferred_association", pooled.evidence
    assert "POOLED" in pooled.context.assay


@check("REAL: assay identity moves potency more than readout type does")
def t_heterogeneity():
    from .real.heterogeneity import analyse
    r = analyse()
    assert r["n_records"] > 1000
    within = r["within_type_across_assays"]["median"]
    across = r["across_readout_types"]["median"]
    assert within is not None and across is not None
    assert within >= across, (
        f"expected within-type/across-assay spread ({within:.2f}) to be at "
        f"least the across-type spread ({across:.2f})")
    assert within >= 0.5, f"spread collapsed to {within:.2f}; snapshot changed?"


def main() -> int:
    width = max(len(n) for n, _, _ in _results)
    failed = 0
    for name, ok, msg in _results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}")
        if not ok:
            failed += 1
            print("        " + msg.replace("\n", "\n        "))
    print(f"\n{len(_results) - failed}/{len(_results)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
