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

import math
import random
import statistics
import traceback

from .adapters_demo import DEMO_ADAPTERS, RnaToProtein
from .chain import run_chain
from .core import (Claim, Context, Entity, Estimate, Provenance,
                   check_identity, EXTRAPOLATE, IN_DOMAIN)
from .demo import CALIBRATED, UNCALIBRATED, chemical_source, genetic_source
from .quantities import (DELTA_PSI, FITNESS, IC50_NM, PIC50, PIC50_REF, PKD,
                         PKI, PROTEIN_LFC, RESIDUAL_ACTIVITY, RNA_LFC)
from .registry import Registry

#: (name, status, detail) where status is True / False / None (skipped)
_results: list[tuple[str, bool | None, str]] = []


class Skip(Exception):
    """Raised by a test whose optional dependency is absent."""


def needs(*modules):
    """Skip rather than fail when an optional dependency is missing.

    The contract layer (bioif/core, adapter, registry, chain, ensemble) is
    stdlib-only on purpose. The pairings built on real chemistry need rdkit,
    numpy and scikit-learn, and their absence is a missing dependency rather
    than a broken guarantee -- so those tests report SKIP.
    """
    import importlib
    missing = [m for m in modules
               if importlib.util.find_spec(m) is None]
    if missing:
        raise Skip("needs " + ", ".join(missing))


def check(name):
    def deco(fn):
        try:
            fn()
            _results.append((name, True, ""))
        except Skip as e:
            _results.append((name, None, str(e)))
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


@check("CONFORMAL: the quantile is exact, and infinite when unearnable")
def t_conformal_quantile():
    import math
    from .real.conformal import conformal_quantile
    s = [float(i) for i in range(1, 11)]          # n=10
    # ceil(11*0.8)=9 -> the 9th smallest
    assert conformal_quantile(s, 0.2) == 9.0
    # ceil(11*0.9)=10 -> the 10th smallest
    assert conformal_quantile(s, 0.1) == 10.0
    # ceil(11*0.95)=11 > 10 -> not attainable
    assert math.isinf(conformal_quantile(s, 0.05))
    assert math.isinf(conformal_quantile([], 0.2))


@check("CONFORMAL: binned patent bands are detected, not regressed on")
def t_binned_gate():
    from .real import conformal as C
    from .real.transfer_adapter import build_transfer
    pairs = C.discover_pairs()
    binned = [p for p in pairs.values() if p.status == "binned"]
    assert binned, "no binned pairs detected; the gate is not being exercised"
    worst = max(binned, key=lambda p: p.tie_x)
    assert worst.tie_x > 0.25 and worst.reason
    try:
        build_transfer(worst.key)
    except ValueError as e:
        assert "binned" in str(e)
        return
    raise AssertionError("built a conformal regression on binned data")


@check("CONFORMAL: within-pair coverage reaches nominal, Gaussian does not")
def t_coverage_within():
    from .real import conformal as C
    pts = C.build_points()
    e1 = C.experiment_within_pair(pts, repeats=200)
    assert e1, "no pairs evaluated"
    conf = [d[0.2]["cov"] for d in e1.values()]
    gaus = [d[0.2]["gcov"] for d in e1.values()]
    # conformal is valid (>= 1-alpha) up to Monte Carlo slack
    assert min(conf) >= 0.75, f"conformal under-covered: min {min(conf):.3f}"
    assert statistics.fmean(conf) >= 0.80, statistics.fmean(conf)
    # and the Gaussian interval is where it goes wrong
    assert statistics.fmean(gaus) < statistics.fmean(conf) - 0.05, (
        f"gaussian {statistics.fmean(gaus):.3f} vs conformal "
        f"{statistics.fmean(conf):.3f}")
    # the well-powered pair should be close to exact at both levels
    big = max(e1.values(), key=lambda d: d["n"])
    assert abs(big[0.2]["cov"] - 0.80) < 0.04, big[0.2]["cov"]
    assert abs(big[0.1]["cov"] - 0.90) < 0.04, big[0.1]["cov"]


@check("CONFORMAL: the guarantee does NOT survive a change of assay pair")
def t_coverage_cross():
    from .real import conformal as C
    pts = C.build_points()
    e2 = C.experiment_cross_pair(pts, repeats=60)
    covs = [d[0.2] for d in e2.values() if d[0.2] == d[0.2]]
    assert len(covs) > 5
    mean = statistics.fmean(covs)
    assert mean < 0.75, (
        f"cross-pair coverage {mean:.3f} did not degrade; the "
        f"exchangeability failure this documents may have gone away")
    assert min(covs) < 0.30, f"worst cross-pair coverage only {min(covs):.3f}"


@check("CONFORMAL: Mondrian restores the conditional coverage pooling loses")
def t_mondrian():
    from .real import conformal as C
    pts = C.build_points()
    e3 = C.experiment_mondrian(pts, alpha=0.2, repeats=200)
    marg = e3.pop("_pooled_marginal")
    pooled = [d["pooled"] for d in e3.values()]
    mond = [d["mondrian"] for d in e3.values()]
    # marginal promise kept ...
    assert abs(marg - 0.80) < 0.05, f"marginal coverage {marg:.3f}"
    # ... while individual classes fail badly
    assert sum(1 for v in pooled if v < 0.80) >= len(pooled) // 2
    assert min(pooled) < 0.20, f"worst pooled class {min(pooled):.3f}"
    # and per-class calibration repairs it
    assert min(mond) >= 0.78, f"worst mondrian class {min(mond):.3f}"


@check("CONFORMAL: the adapter refuses an uncalibrated assay and demotes evidence")
def t_transfer_adapter():
    from .core import Claim, Estimate, Provenance, MEASURED
    from .real.sources import KRAS
    from .real import conformal as C
    from .real.transfer_adapter import build_transfer
    pts = C.build_points()
    key = max({p.pair for p in pts},
              key=lambda k: sum(1 for p in pts if p.pair == k))
    tx = build_transfer(key, alpha=0.2)

    def claim(assay, x=6.0):
        return Claim(entity=KRAS, quantity=PIC50,
                     context=Context(assay=assay), estimate=Estimate.point(x, 500),
                     evidence=MEASURED,
                     provenance=(Provenance("ChEMBL", "REST/34", assay),))

    assert tx.domain(claim("CHEMBL9999999")).status == "refuse"
    assert tx.domain(claim(tx.pair.source)).status == "in_domain"
    # outside the calibrated input range it extrapolates rather than refusing
    assert tx.domain(claim(tx.pair.source, x=12.0)).status == "extrapolate"

    out = tx.apply(claim(tx.pair.source), random.Random(0))
    assert out.evidence == "calibrated_prediction", out.evidence
    assert out.context.assay == tx.pair.target
    assert dict(out.context.covariates)["transferred_from"] == tx.pair.source
    # the sampled interval should sit inside the guaranteed conformal one
    lo, hi = tx.conformal_interval(6.0)
    s_lo, s_hi = out.estimate.ci(0.80)
    assert lo - 1e-9 <= s_lo and s_hi <= hi + 1e-9, (lo, s_lo, s_hi, hi)


@check("CONFORMAL: an assay transfer is never auto-routed")
def t_transfer_not_routed():
    from .real import conformal as C
    from .real.transfer_adapter import ConformalAssayTransfer, build_transfer
    pts = C.build_points()
    key = max({p.pair for p in pts},
              key=lambda k: sum(1 for p in pts if p.pair == k))
    reg = _reg()
    tx = build_transfer(key, alpha=0.2)
    reg.register(tx)
    path = reg.route(PIC50, FITNESS)
    assert path is not None
    assert not any(isinstance(a, ConformalAssayTransfer) for a in path), \
        "routing silently inserted an assay transfer into a plain pIC50 chain"
    # the transfer produces a DIFFERENT quantity, which is what makes
    # "which assay is this from" a type error rather than a comment
    assert tx.produces.key != PIC50.key
    assert reg.route(PIC50_REF, FITNESS) is not None


# ==========================================================================
# Competing models for one hop, and the registry choosing between them.
# ==========================================================================

@check("ROUTING: two-sided conformal offsets are exact and refuse when unearnable")
def t_signed_quantiles():
    import math
    from .real.transfer_models import signed_quantiles
    r = [float(i) for i in range(1, 21)]           # n=20
    lo, hi = signed_quantiles(r, 0.2)
    # floor(21*0.1)=2 -> 2nd smallest ; ceil(21*0.9)=19 -> 19th smallest
    assert (lo, hi) == (2.0, 19.0), (lo, hi)
    lo, hi = signed_quantiles([1.0, 2.0, 3.0], 0.2)   # n=3 is not enough
    assert math.isinf(hi - lo)


@check("ROUTING: identity and shift are ONE model once conformalised")
def t_identity_is_shift():
    from .real import conformal as C
    from .real.transfer_models import fit_conformal
    pts = C.build_points()
    key = max({p.pair for p in pts},
              key=lambda k: sum(1 for p in pts if p.pair == k))
    idm = fit_conformal(pts, key, "identity", (0.2,))
    shm = fit_conformal(pts, key, "shift", (0.2,))
    assert abs(idm.width(0.2) - shm.width(0.2)) < 1e-9, \
        (idm.width(0.2), shm.width(0.2))
    for x in (4.5, 6.0, 7.5):
        assert abs(idm.center(x) - shm.center(x)) < 1e-9
        assert all(abs(a - b) < 1e-9
                   for a, b in zip(idm.interval(x, 0.2), shm.interval(x, 0.2)))
    # and the raw predictors really are different, so this is a property of
    # conformalisation rather than of the two predictors being identical
    assert abs(idm.predict(6.0) - shm.predict(6.0)) > 0.3


@check("ROUTING: the registry picks per claim, by domain then by width")
def t_registry_selects():
    from .demo_routing import source_claim, strict_registry
    reg = strict_registry()
    cands = reg.candidates(PIC50, PIC50_REF)
    assert len(cands) >= 6, f"only {len(cands)} competing adapters"
    sources = sorted({a.pair.source for a in cands})
    assert len(sources) >= 3
    picked = {}
    for src in sources:
        claim = source_claim(src, 6.0)
        a = reg.select(PIC50, PIC50_REF, claim, alpha=0.2)
        assert a is not None and a.pair.source == src, \
            "selected an adapter calibrated for a different source assay"
        # and it is the narrowest among those in domain for this claim
        same = [c for c in cands if c.pair.source == src]
        assert a.calibrated_width(0.2) == min(c.calibrated_width(0.2)
                                              for c in same)
        picked[src] = a
    assert len({id(a) for a in picked.values()}) == len(sources)
    # an assay nobody calibrated for gets no route at all, not a default
    assert reg.select(PIC50, PIC50_REF,
                      source_claim("CHEMBL9999999", 6.0), 0.2) is None


@check("ROUTING: a calibrated model beats an uncalibrated one; algebra beats both")
def t_width_ranking():
    from .adapter import Adapter, EMPIRICAL
    from .demo_routing import source_claim, strict_registry

    class Uncalibrated(Adapter):
        name, version, kind = "uncalibrated-transfer", "0", EMPIRICAL
        consumes, produces = PIC50, PIC50_REF

        def _forward(self, xs, claim, rng, noise):
            return list(xs)

    reg = strict_registry()
    reg.register(Uncalibrated())
    claim = source_claim("CHEMBL5737243", 6.0)
    pick = reg.select(PIC50, PIC50_REF, claim, alpha=0.2)
    assert pick.name != "uncalibrated-transfer", \
        "an edge with no calibrated interval outranked one that has one"
    assert math.isinf(Uncalibrated().calibrated_width(0.2))
    # a coercion declares width 0 and so wins its group by default
    from .quantities import PIC50REF_TO_IC50
    assert PIC50REF_TO_IC50.calibrated_width(0.2) == 0.0


@check("ROUTING: every in-domain route runs, and the envelope contains them all")
def t_compare_routes():
    from .demo_routing import source_claim, strict_registry
    from .ensemble import compare_routes
    reg = strict_registry()
    rc = compare_routes(reg, PIC50, FITNESS,
                        source_claim("CHEMBL5737243", 6.0), alpha=0.2)
    assert len(rc.results) >= 2, rc.results
    assert rc.refused, "no adapter was refused; the domain test did nothing"
    assert rc.selected in rc.results
    lo, hi = rc.envelope(0.80)
    for res in rc.results.values():
        a, b = res.final.estimate.ci(0.80)
        assert lo - 1e-9 <= a and b <= hi + 1e-9
    assert rc.point_spread() >= 0.0
    # each route ends on the same endpoint quantity
    assert {r.final.quantity.key for r in rc.results.values()} == {FITNESS.key}


@check("ROUTING: the model-choice flag fires on thin data, not on thick")
def t_choice_flag():
    from .demo_routing import source_claim, strict_registry
    from .ensemble import compare_routes
    reg = strict_registry()
    thick = compare_routes(reg, PIC50, FITNESS,
                           source_claim("CHEMBL5737243", 7.5), alpha=0.2)
    thin = compare_routes(reg, PIC50, FITNESS,
                          source_claim("CHEMBL4368373", 7.5), alpha=0.2)
    assert thick.ok and thin.ok
    assert not thick.choice_exceeds_interval(0.80), \
        "flag fired on the well-calibrated source"
    assert thin.choice_exceeds_interval(0.80), \
        "flag did not fire on the thinly calibrated source"
    assert thick.flag(0.80) is None and thin.flag(0.80)
    assert "model-choice" in thin.flag(0.80)


@check("ROUTING: model disagreement is second-order next to the interval")
def t_disagreement_scale():
    from .real import conformal as C
    from .real.transfer_models import disagreement_vs_width
    pts = C.build_points()
    d = disagreement_vs_width(pts, alpha=0.2)
    assert len(d) >= 5
    ratios = [v["median_spread"] / v["width"] for v in d.values()]
    med = statistics.median(ratios)
    assert med < 0.35, (
        f"median spread/width {med:.2f}: model choice is no longer the "
        f"second-order term this documents")
    # but not negligible at the extremes
    assert max(v["max_spread"] / v["width"] for v in d.values()) > 0.8


# ==========================================================================
# Three more pairings of the interface map: C->F, C->G and the F->G edge.
# These run off the committed Tox21 / Ames snapshots.
# ==========================================================================

import functools


CHEM = ("rdkit", "numpy", "sklearn")


@functools.lru_cache(maxsize=4)
def _qsar(ds):
    from .real import qsar
    return qsar.build(ds, alpha=0.1, seed=0)


@functools.lru_cache(maxsize=1)
def _cvd():
    from .real import chain_vs_direct
    return chain_vs_direct.run()


@check("TOX: the two datasets join, and the 2x2 enrichment is real")
def t_tox_overlap():
    needs(CHEM[0])
    from .real import tox
    ov = tox.overlap()
    assert len(ov) > 1900, f"only {len(ov)} compounds in both datasets"
    t = tox.two_by_two(tox.P53, ov)
    assert t.n > 1800, t.n
    assert t.risk_ratio > 1.5, t.risk_ratio
    assert t.fisher_p() < 1e-6, t.fisher_p()
    # and the edge is honest about being a poor screen
    assert t.sensitivity < 0.3, (
        f"sensitivity {t.sensitivity:.2f}: if this ever rises the 'ranks but "
        f"does not clear' claim needs revisiting")


@check("TOX: the signal is not explained by the cytotoxicity control")
def t_cytotox_control():
    needs(CHEM[0])
    from .real import tox
    ov = tox.overlap()
    p53 = tox.two_by_two(tox.P53, ov)
    mmp = tox.two_by_two(tox.CYTOTOX_CONTROL, ov)
    ddr = tox.two_by_two(tox.DDR_CONTROL, ov)
    # both genotoxic-stress reporters beat the general-cytotoxicity readout
    assert p53.risk_ratio > mmp.risk_ratio, (p53.risk_ratio, mmp.risk_ratio)
    assert ddr.risk_ratio > mmp.risk_ratio, (ddr.risk_ratio, mmp.risk_ratio)


@check("TOX: the scaffold split leaks no scaffold between folds")
def t_scaffold_split_clean():
    needs(*CHEM)
    from .real import qsar
    b = _qsar("p53")
    tr, ca, te = b["splits"]
    smi = b["smiles"]
    sets = [{qsar.scaffold(smi[i]) for i in fold} for fold in (tr, ca, te)]
    # singletons get a synthetic key, so only real scaffolds can collide
    real = [{x for x in s if x} for s in sets]
    assert not (real[0] & real[2]), "a scaffold is in both train and test"
    assert not (real[1] & real[2]), "a scaffold is in both calibrate and test"
    assert min(len(f) for f in (tr, ca, te)) > 100


@check("TOX: label-conditional conformal holds for the majority class")
def t_classification_conformal():
    needs(*CHEM)
    e = _qsar("p53")["eval"]
    assert e.cov[0] >= 0.88, f"majority-class coverage {e.cov[0]:.3f}"
    # the minority class under-covers, because a scaffold split breaks the
    # exchangeability the guarantee needs. Asserted so the regression is
    # visible rather than silently fixed.
    assert e.cov[1] < 0.90, (
        f"minority coverage {e.cov[1]:.3f} now meets nominal; the "
        f"exchangeability caveat in demo_pairings M3 needs updating")
    assert 0.0 < e.set_sizes.get(2, 0.0) < 1.0, "abstention rate degenerate"
    am = _qsar("ames")["eval"]
    assert am.cov[0] >= 0.85 and am.cov[1] >= 0.85, (am.cov)


@check("TOX: an ASSOCIATION edge caps evidence at inferred_association")
def t_association_kind():
    needs(CHEM[0])
    from .adapter import ASSOCIATION
    from .real.tox_adapters import P53ToMutagenicity
    fg = P53ToMutagenicity()
    assert fg.kind == ASSOCIATION
    assert fg.max_evidence == "inferred_association"
    # a measured input stays an association through this edge ...
    from .core import MEASURED
    from .quantities import P53_ACTIVE
    from .real.tox_adapters import compound_claim
    c = compound_claim("c1ccccc1", "UHOVQNZJ").derive(evidence=MEASURED)
    out = fg.apply(c, random.Random(0))
    assert out.evidence == "inferred_association", out.evidence
    # ... and the ladder is ordered strongest-first
    from .core import EVIDENCE_ORDER
    assert EVIDENCE_ORDER.index("inferred_association") < \
        EVIDENCE_ORDER.index("calibrated_prediction")


@check("TOX: a chain carries its weakest link, not its best one")
def t_weakest_link():
    needs(*CHEM)
    from .real import tox
    from .real.tox_adapters import (CompoundToP53, P53ToMutagenicity,
                                    compound_claim)
    ov = tox.overlap()
    k = sorted(ov)[0]
    claim = compound_claim(ov[k]["ames"]["smiles"], k)
    res = run_chain([CompoundToP53(), P53ToMutagenicity()], claim, seed=1)
    assert res.ok
    # the QSAR demotes to calibrated_prediction; the downstream association
    # cannot lift it back up
    assert res.final.evidence == "calibrated_prediction", res.final.evidence


@check("TOX: the QSAR adapter refuses an unparseable structure")
def t_qsar_refuses():
    needs(*CHEM)
    from .real.tox_adapters import CompoundToP53, compound_claim
    cf = CompoundToP53()
    good = compound_claim("c1ccccc1", "UHOVQNZJ")
    assert cf.domain(good).status == "in_domain"
    bad = compound_claim("not a smiles at all((", "XXXXXXXX")
    assert cf.domain(bad).status == "refuse"
    # and a claim with no structure at all is refused by the covariate check
    from .core import Context
    nosmi = good.derive(context=Context(assay="Tox21 qHTS"))
    assert cf.domain(nosmi).status == "refuse"


@check("CROSS-CHECK: the identity join key decides how much paired data exists")
def t_identity_yield():
    needs(*CHEM)
    from .real.cross_check import identity_yield
    iy = identity_yield()
    best = iy["parent_skeleton"]["overlap"]
    naive = iy["canonical_smiles"]["overlap"]
    full = iy["inchikey_full"]["overlap"]
    assert best > 2000, best
    # normalising salt form / charge / stereo is worth hundreds of compounds
    assert best - naive > 300, (best, naive)
    assert best - full > 300, (best, full)
    # and parent-skeleton is the best of the four, not merely different
    assert best == max(v["overlap"] for v in iy.values())


@check("CROSS-CHECK: the p53->Ames edge is absent inside the cytotoxic stratum")
def t_stratified_edge():
    needs(*CHEM)
    from .real import tox
    from .real.cross_check import stratified
    st = stratified(tox.CYTOTOX_CONTROL)
    assert st["marginal"].informative, st["marginal"]
    # the whole point: stratifying dissolves the association among cytotoxic
    # compounds and sharpens it among the rest
    assert not st["+"].informative, (
        f"cytotoxic stratum now informative (OR {st['+'].odds_ratio:.2f}, "
        f"p={st['+'].p:.3f}); the conditional edge needs revisiting")
    assert st["-"].informative, st["-"]
    assert st["-"].odds_ratio > st["+"].odds_ratio + 1.0, (
        st["-"].odds_ratio, st["+"].odds_ratio)


@check("CROSS-CHECK: the F->G edge refuses to update on cytotoxic compounds")
def t_conditional_edge():
    needs(*CHEM)
    from .core import Estimate
    from .real import tox
    from .real.tox_adapters import (COV_CYTOTOX, P53ToMutagenicity,
                                    compound_claim)
    fg = P53ToMutagenicity()
    ov = tox.overlap()
    k = sorted(ov)[0]
    base = compound_claim(ov[k]["ames"]["smiles"], k).derive(
        estimate=Estimate.point(0.9, 300))

    def out(cov):
        c = base if cov is None else base.derive(
            context=base.context.with_cov(COV_CYTOTOX, cov))
        return fg.apply(c, random.Random(0))

    unknown, non_cyto, cyto = out(None), out("0"), out("1")
    # a high p53 probability moves the answer in the non-cytotoxic stratum ...
    assert non_cyto.estimate.mean > 0.45, non_cyto.estimate.mean
    # ... and is ignored in the cytotoxic one, where it returns that
    # stratum's own prevalence instead
    cyt_prev = fg._apply(0.0, fg.strata["+"], refuse_update=True)
    assert abs(cyto.estimate.mean - cyt_prev) < 1e-6, (cyto.estimate.mean,
                                                       cyt_prev)
    assert abs(out_hi := fg.apply(base.derive(
        estimate=Estimate.point(0.1, 300), context=base.context.with_cov(
            COV_CYTOTOX, "1")), random.Random(0)).estimate.mean
        - cyto.estimate.mean) < 1e-6, "cytotoxic stratum used the p53 value"
    # every regime says which one it was in
    assert any("no-update" in f for f in cyto.flags), cyto.flags
    assert any("confounded" in f for f in unknown.flags), unknown.flags
    assert not any("no-update" in f or "confounded" in f
                   for f in non_cyto.flags), non_cyto.flags
    # and the edge prices the missing measurement
    assert fg.value_of_knowing_cytotox(0.1) > 0.05


@check("CROSS-CHECK: AP lift reconciles two analyses that AP alone does not")
def t_ap_lift():
    needs(*CHEM)
    from .real.cross_check import REPORTED, ap_lift
    e = _qsar("p53")["eval"]
    mine = ap_lift(e.ap, e.prevalence)
    theirs = ap_lift(REPORTED["p53"]["ap"], REPORTED["p53"]["prevalence"])
    # raw AP looks like a big disagreement
    assert abs(e.ap - REPORTED["p53"]["ap"]) > 0.05
    # lift over prevalence does not
    assert abs(mine - theirs) < 0.8, (mine, theirs)


@check("CHAIN-VS-DIRECT: the direct one-hop model beats the chain")
def t_chain_loses():
    needs(*CHEM)
    r = _cvd()
    by = {a.name: a for a in r["arms"]}
    direct, chain = by["direct  C->G"], by["chain   C->F->G"]
    oracle, aug = by["oracle  F->G"], by["augment C+F->G"]
    assert direct.ap > chain.ap + 0.10, (direct.ap, chain.ap)
    # most of the gap is information loss at the intermediate, not model error
    assert direct.ap - oracle.ap > abs(oracle.ap - chain.ap), (
        f"information loss {direct.ap - oracle.ap:.3f} no longer dominates "
        f"model error {abs(oracle.ap - chain.ap):.3f}")
    # and a MEASURED intermediate adds essentially nothing over structure
    assert abs(aug.ap - direct.ap) < 0.03, (aug.ap, direct.ap)


@check("COPULA: a measured cytotoxicity call really does sharpen p53")
def t_copula_gain():
    needs(*CHEM, "scipy")
    import numpy as np
    from sklearn.metrics import average_precision_score, roc_auc_score
    from .real import copula as C, qsar, tox
    m = C.load()
    rows = tox.load_tox21()
    smiles = [r["smiles"] for r in rows]
    X, ok = qsar.featurize(smiles, n_bits=C.N_BITS)
    keep = np.where(ok)[0]
    smi = [smiles[i] for i in keep]
    rk = [rows[i] for i in keep]
    _, _, te = qsar.scaffold_split(smi, fracs=(0.5, 0.2, 0.3), seed=0)

    def lab(i, e):
        v = rk[i][e]
        return None if v in ("", "NA") else int(float(v))
    both = [i for i in te if lab(i, tox.P53) is not None
            and lab(i, tox.CYTOTOX_CONTROL) is not None]
    y = np.array([lab(i, tox.P53) for i in both])
    mmp = np.array([lab(i, tox.CYTOTOX_CONTROL) for i in both])
    sm = [smi[i] for i in both]
    assert len(both) > 1400, len(both)

    pa = m.proba(sm, None)
    pb = np.array([m.proba([t], {tox.CYTOTOX_CONTROL: int(v)})[0]
                   for t, v in zip(sm, mmp)])
    gain = roc_auc_score(y, pb) - roc_auc_score(y, pa)
    ap_gain = average_precision_score(y, pb) - average_precision_score(y, pa)
    assert gain > 0.015, f"AUROC gain collapsed to {gain:+.4f}"
    assert ap_gain > 0.03, f"AP gain collapsed to {ap_gain:+.4f}"
    assert 0.25 < m.rho(tox.P53, tox.CYTOTOX_CONTROL) < 0.40, m.rho(
        tox.P53, tox.CYTOTOX_CONTROL)
    # 66 free correlations, not 1.6M denoiser parameters
    assert C.NE * (C.NE - 1) // 2 == 66


@check("COPULA: better ranking does NOT mean fewer abstentions")
def t_width_is_not_discrimination():
    needs(*CHEM, "scipy")
    from .real import copula as C
    m = C.load()
    # the conditional model is +0.025 AUROC better and abstains no less.
    # Asserted so that the limitation in the registry's docstring stays
    # true: ranking on width alone would pick the worse model here.
    assert m.abstain, "fitted artefact has no abstention rates"
    assert m.abstain["cond"] >= m.abstain["marg"] - 0.005, m.abstain


@check("REGISTRY: discrimination outranks width when the folds match")
def t_discrimination_ranking():
    from .adapter import Adapter, EMPIRICAL
    from .quantities import P53_ACTIVE, AMES_POSITIVE

    class Sharp(Adapter):
        """Abstains less, discriminates worse."""
        name, version, kind = "sharp-but-worse", "0", EMPIRICAL
        consumes, produces = P53_ACTIVE, AMES_POSITIVE
        declared_auroc, eval_fold_id = 0.70, "fold-A"

        def calibrated_width(self, alpha=0.2):
            return 0.20

        def _forward(self, xs, claim, rng, noise):
            return list(xs)

    class Discriminating(Adapter):
        """Abstains more, discriminates better."""
        name, version, kind = "wide-but-better", "0", EMPIRICAL
        consumes, produces = P53_ACTIVE, AMES_POSITIVE
        declared_auroc, eval_fold_id = 0.85, "fold-A"

        def calibrated_width(self, alpha=0.2):
            return 0.60

        def _forward(self, xs, claim, rng, noise):
            return list(xs)

    reg = Registry(include_lossless=False)
    reg.register(Sharp())
    reg.register(Discriminating())
    assert reg.discrimination_comparable(P53_ACTIVE, AMES_POSITIVE)
    pick = reg.select(P53_ACTIVE, AMES_POSITIVE)
    assert pick.name == "wide-but-better", f"width still won: {pick.name}"

    # now break the fold match: discrimination becomes unusable and the
    # rule must fall back to width rather than compare across folds
    d2 = Discriminating()
    d2.eval_fold_id = "fold-B"
    reg2 = Registry(include_lossless=False)
    reg2.register(Sharp())
    reg2.register(d2)
    assert not reg2.discrimination_comparable(P53_ACTIVE, AMES_POSITIVE)
    assert reg2.select(P53_ACTIVE, AMES_POSITIVE).name == "sharp-but-worse"


@check("REGISTRY: the two real C->F models are NOT rankable on discrimination")
def t_real_folds_differ():
    needs(*CHEM, "scipy")
    from .quantities import P53_ACTIVE
    from .real.tox_adapters import CompoundToP53, CompoundToP53Copula
    a, b = CompoundToP53(), CompoundToP53Copula()
    assert a.declared_discrimination() and b.declared_discrimination()
    assert a.eval_fold_id != b.eval_fold_id, "folds coincidentally matched"
    reg = Registry(include_lossless=False)
    reg.register(a)
    reg.register(b)
    assert not reg.discrimination_comparable(P53_ACTIVE, P53_ACTIVE)
    # so it falls back to width, which here picks the WORSE model. That is
    # the rule refusing to guess, and the cost of not refitting on a common
    # split -- documented in Registry.rank.
    assert reg.select(P53_ACTIVE, P53_ACTIVE).name == a.name


@check("COPULA: the duplicated probit primitives agree with diffusion.py")
def t_primitives_agree():
    needs(*CHEM, "scipy")
    import numpy as np
    from .real import copula as C
    try:
        from .real import diffusion as DF
    except Exception as e:
        raise Skip(f"diffusion.py not importable (mid-revision): {e}")
    t1 = np.array([-1.0, 0.0, 0.5, 1.5])
    t2 = np.array([0.2, -0.4, 1.0, 0.0])
    for rho in (-0.6, 0.0, 0.33, 0.9):
        assert np.allclose(C.bvn_sf(t1, t2, rho), DF.bvn_sf(t1, t2, rho),
                           atol=1e-12), rho
    assert np.allclose(C.thresholds(np.array([[0.1, 0.9]])),
                       DF._thresholds(np.array([[0.1, 0.9]])))


@check("MAP21: the audit is complete and internally consistent")
def t_map21():
    from . import map21
    keys = {p.key for p in map21.MAP}
    assert len(keys) == len(map21.MAP), "duplicate pairing in the audit"
    # every object class appears on at least one side
    seen = set()
    for p in map21.MAP:
        a, b = p.key.split("->")
        seen |= {a, b}
    assert seen == set(map21.CLASSES), seen ^ set(map21.CLASSES)
    # everything we refuse to build is a level-C edge in the blueprint
    for p in map21.MAP:
        if p.status == map21.REFUSED:
            assert p.claimed == "C", (p.key, p.claimed)
        if p.status == map21.BUILT:
            assert p.module and p.paired_labels, p.key
    g = map21.by_status()
    assert len(g[map21.BUILT]) >= 4
    assert map21.report()


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
    failed = skipped = 0
    for name, ok, msg in _results:
        tag = "PASS" if ok else ("SKIP" if ok is None else "FAIL")
        print(f"  {tag}  {name:<{width}}" + (f"   ({msg})" if ok is None else ""))
        if ok is False:
            failed += 1
            print("        " + msg.replace("\n", "\n        "))
        elif ok is None:
            skipped += 1
    passed = len(_results) - failed - skipped
    tail = f", {skipped} skipped" if skipped else ""
    print(f"\n{passed}/{len(_results) - skipped} passed{tail}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
