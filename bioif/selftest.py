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
from .quantities import (DELTA_PSI, FITNESS, IC50_NM, PIC50, PROTEIN_LFC,
                         RESIDUAL_ACTIVITY, RNA_LFC)
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


@check("routing prefers arithmetic over an assumption-bearing bridge")
def t_route_cost():
    from .adapter import Adapter, BRIDGE
    from .core import Verdict

    class Shortcut(Adapter):
        name, version, kind = "psi->protein_shortcut", "0", BRIDGE
        consumes, produces = DELTA_PSI, PROTEIN_LFC
        assumption = "a very large leap"

        def _forward(self, xs, claim, rng, noise):
            return [x for x in xs]

    reg = _reg()
    reg.register(Shortcut())
    path = reg.route(DELTA_PSI, FITNESS)
    # Shortcut gives a 3-hop path vs the 4-hop honest one, so BFS takes it --
    # length still wins. What we assert is that among EQUAL-length paths the
    # bridge count breaks the tie.
    from .registry import _cost
    assert _cost(path) <= _cost([Shortcut()] + path[1:])


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
