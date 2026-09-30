"""
bioif.demo_real -- the same chain, with the affinity source replaced by a
REAL data source.

Run:  python3 -m bioif.demo_real

The stub `affinity-model@0.1-stub` is gone. In its place: measured
bioactivity from the EMBL-EBI ChEMBL REST API, committed as a snapshot under
`bioif/real/_snapshot/` so this is reproducible offline. Everything
downstream of the source (occupancy, residual activity, fitness) is still the
same stub set from `adapters_demo.py`, and is still labelled as such.

The question this demo answers is not "is the prediction better". It is:
**what does the interface do differently when the data is real?** Four
things, and every one of them is a behaviour the stub demo could not
exercise because a stub emits exactly one tidy number.
"""
from __future__ import annotations

import statistics

from .adapters_demo import (ActivityToFitness, OccupancyToResidualActivity,
                            ProteinToResidualActivity)
from .chain import run_chain
from .core import MEASURED, Context
from .quantities import FITNESS, PIC50, PKD, PKI, RESIDUAL_ACTIVITY
from .real import chembl, sources
from .real.heterogeneity import report as heterogeneity_report
from .registry import Registry

SEED = 7


def build_registry() -> Registry:
    reg = Registry()
    occ = OccupancyToResidualActivity()
    occ.target = sources.KRAS                 # both arms converge on KRAS
    reg.register(occ)
    reg.register(ProteinToResidualActivity())
    reg.register(ActivityToFitness())
    reg.register(sources.KiToIC50())
    return reg


def hr(t):
    print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def main() -> None:
    reg = build_registry()

    # -- R1 identity resolution is allowed to refuse ----------------------
    hr("R1  resolving a target name against ChEMBL")
    amb = chembl.resolve_target("KRAS", target_type="")
    print(f"  resolve_target('KRAS')  with no target_type constraint "
          f"->  ok={amb.ok}")
    print(f"  ! {amb.reason}")
    print(f"  the {len(amb.candidates)} candidates ChEMBL actually returns:")
    for h in amb.candidates:
        print(f"      {h.chembl_id:<16}{h.target_type:<28}{h.pref_name}")
    res = chembl.resolve_target("KRAS")
    print(f"\n  with organism + target_type constrained -> "
          f"{res.resolved.chembl_id} ({res.resolved.pref_name})")
    print("  A gene symbol is not an identifier. Taking the first row here is\n"
          "  how a compound measured against a protein-protein interaction\n"
          "  ends up scored as a direct inhibitor of the single protein.")

    # -- R2 what a real source actually returns ---------------------------
    hr("R2  the real source returns a SET of claims, not a number")
    claims = sources.measured_affinity_claims()
    print(f"  compound {chembl.DEMO_MOLECULE} on {chembl.DEMO_TARGET}")
    print(f"  {len(claims)} measured records, "
          f"{len(sources.group_by_assay(claims))} distinct assays")
    for c in claims:
        cov = dict(c.context.covariates)
        print(f"    {c.context.assay:<16}{cov['standard_type']:<6}"
              f"{c.quantity.unit:<7}{c.estimate.mean:5.2f}   "
              f"{c.provenance[0].note[:62]}")
    print("\n  Note the quantity column: these are NOT all the same quantity.")

    # -- R3 the contract refuses the conversion everyone reaches for ------
    hr("R3  pKi -> pIC50 looks like arithmetic; the contract says no")
    ic50 = [c for c in claims if c.quantity.key == PIC50.key]
    kd = [c for c in claims if c.quantity.key == PKD.key]
    print(f"  {len(ic50)} pIC50 records route into the chain.")
    print(f"  {len(kd)} pKd records do not: no edge from pKd to pIC50 is "
          f"registered,")
    print(f"     and none should be without the mechanism and [S]/Km.")
    print(f"  route(pKd -> fitness): {reg.route(PKD, FITNESS)}")
    if kd:
        k = kd[0].derive(quantity=PKI)      # pretend it were a Ki, to try
        v = sources.KiToIC50().domain(k)
        print(f"\n  and the Cheng-Prusoff bridge, offered a real record:")
        print(f"    {v.status.upper()}: {v.reason}")
        supplied = k.derive(context=k.context.with_cov("substrate_over_km", "10")
                            .with_cov("mechanism", "competitive"))
        print(f"    with [S]/Km and mechanism supplied by hand: "
              f"{sources.KiToIC50().domain(supplied).status}")
        print("    -- which is the point: the conversion is a scientific claim\n"
              "       about the assay, so somebody has to make it on the record.")

    # -- R4 run the real chain, once per assay ----------------------------
    hr("R4  the same chain, run once per assay that measured this compound")
    path = reg.route(PIC50, FITNESS)
    print("  routed: " + " -> ".join(a.name for a in path))
    ctx_add = dict(dose_uM=0.30, system="depmap:ACH-000019", time_h=120.0)
    ends = []
    print(f"\n  {'assay':<16}{'pIC50':>7}   {'fitness_effect (90% CI)':<28}evidence")
    for c in sorted(ic50, key=lambda x: x.estimate.mean):
        src = c.derive(context=Context(
            assay=c.context.assay, covariates=c.context.covariates, **ctx_add))
        r = run_chain(path, src, seed=SEED)
        if not r.ok:
            continue
        ends.append(r.final.estimate.mean)
        print(f"  {c.context.assay:<16}{c.estimate.mean:7.2f}   "
              f"{r.final.estimate.summary():<28}{r.final.evidence}")
    if len(ends) > 1:
        vals = [c.estimate.mean for c in ic50]
        in_span, out_span = max(vals) - min(vals), max(ends) - min(ends)
        print(f"\n  Same compound. Same target. Same chain. Same dose.")
        print(f"  Endpoint ranges {min(ends):+.3f} to {max(ends):+.3f} "
              f"({out_span:.3f} gene-effect units) purely from\n"
              f"  which assay's number you fed in.")
        print(f"\n  But note the compression: {in_span:.2f} log units of input "
              f"({10 ** in_span:.0f}x in potency)\n"
              f"  became {out_span:.3f} units of output. At {ctx_add['dose_uM']} uM "
              f"every one of these\n  affinities is already near-saturating, so "
              f"the endpoint barely moves.\n"
              f"  Read that as the chain telling you where NOT to spend money:\n"
              f"  at this dose, a better affinity number buys almost nothing.")

    # -- R5 what pooling hides -------------------------------------------
    hr("R5  what the naive 'just ask the database for the affinity' does")
    pooled = sources.naive_pooled_claim(claims)
    psrc = pooled.derive(context=Context(
        assay=pooled.context.assay, covariates=pooled.context.covariates,
        **ctx_add))
    pr = run_chain(path, psrc, seed=SEED)
    print(f"  pooled source : pIC50 {pooled.estimate.summary()}  "
          f"[{pooled.evidence}]")
    print(f"  pooled endpoint: fitness {pr.final.estimate.summary()}  "
          f"[{pr.final.evidence}]")
    print(f"  per-assay endpoints spanned {max(ends) - min(ends):.3f} units; "
          f"the pooled\n  run reports a single interval of "
          f"{pr.final.estimate.ci()[1] - pr.final.estimate.ci()[0]:.3f}.")
    for f in pooled.flags:
        print(f"    flag: {f}")

    # -- R6 evidence can only degrade -------------------------------------
    hr("R6  evidence level along the chain")
    src = ic50[0].derive(context=Context(
        assay=ic50[0].context.assay, covariates=ic50[0].context.covariates,
        **ctx_add))
    r = run_chain(path, src, seed=SEED)
    print(f"  {'source':<36}{src.evidence}")
    for t in r.traces:
        print(f"  after {t.adapter.name:<30}{t.out.evidence}")
    print("\n  The source is a real measurement. The endpoint is a mechanistic\n"
          "  hypothesis, because a bridge sits between them. No downstream\n"
          "  machinery can restore the label, which is the intended behaviour.")

    # -- R7 the real number that justifies all of the above ---------------
    hr("R7  does assay identity actually matter? (measured, not asserted)")
    print(heterogeneity_report())

    print("\n" + "=" * 74)
    print("✅ source: measured ChEMBL bioactivity (snapshot committed).\n"
          "⚠️  occupancy / residual-activity / fitness adapters are still the\n"
          "    illustrative stubs from adapters_demo.py. The endpoint numbers\n"
          "    above are therefore NOT predictions; the spread between them is\n"
          "    the result worth reading.")
    print("=" * 74)


if __name__ == "__main__":
    main()
