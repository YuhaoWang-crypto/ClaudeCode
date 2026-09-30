"""
bioif.demo -- the minimal end-to-end demonstration.

Run:  python3 -m bioif.demo

Six scenarios, each isolating one failure mode that breaks long causal
chains in practice:

  S1  a four-hop genetic chain runs, and its uncertainty is attributed to
      the individual hops
  S2  the same chain asked about a cell system outside the calibration set
      -- extrapolates loudly rather than answering quietly
  S3  a chemical chain reaches the same endpoint through a *different* set
      of adapters, found automatically by routing, with the unit algebra
      inserted for free
  S4  the same chemical chain with the dose omitted -- refuses, because an
      occupancy without a dose is not a number
  S5  naive point-estimate composition vs the typed chain, on identical
      inputs
  S6  the actual deliverable: which single measurement most tightens the
      endpoint per unit cost

⚠️  All adapters are illustrative stubs (see adapters_demo.py). The outputs
below demonstrate interface behaviour, not biology.
"""
from __future__ import annotations

import random

from .adapters_demo import DEMO_ADAPTERS
from .chain import Experiment, rank_experiments, run_chain
from .core import Claim, Context, Entity, Estimate, Provenance, check_identity
from .quantities import DELTA_PSI, FITNESS, PIC50
from .registry import Registry

N = 6000
SEED = 7

CALIBRATED = Context(system="depmap:ACH-000019", assay="CRISPR-KO fitness",
                     time_h=120.0)
UNCALIBRATED = CALIBRATED.with_(system="depmap:ACH-000999")


def build_registry() -> Registry:
    reg = Registry()                       # lossless unit algebra pre-loaded
    for a in DEMO_ADAPTERS:
        reg.register(a)
    return reg


def genetic_source(ctx: Context) -> Claim:
    rng = random.Random(SEED)
    return Claim(
        entity=Entity("hgvs_g", "NC_000017.11:g.7676154G>A", build="GRCh38"),
        quantity=DELTA_PSI,
        context=ctx.with_(assay="splice model, in silico"),
        estimate=Estimate.normal(-0.62, 0.12, N, rng),
        provenance=(Provenance("splice-model", "0.1-stub",
                               "illustrative, not fitted"),),
    )


def chemical_source(ctx: Context) -> Claim:
    rng = random.Random(SEED + 1)
    return Claim(
        entity=Entity("inchikey", "RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
        quantity=PIC50,
        context=ctx.with_(assay="affinity model, in silico"),
        estimate=Estimate.normal(7.40, 0.60, N, rng),
        provenance=(Provenance("affinity-model", "0.1-stub",
                               "illustrative, not fitted"),),
    )


def naive_point_chain(path, source) -> float:
    """
    What the same wiring produces when every hop is a point estimate and no
    hop is allowed to refuse: one confident number, no error bar, no flags,
    no record that a cross-level assumption was made four times.
    """
    rng = random.Random(0)
    claim = source.derive(estimate=Estimate.point(source.estimate.mean))
    for a in path:
        claim = a.apply(claim, rng, noise=False)
    return claim.estimate.mean


def hr(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def main() -> None:
    reg = build_registry()

    hr("REGISTRY")
    print(reg.describe())

    # ---- identifier lint, the cheapest bug class ------------------------
    hr("S0  identifier lint")
    bad = Entity("ensembl_gene", "ENSG00000141510")          # no build
    for p in check_identity(bad):
        print(f"  REJECTED: {p}")
    good = Entity("ensembl_gene", "ENSG00000141510", "GRCh38.p14/r112")
    print(f"  accepted: {good}")

    # ---- S1 genetic arm, in domain --------------------------------------
    hr("S1  genetic chain, calibrated system")
    gpath = reg.route(DELTA_PSI, FITNESS)
    print("routed: " + " -> ".join(a.name for a in gpath) + "\n")
    g1 = run_chain(gpath, genetic_source(CALIBRATED), seed=SEED)
    print(g1.report())

    # ---- S2 genetic arm, out of domain ----------------------------------
    hr("S2  same chain, system outside the calibration set")
    g2 = run_chain(gpath, genetic_source(UNCALIBRATED), seed=SEED)
    print(g2.report())
    if g1.ok and g2.ok:
        w1 = g1.final.estimate.ci()[1] - g1.final.estimate.ci()[0]
        w2 = g2.final.estimate.ci()[1] - g2.final.estimate.ci()[0]
        print(f"\n  90% interval widens {w1:.3f} -> {w2:.3f} "
              f"({w2 / w1:.2f}x) purely from the declared context mismatch.")

    # ---- S3 chemical arm, routed automatically --------------------------
    hr("S3  chemical chain to the same endpoint (auto-routed)")
    cpath = reg.route(PIC50, FITNESS)
    print("routed: " + " -> ".join(a.name for a in cpath))
    print("  (the two unit hops were inserted by the registry, not by hand)\n")
    dosed = CALIBRATED.with_(dose_uM=0.30)
    c1 = run_chain(cpath, chemical_source(dosed), seed=SEED)
    print(c1.report())

    # ---- the comparison the convergence node buys -----------------------
    if g1.ok and c1.ok:
        hr("S3b  what the shared quantity actually buys")
        g_res = [t.out for t in g1.traces
                 if t.out.quantity.name == "residual_target_activity"][0]
        c_res = [t.out for t in c1.traces
                 if t.out.quantity.name == "residual_target_activity"][0]
        print(f"  genetic arm   {g_res.entity}  {g_res.estimate.summary()}")
        print(f"                {g_res.context.describe()}")
        print(f"  chemical arm  {c_res.entity}  {c_res.estimate.summary()}")
        print(f"                {c_res.context.describe()}")
        same = g_res.entity == c_res.entity and \
            g_res.quantity.key == c_res.quantity.key
        print(f"\n  same entity AND same quantity: {same}")
        print("\n  That is the whole convergence condition. Both arms end on the\n"
              "  same target entity with the same reference state, and the thing\n"
              "  that differed between them -- a variant on one side, a compound\n"
              "  at a dose on the other -- was pushed into context rather than\n"
              "  left in the identifier. Only now can a model fitted on genetic\n"
              "  perturbation score a molecule. The bridge flags on both numbers\n"
              "  are the price of doing it, and they stay attached.")

    # ---- S4 refusal -----------------------------------------------------
    hr("S4  same chemical chain, dose omitted")
    c2 = run_chain(cpath, chemical_source(CALIBRATED), seed=SEED)
    print(c2.report())

    # ---- S5 naive vs typed ----------------------------------------------
    hr("S5  naive point composition vs the typed chain")
    naive = naive_point_chain(gpath, genetic_source(UNCALIBRATED))
    print(f"  naive pipeline  : fitness_effect = {naive:+.3f}")
    print("                    (no interval, no flags, no domain check,")
    print("                     and it answered a question about a cell")
    print("                     system no link was ever fitted on)")
    if g2.ok:
        lo, hi = g2.final.estimate.ci()
        print(f"  typed chain     : fitness_effect = "
              f"{g2.final.estimate.summary()}")
        print(f"                    {len(g2.final.flags)} assumption flags, "
              f"1 extrapolation warning")
        print(f"  The point estimates agree to ~{abs(naive - g2.final.estimate.mean):.2f}. "
              f"The difference is\n  that only one of them tells you the answer "
              f"spans {hi - lo:.2f} units and\n  rests on three unvalidated "
              f"cross-level bridges.")

    # ---- S6 value of information ----------------------------------------
    hr("S6  which measurement to actually run")
    catalogue = {
        "psi->transcript_dosage": Experiment(
            "psi->transcript_dosage",
            "RT-PCR / targeted RNA-seq of isoform ratio in this system", 4.0),
        "rna->protein": Experiment(
            "rna->protein",
            "targeted proteomics (PRM) of the gene product, matched samples",
            9.0),
        "protein->residual_activity": Experiment(
            "protein->residual_activity",
            "pathway activity reporter across a knockdown titration", 12.0),
        "residual_activity->fitness": Experiment(
            "residual_activity->fitness",
            "CRISPR-KO fitness arm in THIS cell system", 20.0),
    }
    rows = rank_experiments(g2, catalogue)
    if rows:
        print(f"  {'rank':<5}{'variance removed / cost':<26}{'share':<8}"
              f"experiment")
        for i, (eff, share, exp) in enumerate(rows, 1):
            print(f"  {i:<5}{eff:<26.4f}{share * 100:<8.1f}{exp.description}")
        print("\n  Read this as the output of the pipeline. The endpoint number\n"
              "  is a hypothesis; the ranking above is the decision it supports.")

    print("\n" + "=" * 74)
    print("⚠️  Every constant above is an illustrative stub. What is being\n"
          "    demonstrated is interface behaviour: typed seams, automatic\n"
          "    routing, honest widening, refusal, and variance attribution.")
    print("=" * 74)


if __name__ == "__main__":
    main()
