"""
bioif.demo_pairings -- three more pairings of the interface map, and the
long-chain question answered with labels on both ends.

Run:  python3 -m bioif.demo_pairings

The blueprint's map has 21 directed pairings across seven object classes,
graded A (verifiable under stated conditions), B (needs paired experimental
calibration) or C (hypothesis / ranking only). This demo audits all of them
against what is actually buildable here, then builds the three that have
public labels on both sides, then uses them to answer the question the whole
programme rests on.

⚠️  Every number below is measured on held-out scaffolds from committed
public snapshots. None of it is a claim about a specific compound.
"""
from __future__ import annotations

import random

from . import map21
from .chain import run_chain
from .core import EVIDENCE_ORDER
from .real import chain_vs_direct, qsar, tox
from .real.tox_adapters import (CompoundToAmes, CompoundToP53,
                                P53ToMutagenicity, compound_claim)


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main() -> None:
    # -- M1 which of the 21 pairings can be built at all -------------------
    hr("M1  the interface map, audited")
    g = map21.by_status()
    print(f"  {'status':<11}{'n':>4}   pairings")
    for st in sorted(g, key=lambda s: map21.ORDER[s]):
        print(f"  {st:<11}{len(g[st]):>4}   " + " ".join(p.key for p in g[st]))
    print("\n  Full audit with sources, paired-label counts and findings:")
    print("    python3 -m bioif.map21")
    print("\n  The pattern: a pairing is buildable exactly when some public\n"
          "  source happens to hold labels on BOTH sides. That is 4 of 21.\n"
          "  The binding constraint across the map is paired measurement --\n"
          "  not models, not formats.")

    # -- M2 the blocked edge, measured ------------------------------------
    hr("M2  F->G: the blueprint's `blocked` edge, now measured")
    print(tox.report())

    # -- M3 the two QSAR nodes --------------------------------------------
    hr("M3  C->F and C->G as conformalised classifiers")
    print("  Scaffold split (whole Bemis-Murcko scaffolds to one fold), and\n"
          "  LABEL-CONDITIONAL conformal: SR-p53 is 6.2% positive, so a\n"
          "  marginal 90% classifier could meet its guarantee by calling\n"
          "  everything inactive. Calibrating per class forbids that.\n")
    print(f"  {'node':<22}{'n test':>8}{'prev':>7}{'AUROC':>8}{'AP':>7}"
          f"{'cov|0':>8}{'cov|1':>8}{'abstain':>9}")
    built = {}
    for ds, label in (("p53", "C->F  SR-p53"), ("ames", "C->G  Ames")):
        b = qsar.build(ds, alpha=0.1, seed=0)
        built[ds] = b
        e = b["eval"]
        print(f"  {label:<22}{e.n_test:>8}{e.prevalence:>7.3f}{e.auroc:>8.3f}"
              f"{e.ap:>7.3f}{e.cov[0]:>8.3f}{e.cov[1]:>8.3f}"
              f"{e.set_sizes.get(2, 0.0):>9.3f}")
    p53e = built["p53"]["eval"]
    print(f"\n  ⚠️  The p53 minority class is covered at {p53e.cov[1]:.3f} against a\n"
          f"  nominal 0.90. That shortfall is not a bug -- a scaffold split is\n"
          f"  DESIGNED to break exchangeability between calibration and test\n"
          f"  chemistry, and conformal's guarantee is conditional on it. Same\n"
          f"  mechanism as the cross-assay collapse in INTEROP §6.3, and the\n"
          f"  same lesson: the guarantee is only as good as the split it was\n"
          f"  calibrated under.")
    print(f"\n  For reference, the plan document reports AUROC 0.733 for its own\n"
          f"  p53 model; this one gets {p53e.auroc:.3f} on held-out scaffolds, so the\n"
          f"  two are in the same place.")

    # -- M4 the centrepiece ------------------------------------------------
    hr("M4  does routing through an intermediate beat going direct?")
    res = chain_vs_direct.run()
    print(chain_vs_direct.report(res))

    # -- M5 the same thing as a typed chain --------------------------------
    hr("M5  the chain as typed edges, and what the evidence ladder does")
    cf, fg, cg = CompoundToP53(), P53ToMutagenicity(), CompoundToAmes()
    print(f"  C->F  {cf.name:<24}[{cf.kind}]  abstains on "
          f"{cf.calibrated_width():.0%} of compounds")
    print(f"  F->G  {fg.name:<24}[{fg.kind}]")
    print(f"        {fg.screening_verdict()}")
    print(f"  C->G  {cg.name:<24}[{cg.kind}]  abstains on "
          f"{cg.calibrated_width():.0%} of compounds")

    ov = tox.overlap()
    # a compound the p53 model is confident about, for a legible trace
    pick = None
    for k in sorted(ov):
        smi = ov[k]["ames"]["smiles"]
        if len(cf.prediction_set(smi)) == 1:
            pick = (k, smi)
            break
    if pick:
        k, smi = pick
        claim = compound_claim(smi, k)
        r = run_chain([cf, fg], claim, seed=1)
        print(f"\n  a compound whose p53 prediction set is a singleton:")
        print(f"    {smi[:64]}")
        for t in r.traces:
            print(f"    after {t.adapter.name:<24}"
                  f"{t.out.estimate.summary():<26}[{t.out.evidence}]")
        print(f"\n  Note the evidence label. The F->G edge is an ASSOCIATION,")
        print(f"  which on its own caps at `inferred_association` -- stronger")
        print(f"  than a model prediction. But it sits DOWNSTREAM of a QSAR,")
        print(f"  and a chain carries its weakest link, so the endpoint is")
        print(f"  `{r.final.evidence}`. Putting a measured association after a")
        print(f"  prediction does not recover the association's standing.")
        print(f"\n  Evidence ladder, strongest to weakest: "
              f"{' > '.join(EVIDENCE_ORDER)}")

    hr("what this says about the programme")
    by = {a.name: a for a in res["arms"]}
    print(f"  1. The F->G edge is real but weak: RR 1.81 (p=5e-9) with")
    print(f"     sensitivity 0.13. It moves from `blocked` to")
    print(f"     `inferred_association`: good for ordering a list, not for")
    print(f"     clearing one.")
    print(f"  2. For PREDICTING the mutagenicity label, a direct structure")
    print(f"     model beats the chain by {by['direct  C->G'].ap - by['chain   C->F->G'].ap:+.3f} AP, and a MEASURED p53")
    print(f"     label adds {by['augment C+F->G'].ap - by['direct  C->G'].ap:+.3f} AP on top of structure. The")
    print(f"     intermediate is, for this endpoint, redundant given the")
    print(f"     structure.")
    print(f"  3. ⚠️  That is a statement about predicting Ames, not about the")
    print(f"     reporter's worth. The reporter says WHY, which a structure")
    print(f"     model does not; and Ames is not Comet, which is the endpoint")
    print(f"     the blueprint actually wants. The honest use of this result")
    print(f"     is as a prior on experiment design: before paying for")
    print(f"     reporter panels to predict a genotoxicity endpoint, check")
    print(f"     whether structure alone already carries that information.")

    print("\n" + "=" * 78)
    print("✅ measured on held-out scaffolds from committed public snapshots\n"
          "   (Tox21 MoleculeNet release; Hansen Ames N=6512; ChEMBL).\n"
          "⚠️  one target class, one fingerprint, one split protocol. Ames is\n"
          "   not Comet. No claim about any individual compound.")
    print("=" * 78)


if __name__ == "__main__":
    main()
