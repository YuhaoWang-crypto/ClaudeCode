"""Command-line entry point.

    python -m mri_cd8_til inventory      # what public data exists, verified live
    python -m mri_cd8_til intervals      # confidence intervals on the published AUCs
    python -m mri_cd8_til plan           # the multi-centre plan, with measured numbers
"""

from __future__ import annotations

import argparse

from .cohorts import REFERENCE_STUDY, REGISTRY, LabelTier, by_tier, total_public_mri
from .labels.schema import LabelTask
from .modeling.metrics import auc_interval_from_reported


def cmd_inventory(_: argparse.Namespace) -> int:
    print(
        f"Public breast-MRI patients: {total_public_mri()} - "
        f"{total_public_mri(count_overlapping=True)}\n"
        "  (a range, because ACRIN-6698 is the I-SPY2 DWI sub-study and its 385\n"
        "   patients cannot be matched to the ISPY2 IDs to confirm or rule out overlap)\n"
    )
    for tier, heading in (
        (LabelTier.SPATIAL_TIL, "TIER A -- spatially resolved TIL label + imaging"),
        (LabelTier.TRANSCRIPTOMIC, "TIER B -- transcriptomic immune label + imaging"),
        (LabelTier.IMAGING_ONLY, "TIER C -- imaging only (segmentation / harmonisation)"),
    ):
        print(f"{heading}")
        print("-" * len(heading))
        for cohort in by_tier(tier):
            paired = f"{cohort.n_paired} paired" if cohort.n_paired else "no paired label"
            print(f"  {cohort.name}")
            print(f"      MRI n={cohort.n_with_mri:<5} {paired:<18} centres: {cohort.centers}")
            print(f"      link: {cohort.identifier_link}")
            for caveat in cohort.caveats:
                print(f"      ! {caveat}")
        print()
    print("Re-measure everything with: python scripts/refresh_inventory.py")
    return 0


def cmd_intervals(_: argparse.Namespace) -> int:
    print(f"Source: {REFERENCE_STUDY['citation']}")
    print(
        f"Design: n={REFERENCE_STUDY['n_total']} "
        f"({REFERENCE_STUDY['n_train']}/{REFERENCE_STUDY['n_validation']}), "
        f"{REFERENCE_STUDY['centers']} centre, "
        f"external validation: {REFERENCE_STUDY['external_validation']}\n"
    )
    for name, auc, n_pos, n_neg in (
        ("RM-whole  training  ", 0.973, 46, 91),
        ("RM-whole  validation", 0.985, 15, 30),
        ("RM-peri   training  ", 0.993, 46, 45),
        ("RM-peri   validation", 0.984, 15, 15),
    ):
        print(f"  {name}  {auc_interval_from_reported(auc, n_pos, n_neg).describe()}")
    print(
        "\nClass counts are assumed roughly balanced (not published per split);"
        "\nthe widths depend on that assumption, the conclusion does not."
    )
    return 0


def cmd_plan(_: argparse.Namespace) -> int:
    tier_a = REGISTRY["tcga_brca"]
    tier_b = REGISTRY["ispy2"]
    print("MULTI-CENTRE PLAN\n")
    print(f"  Train  : {tier_b.name}")
    print(f"           n={tier_b.n_paired}, {tier_b.centers}")
    print(f"           label: {tier_b.label_source}")
    print(f"           supports: {[t.value for t in LabelTask if t.value != LabelTask.PERI.value]}")
    print(f"\n  Validate: {tier_a.name}")
    print(f"           n={tier_a.n_paired}, {tier_a.centers}")
    print(f"           label: {tier_a.label_source}")
    print(f"           supports: {[t.value for t in LabelTask]}")
    print(
        "\n  The binding constraint is the label, not the imaging: no public cohort"
        "\n  pairs CD8 IHC spatial phenotype with MRI. Tier B cannot evaluate"
        "\n  desert-vs-excluded at all, because bulk expression has no way to encode"
        "\n  where the lymphocytes are."
    )
    print("\n  See docs/DESIGN.md for the full protocol.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="mri_cd8_til", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("inventory", help="public cohorts by label tier").set_defaults(func=cmd_inventory)
    sub.add_parser("intervals", help="CIs on the published AUCs").set_defaults(func=cmd_intervals)
    sub.add_parser("plan", help="the multi-centre plan").set_defaults(func=cmd_plan)
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
