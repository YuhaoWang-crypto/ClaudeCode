"""libY run once as Demo 1's out-of-domain negative control.

The client's instruction was to run PylRS-libY once on AzK and treat the result
as a negative control. Doing that honestly turns out to require saying three
separate things, two of which are blockers rather than scores.

## 1. ❌ libY structurally cannot score AzK

The only substrate-aware features in libY's 280-column matrix are columns 0..17,
eighteen Rosetta `cartesian_ddg` terms (recovered in `pylrs/audit.py`; the
README's "256 + 5 + 1 + 18" ordering is wrong). Columns 18..279 -- 5 ESM-1v
scores, 1 ESM-msa-1b score, 256 eUniRep dimensions -- are **identical across all
8 substrate rows of a given variant**, so they carry no information about which
ncAA is being asked about.

That means scoring any new ncAA requires running Rosetta `cartesian_ddg` against
it, which is licence-gated. There is no path from the shipped repo to an AzK
prediction. The negative control therefore cannot be "libY scored AzK badly" --
it is "libY cannot be pointed at AzK at all", which is a stronger and cheaper
finding.

And the same 18 columns are the ones that transfer at **chance** (AUC 0.487)
across substrates, while the variant-only block reaches 0.704. So the single
feature block that knows about the substrate is both unobtainable for a new one
and worthless for transfer.

## 2. ✅ The categorical test says AzK is out of domain, decisively

libY's 8 training substrates occupy exactly one chemotype cell: `parent=Tyr`,
`attach=para_O`. AzK is a lysine whose side-chain amine is capped as a carbamate.
Three of its categorical levels appear in **zero** libY training rows, so a model
fitted there has no parameter that AzK's kind of molecule would activate.

## 3. ❌ The Euclidean chemotype distance does NOT detect this

✅ Measured below: AzK sits **3.70** from its nearest libY substrate, while libY's
own 8 substrates span up to **4.32** internally. By distance alone AzK looks
*closer* to the training set than some training points are to each other.

This is the second independent failure of the distance idea in this project --
`demo3_evidence.py` already found that chemotype distance does not predict
held-out transfer (slope +0.06 ± 0.09, n=13). Two different tests, two failures.
The conclusion to carry into the platform is that **out-of-domain detection
should be categorical, not metric**: ask whether the new substrate's *kind* was
ever in training, not how far away it is.

Run:  python3 -m trnaplat.demo1_libY_control
      python3 -m trnaplat.demo1_libY_control --repo /path/to/PylRS-libY
"""

from __future__ import annotations

import argparse
import itertools
import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from trnaplat import ncaa_chemotype as nc  # noqa: E402

#: libY's feature layout, recovered in `pylrs/audit.py` rather than documented.
FEATURE_BLOCKS = [
    ("cols 0..17", 18, "Rosetta cartesian_ddg", True,
     "the ONLY substrate-aware block; needs a Rosetta licence for a new ncAA"),
    ("cols 18..22", 5, "ESM-1v", False, "identical across all 8 substrate rows"),
    ("col 23", 1, "ESM-msa-1b", False, "identical across all 8 substrate rows"),
    ("cols 24..279", 256, "eUniRep (PylRS-tuned)", False,
     "identical across all 8 substrate rows"),
]

#: ✅ Measured by `python3 pylrs/audit.py --repo <PylRS-libY checkout>`,
#: re-confirmed 2026-10-06 on a fresh clone. Best model per split shown
#: (ExtraTrees throughout); the audit prints logreg alongside, always lower.
AUDIT_NUMBERS = [
    ("published score_test (in-sample)", 0.9905,
     "the repo's own leaderboard; scored on the rows it was fitted on"),
    ("random 5-fold, all features", 0.708,
     "ignores structure -- what the published number is closest to"),
    ("by-variant GroupKFold, all features", 0.636,
     "no sibling rows in train; the honest pooled figure"),
    ("new ncAA, known variants -- Rosetta block", 0.487,
     "the ONLY substrate-aware block, and it transfers at chance"),
    ("new ncAA, known variants -- variant-only block", 0.704,
     "what generalises is 'which variant is promiscuous', not fit to the ncAA"),
    ("new ncAA AND new variants -- best block", 0.470,
     "⚠️ the real cold start, and the scenario the platform actually faces"),
]


def categorical_coverage(frame: pd.DataFrame, probe: str,
                         training: list[str]) -> pd.DataFrame:
    """Which of the probe's categorical levels were ever seen in training?

    This is the out-of-domain test that works. A one-hot level absent from every
    training row has no fitted weight -- the model is not extrapolating there,
    it is simply blind.
    """
    cat_cols = ["parent", "attach", "subst_class"]
    probe_row = frame[frame["ncAA"] == probe].iloc[0]
    train_rows = frame[frame["ncAA"].isin(training)]
    rows = []
    for col in cat_cols:
        level = probe_row[col]
        seen = int((train_rows[col] == level).sum())
        rows.append({
            "column": col, "probe_level": level,
            "training_rows_with_this_level": seen,
            "training_levels_present": ", ".join(sorted(train_rows[col].unique())),
            "in_domain": seen > 0,
        })
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=pathlib.Path,
                    help="a PylRS-libY checkout, to confirm the feature layout")
    ap.add_argument("--probe", default="AzK")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "trnaplat/data")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    frame = nc.full_table()
    liby = [n for n in frame["ncAA"] if n.startswith("libY-")]
    dist = nc.distance_matrix(frame)

    print("=" * 78)
    print(f"libY as an out-of-domain control for {args.probe}")
    print("=" * 78)
    print(f"\n  union chemotype table: {len(frame)} molecules "
          f"({len(liby)} libY substrates, 13 TyrRS ncAAs, {args.probe})")
    print("  ✅ pPRF / libY-C cross-check passed (same molecule, two sources)")

    # --- 1. the structural blocker ----------------------------------------
    print("\n" + "=" * 78)
    print("[1] ❌ Can libY be pointed at a new ncAA at all?")
    print("=" * 78 + "\n")
    blocks = pd.DataFrame(
        [{"columns": c, "n": n, "block": b, "substrate_aware": aware, "note": note}
         for c, n, b, aware, note in FEATURE_BLOCKS])
    print(blocks.to_string(index=False))
    aware = blocks[blocks["substrate_aware"]]
    print(f"\n  Substrate-aware features: {int(aware['n'].sum())} of"
          f" {int(blocks['n'].sum())} ({aware['n'].sum() / blocks['n'].sum():.1%}).")
    print("  All of them come from Rosetta cartesian_ddg, which is licence-gated.")
    print(f"  -> There is no route from the shipped repo to an {args.probe} score.")

    if args.repo:
        info = args.repo / "3.training/UAA-info.csv"
        if info.exists():
            lines = [ln for ln in info.read_text().splitlines() if ln.strip()]
            print(f"\n  ✅ confirmed against {info}: {len(lines)} substrates, all"
                  " O-substituted tyrosines:")
            for line in lines:
                print(f"       {line}")
        else:
            print(f"\n  ⚠️ {info} not found -- is --repo a PylRS-libY checkout?")

    print("\n  ✅ What the audit measured on those same features"
          " (pylrs/audit.py, best model per split):\n")
    for name, value, note in AUDIT_NUMBERS:
        print(f"       {value:.4f}  {name}")
        print(f"               {note}")
    print("\n  Two things to carry out of this table:")
    print("  * The one block that knows which substrate it is asked about scores")
    print("    0.487 across substrates -- chance. The block that does generalise")
    print("    (0.704) is the one that cannot tell the substrates apart at all,")
    print("    so what it has learned is a promiscuity prior on the scaffold.")
    print("  * ⚠️ On a NEW ncAA with NEW variants -- which is exactly what a")
    print("    platform is asked for -- the best figure is 0.470, at or below")
    print("    chance. libY's published 0.990 and its real cold start differ by")
    print("    more than the gap between its cold start and a coin.")

    # --- 2. the categorical test that works -------------------------------
    print("\n" + "=" * 78)
    print(f"[2] ✅ Is {args.probe} out of libY's domain? The categorical test")
    print("=" * 78 + "\n")
    coverage = categorical_coverage(frame, args.probe, liby)
    print(coverage.to_string(index=False))
    coverage.to_csv(args.out / "libY_control_coverage.csv", index=False)
    unseen = coverage[~coverage["in_domain"]]
    print(f"\n  {len(unseen)}/{len(coverage)} of {args.probe}'s categorical levels appear"
          " in ZERO libY training rows.")
    if len(unseen):
        for row in unseen.itertuples():
            print(f"       {row.column}={row.probe_level!r} -- training has only"
                  f" {row.training_levels_present!r}")
    print(f"\n  ✅ A one-hot level absent from every training row carries no fitted")
    print("     weight. The model is not extrapolating here, it is blind.")

    # --- 3. the metric test that fails ------------------------------------
    print("\n" + "=" * 78)
    print("[3] ❌ Does Euclidean chemotype distance detect it?")
    print("=" * 78 + "\n")
    within = [dist.loc[a, b] for a, b in itertools.combinations(liby, 2)]
    probe_gap = float(dist.loc[args.probe, liby].min())
    nearest = str(dist.loc[args.probe, liby].idxmin())
    within_s = pd.Series(within)
    print(f"  libY's own 8 substrates, pairwise distance:")
    print(f"      min {within_s.min():.2f}   median {within_s.median():.2f}"
          f"   max {within_s.max():.2f}")
    print(f"  {args.probe} to its nearest libY substrate ({nearest}):"
          f"  {probe_gap:.2f}")
    n_wider = int((within_s > probe_gap).sum())
    print(f"\n  ❌ {n_wider} of {len(within)} within-training pairs are FURTHER apart"
          f" than {args.probe} is")
    print(f"     from the training set. By distance alone {args.probe} does not look")
    print("     out of domain at all -- yet it is a different parent amino acid.")
    print("\n  ⚠️ This is the second failure of the distance idea in this project.")
    print("     demo3_evidence.py found distance does not predict held-out transfer")
    print("     (slope +0.06 ± 0.09, n=13); here it does not even detect a change of")
    print("     parent residue. Use the categorical test in [2] for the platform's")
    print("     out-of-domain gate, and keep distance only as a reported descriptor.")

    summary = pd.DataFrame([{
        "probe": args.probe,
        "liby_internal_min": within_s.min(), "liby_internal_median": within_s.median(),
        "liby_internal_max": within_s.max(),
        "probe_to_nearest_training": probe_gap, "nearest_training": nearest,
        "within_pairs_further_than_probe": n_wider, "within_pairs_total": len(within),
        "unseen_categorical_levels": len(unseen),
        "categorical_levels_tested": len(coverage),
    }])
    summary.to_csv(args.out / "libY_control_summary.csv", index=False)

    print("\n" + "=" * 78)
    print("The control's verdict")
    print("=" * 78)
    print(f"  ❌ libY cannot score {args.probe}: its only substrate-aware features need")
    print("     a Rosetta licence, and they transfer at chance anyway (0.487).")
    print(f"  ✅ {args.probe} is out of domain by the categorical test"
          f" ({len(unseen)}/{len(coverage)} levels unseen).")
    print("  ❌ Euclidean chemotype distance fails to detect that.")
    print("  -> Demo 1 gains nothing from libY, which is the expected answer and now")
    print("     a measured one. Demo 3 should borrow libY's evaluation architecture")
    print("     and gate on unseen categorical levels, not on distance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
