"""McTyrRS as a second screenable scaffold, seeded by transfer from MjTyrRS.

Why this is worth doing: the cross-species comparison (`crossspecies.py`) found
that on the one ncAA measured on both scaffolds, the *Methanosaeta concilii*
enzyme discriminates better than the *Methanococcus jannaschii* one
(F+/F- 3.69 vs 2.79), and that the two species converged independently on the
same two substitutions. McTyrRS is therefore not a curiosity, it is a real
second starting point -- but it has only a handful of measured variants of its
own, so on its own there is nothing to train on.

Transfer fixes that. The 76 MjTyrRS clones are rewritten into McTyrRS numbering
through the verified pairwise alignment (54.1% identity, every engineered site
maps 1:1), giving McTyrRS a 76-clone prior it did not have.

⚠️ What transfer does and does not license
------------------------------------------
Rewriting coordinates is not the same as transferring function. A substitution
that works in MjTyrRS may not work at the equivalent McTyrRS position: the two
pockets differ at 46% of positions, including E107/D111 and F108/Y112 inside the
site. The transferred clones are a **prior over where and what to mutate**, not
predicted McTyrRS hits. `report_transfer_risk` quantifies the part that is
checkable -- how often the native residue differs between the scaffolds at the
positions being transferred.

✅ Sheet 3 also contains the only measured non-hits anywhere in this dataset.
Mut4, Mut5 and Mut7 were built and the paper states Mut6 "is selected because
other ones are not good". They are labelled `weak` rather than 0 because no
number is given for them, but they are measured rejections, not invented decoys,
and they are the first real negatives this project has had.

Run:  python3 pylrs/mc_scaffold.py
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from crossspecies import MC_VARIANTS, align_positions, scaffolds  # noqa: E402

#: Outcome per McTyrRS variant, read from sheet 3's own wording.
MC_OUTCOME = {
    "Mc Mut4": "weak",      # "other ones are not good"
    "Mc Mut5": "weak",
    "Mc Mut7": "weak",
    "Mc Mut6": "hit",       # F+/F- = 1.11 -- selected, but barely discriminating
    "Mc Mut6+RM": "hit",    # F+/F- = 3.69 -- the best McTyrRS variant reported
}


def transfer(mj_variants: pd.DataFrame, mapping: dict[int, int],
             mj: str, mc: str) -> pd.DataFrame:
    """Rewrite MjTyrRS clone mutations into McTyrRS numbering.

    The substituted residue is carried across; the 'from' letter becomes
    McTyrRS's own native residue, which is what makes the result a valid
    McTyrRS variant rather than a mislabelled Mj one.
    """
    reverse = {v: k for k, v in mapping.items()}       # Mj position -> Mc position
    rows = []
    for row in mj_variants.itertuples():
        out, dropped, divergent = [], [], 0
        for mut in str(row.mutations).split("/"):
            if not mut:
                continue
            pos, new = int(mut[1:-1]), mut[-1]
            mc_pos = reverse.get(pos)
            if mc_pos is None:
                dropped.append(mut)
                continue
            if mc[mc_pos - 1] != mj[pos - 1]:
                divergent += 1
            if mc[mc_pos - 1] == new:
                continue      # already the McTyrRS native residue: not a mutation
            out.append(f"{mc[mc_pos - 1]}{mc_pos}{new}")
        if not out:
            continue
        rows.append({
            "clone": f"{row.clone}->Mc", "ncAA": row.ncAA,
            "mutations": "/".join(sorted(out, key=lambda m: int(m[1:-1]))),
            "n_mut": len(out), "origin": "transferred from MjTyrRS",
            "dropped": "/".join(dropped),
            "divergent_sites": divergent,
            "label": 1, "outcome": "hit-on-MjTyrRS",
        })
    return pd.DataFrame(rows)


def native_mc_variants(mapping, mj, mc) -> pd.DataFrame:
    rows = []
    for variant in MC_VARIANTS:
        outcome = MC_OUTCOME[variant["name"]]
        rows.append({
            "clone": variant["name"], "ncAA": "pAzF",
            "mutations": "/".join(variant["mutations"]),
            "n_mut": len(variant["mutations"]),
            "origin": "measured on McTyrRS (sheet 3)",
            "dropped": "", "divergent_sites": 0,
            "label": 1 if outcome == "hit" else 0,
            "outcome": outcome,
            "F_plus": variant.get("f_plus"), "F_minus": variant.get("f_minus"),
        })
    return pd.DataFrame(rows)


def report_transfer_risk(mj_variants, mapping, mj, mc) -> pd.DataFrame:
    """Per position: how often the two scaffolds disagree on the native residue."""
    reverse = {v: k for k, v in mapping.items()}
    counts: dict[int, int] = {}
    for muts in mj_variants["mutations"]:
        for mut in str(muts).split("/"):
            if mut:
                counts[int(mut[1:-1])] = counts.get(int(mut[1:-1]), 0) + 1
    rows = []
    for mj_pos, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        mc_pos = reverse.get(mj_pos)
        rows.append({
            "Mj_pos": mj_pos, "Mj_native": mj[mj_pos - 1],
            "Mc_pos": mc_pos if mc_pos else None,
            "Mc_native": mc[mc_pos - 1] if mc_pos else "-",
            "same": (mc_pos is not None and mc[mc_pos - 1] == mj[mj_pos - 1]),
            "clones_touching": n,
        })
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaigns", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/literature_campaigns.csv"))
    ap.add_argument("--variants", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/tyrrs_variants.csv"))
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("pylrs/data"))
    args = ap.parse_args(argv)

    seqs = scaffolds(args.campaigns)
    mj, mc = seqs["MjTyrRS"], seqs["McTyrRS"]
    mapping = align_positions(mc, mj)        # Mc -> Mj

    mj_pos = pd.read_csv(args.variants).query("label == 1")
    mj_pos = mj_pos[~mj_pos.duplicated(subset=["mutations"])]

    transferred = transfer(mj_pos, mapping, mj, mc)
    native = native_mc_variants(mapping, mj, mc)
    table = pd.concat([native, transferred], ignore_index=True)

    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / "mc_variants.csv"
    table.to_csv(target, index=False)

    print(f"McTyrRS {len(mc)} aa | MjTyrRS {len(mj)} aa | {len(mapping)} aligned positions")
    print(f"transferred {len(transferred)} of {len(mj_pos)} MjTyrRS clones")
    print(f"native McTyrRS variants from sheet 3: {len(native)}")
    print(f"-> {target}\n")

    print("=" * 74)
    print("Measured McTyrRS variants -- the only measured non-hits in this dataset")
    print("=" * 74)
    cols = ["clone", "mutations", "outcome", "F_plus", "F_minus"]
    print(native[cols].to_string(index=False))
    n_neg = int((native["label"] == 0).sum())
    print(f"\n  {n_neg} measured rejections (Mut4/5/7). ⚠️ The paper gives no numbers")
    print("  for them, only that Mut6 was chosen because the others were not good,")
    print("  so they are labelled `weak`, not a hard zero. They are still the first")
    print("  negatives here that someone actually built.")

    print("\n" + "=" * 74)
    print("Transfer risk: does McTyrRS have the same native residue to mutate from?")
    print("=" * 74)
    risk = report_transfer_risk(mj_pos, mapping, mj, mc)
    print(risk.to_string(index=False))
    same = risk[risk["same"]]["clones_touching"].sum()
    total = risk["clones_touching"].sum()
    print(f"\n  {same}/{total} ({same / total:.0%}) of transferred substitutions start")
    print("  from the same native residue in both scaffolds.")
    diverging = risk[~risk["same"]]
    if len(diverging):
        print("  ⚠️ Diverging sites, where the transfer changes meaning:")
        for row in diverging.itertuples():
            print(f"     Mj {row.Mj_native}{row.Mj_pos} vs Mc {row.Mc_native}{row.Mc_pos}"
                  f"  ({row.clones_touching} clones)")

    drops = transferred[transferred["dropped"] != ""]
    if len(drops):
        print(f"\n  {len(drops)} clones lost a substitution with no McTyrRS equivalent:")
        for row in drops.head(5).itertuples():
            print(f"     {row.clone}: dropped {row.dropped}")

    print("\n" + "=" * 74)
    print("What this scaffold is now ready for")
    print("=" * 74)
    pazf = table[table["ncAA"] == "pAzF"]
    print(f"  pAzF rows on McTyrRS: {len(pazf)} "
          f"({int((pazf['label'] == 1).sum())} hits, {int((pazf['label'] == 0).sum())} rejected)")
    print("  Screen it the same way as MjTyrRS:")
    print("    python3 pylrs/screen_pazf.py --variants pylrs/data/mc_variants.csv \\")
    print("        --scaffold McTyrRS --positions 33,111,162,163,166")
    print("\n  ⚠️ Calibrate before trusting it. On MjTyrRS the same model recovered")
    print("  only 4 of 7 known pAzF clones; on a transferred prior expect worse, and")
    print("  read pazf_calibration.csv for the Mc run before ordering anything.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
