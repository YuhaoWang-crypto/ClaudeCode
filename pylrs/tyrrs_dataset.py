"""A labelled MjTyrRS variant x ncAA dataset, mined from sheet 2 of the curation file.

Where the data actually lives
-----------------------------
The campaign rows of `dataset_from_paper_v2.0.xlsx` hold prose, not variants. The
clone tables are **embedded images** (21 of them across the 3 sheets), anchored
below a "figureN" label cell. So this module carries a hand transcription of
those tables, with the source image recorded per campaign. ⚠️ Transcribed by eye
from published figures: treat `provenance` as the thing to re-check before anyone
orders DNA, and re-read the image rather than trusting this file blindly.

Two campaigns are deliberately **not** transcribed:

* 2.9 (OCF3Phe, image13) is a dense 14-clone sequence alignment. Reading 14x9
  residues off an alignment image is error-prone in a way a clone table is not,
  and a wrong residue here is silent. Left out rather than guessed.
* 2.15 (Nal, image18) reports only IC50 values, with no clone sequences.

Numbering is checked, not assumed: `literature.verify_numbering` confirms all 17
cited positions (Y32, L65, A67, H70, Y102, V103, E107, F108, Q109, L110, Y114,
Q155, D158, I159, H160, Y161, L162) against the 306-residue WT sequence. The raw
cell prefixes the sequence with "WT:", which if not stripped shifts every residue
by two and quietly invalidates every mutation below.

⚠️ The negatives are not measured
---------------------------------
These papers report **survivors only**. There is no published set of variants
shown not to work, so a classifier trained here is positive-unlabelled. `decoys()`
samples random members of each campaign's own randomised library as presumed
negatives, which is defensible — the libraries carry ~10^9 members and a handful
survive selection, so a random member is almost certainly non-functional — but it
is an assumption, not an observation, and it is the single biggest caveat on any
number this dataset produces.

Run:  python3 pylrs/tyrrs_dataset.py --campaigns pylrs/data/literature_campaigns.csv
"""

from __future__ import annotations

import argparse
import pathlib

import numpy as np
import pandas as pd

AA20 = "ACDEFGHIKLMNPQRSTVWY"

#: Each campaign: the ncAA, the positions its clone table reports, and the clones.
#: `clones` maps a clone name to the residues at `positions`, in order; "." means
#: the published table marked that residue unchanged from wild type. `extra` holds
#: substitutions reported in a footnote rather than a column. `freq` is the clone
#: count behind a unique sequence where the paper gives it.
CAMPAIGNS = [
    dict(record="2.1", ncAA="pBpa", ncAA_full="p-benzoyl-L-phenylalanine",
         provenance="image7.png (sheet2 figure1)",
         positions=[32, 107, 158, 159, 162],
         clones={"MjBpaRS-1": "GSTSL", "MjBpaRS-2": "ALTAL", "MjBpaRS-3": "GPTSL",
                 "MjBpaRS-4": "WSTSR", "MjBpaRS-5": "ASTSA", "MjBpaRS-6": "GPGTL"}),

    dict(record="2.3", ncAA="pAzF", ncAA_full="p-azido-L-phenylalanine",
         provenance="image8.png (sheet2 Figure2, 'Table 1 AzPheRS Sequences')",
         positions=[32, 107, 158, 159, 162],
         clones={"AzPheRS-1": "TNPLQ", "AzPheRS-2": "TSPSQ", "AzPheRS-3": "TSPLQ",
                 "AzPheRS-4": "LTPVQ", "AzPheRS-5": "ARVID",
                 "AzPheRS-6": "GTTYL", "AzPheRS-7": "LPQIS"},
         extra={"AzPheRS-6": [(160, "Y")], "AzPheRS-7": [(161, "S")]}),

    dict(record="2.4a", ncAA="pIF", ncAA_full="p-isopropyl-phenylalanine",
         provenance="image9.png (sheet2 figure3)",
         positions=[32, 102, 103, 107, 158, 159, 162],
         clones={"pIF-RS": "GCAPGY."}),

    dict(record="2.4b", ncAA="pAF", ncAA_full="p-amino-phenylalanine",
         provenance="image9.png (sheet2 figure3)",
         positions=[32, 102, 103, 107, 158, 159, 162],
         clones={"pAF-RS": "T..TPLA"}),

    dict(record="2.4c", ncAA="OAY", ncAA_full="O-allyl-tyrosine",
         provenance="image10.png (sheet2 figure4); OAY-RS(1) agrees with figure3",
         positions=[32, 107, 158, 159, 162],
         clones={"OAY-RS(1)": "STTYA", "OAY-RS(3)": "THQTE",
                 "OAY-RS(5)": "PMNTG", "OAY-RS(B)": ".ACA."}),

    dict(record="2.5", ncAA="pPRF", ncAA_full="para-propargyloxyphenylalanine",
         provenance="image11.png (sheet2 figure 5)",
         positions=[32, 107, 158, 159, 162],
         clones={"pPR-MjRS-1": "APAIA", "pPR-MjRS-2": "AKAIA", "pPR-MjRS-3": "ARAIP",
                 "pPR-MjRS-4": "HAAIP", "pPR-MjRS-5": "SQAIA", "pPR-MjRS-6": "TSLHP",
                 "pPR-MjRS-7": "AQPGT", "pPR-MjRS-8": "APSLH"},
         extra={"pPR-MjRS-1": [(110, "F")]}),

    dict(record="2.7", ncAA="HQ-Ala", ncAA_full="2-amino-3-(8-hydroxyquinolin-3-yl)propanoic acid",
         provenance="image12.png (sheet2 Figure6)",
         positions=[32, 65, 70, 108, 109, 155, 158, 159, 162],
         clones={"HQ-3D1": "HVMFQLGSL", "HQ-3D2": "LISRENGIP", "HQ-3D3": "EVTREQSSS",
                 "HQ-3D4": "VMTREQSSL", "HQ-3D5": "MCTREQSTH", "HQ-3D6": "LITREQSAR",
                 "HQ-3D7": "LQTREQSYL", "HQ-3D8": "LTNFQQGTL"},
         freq={"HQ-3D1": 1, "HQ-3D2": 1, "HQ-3D3": 1, "HQ-3D4": 4,
               "HQ-3D5": 3, "HQ-3D6": 2, "HQ-3D7": 1, "HQ-3D8": 2}),

    dict(record="2.12a", ncAA="BipAla", ncAA_full="biphenylalanine",
         provenance="image16.png (sheet2 Figure8, Table 1)",
         positions=[32, 65, 70, 108, 109, 155, 158, 162],
         clones={"BipAlaRS1": "HHHWMQGK", "BipAlaRS2": "GVHADQSR", "BipAlaRS3": "HVHSKQGS",
                 "BipAlaRS4": "ASHREQSS", "BipAlaRS5": "HVHRPQGL", "BipAlaRS6": "AVHSHQGE",
                 "BipAlaRS7": "LTSFQQVH"},
         freq={"BipAlaRS1": 1, "BipAlaRS2": 8, "BipAlaRS3": 4, "BipAlaRS4": 3,
               "BipAlaRS5": 2, "BipAlaRS6": 1, "BipAlaRS7": 1}),

    dict(record="2.12b", ncAA="BpyAla", ncAA_full="(2,2'-bipyridin-5-yl)alanine",
         provenance="image16.png (sheet2 Figure8, Table 2)",
         positions=[32, 65, 70, 108, 109, 155, 158, 159, 162],
         clones={"BpyAlaRS1": "GYAFQEGWS", "BpyAlaRS2": "EHHWMQGHH"},
         freq={"BpyAlaRS1": 10, "BpyAlaRS2": 1}),

    dict(record="2.13", ncAA="pBoroPhe", ncAA_full="p-boronophenylalanine",
         provenance="image14.png (sheet2 Figure9, Table S1)",
         positions=[32, 65, 70, 155, 158, 162],
         clones={"1BF6": "SAHQSE", "1BF9": "SAHQSE", "1BE3": "GAHQAE",
                 "1BF10": "SAMQSE", "1BF12": "SAMQSE",
                 "1BG10": "SAMQSE", "1BG11": "SAMQSE"}),

    dict(record="2.14", ncAA="pCMF", ncAA_full="p-carboxymethyl-L-phenylalanine",
         provenance="image17.png (sheet2 Figure10, Table 1)",
         positions=[32, 65, 108, 109, 158, 162],
         clones={"pCMFRS#1": "SAKHGK", "pCMFRS#2": "SAAQGK", "pCMFRS#3": "ASRNSH",
                 "pCMFRS#4": "SANYGK", "pCMFRS#5": "AARQKA"},
         freq={"pCMFRS#1": 11, "pCMFRS#2": 11, "pCMFRS#3": 2,
               "pCMFRS#4": 1, "pCMFRS#5": 1}),

    # Transcribed from a sequence alignment rather than a clone table, so the
    # library consensus rows are carried too and checked against the wild type by
    # `verify_alignment_rows` -- if the transcription or the numbering were off,
    # that check fails instead of silently producing wrong mutations.
    dict(record="2.9", ncAA="OCF3Phe", ncAA_full="p-trifluoromethoxy-L-phenylalanine",
         provenance="image13.png (sheet2 Figure 7, Table S2 alignment)",
         positions=[26, 32, 64, 65, 70, 108, 109, 155, 158, 159, 162],
         clones={
             "OCF3PHE_A6":  "KVIAHQWQAIK", "OCF3PHE_B6":  "KAIAHKWQGIV",
             "OCF3PHE_B7":  "KVILHAWQGIQ", "OCF3PHE_B10": "KAIAHWMQGNL",
             "OCF3PHE_C2":  "KVLGHEWQGIV", "OCF3PHE_D5":  "KVIHHEPQSIS",
             "OCF3PHE_D9":  "KLIPHWMQGAL", "OCF3PHE_E7":  "KVISHTQQAIV",
             "OCF3PHE_F6":  "KAISHQAQAIY", "OCF3PHE_F7":  "KIITHRWQAIS",
             "OCF3PHE_F8":  "KVIQHRESSVH", "OCF3PHE_G2":  "KHIANWMQGAL",
             "OCF3PHE_G5":  "KVITHLGQSIS", "OCF3PHE_H4":  "IVIGHHYQAIH"}),

    dict(record="2.16", ncAA="2-NPA", ncAA_full="2-nitrophenylalanine",
         provenance="image19.png (sheet2 FIgure12, Table S1)",
         positions=[32, 65, 67, 70, 108, 109, 114, 158, 159, 162],
         clones={"2NPA-1": "GHGGLSSTYD", "2NPA-2": "GHGGQLNACD", "2NPA-3": "GHGGYLSAHD",
                 "2NPA-4": "GHGGQLNTYE", "2NPA-5": "GHGGQFGAYD", "2NPA-6": "GHGGECASVE"}),
]


#: The "Library" consensus rows of the OCF3Phe alignment (image13), as three
#: blocks starting at residues 1, 61 and 121. 'X' marks a randomised column.
#: Every non-X residue must equal the wild type, which is what makes the
#: transcription of that alignment checkable rather than trusted.
ALIGNMENT_LIBRARY_ROWS = {
    1: "MDEFEMIKRNTSEIISEEELREVLKKDEKSAXIGFEPSGKIHLGHYLQIKKMIDLQNAGF",
    61: "DIIIXLADLXAYLNQKGELDEIRKIGDYNKKVFEAMGLKAKYVYGSEXXLDKDYTLNVYR",
    121: "LALKTTLKRARRSMELIAREDENPKVAEVIYPIMXVNXXHYXGVDVAVGGMEQRKIHMLA",
}


def verify_alignment_rows(wt: str) -> list[str]:
    """Positions where the transcribed library rows disagree with the wild type."""
    bad = []
    for start, row in ALIGNMENT_LIBRARY_ROWS.items():
        for offset, aa in enumerate(row):
            pos = start + offset
            if aa != "X" and wt[pos - 1] != aa:
                bad.append(f"{pos}: alignment {aa} vs wild type {wt[pos - 1]}")
    return bad


def wild_type(campaigns_csv: pathlib.Path) -> str:
    """The 306-residue MjTyrRS wild-type sequence, numbering-verified."""
    from literature import verify_numbering  # noqa: F401  (same directory)

    table = pd.read_csv(campaigns_csv)
    seq = table.loc[table["scaffold"] == "MjTyrRS", "sequence"].iloc[0]
    bad = verify_numbering(seq)
    if bad:
        raise SystemExit(f"MjTyrRS numbering check failed: {bad}")
    return seq


def mutations_of(campaign: dict, clone: str, wt: str) -> list[str]:
    """Substitutions of one clone, as 'Y32T'-style labels against the wild type."""
    muts = []
    for pos, aa in zip(campaign["positions"], campaign["clones"][clone]):
        native = wt[pos - 1]
        if aa == "." or aa == native:
            continue   # the table marked this position unchanged
        muts.append(f"{native}{pos}{aa}")
    for pos, aa in campaign.get("extra", {}).get(clone, []):
        muts.append(f"{wt[pos - 1]}{pos}{aa}")
    return sorted(muts, key=lambda m: int(m[1:-1]))


def build(campaigns_csv: pathlib.Path) -> tuple[pd.DataFrame, str]:
    """Return the labelled positive set plus the wild-type sequence."""
    wt = wild_type(campaigns_csv)
    rows = []
    for campaign in CAMPAIGNS:
        for clone in campaign["clones"]:
            muts = mutations_of(campaign, clone, wt)
            rows.append({
                "clone": clone,
                "record": campaign["record"],
                "ncAA": campaign["ncAA"],
                "ncAA_full": campaign["ncAA_full"],
                "mutations": "/".join(muts),
                "n_mut": len(muts),
                "positions": ",".join(str(p) for p in campaign["positions"]),
                "frequency": campaign.get("freq", {}).get(clone, pd.NA),
                "label": 1,
                "provenance": campaign["provenance"],
            })
    return pd.DataFrame(rows), wt


def decoys(positives: pd.DataFrame, wt: str, per_positive: int = 20,
           seed: int = 0) -> pd.DataFrame:
    """Presumed-negative library members, sampled from each campaign's own library.

    ⚠️ Assumed, not measured. Each decoy randomises the same positions the campaign
    randomised, which is what the published library actually contained; a random
    member of a ~10^9 library is almost certainly non-functional, but no one
    assayed these. Decoys that collide with a known positive are rejected.
    """
    rng = np.random.default_rng(seed)
    known = set(positives["mutations"])
    rows = []
    for record, group in positives.groupby("record"):
        positions = [int(p) for p in group["positions"].iloc[0].split(",")]
        ncaa = group["ncAA"].iloc[0]
        wanted = per_positive * len(group)
        seen: set[str] = set()
        attempts = 0
        while len(seen) < wanted and attempts < wanted * 50:
            attempts += 1
            muts = []
            for pos in positions:
                aa = AA20[rng.integers(len(AA20))]
                if aa != wt[pos - 1]:
                    muts.append(f"{wt[pos - 1]}{pos}{aa}")
            key = "/".join(muts)
            if not muts or key in known or key in seen:
                continue
            seen.add(key)
            rows.append({
                "clone": f"decoy_{record}_{len(seen)}",
                "record": record, "ncAA": ncaa, "ncAA_full": "",
                "mutations": key, "n_mut": len(muts),
                "positions": ",".join(str(p) for p in positions),
                "frequency": pd.NA, "label": 0,
                "provenance": "sampled from the campaign's randomised positions",
            })
    return pd.DataFrame(rows)


def apply_mutations(wt: str, mutations: str) -> str:
    """Build the full variant sequence, validating each stated wild-type residue."""
    seq = list(wt)
    for mut in filter(None, mutations.split("/")):
        native, pos, new = mut[0], int(mut[1:-1]), mut[-1]
        if seq[pos - 1] != native:
            raise ValueError(f"{mut}: wild type has {seq[pos - 1]} at {pos}")
        seq[pos - 1] = new
    return "".join(seq)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaigns", type=pathlib.Path,
                    default=pathlib.Path("pylrs/data/literature_campaigns.csv"))
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("pylrs/data"))
    ap.add_argument("--decoys-per-positive", type=int, default=20)
    args = ap.parse_args(argv)

    pos, wt = build(args.campaigns)

    bad = verify_alignment_rows(wt)
    if bad:
        raise SystemExit("the transcribed OCF3Phe alignment disagrees with the wild "
                         f"type at: {bad}")
    fixed = sum(len(r) - r.count("X") for r in ALIGNMENT_LIBRARY_ROWS.values())
    print(f"✅ OCF3Phe alignment library rows match the wild type at all {fixed} "
          "non-randomised positions")
    neg = decoys(pos, wt, per_positive=args.decoys_per_positive)
    full = pd.concat([pos, neg], ignore_index=True)

    # Every mutation must apply cleanly to the wild type, or the numbering is wrong.
    for mutations in full["mutations"]:
        apply_mutations(wt, mutations)

    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / "tyrrs_variants.csv"
    full.to_csv(target, index=False)

    print(f"MjTyrRS wild type: {len(wt)} aa, numbering verified")
    print(f"positives: {len(pos)} clones over {pos['ncAA'].nunique()} ncAAs "
          f"({pos['mutations'].nunique()} unique mutation sets)")
    print(f"decoys:    {len(neg)} presumed negatives (⚠️ assumed, not measured)")
    print(f"-> {target}\n")

    summary = (pos.groupby(["record", "ncAA"])
               .agg(clones=("clone", "size"),
                    unique=("mutations", "nunique"),
                    mean_muts=("n_mut", "mean"))
               .reset_index())
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    sites = {}
    for mutations in pos["mutations"]:
        for mut in mutations.split("/"):
            sites[int(mut[1:-1])] = sites.get(int(mut[1:-1]), 0) + 1
    print("\nmutated positions across all positives (position: count):")
    print("  " + ", ".join(f"{p}:{n}" for p, n in sorted(sites.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
