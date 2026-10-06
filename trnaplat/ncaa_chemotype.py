"""Hand-curated chemotype descriptors for the 13 ncAAs in the MjTyrRS dataset.

Why hand-curated rather than computed from SMILES: the properties that decide
whether one synthetase pocket accepts two different ncAAs are *where* the
substituent sits and *what shape* it is, and a generic cheminformatics
fingerprint buries those in hundreds of uninformative bits on a 13-row dataset.
Thirteen rows cannot train a fingerprint; they can carry eight interpretable
columns.

⚠️ Every value below was read off the chemical name in
`pylrs/data/literature_campaigns.csv` by hand. They are structural bookkeeping,
not measurements, and two of them are judgement calls flagged in `NOTES`.

The descriptors exist to answer one question: **how far is a new ncAA from the
nearest one we have data for?** That distance is what sets the evidence tier in
`demo3_evidence.py`, and it is the reason the platform can say "this prediction
is extrapolating" instead of returning a confident number for anything.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

#: parent      the canonical side chain the ncAA is built on
#: attach      where the substituent sits relative to the ring attachment
#: subst_class coarse chemistry of the substituent
#: n_heavy     heavy atoms added versus the parent side chain
#: n_rings     aromatic/heteroaromatic rings added versus the parent
#: hbd/hba     hydrogen-bond donor / acceptor introduced by the substituent
#: charge      formal charge of the substituent at pH 7
#: chelator    can coordinate a transition metal
DESCRIPTORS = [
    # ncAA,      parent, attach,     subst_class,      n_heavy, n_rings, hbd, hba, charge, chelator
    ("pBpa",      "Phe", "para",     "aryl_ketone",          7,       1,   0,   1,      0,        0),
    ("pAzF",      "Phe", "para",     "azide",                3,       0,   0,   1,      0,        0),
    ("pIF",       "Phe", "para",     "alkyl",                3,       0,   0,   0,      0,        0),
    ("pAF",       "Phe", "para",     "amine",                1,       0,   1,   1,      0,        0),
    ("OAY",       "Tyr", "para_O",   "ether_alkene",         3,       0,   0,   0,      0,        0),
    ("pPRF",      "Tyr", "para_O",   "ether_alkyne",         3,       0,   0,   0,      0,        0),
    ("HQ-Ala",    "Ala", "ring_fused", "hydroxyquinoline",   10,       2,   1,   2,      0,        1),
    ("BipAla",    "Phe", "para",     "aryl",                 6,       1,   0,   0,      0,        0),
    ("BpyAla",    "Ala", "ring_fused", "bipyridine",         10,       2,   0,   2,      0,        1),
    ("pBoroPhe",  "Phe", "para",     "boronate",             3,       0,   2,   2,      0,        0),
    ("pCMF",      "Phe", "para",     "carboxylate",          3,       0,   0,   2,     -1,        0),
    # parent Tyr already carries the para oxygen, so the substituent is CF3 = 4
    # heavy atoms, not OCF3 = 5. Caught by the libY-C/pPRF cross-check below.
    ("OCF3Phe",   "Tyr", "para_O",   "fluoroalkyl_ether",    4,       0,   0,   1,      0,        0),
    ("2-NPA",     "Phe", "ortho",    "nitro",                3,       0,   0,   2,      0,        0),
]

COLUMNS = ["ncAA", "parent", "attach", "subst_class", "n_heavy", "n_rings",
           "hbd", "hba", "charge", "chelator"]

#: PylRS-libY's 8 training substrates, from `3.training/UAA-info.csv` as recovered
#: in `pylrs/audit.py`. Encoded in the SAME columns so distances are comparable.
#:
#: ✅ Consistency check built in: substrate C (O-propargyl-Tyr) is the same
#: molecule as `pPRF` in DESCRIPTORS above, encoded independently. The two rows
#: agree, which is the only cross-check this hand-curated table has.
LIBY_SUBSTRATES = [
    ("libY-A OMe-Tyr",        "Tyr", "para_O", "ether_alkyl",        1, 0, 0, 0, 0, 0),
    ("libY-B OEt-Tyr",        "Tyr", "para_O", "ether_alkyl",        2, 0, 0, 0, 0, 0),
    ("libY-C OPropargyl-Tyr", "Tyr", "para_O", "ether_alkyne",       3, 0, 0, 0, 0, 0),
    ("libY-D OCyclohexyl-Tyr", "Tyr", "para_O", "ether_cycloalkyl",  6, 0, 0, 0, 0, 0),
    ("libY-E OPhenyl-Tyr",    "Tyr", "para_O", "ether_aryl",         6, 1, 0, 0, 0, 0),
    ("libY-F OBenzyl-Tyr",    "Tyr", "para_O", "ether_aryl",         7, 1, 0, 0, 0, 0),
    ("libY-G OCyclooctenyl-Tyr", "Tyr", "para_O", "ether_cycloalkene", 8, 0, 0, 0, 0, 0),
    ("libY-H ONitrobenzyl-Tyr", "Tyr", "para_O", "ether_aryl",      10, 1, 0, 2, 0, 0),
]

#: AzK -- Demo 1's substrate, and the out-of-domain probe for libY.
#:
#: ⚠️ Encoded as N-epsilon-((2-azidoethoxy)carbonyl)-L-lysine: a lysine whose
#: side-chain amine is capped as a carbamate carrying an azidoethyl group. The
#: 8 added heavy atoms are C(=O) + O + CH2 + CH2 + N3. If the client's validated
#: construct uses a different AzK regioisomer, change this row before quoting any
#: distance computed from it.
AZK = [
    ("AzK", "Lys", "Ne_carbamate", "carbamate_azide", 8, 0, 1, 3, 0, 0),
]

#: The two entries where the name does not determine the descriptor on its own.
NOTES = {
    "pIF": "the dataset's `ncAA_full` reads p-isopropyl-phenylalanine, while "
           "campaign 2.11 uses the same abbreviation for p-IODO-phenylalanine. "
           "Encoded here as isopropyl, following `ncAA_full`. ⚠️ An iodo "
           "substituent would change n_heavy to 1 and add a halogen column.",
    "HQ-Ala": "8-hydroxyquinolin-3-yl is fused rather than para-substituted, so "
              "`attach` is ring_fused and `parent` is Ala. ⚠️ Encoding it this "
              "way pairs it with BpyAla, the other fused chelator, rather than "
              "with the para-substituted phenylalanines -- which is chemically "
              "right but means its neighbour is one single ncAA, not a family.",
}

#: ncAA groups defined by **what a product needs**, not by chemical class.
#:
#: The multi-handle carrier that Demo 1 targets needs two handles that react
#: with different partners and do not react with each other, so two payloads can
#: be placed on one protein. That requirement, not "all Phe derivatives", is
#: what should define a model's training group -- and `demo3_evidence.py`
#: measures whether that is true.
PRODUCT_GROUPS = {
    "dual_orthogonal_click": {
        "members": ["pAzF", "pPRF"],
        "need": "two mutually bio-orthogonal handles on one carrier: azide "
                "(SPAAC/CuAAC) plus terminal alkyne (CuAAC/thiol-yne). Lets two "
                "different payloads be placed site-specifically on one protein.",
    },
    "photocrosslink_plus_click": {
        "members": ["pBpa", "pAzF"],
        "need": "map an interaction and then conjugate at a second site: "
                "benzophenone photo-crosslinking plus an azide handle.",
    },
    "metal_chelation": {
        "members": ["BpyAla", "HQ-Ala"],
        "need": "a genetically encoded metal site, for radiometal chelation or "
                "artificial metalloenzymes.",
    },
    "fluorine_probe": {
        "members": ["OCF3Phe"],
        "need": "a 19F NMR / MRI reporter with no natural background.",
    },
}

#: Groups defined the way the user said NOT to -- by chemical class -- kept as
#: the comparator that makes the instruction testable rather than assumed.
CHEMISTRY_GROUPS = {
    "phe_para_substituted": {
        "members": ["pBpa", "pAzF", "pIF", "pAF", "BipAla", "pBoroPhe", "pCMF"],
        "need": "(comparator) same parent and same attachment point.",
    },
    "tyr_ether": {
        "members": ["OAY", "pPRF", "OCF3Phe"],
        "need": "(comparator) O-substituted tyrosines -- note this is exactly "
                "the chemotype libY was trained on.",
    },
}


def table() -> pd.DataFrame:
    frame = pd.DataFrame(DESCRIPTORS, columns=COLUMNS)
    if frame["ncAA"].duplicated().any():
        raise ValueError("duplicate ncAA in DESCRIPTORS")
    return frame


def full_table() -> pd.DataFrame:
    """The TyrRS dataset's 13 ncAAs, libY's 8 training substrates, and AzK.

    One frame so that z-scoring -- and therefore every distance -- is computed
    over the same population. Distances taken from `table()` and from this are
    NOT comparable, which is why `demo1_libY_control.py` uses only this one.
    """
    frame = pd.concat([
        table(),
        pd.DataFrame(LIBY_SUBSTRATES, columns=COLUMNS),
        pd.DataFrame(AZK, columns=COLUMNS),
    ], ignore_index=True)
    if frame["ncAA"].duplicated().any():
        dup = frame[frame["ncAA"].duplicated()]["ncAA"].tolist()
        raise ValueError(f"duplicate entries in the union table: {dup}")
    verify_crosscheck(frame)
    return frame


def verify_crosscheck(frame: pd.DataFrame | None = None) -> None:
    """pPRF and libY-C are the same molecule -- make the table prove it.

    O-propargyl-L-tyrosine appears in both datasets, encoded independently from
    two different sources (a figure in the client's spreadsheet, and
    `3.training/UAA-info.csv` in the libY repo). If the two rows disagree, the
    hand-curation is wrong somewhere and every distance below is suspect, so
    this refuses rather than warns.

    ⚠️ One shared molecule is one check. It catches a systematic error in how
    `n_heavy` counts the ether oxygen -- which it did -- and nothing else.
    """
    frame = frame if frame is not None else full_table()
    pairs = [("pPRF", "libY-C OPropargyl-Tyr")]
    compare = [c for c in COLUMNS if c != "ncAA"]
    for a, b in pairs:
        rows = frame[frame["ncAA"].isin([a, b])]
        if len(rows) != 2:
            continue      # one of them is absent from this frame; nothing to check
        left = rows[rows["ncAA"] == a][compare].iloc[0]
        right = rows[rows["ncAA"] == b][compare].iloc[0]
        bad = [c for c in compare if left[c] != right[c]]
        if bad:
            detail = ", ".join(f"{c}: {left[c]!r} vs {right[c]!r}" for c in bad)
            raise ValueError(
                f"{a} and {b} are the same molecule but disagree on {detail}. "
                "Fix DESCRIPTORS / LIBY_SUBSTRATES before trusting any distance."
            )


def feature_matrix(frame: pd.DataFrame | None = None) -> tuple[pd.DataFrame, list[str]]:
    """One-hot the categorical columns, z-score the numeric ones.

    Z-scoring matters: without it `n_heavy` (range 1-10) would dominate every
    distance and the set would collapse to "big substituent vs small".
    """
    frame = table() if frame is None else frame
    cat = pd.get_dummies(frame[["parent", "attach", "subst_class"]],
                         prefix_sep="=", dtype=float)
    num = frame[["n_heavy", "n_rings", "hbd", "hba", "charge", "chelator"]].astype(float)
    num = (num - num.mean()) / num.std(ddof=0).replace(0.0, 1.0)
    out = pd.concat([cat, num], axis=1)
    out.index = frame["ncAA"]
    return out, list(out.columns)


def distance_matrix(frame: pd.DataFrame | None = None) -> pd.DataFrame:
    """Euclidean chemotype distance between every pair of ncAAs.

    ⚠️ This is a distance in a hand-built 8-property space, not a measured
    cross-reactivity. Its only validated use is the one in `demo3_evidence.py`:
    as the x-axis against which held-out transfer accuracy is *measured*. Read
    sideways -- "these two are close so a variant for one will work for the
    other" -- it is an untested assumption.

    Pass `full_table()` to include libY's substrates and AzK; the z-scoring then
    spans that larger population, so those distances are on a different scale
    from the 13-ncAA default. Do not mix the two.
    """
    feats, _ = feature_matrix(frame)
    values = feats.to_numpy(dtype=float)
    n = len(values)
    out = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d = float(np.linalg.norm(values[i] - values[j]))
        out[i, j] = out[j, i] = d
    return pd.DataFrame(out, index=feats.index, columns=feats.index)


def nearest(ncaa: str, exclude: set[str] | None = None) -> tuple[str, float]:
    """Closest other ncAA and its distance, optionally ignoring some."""
    dist = distance_matrix()
    if ncaa not in dist.index:
        raise KeyError(f"{ncaa!r} is not in the descriptor table; add it to DESCRIPTORS")
    row = dist.loc[ncaa].drop(index=[ncaa])
    if exclude:
        row = row.drop(index=[e for e in exclude if e in row.index])
    if row.empty:
        raise ValueError(f"nothing left to compare {ncaa!r} against")
    return str(row.idxmin()), float(row.min())


def main() -> int:
    frame = table()
    print("=" * 74)
    print("ncAA chemotype descriptors (hand-curated, ⚠️ not measured)")
    print("=" * 74)
    print(frame.to_string(index=False))

    print("\n" + "=" * 74)
    print("Nearest neighbour for each ncAA, in the 8-property space")
    print("=" * 74)
    rows = [{"ncAA": n, "nearest": nearest(n)[0], "distance": nearest(n)[1]}
            for n in frame["ncAA"]]
    near = pd.DataFrame(rows).sort_values("distance", ignore_index=True)
    print(near.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\n  most isolated: {near.iloc[-1]['ncAA']} "
          f"(nearest neighbour {near.iloc[-1]['distance']:.2f} away)")
    print(f"  tightest pair: {near.iloc[0]['ncAA']} / {near.iloc[0]['nearest']} "
          f"({near.iloc[0]['distance']:.2f})")

    print("\n" + "=" * 74)
    print("Groups by product need (what a multi-handle carrier requires)")
    print("=" * 74)
    dist = distance_matrix()
    for name, spec in PRODUCT_GROUPS.items():
        members = spec["members"]
        spread = (max(dist.loc[a, b] for a, b in itertools.combinations(members, 2))
                  if len(members) > 1 else 0.0)
        print(f"\n  {name}  ({len(members)} ncAA, widest internal gap {spread:.2f})")
        print(f"    members: {', '.join(members)}")
        print(f"    need:    {spec['need']}")

    print("\n" + "=" * 74)
    print("Comparator groups by chemical class -- the thing NOT to do")
    print("=" * 74)
    for name, spec in CHEMISTRY_GROUPS.items():
        members = spec["members"]
        spread = max(dist.loc[a, b] for a, b in itertools.combinations(members, 2))
        print(f"\n  {name}  ({len(members)} ncAA, widest internal gap {spread:.2f})")
        print(f"    members: {', '.join(members)}")
    print("\n  ⚠️ A wider internal gap is not automatically worse -- it is only worse")
    print("     if accuracy degrades with distance. demo3_evidence.py measures that")
    print("     rather than assuming it.")

    print("\n" + "=" * 74)
    print("Judgement calls in the table")
    print("=" * 74)
    for key, note in NOTES.items():
        print(f"\n  {key}: {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
