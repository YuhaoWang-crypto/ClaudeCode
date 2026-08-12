"""Target definitions for the NPY / PP(PPY) EAB-aptamer campaign.

All peptide sequences are the MATURE, processed hormones taken from UniProt
(signal peptide and the C-terminal GKR amidation signal removed).  All three
family members are exactly 36 aa and align without gaps, so a positional
comparison is the correct alignment.

Honesty labels used throughout this package:
  [MEASURED]  - experimental value from a cited publication / database
  [DERIVED]   - computed deterministically from measured data (e.g. % identity)
  [PREDICTED] - output of a model (ViennaRNA folding, design heuristics)
  [DESIGN]    - a proposed construct; no experimental support yet
"""

# ---------------------------------------------------------------- peptides --
# [MEASURED] UniProt P01303 (NPY_HUMAN),   mature chain 29-64, C-terminally amidated
NPY = "YPSKPDNPGEDAPAEDMARYYSALRHYINLITRQRY"
# [MEASURED] UniProt P01298 (PAHO_HUMAN),  mature chain 30-65, C-terminally amidated
PP = "APLEPVYPGDNATPEQMAQYAADLRRYINMLTRPRY"
# [MEASURED] UniProt P10082 (PYY_HUMAN),   mature chain 29-64, C-terminally amidated
PYY = "YPIKPEAPREDASPEELNRYYASLRHYLNLVTRQRY"

# [MEASURED] DPP-4 removes Tyr1-Pro2 from NPY and PYY in circulation, producing
# NPY(3-36) / PYY(3-36).  An assay for "total NPY" must therefore not depend on
# residues 1-2; an assay that must distinguish the two forms should target them.
NPY_3_36 = NPY[2:]
PYY_3_36 = PYY[2:]

FAMILY = {"NPY": NPY, "PP": PP, "PYY": PYY}

# Amino-acid property classes, used for the "conservative vs radical" call in
# the epitope scan.  Standard Taylor/Livingstone-style grouping.
_CLASS = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic",
    "I": "hydrophobic", "M": "hydrophobic", "F": "aromatic",
    "W": "aromatic", "Y": "aromatic", "P": "proline", "G": "small",
    "S": "polar", "T": "polar", "C": "polar", "N": "polar", "Q": "polar",
    "D": "acidic", "E": "acidic", "K": "basic", "R": "basic", "H": "basic",
}


def identity(a, b):
    """[DERIVED] Fraction of identical positions in a gap-free alignment."""
    assert len(a) == len(b), "family members are all 36 aa; no gaps expected"
    return sum(x == y for x, y in zip(a, b)) / len(a)


def diff_positions(a, b):
    """[DERIVED] 1-based positions where two family members differ."""
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def epitope_scan(on_target, off_targets, window=12, radical_bonus=0.5):
    """[DERIVED] Slide a window along the on-target peptide and score how
    discriminating each window is against every off-target.

    discriminating_score = (# positions differing from EVERY off-target
                            + radical_bonus * # of those that are also a
                              change of amino-acid CLASS) / window

    A position only counts if it differs from *all* off-targets - an aptamer
    epitope that separates NPY from PP but not from PYY is not useful for a
    family assay.  The class bonus rewards windows where the difference is
    chemically radical (D->Y is easier for a binder to read than I->L).

    Returns a list of dicts sorted by start position (1-based, inclusive).
    """
    n = len(on_target)
    out = []
    for start in range(0, n - window + 1):
        end = start + window
        uniq, radical = 0, 0
        uniq_pos = []
        for i in range(start, end):
            res = on_target[i]
            if all(res != off[i] for off in off_targets):
                uniq += 1
                uniq_pos.append(i + 1)
                if any(_CLASS[res] != _CLASS[off[i]] for off in off_targets):
                    radical += 1
        score = (uniq + radical_bonus * radical) / window
        out.append({
            "start": start + 1, "end": end,
            "seq": on_target[start:end],
            "n_unique": uniq,
            "n_radical": radical,
            "unique_positions": uniq_pos,
            "score": round(score, 3),
        })
    return out


def best_epitopes(on_target, off_targets, window=12, top=5):
    """[DERIVED] Top-scoring discriminating windows, best first."""
    scan = epitope_scan(on_target, off_targets, window=window)
    return sorted(scan, key=lambda r: (-r["score"], r["start"]))[:top]


# ------------------------------------------------------------ structure ------
# [MEASURED] Family fold context, relevant to whether a *fragment* is a valid
# SELEX surrogate for the intact hormone:
#   * Pancreatic polypeptide adopts the compact "PP-fold" - an N-terminal
#     polyproline-II helix packed against a C-terminal alpha-helix through a
#     hydrophobic core (the classic avian-PP crystal structure, PDB 1PPT).
#     A linear N-terminal fragment therefore does NOT reproduce the surface
#     presented by intact PP.
#   * NPY is substantially less ordered in plain aqueous buffer and only
#     becomes helical on membrane / micelle binding, so a linear NPY fragment
#     is a much better surrogate for the free hormone than a PP fragment is
#     for free PP.
# Practical consequence, encoded in the design: fragment-based selection is
# lower risk for NPY than for PP.  For PP, select against the intact hormone
# and use the fragment only for counter-/negative-selection framing.
FOLD_NOTE = {
    "NPY": "largely disordered/partially helical in aqueous buffer; "
           "linear fragment is a reasonable SELEX surrogate",
    "PP": "compact PP-fold (PPII helix packed on alpha-helix, cf. PDB 1PPT); "
          "a linear fragment is NOT a faithful surrogate of the folded surface",
}
