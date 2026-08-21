"""Step 5b: peptide-level features that drive the selection score.

Nine features, each mapped to [0, 1] with "higher = better neoantigen".
Every one is traceable to a published observation; none is a vendor formula.
See reference/scoring.md for the citations and the exact rationale.

  presentation    MHC-I eluted-ligand %rank                     (NetMHCpan-4.1)
  agretopicity    log10(WT %rank / MUT %rank)                    Duan 2014, Ghorani 2018
  expression      tumor RNA abundance of the source gene         Moderna workflow step 3
  clonality       cancer-cell fraction of the variant            McGranahan 2016
  dissimilarity   BLOSUM distance of mutant from its self peptide  Richman 2019
  tcr_prior       ungapped-alignment similarity to IEDB positives  Luksza 2017 (R score)
  hydrophobicity  hydrophobicity of the TCR-facing residues      Chowell 2015, Wells 2020
  mhc2_support    a CD4 helper epitope inside the same 25mer     Kreiter 2015, Ott 2017
  anchor_penalty  mutation sits only at an MHC anchor position   Duan 2014 (modifier)
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

AA = "ARNDCQEGHILKMFPSTWYV"
AA_IDX = {a: i for i, a in enumerate(AA)}

# BLOSUM62, row/col order = AA above.
_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62 = np.array([[int(x) for x in row.split()]
                     for row in _B62.strip().splitlines()], dtype=np.int16)

# Kyte-Doolittle hydropathy
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Luksza et al. 2017 published constants for the TCR-recognition probability R
LUKSZA_A = 26.0
LUKSZA_K = 4.87


def encode(peptides: Sequence[str]) -> np.ndarray:
    """Peptides of equal length -> int matrix; unknown residues -> -1."""
    L = len(peptides[0])
    M = np.full((len(peptides), L), -1, dtype=np.int8)
    for i, p in enumerate(peptides):
        for j, c in enumerate(p[:L]):
            M[i, j] = AA_IDX.get(c, -1)
    return M


# --------------------------------------------------------------------------
# Individual features
# --------------------------------------------------------------------------

def f_presentation(rank: float, midpoint: float = 0.5, steep: float = 1.5) -> float:
    """%rank -> [0,1]; 0.5% rank maps to 0.5, 0.05% -> ~0.97, 2% -> ~0.11."""
    if rank is None or (isinstance(rank, float) and math.isnan(rank)) or rank <= 0:
        return 0.0
    return 1.0 / (1.0 + (rank / midpoint) ** steep)


def f_agretopicity(mut_rank: float, wt_rank: float, cap: float = 2.0) -> float:
    """log10(WT rank / MUT rank), scaled. >0 means the mutation *created* binding.

    Returns 0.5 when unknown, so a missing WT never masquerades as evidence.
    """
    if any(v is None or (isinstance(v, float) and math.isnan(v)) or v <= 0
           for v in (mut_rank, wt_rank)):
        return 0.5
    lr = math.log10(wt_rank / mut_rank)
    return float(np.clip((lr + cap) / (2 * cap), 0.0, 1.0))


def f_expression(tpm: float, sat: float = 100.0) -> float:
    if tpm is None or (isinstance(tpm, float) and math.isnan(tpm)) or tpm <= 0:
        return 0.0
    return float(np.clip(math.log10(tpm + 1) / math.log10(sat + 1), 0.0, 1.0))


def f_clonality(ccf: float) -> float:
    if ccf is None or (isinstance(ccf, float) and math.isnan(ccf)):
        return 0.5
    return float(np.clip(ccf, 0.0, 1.0))


def f_dissimilarity(mut: str, wt: Optional[str]) -> float:
    """How chemically different the mutant is from the self peptide it came from.

    A conservative substitution (I->V) is close to self and more likely to hit a
    tolerized repertoire; a radical one (G->W) is more foreign. Normalized from
    the BLOSUM62 score of the substituted position(s).
    """
    if not wt or len(wt) != len(mut):
        return 0.75            # neo-ORF peptide: foreign by construction, but unverified
    diffs = [(a, b) for a, b in zip(mut, wt) if a != b]
    if not diffs:
        return 0.0
    scores = [BLOSUM62[AA_IDX.get(b, 0), AA_IDX.get(a, 0)] for a, b in diffs]
    mean = float(np.mean(scores))          # -4 (radical) .. +11 (identical-ish)
    return float(np.clip((3.0 - mean) / 7.0, 0.0, 1.0))


def f_hydrophobicity(pep: str, tcr_positions: Optional[Sequence[int]] = None) -> float:
    """Mean Kyte-Doolittle over the TCR-facing residues, mapped to [0,1]."""
    L = len(pep)
    pos = tcr_positions or [i for i in range(2, L - 1)]
    vals = [KD.get(pep[i], 0.0) for i in pos if i < L]
    if not vals:
        return 0.5
    return float(np.clip((np.mean(vals) + 4.5) / 9.0, 0.0, 1.0))


def tcr_prior_scores(peptides: Sequence[str], reference: Sequence[str],
                     chunk: int = 256) -> np.ndarray:
    """Luksza-style R = Z/(1+Z) with Z = sum_j exp(-k*(a - s_j)).

    s_j is an *ungapped* BLOSUM62 alignment score against every reference
    epitope of the same length -- a fast approximation of the Smith-Waterman
    score used in the original paper (see reference/scoring.md, labelled).
    """
    if not len(peptides) or not len(reference):
        return np.zeros(len(peptides))
    L = len(peptides[0])
    ref = [r for r in reference if len(r) == L]
    if not ref:
        return np.zeros(len(peptides))
    R = encode(ref)                                   # (nref, L)
    out = np.zeros(len(peptides))
    for i in range(0, len(peptides), chunk):
        P = encode(list(peptides[i:i + chunk]))       # (nb, L)
        # (nb, nref, L) BLOSUM lookup, summed over L
        s = BLOSUM62[P[:, None, :], R[None, :, :]].sum(axis=2).astype(np.float64)
        z = np.exp(-LUKSZA_K * (LUKSZA_A - s))
        z = np.clip(z, 0, 1e12).sum(axis=1)
        out[i:i + chunk] = z / (1.0 + z)
    return out


# --------------------------------------------------------------------------
# Table-level driver
# --------------------------------------------------------------------------

def compute_features(cand: pd.DataFrame, variants: pd.DataFrame,
                     iedb_positive: Sequence[str] = (),
                     self_kmers: Optional[Dict[int, set]] = None) -> pd.DataFrame:
    """cand: joined peptide x allele table with mut_rank / wt_rank.
    variants: annotated variant table (tpm, ccf) keyed by var_id.
    """
    d = cand.copy()
    vmeta = variants.copy()
    if "var_id" not in vmeta.columns:
        vmeta["var_id"] = vmeta["gene"].astype(str) + ":" + vmeta["protein_change"].astype(str)
    vmeta = vmeta.drop_duplicates("var_id").set_index("var_id")
    for col in ("tpm", "ccf", "dna_vaf", "is_clonal"):
        if col in vmeta.columns:
            d[col] = d["var_id"].map(vmeta[col])

    d["feat_presentation"] = d["mut_rank"].map(f_presentation)
    d["feat_agretopicity"] = [f_agretopicity(m, w) for m, w in
                              zip(d["mut_rank"], d.get("wt_rank", pd.Series([None] * len(d))))]
    d["feat_expression"] = d["tpm"].map(f_expression) if "tpm" in d else 0.0
    d["feat_clonality"] = d["ccf"].map(f_clonality) if "ccf" in d else 0.5
    d["feat_dissimilarity"] = [f_dissimilarity(m, w) for m, w in
                               zip(d["mut_peptide"], d.get("wt_peptide", pd.Series([None] * len(d))))]
    d["feat_hydrophobicity"] = d["mut_peptide"].map(f_hydrophobicity)

    d["feat_tcr_prior"] = 0.0
    if len(iedb_positive):
        for L, grp in d.groupby(d["mut_peptide"].str.len()):
            peps = list(grp["mut_peptide"])
            vals = tcr_prior_scores(peps, [r for r in iedb_positive if len(r) == L])
            d.loc[grp.index, "feat_tcr_prior"] = vals

    # novelty gate: the mutant k-mer must not exist anywhere in the self proteome
    if self_kmers:
        def _novel(p):
            s = self_kmers.get(len(p))
            return True if s is None else (p not in s)
        d["is_novel_vs_self"] = d["mut_peptide"].map(_novel)
    else:
        d["is_novel_vs_self"] = True
    return d


def add_mhc2_support(d: pd.DataFrame, class2_hits: pd.DataFrame) -> pd.DataFrame:
    """feat_mhc2_support = 1 when the same variant also yields a class-II binder."""
    d = d.copy()
    if class2_hits is None or class2_hits.empty:
        d["feat_mhc2_support"] = 0.0
        return d
    good = set(class2_hits.loc[class2_hits["mut_rank"] <= 10.0, "var_id"])
    strong = set(class2_hits.loc[class2_hits["mut_rank"] <= 2.0, "var_id"])
    d["feat_mhc2_support"] = d["var_id"].map(
        lambda v: 1.0 if v in strong else (0.5 if v in good else 0.0))
    return d
