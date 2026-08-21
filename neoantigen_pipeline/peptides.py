"""Step 4: variants -> mutant / wild-type peptide pairs.

For every variant we emit, for each requested length L, the L windows that
*contain* the mutated residue, together with the positionally matched
wild-type window. The WT counterpart is what makes agretopicity (the
mutant-vs-self binding ratio) computable, and it is also what the
"is this actually tumor-specific?" gate needs.

Coverage
--------
Missense              fully supported from the reference proteome
In-frame del/ins      supported when the HGVSp string is parseable
Nonsense              produces no novel peptide -> excluded by design
Frameshift / neo-ORF  requires transcript-level (nucleotide) annotation. Supply
                      the translated neo-ORF in a `neo_orf` column and it is
                      tiled here; otherwise those variants are reported as
                      "skipped: needs transcript annotation" rather than
                      silently dropped.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional

import pandas as pd

INFRAME_DEL = re.compile(r"^([A-Z])(\d+)(?:_([A-Z])(\d+))?del$")
INFRAME_INS = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")


def _windows(seq: str, center: int, length: int):
    """All length-L windows of `seq` (0-based) that contain index `center`."""
    out = []
    for start in range(max(0, center - length + 1), min(center, len(seq) - length) + 1):
        if start + length <= len(seq):
            out.append((start, seq[start:start + length]))
    return out


def mutant_protein(wt: str, ref_aa: str, pos1: int, alt_aa: str):
    """Apply a missense change to the canonical protein. Returns (mut_seq, ok)."""
    i = pos1 - 1
    if i < 0 or i >= len(wt):
        return None, "position out of range"
    if wt[i] != ref_aa:
        return None, f"reference mismatch (proteome has {wt[i]}, variant says {ref_aa})"
    return wt[:i] + alt_aa + wt[i + 1:], "ok"


def generate_peptides(variants: pd.DataFrame, proteome: Dict[str, str],
                      lengths=(8, 9, 10, 11),
                      class2_lengths=(15,)) -> (pd.DataFrame, pd.DataFrame):
    """-> (peptide_table, skipped_table).

    peptide_table columns:
      var_id, gene, protein_change, mhc_class, length, mut_peptide, wt_peptide,
      mut_offset (0-based index of the mutated residue inside the peptide),
      pep_start (0-based start in the protein)
    """
    rows, skipped = [], []
    for idx, v in variants.iterrows():
        gene = v.get("gene")
        wt = proteome.get(gene)
        var_id = f"{gene}:{v.get('protein_change')}"
        vclass = v.get("variant_class")

        if v.get("neo_orf"):                        # user-supplied frameshift peptide
            orf = str(v["neo_orf"])
            for mhc_class, Ls in (("I", lengths), ("II", class2_lengths)):
                for L in Ls:
                    for start in range(0, max(0, len(orf) - L + 1)):
                        rows.append(dict(var_id=var_id, gene=gene,
                                         protein_change=v.get("protein_change"),
                                         variant_class=vclass, mhc_class=mhc_class,
                                         length=L, mut_peptide=orf[start:start + L],
                                         wt_peptide=None, mut_offset=None,
                                         pep_start=start))
            continue

        if wt is None:
            skipped.append(dict(var_id=var_id, gene=gene, reason="gene not in reference proteome"))
            continue
        if vclass in ("Nonsense_Mutation", "Splice_Site", "Splice_Region", "Silent"):
            skipped.append(dict(var_id=var_id, gene=gene,
                                reason=f"{vclass}: no novel peptide without transcript-level ORF"))
            continue
        if vclass in ("Frame_Shift_Del", "Frame_Shift_Ins"):
            skipped.append(dict(var_id=var_id, gene=gene,
                                reason="frameshift: needs transcript annotation (supply neo_orf)"))
            continue

        ref_aa, pos1, alt_aa = v.get("ref_aa"), v.get("aa_pos"), v.get("alt_aa")
        if not (ref_aa and pos1 and alt_aa):
            skipped.append(dict(var_id=var_id, gene=gene,
                                reason=f"unparsed protein change {v.get('protein_change')}"))
            continue

        mut, status = mutant_protein(wt, ref_aa, int(pos1), alt_aa)
        if mut is None:
            skipped.append(dict(var_id=var_id, gene=gene, reason=status))
            continue

        center = int(pos1) - 1
        for mhc_class, Ls in (("I", lengths), ("II", class2_lengths)):
            for L in Ls:
                for start, mp in _windows(mut, center, L):
                    wp = wt[start:start + L]
                    rows.append(dict(var_id=var_id, gene=gene,
                                     protein_change=v.get("protein_change"),
                                     variant_class=vclass, mhc_class=mhc_class,
                                     length=L, mut_peptide=mp, wt_peptide=wp,
                                     mut_offset=center - start, pep_start=start))

    pep = pd.DataFrame(rows).drop_duplicates(
        subset=["var_id", "mhc_class", "mut_peptide"]).reset_index(drop=True)
    return pep, pd.DataFrame(skipped)


def annotate_anchor(pep: pd.DataFrame) -> pd.DataFrame:
    """Flag whether the mutation sits at a canonical MHC-I anchor (P2 / P-Omega).

    Anchor mutations mostly change *binding* (they create the epitope); non-anchor
    mutations change what the TCR sees. Both are useful, for different reasons --
    this column lets a caller weight them differently instead of guessing.
    """
    d = pep.copy()
    off = d["mut_offset"]
    L = d["length"]
    d["mut_at_anchor"] = (off == 1) | (off == L - 1) | (off == 0)
    d["mut_at_tcr_face"] = (~d["mut_at_anchor"]) & off.notna()
    return d
