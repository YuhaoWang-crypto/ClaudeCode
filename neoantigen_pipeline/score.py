"""Step 6a: features -> one composite neoantigen score, and the gates that
decide what is even eligible.

Two deliberately separate layers:

  gates   binary, biology-driven, non-negotiable ("must be expressed",
          "must be presented", "must not exist in the self proteome")
  score   a weighted, ordered preference among whatever survived the gates

Keeping them apart is what makes the output auditable: a candidate that fails
is reported with the gate that killed it, not silently ranked last.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

FEATURE_COLS = {
    "presentation": "feat_presentation",
    "agretopicity": "feat_agretopicity",
    "expression": "feat_expression",
    "clonality": "feat_clonality",
    "dissimilarity": "feat_dissimilarity",
    "tcr_prior": "feat_tcr_prior",
    "hydrophobicity": "feat_hydrophobicity",
    "mhc2_support": "feat_mhc2_support",
}


def apply_peptide_gates(d: pd.DataFrame, gates) -> pd.DataFrame:
    d = d.copy()
    d["gate_presented"] = d["mut_rank"].fillna(99) <= gates.max_rank_mhc1
    d["gate_novel"] = d["is_novel_vs_self"] if gates.require_novel_vs_self else True
    d["gate_expressed"] = d["tpm"].fillna(0) >= gates.min_tpm if "tpm" in d else True
    d["gate_clonal"] = d["ccf"].fillna(0) >= gates.min_ccf if "ccf" in d else True
    if gates.drop_anchor_only and "mut_at_anchor" in d:
        d["gate_tcr_face"] = ~d["mut_at_anchor"].fillna(False)
    else:
        d["gate_tcr_face"] = True
    cols = [c for c in d.columns if c.startswith("gate_")]
    d["passes"] = d[cols].fillna(False).all(axis=1)
    return d


def composite_score(d: pd.DataFrame, weights: Dict[str, float],
                    anchor_penalty: float = 0.85) -> pd.DataFrame:
    """Weighted sum of the [0,1] features, with an optional anchor-only modifier.

    A mutation that only changes an MHC anchor residue makes the peptide bind
    but leaves the TCR-facing surface identical to self; it is down-weighted,
    not excluded (Duan 2014 show both classes can be immunogenic).
    """
    d = d.copy()
    s = np.zeros(len(d))
    used = {}
    for name, w in weights.items():
        col = FEATURE_COLS.get(name)
        if col is None or col not in d.columns:
            continue
        s = s + w * d[col].fillna(0).to_numpy(dtype=float)
        used[name] = w
    d["_weights_used"] = [used] * len(d)
    if "mut_at_anchor" in d.columns:
        pen = np.where(d["mut_at_anchor"].fillna(False).to_numpy(), anchor_penalty, 1.0)
        s = s * pen
    d["neo_score"] = s
    return d


def best_per_variant(d: pd.DataFrame) -> pd.DataFrame:
    """A vaccine encodes one minigene per *mutation*, so collapse peptide x allele
    to the single best-scoring representative epitope of each variant, while
    keeping the supporting counts."""
    if d.empty:
        return d
    grp = d.sort_values("neo_score", ascending=False).groupby("var_id", as_index=False)
    best = grp.head(1).copy()
    counts = (d.groupby("var_id")
                .agg(n_epitopes=("mut_peptide", "nunique"),
                     n_alleles=("allele", "nunique"),
                     n_strong=("binder", lambda x: int((x == "strong").sum())),
                     best_rank=("mut_rank", "min"))
                .reset_index())
    best = best.merge(counts, on="var_id", how="left")
    return best.sort_values("neo_score", ascending=False).reset_index(drop=True)


def rank_table(best: pd.DataFrame) -> pd.DataFrame:
    cols = ["var_id", "gene", "protein_change", "allele", "mut_peptide", "wt_peptide",
            "length", "mut_offset", "mut_rank", "wt_rank", "binder", "tpm", "dna_vaf",
            "ccf", "n_epitopes", "n_alleles", "n_strong",
            "feat_presentation", "feat_agretopicity", "feat_expression",
            "feat_clonality", "feat_dissimilarity", "feat_tcr_prior",
            "feat_hydrophobicity", "feat_mhc2_support", "neo_score"]
    have = [c for c in cols if c in best.columns]
    out = best[have].copy()
    out.insert(0, "rank", range(1, len(out) + 1))
    return out
