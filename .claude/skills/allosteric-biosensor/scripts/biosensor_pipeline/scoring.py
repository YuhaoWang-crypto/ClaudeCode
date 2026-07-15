"""
scoring.py -- in-silico validation metrics.

Two strictly separated tiers (this separation is the whole point):

  ✅ RIGOROUS geometric measurements read off a predicted 3D model:
       * catalytic constellation compactness (Ca-Ca distances among the
         reporter's catalytic residues) -- a real number computed from atoms.
       * per-domain confidence pulled from the predictor (pLDDT/pTM/ipTM).

  ⚠️ HYPOTHESIS-level interpretation:
       * the "switch proxy" that combines those numbers into a single
         switchability ranking.  Structure prediction does NOT compute an
         allosteric dynamic range; this proxy is an ILLUSTRATIVE heuristic,
         never a predicted DR.  A wet-lab kobs(+L)/kobs(-L) titration remains
         the only ground truth (paper, Methods "Quantification of biosensor
         performance parameters").
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import math


# ---------------------------------------------------------------------------
# index mapping: reporter residue -> chimera residue
# ---------------------------------------------------------------------------
def map_reporter_index(chimera, reporter_idx: int) -> int:
    """Map a 0-idx reporter residue to its 0-idx position in the chimera."""
    q = chimera.insertion_index
    insert_len = 2 * len(chimera.flank_linker) + chimera.meta["cp_len"]
    return reporter_idx if reporter_idx < q else reporter_idx + insert_len


# ---------------------------------------------------------------------------
# ✅ geometric measurement from a predicted model (CIF/PDB)
# ---------------------------------------------------------------------------
def catalytic_constellation(cif_path: str, chimera, reporter) -> dict:
    """Measure Ca-Ca distances among the reporter catalytic residues in a
    predicted chimera model.  Returns real distances (Angstrom) + a
    compactness summary.  ✅ (measurement); interpretation is ⚠️.

    For a catalytically competent class-A beta-lactamase the S70/K73/S130/E166
    Ca atoms sit within a compact ~8-13 A shell.  Gross expansion of this
    shell in a model is evidence the active site was disrupted by the
    insertion; preservation is necessary (not sufficient) for the ON state.
    """
    import warnings; warnings.filterwarnings("ignore")
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    import numpy as np

    cif = pdbx.CIFFile.read(cif_path)
    arr = pdbx.get_structure(cif, model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    # take the protein chain (largest)
    chains = {c: (arr.chain_id == c).sum() for c in set(arr.chain_id)}
    prot_chain = max(chains, key=chains.get)
    arr = arr[arr.chain_id == prot_chain]
    ca = arr[arr.atom_name == "CA"]

    # chimera residues are renumbered 1..N along the chain
    coords = {}
    for label, ridx in reporter.catalytic.items():
        cidx = map_reporter_index(chimera, ridx)          # 0-idx in chimera
        want = cidx + 1                                    # CA res_id (1-based)
        sel = ca[ca.res_id == want]
        if len(sel) == 1:
            coords[label] = sel.coord[0]

    labels = list(coords)
    dists = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            d = float(np.linalg.norm(coords[a] - coords[b]))
            dists[f"{a}-{b}"] = round(d, 2)

    vals = list(dists.values())
    return {
        "measured_pairs": dists,
        "n_catalytic_found": len(coords),
        "max_pair_A": round(max(vals), 2) if vals else None,
        "mean_pair_A": round(sum(vals) / len(vals), 2) if vals else None,
    }


# reference constellation for intact TEM-1 (measured from 3GMW_A CA coords).
# Populated by structure.reference_constellation(); frozen here for offline use.
TEM1_REF_CONSTELLATION = {
    "n_catalytic_found": 6,
    # filled in at runtime if the reference is recomputed; used only for a
    # relative "expansion" readout, never as a hard threshold.
}


# ---------------------------------------------------------------------------
# ⚠️ switch proxy -- an illustrative ranking, NOT a predicted dynamic range
# ---------------------------------------------------------------------------
@dataclass
class SwitchFeatures:
    """Inputs to the proxy, each sourced from a Boltz prediction."""
    fold_confidence: Optional[float] = None      # chimera complex pTM/confidence [0,1]
    receptor_plddt: Optional[float] = None       # mean pLDDT of receptor domain [0,100]
    ligand_binding: Optional[float] = None       # Boltz binding probability / affinity proxy [0,1]
    active_site_ok: Optional[bool] = None        # constellation not grossly expanded
    apo_holo_shift: Optional[float] = None       # |apo-holo| confidence change (allostery hint)
    notes: str = ""


def switch_proxy(f: SwitchFeatures) -> dict:
    """Combine features into a single [0,1] switchability proxy.

    ⚠️ HYPOTHESIS.  Weights are chosen for illustration and are NOT calibrated
    against measured dynamic ranges.  The purpose is only to give a
    *reproducible, transparent* ranking of a small design library so that the
    few constructs worth wet-lab testing can be prioritized -- exactly the
    "fewer than ten variants" triage the paper performs at the bench.
    """
    terms = {}
    if f.fold_confidence is not None:
        terms["fold"] = (0.30, _clip01(f.fold_confidence))
    if f.receptor_plddt is not None:
        terms["receptor_fold"] = (0.20, _clip01(f.receptor_plddt / 100.0))
    if f.ligand_binding is not None:
        terms["binding"] = (0.30, _clip01(f.ligand_binding))
    if f.active_site_ok is not None:
        terms["active_site"] = (0.10, 1.0 if f.active_site_ok else 0.0)
    if f.apo_holo_shift is not None:
        # a modest apo->holo change is the allosteric-coupling signature;
        # saturating at 0.2 so it rewards *some* change, not chaos.
        terms["allostery_hint"] = (0.10, _clip01(f.apo_holo_shift / 0.2))

    if not terms:
        return {"proxy": None, "components": {}, "label": "no features"}

    wsum = sum(w for w, _ in terms.values())
    score = sum(w * v for w, v in terms.values()) / wsum
    return {
        "proxy": round(score, 3),
        "components": {k: round(v, 3) for k, (_, v) in terms.items()},
        "weights": {k: w for k, (w, _) in terms.items()},
        "label": "⚠️ hypothesis: illustrative switchability proxy, NOT a predicted DR",
    }


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))
