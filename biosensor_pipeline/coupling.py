"""
coupling.py -- allosteric-coupling readout from apo vs holo predicted models.

The dynamic range of these switches comes from ligand binding (at the receptor
domain) allosterically ordering the *reporter* active site. Boltz writes its
per-residue confidence (pLDDT) into the B-factor column of the output CIF, so we
can read, for free, whether the reporter's catalytic residues become MORE
ordered when the ligand is present -- a structural echo of the paper's
entropic ON-state mechanism.

  ✅ measurement: per-residue pLDDT (from the model) + catalytic constellation.
  ⚠️ interpretation: "active site more ordered in holo ⇒ ON-coupling". pLDDT is
     model confidence, not a physical order parameter -- MD (md_entropy.py) is
     the rigorous route. This is a fast triage signal, never a dynamic range.
"""

from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
import os

from .systems import TEM1, get_system
from .screen import build_library
from .scoring import map_reporter_index

OUT = os.path.join(os.getcwd(), "biosensor_out")


def _ca(cif_path):
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    f = pdbx.CIFFile.read(cif_path)
    arr = pdbx.get_structure(f, model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    chain = max(set(arr.chain_id), key=lambda c: (arr.chain_id == c).sum())
    arr = arr[arr.chain_id == chain]
    return arr[arr.atom_name == "CA"]


def active_site_order(cif_path, chimera) -> dict:
    """Mean/per-residue pLDDT at the reporter catalytic residues + global mean."""
    import numpy as np
    ca = _ca(cif_path)
    want = {lab: map_reporter_index(chimera, idx) + 1 for lab, idx in TEM1.catalytic.items()}
    per = {}
    for lab, rid in want.items():
        sel = ca[ca.res_id == rid]
        if len(sel):
            per[lab] = round(float(sel.b_factor[0]), 2)
    vals = list(per.values())
    return {
        "active_site_mean_plddt": round(sum(vals) / len(vals), 2) if vals else None,
        "global_mean_plddt": round(float(np.mean(ca.b_factor)), 2),
        "per_residue": per,
    }


def coupling(system_key: str, apo_cif: str, holo_cif: str) -> dict:
    """Compare apo vs holo active-site ordering for one system."""
    s = get_system(system_key)
    ch = build_library(s)[len(s.receptor.loop_sites) // 2]
    apo = active_site_order(apo_cif, ch)
    holo = active_site_order(holo_cif, ch)
    d_as = round(holo["active_site_mean_plddt"] - apo["active_site_mean_plddt"], 2)
    d_gl = round(holo["global_mean_plddt"] - apo["global_mean_plddt"], 2)
    return {
        "system": system_key, "analyte": s.receptor.ligand_name,
        "apo": apo, "holo": holo,
        "delta_active_site_plddt": d_as,
        "delta_global_plddt": d_gl,
        "on_coupling_hint": d_as > 0,   # ⚠️ hypothesis
        "label": "⚠️ pLDDT ordering is a coupling *hint*, not a dynamic range; "
                 "confirm with MD (md_entropy.py) and, ultimately, a wet-lab kobs assay.",
    }


def main():
    import json
    # dfhbi is the system with both apo and holo CIFs downloaded
    pairs = [("dfhbi", "dfhbi_apo.cif", "dfhbi_holo.cif")]
    results = []
    for key, apo, holo in pairs:
        ap, ho = os.path.join(OUT, apo), os.path.join(OUT, holo)
        if not (os.path.exists(ap) and os.path.exists(ho)):
            print(f"[skip] {key}: need {apo} and {holo}")
            continue
        r = coupling(key, ap, ho)
        results.append(r)
        print(f"=== {key} ({r['analyte']}) allosteric-coupling readout ===")
        print(f"  apo : active-site pLDDT {r['apo']['active_site_mean_plddt']}  global {r['apo']['global_mean_plddt']}")
        print(f"  holo: active-site pLDDT {r['holo']['active_site_mean_plddt']}  global {r['holo']['global_mean_plddt']}")
        print(f"  Δ active-site = {r['delta_active_site_plddt']:+.2f}  Δ global = {r['delta_global_plddt']:+.2f}  "
              f"ON-coupling hint = {r['on_coupling_hint']}")
        print("  " + r["label"])
    if results:
        with open(os.path.join(OUT, "coupling_readout.json"), "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"\nwrote {OUT}/coupling_readout.json")


if __name__ == "__main__":
    main()
