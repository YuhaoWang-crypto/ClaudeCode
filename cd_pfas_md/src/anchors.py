"""Automated APR anchor-atom resolution.

APR needs a small set of host atoms to define the pulling axis plus a guest
anchor atom. Hand-picking Amber atom names is the usual friction point; this
module resolves them heuristically from the solvated topology so the Modal run
can proceed unattended, and writes the resolved 0-based atom INDICES back into
apr_manifest.json for run_apr to consume.

Heuristic (documented so you can override):
  * host axis   = the 3 host heavy atoms closest to the host centroid projected
                  onto its principal axis, i.e. atoms framing the cavity rim.
  * guest anchor = the guest's formally charged head atom (S of a sulfonate, C of
                  a carboxylate) — the group that sits in/near the portal.

These are REASONABLE defaults for a β-CD-like macrocycle + a head-charged guest,
not a substitute for chemical judgement. Inspect the indices it writes before a
production run; override by editing apr_manifest.json.
"""
from __future__ import annotations

import json
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.anchors")


def _residue_atom_indices(universe, resname_startswith: str) -> list[int]:
    sel = universe.select_atoms(f"resname {resname_startswith}*")
    return list(sel.atoms.indices)


def resolve_into_manifest(work: Path) -> dict:
    """Populate host_axis_index / guest_anchor_index in the window manifest."""
    work = Path(work)
    manifest = json.loads((work / "apr_manifest.json").read_text())
    try:
        import MDAnalysis as mda
        import numpy as np
    except ImportError as exc:  # keep importable without MDAnalysis
        raise ImportError("MDAnalysis required for anchor resolution") from exc

    u = mda.Universe(manifest["prmtop"], manifest["inpcrd"], format="INPCRD",
                     topology_format="PRMTOP")

    # host = first non-water, non-ion residue block (the macrocycle); guest = next.
    solvent = {"WAT", "HOH", "Na+", "Cl-", "NA", "CL"}
    residues = [r for r in u.residues if r.resname not in solvent]
    if len(residues) < 2:
        raise RuntimeError("expected >=2 solute residues (host + guest) in topology")
    host_res, guest_res = residues[0], residues[1]

    host = u.atoms[[a.index for a in host_res.atoms if a.mass > 2.0]]  # heavy atoms
    centroid = host.center_of_mass()
    # principal axis of the host
    shifted = host.positions - centroid
    cov = np.cov(shifted.T)
    _, vecs = np.linalg.eigh(cov)
    axis = vecs[:, -1]
    proj = shifted @ axis
    # 3 atoms nearest the rim plane (|projection| small) frame the axis origin
    order = np.argsort(np.abs(proj))
    host_axis_idx = [int(host.atoms[i].index) for i in order[:3]]

    # guest anchor = most negatively charged heavy atom (head group)
    charges = guest_res.atoms.charges
    heavy = [a for a in guest_res.atoms if a.mass > 2.0]
    guest_anchor = min(heavy, key=lambda a: a.charge)
    guest_anchor_idx = int(guest_anchor.index)

    manifest["host_axis_index"] = host_axis_idx
    manifest["guest_anchor_index"] = guest_anchor_idx
    manifest["anchor_method"] = "auto (principal-axis rim + head charge); review before prod"
    (work / "apr_manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info("resolved anchors: host_axis=%s guest_anchor=%d",
             host_axis_idx, guest_anchor_idx)
    log.warning("anchor indices are HEURISTIC — inspect apr_manifest.json before a "
                "production APR run.")
    return manifest
