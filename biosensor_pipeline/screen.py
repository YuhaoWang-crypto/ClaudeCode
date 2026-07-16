"""
screen.py -- build a focused chimera library and orchestrate scoring.

Mirrors the paper's bench workflow: because ML-designed binders are small,
the circular-permutation library is tiny ("fewer than ten variants"), so we
enumerate one chimera per structure-derived loop site and rank them.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from .design import build_chimera, build_terminal_fusion, verify_chimera, Chimera
from .systems import System


def pocket_safe_sites(loop_sites, pocket_indices, guard: int = 2):
    """Drop permutation sites within `guard` residues of any pocket residue.

    Circularly permuting *inside* the ligand pocket would destroy binding, so
    pocket-adjacent loop residues are excluded. (Integrated from the
    biosensor-chimera-design skill.) If no pocket is known, returns all sites.
    """
    if not pocket_indices:
        return list(loop_sites), []
    pocket = set(pocket_indices)
    safe, dropped = [], []
    for s in loop_sites:
        if any(abs(s - p) <= guard for p in pocket):
            dropped.append(s)
        else:
            safe.append(s)
    return safe, dropped


def gs_linker_for(n_gap_residues: int = 4) -> str:
    """A flexible Gly/Ser linker 'of sufficient length' (paper's phrasing).

    Default GS pattern; length is a design knob validated in silico.  The
    paper shows linker length trades dynamic range against kcat, so this is a
    parameter, not a fixed constant.
    """
    unit = "GGS"
    reps = max(1, (n_gap_residues + len(unit) - 1) // len(unit))
    return (unit * reps)[: max(len(unit), n_gap_residues)]


def build_library(system: System, gs_linker: str = "GGSGGSGGS",
                  pocket_indices=None, orientations=("N", "C")) -> list[Chimera]:
    """The focused chimera library, dispatched on the reporter topology.

    - insertion (TEM-1, PQQ-GDH): one chimera per candidate permutation site.
    - cp_reporter_terminal (NanoLuc): reporter permuted once; binder permuted at
      each site and fused at a new terminus in each requested orientation.

    pocket_indices (0-based receptor residues) are excluded from permutation
    (pocket-safe), if provided.
    """
    rc, rp = system.receptor, system.reporter
    sites, dropped = pocket_safe_sites(rc.loop_sites, pocket_indices or getattr(rc, "pocket_indices", None))

    lib = []
    if rp.mode == "cp_reporter_terminal":
        for site in sites:
            for orient in orientations:
                ch = build_terminal_fusion(
                    name=f"cp{rc.name}-{site}_{rp.name}-cp{rp.cp_site + 1}-{orient}",
                    reporter_name=rp.name, reporter_seq=rp.seq, reporter_cp_site=rp.cp_site,
                    receptor_name=rc.name, receptor_seq=rc.seq, permutation_site=site,
                    gs_linker=gs_linker, terminal_linker=rp.terminal_linker,
                    reporter_cp_linker=rp.reporter_cp_linker, orientation=orient,
                )
                _annotate(ch, rp, rc, dropped)
                lib.append(ch)
    else:
        q = rp.insertion_sites[system.primary_insertion]
        for site in sites:
            ch = build_chimera(
                name=f"cp{rc.name}-{site}_{rp.name}-{system.primary_insertion}",
                reporter_name=rp.name, reporter_seq=rp.seq,
                receptor_name=rc.name, receptor_seq=rc.seq,
                insertion_index=q, permutation_site=site, gs_linker=gs_linker,
            )
            _annotate(ch, rp, rc, dropped)
            lib.append(ch)
    return lib


def _annotate(ch, rp, rc, dropped):
    ch.meta["verify"] = verify_chimera(ch, rp.seq, rc.seq)
    ch.meta["ligand"] = {"name": rc.ligand_name, "smiles": rc.ligand_smiles}
    ch.meta["readout"] = rp.readout
    ch.meta["pocket_dropped_sites"] = dropped


def library_summary(system: System, lib: list[Chimera]) -> dict:
    return {
        "system": system.key,
        "role": system.role,
        "receptor": system.receptor.name,
        "reporter": system.reporter.name,
        "analyte": system.receptor.ligand_name,
        "insertion": f"{system.reporter.name}-{system.primary_insertion}",
        "n_variants": len(lib),
        "all_constructions_valid": all(c.meta["verify"]["all_ok"] for c in lib),
        "variants": [
            {
                "name": c.name,
                "permutation_site": c.permutation_site,
                "removed_residue": c.meta["verify"]["removed_residue"],
                "length": c.length,
                "receptor_domain_span": c.receptor_domain_span,
                "valid": c.meta["verify"]["all_ok"],
            }
            for c in lib
        ],
    }
