"""
screen.py -- build a focused chimera library and orchestrate scoring.

Mirrors the paper's bench workflow: because ML-designed binders are small,
the circular-permutation library is tiny ("fewer than ten variants"), so we
enumerate one chimera per structure-derived loop site and rank them.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from .design import build_chimera, verify_chimera, Chimera
from .systems import System


def gs_linker_for(n_gap_residues: int = 4) -> str:
    """A flexible Gly/Ser linker 'of sufficient length' (paper's phrasing).

    Default GS pattern; length is a design knob validated in silico.  The
    paper shows linker length trades dynamic range against kcat, so this is a
    parameter, not a fixed constant.
    """
    unit = "GGS"
    reps = max(1, (n_gap_residues + len(unit) - 1) // len(unit))
    return (unit * reps)[: max(len(unit), n_gap_residues)]


def build_library(system: System, gs_linker: str = "GGSGGSGGS") -> list[Chimera]:
    """One chimera per candidate permutation site (the whole focused library)."""
    rc = system.receptor
    rp = system.reporter
    q = rp.insertion_sites[system.primary_insertion]

    lib = []
    for site in rc.loop_sites:
        ch = build_chimera(
            name=f"cp{rc.name}-{site}_{rp.name}-{system.primary_insertion}",
            reporter_name=rp.name,
            reporter_seq=rp.seq,
            receptor_name=rc.name,
            receptor_seq=rc.seq,
            insertion_index=q,
            permutation_site=site,
            gs_linker=gs_linker,
        )
        chk = verify_chimera(ch, rp.seq, rc.seq)
        ch.meta["verify"] = chk
        ch.meta["ligand"] = {"name": rc.ligand_name, "smiles": rc.ligand_smiles}
        lib.append(ch)
    return lib


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
