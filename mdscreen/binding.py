"""End-state MM/GBSA binding free-energy estimation for ranking ligands.

Single-trajectory approximation (the standard, cheap MM/GBSA variant):

    ΔG_bind ≈ <E_complex> − <E_receptor> − <E_ligand>

evaluated on the *dry* (solvent-stripped) frames of the complex trajectory,
with solvation modelled implicitly (GBn2). Configurational entropy (−TΔS) is
omitted — this is appropriate for *ranking* a congeneric series, not for
absolute affinities. The more negative the score, the stronger the predicted
binder.

Atom-ordering assumption: complex systems built by ``prepare.build_complex``
place protein atoms first, ligand atoms second, solvent last. We therefore map
trajectory coordinates by index: [0:n_prot] = receptor, [n_prot:n_prot+n_lig]
= ligand, and the concatenation = complex.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .config import SimConfig
from .prepare import fix_protein, load_ligand
from .utils import log, ensure_dir, write_json


@dataclass
class BindingResult:
    dg_bind_kcal: float
    dg_std_kcal: float
    e_complex: float
    e_receptor: float
    e_ligand: float
    n_frames_used: int


def _implicit_generator(cfg: SimConfig, molecules):
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator
    ff_kwargs = {
        "constraints": None,
        "rigidWater": False,
        "removeCMMotion": False,
    }
    return SystemGenerator(
        forcefields=[cfg.protein_forcefield, "implicit/gbn2.xml"],
        small_molecule_forcefield=cfg.small_molecule_forcefield,
        molecules=molecules if molecules else None,
        forcefield_kwargs=ff_kwargs,
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
        cache=None,
    )


def _energy(system, topology, positions_nm) -> float:
    """Potential energy (kcal/mol) of `system` at the given positions."""
    from openmm import Context, VerletIntegrator, unit
    from openmm import Platform
    integ = VerletIntegrator(1.0 * unit.femtoseconds)
    try:
        ctx = Context(system, integ, Platform.getPlatformByName("CPU"))
    except Exception:
        ctx = Context(system, integ)
    ctx.setPositions(positions_nm * unit.nanometer)
    e = ctx.getState(getEnergy=True).getPotentialEnergy()
    val = e.value_in_unit(unit.kilocalorie_per_mole)
    del ctx, integ
    return val


def compute_mmgbsa(topology_pdb: str, trajectory_dcd: str,
                   receptor_source: str, ligand_source: str,
                   cfg: SimConfig, outdir: str, ligand_name: str = "LIG",
                   stride: int = 5, max_frames: int = 50) -> BindingResult:
    import MDAnalysis as mda
    from openmm import app

    ensure_dir(outdir)

    # Rebuild dry topologies (must reproduce build_complex atom ordering).
    prot_top, prot_pos = fix_protein(receptor_source)
    ligand = load_ligand(ligand_source, ligand_name, cfg.ligand_charge_method)

    # complex = protein first, ligand second (matches prepare.build_complex)
    modeller = app.Modeller(prot_top, prot_pos)
    modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())
    complex_top = modeller.topology

    gen_complex = _implicit_generator(cfg, [ligand])
    sys_complex = gen_complex.create_system(complex_top)
    gen_rec = _implicit_generator(cfg, [ligand])
    sys_rec = gen_rec.create_system(prot_top)
    lig_top = ligand.to_topology().to_openmm()
    gen_lig = _implicit_generator(cfg, [ligand])
    sys_lig = gen_lig.create_system(lig_top)

    n_prot = prot_top.getNumAtoms()
    n_lig = lig_top.getNumAtoms()
    n_complex = n_prot + n_lig
    log.info("MM/GBSA: n_prot=%d n_lig=%d", n_prot, n_lig)

    u = mda.Universe(topology_pdb, trajectory_dcd)
    frames = list(range(0, len(u.trajectory), stride))[:max_frames]

    e_c, e_r, e_l = [], [], []
    for fi in frames:
        u.trajectory[fi]
        # positions of the first n_complex atoms, Å -> nm
        coords = u.atoms.positions[:n_complex] / 10.0
        rec_xyz = coords[:n_prot]
        lig_xyz = coords[n_prot:n_complex]
        e_c.append(_energy(sys_complex, complex_top, coords))
        e_r.append(_energy(sys_rec, prot_top, rec_xyz))
        e_l.append(_energy(sys_lig, lig_top, lig_xyz))

    e_c, e_r, e_l = map(np.array, (e_c, e_r, e_l))
    dg = e_c - e_r - e_l
    result = BindingResult(
        dg_bind_kcal=float(np.mean(dg)),
        dg_std_kcal=float(np.std(dg)),
        e_complex=float(np.mean(e_c)),
        e_receptor=float(np.mean(e_r)),
        e_ligand=float(np.mean(e_l)),
        n_frames_used=len(frames),
    )
    write_json(os.path.join(outdir, "mmgbsa.json"), result.__dict__)
    log.info("MM/GBSA ΔG_bind = %.2f ± %.2f kcal/mol (%d frames)",
             result.dg_bind_kcal, result.dg_std_kcal, result.n_frames_used)
    return result
