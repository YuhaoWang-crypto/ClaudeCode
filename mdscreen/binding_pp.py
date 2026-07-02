"""Protein–protein / protein–peptide end-state MM/GBSA + interface stability.

The stock ``binding.py`` assumes a protein *receptor* + a small-molecule *ligand*
(OpenFF/GAFF2). For an antibody/peptide *binder* against a peptide *target*
(both proteins) we instead:

  * clean the docked complex with PDBFixer (ff14SB-compatible; non-standard
    residues such as Aib are replaced by their standard parent),
  * define the two components purely by *atom index* in the fixed complex
    (chain-0 atoms = binder, chain-1 atoms = target — the order they appear in
    the input PDB), which is robust to chain-ID renaming by OpenMM,
  * compute single-trajectory end-state MM/GBSA:

        ΔG_bind ≈ <E_complex> − <E_binder> − <E_target>

    on solvent-stripped frames with implicit GBn2 solvent (entropy omitted →
    a *ranking* score, not an absolute affinity),
  * and measure pose stability: binder Cα-RMSD after superposing on the target,
    plus the inter-chain heavy-atom contact count over the trajectory (does the
    Boltz-predicted pose stay bound?).

Atom-order contract: matches ``prepare.build_protein`` — protein atoms first
(binder then target, in input-PDB order), solvent last. We only ever touch the
first ``n_prot`` atoms of each trajectory frame.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, asdict
from typing import List, Tuple

import numpy as np

from .config import SimConfig
from .prepare import fix_protein
from .utils import log, ensure_dir, write_json


@dataclass
class PPBindingResult:
    dg_bind_kcal: float
    dg_std_kcal: float
    e_complex: float
    e_binder: float
    e_target: float
    n_frames_used: int
    # pose stability
    binder_rmsd_mean_nm: float
    binder_rmsd_max_nm: float
    interface_contacts_t0: int
    interface_contacts_mean: float
    contact_retention: float          # mean/t0 — 1.0 = fully retained, ~0 = unbound
    n_binder_atoms: int
    n_target_atoms: int


def _chain_atom_counts(top) -> List[int]:
    return [sum(1 for _ in ch.atoms()) for ch in top.chains()]


def _implicit_system(topology):
    """ff14SB + GBn2 implicit-solvent System for a (sub)topology of protein atoms."""
    from openmm import app
    ff = app.ForceField("amber/ff14SB.xml", "implicit/gbn2.xml")
    return ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                           constraints=None, rigidWater=False, removeCMMotion=False)


def _energy(system, positions_nm) -> float:
    from openmm import Context, VerletIntegrator, Platform, unit
    integ = VerletIntegrator(1.0 * unit.femtoseconds)
    try:
        ctx = Context(system, integ, Platform.getPlatformByName("CPU"))
    except Exception:
        ctx = Context(system, integ)
    ctx.setPositions(positions_nm * unit.nanometer)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    del ctx, integ
    return e


def _submodeller(fixed_top, fixed_pos, keep_range: Tuple[int, int]):
    """Return an OpenMM Modeller keeping only atoms in [lo, hi)."""
    from openmm import app
    m = app.Modeller(fixed_top, fixed_pos)
    atoms = list(m.topology.atoms())
    lo, hi = keep_range
    to_delete = [a for i, a in enumerate(atoms) if not (lo <= i < hi)]
    m.delete(to_delete)
    return m


def compute_mmgbsa_pp(complex_pdb: str, topology_pdb: str, trajectory_dcd: str,
                      cfg: SimConfig, outdir: str,
                      stride: int = 5, max_frames: int = 50,
                      contact_cutoff_ang: float = 4.5) -> PPBindingResult:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms

    ensure_dir(outdir)

    # 1) Rebuild the dry, fixed complex (same recipe as build_protein) and get
    #    the per-chain atom split by *order*: chain0 = binder, chain1 = target.
    fixed_top, fixed_pos = fix_protein(complex_pdb)
    counts = _chain_atom_counts(fixed_top)
    if len(counts) < 2:
        raise ValueError(f"expected >=2 protein chains, got {len(counts)}")
    n_binder = counts[0]
    n_target = counts[1]
    n_prot = n_binder + n_target
    log.info("MM/GBSA(pp): n_binder=%d n_target=%d n_prot=%d", n_binder, n_target, n_prot)

    # 2) Implicit-solvent systems for complex / binder / target.
    sys_complex = _implicit_system(fixed_top)
    binder_m = _submodeller(fixed_top, fixed_pos, (0, n_binder))
    target_m = _submodeller(fixed_top, fixed_pos, (n_binder, n_prot))
    sys_binder = _implicit_system(binder_m.topology)
    sys_target = _implicit_system(target_m.topology)

    # 3) Load trajectory, keep protein molecules WHOLE across the periodic box.
    #    DCD frames are wrapped, so a naive read splits chains across the box
    #    (interface contacts -> 0, MM/GBSA -> ~0). We unwrap the protein by its
    #    bond graph (guessed from the compact first frame), then per frame bring
    #    the target chain into the binder's periodic image (orthorhombic min-image
    #    on chain centres of mass). All energetics/contacts use these coords.
    from MDAnalysis.analysis import distances as _dist
    u = mda.Universe(topology_pdb, trajectory_dcd)
    prot = u.atoms[:n_prot]
    u.trajectory[0]
    try:
        prot.guess_bonds()
        from MDAnalysis import transformations as _trans
        u.trajectory.add_transformations(_trans.unwrap(prot))
    except Exception as exc:
        log.warning("unwrap setup failed (%s); proceeding without", exc)

    names = np.array([a.name for a in prot.atoms])
    is_ca = names == "CA"
    is_heavy = ~np.char.startswith(names, "H")
    b_ca = np.where(is_ca[:n_binder])[0]
    t_ca = np.where(is_ca[n_binder:n_prot])[0] + n_binder
    b_hv = np.where(is_heavy[:n_binder])[0]
    t_hv = np.where(is_heavy[n_binder:n_prot])[0] + n_binder

    def imaged(fi):
        """Protein coords (Å) with target min-imaged into the binder's box image."""
        u.trajectory[fi]
        xyz = prot.positions.copy()
        box = u.dimensions[:3]
        if box is not None and np.all(box > 0):
            bcom = xyz[:n_binder].mean(0)
            tcom = xyz[n_binder:n_prot].mean(0)
            shift = np.round((tcom - bcom) / box) * box
            xyz[n_binder:n_prot] -= shift
        return xyz

    frames = list(range(0, len(u.trajectory), stride))[:max_frames]
    xyz0 = imaged(frames[0])
    ref_b_ca = xyz0[b_ca]; ref_t_ca = xyz0[t_ca]
    ref_t_com = ref_t_ca.mean(0)

    from MDAnalysis.analysis.align import rotation_matrix
    e_c, e_b, e_t, rmsds, contacts = [], [], [], [], []
    for fi in frames:
        xyz = imaged(fi)
        nm = xyz / 10.0
        e_c.append(_energy(sys_complex, nm))
        e_b.append(_energy(sys_binder, nm[:n_binder]))
        e_t.append(_energy(sys_target, nm[n_binder:n_prot]))
        # interface heavy-atom contacts
        dm = _dist.distance_array(xyz[b_hv], xyz[t_hv])
        contacts.append(int((dm < contact_cutoff_ang).sum()))
        # binder drift after superposing target onto frame-0 target
        t_ca_f = xyz[t_ca]; t_com_f = t_ca_f.mean(0)
        R, _ = rotation_matrix(t_ca_f - t_com_f, ref_t_ca - ref_t_com)
        moved = (xyz[b_ca] - t_com_f) @ R.T + ref_t_com
        d = moved - ref_b_ca
        rmsds.append(float(np.sqrt((d * d).sum(1).mean())) / 10.0)   # nm
    e_c, e_b, e_t = map(np.array, (e_c, e_b, e_t))
    dg = e_c - e_b - e_t
    rmsds = np.array(rmsds); contacts = np.array(contacts, dtype=float)
    c0 = float(contacts[0]) if len(contacts) else 0.0
    cmean = float(contacts.mean()) if len(contacts) else 0.0

    res = PPBindingResult(
        dg_bind_kcal=float(dg.mean()), dg_std_kcal=float(dg.std()),
        e_complex=float(e_c.mean()), e_binder=float(e_b.mean()), e_target=float(e_t.mean()),
        n_frames_used=len(frames),
        binder_rmsd_mean_nm=float(rmsds.mean()), binder_rmsd_max_nm=float(rmsds.max()),
        interface_contacts_t0=int(c0), interface_contacts_mean=cmean,
        contact_retention=float(cmean / c0) if c0 > 0 else 0.0,
        n_binder_atoms=n_binder, n_target_atoms=n_target,
    )
    write_json(os.path.join(outdir, "mmgbsa_pp.json"), asdict(res))
    log.info("MM/GBSA(pp) dG=%.2f+/-%.2f kcal/mol | binderRMSD=%.2f nm | contacts %.0f->%.0f (ret %.2f)",
             res.dg_bind_kcal, res.dg_std_kcal, res.binder_rmsd_mean_nm,
             res.interface_contacts_t0, res.interface_contacts_mean, res.contact_retention)
    return res
