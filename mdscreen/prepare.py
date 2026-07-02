"""System preparation: protein cleanup, ligand parameterisation, solvation.

This replaces the AmberTools/tleap machinery used by the making-it-rain
notebooks with the modern OpenMM/OpenFF stack:

  * PDBFixer          -> add missing atoms/residues, protonate, cap
  * OpenFF Molecule   -> read ligand (SMILES or SDF), assign charges/conformers
  * SystemGenerator   -> GAFF2 (or OpenFF) ligand params + ff14SB protein params
  * Modeller          -> merge protein+ligand, add solvent + ions

Everything here runs on the GPU worker but only needs a CPU; the heavy MD is
in ``simulate.py``.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Tuple

from .config import SimConfig
from .utils import log, is_smiles


@dataclass
class PreparedSystem:
    """A fully parameterised, solvated OpenMM system ready to integrate."""
    system: object                 # openmm.System
    topology: object               # openmm.app.Topology
    positions: object              # openmm.unit.Quantity (Nx3)
    ligand_resname: Optional[str]  # residue name of the ligand, if any
    n_atoms: int


# ---------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------
def fix_protein(pdb_input: str, ph: float = 7.0, keep_water: bool = False):
    """Clean a raw PDB with PDBFixer and return (topology, positions).

    ``pdb_input`` may be a local path or a 4-letter PDB id (downloaded).
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile  # noqa: F401 (imported for side effects/typing)

    if len(pdb_input) == 4 and not os.path.exists(pdb_input):
        log.info("Fetching PDB %s from RCSB", pdb_input)
        fixer = PDBFixer(pdbid=pdb_input)
    else:
        fixer = PDBFixer(filename=pdb_input)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    if not keep_water:
        fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    log.info("Protein prepared: %d chains", fixer.topology.getNumChains())
    return fixer.topology, fixer.positions


# ---------------------------------------------------------------------------
# Ligand
# ---------------------------------------------------------------------------
def load_ligand(source: str, name: str = "LIG", charge_method: str = "am1bcc"):
    """Load a ligand from SMILES or a file into an OpenFF Molecule.

    A 3D conformer and partial charges are generated so the molecule can be
    parameterised by the SystemGenerator.
    """
    from openff.toolkit import Molecule

    if is_smiles(source):
        log.info("Parsing ligand from SMILES: %s", source)
        mol = Molecule.from_smiles(source, allow_undefined_stereo=True)
        mol.generate_conformers(n_conformers=1)
    else:
        log.info("Loading ligand from file: %s", source)
        mol = Molecule.from_file(source, allow_undefined_stereo=True)
        if isinstance(mol, list):
            mol = mol[0]
        if mol.n_conformers == 0:
            mol.generate_conformers(n_conformers=1)

    try:
        if charge_method == "am1bcc":
            mol.assign_partial_charges("am1bcc")
        else:
            mol.assign_partial_charges("gasteiger")
    except Exception as exc:  # pragma: no cover - depends on toolkit backend
        log.warning("am1bcc charge assignment failed (%s); using gasteiger", exc)
        mol.assign_partial_charges("gasteiger")

    mol.name = name
    return mol


def _make_system_generator(cfg: SimConfig, molecules, solvate: bool):
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator

    constraints = {
        "HBonds": app.HBonds,
        "AllBonds": app.AllBonds,
        "None": None,
        None: None,
    }[cfg.constraints]

    forcefields = [cfg.protein_forcefield]
    if solvate:
        forcefields.append(cfg.water_forcefield)

    ff_kwargs = {
        "constraints": constraints,
        "rigidWater": True,
        "removeCMMotion": True,
        "hydrogenMass": cfg.hydrogen_mass_amu * unit.amu,
    }
    periodic_kwargs = {
        "nonbondedMethod": app.PME if solvate else app.NoCutoff,
    }
    if solvate:
        periodic_kwargs["nonbondedCutoff"] = cfg.nonbonded_cutoff_nm * unit.nanometer

    return SystemGenerator(
        forcefields=forcefields,
        small_molecule_forcefield=cfg.small_molecule_forcefield,
        molecules=molecules if molecules else None,
        forcefield_kwargs=ff_kwargs,
        periodic_forcefield_kwargs=periodic_kwargs if solvate else None,
        nonperiodic_forcefield_kwargs=None if solvate else periodic_kwargs,
        cache=None,
    )


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def build_complex(pdb_input: str, ligand_source: str, cfg: SimConfig,
                  ligand_name: str = "LIG") -> PreparedSystem:
    """Protein + ligand solvated system (the flagship screening workflow)."""
    from openmm import app, unit
    from openff.toolkit import Molecule  # noqa: F401

    prot_top, prot_pos = fix_protein(pdb_input)
    ligand = load_ligand(ligand_source, ligand_name, cfg.ligand_charge_method)

    generator = _make_system_generator(cfg, [ligand], solvate=True)

    # Merge protein and ligand into one Modeller object.
    modeller = app.Modeller(prot_top, prot_pos)
    lig_top = ligand.to_topology().to_openmm()
    lig_pos = ligand.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)

    _solvate(modeller, generator, cfg)
    system = generator.create_system(modeller.topology)
    _add_barostat(system, cfg)

    return PreparedSystem(system, modeller.topology, modeller.positions,
                          ligand_resname=_ligand_resname(modeller.topology, ligand_name),
                          n_atoms=modeller.topology.getNumAtoms())


def build_protein(pdb_input: str, cfg: SimConfig) -> PreparedSystem:
    """Apo protein solvated system (protein-only MD / apo reference)."""
    from openmm import app

    prot_top, prot_pos = fix_protein(pdb_input)
    generator = _make_system_generator(cfg, [], solvate=True)
    modeller = app.Modeller(prot_top, prot_pos)
    _solvate(modeller, generator, cfg)
    system = generator.create_system(modeller.topology)
    _add_barostat(system, cfg)
    return PreparedSystem(system, modeller.topology, modeller.positions,
                          ligand_resname=None,
                          n_atoms=modeller.topology.getNumAtoms())


def build_ligand_only(ligand_source: str, cfg: SimConfig,
                      ligand_name: str = "LIG", solvate: bool = True) -> PreparedSystem:
    """Small-molecule-only MD (e.g. conformational stability of a candidate)."""
    from openmm import app

    ligand = load_ligand(ligand_source, ligand_name, cfg.ligand_charge_method)
    generator = _make_system_generator(cfg, [ligand], solvate=solvate)
    modeller = app.Modeller(ligand.to_topology().to_openmm(),
                            ligand.conformers[0].to_openmm())
    if solvate:
        _solvate(modeller, generator, cfg)
    system = generator.create_system(modeller.topology)
    if solvate:
        _add_barostat(system, cfg)
    return PreparedSystem(system, modeller.topology, modeller.positions,
                          ligand_resname=_ligand_resname(modeller.topology, ligand_name),
                          n_atoms=modeller.topology.getNumAtoms())


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _solvate(modeller, generator, cfg: SimConfig) -> None:
    from openmm import unit
    with_ff = generator.forcefield
    modeller.addSolvent(
        with_ff,
        model=cfg.water_model,
        padding=cfg.solvent_padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_molar * unit.molar,
        positiveIon=cfg.positive_ion,
        negativeIon=cfg.negative_ion,
        neutralize=True,
    )
    log.info("Solvated system: %d atoms", modeller.topology.getNumAtoms())


def _add_barostat(system, cfg: SimConfig) -> None:
    from openmm import MonteCarloBarostat, unit
    system.addForce(
        MonteCarloBarostat(cfg.pressure_bar * unit.bar,
                           cfg.temperature_k * unit.kelvin,
                           cfg.barostat_interval)
    )


def _ligand_resname(topology, ligand_name: str) -> Optional[str]:
    # OpenFF assigns the residue name from mol.name (upper-cased, truncated).
    target = ligand_name.upper()[:3]
    for res in topology.residues():
        if res.name.upper()[:3] == target:
            return res.name
    # Fall back: the non-standard residue that is not water/ion/protein.
    standard = {"HOH", "WAT", "NA", "CL", "NA+", "CL-"}
    for res in topology.residues():
        if res.name not in standard and not _is_protein_residue(res.name):
            return res.name
    return None


_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}


def _is_protein_residue(name: str) -> bool:
    return name.upper() in _AA
