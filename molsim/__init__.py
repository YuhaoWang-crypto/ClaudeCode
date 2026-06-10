"""
molsim — Lightweight Molecular Simulation & Docking Toolkit
============================================================

Quick start
-----------
>>> from molsim import read_pdb, DockingEngine, analyze
>>> receptor = read_pdb("receptor.pdb")
>>> ligand = read_pdb("ligand.pdb")
>>> engine = DockingEngine(receptor, ligand)
>>> engine.set_binding_site(center=[x, y, z], radius=12)
>>> results = engine.dock(n_random=2000, n_refine=20, verbose=True)

Key modules
-----------
- io         : PDB / MOL2 / SDF file I/O
- forcefield : Atom type assignment & Gasteiger partial charges
- energy     : Non-bonded energy calculation & L-BFGS-B minimization
- docking    : Monte Carlo docking with Vina-style scoring
- analysis   : RMSD, RMSF, Rg, Ramachandran, secondary structure, H-bonds
- viz        : Matplotlib-based plots
"""

from .core import Atom, Bond, Molecule
from .io import read_pdb, write_pdb, read_mol2, read_sdf
from .forcefield import assign_atom_types, gasteiger_charges
from .energy import EnergyCalculator, interaction_energy, minimize_energy
from .docking import DockingEngine, DockingResult
from .analysis import (
    rmsd, aligned_rmsd, rmsf, radius_of_gyration,
    ramachandran_angles, secondary_structure, distance_matrix,
    find_hbonds, find_hydrophobic_contacts, residue_contacts,
    Trajectory,
)
from . import viz

__version__ = "0.1.0"

__all__ = [
    # Core
    "Atom", "Bond", "Molecule",
    # I/O
    "read_pdb", "write_pdb", "read_mol2", "read_sdf",
    # Force field
    "assign_atom_types", "gasteiger_charges",
    # Energy
    "EnergyCalculator", "interaction_energy", "minimize_energy",
    # Docking
    "DockingEngine", "DockingResult",
    # Analysis
    "rmsd", "aligned_rmsd", "rmsf", "radius_of_gyration",
    "ramachandran_angles", "secondary_structure", "distance_matrix",
    "find_hbonds", "find_hydrophobic_contacts", "residue_contacts",
    "Trajectory",
    # Viz
    "viz",
]
