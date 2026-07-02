"""Simulation and workflow configuration for mdscreen.

All tunable parameters live here as dataclasses so that a run can be fully
described by a single serialisable object (handy for reproducibility and for
shipping the config across the wire to a Modal GPU worker).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class SimConfig:
    """Physical + engine parameters for an OpenMM MD run.

    Defaults mirror the choices made by the making-it-rain notebooks
    (ff14SB protein force field, GAFF2 for small molecules, TIP3P water,
    Langevin thermostat, PME electrostatics) but are all overridable.
    """

    # --- force fields -----------------------------------------------------
    protein_forcefield: str = "amber/ff14SB.xml"
    water_forcefield: str = "amber/tip3p_standard.xml"
    water_model: str = "tip3p"                      # tip3p | tip4pew | spce
    small_molecule_forcefield: str = "gaff-2.11"    # gaff-2.11 | openff-2.1.0
    ligand_charge_method: str = "am1bcc"            # am1bcc (AmberTools) | gasteiger

    # --- solvation / box --------------------------------------------------
    solvent_padding_nm: float = 1.0                 # min distance solute->box wall
    ionic_strength_molar: float = 0.15              # ~physiological NaCl
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"

    # --- integrator / thermodynamics -------------------------------------
    temperature_k: float = 300.0
    friction_per_ps: float = 1.0
    timestep_fs: float = 4.0                         # 4 fs is safe with HMR
    hydrogen_mass_amu: float = 1.5                   # hydrogen mass repartitioning
    pressure_bar: float = 1.0
    barostat_interval: int = 25

    # --- nonbonded --------------------------------------------------------
    nonbonded_cutoff_nm: float = 1.0
    constraints: str = "HBonds"                      # HBonds | AllBonds | None

    # --- protocol lengths -------------------------------------------------
    minimize_max_iterations: int = 0                 # 0 = until convergence
    nvt_equil_ps: float = 200.0                      # heating / NVT equilibration
    npt_equil_ps: float = 200.0                      # density equilibration
    production_ns: float = 5.0                        # production trajectory length

    # --- output -----------------------------------------------------------
    report_interval_ps: float = 10.0                 # trajectory + log cadence
    platform: str = "CUDA"                           # CUDA | OpenCL | CPU
    precision: str = "mixed"                          # mixed | single | double
    seed: int = 0                                     # 0 = random

    # ---------------------------------------------------------------------
    def steps(self, picoseconds: float) -> int:
        """Convert a duration in ps to an integer number of MD steps."""
        return max(1, int(round(picoseconds * 1000.0 / self.timestep_fs)))

    @property
    def report_steps(self) -> int:
        return self.steps(self.report_interval_ps)

    @property
    def production_steps(self) -> int:
        return self.steps(self.production_ns * 1000.0)

    @property
    def nvt_steps(self) -> int:
        return self.steps(self.nvt_equil_ps)

    @property
    def npt_steps(self) -> int:
        return self.steps(self.npt_equil_ps)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SimConfig":
        fields = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in fields})


@dataclass
class ScreenConfig:
    """High-level configuration for a virtual-screening / activity campaign."""

    # Receptor (a PDB path/id) shared by every ligand in the screen.
    receptor: str = ""
    # List of ligands: each entry is a SMILES string or a path to an SDF/MOL2.
    ligands: List[str] = field(default_factory=list)
    ligand_names: List[str] = field(default_factory=list)

    # Run an apo (ligand-free) trajectory so holo-vs-apo dynamics can be
    # compared -> this is what lets us flag allosteric perturbation.
    run_apo_reference: bool = True
    # Compute MM/GBSA end-state binding free energy for ranking.
    compute_binding: bool = True
    # Run dynamic cross-correlation / network allostery analysis.
    compute_allostery: bool = True

    # Screening runs are usually short; override the sim defaults for speed.
    sim: SimConfig = field(default_factory=lambda: SimConfig(production_ns=2.0))

    # Decision thresholds used to classify a ligand (kcal/mol).
    active_dg_threshold: float = -6.0      # <= this -> likely active binder
    strong_dg_threshold: float = -9.0      # <= this -> strong binder
    pose_rmsd_stable_nm: float = 0.25      # ligand RMSD below this -> stable pose

    def to_dict(self) -> dict:
        d = dataclasses.asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "ScreenConfig":
        d = dict(d)
        sim = d.pop("sim", None)
        cfg = cls(**{k: v for k, v in d.items()
                     if k in {f.name for f in dataclasses.fields(cls)}})
        if sim:
            cfg.sim = SimConfig.from_dict(sim)
        return cfg
