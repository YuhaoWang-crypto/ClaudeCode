"""OpenMM MD engine: minimise -> NVT -> NPT -> production.

Produces a topology PDB and a DCD trajectory plus a CSV state log, exactly
the artifacts the analysis layer consumes.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

from .config import SimConfig
from .prepare import PreparedSystem
from .utils import log, timer, ensure_dir


@dataclass
class SimulationResult:
    topology_pdb: str
    trajectory_dcd: str
    state_log_csv: str
    final_pdb: str
    ligand_resname: Optional[str]
    production_ns: float


def _make_integrator(cfg: SimConfig):
    from openmm import LangevinMiddleIntegrator, unit
    integ = LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        cfg.friction_per_ps / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    if cfg.seed:
        integ.setRandomNumberSeed(cfg.seed)
    return integ


def _make_simulation(prepared: PreparedSystem, cfg: SimConfig):
    from openmm import app, Platform

    integrator = _make_integrator(cfg)
    try:
        platform = Platform.getPlatformByName(cfg.platform)
        props = {}
        if cfg.platform == "CUDA":
            props = {"Precision": cfg.precision}
        simulation = app.Simulation(prepared.topology, prepared.system,
                                    integrator, platform, props)
    except Exception as exc:  # pragma: no cover - platform fallback
        log.warning("Platform %s unavailable (%s); falling back to CPU",
                    cfg.platform, exc)
        platform = Platform.getPlatformByName("CPU")
        simulation = app.Simulation(prepared.topology, prepared.system,
                                    integrator, platform)
    simulation.context.setPositions(prepared.positions)
    return simulation


def run(prepared: PreparedSystem, cfg: SimConfig, outdir: str,
        tag: str = "sim") -> SimulationResult:
    """Full protocol. Returns paths to the produced artifacts."""
    from openmm import app, unit
    from openmm import MonteCarloBarostat

    ensure_dir(outdir)
    top_pdb = os.path.join(outdir, f"{tag}_topology.pdb")
    traj_dcd = os.path.join(outdir, f"{tag}_production.dcd")
    log_csv = os.path.join(outdir, f"{tag}_state.csv")
    final_pdb = os.path.join(outdir, f"{tag}_final.pdb")

    simulation = _make_simulation(prepared, cfg)

    # Save the topology (with initial positions) for the analysis layer.
    with open(top_pdb, "w") as fh:
        app.PDBFile.writeFile(prepared.topology, prepared.positions, fh)

    with timer("energy minimisation"):
        simulation.minimizeEnergy(
            maxIterations=cfg.minimize_max_iterations or 0)

    # --- NVT: heat to target temperature (barostat disabled) --------------
    _set_barostat_frequency(simulation, 0)
    simulation.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin)
    with timer(f"NVT equilibration ({cfg.nvt_equil_ps} ps)"):
        simulation.step(cfg.nvt_steps)

    # --- NPT: equilibrate density (barostat on) ---------------------------
    _set_barostat_frequency(simulation, cfg.barostat_interval)
    with timer(f"NPT equilibration ({cfg.npt_equil_ps} ps)"):
        simulation.step(cfg.npt_steps)

    # --- Production -------------------------------------------------------
    simulation.reporters.append(app.DCDReporter(traj_dcd, cfg.report_steps))
    simulation.reporters.append(app.StateDataReporter(
        log_csv, cfg.report_steps, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True, speed=True,
    ))
    with timer(f"production ({cfg.production_ns} ns)"):
        simulation.step(cfg.production_steps)

    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, "w") as fh:
        app.PDBFile.writeFile(prepared.topology, state.getPositions(), fh)

    log.info("Simulation artifacts written to %s", outdir)
    return SimulationResult(
        topology_pdb=top_pdb, trajectory_dcd=traj_dcd, state_log_csv=log_csv,
        final_pdb=final_pdb, ligand_resname=prepared.ligand_resname,
        production_ns=cfg.production_ns,
    )


def _set_barostat_frequency(simulation, frequency: int) -> None:
    """Toggle the MonteCarloBarostat (frequency 0 == effectively off)."""
    from openmm import MonteCarloBarostat
    system = simulation.system
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            force.setFrequency(frequency)
            simulation.context.reinitialize(preserveState=True)
            return
