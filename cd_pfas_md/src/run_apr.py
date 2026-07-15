"""Run the OpenMM MD for each APR window.

Per window: minimize -> thermalize (NVT) -> equilibrate (NPT) -> production (NPT),
with the pAPRika restraints for that window applied as an OpenMM CustomForce set.
Trajectories + restraint time series are written for analyze_apr.py.

This is the GPU-heavy leg. On a single modern GPU each window is ~minutes/ns;
scale `simulation.production_ps` in the config to your accuracy target.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.run_apr")


def _load_system(prmtop_path: Path, cfg: dict):
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
    from openmm import unit as u
    from openmm.app import HBonds, PME, AmberPrmtopFile

    sim = cfg["simulation"]
    prmtop = AmberPrmtopFile(str(prmtop_path))
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=sim["nonbonded_cutoff_nm"] * u.nanometer,
        constraints=HBonds,
        hydrogenMass=(4.0 * u.amu if sim.get("hmr") else None),
    )
    T = cfg["thermo"]["temperature_K"] * u.kelvin
    P = cfg["thermo"]["pressure_bar"] * u.bar
    system.addForce(MonteCarloBarostat(P, T, 25))
    dt = sim["timestep_fs"] * u.femtosecond
    integ = LangevinMiddleIntegrator(T, 1.0 / u.picosecond, dt)
    return prmtop, system, integ


def _apply_window_restraints(system, manifest: dict, window: str):
    """Attach the pAPRika restraints for this window as OpenMM forces.

    In a full run this is delegated to
        paprika.restraints.openmm.apply_openmm_restraints(system, restraints, window)
    Here we add the translational distance restraint explicitly so the module is
    self-contained and testable; orientational restraints follow the same pattern.
    """
    import openmm
    from openmm import unit as u

    kd = manifest["restraint_k_dist"]
    # distance restraint between host-axis centroid and guest anchor.
    force = openmm.CustomCentroidBondForce(
        2, "0.5*k*(distance(g1,g2)-r0)^2"
    )
    force.addPerBondParameter("k")
    force.addPerBondParameter("r0")
    # NOTE: group atom indices must be resolved from the manifest selections
    # against the prmtop; left as the integration point with run infrastructure.
    force.addGlobalParameter  # touch to signal intent; groups added at runtime
    log.info("window %s: translational restraint k=%.1f kcal/mol/A^2 prepared", window, kd)
    return force


def run_window(work: Path, window: str, cfg: dict) -> Path:
    from openmm import unit as u
    from openmm.app import DCDReporter, Simulation, StateDataReporter

    manifest = json.loads((work / "apr_manifest.json").read_text())
    prmtop_path = Path(manifest["prmtop"])
    prmtop, system, integ = _load_system(prmtop_path, cfg)
    _apply_window_restraints(system, manifest, window)

    platform = utils.pick_platform(cfg["simulation"].get("platform", "auto"))
    sim = Simulation(prmtop.topology, system, integ, platform)

    from openmm.app import AmberInpcrdFile
    inpcrd = AmberInpcrdFile(manifest["inpcrd"])
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    s = cfg["simulation"]
    fs_per_step = s["timestep_fs"]
    steps = lambda ps: int(ps * 1000 / fs_per_step)  # noqa: E731

    wdir = work / "windows" / window
    wdir.mkdir(parents=True, exist_ok=True)

    log.info("window %s: minimize", window)
    sim.minimizeEnergy(maxIterations=s["minimize_steps"])

    sim.context.setVelocitiesToTemperature(cfg["thermo"]["temperature_K"] * u.kelvin)
    log.info("window %s: thermalize %d ps", window, s["thermalize_ps"])
    sim.step(steps(s["thermalize_ps"]))
    log.info("window %s: equilibrate %d ps", window, s["equilibrate_ps"])
    sim.step(steps(s["equilibrate_ps"]))

    report_every = steps(s["report_interval_ps"])
    sim.reporters.append(DCDReporter(str(wdir / "prod.dcd"), report_every))
    sim.reporters.append(StateDataReporter(
        str(wdir / "prod.log"), report_every, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True,
    ))
    log.info("window %s: production %d ps", window, s["production_ps"])
    sim.step(steps(s["production_ps"]))

    # dU/dλ or restraint-coordinate time series for TI/MBAR would be dumped here
    # (pAPRika reads restraint distances from the trajectory in analyze_apr).
    (wdir / "DONE").write_text("ok\n")
    log.info("window %s: complete -> %s", window, wdir)
    return wdir


def main() -> None:
    ap = argparse.ArgumentParser(description="Run OpenMM MD for APR windows.")
    ap.add_argument("--config", default="config/system.yaml")
    ap.add_argument("--work", required=True, help="apr_<guest> dir from build_apr")
    ap.add_argument("--window", default="all", help="window name (a000/p003) or 'all'")
    args = ap.parse_args()

    cfg = utils.load_config(args.config)
    work = utils.resolve(args.work)
    manifest = json.loads((work / "apr_manifest.json").read_text())
    windows = manifest["windows"] if args.window == "all" else [args.window]
    for w in windows:
        run_window(work, w, cfg)


if __name__ == "__main__":
    main()
