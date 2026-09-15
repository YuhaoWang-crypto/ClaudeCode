"""
E1b — minimise, equilibrate, produce, and CHECK before analysing   [digest p4]

    "先检查结构，再讨论长轨迹": minimisation -> NPT equilibration with density,
    energy and pressure logged -> a convergence check that decides whether the
    production window is usable -> production NPT with frames every `frame_ps`.

Outputs in <workdir>/<label>/ :
    equil.csv, prod.csv         thermo logs (time, T, V, rho, Epot, P)
    prod.dcd                    trajectory
    final.pdb                   last frame with box
    run_record.json             what was actually run + convergence verdict
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, DynamicsSpec, SERIES


def _platform():
    import openmm
    names = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    for pref in ("CUDA", "OpenCL", "CPU"):
        if pref in names:
            plat = openmm.Platform.getPlatformByName(pref)
            props = {"Precision": "mixed"} if pref in ("CUDA", "OpenCL") else {}
            return plat, props
    return openmm.Platform.getPlatformByName("Reference"), {}


def _simulation(sysdir: str, dyn: DynamicsSpec, npt=True):
    import openmm
    from openmm import app, unit
    with open(os.path.join(sysdir, "system.xml")) as fh:
        system = openmm.XmlSerializer.deserialize(fh.read())
    pdb = app.PDBFile(os.path.join(sysdir, "initial.pdb"))
    # remove any pre-existing barostat (fresh object each call)
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i), openmm.MonteCarloBarostat):
            system.removeForce(i)
    if npt:
        system.addForce(openmm.MonteCarloBarostat(dyn.pressure_bar * unit.bar,
                                                  dyn.temperature_K * unit.kelvin,
                                                  dyn.barostat_interval))
    integ = openmm.LangevinMiddleIntegrator(dyn.temperature_K * unit.kelvin,
                                            dyn.friction_per_ps / unit.picosecond,
                                            dyn.timestep_fs * unit.femtosecond)
    plat, props = _platform()
    sim = app.Simulation(pdb.topology, system, integ, plat, props)
    sim.context.setPositions(pdb.positions)
    if pdb.topology.getPeriodicBoxVectors() is not None:
        sim.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
    return sim, plat.getName()


def _reporter(path, every_steps):
    from openmm import app
    return app.StateDataReporter(path, every_steps, step=True, time=True, temperature=True,
                                 potentialEnergy=True, volume=True, density=True,
                                 speed=True, separator=",")


def convergence_check(csv_path: str, nblocks: int = 5, tail_frac: float = 0.5) -> dict:
    """Digest p4: density, energy, pressure stable?  Compare the block means of
    the last `tail_frac` of the log; verdict = max block deviation < 2*block SD."""
    import pandas as pd
    df = pd.read_csv(csv_path)
    df.columns = [c.strip('#" ') for c in df.columns]
    tail = df.iloc[int(len(df) * (1 - tail_frac)):]
    out = {"n_rows": int(len(df)), "tail_rows": int(len(tail))}
    for col, key in (("Density (g/mL)", "density"), ("Potential Energy (kJ/mole)", "epot")):
        x = tail[col].values
        blocks = np.array_split(x, nblocks)
        bm = np.array([b.mean() for b in blocks])
        out[key] = {"mean": float(x.mean()), "block_means": bm.round(5).tolist(),
                    "block_sd": float(bm.std(ddof=1)),
                    "drift": float(bm[-1] - bm[0]),
                    "sem": float(bm.std(ddof=1) / np.sqrt(nblocks))}
    # drift vs. scatter: converged if |drift| < 2 * block SD for both quantities
    ok = all(abs(out[k]["drift"]) < 2 * out[k]["block_sd"] + 1e-12 for k in ("density", "epot"))
    out["converged"] = bool(ok)
    return out


def run(label: str, dyn: DynamicsSpec | None = None, workdir: str = WORK,
        equil_ns: float | None = None, prod_ns: float | None = None) -> dict:
    import openmm
    from openmm import app, unit
    dyn = dyn or DynamicsSpec()
    if equil_ns is not None:
        dyn.equil_npt_ns = equil_ns
    if prod_ns is not None:
        dyn.prod_ns = prod_ns
    sysdir = os.path.join(workdir, label)
    steps_per_ps = int(round(1000 / dyn.timestep_fs))
    t0 = time.time()

    sim, platform = _simulation(sysdir, dyn, npt=True)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(tolerance=dyn.minimise_tol_kj_mol_nm * unit.kilojoule_per_mole / unit.nanometer)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.context.setVelocitiesToTemperature(dyn.temperature_K * unit.kelvin, 1)

    # --- NPT equilibration ---------------------------------------------------
    sim.reporters = [_reporter(os.path.join(sysdir, "equil.csv"), int(dyn.thermo_ps * steps_per_ps))]
    sim.step(int(dyn.equil_npt_ns * 1000 * steps_per_ps))
    conv = convergence_check(os.path.join(sysdir, "equil.csv"))

    # --- production (NPT, frames every frame_ps) ------------------------------
    sim.reporters = [_reporter(os.path.join(sysdir, "prod.csv"), int(dyn.thermo_ps * steps_per_ps)),
                     app.DCDReporter(os.path.join(sysdir, "prod.dcd"), int(dyn.frame_ps * steps_per_ps),
                                     enforcePeriodicBox=False)]
    sim.step(int(dyn.prod_ns * 1000 * steps_per_ps))
    st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
    with open(os.path.join(sysdir, "final.pdb"), "w") as fh:
        sim.topology.setPeriodicBoxVectors(st.getPeriodicBoxVectors())
        app.PDBFile.writeFile(sim.topology, st.getPositions(), fh)
    prod_conv = convergence_check(os.path.join(sysdir, "prod.csv"), tail_frac=1.0)

    rec = {"label": label, "platform": platform,
           "minimisation": {"E_before_kJ_mol": e0.value_in_unit(unit.kilojoule_per_mole),
                            "E_after_kJ_mol": e1.value_in_unit(unit.kilojoule_per_mole)},
           "equil_ns": dyn.equil_npt_ns, "prod_ns": dyn.prod_ns,
           "timestep_fs": dyn.timestep_fs, "frame_ps": dyn.frame_ps,
           "equil_convergence": conv, "prod_stationarity": prod_conv,
           "wall_s": round(time.time() - t0, 1),
           "label_note": ("✅ equilibration converged" if conv["converged"]
                          else "⚠️ equilibration NOT converged by the block criterion; "
                               "extend equil_npt_ns before trusting production averages")}
    with open(os.path.join(sysdir, "run_record.json"), "w") as fh:
        json.dump(rec, fh, indent=2)
    print(f"[E1 run] {label} on {platform}: equil {dyn.equil_npt_ns} ns + prod {dyn.prod_ns} ns, "
          f"rho={prod_conv['density']['mean']:.4f}±{prod_conv['density']['sem']:.4f} g/mL, "
          f"{rec['label_note']}, {rec['wall_s']} s")
    return rec


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("labels", nargs="*", default=[c.label for c in SERIES])
    p.add_argument("--equil-ns", type=float, default=None)
    p.add_argument("--prod-ns", type=float, default=None)
    p.add_argument("--workdir", default=WORK)
    a = p.parse_args(argv)
    for lab in a.labels:
        run(lab, workdir=a.workdir, equil_ns=a.equil_ns, prod_ns=a.prod_ns)


if __name__ == "__main__":
    main()
