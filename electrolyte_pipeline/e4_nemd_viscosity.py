"""
E4b — shear viscosity by the periodic-perturbation NEMD method (Hess 2002)

OpenMM does not expose the pressure tensor, so Green–Kubo is not available
here; instead a spatially periodic acceleration  a_x(z) = A cos(k z),
k = 2π/L_z, is applied and the steady-state velocity-profile amplitude V gives

        η = A ρ / (V k²)

Implementation notes (the digest asks for "驱动条件" to be reported):
  * force on atom i  = m_i A cos(k z_i), set as a per-particle parameter of a
    CustomExternalForce and refreshed every `refresh_steps` from the current z
    (a z-dependent energy expression would add spurious z-forces);
  * NVT at the E1 equilibrated box, Nosé–Hoover thermostat (a Langevin
    thermostat would damp the flow itself and bias η upward);
  * two amplitudes are run to check linear response; V is measured from the
    mass-weighted projection  V = 2 Σ m_i v_x,i cos(k z_i) / Σ m_i, averaged
    after a settling period, with block errors.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, DynamicsSpec

AMU = 1.66053906660e-27


def run(label: str, workdir: str = WORK, amplitudes=(0.01, 0.02), ns: float = 2.0,
        refresh_steps: int = 20, settle_frac: float = 0.25, nblocks: int = 5) -> dict:
    """amplitudes in nm/ps² (0.01 nm/ps² = 1e10 m/s², typical for these boxes)."""
    import openmm
    from openmm import app, unit
    dyn = DynamicsSpec()
    sysdir = os.path.join(workdir, label)
    with open(os.path.join(sysdir, "system.xml")) as fh:
        system = openmm.XmlSerializer.deserialize(fh.read())
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i), openmm.MonteCarloBarostat):
            system.removeForce(i)
    pdb = app.PDBFile(os.path.join(sysdir, "final.pdb"))
    n = system.getNumParticles()
    mass = np.array([system.getParticleMass(i).value_in_unit(unit.amu) for i in range(n)])
    dof = 3 * n - system.getNumConstraints() - 3          # kinetic temperature needs the constrained DOF count
    ext = openmm.CustomExternalForce("-fx*x")
    ext.addPerParticleParameter("fx")
    for i in range(n):
        ext.addParticle(i, [0.0])
    system.addForce(ext)
    from electrolyte_pipeline.e1_run import _platform
    plat, props = _platform()
    results = {"label": label, "method": "periodic-perturbation NEMD, Nosé–Hoover NVT",
               "refresh_steps": refresh_steps, "ns_per_amplitude": ns, "runs": []}
    for A in amplitudes:
        integ = openmm.NoseHooverIntegrator(dyn.temperature_K * unit.kelvin, 25 / unit.picosecond,
                                            dyn.timestep_fs * unit.femtosecond)
        sim = app.Simulation(pdb.topology, system, integ, plat, props)
        sim.context.setPositions(pdb.positions)
        sim.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
        sim.context.setVelocitiesToTemperature(dyn.temperature_K * unit.kelvin, 7)
        Lz = pdb.topology.getPeriodicBoxVectors()[2][2].value_in_unit(unit.nanometer)
        k = 2 * np.pi / Lz
        nsteps = int(ns * 1000 * 1000 / dyn.timestep_fs)
        ncycles = nsteps // refresh_steps
        Vs, Ts = [], []
        t0 = time.time()
        for c in range(ncycles):
            st = sim.context.getState(getPositions=True, getVelocities=True)
            pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            vel = st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            cz = np.cos(k * pos[:, 2])
            fx = mass * A * cz                          # amu·nm/ps² = kJ/mol/nm
            for i in range(n):
                ext.setParticleParameters(i, i, [fx[i]])
            ext.updateParametersInContext(sim.context)
            Vs.append(2 * np.sum(mass * vel[:, 0] * cz) / mass.sum())
            Ts.append(np.sum(mass * np.sum(vel**2, axis=1)) / dof * AMU * 1e6 / 1.380649e-23)
            sim.step(refresh_steps)
        Vs = np.array(Vs); Ts = np.array(Ts)
        tail = Vs[int(len(Vs) * settle_frac):]
        bm = np.array([b.mean() for b in np.array_split(tail, nblocks)])
        V, Vsem = bm.mean(), bm.std(ddof=1) / np.sqrt(nblocks)
        vol = np.prod([pdb.topology.getPeriodicBoxVectors()[i][i].value_in_unit(unit.nanometer) for i in range(3)])
        rho = mass.sum() * AMU / (vol * 1e-27)                                   # kg/m³
        eta = A * 1e9 / 1e-24 * rho / (V * 1e9 / 1e-12 * (k * 1e9) ** 2)          # Pa·s   (A: nm/ps² -> m/s²; V: nm/ps -> m/s; k: 1/nm -> 1/m)
        eta_sem = eta * Vsem / abs(V)
        results["runs"].append({"A_nm_ps2": A, "k_per_nm": k, "V_nm_ps": float(V), "V_sem": float(Vsem),
                                "T_mean_K": float(Ts[int(len(Ts) * settle_frac):].mean()),
                                "eta_mPa_s": eta * 1e3, "eta_sem_mPa_s": eta_sem * 1e3,
                                "wall_s": round(time.time() - t0, 1),
                                "profile_trace": Vs[:: max(1, len(Vs) // 400)].round(5).tolist()})
        del sim
    etas = np.array([r["eta_mPa_s"] for r in results["runs"]])
    sems = np.array([r["eta_sem_mPa_s"] for r in results["runs"]])
    spread = float(abs(etas[0] - etas[-1]) / etas.mean()) if len(etas) > 1 else 0.0
    results["eta_mPa_s"] = float(etas.mean())
    results["eta_sem_mPa_s"] = float(np.sqrt(np.sum(sems**2)) / len(sems))
    results["linearity_spread"] = spread
    results["label_note"] = ("✅ linear response (amplitudes agree within 15%)" if spread < 0.15
                             else f"⚠️ amplitudes differ by {spread*100:.0f}% — non-linear or unconverged; "
                                  f"run longer / smaller A")
    with open(os.path.join(sysdir, "e4_nemd.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"[E4 NEMD] {label}: eta = {results['eta_mPa_s']:.3f} ± {results['eta_sem_mPa_s']:.3f} mPa·s "
          f"({results['label_note']})")
    return results


if __name__ == "__main__":
    import sys
    run(sys.argv[1], ns=float(sys.argv[2]) if len(sys.argv) > 2 else 2.0)
