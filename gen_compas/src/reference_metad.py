"""Ground-truth free-energy landscape of NANMA by well-tempered metadynamics.

This is the *reference* the Gen-COMPAS result is validated against.  It biases
(phi, psi) explicitly -- exactly the kind of predefined-CV method Gen-COMPAS is
meant to avoid -- and is used here only to provide a converged answer.

Usage: python3 src/reference_metad.py <walker_id> <n_walkers> <ns>
"""

import os
import sys

import numpy as np
import openmm
from openmm import unit
from openmm.app import Simulation, PDBFile
from openmm.app.metadynamics import Metadynamics, BiasVariable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def build_torsion_force(indices):
    f = openmm.CustomTorsionForce("theta")
    f.addTorsion(*indices)
    return f


def main():
    walker = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    n_walkers = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    ns = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0

    top, system, pos, _ = common.get_system()
    m = common.atom_index_map()
    phi_idx = (m["C_0"], m["N_1"], m["CA_1"], m["C_1"])
    psi_idx = (m["N_1"], m["CA_1"], m["C_1"], m["N_2"])

    width = 0.35  # rad, gaussian width
    bv_phi = BiasVariable(build_torsion_force(phi_idx), -np.pi, np.pi,
                          width, True, gridWidth=72)
    bv_psi = BiasVariable(build_torsion_force(psi_idx), -np.pi, np.pi,
                          width, True, gridWidth=72)

    bias_dir = os.path.join(common.RESULTS, "metad_bias")
    os.makedirs(bias_dir, exist_ok=True)

    meta = Metadynamics(
        system,
        [bv_phi, bv_psi],
        common.TEMPERATURE,
        biasFactor=8.0,
        height=1.0 * unit.kilojoules_per_mole,
        frequency=250,
        saveFrequency=5000,
        biasDir=bias_dir,
    )

    integ = openmm.LangevinMiddleIntegrator(
        common.TEMPERATURE, common.FRICTION, common.TIMESTEP
    )
    plat = openmm.Platform.getPlatformByName("CPU")
    sim = Simulation(top, system, integ, plat, {"Threads": "1"})
    sim.context.setPositions(pos)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(common.TEMPERATURE, 1234 + walker)

    n_steps = int(ns * 1000 / 0.002)
    chunk = 25000
    done = 0
    while done < n_steps:
        meta.step(sim, min(chunk, n_steps - done))
        done += chunk
        fe = meta.getFreeEnergy().value_in_unit(unit.kilocalorie_per_mole)
        print(f"[metad w{walker}] {done * 0.002 / 1000:7.2f} ns  "
              f"span={fe.max() - fe.min():6.2f} kcal/mol", flush=True)

    fe = meta.getFreeEnergy().value_in_unit(unit.kilocalorie_per_mole)
    np.save(os.path.join(common.RESULTS, f"metad_fes_w{walker}.npy"), fe)
    print("[metad] done", flush=True)


if __name__ == "__main__":
    main()
