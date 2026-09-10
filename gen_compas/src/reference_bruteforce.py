"""Brute-force unbiased MD reference for NANMA.

Two purposes:
  1. an unbiased (but slowly converging) FEL to cross-check metadynamics;
  2. the cost baseline -- how much plain MD is needed to observe the same
     number of C7eq <-> C7ax transitions that Gen-COMPAS produces.

Usage: python3 src/reference_bruteforce.py <replica_id> <ns>
"""

import os
import sys

import numpy as np
import openmm
from openmm import unit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def main():
    rep = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    ns = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0

    sim = common.make_simulation("CPU", threads=1)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(common.TEMPERATURE, 777 + rep)
    sim.step(50000)  # 100 ps equilibration

    save_every = 500  # 1 ps
    n_steps = int(ns * 1000 / 0.002)
    n_frames = n_steps // save_every

    phis = np.empty(n_frames, dtype=np.float32)
    psis = np.empty(n_frames, dtype=np.float32)
    # store a subset of full coordinates for later analysis (every 20 ps)
    coord_stride = 20
    coords = np.empty((n_frames // coord_stride + 1, 22, 3), dtype=np.float32)

    out = os.path.join(common.RESULTS, f"bruteforce_r{rep}")
    for i in range(n_frames):
        sim.step(save_every)
        x = sim.context.getState(getPositions=True).getPositions(
            asNumpy=True).value_in_unit(unit.nanometer)
        ph, ps = common.phi_psi(x)
        phis[i], psis[i] = ph, ps
        if i % coord_stride == 0:
            coords[i // coord_stride] = x
        if (i + 1) % 5000 == 0:
            np.savez_compressed(out + ".npz", phi=phis[:i + 1], psi=psis[:i + 1],
                                coords=coords[:i // coord_stride + 1])
            st = common.which_basin(phis[:i + 1])
            st = st[st >= 0]
            trans = int(np.sum(np.diff(st) != 0)) if len(st) > 1 else 0
            print(f"[bf r{rep}] {(i + 1) * save_every * 0.002 / 1000:7.2f} ns  "
                  f"transitions={trans}", flush=True)

    np.savez_compressed(out + ".npz", phi=phis, psi=psis, coords=coords)
    print("[bf] done", flush=True)


if __name__ == "__main__":
    main()
