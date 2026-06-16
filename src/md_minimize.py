"""
md_minimize.py
--------------
In-environment sanity check (CPU, no GPU/AmberTools): build the protein-only
WT and His-pair mutant with Amber14 + GBn2 implicit solvent, energy-minimise,
and report the potential energy + max force. This does NOT compute substrate
binding (glucose/Zn need GAFF/ZAFF params we can't generate here) -- it only
confirms the engineered His-pair model is geometrically sane (no clashes, the
minimiser converges), which is the prerequisite before the full Colab MD.

Usage:
    python src/md_minimize.py [maxiter]   # default 500
"""
from __future__ import annotations
import sys
import time
import numpy as np
from openmm.app import ForceField, Modeller, NoCutoff, HBonds, Simulation
from openmm import LangevinMiddleIntegrator, Platform
from openmm.unit import kelvin, picosecond, femtosecond, kilojoule_per_mole, nanometer
import build_md_models as B

KJ = kilojoule_per_mole


def minimise(apply_mut: bool, label: str, maxiter: int):
    t0 = time.time()
    fixer = B._protein_fixer(apply_mut=apply_mut)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)
    sim = Simulation(modeller.topology, system, integ, Platform.getPlatformByName("CPU"))
    sim.context.setPositions(modeller.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    print(f"[{label}] built {system.getNumParticles()} atoms, E_initial={e0:.0f} kJ/mol "
          f"({time.time()-t0:.0f}s) -- minimising...", flush=True)
    sim.minimizeEnergy(maxIterations=maxiter)
    st = sim.context.getState(getEnergy=True, getForces=True)
    e1 = st.getPotentialEnergy().value_in_unit(KJ)
    f = st.getForces(asNumpy=True).value_in_unit(KJ / nanometer)
    fmax = float(np.linalg.norm(f, axis=1).max())
    print(f"[{label}] E_min={e1:.0f} kJ/mol  Fmax={fmax:.1f} kJ/mol/nm  "
          f"({time.time()-t0:.0f}s total)", flush=True)
    return e1


if __name__ == "__main__":
    sys.path.insert(0, B.HERE)
    maxiter = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    ewt = minimise(False, "WT     ", maxiter)
    emut = minimise(True, "mutHis ", maxiter)
    print(f"\nDelta_E(min, mutant - WT) = {emut - ewt:.0f} kJ/mol  "
          f"(whole-protein implicit-solvent energy; clash/feasibility check only,\n"
          f" NOT the substrate-binding ddG -- that needs the Colab MM-GBSA run).", flush=True)
