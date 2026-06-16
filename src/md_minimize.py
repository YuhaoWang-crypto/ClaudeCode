"""
md_minimize.py
--------------
In-environment sanity check (CPU, no GPU/AmberTools): build the protein-only
WT and His-pair mutant with Amber14 + GBn2 implicit solvent, energy-minimise,
and report the potential energy + max force. This does NOT compute substrate
binding (glucose/Zn need GAFF/ZAFF params we can't generate here) -- it only
confirms the engineered His-pair model is geometrically sane (no clashes, the
minimiser converges), which is the prerequisite before the full Colab MD.
"""
from __future__ import annotations
import sys
from pdbfixer import PDBFixer
from openmm.app import ForceField, Modeller, NoCutoff, HBonds
from openmm import LangevinMiddleIntegrator, Platform
from openmm.unit import kelvin, picosecond, femtosecond, kilojoule_per_mole, nanometer
import build_md_models as B


def minimise(apply_mut: bool, label: str):
    fixer = B._protein_fixer(apply_mut=apply_mut)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)
    from openmm.app import Simulation
    sim = Simulation(modeller.topology, system, integ,
                     Platform.getPlatformByName("CPU"))
    sim.context.setPositions(modeller.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=2000)
    st = sim.context.getState(getEnergy=True, getForces=True)
    e1 = st.getPotentialEnergy()
    import numpy as np
    f = st.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole / nanometer)
    fmax = float(np.linalg.norm(f, axis=1).max())
    n = system.getNumParticles()
    print(f"[{label}] atoms={n}  E_initial={e0.value_in_unit(kilojoule_per_mole):.0f}  "
          f"E_min={e1.value_in_unit(kilojoule_per_mole):.0f} kJ/mol  Fmax={fmax:.1f} kJ/mol/nm")
    return e1.value_in_unit(kilojoule_per_mole)


if __name__ == "__main__":
    sys.path.insert(0, B.HERE)
    ewt = minimise(False, "WT      ")
    emut = minimise(True, "mutHis  ")
    print(f"\nΔE(min, mutant - WT) = {emut - ewt:.0f} kJ/mol  "
          f"(whole-protein implicit-solvent energy; clash/feasibility check only,\n"
          f" NOT the substrate-binding ΔΔG — that needs the Colab MM-GBSA run).")
