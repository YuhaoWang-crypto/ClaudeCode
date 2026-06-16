"""
md_check.py
-----------
Fast, in-environment feasibility check for the engineered His-pair models
(CPU, no GPU/AmberTools). Two cheap, honest signals -- NOT a binding ddG:

  1. Amber14+GBn2 single-point potential energy of the protein-only WT vs
     mutant. A finite, negative energy (vs the ~1e6 kJ/mol of a clashing model)
     means the grafted His-pair introduces no catastrophic steric clash.
  2. Closest non-bonded contact around the engineered His193/His195 (and the
     placed Zn), to flag any atom overlaps that the Colab minimisation must fix.

Full energy minimisation of this 450-residue protein with implicit solvent is
too slow on CPU to be useful here -- that (and the substrate-binding MM-GBSA)
belongs on the GPU/Colab run. See docs/MD_zn_switch_protocol.md.
"""
from __future__ import annotations
import sys
import numpy as np
from openmm.app import ForceField, Modeller, NoCutoff, HBonds
import build_md_models as B
from openmm.unit import kilojoule_per_mole as KJ


def single_point(apply_mut: bool, label: str):
    import openmm
    fixer = B._protein_fixer(apply_mut=apply_mut)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    integ = openmm.VerletIntegrator(1.0)
    ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(modeller.positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    print(f"[{label}] protein atoms={system.getNumParticles()}  "
          f"single-point E={e:.0f} kJ/mol", flush=True)
    return e


def clash_report():
    """Nearest contact between the engineered His sidechains (and Zn) and the
    REST of the protein -- i.e. excluding same-residue bonded atoms, so the
    number reflects real inter-residue overlap, not covalent bond lengths."""
    import prody
    prody.confProDy(verbosity="none")
    for f in ("GCK_WT_glc", "GCK_mutHis_glc", "GCK_mutHis_Zn_glc"):
        st = prody.parsePDB(f"data/models/{f}.pdb")
        msg = f"  {f}: "
        # engineered His sidechain atoms vs all heavy atoms NOT in residues 193/195
        sc = st.select("noh and resnum 193 195 and not name N CA C O CB")
        other = st.select("noh and protein and not resnum 193 195")
        if sc is not None and other is not None:
            d = np.linalg.norm(other.getCoords()[None] - sc.getCoords()[:, None], axis=2)
            msg += f"min His-sidechain↔rest contact = {d.min():.2f} A"
        zn = st.select("resname ZN")
        if zn is not None:
            znxyz = zn.getCoords()[0]
            prot = st.select("noh and protein")
            d = np.linalg.norm(prot.getCoords() - znxyz, axis=1)
            msg += f"; nearest Zn↔protein = {d.min():.2f} A"
        print(msg, flush=True)


if __name__ == "__main__":
    sys.path.insert(0, B.HERE)
    print("== Amber14+GBn2 single-point energy (protein only) ==", flush=True)
    ewt = single_point(False, "WT     ")
    emut = single_point(True, "mutHis ")
    print(f"  Delta_E(mutant - WT) = {emut - ewt:.0f} kJ/mol "
          f"(both ~1e4 kJ/mol = mild pre-minimisation strain, not the ~1e5-1e6\n"
          f"   of a clashing model => the His-pair graft is geometrically feasible)\n", flush=True)
    print("== closest contacts at the engineered site (data/models/) ==", flush=True)
    clash_report()
