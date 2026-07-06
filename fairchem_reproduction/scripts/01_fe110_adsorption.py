"""
01_fe110_adsorption.py
Reproduce the VASP adsorption-energy deliverable with fairchem UMA (OC20 task).

  Eads = E(slab+adsorbate) - E(slab) - E(adsorbate)

Original report value (VASP, PBE): Eads(C13H10N4S / Fe(110)) = -2.88 eV.

All three energies are computed with the SAME UMA task ("oc20") so they share one
energy reference -> the difference is a physically meaningful adsorption energy.
The isolated molecule is evaluated in the same periodic cell (vacuum) as the slab
to keep the reference identical.
"""
import os, json, time
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from fairchem.core import pretrained_mlip, FAIRChemCalculator

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
os.makedirs(RDIR, exist_ok=True)

FMAX = 0.05
MAXSTEPS = 250

print("Loading UMA (uma-s-1p1, oc20 task)...")
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
def oc20_calc():
    return FAIRChemCalculator(predictor, task_name="oc20")


def relax(atoms, label, fix_mask=None):
    atoms = atoms.copy()
    if fix_mask is not None:
        atoms.set_constraint(FixAtoms(mask=fix_mask))
    atoms.calc = oc20_calc()
    t0 = time.time()
    e0 = atoms.get_potential_energy()
    opt = LBFGS(atoms, logfile=os.path.join(RDIR, f"opt_{label}.log"))
    opt.run(fmax=FMAX, steps=MAXSTEPS)
    e = atoms.get_potential_energy()
    # FixAtoms zeroes forces on constrained atoms, so a global max is correct
    fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1)).max())
    dt = time.time() - t0
    write(os.path.join(SDIR, f"relaxed_{label}.xyz"), atoms)
    print(f"  [{label}] E0={e0:.3f} -> E={e:.3f} eV, steps={opt.get_number_of_steps()}, "
          f"fmax={fmax:.3f}, {dt:.0f}s ({len(atoms)} atoms)")
    return atoms, float(e), int(opt.get_number_of_steps()), float(fmax), float(dt)


def main():
    res = {"reference_report_eV": -2.88, "task": "oc20", "model": "uma-s-1p1",
           "fmax": FMAX, "maxsteps": MAXSTEPS}

    # ---- 1. bare Fe(110) slab (bottom layers already frozen in the .traj) ----
    slab = read(os.path.join(SDIR, "fe110_slab.traj"))
    fix_mask = np.zeros(len(slab), bool)
    if slab.constraints:
        fix_mask[slab.constraints[0].get_indices()] = True
    slab_r, E_slab, s1, f1, t1 = relax(slab, "fe110_slab", fix_mask=fix_mask)
    res["E_slab"] = E_slab

    # ---- 2. isolated adsorbate in the SAME cell (consistent reference) -------
    mol = read(os.path.join(SDIR, "inhibitor_C13H10N4S.xyz"))
    mol.set_cell(slab.cell); mol.set_pbc(True)
    mol.center()
    mol_r, E_mol, s2, f2, t2 = relax(mol, "inhibitor_in_slabcell")
    res["E_adsorbate"] = E_mol

    # ---- 3. adsorption complex: place molecule flat above top Fe layer ------
    top_z = slab_r.positions[~fix_mask][:, 2].max()
    ad = mol_r.copy()
    # lay the molecule roughly flat (its principal plane parallel to surface)
    ad.translate(-ad.get_center_of_mass())
    # orient: align smallest principal axis with z
    inertia = np.cov((ad.positions).T)
    evals, evecs = np.linalg.eigh(inertia)
    R = evecs[:, ::-1]  # largest variance -> x,y ; smallest -> z
    ad.set_positions(ad.positions @ R)
    # drop S atom closest to surface
    ad.translate([slab_r.cell[0][0]/2 + slab_r.cell[1][0]/2,
                  slab_r.cell[1][1]/2, 0])
    zmin = ad.positions[:, 2].min()
    height = 2.3
    ad.translate([0, 0, top_z + height - zmin])

    complex_atoms = slab_r.copy()
    complex_atoms += ad
    cmask = np.zeros(len(complex_atoms), bool)
    cmask[:len(slab_r)] = fix_mask  # keep same frozen bottom layers
    comp_r, E_comp, s3, f3, t3 = relax(complex_atoms, "fe110_inhibitor", fix_mask=cmask)
    res["E_complex"] = E_comp

    Eads = E_comp - E_slab - E_mol
    res["Eads_eV"] = float(Eads)
    res["Eads_diff_vs_report_eV"] = float(Eads - (-2.88))
    res["timing_s"] = {"slab": t1, "mol": t2, "complex": t3}
    res["steps"] = {"slab": s1, "mol": s2, "complex": s3}

    print("\n=== Fe(110) adsorption ===")
    print(f"  E_slab      = {E_slab:.3f} eV")
    print(f"  E_adsorbate = {E_mol:.3f} eV")
    print(f"  E_complex   = {E_comp:.3f} eV")
    print(f"  Eads (UMA)  = {Eads:.3f} eV   |   report (VASP) = -2.88 eV")

    with open(os.path.join(RDIR, "fe110_adsorption.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print("  -> results/fe110_adsorption.json")


if __name__ == "__main__":
    main()
