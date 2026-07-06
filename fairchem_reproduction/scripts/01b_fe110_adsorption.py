"""
01b_fe110_adsorption.py  (corrected)
Adsorption energy of C13H10N4S on Fe(110) with UMA (oc20), done properly:

Fixes vs first attempt:
  * bare-slab reference is relaxed to the TRUE minimum (-684.85 eV); the naive
    6-step relaxation was stuck 11.7 eV too high, which spuriously gave -10.8 eV.
  * the adsorbate is docked in several binding orientations (S-down, N-down,
    flat) over different surface sites so it actually chemisorbs instead of
    desorbing; the lowest-energy configuration is reported.

  Eads = E(slab+ads) - E(slab) - E(gas molecule)      (all UMA/oc20)
Report (VASP/PBE): Eads = -2.88 eV.
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
FMAX, MAXSTEPS = 0.05, 250

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
def cal(): return FAIRChemCalculator(predictor, task_name="oc20")

def relax(atoms, mask=None, steps=MAXSTEPS, label="x", logdir=RDIR):
    atoms = atoms.copy()
    if mask is not None: atoms.set_constraint(FixAtoms(mask=mask))
    atoms.calc = cal()
    opt = LBFGS(atoms, logfile=os.path.join(logdir, f"opt_{label}.log"))
    opt.run(fmax=FMAX, steps=steps)
    e = float(atoms.get_potential_energy())
    fmax = float(np.sqrt((atoms.get_forces()**2).sum(1)).max())
    return atoms, e, opt.get_number_of_steps(), fmax

def orient_atom_down(mol, anchor_idx):
    """Rotate molecule so the vector (COM->anchor) points in -z (anchor at bottom)."""
    mol = mol.copy()
    v = mol.positions[anchor_idx] - mol.get_center_of_mass()
    v = v/np.linalg.norm(v)
    target = np.array([0,0,-1.0])
    axis = np.cross(v, target); s = np.linalg.norm(axis)
    if s < 1e-6:
        return mol
    axis /= s; ang = np.degrees(np.arccos(np.clip(v@target,-1,1)))
    mol.rotate(ang, axis, center=mol.get_center_of_mass())
    return mol

def dock(mol, slab, anchor_idx, site_xy, height=2.2):
    ad = orient_atom_down(mol, anchor_idx)
    top_z = slab.positions[:,2].max()
    ad.translate([site_xy[0]-ad.positions[anchor_idx,0],
                  site_xy[1]-ad.positions[anchor_idx,1],
                  top_z + height - ad.positions[anchor_idx,2]])
    return slab + ad

def main():
    res = {"reference_report_eV": -2.88, "task":"oc20", "model":"uma-s-1p1"}

    # references
    slab = read(os.path.join(SDIR, "fe110_slab_relaxed_deep.xyz"))
    slab.pbc = True
    mask = np.zeros(len(slab), bool)
    zmin = slab.positions[:,2].min()
    mask[slab.positions[:,2] < zmin + 2.3] = True   # freeze bottom ~2 layers
    slab_r, E_slab, _, _ = relax(slab, mask=mask, label="slab_ref")
    res["E_slab"] = E_slab

    mol = read(os.path.join(SDIR, "relaxed_inhibitor_in_slabcell.xyz"))
    mol.set_cell(slab.cell); mol.pbc = True; mol.center()
    mol_r, E_mol, _, _ = relax(mol, label="mol_ref")
    res["E_adsorbate"] = E_mol
    print(f"references: E_slab={E_slab:.3f}  E_mol={E_mol:.3f}")

    # binding-atom anchors on the molecule
    sym = list(mol_r.symbols)
    S = sym.index("S")
    Ns = [i for i,s in enumerate(sym) if s=="N"]
    # a top Fe site near the cell centre
    top = slab_r.positions[slab_r.positions[:,2] > slab_r.positions[:,2].max()-0.3]
    centre = slab_r.cell.diagonal()[:2]/2
    site = top[np.argmin(np.linalg.norm(top[:,:2]-centre,axis=1))][:2]
    # a second site shifted by ~1.4 A (bridge-ish)
    site2 = site + np.array([1.2, 0.0])

    configs = [
        ("S_down_top",  S,        site,  2.2),
        ("N_down_top",  Ns[0],    site,  2.2),
        ("N2_down_brg", Ns[1],    site2, 2.2),
        ("S_down_close",S,        site,  1.9),
    ]
    best = None
    res["configs"] = {}
    for name, anchor, s_xy, h in configs:
        comp = dock(mol_r, slab_r, anchor, s_xy, h)
        cmask = np.zeros(len(comp), bool); cmask[:len(slab_r)] = mask
        t0 = time.time()
        comp_r, E_comp, steps, fm = relax(comp, mask=cmask, label=f"cfg_{name}")
        Eads = E_comp - E_slab - E_mol
        # min adsorbate-surface distance
        fe = comp_r[[i for i,s in enumerate(comp_r.symbols) if s=="Fe"]]
        ads = comp_r[[i for i,s in enumerate(comp_r.symbols) if s!="Fe"]]
        dmin = min(np.linalg.norm(fe.positions-ads.positions[j],axis=1).min()
                   for j,s in enumerate(ads.symbols) if s in ("S","N"))
        res["configs"][name] = {"Eads_eV": float(Eads), "steps": steps,
                                "fmax": fm, "min_SN_Fe_dist": float(dmin),
                                "time_s": round(time.time()-t0)}
        print(f"  [{name}] Eads={Eads:.3f} eV  min(S/N-Fe)={dmin:.2f} A  "
              f"steps={steps} {time.time()-t0:.0f}s")
        if best is None or Eads < best[1]:
            best = (name, Eads, comp_r)
    write(os.path.join(SDIR, "relaxed_fe110_best.xyz"), best[2])
    res["best_config"] = best[0]
    res["Eads_eV"] = float(best[1])
    res["Eads_diff_vs_report_eV"] = float(best[1] - (-2.88))
    print(f"\nBEST: {best[0]}  Eads = {best[1]:.3f} eV  (report -2.88 eV)")
    with open(os.path.join(RDIR, "fe110_adsorption.json"), "w") as fh:
        json.dump(res, fh, indent=2)

if __name__ == "__main__":
    main()
