"""
E5 — quantum chemistry on ion–solvent clusters   [digest p2]

The four things the digest says a QC model must state, made explicit here:
  ① species: Li+ + one EC / one DME (gas-phase, 1:1)  and one Li+ shell cut from
     the MD liquid (E2 representative configuration);
  ② total charge +1, singlet;
  ③ several *starting placements* (carbonyl O vs ether O for EC; monodentate vs
     bidentate for DME) optimised at the same level — the lowest is reported
     and the others kept, because "初始位置会影响优化所得极小值";
  ④ environment: gas phase (no continuum); this is stated on every number.

Binding energy  ΔE_bind = E(complex) − E(Li+) − E(solvent)   (digest p2, D)
   * solvent fragment at its *own* optimised geometry ("relaxed" ΔE) and at the
     complex geometry ("vertical" / interaction energy) — both are given;
   * counterpoise (Boys–Bernardi) BSSE correction on the interaction energy.

Level: B3LYP/def2-SVP optimisation (geomeTRIC), def2-TZVP single points.
GFN2-xTB values are computed alongside for E8 (potential validation).
Electronic energies only — no ZPE / thermal corrections, no solvent.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG

HARTREE_KCAL = 627.509474
OPT_BASIS, SP_BASIS, XC = "def2-svp", "def2-tzvp", "b3lyp"
OUT = os.path.join(WORK, "e5_qc")
os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------
# starting geometries from RDKit conformers + Li placement rules
# --------------------------------------------------------------------------
def _rdkit_xyz(smiles, seed=1):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                    for i in range(m.GetNumAtoms())])
    return m, sym, xyz


def starting_structures():
    """Return {name: (symbols, coords)} for Li+–EC and Li+–DME with different Li placements."""
    starts = {}
    # --- EC: carbonyl O (index of O double-bonded to C) vs ring ether O
    m, sym, xyz = _rdkit_xyz("O=C1OCCO1")
    c_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C" and
             any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds())][0]
    o_carb = [n.GetIdx() for n in m.GetAtomWithIdx(c_idx).GetNeighbors() if n.GetSymbol() == "O" and
              m.GetBondBetweenAtoms(c_idx, n.GetIdx()).GetBondTypeAsDouble() == 2][0]
    o_eth = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetIdx() != o_carb][0]
    v = xyz[o_carb] - xyz[c_idx]; v /= np.linalg.norm(v)
    starts["EC_carbonylO"] = (sym + ["Li"], np.vstack([xyz, xyz[o_carb] + 1.9 * v]))
    v2 = xyz[o_eth] - xyz[c_idx]; v2 /= np.linalg.norm(v2)
    starts["EC_etherO"] = (sym + ["Li"], np.vstack([xyz, xyz[o_eth] + 1.9 * v2]))
    # --- DME: bidentate (Li above midpoint of the two O) vs monodentate at one O
    m, sym, xyz = _rdkit_xyz("COCCOC")
    o_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O"]
    c_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C"]
    # make a gauche (chelating) conformer: rotate so that O...O ≈ 2.8 Å by placing Li at midpoint + offset
    mid = xyz[o_idx].mean(axis=0)
    centroid = xyz[c_idx].mean(axis=0)
    away = mid - centroid; away /= (np.linalg.norm(away) + 1e-9)
    starts["DME_bidentate"] = (sym + ["Li"], np.vstack([xyz, mid + 1.6 * away]))
    o0 = o_idx[0]
    nbr = [n.GetIdx() for n in m.GetAtomWithIdx(o0).GetNeighbors()]
    v = xyz[o0] - xyz[nbr].mean(axis=0); v /= np.linalg.norm(v)
    starts["DME_monodentate"] = (sym + ["Li"], np.vstack([xyz, xyz[o0] + 1.9 * v]))
    return starts


# --------------------------------------------------------------------------
# PySCF helpers
# --------------------------------------------------------------------------
def _mol(sym, xyz, charge, basis, ghost_mask=None):
    from pyscf import gto
    atoms = []
    for i, (s, r) in enumerate(zip(sym, xyz)):
        name = f"ghost-{s}" if (ghost_mask is not None and ghost_mask[i]) else s
        atoms.append([name, tuple(r)])
    return gto.M(atom=atoms, charge=charge, spin=0, basis=basis, verbose=0)


def _energy(sym, xyz, charge, basis=SP_BASIS, ghost_mask=None):
    from pyscf import dft
    mol = _mol(sym, xyz, charge, basis, ghost_mask)
    mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = 3
    return mf.kernel()


def _optimize(sym, xyz, charge, basis=OPT_BASIS, maxsteps=150):
    from pyscf import dft
    from pyscf.geomopt.geometric_solver import optimize
    mol = _mol(sym, xyz, charge, basis)
    mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = 3
    mol_eq = optimize(mf, maxsteps=maxsteps, convergence_energy=1e-5, convergence_grms=5e-4,
                      convergence_gmax=1e-3, convergence_drms=2e-3, convergence_dmax=4e-3)
    return mol_eq.atom_coords(unit="Angstrom")


def _xtb_energy(sym, xyz, charge):
    """GFN2-xTB total energy (Hartree) via xtb-python; None if unavailable."""
    try:
        from xtb.interface import Calculator, Param
        from xtb.utils import get_method
        Z = {"H": 1, "Li": 3, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16}
        calc = Calculator(Param.GFN2xTB, np.array([Z[s] for s in sym]), np.array(xyz) / 0.52917721092,
                          charge=float(charge))
        calc.set_verbosity(0)
        return float(calc.singlepoint().get_energy())
    except Exception as e:
        return None


def binding(sym, xyz, li_index, charge=+1):
    """ΔE (relaxed and vertical) + counterpoise on the vertical interaction energy, def2-TZVP."""
    mask_li = np.zeros(len(sym), bool); mask_li[li_index] = True
    E_cplx = _energy(sym, xyz, charge)
    # fragments at complex geometry (vertical), with and without ghost functions
    E_solv_v = _energy([s for s, m in zip(sym, mask_li) if not m], xyz[~mask_li], 0)
    E_li_v = _energy(["Li"], xyz[mask_li], +1)
    E_solv_cp = _energy(sym, xyz, 0, ghost_mask=mask_li)             # solvent in full basis
    E_li_cp = _energy(sym, xyz, +1, ghost_mask=~mask_li)             # Li in full basis
    # relaxed solvent: optimise the free solvent from the complex geometry
    solv_sym = [s for s, m in zip(sym, mask_li) if not m]
    xyz_solv_opt = _optimize(solv_sym, xyz[~mask_li], 0)
    E_solv_r = _energy(solv_sym, xyz_solv_opt, 0)
    dE_vert = (E_cplx - E_li_v - E_solv_v) * HARTREE_KCAL
    dE_cp = (E_cplx - E_li_cp - E_solv_cp) * HARTREE_KCAL
    bsse = dE_cp - dE_vert
    dE_relax = (E_cplx - E_li_v - E_solv_r) * HARTREE_KCAL
    deform = (E_solv_v - E_solv_r) * HARTREE_KCAL
    xtb_c, xtb_s, xtb_li = _xtb_energy(sym, xyz, charge), _xtb_energy(solv_sym, xyz[~mask_li], 0), _xtb_energy(["Li"], xyz[mask_li], 1)
    return {"E_complex_Ha": E_cplx, "dE_vertical_kcal": dE_vert, "dE_vertical_CP_kcal": dE_cp,
            "BSSE_kcal": bsse, "dE_relaxed_kcal": dE_relax, "solvent_deformation_kcal": deform,
            "dE_relaxed_CP_kcal": dE_relax + bsse,
            "xtb_dE_vertical_kcal": ((xtb_c - xtb_s - xtb_li) * HARTREE_KCAL) if None not in (xtb_c, xtb_s, xtb_li) else None}


def _write_xyz(path, sym, xyz, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(sym)}\n{comment}\n")
        for s, r in zip(sym, xyz):
            fh.write(f"{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")


def _li_o_distances(sym, xyz, li):
    o = [i for i, s in enumerate(sym) if s == "O"]
    return sorted(float(np.linalg.norm(xyz[i] - xyz[li])) for i in o)


def run_pairs(force=False) -> dict:
    """Li+ with one EC / one DME from several starting placements."""
    cache = os.path.join(OUT, "e5_pairs.json")
    if os.path.exists(cache) and not force:
        return json.load(open(cache))
    res = {"level": f"{XC}/{OPT_BASIS} opt, {XC}/{SP_BASIS} single points; gas phase; electronic energies",
           "structures": {}}
    for name, (sym, xyz) in starting_structures().items():
        t0 = time.time()
        li = len(sym) - 1
        xyz_opt = _optimize(sym, xyz, +1)
        b = binding(sym, xyz_opt, li)
        b["Li_O_distances_A"] = _li_o_distances(sym, xyz_opt, li)
        b["start_Li_O_distances_A"] = _li_o_distances(sym, xyz, li)
        b["wall_s"] = round(time.time() - t0, 1)
        _write_xyz(os.path.join(OUT, f"{name}_opt.xyz"), sym, xyz_opt, f"{name} {XC}/{OPT_BASIS}")
        res["structures"][name] = b
        print(f"[E5] {name:16s} ΔE_relaxed(CP) = {b['dE_relaxed_CP_kcal']:7.2f} kcal/mol  "
              f"(vertical {b['dE_vertical_kcal']:7.2f}, BSSE {b['BSSE_kcal']:+.2f}, deform {b['solvent_deformation_kcal']:+.2f})  "
              f"Li–O {['%.2f' % d for d in b['Li_O_distances_A'][:2]]} Å   {b['wall_s']} s")
    for solv in ("EC", "DME"):
        names = [n for n in res["structures"] if n.startswith(solv)]
        best = min(names, key=lambda n: res["structures"][n]["E_complex_Ha"])
        res[f"{solv}_lowest"] = best
        res[f"{solv}_spread_kcal"] = float((max(res["structures"][n]["E_complex_Ha"] for n in names)
                                            - res["structures"][best]["E_complex_Ha"]) * HARTREE_KCAL)
    json.dump(res, open(cache, "w"), indent=1)
    return res


def run_md_shell(label="LiTFSI_r10", force=False) -> dict | None:
    """Single-point binding of a Li+ first shell cut from the MD liquid (E2 output),
    after a GFN2-xTB relaxation (DFT optimisation of ~50 atoms is skipped on CPU)."""
    import MDAnalysis as mda
    pdb = os.path.join(FIG, f"e2_shell_{label}.pdb")
    cache = os.path.join(OUT, f"e5_shell_{label}.json")
    if not os.path.exists(pdb):
        return None
    if os.path.exists(cache) and not force:
        return json.load(open(cache))
    u = mda.Universe(pdb)
    sym = [a.element if a.element else a.name[0] for a in u.atoms]
    sym = [s.capitalize() for s in sym]
    xyz = u.atoms.positions.astype(float)
    li = [i for i, s in enumerate(sym) if s == "Li"][0]
    n_tfsi = sum(1 for s in sym if s == "S") // 2
    charge = 1 - n_tfsi
    t0 = time.time()
    # xtb relaxation
    try:
        from xtb.interface import Calculator, Param
        from xtb.ase_calculator import XTB
        from ase import Atoms
        from ase.optimize import BFGS
        atoms = Atoms(sym, positions=xyz); atoms.set_initial_charges([0] * len(sym))
        atoms.calc = XTB(method="GFN2-xTB", charge=charge)
        BFGS(atoms, logfile=None).run(fmax=0.05, steps=300)
        xyz = atoms.get_positions()
        relaxed = "GFN2-xTB (ASE BFGS, fmax 0.05 eV/Å)"
    except Exception as e:
        relaxed = f"none ({type(e).__name__})"
    # interaction energy of Li+ with the whole shell (shell fragment charge = -n_tfsi), def2-SVP to keep CPU time sane
    mask = np.zeros(len(sym), bool); mask[li] = True
    shell_sym = [s for s, m in zip(sym, mask) if not m]
    E_c = _energy(sym, xyz, charge, basis=OPT_BASIS)
    E_s = _energy(shell_sym, xyz[~mask], charge - 1, basis=OPT_BASIS)
    E_li = _energy(["Li"], xyz[mask], 1, basis=OPT_BASIS)
    E_s_cp = _energy(sym, xyz, charge - 1, basis=OPT_BASIS, ghost_mask=mask)
    E_li_cp = _energy(sym, xyz, 1, basis=OPT_BASIS, ghost_mask=~mask)
    out = {"label": label, "n_atoms": len(sym), "n_DME": sum(1 for s in sym if s == "C") // 4 if n_tfsi == 0 else None,
           "n_TFSI": n_tfsi, "cluster_charge": charge, "relaxation": relaxed,
           "level": f"{XC}/{OPT_BASIS} single points",
           "dE_Li_shell_vertical_kcal": (E_c - E_s - E_li) * HARTREE_KCAL,
           "dE_Li_shell_vertical_CP_kcal": (E_c - E_s_cp - E_li_cp) * HARTREE_KCAL,
           "Li_O_distances_A": _li_o_distances(sym, xyz, li)[:6], "wall_s": round(time.time() - t0, 1),
           "note": "⚠️ one MD snapshot, xtb-relaxed, no thermal/solvent correction: illustrates that the "
                   "Li+–shell interaction is NOT the sum of 1:1 binding energies; not a free energy"}
    _write_xyz(os.path.join(OUT, f"shell_{label}_xtbopt.xyz"), sym, xyz)
    json.dump(out, open(cache, "w"), indent=1)
    print(f"[E5] MD shell {label}: {len(sym)} atoms, charge {charge}, ΔE(Li+…shell, CP) = "
          f"{out['dE_Li_shell_vertical_CP_kcal']:.1f} kcal/mol   {out['wall_s']} s")
    return out


def report(force=False) -> dict:
    print("E5  Li+–solvent clusters, gas phase, B3LYP/def2-TZVP//def2-SVP (electronic energies only)")
    pairs = run_pairs(force)
    for solv in ("EC", "DME"):
        print(f"    {solv}: lowest start = {pairs[f'{solv}_lowest']}, spread between starts = "
              f"{pairs[f'{solv}_spread_kcal']:.2f} kcal/mol")
    shell = run_md_shell(force=force)
    _plot(pairs, shell)
    return {"pairs": pairs, "shell": shell}


def _plot(pairs, shell):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(pairs["structures"])
    fig, ax = plt.subplots(figsize=(7, 3.8))
    x = np.arange(len(names))
    for j, (key, lab) in enumerate((("dE_vertical_kcal", "vertical"), ("dE_vertical_CP_kcal", "vertical + CP"),
                                     ("dE_relaxed_CP_kcal", "relaxed + CP"), ("xtb_dE_vertical_kcal", "GFN2-xTB vertical"))):
        y = [pairs["structures"][n].get(key) or np.nan for n in names]
        ax.bar(x + (j - 1.5) * 0.2, y, width=0.2, label=lab)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8); ax.set_ylabel("ΔE_bind (kcal/mol)")
    ax.axhline(0, color="k", lw=0.5); ax.legend(fontsize=7)
    t = "E5  Li+–solvent binding: definition matters (vertical / CP / relaxed) and so does the starting placement"
    if shell:
        t += f"\nMD shell ({shell['n_atoms']} atoms, {shell['n_TFSI']} TFSI): ΔE(Li+…shell,CP) = {shell['dE_Li_shell_vertical_CP_kcal']:.0f} kcal/mol"
    ax.set_title(t, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e5_qc_binding.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(force="--force" in sys.argv)
