"""
00_build_structures.py
Build all molecular and surface structures needed to reproduce the project
"CCB012626YW01" (RSV / C13H10N4S corrosion & binding study) with fairchem UMA.

Systems built:
  * C13H10N4S  = 3-mercapto-4-phenyl-5-(pyridin-4-yl)-4H-1,2,4-triazole
                 (the corrosion-inhibitor molecule adsorbed on Fe(110))
  * RSV        = Rosuvastatin (C22H28FN3O6S), the statin bound to TiOx clusters
  * Fe(110)    = BCC iron (110) slab
  * TiO2       = n-type stoichiometric titania cluster (H/OH saturated)
  * TiO        = p-type (oxygen-poor) titania cluster (H/OH saturated)

Geometries are written to ../structures/*.xyz (with cell info for periodic ones
saved as .traj / .json where relevant).
"""
import os
import json
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import bcc110, add_adsorbate

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
os.makedirs(SDIR, exist_ok=True)


# -----------------------------------------------------------------------------
# Helpers: build a 3D molecule from SMILES with RDKit, return an ASE Atoms
# -----------------------------------------------------------------------------
def mol_from_smiles(smiles, name, seed=42):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"RDKit could not parse SMILES for {name}: {smiles}"
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    assert AllChem.EmbedMolecule(m, params) == 0, f"embed failed for {name}"
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    atoms = Atoms(symbols=syms, positions=pos)
    formula = atoms.get_chemical_formula()
    print(f"  {name}: {formula}  ({len(atoms)} atoms)")
    return atoms, formula


def build_molecules():
    print("[molecules]")
    out = {}
    # C13H10N4S : 3-mercapto-4-phenyl-5-(pyridin-4-yl)-4H-1,2,4-triazole
    inhib_smiles = "Sc1nnc(-c2ccncc2)n1-c1ccccc1"
    inhib, f1 = mol_from_smiles(inhib_smiles, "C13H10N4S (triazole-thiol inhibitor)")
    assert f1 == "C13H10N4S", f"formula mismatch: got {f1}, expected C13H10N4S"
    write(os.path.join(SDIR, "inhibitor_C13H10N4S.xyz"), inhib)
    out["inhibitor"] = {"smiles": inhib_smiles, "formula": f1, "natoms": len(inhib),
                        "charge": 0, "spin": 1}

    # Rosuvastatin (RSV)
    rsv_smiles = ("CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1"
                  "/C=C/[C@@H](O)C[C@@H](O)CC(=O)O")
    rsv, f2 = mol_from_smiles(rsv_smiles, "Rosuvastatin (RSV)")
    write(os.path.join(SDIR, "rosuvastatin_RSV.xyz"), rsv)
    out["rsv"] = {"smiles": rsv_smiles, "formula": f2, "natoms": len(rsv),
                  "charge": 0, "spin": 1}
    return out


# -----------------------------------------------------------------------------
# Fe(110) BCC slab
# -----------------------------------------------------------------------------
def build_fe110(size=(5, 5, 4), a=2.866, vacuum=8.0):
    print("[Fe(110) slab]")
    slab = bcc110("Fe", size=size, a=a, vacuum=vacuum)
    slab.pbc = (True, True, True)
    # freeze the bottom two layers (typical surface-science protocol)
    zs = slab.positions[:, 2]
    layer_tol = 0.5
    unique_layers = sorted(set(np.round(zs / layer_tol) * layer_tol))
    bottom_cut = sorted(set(np.round(zs, 1)))[1] + 0.1  # bottom 2 layers frozen
    from ase.constraints import FixAtoms
    mask = slab.positions[:, 2] < bottom_cut
    slab.set_constraint(FixAtoms(mask=mask))
    print(f"  Fe(110): {len(slab)} atoms, cell diag {np.round(slab.cell.lengths(),2)}, "
          f"{int(mask.sum())} atoms fixed")
    write(os.path.join(SDIR, "fe110_slab.traj"), slab)
    write(os.path.join(SDIR, "fe110_slab.xyz"), slab)
    return {"natoms": len(slab), "size": list(size), "a": a, "vacuum": vacuum,
            "n_fixed": int(mask.sum())}


# -----------------------------------------------------------------------------
# TiOx clusters (cut from bulk, H/OH saturated) -> finite molecular models
# -----------------------------------------------------------------------------
def _saturate_and_neutralize(atoms):
    """Cap under-coordinated edge O with H (mimics the report's
    'hydrogen-saturated edge dangling bonds')."""
    from ase.neighborlist import natural_cutoffs, NeighborList
    cut = natural_cutoffs(atoms, mult=1.1)
    nl = NeighborList(cut, self_interaction=False, bothways=True)
    nl.update(atoms)
    add = []
    for i, at in enumerate(atoms):
        if at.symbol != "O":
            continue
        nbrs, _ = nl.get_neighbors(i)
        if len(nbrs) < 2:  # dangling / terminal oxygen -> add H
            # place H ~0.98 A away, pointing outward from cluster centroid
            v = atoms.positions[i] - atoms.get_center_of_mass()
            v = v / (np.linalg.norm(v) + 1e-9)
            add.append(atoms.positions[i] + 0.98 * v)
    from ase import Atom
    for p in add:
        atoms.append(Atom("H", position=p))
    return atoms


def build_tio2_cluster():
    print("[TiO2 cluster (n-type)]")
    from ase.spacegroup import crystal
    # anatase TiO2
    a, c = 3.7845, 9.5143
    anatase = crystal(["Ti", "O"], basis=[(0, 0, 0), (0, 0, 0.208)],
                      spacegroup=141, cellpar=[a, a, c, 90, 90, 90])
    sc = anatase * (2, 2, 1)
    com = sc.get_center_of_mass()
    d = np.linalg.norm(sc.positions - com, axis=1)
    keep = d < 4.2
    clu = sc[keep]
    clu = _saturate_and_neutralize(clu)
    clu.set_pbc(False)
    clu.center(vacuum=8.0)
    f = clu.get_chemical_formula()
    nTi = sum(1 for s in clu.symbols if s == "Ti")
    nO = sum(1 for s in clu.symbols if s == "O")
    print(f"  TiO2 cluster: {f}  Ti={nTi} O={nO} ({len(clu)} atoms)")
    write(os.path.join(SDIR, "tio2_cluster_ntype.xyz"), clu)
    return {"formula": f, "natoms": len(clu), "nTi": nTi, "nO": nO,
            "charge": 0, "spin": 1}


def build_tio_cluster():
    print("[TiO cluster (p-type)]")
    from ase.build import bulk
    # rocksalt TiO (Ti:O = 1:1, oxygen-poor relative to TiO2)
    tio = bulk("TiO", crystalstructure="rocksalt", a=4.177)
    sc = tio * (3, 3, 3)
    com = sc.get_center_of_mass()
    d = np.linalg.norm(sc.positions - com, axis=1)
    keep = d < 4.0
    clu = sc[keep]
    clu = _saturate_and_neutralize(clu)
    clu.set_pbc(False)
    clu.center(vacuum=8.0)
    f = clu.get_chemical_formula()
    nTi = sum(1 for s in clu.symbols if s == "Ti")
    nO = sum(1 for s in clu.symbols if s == "O")
    print(f"  TiO cluster: {f}  Ti={nTi} O={nO} ({len(clu)} atoms)")
    write(os.path.join(SDIR, "tio_cluster_ptype.xyz"), clu)
    return {"formula": f, "natoms": len(clu), "nTi": nTi, "nO": nO,
            "charge": 0, "spin": 1}


if __name__ == "__main__":
    meta = {}
    meta["molecules"] = build_molecules()
    meta["fe110"] = build_fe110()
    meta["tio2_ntype"] = build_tio2_cluster()
    meta["tio_ptype"] = build_tio_cluster()
    with open(os.path.join(SDIR, "structures_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print("\nAll structures written to structures/. Meta saved.")
