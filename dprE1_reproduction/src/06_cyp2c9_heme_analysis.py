#!/usr/bin/env python3
"""
Part 2b - CYP2C9 heme-interaction analysis (tests the paper's SAFETY claim).

The paper's central selectivity argument is *mechanistic*, not score-based:
TCA007 sits ~6.6 A from the heme iron (a water bridges to the Fe), so it inhibits
CYP2C9 reversibly WITHOUT coordinating the heme; and the de novo GTD compounds
"do not interact with the porphyrin ring". A raw docking score cannot show this,
so here we measure the geometry directly.

For each ligand docked into CYP2C9 (5W0C) we report:
  * min distance from any ligand atom to the heme Fe        (coordination test)
  * min distance from any ligand atom to a porphyrin-core atom
  * whether any atom is within 2.5 A of Fe  (= direct heme coordination)

Validation anchor: TCA007 is the native co-crystal ligand (9W6) in 5W0C, so its
CRYSTAL pose gives a ground-truth Fe distance to check against the paper's 6.6 A.
"""
import sys
import json
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "src"))
from importlib import import_module
d = import_module("04_docking")
OUT = BASE / "results"
COORD_CUTOFF = 2.5   # A: metal-coordination threshold


def heme_atoms(pdb="5W0C"):
    """Return Fe coord and porphyrin-core atom coords (pyrrole N + meso/core C)."""
    fe, core = None, []
    for line in open(BASE / "data" / "pdb" / f"{pdb}.pdb"):
        if line[:6] not in ("ATOM  ", "HETATM") or line[21] != "A":
            continue
        if line[17:20].strip() != "HEM":
            continue
        atom = line[12:16].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if atom == "FE":
            fe = xyz
        # porphyrin core: 4 pyrrole N + the alpha/meso carbons (exclude tails)
        elif atom.startswith(("NA", "NB", "NC", "ND")) or atom.startswith(
                ("C1", "C2", "C3", "C4", "CH")):
            core.append(xyz)
    return fe, np.array(core)


def native_ligand_coords(pdb="5W0C", resn="9W6"):
    xyz = []
    for line in open(BASE / "data" / "pdb" / f"{pdb}.pdb"):
        if line[:6] in ("ATOM  ", "HETATM") and line[21] == "A" \
                and line[17:20].strip() == resn:
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def pose_coords(pose_pdbqt):
    xyz = []
    for line in open(pose_pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def geom(lig_xyz, fe, core):
    d_fe = float(np.min(np.linalg.norm(lig_xyz - fe, axis=1)))
    d_core = float(np.min(np.linalg.norm(
        lig_xyz[:, None, :] - core[None, :, :], axis=2)))
    coordinates_heme = d_fe <= COORD_CUTOFF
    return d_fe, d_core, coordinates_heme


def main():
    from vina import Vina
    d.clean_structures()
    receptor_pdbqt, center = d.prep_receptor("5W0C")
    fe, core = heme_atoms("5W0C")

    results = {}

    # ground-truth: TCA007 native crystal pose
    d_fe, d_core, coord = geom(native_ligand_coords(), fe, core)
    results["TCA007_native(9W6)"] = dict(
        min_dist_to_Fe=round(d_fe, 2), min_dist_to_porphyrin=round(d_core, 2),
        coordinates_heme=coord, source="crystal")
    print(f"TCA007 (crystal 9W6): Fe {d_fe:.2f} A  porphyrin {d_core:.2f} A  "
          f"coordinates_heme={coord}   [paper: ~6.6 A, no coordination]")

    # docked GTD series + TCA1
    for name, smi in d.CAND.items():
        lp = d.PREP / f"5W0C_{name}_heme.pdbqt"
        if not d.prep_ligand(name, smi, lp):
            print(f"{name}: ligand prep failed"); continue
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(str(lp))
        v.compute_vina_maps(center=list(center), box_size=[22, 22, 22])
        v.dock(exhaustiveness=10, n_poses=5)
        pose = str(lp).replace(".pdbqt", "_pose.pdbqt")
        v.write_poses(pose, n_poses=1, overwrite=True)
        d_fe, d_core, coord = geom(pose_coords(pose), fe, core)
        results[name] = dict(min_dist_to_Fe=round(d_fe, 2),
                             min_dist_to_porphyrin=round(d_core, 2),
                             coordinates_heme=coord, source="docked")
        print(f"{name:9}: Fe {d_fe:5.2f} A  porphyrin {d_core:5.2f} A  "
              f"coordinates_heme={coord}")

    # verdict
    gtd = {k: v for k, v in results.items() if k.startswith("GTD")}
    none_coord = all(not v["coordinates_heme"] for v in gtd.values())
    min_fe = min(v["min_dist_to_Fe"] for v in gtd.values())
    print("-" * 60)
    print(f"GTD candidates coordinating heme Fe (<{COORD_CUTOFF} A): "
          f"{sum(v['coordinates_heme'] for v in gtd.values())}/{len(gtd)}")
    print(f"Closest any-GTD approach to Fe: {min_fe:.2f} A")
    print(f"Paper claim reproduced (no GTD heme coordination): {none_coord}")

    json.dump(results, open(OUT / "cyp2c9_heme_analysis.json", "w"), indent=2)
    print(f"Wrote {OUT / 'cyp2c9_heme_analysis.json'}")


if __name__ == "__main__":
    main()
