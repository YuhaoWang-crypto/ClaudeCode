#!/usr/bin/env python3
"""
Prepare MD starting complexes for the Modal pipeline (Part 3).

For each ligand we dock into the target with Vina (reusing Part 2), export the
top pose to SDF, and copy the cleaned protein PDB. This mirrors the paper, which
starts each 500 ns MD from the docked protein-ligand complex.

Output layout consumed by app.py:
    modal_md/inputs/<system>/protein.pdb
    modal_md/inputs/<system>/ligand.sdf     # top docked pose, 3D

Usage:
    python modal_md/prepare_inputs.py                 # all GTD + TCA1 into DprE1_WT
    python modal_md/prepare_inputs.py --target DprE1_Y314C --only GTD_9.7 GTD_9.4
"""
import sys
import argparse
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "src"))
from importlib import import_module
d = import_module("04_docking")

INPUTS = Path(__file__).resolve().parent / "inputs"


def export_pose(receptor_pdbqt, center, name, smiles, target_pdb, out_sdf):
    """Dock and write the best pose as SDF."""
    from vina import Vina
    lp = d.PREP / f"{target_pdb}_{name}_mdprep.pdbqt"
    if not d.prep_ligand(name, smiles, lp):
        return False
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(str(lp))
    v.compute_vina_maps(center=list(center), box_size=[22, 22, 22])
    v.dock(exhaustiveness=16, n_poses=10)
    pose_pdbqt = str(out_sdf).replace(".sdf", "_pose.pdbqt")
    v.write_poses(pose_pdbqt, n_poses=1, overwrite=True)
    # pdbqt -> sdf (adds Hs, 3D preserved)
    subprocess.run(["obabel", pose_pdbqt, "-osdf", "-O", str(out_sdf), "-h"],
                   capture_output=True, text=True)
    return out_sdf.exists()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="DprE1_WT", choices=list(d.TARGETS))
    ap.add_argument("--only", nargs="*", default=None,
                    help="subset of ligand names; default = all")
    args = ap.parse_args()

    d.clean_structures()
    t = d.TARGETS[args.target]
    receptor_pdbqt, center = d.prep_receptor(t["pdb"])

    ligands = dict(d.CAND)
    if args.only:
        ligands = {k: v for k, v in ligands.items() if k in args.only}

    for name, smi in ligands.items():
        sysname = f"{args.target}__{name}"
        outdir = INPUTS / sysname
        outdir.mkdir(parents=True, exist_ok=True)
        # protein (chain A, no cofactor -> parametrised by ff14SB in tleap;
        # see README for adding FAD/HEM cofactor parameters)
        (outdir / "protein.pdb").write_bytes(
            (d.PREP / f"{t['pdb']}_protein.pdb").read_bytes())
        ok = export_pose(receptor_pdbqt, center, name, smi, t["pdb"],
                         outdir / "ligand.sdf")
        print(f"{sysname}: {'OK' if ok else 'FAILED'}")

    print(f"\nInputs written under {INPUTS}")
    print("Next:  modal run modal_md/app.py --all --ns 500   (or --system <name>)")


if __name__ == "__main__":
    main()
