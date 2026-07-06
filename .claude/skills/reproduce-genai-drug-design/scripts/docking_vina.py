#!/usr/bin/env python3
"""
Part 2 - Open-source molecular docking (AutoDock Vina) reproduction.

Reproduces the docking experiments in Table 3 of Chikhale et al. using Vina
instead of GOLD/ChemPLP:
  * DprE1 wild-type      (PDB 4KW5, native ligand 1W6 = TCA1)
  * DprE1 mutant Y314C   (PDB 5OEL, native ligand 1W6 = TCA1)
  * CYP2C9 (off-target)  (PDB 5W0C, native ligand 9W6 = TCA007, cofactor HEM)

Cofactors (FAD / HEM) are kept as part of the rigid receptor so ligand pi-pi
stacking with them is represented. The search box is centred on the native
co-crystallised ligand. Re-docking the native ligand gives an RMSD-to-crystal
validation (paper reports 0.937 A for TCA1 in 4KW5, 0.590 A for TCA007 in 5W0C).

Scores are Vina binding affinities (kcal/mol, more negative = better) and are
NOT on the same scale as GOLD ChemPLP; we compare RANKINGS/trends, not absolute
values.
"""
import os
import json
import subprocess
import numpy as np
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
BASE = Path(__file__).resolve().parent.parent
PREP = BASE / "data" / "prep"
PDB = BASE / "data" / "pdb"
OUT = BASE / "results"
PREP.mkdir(exist_ok=True, parents=True)

# native-ligand code and cofactor per target
TARGETS = {
    "DprE1_WT":    dict(pdb="4KW5", cofactor="FAD", native="1W6"),
    "DprE1_Y314C": dict(pdb="5OEL", cofactor="FAD", native="1W6"),
    "CYP2C9":      dict(pdb="5W0C", cofactor="HEM", native="9W6"),
}

# ligands to dock (SMILES from Appendix-2); TCA derivatives left out unless
# a structure is available - GTD series + TCA1 is the core comparison
from importlib import import_module
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
CAND = import_module("candidates").CANDIDATES


# ------------------------------------------------------------------ cleaning --
def clean_structures():
    """Write protein-only (altloc-stripped) receptor, cofactor, native ligand."""
    for t in TARGETS.values():
        pdb, cof, lig = t["pdb"], t["cofactor"], t["native"]
        prot, cofactor, native = [], [], []
        seen = set()
        for line in open(PDB / f"{pdb}.pdb"):
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[21] != "A":          # chain A only
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):  # keep a single conformer
                continue
            line = line[:16] + " " + line[17:]  # blank the altloc column
            resn = line[17:20].strip()
            if rec == "ATOM":
                prot.append(line)
            elif resn == cof:
                cofactor.append(line)
            elif resn == lig:
                native.append(line)
        (PREP / f"{pdb}_protein.pdb").write_text("".join(prot) + "END\n")
        (PREP / f"{pdb}_cofactor.pdb").write_text("".join(cofactor) + "END\n")
        (PREP / f"{pdb}_native.pdb").write_text("".join(native) + "END\n")
        print(f"{pdb}: protein={len(prot)} {cof}={len(cofactor)} native({lig})={len(native)}")


# --------------------------------------------------------------- receptor prep --
def prep_receptor(pdb):
    """protein -> PDBQT (meeko) then append cofactor PDBQT (obabel) = rigid receptor."""
    base = PREP / f"{pdb}_protein"
    subprocess.run(
        ["mk_prepare_receptor.py", "--read_pdb", str(PREP / f"{pdb}_protein.pdb"),
         "-o", str(base), "-p", "--allow_bad_res", "--default_altloc", "A",
         "--charge_model", "gasteiger", "--box_enveloping",
         str(PREP / f"{pdb}_native.pdb"), "--padding", "6"],
        capture_output=True, text=True,
    )
    prot_pdbqt = Path(str(base) + ".pdbqt")
    if not prot_pdbqt.exists():
        raise RuntimeError(f"receptor prep failed for {pdb}")

    # cofactor -> pdbqt (rigid) via obabel, then merge ATOM/HETATM lines
    cof_pdbqt = PREP / f"{pdb}_cofactor.pdbqt"
    subprocess.run(["obabel", str(PREP / f"{pdb}_cofactor.pdb"),
                    "-opdbqt", "-O", str(cof_pdbqt), "-xr", "-p", "7.4"],
                   capture_output=True, text=True)

    merged = PREP / f"{pdb}_receptor.pdbqt"
    lines = [l for l in prot_pdbqt.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]
    if cof_pdbqt.exists():
        lines += [l for l in cof_pdbqt.read_text().splitlines()
                  if l.startswith(("ATOM", "HETATM"))]
    merged.write_text("\n".join(lines) + "\n")

    # box center from native ligand centroid
    coords = []
    for line in open(PREP / f"{pdb}_native.pdb"):
        if line.startswith(("ATOM", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    c = np.mean(coords, axis=0)
    return str(merged), tuple(c)


# ---------------------------------------------------------------- ligand prep --
def prep_ligand(name, smiles, out_pdbqt):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=1, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol)
    sdf = str(out_pdbqt).replace(".pdbqt", ".sdf")
    Chem.SDWriter(sdf).write(mol)
    r = subprocess.run(["mk_prepare_ligand.py", "-i", sdf, "-o", str(out_pdbqt)],
                       capture_output=True, text=True)
    return out_pdbqt.exists()


def native_pdbqt(pdb):
    """native crystal ligand -> pdbqt for RMSD reference (obabel, keep coords)."""
    out = PREP / f"{pdb}_native.pdbqt"
    subprocess.run(["obabel", str(PREP / f"{pdb}_native.pdb"),
                    "-opdbqt", "-O", str(out), "-p", "7.4"],
                   capture_output=True, text=True)
    return out


# ---------------------------------------------------------------------- dock --
def dock(receptor_pdbqt, center, ligand_pdbqt, box=22.0, exhaustiveness=10):
    from vina import Vina
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=[box, box, box])
    v.dock(exhaustiveness=exhaustiveness, n_poses=10)
    energies = v.energies(n_poses=1)
    return float(energies[0][0])  # best affinity kcal/mol


def main():
    clean_structures()
    receptors = {}
    for name, t in TARGETS.items():
        rp, center = prep_receptor(t["pdb"])
        receptors[name] = (rp, center)
        print(f"receptor ready: {name} ({t['pdb']})  box center {np.round(center,1)}")

    results = {name: {} for name in TARGETS}
    ligands = dict(CAND)  # TCA1 + GTD_9.1..9.10

    for name, (rp, center) in receptors.items():
        print(f"\n=== docking into {name} ===")
        for lig, smi in ligands.items():
            lp = PREP / f"{TARGETS[name]['pdb']}_{lig}.pdbqt"
            if not prep_ligand(lig, smi, lp):
                print(f"  {lig}: ligand prep FAILED"); continue
            try:
                score = dock(rp, center, lp)
            except Exception as e:
                print(f"  {lig}: dock error {e}"); continue
            results[name][lig] = round(score, 2)
            print(f"  {lig:9}  Vina affinity = {score:6.2f} kcal/mol")

    json.dump(results, open(OUT / "docking_results.json", "w"), indent=2)
    print(f"\nWrote {OUT / 'docking_results.json'}")

    # quick trend summary vs paper (GTD series should beat TCA1 on DprE1)
    for name in ("DprE1_WT", "DprE1_Y314C"):
        r = results[name]
        if "TCA1" in r:
            gtd = [v for k, v in r.items() if k.startswith("GTD")]
            better = sum(1 for v in gtd if v < r["TCA1"])
            print(f"{name}: {better}/{len(gtd)} GTD candidates dock better than TCA1 "
                  f"(TCA1={r['TCA1']})")


if __name__ == "__main__":
    main()
