#!/usr/bin/env python3
"""Convert all ChemDraw .cdx ligand files to SMILES via the OpenBabel CLI,
running each file in an isolated subprocess so a hard crash (SIGABRT) in the
C++ CDX reader only loses that one file instead of the whole run.

Key format:  "<doi-folder>||<file-stem>"  ->  canonical SMILES
"""
import os
import glob
import json
import subprocess
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = "data/ligands_extracted/Structures_of_ligands"
OUT = "data/cdx_smiles.json"


def convert(path):
    try:
        r = subprocess.run(
            ["obabel", path, "-osmi"],
            capture_output=True, text=True, timeout=20,
        )
    except (subprocess.TimeoutExpired, Exception):
        return None
    if r.returncode != 0 or not r.stdout.strip():
        return None
    smi = r.stdout.strip().split("\t")[0].split()[0].strip()
    if not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def main():
    files = sorted(glob.glob(os.path.join(ROOT, "*", "*.cdx")))
    out = {}
    ok = fail = 0
    for i, path in enumerate(files):
        folder = os.path.basename(os.path.dirname(path))
        stem = os.path.splitext(os.path.basename(path))[0]
        smi = convert(path)
        if smi:
            out[folder + "||" + stem] = smi
            ok += 1
        else:
            fail += 1
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(files)}  ok={ok} fail={fail}", flush=True)
    json.dump(out, open(OUT, "w"))
    print(f"DONE  converted={ok}  failed={fail}  total_files={len(files)}")


if __name__ == "__main__":
    main()
