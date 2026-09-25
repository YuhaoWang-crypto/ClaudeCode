"""Stage 4 - pick the screening subset and build 3D ligands.

Docking the full 2,200-compound set into 6 receptors is ~13,200 docking runs,
which does not fit the compute budget here. The subset is therefore chosen to
preserve the two properties the benchmark depends on:

* intact analogue series - whole Bemis-Murcko scaffold groups are taken, never
  sampled across. A first attempt sampled one compound per scaffold round
  robin; with 935 scaffolds available the 800-compound quota filled on the
  first pass, so every selected compound had a unique scaffold. That makes a
  scaffold split indistinguishable from a random split and erases the
  within-series SAR that QSAR is supposed to learn - the benchmark would have
  been measuring nothing.
* a cap of MAX_PER_SCAFFOLD on any one group, so the largest series (76
  compounds) cannot eat a tenth of the budget.

Scaffolds are drawn in random order, so the mix of large series and singletons
follows the source distribution rather than being engineered.

Potency is *not* used to select compounds. Selecting on the label would bias
the benchmark before a single model is fitted.

3D embedding uses ETKDGv3 with a fixed seed and MMFF94s minimisation; the
lowest-energy of N_CONFS conformers is written out. Docking re-samples torsions
anyway, so this only has to be a sane starting geometry.
"""

from __future__ import annotations

import json
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
LIGDIR = ROOT / "data" / "ligands"

N_SELECT = 800
MAX_PER_SCAFFOLD = 25
SEED = 42
N_CONFS = 8


def select_subset(df: pd.DataFrame) -> pd.DataFrame:
    """Take whole scaffold groups in random order until the quota is met."""
    rng = np.random.default_rng(SEED)

    by_scaffold: dict[str, list[int]] = defaultdict(list)
    for idx, scaf in df["scaffold"].items():
        by_scaffold[scaf].append(idx)

    scaffolds = list(by_scaffold)
    rng.shuffle(scaffolds)

    chosen: list[int] = []
    for scaf in scaffolds:
        members = by_scaffold[scaf]
        if len(members) > MAX_PER_SCAFFOLD:
            members = list(rng.choice(members, MAX_PER_SCAFFOLD, replace=False))
        # Keep the group intact: stop before a group that would overshoot,
        # rather than splitting it across the selection boundary.
        if len(chosen) + len(members) > N_SELECT:
            if len(chosen) >= N_SELECT * 0.97:
                break
            continue
        chosen.extend(members)
        if len(chosen) >= N_SELECT:
            break
    return df.loc[chosen]


def embed(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=N_CONFS, params=params)
    if not cids:
        # Fall back to a random-coordinate embed for awkward macrocycles.
        params.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=N_CONFS, params=params)
        if not cids:
            return None
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(
            mol, mmffVariant="MMFF94s", maxIters=500
        )
        energies = [e for _conv, e in res]
        best = int(np.argmin(energies))
    except Exception:
        best = 0
    keep = Chem.Mol(mol)
    keep.RemoveAllConformers()
    keep.AddConformer(mol.GetConformer(cids[best]), assignId=True)
    return keep


def main() -> None:
    LIGDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA / "cdk2_dataset.csv")

    subset = select_subset(df).reset_index(drop=True)
    print(
        f"selected {len(subset)} of {len(df)} compounds across "
        f"{subset['scaffold'].nunique()} scaffolds "
        f"(source pool had {df['scaffold'].nunique()})"
    )
    print(
        f"  pIC50 selected: mean {subset['pchembl'].mean():.2f} "
        f"sd {subset['pchembl'].std():.2f} "
        f"range {subset['pchembl'].min():.2f}-{subset['pchembl'].max():.2f}"
    )
    print(
        f"  pIC50 full pool: mean {df['pchembl'].mean():.2f} "
        f"sd {df['pchembl'].std():.2f}   (selection is label-blind, so these "
        f"should agree)"
    )

    records = []
    failures = []
    n_cached = 0
    for i, row in subset.iterrows():
        lig_id = f"L{i:04d}"
        sdf = LIGDIR / f"{lig_id}.sdf"
        pdbqt = LIGDIR / f"{lig_id}.pdbqt"

        # Resume support: embedding is deterministic given SEED, so anything
        # already on disk is exactly what this run would produce.
        if sdf.exists() and pdbqt.exists() and pdbqt.stat().st_size > 0:
            n_cached += 1
        else:
            mol = embed(row["smiles"])
            if mol is None:
                failures.append(row["molecule_chembl_id"])
                continue
            mol.SetProp("_Name", lig_id)
            Chem.MolToMolFile(mol, str(sdf))
            proc = subprocess.run(
                ["obabel", str(sdf), "-O", str(pdbqt)],
                capture_output=True,
            )
            if proc.returncode != 0 or not pdbqt.exists() or pdbqt.stat().st_size == 0:
                failures.append(row["molecule_chembl_id"])
                continue
        records.append(
            {
                "ligand_id": lig_id,
                "molecule_chembl_id": row["molecule_chembl_id"],
                "smiles": row["smiles"],
                "pchembl": row["pchembl"],
                "scaffold": row["scaffold"],
                "n_meas": row["n_meas"],
                "sdf": sdf.name,
                "pdbqt": pdbqt.name,
            }
        )
        if (i + 1) % 100 == 0:
            print(f"  embedded {i + 1}/{len(subset)}")

    out = pd.DataFrame(records)
    out.to_csv(DATA / "screening_set.csv", index=False)
    print(
        f"\n3D-ready ligands: {len(out)}  "
        f"({n_cached} reused from a previous run, failed: {len(failures)})"
    )

    (DATA / "ligand_prep_summary.json").write_text(
        json.dumps(
            {
                "n_requested": N_SELECT,
                "max_per_scaffold": MAX_PER_SCAFFOLD,
                "n_selected": int(len(subset)),
                "n_embedded": int(len(out)),
                "n_failed": len(failures),
                "n_reused_from_disk": n_cached,
                "n_scaffolds": int(out["scaffold"].nunique()),
                "seed": SEED,
                "failed_ids": failures[:50],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
