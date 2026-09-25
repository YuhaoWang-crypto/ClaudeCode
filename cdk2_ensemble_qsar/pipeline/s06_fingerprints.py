"""Stage 6 - per-residue interaction fingerprints from the docked poses.

For every (ligand, slice) pose this records, for each of the pocket
contact-shell residues, four geometric channels:

  mind  minimum heavy-atom distance from the residue to the ligand (capped)
  ncon  ligand heavy atoms within CONTACT_CUTOFF of the residue
  apol  apolar atom pairs (C/S/halogen on both sides) within CONTACT_CUTOFF
  polr  N/O-to-N/O pairs within POLAR_CUTOFF

Honesty note on what these are NOT: `polr` is a distance-based proxy for a
polar contact, not a hydrogen bond. The receptors carry no explicit hydrogens
and no protonation states were assigned per residue, so donor/acceptor roles
and H-bond geometry (angle, directionality) are simply not determined here.
Calling these H-bonds would overclaim. The same applies to aromatic stacking:
a ring-centroid distance channel was considered and left out rather than
reported without the ring-normal angle that makes it meaningful.

The value of the representation is that it is *positional* - it says which
part of the pocket a ligand touches, which a scalar docking score throws away.
That is the hypothesis the benchmark tests, and a distance-geometric channel
set is enough to test it.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from rdkit import Chem, RDLogger
from scipy.spatial import cKDTree

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "structures"
RESULTS = ROOT / "results"
POSES = RESULTS / "poses"

CONTACT_CUTOFF = 4.5  # Angstrom, apolar/general contact
POLAR_CUTOFF = 3.5  # Angstrom, N/O-to-N/O proximity
DIST_CAP = 10.0  # mind is capped here so "far" is one value, not noise

APOLAR = {"C", "S", "SE", "F", "CL", "BR", "I"}
POLAR = {"N", "O"}


def residue_atoms(pdb_path: Path, keep: set[int]) -> dict[int, dict]:
    """Heavy atoms of the wanted residues, grouped by residue number."""
    parser = PDBParser(QUIET=True)
    chain = next(
        next(parser.get_structure("r", str(pdb_path)).get_models()).get_chains()
    )
    out: dict[int, dict] = {}
    for res in chain:
        if res.id[0] != " " or res.id[1] not in keep:
            continue
        coords, elems = [], []
        for atom in res:
            el = atom.element.strip().upper()
            if el == "H":
                continue
            coords.append(atom.get_coord())
            elems.append(el)
        if coords:
            out[res.id[1]] = {
                "coords": np.array(coords, dtype=float),
                "elems": np.array(elems),
                "resname": res.get_resname().strip(),
            }
    return out


def pose_features(
    lig_xyz: np.ndarray, lig_elems: np.ndarray, residues: dict[int, dict],
    pocket: list[int],
) -> dict[str, float]:
    """Four channels per pocket residue plus whole-pose summaries."""
    lig_apolar = np.array([e in APOLAR for e in lig_elems])
    lig_polar = np.array([e in POLAR for e in lig_elems])

    feats: dict[str, float] = {}
    contacted_any = np.zeros(len(lig_xyz), dtype=bool)
    n_res_contacted = 0
    total_contacts = 0

    for rnum in pocket:
        info = residues.get(rnum)
        if info is None:
            # Residue unresolved in this crystal form. Encoded as "maximally
            # far / no contact" and flagged, so a missing loop never looks
            # like a measured zero-distance contact.
            feats[f"r{rnum}_mind"] = DIST_CAP
            feats[f"r{rnum}_ncon"] = 0.0
            feats[f"r{rnum}_apol"] = 0.0
            feats[f"r{rnum}_polr"] = 0.0
            feats[f"r{rnum}_missing"] = 1.0
            continue
        feats[f"r{rnum}_missing"] = 0.0

        d = np.linalg.norm(
            lig_xyz[:, None, :] - info["coords"][None, :, :], axis=2
        )  # (n_lig, n_res)
        mind = float(d.min())
        feats[f"r{rnum}_mind"] = min(mind, DIST_CAP)

        near = d <= CONTACT_CUTOFF
        lig_near = near.any(axis=1)
        ncon = int(lig_near.sum())
        feats[f"r{rnum}_ncon"] = float(ncon)
        contacted_any |= lig_near
        if ncon:
            n_res_contacted += 1
            total_contacts += ncon

        res_apolar = np.array([e in APOLAR for e in info["elems"]])
        res_polar = np.array([e in POLAR for e in info["elems"]])
        feats[f"r{rnum}_apol"] = float(
            (near & lig_apolar[:, None] & res_apolar[None, :]).sum()
        )
        feats[f"r{rnum}_polr"] = float(
            ((d <= POLAR_CUTOFF) & lig_polar[:, None] & res_polar[None, :]).sum()
        )

    feats["g_n_res_contacted"] = float(n_res_contacted)
    feats["g_total_contacts"] = float(total_contacts)
    feats["g_lig_heavy"] = float(len(lig_xyz))
    feats["g_buried_frac"] = float(contacted_any.mean()) if len(lig_xyz) else 0.0
    return feats


def main() -> None:
    manifest = json.loads((STRUCT / "ensemble_manifest.json").read_text())
    pocket = manifest["pocket_residues"]
    print(f"pocket residues ({len(pocket)}): {pocket}")

    scores = pd.read_csv(RESULTS / "docking_scores.csv", index_col=0)

    frames = []
    for sl in manifest["slices"]:
        pid = sl["pdb_id"]
        pose_file = POSES / f"{pid}.sdf.gz"
        if not pose_file.exists():
            print(f"{pid}: no poses, skipped")
            continue
        residues = residue_atoms(STRUCT / sl["receptor_pdb"], set(pocket))
        missing = sorted(set(pocket) - set(residues))
        print(
            f"{pid}: {len(residues)}/{len(pocket)} pocket residues resolved"
            + (f"  missing {missing}" if missing else "")
        )

        rows = []
        with gzip.open(pose_file, "rb") as fh:
            supplier = Chem.ForwardSDMolSupplier(fh, removeHs=True, sanitize=False)
            for mol in supplier:
                if mol is None:
                    continue
                lig_id = (
                    mol.GetProp("LIGAND_ID")
                    if mol.HasProp("LIGAND_ID")
                    else mol.GetProp("_Name")
                )
                conf = mol.GetConformer()
                xyz, elems = [], []
                for i, atom in enumerate(mol.GetAtoms()):
                    if atom.GetAtomicNum() <= 1:
                        continue
                    p = conf.GetAtomPosition(i)
                    xyz.append([p.x, p.y, p.z])
                    elems.append(atom.GetSymbol().upper())
                if not xyz:
                    continue
                feats = pose_features(
                    np.array(xyz), np.array(elems), residues, pocket
                )
                feats["ligand_id"] = lig_id
                feats["receptor"] = pid
                feats["score"] = (
                    float(scores.loc[lig_id, pid])
                    if lig_id in scores.index and pid in scores.columns
                    else np.nan
                )
                rows.append(feats)
        df = pd.DataFrame(rows)
        print(f"  {len(df)} poses featurised")
        frames.append(df)

    long = pd.concat(frames, ignore_index=True)
    long.to_csv(RESULTS / "ifp_long.csv.gz", index=False)
    print(f"\nwrote ifp_long.csv.gz  shape={long.shape}")
    print(
        f"  {long['ligand_id'].nunique()} ligands x "
        f"{long['receptor'].nunique()} slices"
    )


if __name__ == "__main__":
    main()
