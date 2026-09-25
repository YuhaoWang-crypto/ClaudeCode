"""Stage 2 - turn CDK2 crystal entries into a comparable conformational ensemble.

The whole premise of the benchmark is that each structure is a *slice* of the
same conformational landscape, so the slices have to live in one frame:

* every entry is superposed onto the reference (1HCK chain A) on the CA atoms
  of the *pocket* residues, not on all CA atoms. Superposing globally puts
  1FIN - the cyclin-A-bound active conformation - 3.9 A from the reference,
  which walks its ATP site out of the shared box even though the pocket itself
  overlays well. Aligning on the pocket keeps the thing being docked into
  registered across slices; the global RMSD is still reported, as the honest
  measure of how far apart the slices are.
* the ATP-site ligand is found geometrically - the HETATM group closest to the
  hinge residue Leu83 - rather than by a hard-coded three-letter code.
* only hinge-proximal ligands size the box. 3PXF carries two copies of ANS but
  both sit in the allosteric pocket under the C-helix and its ATP site is
  empty, so 3PXF contributes a receptor slice but no box geometry. That slice
  is deliberately kept: it is the negative control for "does an inappropriate
  conformer degrade the ensemble?".
* one shared box is used for every slice. Per-slice boxes would make docking
  scores incomparable, which is exactly the confound the benchmark tests for.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "structures"

REFERENCE = "1HCK"
SLICES = ["1HCK", "1AQ1", "2VTA", "1FIN", "1KE5", "3PXF"]
CDK2_CHAIN = "A"

HINGE_RESIDUE = 83  # Leu83 - the CDK2 hinge NH/C=O that ATP-site ligands read
# The alignment frame is the ligand's first contact shell, NOT a fat sphere
# around the pocket centroid. A 15 A sphere sweeps in the activation/T-loop
# (res ~145-165), which cyclin A remodels by 8-15 A in 1FIN - real biology, but
# it drags the fit and pushes the ATP site out of register. Contact-shell
# residues are the ones the docking box actually has to line up.
POCKET_CONTACT_CUTOFF = 6.0  # heavy-atom distance to the reference ligand
ATP_SITE_MAX_HINGE_DIST = 8.0  # beyond this a ligand is not in the ATP site
BOX_PADDING = 6.0  # Angstrom beyond the union of native ligand atoms
MIN_BOX = 20.0

# Cryo/buffer junk that is never the ligand of interest.
SOLVENT = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "MPD", "PEG", "ACT", "CL",
    "NA", "K", "MG", "MN", "ZN", "CA", "IOD", "BR", "NO3", "TRS", "DMS",
}

parser = PDBParser(QUIET=True)


class ProteinOnly(Select):
    """Chain-A standard residues, highest-occupancy altloc, no hydrogens."""

    def accept_chain(self, chain):
        return chain.id == CDK2_CHAIN

    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        if atom.element == "H":
            return False
        return atom.get_altloc() in (" ", "A")


class OneLigand(Select):
    """Exactly one hetero residue, heavy atoms, highest-occupancy altloc."""

    def __init__(self, residue):
        self.target = residue

    def accept_chain(self, chain):
        return chain.id == CDK2_CHAIN

    def accept_residue(self, residue):
        return residue is self.target

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def load_chain(pdb_id: str):
    structure = parser.get_structure(pdb_id, STRUCT / f"{pdb_id}.pdb")
    model = next(structure.get_models())
    return structure, model[CDK2_CHAIN]


def ca_map(chain) -> dict[int, object]:
    return {
        res.id[1]: res["CA"]
        for res in chain
        if res.id[0] == " " and "CA" in res
    }


def find_site_ligand(chain):
    """Return the hetero residue whose atoms come closest to the Leu83 CA."""
    cas = ca_map(chain)
    if HINGE_RESIDUE not in cas:
        raise RuntimeError(f"hinge residue {HINGE_RESIDUE} missing")
    hinge = cas[HINGE_RESIDUE].get_coord()

    best, best_d = None, np.inf
    for res in chain:
        het = res.id[0]
        if het == " " or res.get_resname().strip() in SOLVENT:
            continue
        coords = np.array([a.get_coord() for a in res if a.element != "H"])
        if len(coords) < 6:  # ions and tiny fragments are not the site ligand
            continue
        d = np.linalg.norm(coords - hinge, axis=1).min()
        if d < best_d:
            best, best_d = res, d
    if best is None:
        raise RuntimeError("no ATP-site ligand found")
    return best, best_d


def rmsd_on(ref_cas: dict, cas: dict, residues: list[int]) -> float:
    """CA RMSD over the given residue numbers, with no further fitting."""
    a = np.array([ref_cas[i].get_coord() for i in residues])
    b = np.array([cas[i].get_coord() for i in residues])
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def reference_pocket_residues() -> tuple[list[int], str]:
    """Residue numbers of the reference ligand's first contact shell.

    Shared with the fingerprint stage so both use exactly one definition of
    "the pocket".
    """
    _, ref_chain = load_chain(REFERENCE)
    ref_ligand, _ = find_site_ligand(ref_chain)
    ref_lig_atoms = np.array(
        [a.get_coord() for a in ref_ligand if a.element != "H"]
    )
    residues = []
    for res in ref_chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        heavy = np.array([a.get_coord() for a in res if a.element != "H"])
        d = np.linalg.norm(
            heavy[:, None, :] - ref_lig_atoms[None, :, :], axis=2
        ).min()
        if d <= POCKET_CONTACT_CUTOFF:
            residues.append(res.id[1])
    return residues, ref_ligand.get_resname().strip()


def main() -> None:
    _, ref_chain = load_chain(REFERENCE)
    ref_cas = ca_map(ref_chain)

    pocket_residues, ref_lig_name = reference_pocket_residues()
    ref_ligand, _ = find_site_ligand(ref_chain)
    print(
        f"reference {REFERENCE}: {len(pocket_residues)} contact-shell residues "
        f"within {POCKET_CONTACT_CUTOFF:.0f} A of {ref_lig_name}"
    )
    print(f"  {pocket_residues}\n")

    manifest: list[dict] = []
    ligand_clouds: list[np.ndarray] = []

    for pdb_id in SLICES:
        structure, chain = load_chain(pdb_id)
        cas = ca_map(chain)

        # --- superpose onto the reference on shared POCKET CA atoms ---------
        shared_pocket = sorted(set(pocket_residues) & set(cas))
        shared_all = sorted(set(ref_cas) & set(cas))
        if pdb_id == REFERENCE:
            pocket_rmsd, global_rmsd = 0.0, 0.0
        else:
            sup = Superimposer()
            sup.set_atoms(
                [ref_cas[i] for i in shared_pocket], [cas[i] for i in shared_pocket]
            )
            # Apply to every atom of the entry so the ligand travels with it.
            sup.apply(list(structure.get_atoms()))
            pocket_rmsd = float(sup.rms)
            # Global RMSD is measured *in the pocket-fitted frame* - it is the
            # honest statement of how differently the rest of the kinase sits
            # once the pockets are registered.
            global_rmsd = rmsd_on(ref_cas, ca_map(chain), shared_all)

        ligand, hinge_dist = find_site_ligand(chain)
        lig_coords = np.array(
            [a.get_coord() for a in ligand if a.element != "H"]
        )
        in_atp_site = hinge_dist <= ATP_SITE_MAX_HINGE_DIST
        if in_atp_site:
            ligand_clouds.append(lig_coords)

        # --- write the aligned protein-only receptor ------------------------
        io = PDBIO()
        io.set_structure(structure)
        rec_pdb = STRUCT / f"{pdb_id}_receptor.pdb"
        io.save(str(rec_pdb), ProteinOnly())

        # --- write the aligned native ligand (for the cross-docking control) -
        # Written through PDBIO rather than by hand: hand-formatted HETATM
        # records silently break RDKit's conformer parsing for ligands with
        # four-character atom names (STU, LS1).
        lig_pdb = STRUCT / f"{pdb_id}_native_ligand.pdb"
        io.save(str(lig_pdb), OneLigand(ligand))

        n_res = sum(1 for r in chain if r.id[0] == " ")
        manifest.append(
            {
                "pdb_id": pdb_id,
                "native_ligand": ligand.get_resname().strip(),
                "hinge_min_dist_A": round(float(hinge_dist), 2),
                "ligand_in_atp_site": bool(in_atp_site),
                "pocket_rmsd_to_ref_A": round(pocket_rmsd, 3),
                "global_rmsd_in_pocket_frame_A": round(global_rmsd, 3),
                "n_pocket_ca": len(shared_pocket),
                "n_residues": n_res,
                "n_ligand_heavy_atoms": int(len(lig_coords)),
                "receptor_pdb": rec_pdb.name,
                "native_ligand_pdb": lig_pdb.name,
            }
        )
        flag = "ATP-site" if in_atp_site else "NOT-ATP-site (box excluded)"
        print(
            f"{pdb_id}: ligand={ligand.get_resname():>4s} hinge_d={hinge_dist:5.2f}A "
            f"pocket_rmsd={pocket_rmsd:5.3f}A global_rmsd={global_rmsd:5.3f}A "
            f"res={n_res}  [{flag}]"
        )

    # --- one shared docking box over the union of native ligands ------------
    allcoords = np.vstack(ligand_clouds)
    center = (allcoords.min(0) + allcoords.max(0)) / 2.0
    extent = allcoords.max(0) - allcoords.min(0) + 2 * BOX_PADDING
    size = np.maximum(extent, MIN_BOX)

    box = {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": round(float(size[0]), 3),
        "size_y": round(float(size[1]), 3),
        "size_z": round(float(size[2]), 3),
    }
    print(f"\nshared box: center={center.round(2)} size={size.round(2)}")

    # --- convert receptors to pdbqt (protonated at pH 7.4) ------------------
    for entry in manifest:
        src = STRUCT / entry["receptor_pdb"]
        dst = src.with_suffix(".pdbqt")
        subprocess.run(
            ["obabel", str(src), "-O", str(dst), "-xr", "-p", "7.4"],
            check=True,
            capture_output=True,
        )
        entry["receptor_pdbqt"] = dst.name
        print(f"  pdbqt {dst.name} ({dst.stat().st_size // 1024} KB)")

    (STRUCT / "ensemble_manifest.json").write_text(
        json.dumps(
            {
                "reference": REFERENCE,
                "box": box,
                "pocket_residues": pocket_residues,
                "pocket_contact_cutoff_A": POCKET_CONTACT_CUTOFF,
                "slices": manifest,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
