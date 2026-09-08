"""Structure input: fetching, minimal PDB parsing, and ligand-contact ground truth.

Deliberately dependency-free (numpy only) so the analysis side of the pipeline
runs anywhere. AF2BIND itself consumes a PDB file, so that is the only format
handled here; convert mmCIF upstream if needed.
"""

from __future__ import annotations

import os
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

RCSB_URL = "https://files.rcsb.org/view/{code}.pdb"
AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.pdb"

#: Common crystallisation additives / cryoprotectants / buffer components. These
#: are HETATMs but are not the "small-molecule ligand" a binding-site predictor
#: is being judged on, so they are excluded from ground truth by default.
CRYSTALLIZATION_ADDITIVES = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "MPD",
    "ACT", "ACY", "FMT", "CIT", "TRS", "MES", "EPE", "IMD", "DMS", "IPA", "MOH",
    "NO3", "CL", "BR", "IOD", "NA", "K", "MG", "CA", "ZN", "MN", "CD", "NI",
    "CO", "CU", "FE", "HG", "AZI", "SCN", "BME", "DTT", "TLA", "MLI", "OXL",
    "PGO", "P6G", "12P", "15P", "2PE", "XPE", "SIN", "BEN", "URE", "GLC",
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


@dataclass
class Residue:
    chain: str
    resi: int
    icode: str
    resn: str
    hetero: bool
    atoms: dict = field(default_factory=dict)  # name -> (3,) xyz
    bfactors: dict = field(default_factory=dict)  # name -> float

    @property
    def key(self) -> tuple[str, int]:
        return (self.chain, self.resi)

    @property
    def one_letter(self) -> str:
        return AA3_TO_1.get(self.resn, "X")

    def coords(self, heavy_only: bool = True) -> np.ndarray:
        names = [n for n in self.atoms if not (heavy_only and n.startswith("H"))]
        if not names:
            return np.zeros((0, 3))
        return np.stack([self.atoms[n] for n in names])

    def mean_bfactor(self) -> float:
        vals = list(self.bfactors.values())
        return float(np.mean(vals)) if vals else float("nan")


def fetch_structure(target: str, out_dir: str | Path, afdb_version: int = 4) -> Path:
    """Resolve `target` to a local PDB file.

    Accepts an existing file path, a 4-character PDB ID, or anything else, which
    is treated as a UniProt accession and pulled from the AlphaFold database.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if os.path.isfile(target):
        return Path(target)
    if len(target) == 4 and target[0].isdigit():
        url = RCSB_URL.format(code=target.upper())
        dest = out_dir / f"{target.upper()}.pdb"
    else:
        url = AFDB_URL.format(acc=target.upper(), ver=afdb_version)
        dest = out_dir / f"AF-{target.upper()}-F1-model_v{afdb_version}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            body = r.read()
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"could not download {url}: {exc}") from exc
    if b"ATOM" not in body[:200000]:
        raise RuntimeError(f"{url} did not return a PDB file")
    dest.write_bytes(body)
    return dest


def parse_pdb(path: str | Path, model: int = 1) -> list[Residue]:
    """Parse ATOM/HETATM records of one model into per-residue records."""
    residues: dict[tuple, Residue] = {}
    order: list[tuple] = []
    current_model = 1
    for line in Path(path).read_text().splitlines():
        rec = line[:6]
        if rec == "MODEL ":
            current_model = int(line[10:14])
            continue
        if rec == "ENDMDL":
            if current_model >= model:
                break
            continue
        if rec not in ("ATOM  ", "HETATM"):
            continue
        if current_model != model:
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip() or "A"
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        xyz = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
        try:
            bfac = float(line[60:66])
        except ValueError:
            bfac = 0.0
        key = (chain, resi, icode, resn)
        if key not in residues:
            residues[key] = Residue(
                chain=chain, resi=resi, icode=icode, resn=resn,
                hetero=(rec == "HETATM"),
            )
            order.append(key)
        residues[key].atoms.setdefault(name, xyz)
        residues[key].bfactors.setdefault(name, bfac)
    return [residues[k] for k in order]


def protein_residues(residues: list[Residue], chain: str | None = None) -> list[Residue]:
    out = [r for r in residues if r.resn in AA3_TO_1 and "CA" in r.atoms]
    if chain is not None:
        out = [r for r in out if r.chain == chain]
    return out


def ligands(
    residues: list[Residue],
    min_heavy_atoms: int = 6,
    exclude: set[str] | None = None,
) -> list[Residue]:
    """HETATM groups that plausibly are the small-molecule ligand of interest."""
    exclude = CRYSTALLIZATION_ADDITIVES if exclude is None else exclude
    out = []
    for r in residues:
        if not r.hetero or r.resn in AA3_TO_1 or r.resn in exclude:
            continue
        if len(r.coords()) >= min_heavy_atoms:
            out.append(r)
    return out


def assign_ligands_to_chain(
    residues: list[Residue],
    ligand: list[Residue],
    chain: str,
) -> list[Residue]:
    """Keep only ligands whose nearest protein chain is `chain`.

    Crystal forms with several copies in the asymmetric unit contain one ligand
    per copy. Scoring the target chain against all of them turns a neighbouring
    copy's pocket into spurious ground truth.
    """
    prot = protein_residues(residues)
    if not prot:
        return []
    keep = []
    for lig in ligand:
        lig_xyz = lig.coords()
        best_d, best_chain = np.inf, None
        for r in prot:
            xyz = r.coords()
            if len(xyz) == 0:
                continue
            d = ((xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1).min()
            if d < best_d:
                best_d, best_chain = d, r.chain
        if best_chain == chain:
            keep.append(lig)
    return keep


def contact_residues(
    targets: list[Residue],
    ligand: list[Residue],
    cutoff: float = 5.0,
) -> set[tuple[str, int]]:
    """Protein residues with any heavy atom within `cutoff` A of any ligand atom."""
    if not ligand:
        return set()
    lig_xyz = np.concatenate([r.coords() for r in ligand])
    hits = set()
    for r in targets:
        xyz = r.coords()
        if len(xyz) == 0:
            continue
        d2 = ((xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1)
        if d2.min() <= cutoff * cutoff:
            hits.add(r.key)
    return hits


def strip_low_plddt(
    in_path: str | Path, out_path: str | Path, min_plddt: float = 70.0
) -> int:
    """Write a copy of a predicted model with low-pLDDT residues removed.

    AF2BIND is trained on ordered structures; disordered AlphaFold tails produce
    spurious high scores. Returns the number of residues kept.
    """
    residues = parse_pdb(in_path)
    keep = {
        r.key for r in protein_residues(residues) if r.mean_bfactor() >= min_plddt
    }
    lines = []
    for line in Path(in_path).read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            lines.append(line)
            continue
        chain = line[21].strip() or "A"
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        if (chain, resi) in keep:
            lines.append(line)
    Path(out_path).write_text("\n".join(lines) + "\n")
    return len(keep)


def write_bfactor_pdb(
    in_path: str | Path,
    out_path: str | Path,
    values: dict[tuple[str, int], float],
    scale: float = 100.0,
    default: float = 0.0,
) -> None:
    """Copy a PDB, replacing the B-factor column with per-residue scores.

    Load in PyMOL and `spectrum b, blue_white_red` to colour by p(bind).
    """
    lines = []
    for line in Path(in_path).read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            lines.append(line)
            continue
        chain = line[21].strip() or "A"
        try:
            resi = int(line[22:26])
        except ValueError:
            lines.append(line)
            continue
        v = values.get((chain, resi), default) * scale
        lines.append(f"{line[:60]}{v:6.2f}{line[66:]}")
    Path(out_path).write_text("\n".join(lines) + "\n")
