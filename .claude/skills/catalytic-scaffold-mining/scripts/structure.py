"""Minimal structure I/O and superposition for catalytic-scaffold mining.

Deliberately dependency-light (numpy only) so the pipeline runs anywhere:
AlphaFold DB models are plain, well-formed PDB/mmCIF and do not need a full
parser stack.

In AlphaFold models the PDB B-factor column carries **pLDDT**, not a thermal
factor. `Structure.plddt()` exposes that; never average it over the whole
chain when judging an active site -- score the catalytic shell separately.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

# Atom that carries the carboxylate carbon, used for acidic-cluster geometry.
CARBOXYLATE_ATOM = {"ASP": "CG", "GLU": "CD"}


@dataclass
class Residue:
    resnum: int           # author residue number as written in the file
    resname: str          # three-letter code
    seq_index: int        # 1-based position in the chain's ordered residue list
    atoms: dict = field(default_factory=dict)   # atom name -> xyz
    bfactors: dict = field(default_factory=dict)

    @property
    def one(self) -> str:
        return THREE_TO_ONE.get(self.resname, "X")

    @property
    def ca(self):
        return self.atoms.get("CA")

    @property
    def plddt(self):
        """AlphaFold stores pLDDT in the B-factor column; constant per residue."""
        return self.bfactors.get("CA")

    def carboxylate_c(self):
        """Side-chain carboxylate carbon for Asp/Glu, else None."""
        name = CARBOXYLATE_ATOM.get(self.resname)
        return self.atoms.get(name) if name else None


class Structure:
    """One chain of one model, in file order."""

    def __init__(self, residues, name="query", chain="A"):
        self.residues = residues
        self.name = name
        self.chain = chain
        self._by_resnum = {r.resnum: r for r in residues}
        self._by_index = {r.seq_index: r for r in residues}

    def __len__(self):
        return len(self.residues)

    @property
    def sequence(self) -> str:
        return "".join(r.one for r in self.residues)

    def by_resnum(self, n):
        return self._by_resnum.get(n)

    def by_index(self, i):
        """1-based position in the chain -- this is what Foldseek aligns on."""
        return self._by_index.get(i)

    def resnum_to_index(self, n):
        r = self._by_resnum.get(n)
        return r.seq_index if r else None

    def ca_array(self, indices):
        return np.array([self._by_index[i].ca for i in indices], dtype=float)

    def plddt(self, indices=None):
        """Mean pLDDT over the given 1-based indices (all residues if None)."""
        rs = self.residues if indices is None else [
            self._by_index[i] for i in indices if i in self._by_index
        ]
        vals = [r.plddt for r in rs if r.plddt is not None]
        return float(np.mean(vals)) if vals else None


def _open_text(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_pdb(path, chain=None, model=1) -> Structure:
    """Parse the first model of a PDB file (plain or .gz)."""
    residues, seen, idx = [], {}, 0
    current_model = 1
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                current_model = int(line.split()[1])
            elif line.startswith("ENDMDL"):
                if current_model == model:
                    break
            if current_model != model:
                continue
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname not in THREE_TO_ONE:
                continue                        # skip ligands/waters
            ch = line[21].strip() or "A"
            if chain is not None and ch != chain:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            resnum = int(line[22:26])
            icode = line[26].strip()
            key = (ch, resnum, icode)
            if key not in seen:
                idx += 1
                res = Residue(resnum=resnum, resname=resname, seq_index=idx)
                seen[key] = res
                residues.append(res)
            res = seen[key]
            atom = line[12:16].strip()
            res.atoms[atom] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
            try:
                res.bfactors[atom] = float(line[60:66])
            except ValueError:
                pass
    if chain is None and residues:
        # keep only the first chain encountered, re-indexed
        first = None
        with _open_text(path) as fh:
            for line in fh:
                if line.startswith("ATOM  "):
                    first = line[21].strip() or "A"
                    break
        if first is not None:
            return read_pdb(path, chain=first, model=model)
    return Structure(residues, name=str(path), chain=chain or "A")


def read_cif(path, chain=None) -> Structure:
    """Parse the atom_site loop of an mmCIF file (plain or .gz).

    Handles the AlphaFold DB flavour: a single `loop_` with named columns.
    """
    cols, rows, in_loop, header = {}, [], False, []
    with _open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("loop_"):
                in_loop, header = True, []
                continue
            if in_loop and s.startswith("_atom_site."):
                header.append(s.split(".", 1)[1])
                continue
            if header and (s.startswith("ATOM") or s.startswith("HETATM")):
                in_loop = False
                if not cols:
                    cols = {n: i for i, n in enumerate(header)}
                rows.append(s.split())
                continue
            if cols and (s.startswith("ATOM") or s.startswith("HETATM")):
                rows.append(s.split())
            elif cols and (s == "#" or s.startswith("loop_")):
                if rows:
                    break
    if not cols:
        raise ValueError(f"no atom_site loop found in {path}")

    def col(row, name, default=None):
        i = cols.get(name)
        return row[i] if i is not None and i < len(row) else default

    residues, seen, idx = [], {}, 0
    model_target = None
    for row in rows:
        resname = col(row, "label_comp_id", "")
        if resname not in THREE_TO_ONE:
            continue
        mdl = col(row, "pdbx_PDB_model_num", "1")
        if model_target is None:
            model_target = mdl
        if mdl != model_target:
            continue
        ch = col(row, "auth_asym_id") or col(row, "label_asym_id") or "A"
        if chain is not None and ch != chain:
            continue
        alt = col(row, "label_alt_id", ".")
        if alt not in (".", "?", "A"):
            continue
        rn = col(row, "auth_seq_id") or col(row, "label_seq_id")
        try:
            resnum = int(rn)
        except (TypeError, ValueError):
            continue
        key = (ch, resnum)
        if key not in seen:
            idx += 1
            res = Residue(resnum=resnum, resname=resname, seq_index=idx)
            seen[key] = res
            residues.append(res)
        res = seen[key]
        atom = col(row, "label_atom_id", "").strip('"')
        try:
            res.atoms[atom] = np.array(
                [float(col(row, "Cartn_x")), float(col(row, "Cartn_y")),
                 float(col(row, "Cartn_z"))]
            )
            res.bfactors[atom] = float(col(row, "B_iso_or_equiv", "0") or 0)
        except (TypeError, ValueError):
            continue
    return Structure(residues, name=str(path), chain=chain or "A")


def read_structure(path, chain=None) -> Structure:
    p = str(path)
    base = p[:-3] if p.endswith(".gz") else p
    if base.endswith((".cif", ".mmcif")):
        return read_cif(path, chain=chain)
    return read_pdb(path, chain=chain)


def superpose(P, Q):
    """Optimally superpose P onto Q. Returns (rmsd, per-pair deviations)."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or len(P) < 3:
        return float("nan"), np.array([])
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    rot = V @ np.diag([1.0, 1.0, d]) @ Wt
    dev = np.linalg.norm((Pc @ rot) - Qc, axis=1)
    return float(np.sqrt((dev ** 2).mean())), dev


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal rigid superposition. Returns (rmsd, n)."""
    rmsd, dev = superpose(P, Q)
    return rmsd, len(dev) if len(dev) else len(P)


def core_rmsd(P, Q, cutoff=4.0, min_frac=0.5, max_iter=10):
    """Outlier-trimmed RMSD -- the metric to judge an active site by.

    A flexible loop next to the catalytic residues (activation loops, lid
    segments, crystal-vs-model differences) can drag a whole-shell RMSD from
    ~1 A to >3 A while the catalytic core itself superimposes perfectly. This
    iteratively drops pairs deviating by more than `cutoff` and re-fits, which
    is what "RMSD over the conserved catalytic region" means in practice.

    Returns (rmsd, n_core, fraction_kept).
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    n0 = len(P)
    if n0 < 3:
        return float("nan"), n0, 0.0
    keep = np.arange(n0)
    rmsd, dev = superpose(P, Q)
    for _ in range(max_iter):
        new = keep[dev <= cutoff]
        if len(new) < max(3, int(min_frac * n0)) or len(new) == len(keep):
            break
        keep = new
        rmsd, dev = superpose(P[keep], Q[keep])
    return rmsd, len(keep), len(keep) / n0
