r"""Download and parse real antibody-antigen structures into training examples.

Pulls mmCIF from RCSB, reads one CA per residue per auth chain, classifies each
chain as heavy / light / antigen (ANARCI when available, else J-motif + length
heuristics), and annotates CDRs. The result feeds abproto.train.dataset.

Structures are the real 3D data the IPA model needs (the GraphQL pipeline only
had sequences). Files are cached under ``cache_dir``.
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from .anarci import annotate_chain

_THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
CIF_URL = "https://files.rcsb.org/download/{}.cif"


@dataclass
class Chain:
    chain_id: str
    seq: str
    ca: np.ndarray                  # (N, 3)
    role: str = "?"                 # "heavy" | "light" | "antigen"
    cdr_spans: dict = field(default_factory=dict)


@dataclass
class ComplexStructure:
    pdb_id: str
    heavy: Optional[Chain] = None
    light: Optional[Chain] = None
    antigen: list[Chain] = field(default_factory=list)

    def ok(self) -> bool:
        return (self.heavy is not None or self.light is not None) and len(self.antigen) > 0


def fetch_cif(pdb_id: str, cache_dir: Path, timeout: float = 90.0) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb_id.upper()}.cif"
    if not path.exists() or path.stat().st_size == 0:
        urllib.request.urlretrieve(CIF_URL.format(pdb_id.upper()), path)  # noqa: S310
    return path


def parse_ca_chains(cif_path: Path) -> dict[str, tuple[str, np.ndarray]]:
    """auth chain -> (sequence, CA coords). First model, first altloc only."""
    cols: list[str] = []
    in_loop = False
    seqs: dict[str, list[str]] = {}
    xyz: dict[str, list[list[float]]] = {}
    seen: dict[str, set] = {}
    for line in cif_path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1]); in_loop = True; continue
        if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
            f = line.split()
            if len(f) < len(cols):
                continue
            rec = dict(zip(cols, f))
            if rec.get("label_atom_id") != "CA":
                continue
            model = rec.get("pdbx_PDB_model_num", "1")
            if model not in ("1", "."):
                continue
            ch = rec.get("auth_asym_id", rec.get("label_asym_id", "?"))
            key = (ch, rec.get("auth_seq_id"), rec.get("pdbx_PDB_ins_code", "."))
            seen.setdefault(ch, set())
            if key in seen[ch]:
                continue
            seen[ch].add(key)
            try:
                p = [float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"])]
            except (KeyError, ValueError):
                continue
            seqs.setdefault(ch, []).append(_THREE2ONE.get(rec.get("label_comp_id", ""), "X"))
            xyz.setdefault(ch, []).append(p)
        elif in_loop and line.startswith("#") and cols:
            break
    return {ch: ("".join(seqs[ch]), np.asarray(xyz[ch], float)) for ch in seqs}


def _classify_chain(seq: str) -> tuple[str, dict]:
    """Return (role, cdr_spans). role in {heavy, light, antigen}."""
    n = len(seq)
    if 90 <= n <= 260:
        ann_h = annotate_chain(seq, prefer="H")
        ann_l = annotate_chain(seq, prefer="L")
        import re
        has_wg = re.search(r"WG.G", seq) is not None
        has_fg = re.search(r"FG.G", seq) is not None
        if has_wg and ann_h and "3" in ann_h.cdr_spans:
            return "heavy", ann_h.cdr_spans
        if has_fg and ann_l and "3" in ann_l.cdr_spans:
            return "light", ann_l.cdr_spans
    return "antigen", {}


def build_structure(pdb_id: str, cache_dir: Path) -> Optional[ComplexStructure]:
    try:
        cif = fetch_cif(pdb_id, cache_dir)
        chains = parse_ca_chains(cif)
    except Exception:
        return None
    cx = ComplexStructure(pdb_id=pdb_id.upper())
    for ch, (seq, ca) in chains.items():
        if len(seq) < 20:
            continue
        role, spans = _classify_chain(seq)
        c = Chain(ch, seq, ca, role, spans)
        if role == "heavy" and cx.heavy is None:
            cx.heavy = c
        elif role == "light" and cx.light is None:
            cx.light = c
        elif role == "antigen":
            cx.antigen.append(c)
    return cx if cx.ok() else None


def load_structures(pdb_ids: list[str], cache_dir: str) -> tuple[list[ComplexStructure], list[str]]:
    cd = Path(cache_dir)
    good, skipped = [], []
    for pid in pdb_ids:
        cx = build_structure(pid, cd)
        (good if cx is not None else skipped).append(cx if cx is not None else pid)
    good = [c for c in good if isinstance(c, ComplexStructure)]
    return good, skipped
