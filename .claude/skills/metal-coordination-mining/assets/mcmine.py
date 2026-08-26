#!/usr/bin/env python3
"""
mcmine — mechanism-guided metal-coordination mining of predicted structures.

Generalization of the pipeline in Kipouros & Chang, Nature 656, 763 (2026),
"Targeted enzyme discovery using metal-coordination mining"
(doi:10.1038/s41586-026-10716-z; reference code: github.com/yannikipouros/hal-discovery).

The paper's insight in one line: a mechanistically *required* 3D active-site
motif discriminates enzyme function far better than sequence similarity, and a
motif search costs O(N) instead of the O(N^2) of pairwise alignment.

This tool takes that logic and makes it motif-agnostic: the anchor pair, the
probe point, the ligand shell and the classification rules all come from a JSON
motif spec, so a new mechanism becomes a new spec rather than new code.

Subcommands
    mine       classify sites in a folder of structures against a motif spec
    type       enumerate coordination environments (anchor-Xn typing, no rules)
    benchmark  recovery/specificity against known positives + threshold sweeps

Every run writes a JSON summary that records the resolved parameters, so a
result can always be traced back to the thresholds that produced it.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser

# --------------------------------------------------------------------------
# Atom sets
#
# "coordinating": atoms that can actually ligate a metal. Use for coordination
#                 typing and for deciding what IS in the first shell.
# "proximity":    coordinating atoms plus side-chain carbons. Use to decide what
#                 is ABSENT — a carboxylate whose CB/CG is near the probe can
#                 rotate in, so excluding on OD/OE alone under-rejects. This is
#                 the set the published halogenase filter used.
# --------------------------------------------------------------------------

COORDINATING_ATOMS: dict[str, list[str]] = {
    "ASN": ["ND2", "OD1"],
    "GLN": ["NE2", "OE1"],
    "CYS": ["SG"],
    "HIS": ["NE2", "ND1"],
    "MET": ["SD"],
    "TYR": ["OH"],
    "TRP": ["NE1"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "LYS": ["NZ"],
    "ARG": ["NE", "NH1", "NH2"],
}

PROXIMITY_ATOMS: dict[str, list[str]] = {
    "ASN": ["ND2", "OD1"],
    "GLN": ["NE2", "OE1", "OE2"],
    "CYS": ["SG"],
    "HIS": ["NE2", "ND1"],
    "MET": ["SD"],
    "TYR": ["OH"],
    "TRP": ["NE1"],
    "ASP": ["OD1", "OD2", "CB", "CG"],
    "GLU": ["OE1", "OE2", "CB", "CG", "CD"],
    "ALA": ["CB"],
    "GLY": ["CA"],
    "PHE": ["CD1", "CD2", "CG", "CE1", "CE2", "CZ"],
    "ILE": ["CG1", "CG2", "CD1"],
    "LEU": ["CG", "CD1", "CD2"],
    "VAL": ["CG1", "CG2"],
    "LYS": ["CG", "CD", "CE", "NZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["OG"],
    "THR": ["OG1", "CG2"],
    "ARG": ["CZ", "NE", "NH1", "NH2"],
}

ATOM_SETS = {"coordinating": COORDINATING_ATOMS, "proximity": PROXIMITY_ATOMS}


# --------------------------------------------------------------------------
# Motif spec
# --------------------------------------------------------------------------


class MotifSpec:
    """A mechanism encoded as geometry + rules. See reference/motif-spec.md."""

    def __init__(self, raw: dict[str, Any]):
        self.raw = raw
        self.name: str = raw.get("name", "unnamed_motif")
        self.description: str = raw.get("description", "")

        anchor = raw["anchor"]
        self.anchor_residues: list[str] = [r.upper() for r in anchor["residues"]]
        if len(self.anchor_residues) != 2:
            raise ValueError("anchor.residues must list exactly two residue names")
        self.anchor_atoms: dict[str, list[str]] = {
            k.upper(): list(v) for k, v in anchor["anchor_atoms"].items()
        }
        self.anchor_distance_max: float = float(anchor["anchor_distance_max"])
        bb = anchor.get("backbone_hbond_max", None)
        self.backbone_hbond_max: float | None = None if bb is None else float(bb)
        self.same_chain_only: bool = bool(anchor.get("same_chain_only", False))
        self.min_sequence_separation: int = int(anchor.get("min_sequence_separation", 0))

        self.atom_set_name: str = raw.get("atom_set", "proximity")
        if isinstance(self.atom_set_name, dict):  # inline custom map
            self.atom_set: dict[str, list[str]] = {
                k.upper(): list(v) for k, v in self.atom_set_name.items()
            }
            self.atom_set_name = "custom"
        else:
            if self.atom_set_name not in ATOM_SETS:
                raise ValueError(f"atom_set must be one of {list(ATOM_SETS)} or a map")
            self.atom_set = ATOM_SETS[self.atom_set_name]

        self.rules: dict[str, Any] = raw.get("rules", {})
        self.contrast: dict[str, Any] | None = raw.get("contrast_rules")
        self.ligand_shell_cutoff: float = float(raw.get("ligand_shell_cutoff", 4.5))

    @classmethod
    def load(cls, path: str | Path) -> "MotifSpec":
        return cls(json.loads(Path(path).read_text()))

    def resolved(self) -> dict[str, Any]:
        """Everything a reader needs to reproduce the run."""
        return {
            "name": self.name,
            "description": self.description,
            "anchor": {
                "residues": self.anchor_residues,
                "anchor_atoms": self.anchor_atoms,
                "anchor_distance_max": self.anchor_distance_max,
                "backbone_hbond_max": self.backbone_hbond_max,
                "same_chain_only": self.same_chain_only,
                "min_sequence_separation": self.min_sequence_separation,
            },
            "atom_set": self.atom_set_name,
            "ligand_shell_cutoff": self.ligand_shell_cutoff,
            "rules": self.rules,
            "contrast_rules": self.contrast,
        }


# --------------------------------------------------------------------------
# Structure parsing
# --------------------------------------------------------------------------

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)


def _open_maybe_gz(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_structure(path: Path):
    """Parse .pdb/.cif, optionally .gz. Only the first model is used."""
    stem = path.name.lower()
    parser = _CIF_PARSER if (".cif" in stem) else _PDB_PARSER
    with _open_maybe_gz(path) as handle:
        structure = parser.get_structure(path.stem, handle)
    return next(structure.get_models())


class ResidueTable:
    """Flat, vectorized view of one model: residues + their scorable atoms."""

    def __init__(self, model, atom_set: dict[str, list[str]]):
        self.keys: list[tuple[str, int]] = []
        self.resnames: list[str] = []
        self.res_index: dict[tuple[str, int], int] = {}
        # per-residue backbone/anchor atoms
        self.bb_n: dict[int, np.ndarray] = {}
        self.bb_o: dict[int, np.ndarray] = {}
        self.anchor_coord: dict[int, dict[str, np.ndarray]] = defaultdict(dict)

        coords: list[np.ndarray] = []
        atom_res: list[int] = []
        atom_type: list[str] = []

        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # skip hetero/water
                    continue
                key = (chain.id, int(residue.id[1]))
                idx = len(self.keys)
                self.keys.append(key)
                self.res_index[key] = idx
                resname = residue.get_resname().upper()
                self.resnames.append(resname)

                if "N" in residue:
                    self.bb_n[idx] = residue["N"].get_coord()
                if "O" in residue:
                    self.bb_o[idx] = residue["O"].get_coord()
                for atom in residue:
                    name = atom.get_name()
                    self.anchor_coord[idx][name] = atom.get_coord()
                    if name in atom_set.get(resname, ()):
                        coords.append(atom.get_coord())
                        atom_res.append(idx)
                        atom_type.append(resname)

        self.n_residues = len(self.keys)
        self.atom_coords = (
            np.asarray(coords, dtype=float) if coords else np.zeros((0, 3), dtype=float)
        )
        self.atom_res = np.asarray(atom_res, dtype=int)
        self.atom_restype = np.asarray(atom_type, dtype=object)
        self.present_types = sorted(set(atom_type))

    def residues_of(self, resname: str) -> list[int]:
        return [i for i, r in enumerate(self.resnames) if r == resname]

    def resname_at_offset(self, res_idx: int, offset: int) -> str:
        chain, resseq = self.keys[res_idx]
        target = self.res_index.get((chain, resseq + offset))
        return self.resnames[target] if target is not None else "NA"


# --------------------------------------------------------------------------
# Anchor detection
# --------------------------------------------------------------------------


def _min_anchor_pair_distance(
    table: ResidueTable, i: int, j: int, atoms_i: Iterable[str], atoms_j: Iterable[str]
) -> tuple[float, str, str, np.ndarray, np.ndarray]:
    best = (np.inf, "", "", None, None)
    for ai in atoms_i:
        ci = table.anchor_coord[i].get(ai)
        if ci is None:
            continue
        for aj in atoms_j:
            cj = table.anchor_coord[j].get(aj)
            if cj is None:
                continue
            d = float(np.linalg.norm(ci - cj))
            if d < best[0]:
                best = (d, ai, aj, ci, cj)
    return best


def find_anchor_pairs(table: ResidueTable, spec: MotifSpec) -> list[dict[str, Any]]:
    """Pairs of anchor residues satisfying the side-chain and backbone geometry.

    The backbone O(i)-N(j) / N(i)-O(j) pair of constraints is what pins the two
    residues onto adjacent strands; the side-chain distance is what says they
    are preorganized to chelate. Both matter: side-chain distance alone admits
    accidental proximity, backbone alone admits any strand pair.
    """
    ra, rb = spec.anchor_residues
    idx_a = table.residues_of(ra)
    idx_b = table.residues_of(rb) if rb != ra else idx_a
    atoms_a = spec.anchor_atoms.get(ra, [])
    atoms_b = spec.anchor_atoms.get(rb, [])

    pairs: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for i in idx_a:
        for j in idx_b:
            if i == j:
                continue
            key = (min(i, j), max(i, j)) if ra == rb else (i, j)
            if key in seen:
                continue
            chain_i, seq_i = table.keys[i]
            chain_j, seq_j = table.keys[j]
            if spec.same_chain_only and chain_i != chain_j:
                continue
            if (
                spec.min_sequence_separation
                and chain_i == chain_j
                and abs(seq_i - seq_j) < spec.min_sequence_separation
            ):
                continue

            d_anchor, atom_i, atom_j, ci, cj = _min_anchor_pair_distance(
                table, i, j, atoms_a, atoms_b
            )
            if not np.isfinite(d_anchor) or d_anchor >= spec.anchor_distance_max:
                continue

            d_on = d_no = None
            if spec.backbone_hbond_max is not None:
                o_i, n_j = table.bb_o.get(i), table.bb_n.get(j)
                n_i, o_j = table.bb_n.get(i), table.bb_o.get(j)
                if any(x is None for x in (o_i, n_j, n_i, o_j)):
                    continue
                d_on = float(np.linalg.norm(o_i - n_j))
                d_no = float(np.linalg.norm(n_i - o_j))
                if d_on >= spec.backbone_hbond_max or d_no >= spec.backbone_hbond_max:
                    continue

            seen.add(key)
            pairs.append(
                {
                    "i": i,
                    "j": j,
                    "probe": (ci + cj) / 2.0,
                    "anchor_distance": d_anchor,
                    "anchor_atom_i": atom_i,
                    "anchor_atom_j": atom_j,
                    "bb_O_i_N_j": d_on,
                    "bb_N_i_O_j": d_no,
                }
            )
    return pairs


def shell_distances(
    table: ResidueTable, probe: np.ndarray, exclude: set[int]
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    """Min distance from probe to each residue type, plus a per-residue list."""
    if table.atom_coords.shape[0] == 0:
        return {}, []
    d = np.linalg.norm(table.atom_coords - probe, axis=1)
    keep = ~np.isin(table.atom_res, list(exclude)) if exclude else np.ones_like(d, bool)

    mins: dict[str, float] = {}
    per_res: dict[int, tuple[float, str]] = {}
    for k in np.flatnonzero(keep):
        rt = table.atom_restype[k]
        dk = float(d[k])
        if dk < mins.get(rt, np.inf):
            mins[rt] = dk
        r = int(table.atom_res[k])
        if dk < per_res.get(r, (np.inf, ""))[0]:
            per_res[r] = (dk, rt)

    residues = [
        {
            "ligand_chain": table.keys[r][0],
            "ligand_resseq": table.keys[r][1],
            "ligand_resname": rt,
            "ligand_min_dist": dist,
        }
        for r, (dist, rt) in sorted(per_res.items(), key=lambda kv: kv[1][0])
    ]
    return mins, residues


# --------------------------------------------------------------------------
# Site extraction
# --------------------------------------------------------------------------


def sites_in_structure(path: Path, spec: MotifSpec) -> tuple[list[dict], list[dict]]:
    model = parse_structure(path)
    table = ResidueTable(model, spec.atom_set)
    pairs = find_anchor_pairs(table, spec)
    if not pairs:
        return [], []

    ctx = spec.rules.get("sequence_context")
    site_rows: list[dict] = []
    shell_rows: list[dict] = []

    for n, pair in enumerate(pairs):
        i, j = pair["i"], pair["j"]
        mins, residues = shell_distances(table, pair["probe"], {i, j})
        site_id = f"{path.stem}:{n}"
        row: dict[str, Any] = {
            "site_id": site_id,
            "structure": path.stem,
            "chainA": table.keys[i][0],
            "resA": table.keys[i][1],
            "resA_name": table.resnames[i],
            "chainB": table.keys[j][0],
            "resB": table.keys[j][1],
            "resB_name": table.resnames[j],
            "anchor_distance": round(pair["anchor_distance"], 3),
            "bb_O_i_N_j": None if pair["bb_O_i_N_j"] is None else round(pair["bb_O_i_N_j"], 3),
            "bb_N_i_O_j": None if pair["bb_N_i_O_j"] is None else round(pair["bb_N_i_O_j"], 3),
            "protein_length": table.n_residues,
        }
        if ctx:
            anchor_idx = i if ctx.get("from", "anchorA") == "anchorA" else j
            row["context_res"] = table.resname_at_offset(anchor_idx, int(ctx["offset"]))
        for rt, dist in mins.items():
            row[f"closest_{rt}"] = round(dist, 3)
        # first-shell summary, independent of any rule
        shell = sorted(
            r["ligand_resname"] for r in residues if r["ligand_min_dist"] <= spec.ligand_shell_cutoff
        )
        row["shell_residues"] = "+".join(shell) if shell else "none"
        row["coordination_type"] = coordination_label(spec, shell)
        site_rows.append(row)

        for r in residues:
            if r["ligand_min_dist"] <= spec.ligand_shell_cutoff:
                shell_rows.append({"site_id": site_id, "structure": path.stem, **r})

    return site_rows, shell_rows


def coordination_label(spec: MotifSpec, shell: list[str]) -> str:
    """e.g. 2His-1Asp, 2His-2His, 2His-X0 — the anchor plus its extra ligands."""
    a, b = spec.anchor_residues
    anchor = f"2{a.capitalize()}" if a == b else f"1{a.capitalize()}-1{b.capitalize()}"
    if not shell:
        return f"{anchor}-X0"
    counts: dict[str, int] = defaultdict(int)
    for r in shell:
        counts[r] += 1
    extra = "-".join(f"{counts[r]}{r.capitalize()}" for r in sorted(counts))
    return f"{anchor}-{extra}"


# --------------------------------------------------------------------------
# Rule evaluation
# --------------------------------------------------------------------------


def _col(df: pd.DataFrame, resname: str) -> pd.Series:
    """Distance column for a residue type; absent type => +inf (nothing nearby)."""
    col = f"closest_{resname}"
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(np.inf)
    return pd.Series(np.inf, index=df.index)


def apply_rules(df: pd.DataFrame, rules: dict[str, Any]) -> pd.Series:
    """Evaluate a rule block; returns a boolean mask over sites.

    Rule blocks are deliberately small:
      require_absent  — every listed residue type must be FARTHER than min_distance
      require_present — at least one listed type must be WITHIN max_distance
      sequence_context— residue at a fixed sequence offset must be in `allowed`
      min_protein_length
    """
    mask = pd.Series(True, index=df.index)

    for block in rules.get("require_absent", []):
        for res in block["residues"]:
            mask &= _col(df, res.upper()) > float(block["min_distance"])

    for block in rules.get("require_present", []):
        any_present = pd.Series(False, index=df.index)
        for res in block["residues"]:
            any_present |= _col(df, res.upper()) <= float(block["max_distance"])
        mask &= any_present

    ctx = rules.get("sequence_context")
    if ctx and ctx.get("allowed") and "context_res" in df.columns:
        mask &= df["context_res"].isin([r.upper() for r in ctx["allowed"]])

    min_len = rules.get("min_protein_length")
    if min_len is not None:
        mask &= pd.to_numeric(df["protein_length"], errors="coerce").fillna(0) > float(min_len)

    return mask


# --------------------------------------------------------------------------
# Scanning
# --------------------------------------------------------------------------

STRUCTURE_GLOBS = ("*.pdb", "*.pdb.gz", "*.cif", "*.cif.gz")


def list_structures(struct_dir: Path, limit: int | None = None) -> list[Path]:
    files: list[Path] = []
    for pattern in STRUCTURE_GLOBS:
        files.extend(struct_dir.glob(pattern))
    files = sorted(set(files))
    return files[:limit] if limit else files


def scan(
    struct_dir: Path, spec: MotifSpec, limit: int | None, progress_every: int
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], int]:
    files = list_structures(struct_dir, limit)
    if not files:
        raise SystemExit(f"No structures ({', '.join(STRUCTURE_GLOBS)}) in {struct_dir}")

    print(f"Scanning {len(files)} structures in {struct_dir}", flush=True)
    site_rows: list[dict] = []
    shell_rows: list[dict] = []
    no_sites: list[str] = []
    failed: list[str] = []

    for n, path in enumerate(files, 1):
        if n == 1 or n % progress_every == 0 or n == len(files):
            print(f"[{n}/{len(files)}]", flush=True)
        try:
            rows, shells = sites_in_structure(path, spec)
        except Exception as exc:  # malformed file: record, never abort the sweep
            failed.append(f"{path.name}\t{type(exc).__name__}: {exc}")
            continue
        if not rows:
            no_sites.append(path.stem)
            continue
        site_rows.extend(rows)
        shell_rows.extend(shells)

    return pd.DataFrame(site_rows), pd.DataFrame(shell_rows), no_sites, failed, len(files)


def _stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _write_summary(outdir: Path, ts: str, payload: dict[str, Any]) -> None:
    (outdir / f"summary_{ts}.json").write_text(json.dumps(payload, indent=2, default=str))
    print("\n" + json.dumps({k: v for k, v in payload.items() if not k.endswith("_csv")}, indent=2))


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_mine(args: argparse.Namespace) -> None:
    spec = MotifSpec.load(args.motif)
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    ts = _stamp()

    df, shells, no_sites, failed, n_files = scan(
        Path(args.struct_dir).expanduser().resolve(), spec, args.max_structures, args.progress_every
    )

    sites_csv = outdir / f"sites_all_{ts}.csv"
    df.to_csv(sites_csv, index=False)
    (outdir / f"no_sites_{ts}.txt").write_text("\n".join(no_sites))
    if failed:
        (outdir / f"parse_failures_{ts}.txt").write_text("\n".join(failed))

    summary: dict[str, Any] = {
        "motif": spec.name,
        "structures_scanned": n_files,
        "structures_parse_failed": len(failed),
        "structures_without_site": len(no_sites),
        "total_sites": int(len(df)),
        "proteins_with_site": int(df["structure"].nunique()) if len(df) else 0,
        "sites_csv": str(sites_csv),
        "parameters": spec.resolved(),
    }

    if len(df):
        hits = df[apply_rules(df, spec.rules)].copy()
        hits_csv = outdir / f"hits_{spec.name}_{ts}.csv"
        hits.to_csv(hits_csv, index=False)
        ids_txt = outdir / f"hit_accessions_{ts}.txt"
        ids_txt.write_text("\n".join(sorted(hits["structure"].drop_duplicates())))
        summary.update(
            {
                "hits": int(len(hits)),
                "hit_proteins": int(hits["structure"].nunique()),
                "hit_rate_vs_sites": round(len(hits) / len(df), 6),
                "hits_csv": str(hits_csv),
                "hit_accessions": str(ids_txt),
            }
        )
        if spec.contrast:
            contrast = df[apply_rules(df, spec.contrast)]
            cname = spec.contrast.get("name", "contrast")
            contrast_csv = outdir / f"contrast_{cname}_{ts}.csv"
            contrast.to_csv(contrast_csv, index=False)
            summary["contrast_class"] = cname
            summary["contrast_sites"] = int(len(contrast))
            summary["contrast_csv"] = str(contrast_csv)
        if args.emit_shell and len(shells):
            shell_csv = outdir / f"shell_residues_{ts}.csv"
            shells.to_csv(shell_csv, index=False)
            summary["shell_csv"] = str(shell_csv)

    _write_summary(outdir, ts, summary)


def cmd_type(args: argparse.Namespace) -> None:
    """Coordination-environment census: no functional rules, just what's there."""
    spec = MotifSpec.load(args.motif)
    if args.ligand_cutoff is not None:
        spec.ligand_shell_cutoff = args.ligand_cutoff
    spec.atom_set_name = "coordinating"  # typing asks what can ligate, not what is near
    spec.atom_set = COORDINATING_ATOMS

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    ts = _stamp()

    df, shells, no_sites, failed, n_files = scan(
        Path(args.struct_dir).expanduser().resolve(), spec, args.max_structures, args.progress_every
    )
    sites_csv = outdir / f"typed_sites_{ts}.csv"
    df.to_csv(sites_csv, index=False)
    if len(shells):
        shells.to_csv(outdir / f"shell_residues_{ts}.csv", index=False)

    counts = (
        df["coordination_type"].value_counts().rename_axis("coordination_type").reset_index(name="n")
        if len(df)
        else pd.DataFrame(columns=["coordination_type", "n"])
    )
    counts_csv = outdir / f"coordination_type_counts_{ts}.csv"
    counts.to_csv(counts_csv, index=False)

    _write_summary(
        outdir,
        ts,
        {
            "motif": spec.name,
            "mode": "type",
            "structures_scanned": n_files,
            "structures_parse_failed": len(failed),
            "total_sites": int(len(df)),
            "distinct_coordination_types": int(len(counts)),
            "top_types": counts.head(15).to_dict("records"),
            "sites_csv": str(sites_csv),
            "counts_csv": str(counts_csv),
            "parameters": spec.resolved(),
        },
    )


def _read_id_list(value: str) -> set[str]:
    p = Path(value)
    if p.exists():
        return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}
    return {v.strip() for v in value.split(",") if v.strip()}


def _sweep_rules(rules: dict[str, Any], path: str, value: float) -> dict[str, Any]:
    """Override one threshold, addressed like 'require_absent.0.min_distance'."""
    out = json.loads(json.dumps(rules))
    node: Any = out
    parts = path.split(".")
    for part in parts[:-1]:
        node = node[int(part)] if part.isdigit() else node[part]
    node[parts[-1]] = value
    return out


def cmd_benchmark(args: argparse.Namespace) -> None:
    """The step that makes a mining claim believable.

    Never run a motif at database scale before it has recovered known members of
    the target class and rejected known members of the sibling class on a set
    where the answer is already known.
    """
    spec = MotifSpec.load(args.motif)
    positives = _read_id_list(args.positives)
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    ts = _stamp()

    if args.sites_csv:
        df = pd.read_csv(args.sites_csv)
        n_files = df["structure"].nunique()
    else:
        df, _, no_sites, failed, n_files = scan(
            Path(args.struct_dir).expanduser().resolve(),
            spec,
            args.max_structures,
            args.progress_every,
        )
        df.to_csv(outdir / f"sites_all_{ts}.csv", index=False)

    if not len(df):
        raise SystemExit("No sites found — the anchor geometry rejected everything.")

    scanned = set(df["structure"].unique())
    negatives_all = scanned - positives
    missing = positives - scanned
    if missing:
        print(f"WARNING: {len(missing)} positives have no anchor site at all: {sorted(missing)}")

    def evaluate(rules: dict[str, Any]) -> dict[str, Any]:
        hit_ids = set(df[apply_rules(df, rules)]["structure"].unique())
        tp = len(hit_ids & positives)
        fp = len(hit_ids & negatives_all)
        fn = len(positives - hit_ids)
        tn = len(negatives_all - hit_ids)
        return {
            "hits": len(hit_ids),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "TN": tn,
            "recall": round(tp / len(positives), 4) if positives else None,
            "specificity": round(tn / len(negatives_all), 4) if negatives_all else None,
            "precision": round(tp / (tp + fp), 4) if (tp + fp) else None,
            "missed": sorted(positives - hit_ids),
            "false_positives": sorted(hit_ids & negatives_all)[:50],
        }

    baseline = evaluate(spec.rules)
    summary: dict[str, Any] = {
        "motif": spec.name,
        "structures_evaluated": int(n_files),
        "positives_declared": len(positives),
        "positives_without_anchor_site": sorted(missing),
        "baseline": baseline,
        "parameters": spec.resolved(),
    }

    if args.sweep:
        path, lo, hi, step = args.sweep.split(":")
        values = np.arange(float(lo), float(hi) + 1e-9, float(step))
        rows = []
        for v in values:
            res = evaluate(_sweep_rules(spec.rules, path, float(v)))
            rows.append({"parameter": path, "value": float(v), **{
                k: res[k] for k in ("hits", "TP", "FP", "FN", "recall", "specificity", "precision")
            }})
        sweep_df = pd.DataFrame(rows)
        sweep_csv = outdir / f"sweep_{path.replace('.', '_')}_{ts}.csv"
        sweep_df.to_csv(sweep_csv, index=False)
        summary["sweep_csv"] = str(sweep_csv)
        summary["sweep"] = rows
        print("\n" + sweep_df.to_string(index=False))

    _write_summary(outdir, ts, summary)
    if baseline["recall"] is not None and baseline["recall"] < 1.0:
        print(
            "\nRecall < 1.0 — the motif misses known members. Loosen geometry or "
            "revisit the mechanism before mining at scale.",
            file=sys.stderr,
        )


# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mcmine", description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("--motif", required=True, help="Path to motif spec JSON")
        sp.add_argument("--outdir", required=True)
        sp.add_argument("--max-structures", type=int, default=None)
        sp.add_argument("--progress-every", type=int, default=200)

    sp = sub.add_parser("mine", help="Classify sites against a motif spec")
    sp.add_argument("--struct-dir", required=True, help="Folder of <accession>.pdb/.cif files")
    sp.add_argument("--emit-shell", action="store_true", help="Also write per-site ligand table")
    common(sp)
    sp.set_defaults(func=cmd_mine)

    sp = sub.add_parser("type", help="Census of coordination environments (anchor-Xn)")
    sp.add_argument("--struct-dir", required=True)
    sp.add_argument("--ligand-cutoff", type=float, default=None, help="Override shell cutoff (A)")
    common(sp)
    sp.set_defaults(func=cmd_type)

    sp = sub.add_parser("benchmark", help="Recovery/specificity vs known positives")
    sp.add_argument("--struct-dir", help="Folder of structures (omit if --sites-csv given)")
    sp.add_argument("--sites-csv", help="Reuse a sites_all_*.csv from a previous run")
    sp.add_argument("--positives", required=True, help="File or comma list of positive accessions")
    sp.add_argument(
        "--sweep",
        help="Threshold scan, e.g. 'require_absent.0.min_distance:4.0:7.0:0.25'",
    )
    common(sp)
    sp.set_defaults(func=cmd_benchmark)

    return p


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
