#!/usr/bin/env python3
"""Driver for structure-first catalytic-scaffold mining (Stages 1-4).

Two modes, same filters:

  web    -- Foldseek web server against the AlphaFold Database. No downloads.
            Calpha-level geometry only; no pLDDT.
  local  -- local Foldseek binary against a directory of structures (e.g. one
            extracted AlphaFold proteome). Full side chains, so acidic-cluster
            geometry and catalytic-shell pLDDT are available.

Example
-------
    python3 mine.py \
        --query ref/4km5.pdb --chain A --catalytic 225,227,319 \
        --mode web --databases afdb-proteome afdb-swissprot \
        --out out/cgas

    python3 mine.py \
        --query ref/4km5.pdb --chain A --catalytic 225,227,319 \
        --mode local --target-dir proteomes/YEAST --out out/yeast

Every row is a *nominated* candidate. Nothing this script prints is evidence
of catalytic activity -- see SKILL.md and reference/validation.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from catalytic import CatalyticSite, check_site, residue_map   # noqa: E402
from fsweb import (DEFAULT_DATABASES, iter_hits, search,        # noqa: E402
                   target_structure, tm_score)
from structure import read_structure                            # noqa: E402

FOLDSEEK_FIELDS = [
    "query", "target", "qstart", "qend", "tstart", "tend", "qaln", "taln",
    "evalue", "prob", "alntmscore", "qtmscore", "ttmscore", "lddt", "fident",
    "alnlen", "qcov", "tcov", "qlen", "tlen",
]


def run_local(query_path, target_dir, out_prefix, foldseek="foldseek",
              evalue=10.0, max_seqs=10000, threads=8):
    """foldseek easy-search; returns the parsed hit dicts."""
    raw = f"{out_prefix}.foldseek.tsv"
    tmp = f"{out_prefix}.fstmp"
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    cmd = [
        foldseek, "easy-search", query_path, target_dir, raw, tmp,
        "--format-output", ",".join(FOLDSEEK_FIELDS),
        "-e", str(evalue), "--max-seqs", str(max_seqs),
        "--threads", str(threads), "-v", "1",
    ]
    subprocess.run(cmd, check=True)
    hits = []
    with open(raw) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(FOLDSEEK_FIELDS):
                continue
            h = dict(zip(FOLDSEEK_FIELDS, parts))
            for k in ("qstart", "qend", "tstart", "tend", "alnlen", "qlen", "tlen"):
                h[k] = int(float(h[k]))
            for k in ("evalue", "prob", "alntmscore", "qtmscore", "ttmscore",
                      "lddt", "fident", "qcov", "tcov"):
                h[k] = float(h[k])
            hits.append(h)
    return hits


def _accession_from_filename(name):
    base = os.path.basename(name)
    for suffix in (".gz", ".cif", ".pdb", ".mmcif"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    parts = base.split("-")
    if len(parts) >= 2 and parts[0] == "AF":
        return parts[1]
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="query structure (PDB/mmCIF)")
    ap.add_argument("--chain", default=None, help="query chain (default: first)")
    ap.add_argument("--catalytic", required=True,
                    help="comma-separated author residue numbers in the query")
    ap.add_argument("--class", dest="klass", default="acidic",
                    help="required chemical class at each catalytic position")
    ap.add_argument("--mode", choices=["web", "local"], default="web")
    ap.add_argument("--databases", nargs="+", default=list(DEFAULT_DATABASES))
    ap.add_argument("--target-dir", help="local mode: directory of structures")
    ap.add_argument("--foldseek", default="foldseek", help="foldseek binary")
    ap.add_argument("--shell-radius", type=float, default=12.0,
                    help="catalytic shell radius in Angstrom (default 12)")
    ap.add_argument("--window", type=int, default=0,
                    help="allow catalytic counterparts +/- N residues from the "
                         "aligned position (default 0 = strict; 2 is a sane "
                         "setting below ~20%% identity). Offsets are reported.")
    ap.add_argument("--require-all", action="store_true",
                    help="keep only hits conserving every catalytic position")
    ap.add_argument("--min-conserved", type=int, default=0,
                    help="keep hits with at least this many conserved positions")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse an existing raw result file if present")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    query = read_structure(args.query, chain=args.chain)
    resnums = [int(x) for x in args.catalytic.split(",") if x.strip()]
    site = CatalyticSite(resnums=resnums, klass=args.klass)

    missing = [rn for rn in resnums if query.by_resnum(rn) is None]
    if missing:
        sys.exit(f"catalytic residues {missing} not present in {args.query}")
    found = "".join(query.by_resnum(rn).one for rn in resnums)
    print(f"query {args.query} chain={query.chain} len={len(query)} "
          f"catalytic={resnums} ({found}) class={args.klass}")

    rows = []

    if args.mode == "web":
        raw_path = f"{args.out}.raw.json"
        if args.reuse and os.path.exists(raw_path):
            print(f"  reusing {raw_path}")
            raw = json.load(open(raw_path))
        else:
            raw = search(args.query, raw_path, databases=args.databases)
        hits = list(iter_hits(raw))
        print(f"  {len(hits)} raw structural hits")
        for h in hits:
            if not h["qaln"] or not h["taln"]:
                continue
            target = target_structure(h)
            mapping = residue_map(h["qstart"], h["tstart"], h["qaln"], h["taln"])
            chk = check_site(query, target, site, mapping,
                             shell_radius=args.shell_radius, window=args.window)
            rows.append(_row(h, chk, query, target, mapping, mode="web"))
    else:
        if not args.target_dir:
            sys.exit("--target-dir is required in local mode")
        hits = run_local(args.query, args.target_dir, args.out,
                         foldseek=args.foldseek)
        print(f"  {len(hits)} raw structural hits")
        for h in hits:
            path = os.path.join(args.target_dir, h["target"])
            if not os.path.exists(path):
                cand = [p for p in os.listdir(args.target_dir)
                        if p.startswith(h["target"])]
                if not cand:
                    continue
                path = os.path.join(args.target_dir, cand[0])
            try:
                target = read_structure(path)
            except Exception as exc:                      # noqa: BLE001
                print(f"    skip {h['target']}: {exc}")
                continue
            mapping = residue_map(h["qstart"], h["tstart"], h["qaln"], h["taln"])
            chk = check_site(query, target, site, mapping,
                             shell_radius=args.shell_radius, window=args.window)
            h["accession"] = _accession_from_filename(h["target"])
            rows.append(_row(h, chk, query, target, mapping, mode="local"))

    if args.require_all:
        rows = [r for r in rows if r["site_passed"]]
    if args.min_conserved:
        rows = [r for r in rows if r["n_conserved"] >= args.min_conserved]

    rows.sort(key=lambda r: (-(r["n_conserved"]), r["core_rmsd"]
                             if r["core_rmsd"] == r["core_rmsd"] else 9e9))

    tsv = f"{args.out}.candidates.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["accession"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(f"{args.out}.candidates.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    print(f"  wrote {len(rows)} rows -> {tsv}")
    print("  NOTE: every row is a nominated candidate, not a demonstrated "
          "activity (see reference/validation.md)")


def _row(h, chk, query, target, mapping, mode):
    tm = h.get("qtmscore")
    if tm is None:
        tm = tm_score(query, target, mapping, norm_len=len(query))
    return {
        "accession": h.get("accession", h.get("target", "")),
        "description": h.get("description", ""),
        "tax_name": h.get("tax_name", ""),
        "tax_id": h.get("tax_id", ""),
        "db": h.get("db", mode),
        "evalue": h.get("evalue", h.get("eval")),
        "prob": h.get("prob"),
        "seq_identity": h.get("fident", (h.get("seq_id") or 0) / 100.0
                              if h.get("seq_id") is not None else None),
        "tm_score_qnorm": round(float(tm), 4) if tm == tm else None,
        "aln_length": h.get("alnlen", h.get("aln_length")),
        "query_len": h.get("qlen", h.get("qlen")),
        "target_len": h.get("tlen", h.get("tlen")),
        "n_conserved": chk.n_conserved,
        "n_conserved_exact": chk.n_conserved_exact,
        "max_offset": chk.max_offset,
        "n_catalytic": chk.n_total,
        "site_passed": chk.passed,
        "catalytic_residues": ";".join(
            f"{p['query_aa']}{p['query_resnum']}->"
            f"{p['target_aa'] or '-'}{p['target_resnum'] or '-'}"
            + (f"[{p['offset']:+d}]" if p.get("offset") else "")
            for p in chk.positions),
        "core_rmsd": round(chk.core_rmsd, 3) if chk.core_rmsd == chk.core_rmsd
        else float("nan"),
        "n_core": chk.n_core,
        "core_frac": round(chk.core_frac, 2),
        "local_rmsd": round(chk.local_rmsd, 3) if chk.local_rmsd == chk.local_rmsd
        else float("nan"),
        "n_shell": chk.n_shell,
        "shell_plddt": round(chk.shell_plddt, 1) if chk.shell_plddt else None,
        "carboxylate_geometry": ";".join(
            f"{d['pair'][0]}-{d['pair'][1]}:{d['target_dist']}({d['delta']:+})"
            for d in chk.carboxylate_deltas),
        "label": "hypothesis",
    }


if __name__ == "__main__":
    main()
