"""Score a target log2FC profile against the v0.2 library. Python standard library only.

    python -m cell_state_library.query target.json [--method nnls|pair|both]

Input JSON: {"cell_type": "CD4_Naive", "celltype_level": "L2",
             "contrast": "target vs matched control", "log2fc": {"ISG15": 1.2, ...}}
Accepted cell types: the 13 active subsets (author L1/L2 labels are resolved via
cell_type_aliases). Scores for a stimulus that fails the identifiability check
are withheld (value null), not forced to zero.
"""
import argparse
import json
import sqlite3
import sys
from contextlib import closing
from pathlib import Path

from .nnls import score_nnls
from .pair_code import score_pairs

DEFAULT_DB = Path(__file__).resolve().parent / "data" / "cell_state_v0_2.sqlite"
STIM = {"IFNA": "IFNa", "IFNG": "IFNg"}

INTERPRETATION = ("Similarity scores, not probabilities, cytokine concentrations, causal attribution or a validated "
                  "culture recipe. Restricted to author-selected IFN genes; does not assess the full cell state.")


def resolve_cell_type(con, label, level):
    row = con.execute("SELECT cell_type_id FROM cell_type_aliases WHERE alias_label=? AND alias_level=?",
                      (label, level)).fetchone()
    if row is None:
        row = con.execute("SELECT cell_type_id FROM cell_types WHERE cell_type_id=?", (label,)).fetchone()
    if row is None:
        raise ValueError(f"Unknown cell type {label!r} at level {level!r}.")
    return row[0]


def load(con, ct):
    statuses = {c: s for c, s in con.execute(
        "SELECT condition_id, status FROM reference_profiles WHERE cell_type_id=?", (ct,))}
    if not statuses or all(s.startswith("staged") for s in statuses.values()):
        raise ValueError(f"{ct} is staged, not active in v0.2. Active: 12 author L2 subsets + Monocyte.")
    refs = {}
    for cond in ("IFNA", "IFNG"):
        refs[cond] = dict(con.execute(
            "SELECT v.gene, v.log2fc FROM reference_values v JOIN reference_profiles p USING(profile_id) "
            "WHERE p.cell_type_id=? AND p.condition_id=?", (ct, cond)))
    ident = dict(zip(("status", "z_ifna", "z_ifng", "assumed_sigma", "notes"), con.execute(
        "SELECT status, z_ifna, z_ifng, assumed_sigma, notes FROM identifiability WHERE cell_type_id=?", (ct,)).fetchone()))
    panel = [dict(zip(("gene_i", "gene_j", "d_ifna", "d_ifng", "tau", "sigma"), r)) for r in con.execute(
        "SELECT gene_i, gene_j, d_ifna, d_ifng, tau, sigma FROM pair_panels WHERE cell_type_id=? ORDER BY rank", (ct,))]
    identifiable = [c for c, s in statuses.items() if s == "active"]
    return refs, ident, panel, identifiable


def score(query, db=DEFAULT_DB, method="both"):
    if not query.get("contrast") or not isinstance(query.get("log2fc"), dict):
        raise ValueError("Provide contrast and a gene-to-log2FC mapping, not expression counts.")
    with closing(sqlite3.connect(db)) as con:
        ct = resolve_cell_type(con, query.get("cell_type"), query.get("celltype_level"))
        refs, ident, panel, identifiable = load(con, ct)
    out = {"cell_type_id": ct, "identifiability": ident, "identifiable_stimuli": identifiable,
           "interpretation": INTERPRETATION}
    if identifiable == ["IFNA"]:
        out["warning"] = ("IFNg is not identifiable in this cell type: an IFNg exposure gives a small response that "
                          "looks like weak IFNa (in the v0.2 benchmark NNLS called simulated IFNg targets 'IFNa' in "
                          "~99-100% of cases). Read a positive IFNa score here as 'IFN response present, type unresolved'.")
    if method in ("nnls", "both"):
        r = score_nnls(refs, query["log2fc"])
        r["scores"] = {STIM[k]: {"value": v if k in identifiable else None,
                                 "raw_value": v, "identifiable": k in identifiable}
                       for k, v in r["scores"].items()}
        out["nnls"] = r
    if method in ("pair", "both"):
        sigma, tau = panel[0]["sigma"], panel[0]["tau"]
        r = score_pairs(panel, query["log2fc"], sigma, tau, identifiable)
        r["panel_genes"] = sorted({g for p in panel for g in (p["gene_i"], p["gene_j"])})
        r["note"] = ("Log-likelihoods assume iid per-gene noise sigma; they are model-based, not calibrated "
                     "probabilities. Classes limited to identifiable stimuli plus 'none'.")
        out["pair_code"] = r
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("query", type=Path)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--method", choices=("nnls", "pair", "both"), default="both")
    args = ap.parse_args(argv)
    if not args.db.exists():
        ap.exit(2, f"{args.db} not found; run `python -m cell_state_library.build` first.\n")
    try:
        res = score(json.loads(args.query.read_text()), args.db, args.method)
    except ValueError as e:
        ap.exit(2, str(e) + "\n")
    json.dump(res, sys.stdout, ensure_ascii=False, indent=2)
    print()


if __name__ == "__main__":
    main()
