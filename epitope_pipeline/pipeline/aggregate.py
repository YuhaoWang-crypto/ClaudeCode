"""Merge per-model results into one normalized table + a consensus summary."""

from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from typing import Dict, List

from .runner import ModelResult

# %Rank thresholds for calling a "binder" (DTU convention)
SB_RANK = {"mhc_i": 0.5, "mhc_ii": 1.0}
WB_RANK = {"mhc_i": 2.0, "mhc_ii": 5.0}


def all_rows(results: List[ModelResult]) -> List[dict]:
    rows = []
    for res in results:
        for rec in res.records:
            rows.append(rec.as_row())
    return rows


def write_long_table(results: List[ModelResult], path: str) -> str:
    rows = all_rows(results)
    fields = ["model", "category", "seq_id", "position", "allele", "peptide",
              "core", "score", "rank_pct", "affinity_nM", "bind_level"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def _is_binder(category: str, rank) -> str:
    if rank is None:
        return ""
    if category in SB_RANK and rank <= SB_RANK[category]:
        return "SB"
    if category in WB_RANK and rank <= WB_RANK[category]:
        return "WB"
    return ""


def consensus(results: List[ModelResult], top_n: int = 25) -> List[dict]:
    """Group binder-level hits by peptide and count how many models flag it.

    A peptide flagged by several models (e.g. NetMHCpan + NetCTLpan) or across
    categories (MHC-I + B-cell region) is a higher-confidence immunogenic hit.
    """
    agg: Dict[str, dict] = defaultdict(lambda: {
        "peptide": "", "seq_id": "", "models": set(), "categories": set(),
        "alleles": set(), "best_rank": None, "best_affinity": None,
        "positions": set(),
    })
    for res in results:
        for rec in res.records:
            level = rec.bind_level or _is_binder(rec.category, rec.rank)
            if level not in ("SB", "WB") and rec.category not in (
                    "bcell_linear", "bcell_conf"):
                continue
            if rec.category in ("bcell_linear",) and (rec.bind_level != "epitope"):
                continue
            if not rec.peptide:
                continue
            key = f"{rec.seq_id}|{rec.peptide}"
            e = agg[key]
            e["peptide"] = rec.peptide
            e["seq_id"] = rec.seq_id
            e["models"].add(rec.model)
            e["categories"].add(rec.category)
            if rec.allele:
                e["alleles"].add(rec.allele)
            if rec.position is not None:
                e["positions"].add(rec.position)
            if rec.rank is not None:
                e["best_rank"] = (rec.rank if e["best_rank"] is None
                                  else min(e["best_rank"], rec.rank))
            if rec.affinity_nM is not None:
                e["best_affinity"] = (rec.affinity_nM if e["best_affinity"] is None
                                      else min(e["best_affinity"], rec.affinity_nM))

    out = []
    for e in agg.values():
        out.append({
            "peptide": e["peptide"],
            "seq_id": e["seq_id"],
            "n_models": len(e["models"]),
            "models": ",".join(sorted(e["models"])),
            "categories": ",".join(sorted(e["categories"])),
            "n_alleles": len(e["alleles"]),
            "alleles": ",".join(sorted(e["alleles"])),
            "best_rank_pct": e["best_rank"],
            "best_affinity_nM": e["best_affinity"],
        })
    # rank: more models first, then more categories, then stronger rank
    out.sort(key=lambda r: (-r["n_models"], -len(r["categories"].split(",")),
                            r["best_rank_pct"] if r["best_rank_pct"] is not None
                            else 1e9))
    return out[:top_n]


def write_consensus(results: List[ModelResult], path: str, top_n: int = 25) -> str:
    rows = consensus(results, top_n)
    fields = ["peptide", "seq_id", "n_models", "models", "categories",
              "n_alleles", "alleles", "best_rank_pct", "best_affinity_nM"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def write_summary_json(results: List[ModelResult], path: str) -> str:
    summary = {
        "models": [],
        "totals": {"records": 0, "ok": 0, "skipped": 0, "errors": 0},
    }
    for res in results:
        summary["models"].append({
            "model": res.model, "category": res.category, "status": res.status,
            "records": len(res.records), "message": res.message,
            "duration_s": round(res.duration_s, 2), "command": res.command,
        })
        summary["totals"]["records"] += len(res.records)
        if res.status == "ok":
            summary["totals"]["ok"] += 1
        elif res.status in ("not_installed", "needs_structure"):
            summary["totals"]["skipped"] += 1
        else:
            summary["totals"]["errors"] += 1
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)
    return path


def print_report(results: List[ModelResult], cons: List[dict]) -> None:
    print("\n=== Model run summary ===")
    for res in results:
        icon = {"ok": "✅", "not_installed": "⏭️ ", "needs_structure": "⏭️ ",
                "error": "❌"}.get(res.status, "? ")
        print(f"{icon} {res.model:<18} [{res.category:<12}] "
              f"{res.status:<15} {res.message}")
    print("\n=== Top consensus epitopes (flagged by >=1 model) ===")
    if not cons:
        print("  (no binder-level hits — either no tools installed or no binders)")
        return
    print(f"{'peptide':<18}{'seq':<10}{'#mdl':<5}{'cats':<22}"
          f"{'best%Rank':<10}{'models'}")
    for r in cons[:15]:
        rank = "" if r["best_rank_pct"] is None else f"{r['best_rank_pct']:.3f}"
        print(f"{r['peptide']:<18}{r['seq_id'][:9]:<10}{r['n_models']:<5}"
              f"{r['categories'][:21]:<22}{rank:<10}{r['models']}")
