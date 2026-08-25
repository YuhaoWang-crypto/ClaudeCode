#!/usr/bin/env python3
"""
M3 - HLA-DR binding / presentation prediction across the designed panel.

Two orthogonal NetMHCIIpan heads are run over every 15-mer:

  netmhciipan_el  eluted-ligand likelihood - what actually gets presented
  netmhciipan_ba  binding affinity          - peptide-MHC stability

A single-head EL call is the industry default and is what the previous report
used. It over-calls: EL is trained on mass-spec ligandomes and rewards
motif-like peptides even when the measured affinity is poor. Requiring a
*consensus* call (EL %Rank < 1 AND BA %Rank < 10) removes that class of hit
while keeping every peptide that both heads agree on. Both the single-head and
the consensus call are written out so the effect is auditable.

Output: results/m3_binding_long.tsv, one row per (sequence, peptide, allele).
"""
import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, data_path, results_path,  # noqa: E402
                    iedb_mhcii)

MAX_WORKERS = 4          # be a good citizen on the shared IEDB cluster


def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def main():
    cfg = load_config()
    pcfg = cfg["prediction"]
    seqs = read_fasta(data_path("sequences.fasta"))
    with open(results_path("m2_panel_alleles.txt")) as f:
        panel = [l.strip() for l in f if l.strip()]

    # Resume: anything already in the long table is not fetched again, so
    # widening the panel costs only the new allele's calls.
    done = {"netmhciipan_el": set(), "netmhciipan_ba": set()}
    existing = []
    out_path = results_path("m3_binding_long.tsv")
    if os.path.exists(out_path):
        with open(out_path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                existing.append(r)
                if r["el_rank"] not in ("", "None"):
                    done["netmhciipan_el"].add((r["id"], r["allele"]))
                if r["ba_rank"] not in ("", "None"):
                    done["netmhciipan_ba"].add((r["id"], r["allele"]))
        print(f"resuming: {len(existing)} rows already on disk")

    jobs = []
    for method in ("netmhciipan_el", "netmhciipan_ba"):
        for block in chunks(panel, pcfg["alleles_per_request"]):
            for sid, seq in seqs.items():
                todo = [a for a in block if (sid, a) not in done[method]]
                if todo:
                    jobs.append((method, todo, sid, seq))

    print(f"{len(seqs)} sequences x {len(panel)} DR molecules x 2 heads "
          f"-> {len(jobs)} IEDB calls")

    def run(job):
        method, block, sid, seq = job
        rows = iedb_mhcii(pcfg["endpoint"], method, {sid: seq}, block,
                          pcfg["peptide_length"])
        for r in rows:
            r["method"] = method
        print(f"  ok {method:16s} {sid:22s} {len(block)} alleles "
              f"-> {len(rows):5d} rows", flush=True)
        return rows

    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for rows in ex.map(run, jobs):
            results.extend(rows)

    # EL and BA name their score column differently; normalise to score/rank.
    merged = {}
    for r in existing:
        merged[(r["id"], r["allele"], int(r["start"]))] = {
            "id": r["id"], "allele": r["allele"], "start": int(r["start"]),
            "end": int(r["end"]), "peptide": r["peptide"], "core": r["core"],
            "el_rank": None if r["el_rank"] in ("", "None") else float(r["el_rank"]),
            "el_score": None if r["el_score"] in ("", "None") else float(r["el_score"]),
            "ba_rank": None if r["ba_rank"] in ("", "None") else float(r["ba_rank"]),
            "ba_ic50": None if r["ba_ic50"] in ("", "None") else float(r["ba_ic50"]),
        }
    for r in results:
        key = (r["id"], r["allele"], int(r["start"]))
        rec = merged.setdefault(key, {
            "id": r["id"], "allele": r["allele"], "start": int(r["start"]),
            "end": int(r["end"]), "peptide": r["peptide"],
            "core": r["core_peptide"], "el_rank": None, "el_score": None,
            "ba_rank": None, "ba_ic50": None,
        })
        if r["method"] == "netmhciipan_el":
            rec["el_rank"] = float(r["rank"])
            rec["el_score"] = float(r["score"])
            rec["core"] = r["core_peptide"]
        else:
            rec["ba_rank"] = float(r["rank"])
            rec["ba_ic50"] = float(r["ic50"])

    sb, wb = pcfg["sb_rank"], pcfg["wb_rank"]
    bac = pcfg["ba_confirm_rank"]
    out = out_path
    cols = ["id", "allele", "start", "end", "peptide", "core", "el_rank",
            "el_score", "ba_rank", "ba_ic50", "call_el", "call_consensus"]
    n_sb = n_cons = 0
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for rec in sorted(merged.values(), key=lambda r: (r["id"], r["start"], r["allele"])):
            el = rec["el_rank"]
            ba = rec["ba_rank"]
            rec["call_el"] = "SB" if el is not None and el < sb else \
                             "WB" if el is not None and el < wb else "-"
            if rec["call_el"] == "SB" and ba is not None and ba < bac:
                rec["call_consensus"] = "SB"
            elif rec["call_el"] in ("SB", "WB") and ba is not None and ba < bac:
                rec["call_consensus"] = "WB"
            else:
                rec["call_consensus"] = "-"
            n_sb += rec["call_el"] == "SB"
            n_cons += rec["call_consensus"] == "SB"
            w.writerow(rec)

    print(f"\nwrote {len(merged)} peptide-allele rows to {out}")
    print(f"EL-only strong binders : {n_sb}")
    print(f"EL+BA consensus strong : {n_cons}  "
          f"({100*(n_sb-n_cons)/max(n_sb,1):.1f}% of EL calls dropped by BA)")


if __name__ == "__main__":
    main()
