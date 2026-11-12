#!/usr/bin/env python3
"""Deliverable: candidate sequences as FASTA + an EU-numbered edit table.

The edit table is the transferable artefact: it states each design as an
operation on EU-numbered positions, so it applies unchanged to the client's own
IgG1 heavy chain once those sequences are provided. The FASTA files are the
modelled modules (Fv or CH3 homodimer), not full-length IgG1.
"""
import json
import os

cand = json.load(open("results/03_candidates.json"))
os.makedirs("results/sequences", exist_ok=True)

CH3 = cand["scaffold"]["CH3_341_447"]
VH = cand["scaffold"]["VH"]
VL = cand["scaffold"]["VL"]
CH3_LO = 341
SITES = {"CH3-AB": (355, 362), "CH3-CD": (383, 391), "CH3-EF": (413, 422)}
CDRH3 = "WGGDGFYAMDY"

try:
    rank = {d["id"]: d for d in json.load(open("results/10_ranking.json"))["designs_ranked"]}
except FileNotFoundError:
    rank = {}

fa, rows = [], []
for c in cand["candidates"] + cand["controls"]:
    cid = c["id"]
    r = rank.get(cid, {})
    tag = ""
    if r:
        tag = f" rank={r['rank']} score={r['score_total']} {'SELECTED' if r.get('selected_top12') else 'not-selected'}"
    for ch, seq in c["chains"].items():
        fa.append(f">{cid}|chain_{ch}|{c.get('module','')}|{len(seq)}aa{tag}\n" +
                  "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)))

    if cid.startswith("C"):
        site = c["site"]
        if site == "CDRH3":
            where = (f"heavy chain CDRH3 ({CDRH3}); "
                     + ("insert after CDRH3 position 4 (apex)"
                        if c["strategy"] != "S3_replace" else "replace CDRH3 in full"))
        else:
            lo, hi = SITES[site]
            native = CH3[lo - CH3_LO:hi - CH3_LO + 1]
            if c["strategy"] == "S3_replace":
                where = f"replace EU {lo}-{hi} ({native}) in full"
            else:
                h = len(native) // 2
                where = (f"insert after EU {lo + h - 1} inside the {site} loop "
                         f"EU {lo}-{hi} ({native}); native loop retained")
        rows.append({
            "id": cid, "site": site, "strategy": c["strategy"],
            "module_modelled": c["module"],
            "edit": where,
            "inserted_cassette": c["insert"],
            "cassette_length_aa": c["insert_len"],
            "n_junction": c["n_junction"], "payload": c["payload"],
            "c_junction": c["c_junction"],
            "net_length_change_aa": c["length_change_vs_wt"],
            "rank": r.get("rank"), "score_total": r.get("score_total"),
            "selected_top12": r.get("selected_top12"),
        })

open("results/sequences/all_constructs.fasta", "w").write("\n".join(fa) + "\n")
json.dump({"scaffold_assumption": cand["scaffold"]["assumption_note"],
           "reference_scaffold": cand["scaffold"]["antibody"],
           "edits": rows}, open("results/sequences/edit_table.json", "w"), indent=2)

# human-readable TSV
with open("results/sequences/edit_table.tsv", "w") as fh:
    cols = ["id", "site", "strategy", "edit", "n_junction", "payload", "c_junction",
            "cassette_length_aa", "net_length_change_aa", "rank", "score_total",
            "selected_top12"]
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda x: (x["rank"] is None, x["rank"] or 0)):
        fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

print(f"wrote results/sequences/all_constructs.fasta ({len(fa)} records)")
print("wrote results/sequences/edit_table.json and .tsv")
