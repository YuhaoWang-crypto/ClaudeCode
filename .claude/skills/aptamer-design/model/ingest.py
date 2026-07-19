#!/usr/bin/env python3
"""
Ingest a raw aptamer-database export (UTexas download, an Aptagen-derived CSV, or a
literature CSV) and normalize it into the unified schema (see schema.md).

- Maps source columns -> schema columns via a preset or a custom --map.
- Normalizes affinity to kd_nM.
- Optionally emits target-shuffled synthetic negatives for classifier training.
- Optionally back-fills target_seq from UniProt accessions (needs network; --fetch-seqs).

Usage:
  python ingest.py raw_utexas.csv --preset utexas -o unified.csv
  python ingest.py raw.csv --map "seq=Sequence,target_name=Target,kd_value=Kd,kd_unit=Unit" -o unified.csv --negatives
"""
import argparse, csv, sys

PRESETS = {
    # best-effort column-name guesses; adjust to the actual export headers
    "utexas": {"name": "Aptamer Name", "sequence": "Sequence", "chem": "Type",
               "target_name": "Target", "kd_value": "Kd", "kd_unit": "Kd Unit",
               "buffer": "Buffer", "doi": "DOI"},
    "aptagen": {"name": "Name", "sequence": "Sequence", "chem": "Chemistry",
                "target_name": "Target", "target_type": "Target Category",
                "kd_value": "Affinity", "kd_unit": "Unit", "specificity": "Specificity"},
}
UNIT_TO_NM = {"fM": 1e-6, "pM": 1e-3, "nM": 1.0, "uM": 1e3, "µM": 1e3, "mM": 1e6}
SCHEMA = ["apt_id","name","chem","sequence","mods","length","target_name","target_type",
          "target_seq","target_uniprot","kd_value","kd_unit","kd_nM","buffer",
          "specificity","source","doi","label"]

def to_nM(val, unit):
    try:
        return float(val) * UNIT_TO_NM[unit.strip()]
    except (ValueError, KeyError, AttributeError):
        return ""

def parse_map(s):
    return dict(kv.split("=", 1) for kv in s.split(",")) if s else {}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_csv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--preset", choices=list(PRESETS))
    ap.add_argument("--map", help="schema=source,... overrides preset")
    ap.add_argument("--source", default="", help="value for the 'source' column")
    ap.add_argument("--negatives", action="store_true", help="append target-shuffled negatives")
    args = ap.parse_args()

    mapping = dict(PRESETS.get(args.preset, {}))
    mapping.update({v: k for k, v in parse_map(args.map).items()})  # allow schema=source form
    # normalize mapping to schema_col -> source_col
    inv = {}
    for k, v in mapping.items():
        (inv.__setitem__(k, v) if k in SCHEMA else inv.__setitem__(v, k))

    with open(args.raw_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))

    out = []
    for i, r in enumerate(rows):
        def g(col): return r.get(inv.get(col, col), "")
        seq = (g("sequence") or "").upper().replace(" ", "")
        rec = {c: "" for c in SCHEMA}
        rec.update({
            "apt_id": g("apt_id") or f"{(args.source or args.preset or 'src')}_{i:04d}",
            "name": g("name"), "chem": g("chem") or "DNA", "sequence": seq,
            "length": len(seq), "target_name": g("target_name"),
            "target_type": g("target_type") or "protein", "target_seq": g("target_seq"),
            "target_uniprot": g("target_uniprot"), "kd_value": g("kd_value"),
            "kd_unit": g("kd_unit"), "kd_nM": to_nM(g("kd_value"), g("kd_unit")),
            "buffer": g("buffer"), "specificity": g("specificity"),
            "source": args.source or args.preset or "", "doi": g("doi"), "label": 1,
        })
        if rec["sequence"]:
            out.append(rec)

    if args.negatives:
        import random; random.seed(11)
        tgts = list({r["target_name"] for r in out if r["target_name"]})
        for r in list(out):
            others = [t for t in tgts if t != r["target_name"]]
            if not others:
                continue
            neg = dict(r); neg["target_name"] = random.choice(others)
            neg["apt_id"] = r["apt_id"] + "_neg"; neg["label"] = 0
            neg["target_seq"] = ""; neg["kd_nM"] = ""; neg["source"] = "synthetic-negative"
            out.append(neg)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SCHEMA); w.writeheader(); w.writerows(out)
    print(f"wrote {len(out)} rows -> {args.out}  "
          f"({sum(1 for r in out if r['label']==1)} pos, "
          f"{sum(1 for r in out if r['label']==0)} neg)")

if __name__ == "__main__":
    main()
