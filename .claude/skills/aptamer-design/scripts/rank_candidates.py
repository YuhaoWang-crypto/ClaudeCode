#!/usr/bin/env python3
"""
Rank aptamer candidates from Boltz-2.1 co-folding metrics with a transparent,
use-case-weighted composite score. Standard library only.

Input: a JSON file — a list of candidate objects:
  {
    "id": "DNA-G4",
    "seq": "GGGAGGGTGGGAGGGTTCAGTCGACT",
    "chem": "DNA",                 # "DNA" or "RNA"
    "scaffold": "G4",              # G4 | hairpin | dual-hairpin | pseudoknot | random
    "iptm": 0.918,                 # required
    "ptm": 0.956,                  # optional
    "complex_iplddt": 0.755,       # optional
    "iface_res": 21,               # optional (interface residue count)
    "on_target_frac": 0.67,        # optional (fraction of interface on ligand footprint)
    "offtarget_iptm": 0.60         # optional (best paralog/off-target ipTM)
  }

Usage:
  python rank_candidates.py candidates.json --use-case diagnostic
  python rank_candidates.py candidates.json --use-case therapeutic --csv out.csv

Scoring (all sub-scores 0..1, weighted then scaled to 0..100):
  confidence  = iptm
  fold        = mean(ptm, complex_iplddt) present
  chemistry   = scaffold/chem suitability for the use case
  specificity = (iptm - offtarget_iptm) margin if provided, else neutral 0.5
  on_target   = on_target_frac (therapeutic only; discounted if iface_res small)

Weights:
  diagnostic : confidence .40  fold .20  chemistry .20  specificity .20
  therapeutic: confidence .40  fold .20  chemistry .10  specificity .10  on_target .20
"""
import argparse, csv, json, sys

CHEM_SCORE = {  # suitability by (use_case, chem, scaffold)
    "diagnostic": {"DNA_G4": 1.00, "DNA": 0.85, "RNA_pseudoknot": 0.65, "RNA": 0.55},
    "therapeutic": {"DNA_G4": 0.90, "DNA": 0.85, "RNA_pseudoknot": 0.80, "RNA": 0.75},
}

def chem_key(chem, scaffold):
    chem = (chem or "").upper()
    sca = (scaffold or "").lower().replace("-", "").replace("_", "")
    if chem == "DNA" and "g4" in sca or "quadruplex" in sca:
        return "DNA_G4"
    if chem == "DNA":
        return "DNA"
    if chem == "RNA" and "pseudoknot" in sca or sca == "pk":
        return "RNA_pseudoknot"
    return "RNA"

def clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))

def sub_scores(c, use_case):
    iptm = float(c["iptm"])
    conf = clamp(iptm)

    fold_vals = [c[k] for k in ("ptm", "complex_iplddt") if c.get(k) is not None]
    fold = sum(fold_vals) / len(fold_vals) if fold_vals else conf

    ck = chem_key(c.get("chem"), c.get("scaffold"))
    chem = CHEM_SCORE[use_case].get(ck, 0.5)

    if c.get("offtarget_iptm") is not None:
        # margin over best off-target, scaled: 0.2 gap -> ~full marks
        spec = clamp(0.5 + (iptm - float(c["offtarget_iptm"])) * 2.5)
    else:
        spec = 0.5  # neutral / unknown

    on_t = None
    if c.get("on_target_frac") is not None:
        frac = clamp(float(c["on_target_frac"]))
        # discount small-denominator interfaces (< 10 residues) toward 0.5
        n = c.get("iface_res")
        if n is not None and n < 10:
            frac = 0.5 + (frac - 0.5) * (n / 10.0)
        on_t = frac
    return conf, fold, chem, spec, on_t

def composite(c, use_case):
    conf, fold, chem, spec, on_t = sub_scores(c, use_case)
    if use_case == "diagnostic":
        score = 0.40*conf + 0.20*fold + 0.20*chem + 0.20*spec
    else:  # therapeutic
        ot = on_t if on_t is not None else 0.5
        score = 0.40*conf + 0.20*fold + 0.10*chem + 0.10*spec + 0.20*ot
    return round(score * 100, 1), dict(confidence=round(conf,3), fold=round(fold,3),
                                        chemistry=round(chem,3), specificity=round(spec,3),
                                        on_target=(round(on_t,3) if on_t is not None else None))

def grade(score):
    return ("A+" if score >= 85 else "A" if score >= 78 else "A-" if score >= 74
            else "B+" if score >= 70 else "B" if score >= 64 else "B-" if score >= 60
            else "C+" if score >= 55 else "C" if score >= 48 else "C-")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates_json")
    ap.add_argument("--use-case", choices=["diagnostic", "therapeutic"], default="diagnostic")
    ap.add_argument("--csv", help="also write ranked CSV here")
    args = ap.parse_args()

    with open(args.candidates_json) as fh:
        cands = json.load(fh)

    rows = []
    for c in cands:
        if "iptm" not in c:
            print(f"skip {c.get('id','?')}: no iptm", file=sys.stderr); continue
        score, parts = composite(c, args.use_case)
        rows.append({"id": c.get("id","?"), "chem": c.get("chem",""),
                     "scaffold": c.get("scaffold",""), "iptm": c.get("iptm"),
                     "score": score, "grade": grade(score), "seq": c.get("seq",""),
                     **{f"s_{k}": v for k, v in parts.items()}})
    rows.sort(key=lambda r: r["score"], reverse=True)

    print(f"\nAptamer ranking  (use-case: {args.use_case})")
    print("=" * 78)
    print(f"{'#':>2}  {'ID':<14} {'chem':<4} {'ipTM':>5} {'score':>6} {'grade':>5}  seq")
    print("-" * 78)
    for i, r in enumerate(rows, 1):
        print(f"{i:>2}  {r['id']:<14} {r['chem']:<4} {r['iptm']:>5} {r['score']:>6} "
              f"{r['grade']:>5}  {r['seq']}")
    print("=" * 78)
    print("Reminder: ipTM/score are PREDICTED confidence proxies, not affinity. "
          "Differences < ~3 points are within noise.")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {args.csv}")

if __name__ == "__main__":
    main()
