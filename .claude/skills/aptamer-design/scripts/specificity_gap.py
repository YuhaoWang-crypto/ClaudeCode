#!/usr/bin/env python3
"""
Calibrated specificity gate for aptamer candidates.

Hard lesson (CTLA-4 cross-validation, this repo): absolute ipTM / HDOCK scores are NOT
trustworthy for small-protein x ssRNA — Boltz ipTM saturates ~0.8 and HDOCK docks
polyanionic RNA onto any (especially cationic) protein. A candidate that scores well on
its target can score EVEN BETTER on an unrelated protein (aptamerd6: HDOCK -313.8 on
CTLA-4 vs -382.1 on lysozyme; Boltz 0.51 on CTLA-4 vs 0.81 on lysozyme).

So we NEVER rank by absolute score. A candidate is "specific" only if, on EACH scorer:
  (1) on-target beats its matched scramble decoy on the SAME target (real signal), AND
  (2) on-target beats EVERY off-target (paralog + unrelated) by a margin.
If any off-target beats the on-target on any scorer -> PROMISCUOUS (reject).

Input JSON: list of candidates:
  {
    "id": "aptamerd6",
    "on_target": "CTLA-4",
    "boltz": {"CTLA-4":0.514,"CD28":0.851,"PD-L1":0.848,"Lysozyme":0.814},
    "boltz_decoy_on_target": 0.267,          # scramble vs on-target (optional)
    "hdock": {"CTLA-4":-313.8,"Lysozyme":-382.06}   # more-negative=stronger (optional)
  }

Usage: python specificity_gap.py candidates.json [--margin 0.05]
"""
import argparse, json, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates_json")
    ap.add_argument("--margin", type=float, default=0.05,
                    help="Boltz ipTM margin on-target must beat every off-target by")
    args = ap.parse_args()
    cands = json.load(open(args.candidates_json))

    print(f"\nSpecificity gate  (on-target must beat decoy AND all off-targets)")
    print("=" * 78)
    for c in cands:
        tgt = c["on_target"]; flags = []; verdict = "PASS"
        # --- Boltz ---
        b = c.get("boltz", {})
        if tgt in b:
            on = b[tgt]
            offs = {k: v for k, v in b.items() if k != tgt}
            worst_off = max(offs.values()) if offs else None      # highest off-target ipTM
            dec = c.get("boltz_decoy_on_target")
            if dec is not None:
                gap = on - dec
                flags.append(f"Boltz on-vs-decoy {gap:+.2f}")
                if gap <= 0.03: verdict = "FAIL"; flags.append("no-decoy-gap")
            if worst_off is not None:
                sm = on - worst_off
                topoff = max(offs, key=offs.get)
                flags.append(f"Boltz spec {sm:+.2f} (worst off={topoff}:{worst_off:.2f})")
                if sm < args.margin: verdict = "FAIL"; flags.append("Boltz-PROMISCUOUS")
        # --- HDOCK (more negative = stronger) ---
        h = c.get("hdock", {})
        if tgt in h:
            on = h[tgt]
            offs = {k: v for k, v in h.items() if k != tgt}
            best_off = min(offs.values()) if offs else None       # most-negative off-target
            if best_off is not None:
                sm = best_off - on          # >0 means on-target is MORE negative (better)
                topoff = min(offs, key=offs.get)
                flags.append(f"HDOCK spec {sm:+.1f} (best off={topoff}:{best_off:.0f})")
                if sm <= 0: verdict = "FAIL"; flags.append("HDOCK-PROMISCUOUS")
        mark = "✅" if verdict == "PASS" else "❌"
        print(f"{mark} {c['id']:<14} target={tgt:<8} {verdict}")
        for f in flags: print(f"      - {f}")
    print("=" * 78)
    print("PASS requires positive decoy gap AND on-target beating every off-target on")
    print("every scorer. A single off-target win = PROMISCUOUS -> back to counter-SELEX.")

if __name__ == "__main__":
    main()
