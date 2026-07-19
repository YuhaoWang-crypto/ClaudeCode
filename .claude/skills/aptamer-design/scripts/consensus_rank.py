#!/usr/bin/env python3
"""
Consensus ranking of aptamer candidates from TWO orthogonal 3D scorers:
  - HDOCK (T-SELEX path): docking score (more negative = better) + confidence (0-1)
  - Boltz-2.1: interface ipTM (higher = better)

Combines them by RANK (Borda), because HDOCK score and Boltz ipTM live on different
scales and must not be added directly. A candidate ranked well by BOTH methods is far
more trustworthy than one that only one method likes.

Optional specificity: pass off-target values to require on-target >> off-target in BOTH
scorers (the paralog counter-screen that decides diagnostic suitability).

Input JSON: list of candidates, each:
  {
    "id": "PDL1-RNA-HP", "chem": "RNA",
    "hdock_score": -285.4,          # optional (RNA only; DNA may omit)
    "hdock_confidence": 0.83,       # optional
    "boltz_iptm": 0.771,            # optional but recommended
    "boltz_offtarget_iptm": 0.166,  # optional (specificity)
    "hdock_offtarget_score": -180.0 # optional (specificity)
  }

Usage:
  python consensus_rank.py candidates.json
  python consensus_rank.py candidates.json --decoy-iptm 0.63   # flag on-target below decoy
"""
import argparse, json, sys

def rank_map(items, key, better):
    """Return {id: rank} (1 = best). better='low' or 'high'. Missing values -> None."""
    vals = [(c["id"], c.get(key)) for c in items if c.get(key) is not None]
    vals.sort(key=lambda t: t[1], reverse=(better == "high"))
    return {cid: i + 1 for i, (cid, _) in enumerate(vals)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates_json")
    ap.add_argument("--decoy-iptm", type=float, default=None,
                    help="scramble-decoy on-target ipTM; flags candidates that don't clear it")
    args = ap.parse_args()
    with open(args.candidates_json) as fh:
        cands = json.load(fh)

    hdock_rank = rank_map(cands, "hdock_score", "low")     # more negative = better
    boltz_rank = rank_map(cands, "boltz_iptm", "high")     # higher = better
    n = len(cands)

    rows = []
    for c in cands:
        hr = hdock_rank.get(c["id"]); br = boltz_rank.get(c["id"])
        present = [r for r in (hr, br) if r is not None]
        consensus = sum(present) / len(present) if present else n + 1
        both = hr is not None and br is not None
        agree = both and hr <= n / 2 and br <= n / 2      # both in top half
        # specificity margins (positive = on-target preferred)
        spec = None
        if c.get("boltz_iptm") is not None and c.get("boltz_offtarget_iptm") is not None:
            spec_b = c["boltz_iptm"] - c["boltz_offtarget_iptm"]
        else:
            spec_b = None
        if c.get("hdock_score") is not None and c.get("hdock_offtarget_score") is not None:
            spec_h = c["hdock_offtarget_score"] - c["hdock_score"]   # off higher(worse)=good
        else:
            spec_h = None
        # flags + gate penalty (disqualifiers sink to the bottom)
        flags = []; penalty = 0.0
        if not both:
            flags.append("single-scorer")
        elif agree:
            flags.append("consensus✓")
        else:
            flags.append("disagree")
        if args.decoy_iptm is not None and c.get("boltz_iptm") is not None \
                and c["boltz_iptm"] < args.decoy_iptm:
            flags.append("below-decoy"); penalty += 100
        if spec_b is not None:
            flags.append(f"specB{spec_b:+.2f}")
            if spec_b <= 0.03:
                flags.append("NON-SPECIFIC"); penalty += 100
        rows.append((consensus + penalty, c, hr, br, spec_b, spec_h, flags))

    rows.sort(key=lambda r: r[0])
    print(f"\nConsensus ranking (HDOCK ⊕ Boltz, Borda)   n={n}")
    print("=" * 84)
    print(f"{'#':>2}  {'id':<15} {'chem':<4} {'HDOCKr':>6} {'Boltzr':>6} {'cons':>5}  flags")
    print("-" * 84)
    for i, (cons, c, hr, br, sb, sh, flags) in enumerate(rows, 1):
        print(f"{i:>2}  {c['id']:<15} {c.get('chem',''):<4} "
              f"{('-' if hr is None else hr):>6} {('-' if br is None else br):>6} "
              f"{cons:>5.1f}  {', '.join(flags)}")
    print("=" * 84)
    print("Lead = low consensus rank + consensus✓ + positive specificity + clears decoy.")
    print("Scores are relative confidence proxies, NOT affinities. Confirm by wet-lab SELEX.")

if __name__ == "__main__":
    main()
