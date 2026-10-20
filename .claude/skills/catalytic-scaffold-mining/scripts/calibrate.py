#!/usr/bin/env python3
"""Stage 2: derive score thresholds from known family members, don't import them.

Foldseek reports no DALI Z-score, and published Z cutoffs do not transfer.
The defensible move is to take proteins whose membership in the family is
already established, see where they land on *your* query with *your* tool, and
set the cutoff from that distribution -- then report the distribution, not
just the cutoff.

    python3 calibrate.py out/urt1_proteome.annotated.tsv \
        --positives O64642 A0A0P0VE25 A0A0P0XSC2 \
        --negatives O75072 Q8IYB1

Positives are known members; negatives are proteins that share the fold but
are known to do a different job (the most informative negatives you can pick).
"""

from __future__ import annotations

import argparse
import csv
import statistics as st


def _f(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def describe(label, vals):
    vals = [v for v in vals if v == v]
    if not vals:
        return f"  {label:12s} (none found)"
    return (f"  {label:12s} n={len(vals):3d}  min={min(vals):6.3f}  "
            f"median={st.median(vals):6.3f}  max={max(vals):6.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--positives", nargs="*", default=[])
    ap.add_argument("--negatives", nargs="*", default=[])
    ap.add_argument("--metric", default="core_rmsd",
                    help="column to calibrate (default core_rmsd)")
    ap.add_argument("--lower-is-better", action="store_true", default=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.table), delimiter="\t"))
    by_acc = {}
    for r in rows:
        by_acc.setdefault(r["accession"], r)

    pos = [by_acc[a] for a in args.positives if a in by_acc]
    neg = [by_acc[a] for a in args.negatives if a in by_acc]
    missing = [a for a in args.positives + args.negatives if a not in by_acc]

    print(f"{len(rows)} rows; {len(pos)}/{len(args.positives)} positives, "
          f"{len(neg)}/{len(args.negatives)} negatives located")
    if missing:
        print(f"  NOT IN RESULTS (absence is uninformative -- hit caps, "
              f"clustered DBs): {', '.join(missing)}")

    m = args.metric
    print(f"\n{m}:")
    print(describe("positives", [_f(r[m]) for r in pos]))
    print(describe("negatives", [_f(r[m]) for r in neg]))
    print(describe("all hits", [_f(r[m]) for r in rows]))

    print("\nn_core (how much of the catalytic shell actually aligned):")
    print(describe("positives", [_f(r.get("n_core")) for r in pos]))
    print(describe("negatives", [_f(r.get("n_core")) for r in neg]))

    pv = [_f(r[m]) for r in pos if _f(r[m]) == _f(r[m])]
    nv = [_f(r[m]) for r in neg if _f(r[m]) == _f(r[m])]
    if pv and nv:
        worst_pos, best_neg = max(pv), min(nv)
        print()
        if worst_pos < best_neg:
            mid = (worst_pos + best_neg) / 2
            print(f"  SEPARABLE: worst positive {worst_pos:.3f} < best negative "
                  f"{best_neg:.3f}")
            print(f"  suggested cutoff: {m} <= {mid:.2f} "
                  f"(margin {best_neg - worst_pos:.2f})")
        else:
            print(f"  NOT SEPARABLE on {m} alone: worst positive "
                  f"{worst_pos:.3f} >= best negative {best_neg:.3f}")
            print("  -> add a second criterion (n_core, catalytic-shell pLDDT, "
                  "carboxylate spacing) rather than moving the cutoff")
    else:
        print("\n  need both positives and negatives present to suggest a cutoff")

    print("\nReport the distribution above alongside any cutoff you adopt. "
          "A cutoff without its calibration set is not a result.")


if __name__ == "__main__":
    main()
