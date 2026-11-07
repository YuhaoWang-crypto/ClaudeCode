"""
Template: score a batch of designed glycoproteins and pick what to express.

Copy, point ``DESIGN_DIR`` at a folder of predicted structures (each with its
modelled glycan), and run. Implements stage 4 of
``reference/design-pipeline.md``.

    python3 score_design.py /path/to/predictions --top-frac 0.05

The one rule that makes the output meaningful: **every structure must carry
the same modelled glycan**. Mean per-atom DeltaRASA is size-confounded
(Spearman -0.38 against atom count on real structures), so a campaign mixing
Man6 and Man9 models ranks tree size, not shielding. The script refuses to
rank a mixed batch unless you pass --allow-mixed.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

# Make the repo's glycan_pipeline importable wherever this file is copied to:
# walk up from here (and from the cwd) looking for the package.
def _add_pipeline_to_path() -> None:
    seen = set()
    for start in (os.path.dirname(os.path.abspath(__file__)), os.getcwd()):
        d = start
        while d and d not in seen:
            seen.add(d)
            if os.path.isdir(os.path.join(d, "glycan_pipeline")):
                if d not in sys.path:
                    sys.path.insert(0, d)
                return
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent


_add_pipeline_to_path()

from glycan_pipeline.g1_delta_rasa import site_delta_rasa          # noqa: E402
from glycan_pipeline.glycans import (find_glycan_trees,            # noqa: E402
                                     load_model, sequons_in_structure)

STRUCTURE_EXTS = (".pdb", ".cif", ".mmcif")


def score_one(path: str, context: str = "protein") -> list[dict]:
    """Every glycosylation site in one structure, with its DeltaRASA."""
    model = load_model(path)
    trees = [t for t in find_glycan_trees(model) if t.asn_resseq is not None]
    if not trees:
        return []

    sequons = sequons_in_structure(model)
    rows = site_delta_rasa(model, trees, context=context)
    name = os.path.splitext(os.path.basename(path))[0]
    for row in rows:
        key = (row["chain"], row["resseq"])
        motif = sequons.get(key)
        row["design"] = name
        row["sequon"] = motif or "NONE"
        row["is_nxt"] = bool(motif and motif.endswith("T"))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design_dir", help="folder of predicted structures")
    ap.add_argument("--top-frac", type=float, default=0.05,
                    help="fraction to shortlist, by either statistic (default 0.05)")
    ap.add_argument("--context", default="protein",
                    choices=["protein", "protein+other_glycans"])
    ap.add_argument("--out", default="delta_rasa_scores.csv")
    ap.add_argument("--allow-mixed", action="store_true",
                    help="rank even if designs carry different glycan sizes")
    args = ap.parse_args()

    paths = sorted(
        os.path.join(args.design_dir, f)
        for f in os.listdir(args.design_dir)
        if f.lower().endswith(STRUCTURE_EXTS)
    )
    if not paths:
        print(f"no structures in {args.design_dir}")
        return 1

    rows = []
    for i, path in enumerate(paths, 1):
        try:
            rows.extend(score_one(path, args.context))
        except Exception as exc:
            print(f"  [{i}/{len(paths)}] {os.path.basename(path)}: "
                  f"FAILED ({type(exc).__name__}: {exc})")
            continue
        print(f"  [{i}/{len(paths)}] {os.path.basename(path)}")

    if not rows:
        print("no glycosylated sites found - is the glycan modelled?")
        return 1

    # Guard against the size confound before ranking anything.
    sizes = sorted({r["n_atoms"] for r in rows})
    if len(sizes) > 1 and not args.allow_mixed:
        print("\nREFUSING TO RANK: glycan sizes differ across sites "
              f"(atom counts {sizes[:6]}{'...' if len(sizes) > 6 else ''}).")
        print("Mean per-atom DeltaRASA is size-confounded. Model the SAME")
        print("glycan on every design, or pass --allow-mixed knowingly.")
        print(f"Per-site scores still written to {args.out}.")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    if len(sizes) > 1 and not args.allow_mixed:
        return 2

    # Shortlist on EITHER statistic, as the paper does.
    k = max(1, int(round(len(rows) * args.top_frac)))
    by_mean = sorted(rows, key=lambda r: -r["mean_delta_rasa"])[:k]
    by_count = sorted(rows, key=lambda r: -r["n_contact_atoms"])[:k]
    shortlist = {(r["design"], r["site"]): r for r in by_mean + by_count}

    print(f"\n{len(rows)} sites scored; shortlist = top {args.top_frac:.0%} "
          f"by either statistic -> {len(shortlist)} sites")
    print(f"\n{'design':>22} {'site':>10} {'sequon':>7} "
          f"{'meandRASA':>10} {'contact':>9}")
    for r in sorted(shortlist.values(), key=lambda r: -r["mean_delta_rasa"]):
        print(f"{r['design'][:22]:>22} {r['site']:>10} {r['sequon']:>7} "
              f"{r['mean_delta_rasa']:>9.2f}% "
              f"{r['n_contact_atoms']:>3}/{r['n_atoms']:<5}")

    nxs = [r for r in shortlist.values() if not r["is_nxt"]]
    if nxs:
        print(f"\nnote: {len(nxs)} shortlisted site(s) use NXS, which is "
              "occupied less often than NXT.")
    print("\nNow inspect these by eye: the glycan should sit in a pocket or")
    print("crevice, not merely next to a bulge. Then express, and go to")
    print("reference/wet-validation.md.")
    print(f"\nfull table: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
