"""Would acquiring more source cell lines actually help?

This project has said, repeatedly, that the single biggest lever on its weakest
result is a second source context covering the official panel.  That claim was
resting on the context ablation -- score against number of source lines -- which
shows accuracy improving as lines are added.  An improving curve is not the same
as a curve worth extrapolating, and the honest version needs two things the
ablation did not have.

**A saturating fit.**  Averaging across contexts cancels each source's
context-specific component, and the amount left to cancel shrinks like 1/k.  So
the consensus arm has an asymptote, and the question is not "does it improve"
but "how much is left above where we already are".

**A slope arm.**  The consensus asymptote bounds *averaging*.  It does not bound
an architecture that fits a trend across contexts rather than averaging them --
and that architecture is untestable at k=1 and barely testable at k=3, so a
saturating consensus curve alone cannot conclude that more lines are worthless.

Both arms come from ``transferability-prior-eval``'s kernel, which was written
from measurements on this same class of data; this module only supplies this
project's four lines and prints the result.

    python -m virtualcell.context_curve
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from .data import load_all, shared_perturbations

KERNEL = Path("/root/.claude/skills/synced/transferability-prior-eval")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-perturbations", type=int, default=600)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/context_curve.json"))
    args = ap.parse_args()

    sys.path.insert(0, str(KERNEL))
    from kernel import source_context_curve                      # noqa: E402

    lines, _, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(args.seed)
    perts = list(np.sort(rng.choice(common, min(args.n_perturbations,
                                                len(common)), replace=False)))

    delta = {cl.name: cl.delta for cl in lines}
    index = {cl.name: {n: i for i, n in enumerate(cl.names)} for cl in lines}
    basal = {cl.name: cl.mu for cl in lines}

    out = {}
    for target in lines:
        sources = [cl.name for cl in lines if cl.name != target.name]
        res = source_context_curve(delta, index, sources, target.name, perts,
                                   basal_by_context=basal)
        out[target.name] = res
        print(f"\n--- held out {target.name} "
              f"(sources: {', '.join(sources)}) ---")
        print(f"  {'k sources':>10}{'consensus':>12}{'slope model':>14}")
        for row in res["curve"]:
            s = row.get("slope_model")
            print(f"  {row['k_sources']:>10}{row['consensus']:>12.4f}"
                  + (f"{s:>14.4f}" if s is not None else f"{'-':>14}"))
        if res["fit"]:
            f = res["fit"]
            print(f"  1/k asymptote {f['asymptote']:.4f}; extrapolated "
                  + ", ".join(f"k={k}: {v:.4f}"
                              for k, v in f["extrapolated"].items()))

    # What the two arms say, averaged over folds.
    ks = sorted({r["k_sources"] for v in out.values() for r in v["curve"]})
    print(f"\n{'':<14}" + "".join(f"{'k=' + str(k):>12}" for k in ks))
    for arm in ("consensus", "slope_model"):
        vals = []
        for k in ks:
            per = [r[arm] for v in out.values() for r in v["curve"]
                   if r["k_sources"] == k and arm in r]
            vals.append(np.mean(per) if per else np.nan)
        print(f"{arm:<14}" + "".join(f"{v:>12.4f}" for v in vals))
    asym = np.mean([v["fit"]["asymptote"] for v in out.values() if v["fit"]])
    at3 = np.mean([r["consensus"] for v in out.values() for r in v["curve"]
                   if r["k_sources"] == 3])
    print(f"\nconsensus asymptote (mean over folds): {asym:.4f}")
    print(f"at the 3 sources this project has:      {at3:.4f}")
    print(f"so averaging has at most                {asym - at3:+.4f} left, "
          f"on a scale where perfect shape is 1.0")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
