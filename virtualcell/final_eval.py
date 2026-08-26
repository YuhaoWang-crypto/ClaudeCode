"""Score the configuration that actually ships, against everything it beat.

Every other benchmark in this project answers a design question and uses
whichever settings that question needed.  This one answers a different question
-- *how good is the thing being submitted* -- so it takes the shipping
configuration exactly as :func:`virtualcell.vcc2026.threshold_safe` assembles it
and changes nothing:

* everything except effect scale and the two priors, from nested
  leave-one-line-out cross-validation on the source lines;
* ``beta`` and ``gene_w`` from the magnitude frontier, chosen under the metric
  threshold the challenge enforces rather than the aggregate a tuner prefers;
* ``esm_mix`` from the two-axis benchmark, acting only where the source has no
  measurement at all.

Scored on the single-source genome-wide setting the submission uses -- one K562
screen, out-of-lineage, held-out line's controls only -- because a number from
the four-line atlas would describe a model that cannot answer the official
panel.

    python -m virtualcell.final_eval
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth, evaluate
from .gwps import _tuned_per_fold, frame
from .model import (ContextTransferModel, ControlBaseline, GlobalMeanBaseline,
                    NaiveTransferBaseline, SourceBank)
from .two_axis import load_embeddings
from .vcc2026 import threshold_safe

HEAD = ("discrimination_score_l1", "overlap_at_100", "mae")


def shipping_hyper(fold: str, verbose: bool = False):
    """Exactly what the submission uses, per fold, with nothing re-tuned here."""
    hp = _tuned_per_fold([fold])[fold]
    hp, prior_name = threshold_safe(hp, verbose=verbose)
    return replace(hp, esm_mix=0.25), prior_name


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-eval", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/final_eval.json"))
    args = ap.parse_args()

    source, targets, genes, symbols = frame()
    embed = load_embeddings(symbols)
    bank = SourceBank.build([source], np.ones(1))
    rng = np.random.default_rng(args.seed)

    from .gwps import build_prior

    results: dict[str, dict] = {}
    for target in targets:
        hp, prior_name = shipping_hyper(target.name, verbose=(not results))
        prior = (build_prior(prior_name, symbols, exclude=target.name)
                 if hp.gene_w > 0 else None)

        avail = np.array(sorted(set(target.names) & set(source.names)))
        sel = np.sort(rng.choice(avail, min(args.n_eval, avail.size),
                                 replace=False))
        truth = de_truth(target, sel)

        arms = {
            "control (delta=0)": ControlBaseline(),
            "global mean [challenge baseline]": GlobalMeanBaseline(),
            "naive transfer": NaiveTransferBaseline(),
            "ContextTransfer [shipped]": ContextTransferModel(hp),
        }
        fold = {}
        for name, model in arms.items():
            kw = ({"bank": bank, "gene_prior": prior, "gene_embed": embed}
                  if isinstance(model, ContextTransferModel) else {})
            model.fit([source], target.mu, symbols, **kw)
            fold[name] = evaluate(model.predict(sel), target, sel, truth,
                                  symbols)
        fold["_hyper"] = {"beta": hp.beta, "gene_w": hp.gene_w,
                          "esm_mix": hp.esm_mix, "prior": prior_name}
        results[target.name] = fold

        base = fold["global mean [challenge baseline]"]
        print(f"\n--- held out {target.name} ({sel.size} knockdowns) ---")
        for name, m in fold.items():
            if name.startswith("_"):
                continue
            print(f"  {name:<34}PDS={m['discrimination_score_l1']:.4f}  "
                  f"ovl@100={m['overlap_at_100']:.4f}  MAE={m['mae']:.4f}  "
                  f"score={M.vcc_score(m, base)['avg_score']:+.4f}  "
                  f"balanced={M.vcc_score(m, base, clip=False)['avg_score']:+.4f}",
                  flush=True)

    print("\n" + "=" * 78)
    print("mean over held-out lines -- the shipping configuration")
    names = [n for n in results[targets[0].name] if not n.startswith("_")]
    print(f"\n  {'model':<34}{'PDS':>8}{'ovl@100':>10}{'MAE':>9}"
          f"{'score':>9}{'balanced':>10}{'thresholds':>12}")
    print("  " + "-" * 82)
    for name in names:
        ms = [results[t][name] for t in results]
        bs = [results[t]["global mean [challenge baseline]"] for t in results]
        sc = np.mean([M.vcc_score(m, b)["avg_score"] for m, b in zip(ms, bs)])
        bal = np.mean([M.vcc_score(m, b, clip=False)["avg_score"]
                       for m, b in zip(ms, bs)])
        # does it clear the baseline on every metric, on every fold?
        ok = all(m[k] >= b[k] if k != "mae" else m[k] <= b[k]
                 for m, b in zip(ms, bs) for k in HEAD)
        print(f"  {name:<34}"
              f"{np.mean([m['discrimination_score_l1'] for m in ms]):>8.4f}"
              f"{np.mean([m['overlap_at_100'] for m in ms]):>10.4f}"
              f"{np.mean([m['mae'] for m in ms]):>9.4f}{sc:>9.4f}{bal:>10.4f}"
              f"{('all pass' if ok else 'FAILS'):>12}")

    print("\nthresholds = clears the challenge baseline on discrimination, DE "
          "overlap\nand error, on every fold. The challenge enforces a minimum "
          "on each, so a\nmodel that wins on aggregate while failing one is not "
          "admissible.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
