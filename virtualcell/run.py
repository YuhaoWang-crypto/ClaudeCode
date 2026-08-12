"""Command-line entry point: run the benchmark and write results to JSON."""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import METRIC_ORDER, context_ablation, run
from .model import Hyper


def summarise(payload: dict) -> str:
    """Render the per-cell-line tables and the across-line average."""
    results = payload["results"]
    lines = payload["lines"]
    models = list(next(iter(results.values())).keys())
    baseline = "global mean [challenge baseline]"
    out: list[str] = []

    for name in lines:
        fold = results[name]
        out.append(f"\n### held-out cell line: {name}\n")
        head = f"{'model':<34}" + "".join(f"{m[:13]:>15}" for m in METRIC_ORDER)
        out.append(head + f"{'VCC score':>12}")
        out.append("-" * len(head + f"{'VCC score':>12}"))
        for model in models:
            m = fold[model]
            score = M.vcc_score(m, fold[baseline])["avg_score"]
            row = f"{model:<34}" + "".join(f"{m[k]:>15.4f}" for k in METRIC_ORDER)
            out.append(row + f"{score:>12.4f}")

    out.append("\n### mean across all four held-out cell lines\n")
    head = f"{'model':<34}" + "".join(f"{m[:13]:>15}" for m in METRIC_ORDER)
    out.append(head + f"{'VCC score':>12}")
    out.append("-" * len(head + f"{'VCC score':>12}"))
    for model in models:
        avg = {k: float(np.mean([results[c][model][k] for c in lines]))
               for k in METRIC_ORDER}
        scores = [M.vcc_score(results[c][model], results[c][baseline])["avg_score"]
                  for c in lines]
        row = f"{model:<34}" + "".join(f"{avg[k]:>15.4f}" for k in METRIC_ORDER)
        out.append(row + f"{np.mean(scores):>12.4f}")
    return "\n".join(out)


def _consensus_hyper(path: Path) -> Hyper:
    """Take the most commonly chosen value per switch across the folds.

    The ablation must hold the model fixed while the data varies, so it needs
    one setting rather than four.  Using the modal tuned value keeps it
    representative without re-tuning inside the ablation.
    """
    if not path.exists():
        return Hyper()
    chosen = json.loads(path.read_text())["hyper"]
    fields = {f.name for f in dataclasses.fields(Hyper)}
    out = {}
    for name in fields:
        values = [h[name] for h in chosen.values() if name in h]
        if not values:
            continue
        out[name] = max(set(map(str, values)), key=list(map(str, values)).count)
    coerced = {}
    for name, raw in out.items():
        default = getattr(Hyper(), name)
        if isinstance(default, bool):
            coerced[name] = raw == "True"
        elif isinstance(default, int) and not isinstance(default, bool):
            coerced[name] = int(float(raw))
        else:
            coerced[name] = float(raw)
    return Hyper(**coerced)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regime", default="context", choices=["context", "double"])
    ap.add_argument("--n-eval", type=int, default=None,
                    help="evaluate on a random subset of the shared perturbations")
    ap.add_argument("--tune-eval", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--ablation", action="store_true",
                    help="run the source-context count ablation instead")
    args = ap.parse_args()

    if args.ablation:
        hp = _consensus_hyper(Path("results/context.json"))
        print(f"context ablation with fixed hyperparameters: {hp}\n")
        result = context_ablation(n_eval=args.n_eval or 600, hp=hp, seed=args.seed)
        out = args.out or Path("results/ablation.json")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2, default=str))
        print(f"\nwrote {out}")
        return

    payload = run(regime=args.regime, n_eval=args.n_eval,
                  tune_eval=args.tune_eval, seed=args.seed)
    text = summarise(payload)
    print(text)

    out = args.out or Path(f"results/{args.regime}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "regime": args.regime,
        "lines": payload["lines"],
        "n_eval_perturbations": int(payload["eval_perts"].size),
        "eval_perturbations": payload["eval_perts"].tolist(),
        "results": payload["results"],
        "hyper": {k: dataclasses.asdict(v) for k, v in payload["hyper"].items()},
    }, indent=2, default=str))
    out.with_suffix(".txt").write_text(text)
    print(f"\nwrote {out} and {out.with_suffix('.txt')}")


if __name__ == "__main__":
    main()
