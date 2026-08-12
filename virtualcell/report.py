"""Render the benchmark results as markdown tables.

Hand-transcribing numbers out of JSON into a report is how reports end up
disagreeing with the data they describe.  This emits the tables straight from
``results/*.json`` so the report and the run cannot drift apart.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import metrics as M

RESULTS = Path("results")
HEADLINE = [("discrimination_score_l1", "discrimination"),
            ("overlap_at_100", "DE overlap@100"),
            ("mae", "MAE")]
EXTRA = [("de_direction_match", "DE direction"),
         ("de_spearman_lfc", "DE LFC ρ"),
         ("pearson_delta", "Pearson (effect)"),
         ("mae_delta", "MAE (effect)")]
SHORT = {"control (delta=0)": "control (Δ=0)",
         "global mean [challenge baseline]": "global mean *(baseline)*",
         "naive transfer": "naive transfer",
         "ContextTransfer (ours)": "**ContextTransfer**"}


def _fold_table(payload: dict, cols: list[tuple[str, str]]) -> str:
    results, lines = payload["results"], payload["lines"]
    base = "global mean [challenge baseline]"
    head = ["model"] + [f"{lab}" for _, lab in cols] + ["VCC score", "balanced"]
    out = ["| " + " | ".join(head) + " |",
           "|" + "|".join(["---"] * len(head)) + "|"]
    for model in SHORT:
        avg = {k: float(np.mean([results[c][model][k] for c in lines]))
               for k, _ in cols}
        sc = float(np.mean([M.vcc_score(results[c][model], results[c][base])
                            ["avg_score"] for c in lines]))
        bal = float(np.mean([M.vcc_score(results[c][model], results[c][base],
                                         clip=False)["avg_score"] for c in lines]))
        cells = [f"{avg[k]:.4f}" for k, _ in cols]
        out.append("| " + " | ".join([SHORT[model]] + cells
                                     + [f"{sc:.3f}", f"{bal:+.3f}"]) + " |")
    return "\n".join(out)


def _per_line_table(payload: dict, key: str) -> str:
    results, lines = payload["results"], payload["lines"]
    out = ["| model | " + " | ".join(lines) + " | mean |",
           "|" + "|".join(["---"] * (len(lines) + 2)) + "|"]
    for model in SHORT:
        vals = [results[c][model][key] for c in lines]
        cells = [f"{v:.3f}" for v in vals] + [f"{np.mean(vals):.3f}"]
        out.append("| " + " | ".join([SHORT[model]] + cells) + " |")
    return "\n".join(out)


def _ablation_table(ablation: dict) -> str:
    names = list(ablation)
    sizes = sorted(int(k) for k in ablation[names[0]])
    out = ["| held-out line | " + " | ".join(f"{k} source line{'s' if k > 1 else ''}"
                                             for k in sizes) + " |",
           "|" + "|".join(["---"] * (len(sizes) + 1)) + "|"]
    means = []
    for n in names:
        ys = [float(np.mean([r["vcc_score"] for r in ablation[n][str(k)]]))
              for k in sizes]
        means.append(ys)
        out.append("| " + " | ".join([n] + [f"{y:.3f}" for y in ys]) + " |")
    col = np.array(means).mean(axis=0)
    out.append("| **mean** | " + " | ".join(f"**{y:.3f}**" for y in col) + " |")
    return "\n".join(out)


def main() -> None:
    parts: list[str] = []
    for regime, title in [("context", "Context zero-shot (the challenge's setting)"),
                          ("double", "Double-blind (knockdown unseen everywhere too)")]:
        path = RESULTS / f"{regime}.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text())
        parts.append(f"\n## {title}\n")
        parts.append(f"*{payload['n_eval_perturbations']} knockdowns, "
                     f"mean over {len(payload['lines'])} held-out cell lines*\n")
        parts.append("### Headline metrics\n")
        parts.append(_fold_table(payload, HEADLINE))
        parts.append("\n### Supplementary metrics\n")
        parts.append(_fold_table(payload, EXTRA))
        for key, label in [("discrimination_score_l1", "Discrimination, per cell line"),
                           ("overlap_at_100", "DE overlap@100, per cell line"),
                           ("mae", "MAE, per cell line")]:
            parts.append(f"\n### {label}\n")
            parts.append(_per_line_table(payload, key))
        parts.append("\n### Hyperparameters chosen per fold\n")
        parts.append("```")
        for name, hp in payload["hyper"].items():
            keep = {k: v for k, v in hp.items() if k != "on_target"}
            parts.append(f"{name:8s} " + "  ".join(f"{k}={v}" for k, v in keep.items()))
        parts.append("```")

    abl = RESULTS / "ablation.json"
    if abl.exists():
        parts.append("\n## Context ablation\n")
        parts.append("*VCC score against the number of source cell lines, "
                     "hyperparameters held fixed*\n")
        parts.append(_ablation_table(json.loads(abl.read_text())))

    strat = RESULTS / "stratified.txt"
    if strat.exists():
        parts.append("\n## By measured perturbation strength\n")
        parts.append("```")
        parts.append(strat.read_text().strip())
        parts.append("```")

    text = "\n".join(parts)
    out = RESULTS / "tables.md"
    out.write_text(text)
    print(text)
    print(f"\n\nwrote {out}")


if __name__ == "__main__":
    main()
