"""Post-hoc analysis: where does the model help, and where does it not?

Two questions the headline table cannot answer.

**Does it work on the knockdowns that matter?**  Roughly half of an essential-gene
panel produces almost no transcriptional response at this sequencing depth, and
averages over the whole panel are dominated by those.  Stratifying by measured
effect strength separates "predicts the strong responses" from "predicts that
most things do nothing".

**Which source line does the work?**  The model weights sources by control-profile
similarity to the target.  Printing those weights shows whether the weighting is
finding real structure or just spreading mass evenly.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth, evaluate
from .data import load_all, shared_perturbations
from .model import (ContextTransferModel, GlobalMeanBaseline, Hyper,
                    NaiveTransferBaseline)

STRATA = [(0, 10, "silent (<10 DE genes)"),
          (10, 100, "weak (10-100)"),
          (100, 500, "moderate (100-500)"),
          (500, 10 ** 9, "strong (>500)")]


def stratified(hp: Hyper, n_eval: int = 800, seed: int = 0) -> dict:
    """Score each model within bands of measured perturbation strength."""
    lines, _, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(seed)
    sel = np.sort(rng.choice(common, min(n_eval, common.size), replace=False))

    rows: dict[str, dict[str, list[float]]] = {}
    weights_report: dict[str, dict[str, float]] = {}

    for i, target in enumerate(lines):
        sources = [s for j, s in enumerate(lines) if j != i]
        truth = de_truth(target, sel)
        n_de = truth.sig.sum(axis=1)

        model = ContextTransferModel(hp).fit(sources, target.mu, symbols)
        weights_report[target.name] = dict(zip(
            [s.name for s in sources],
            np.round(model.context_weights(sources, target.mu), 3).tolist()))

        preds = {
            "global mean [challenge baseline]":
                GlobalMeanBaseline().fit(sources, target.mu, symbols).predict(sel),
            "naive transfer":
                NaiveTransferBaseline().fit(sources, target.mu, symbols).predict(sel),
            "ContextTransfer (ours)": model.predict(sel),
        }

        for lo, hi, label in STRATA:
            keep = np.flatnonzero((n_de >= lo) & (n_de < hi))
            if keep.size < 5:
                continue
            sub_truth = type(truth)(sig=truth.sig[keep], lfc=truth.lfc[keep],
                                    names=sel[keep])
            base_m = None
            for name, pred in preds.items():
                m = evaluate(pred[keep], target, sel[keep], sub_truth, symbols)
                if base_m is None:
                    base_m = m
                m["vcc_score"] = M.vcc_score(m, base_m)["avg_score"]
                m["n"] = int(keep.size)
                rows.setdefault(label, {}).setdefault(name, []).append(m)

    return {"strata": rows, "context_weights": weights_report}


def report(result: dict) -> str:
    out = ["\n### performance by measured perturbation strength\n"]
    head = (f"{'stratum':<24}{'model':<34}{'n':>6}{'PDS':>9}{'ovl@100':>9}"
            f"{'dir':>9}{'pearson':>9}{'score':>9}")
    out.append(head)
    out.append("-" * len(head))
    for lo, hi, label in STRATA:
        if label not in result["strata"]:
            continue
        for name, ms in result["strata"][label].items():
            def avg(k):
                return float(np.mean([m[k] for m in ms]))
            out.append(f"{label:<24}{name:<34}{int(np.sum([m['n'] for m in ms])):>6}"
                       f"{avg('discrimination_score_l1'):>9.3f}"
                       f"{avg('overlap_at_100'):>9.3f}"
                       f"{avg('de_direction_match'):>9.3f}"
                       f"{avg('pearson_delta'):>9.3f}"
                       f"{avg('vcc_score'):>9.3f}")
        out.append("")

    out.append("### source-line weights chosen from control profiles alone\n")
    for target, w in result["context_weights"].items():
        pretty = "  ".join(f"{k}={v:.2f}" for k, v in w.items())
        out.append(f"  target {target:<8} {pretty}")
    return "\n".join(out)


def main() -> None:
    from .run import _consensus_hyper
    hp = _consensus_hyper(Path("results/context.json"))
    print(f"stratified analysis with {hp}\n")
    result = stratified(hp)
    text = report(result)
    print(text)
    out = Path("results/stratified.txt")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    Path("results/stratified.json").write_text(json.dumps(result, indent=2,
                                                          default=str))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
