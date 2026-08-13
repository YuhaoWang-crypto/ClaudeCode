"""Score Stack against this project's models on identical folds and knockdowns.

Stack emits one AnnData of predicted *cells* per condition, in raw counts over
its own 15,012-gene universe.  This project's metrics work on effects in log1p
CP10K over the shared 6,642-gene space.  The conversion is the same one applied
to the measured data, so nothing favours either side:

    mean over predicted cells -> CP10K -> log1p -> subtract the predicted
    non-targeting condition

Every model is then scored on the *same* held-out line, the *same* 100
knockdowns, and the *same* genes, so the numbers sit in one table honestly.

What this does and does not test is worth stating plainly.  ``Stack-Large-Aligned``
was post-trained on chemical and cytokine perturbations; here it is asked to do
in-context transfer for genetic knockdowns.  The paper claims in-context learning
over arbitrary conditions, so this is within its stated scope, but it is not the
modality it was aligned on, and a weak score should be read accordingly.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np

from . import metrics as M
from .benchmark import DETruth, de_truth, evaluate
from .data import load_all, shared_perturbations
from .model import (ContextTransferModel, ControlBaseline, GlobalMeanBaseline,
                    Hyper, NaiveTransferBaseline)

BENCH = Path("/home/user/stack_bench")
CONTROL = "non-targeting"
NAME = {"k562": "K562", "rpe1": "RPE1", "hepg2": "HepG2", "jurkat": "Jurkat"}


def _profile(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Mean log1p CP10K profile of one predicted condition, plus its symbols."""
    a = ad.read_h5ad(path)
    X = a.X
    counts = np.asarray(X.sum(axis=0)).ravel() if hasattr(X, "sum") else X.sum(0)
    total = counts.sum()
    if total <= 0:
        total = 1.0
    return np.log1p(counts / total * 1e4), np.asarray(a.var_names, dtype=str)


def stack_delta(target: str, perts: np.ndarray, symbols: np.ndarray
                ) -> tuple[np.ndarray, np.ndarray] | None:
    """Predicted effects for one fold, projected onto this project's gene axis.

    Also returns the mask of genes Stack actually covers.  Its 15,012-gene
    universe misses roughly a third of this project's 6,642, and scoring the
    remainder as zero effect would penalise Stack for a gene-space mismatch
    rather than for prediction quality -- so every model gets restricted to the
    same genes.
    """
    out = BENCH / f"pred_{target}"
    ctrl_path = out / f"{CONTROL}.h5ad"
    if not ctrl_path.exists():
        return None

    ctrl, sym = _profile(ctrl_path)
    # first occurrence wins, matching how the inputs were built
    lookup: dict[str, int] = {}
    for i, s in enumerate(sym):
        lookup.setdefault(s, i)
    take = np.array([lookup.get(s, -1) for s in symbols])
    have = take >= 0
    print(f"  {have.sum()}/{symbols.size} genes recovered from Stack output")

    delta = np.zeros((perts.size, symbols.size))
    missing = []
    for i, p in enumerate(perts):
        f = out / f"{p}.h5ad"
        if not f.exists():
            missing.append(p)
            continue
        prof, _ = _profile(f)
        d = prof - ctrl
        delta[i, have] = d[take[have]]
    if missing:
        print(f"  WARNING: {len(missing)} conditions absent from Stack output")
    return delta, have


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results/stack_compare.json"))
    args = ap.parse_args()

    lines, genes, symbols = load_all()
    by_name = {cl.name: cl for cl in lines}
    common = shared_perturbations(lines)
    hyper = {}
    ctx = Path("results/context.json")
    if ctx.exists():
        for name, h in json.loads(ctx.read_text())["hyper"].items():
            clean = {}
            for k, v in h.items():
                default = getattr(Hyper(), k)
                clean[k] = (bool(v) if isinstance(default, bool)
                            else int(float(v)) if isinstance(default, int)
                            else float(v))
            hyper[name] = Hyper(**clean)

    results: dict[str, dict[str, dict[str, float]]] = {}
    restricted: dict[str, dict[str, dict[str, float]]] = {}
    for tag, pretty in NAME.items():
        pert_file = BENCH / f"perts_{tag}.txt"
        if not pert_file.exists():
            continue
        perts = np.loadtxt(pert_file, dtype=str)
        perts = np.array([p for p in perts if p in set(common)])
        target = by_name[pretty]
        sources = [cl for cl in lines if cl.name != pretty]
        print(f"\n=== {pretty}: {perts.size} knockdowns ===")

        truth = de_truth(target, perts)
        fold: dict[str, dict[str, float]] = {}

        models = {
            "control (delta=0)": ControlBaseline(),
            "global mean [challenge baseline]": GlobalMeanBaseline(),
            "naive transfer": NaiveTransferBaseline(),
            "ContextTransfer (ours)": ContextTransferModel(
                hyper.get(pretty, Hyper())),
        }
        for name, model in models.items():
            model.fit(sources, target.mu, symbols)
            fold[name] = evaluate(model.predict(perts), target, perts, truth,
                                  symbols)

        got = stack_delta(tag, perts, symbols)
        if got is not None:
            d, have = got
            # Restrict every model to the genes Stack covers, so the comparison
            # measures prediction quality rather than gene-space overlap.
            sub = replace(target, ctrl=target.ctrl[:, have],
                          pert=target.pert[:, have], genes=target.genes[have],
                          symbols=target.symbols[have], _mu=None, _delta=None)
            sub_sources = [replace(s_, ctrl=s_.ctrl[:, have],
                                   pert=s_.pert[:, have], genes=s_.genes[have],
                                   symbols=s_.symbols[have], _mu=None,
                                   _delta=None) for s_ in sources]
            sub_sym = symbols[have]
            sub_truth = de_truth(sub, perts)
            fold_r: dict[str, dict[str, float]] = {}
            for name, model in models.items():
                model.fit(sub_sources, sub.mu, sub_sym)
                fold_r[name] = evaluate(model.predict(perts), sub, perts,
                                        sub_truth, sub_sym)
            fold_r["Stack (Arc, in-context)"] = evaluate(
                d[:, have], sub, perts, sub_truth, sub_sym)
            restricted[pretty] = fold_r
        results[pretty] = fold

        for label, table in (("all 6,642 genes", fold),
                             (f"restricted to genes Stack covers", 
                              restricted.get(pretty))):
            if not table:
                continue
            print(f"  -- {label} --")
            base = table["global mean [challenge baseline]"]
            for name, m in table.items():
                print(f"    {name:<34} PDS={m['discrimination_score_l1']:.3f}  "
                      f"ovl@100={m['overlap_at_100']:.3f}  MAE={m['mae']:.4f}  "
                      f"score={M.vcc_score(m, base)['avg_score']:.3f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"results": results,
                                    "restricted": restricted,
                                    "lines": [NAME[t] for t in NAME
                                              if NAME[t] in results]},
                                   indent=2, default=str))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
