"""Leave-one-cell-line-out benchmark for zero-shot perturbation prediction.

Each fold hides one cell line completely: the model sees that line's
non-targeting controls and nothing else from it, trains on the remaining lines,
and predicts knockdown responses that are then scored against the held-out
measurements.  That is the Virtual Cell Challenge 2 setting reproduced on public
data, with four contexts instead of one.

Two regimes are run:

``context``
    The knocked-down gene was measured in the source lines, just never in the
    target.  This is the challenge's own setting.

``double``
    The gene is additionally deleted from every source line, so neither the
    context nor the perturbation has ever been seen.  Nothing in the training
    data names this gene; the model has to reach it through perturbations whose
    effects resemble it.

Hyperparameters are chosen by a nested leave-one-out over the *source* lines
only, so the held-out line never influences model selection.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, replace

import numpy as np

from . import metrics as M
from .data import CellLine, load_all, shared_perturbations
from .model import (ContextTransferModel, ControlBaseline, GlobalMeanBaseline,
                    Hyper, NaiveTransferBaseline, SourceBank)

METRIC_ORDER = ["discrimination_score_l1", "overlap_at_100", "overlap_at_50",
                "de_direction_match", "de_spearman_lfc", "pearson_delta",
                "mae", "mae_delta"]


@dataclass
class DETruth:
    """Measured differential expression for one cell line's perturbations."""

    sig: np.ndarray      # (n_pert, G) bool, FDR < 0.05
    lfc: np.ndarray      # (n_pert, G) log2 fold change
    names: np.ndarray


def de_truth(line: CellLine, perts: np.ndarray) -> DETruth:
    """Call DE for the given perturbations using the control-replicate null."""
    batch, samp = M.control_variance_model(line.ctrl, line.ctrl_ncells)
    index = {n: i for i, n in enumerate(line.names)}
    rows = np.array([index[p] for p in perts])
    sig, lfc = M.de_significance(line.delta[rows], line.ncells[rows],
                                 batch, samp, line.mu)
    return DETruth(sig=sig, lfc=lfc, names=perts)


def evaluate(pred_delta: np.ndarray, target: CellLine, perts: np.ndarray,
             truth: DETruth, symbols: np.ndarray) -> dict[str, float]:
    """Score a prediction with the Virtual Cell Challenge metric suite."""
    index = {n: i for i, n in enumerate(target.names)}
    rows = np.array([index[p] for p in perts])
    real_delta = target.delta[rows]
    real_expr = target.pert[rows]
    pred_expr = target.mu[None, :] + pred_delta
    pred_lfc = pred_delta / M.LOG2

    return {
        "discrimination_score_l1": float(np.mean(
            M.discrimination_score(pred_delta, real_delta, perts, symbols))),
        "overlap_at_100": float(np.mean(
            M.overlap_at_k(pred_lfc, truth.lfc, truth.sig, k=100))),
        "overlap_at_50": float(np.mean(
            M.overlap_at_k(pred_lfc, truth.lfc, truth.sig, k=50))),
        "de_direction_match": float(np.nanmean(
            M.de_direction_match(pred_lfc, truth.lfc, truth.sig))),
        "de_spearman_lfc": float(np.nanmean(
            M.de_spearman_lfc(pred_lfc, truth.lfc, truth.sig))),
        "pearson_delta": float(np.mean(M.pearson_delta(pred_delta, real_delta))),
        "mae": float(np.mean(M.mae(pred_expr, real_expr))),
        "mae_delta": float(np.mean(M.mae_delta(pred_delta, real_delta))),
    }


def drop_perturbations(line: CellLine, drop: set[str]) -> CellLine:
    """Return a copy of a cell line with the named knockdowns removed."""
    keep = np.array([n not in drop for n in line.names])
    # _delta is cached against the old row set, so it must be invalidated; _mu
    # depends only on the controls and survives.
    return replace(line, pert=line.pert[keep], names=line.names[keep],
                   ncells=line.ncells[keep], _delta=None)


# --------------------------------------------------------------------------
# hyperparameter selection, source lines only
# --------------------------------------------------------------------------

SEARCH: dict[str, list] = {
    "gamma":       [0.0, 0.25, 0.5, 0.75, 1.0],
    "beta":        [0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0],
    "shrink":      [0.0, 0.5, 1.0, 2.0],
    "smooth":      [0.0, 0.15, 0.3, 0.5, 0.7],
    "use_global":  [0.0, 0.1, 0.2, 0.35],
    "temp":        [0.1, 0.25, 0.6, np.inf],
    "n_neighbors": [5, 15, 30],
    "mod_clip":    [2.0, 3.0, 5.0],
    "rank_mix":    [0.0, 0.25, 0.5, 0.75, 1.0],
    "rank":        [20, 40, 80, 160, 320],
    "unseen_k":    [5, 15, 25, 50],
}


def _inner_score(hp: Hyper, sources: list[CellLine], symbols: np.ndarray,
                 n_eval: int, rng: np.random.Generator,
                 cache: dict, regime: str = "context") -> float:
    """Mean leaderboard score across held-out *source* lines."""
    if len(sources) < 2:
        return 0.0
    scores = []
    for i, held in enumerate(sources):
        rest = [s for j, s in enumerate(sources) if j != i]

        key = ("fold", i)
        if key not in cache:
            avail = sorted(set(held.names) & {n for s in rest for n in s.names})
            sel = np.array(avail)
            if sel.size > n_eval:
                sel = rng.choice(sel, n_eval, replace=False)
                sel.sort()
            if regime == "double":
                # tune under the same handicap the model will face: the inner
                # sources must not have seen these knockdowns either
                rest = [drop_perturbations(s, set(sel)) for s in rest]
            truth = de_truth(held, sel)
            # the reference model does not depend on hyperparameters
            base = GlobalMeanBaseline().fit(rest, held.mu, symbols)
            base_m = evaluate(base.predict(sel), held, sel, truth, symbols)
            cache[key] = (sel, truth, base_m, rest)
        sel, truth, base_m, rest = cache[key]

        # The bank depends on the source set and on the context weights, which
        # only ``temp`` moves; caching on that pair keeps the search cheap
        # without ever approximating the model being scored.
        bkey = ("bank", i, float(hp.temp))
        if bkey not in cache:
            w = ContextTransferModel(hp).context_weights(rest, held.mu)
            cache[bkey] = SourceBank.build(rest, w)

        model = ContextTransferModel(hp).fit(rest, held.mu, symbols,
                                             bank=cache[bkey])
        user_m = evaluate(model.predict(sel), held, sel, truth, symbols)
        scores.append(M.vcc_score(user_m, base_m)["avg_score"])
    return float(np.mean(scores))


def tune(sources: list[CellLine], symbols: np.ndarray, n_eval: int = 400,
         passes: int = 2, seed: int = 0, verbose: bool = True,
         regime: str = "context") -> Hyper:
    """Coordinate ascent over the switches, scored on source lines only."""
    rng = np.random.default_rng(seed)
    cache: dict = {}
    hp = Hyper()
    best = _inner_score(hp, sources, symbols, n_eval, rng, cache, regime)
    if verbose:
        print(f"      start avg_score={best:.4f}")

    for p in range(passes):
        improved = False
        for key, values in SEARCH.items():
            current = getattr(hp, key)
            for v in values:
                if v == current:
                    continue
                cand = replace(hp, **{key: v})
                score = _inner_score(cand, sources, symbols, n_eval, rng, cache,
                                     regime)
                if score > best + 1e-6:
                    best, hp, improved = score, cand, True
            if verbose:
                print(f"      pass{p} {key:<12} -> {getattr(hp, key)!r:>8} "
                      f"(avg_score={best:.4f})", flush=True)
        if not improved:
            break
    return hp


# --------------------------------------------------------------------------
# the benchmark itself
# --------------------------------------------------------------------------

def run(regime: str = "context", n_eval: int | None = None, tune_eval: int = 400,
        seed: int = 0, verbose: bool = True) -> dict:
    lines, genes, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(seed)

    if n_eval is not None and common.size > n_eval:
        eval_perts = np.sort(rng.choice(common, n_eval, replace=False))
    else:
        eval_perts = common

    if verbose:
        print(f"\n{'=' * 78}\nregime={regime}  "
              f"{len(lines)} cell lines, {genes.size} genes, "
              f"{eval_perts.size} evaluated perturbations\n{'=' * 78}")

    results: dict[str, dict[str, dict[str, float]]] = {}
    chosen: dict[str, Hyper] = {}

    for i, target in enumerate(lines):
        t0 = time.time()
        sources = [s for j, s in enumerate(lines) if j != i]
        if verbose:
            print(f"\n--- held-out cell line: {target.name} "
                  f"(sources: {', '.join(s.name for s in sources)}) ---")

        if regime == "double":
            # the perturbation is unknown everywhere, not just in the target
            sources = [drop_perturbations(s, set(eval_perts)) for s in sources]

        truth = de_truth(target, eval_perts)
        if verbose:
            print(f"    measured DE genes/perturbation: "
                  f"median {np.median(truth.sig.sum(1)):.0f}")

        if verbose:
            print("    tuning on source lines only...")
        hp = tune(sources, symbols, n_eval=tune_eval, seed=seed, verbose=verbose,
                  regime=regime)
        chosen[target.name] = hp

        models = {
            "control (delta=0)": ControlBaseline(),
            "global mean [challenge baseline]": GlobalMeanBaseline(),
            "naive transfer": NaiveTransferBaseline(),
            "ContextTransfer (ours)": ContextTransferModel(hp),
        }
        fold: dict[str, dict[str, float]] = {}
        for name, model in models.items():
            model.fit(sources, target.mu, symbols)
            fold[name] = evaluate(model.predict(eval_perts), target,
                                  eval_perts, truth, symbols)
        results[target.name] = fold
        if verbose:
            base = fold["global mean [challenge baseline]"]
            for name, m in fold.items():
                s = M.vcc_score(m, base)["avg_score"]
                print(f"    {name:<34} PDS={m['discrimination_score_l1']:.3f}  "
                      f"ovl@100={m['overlap_at_100']:.3f}  "
                      f"MAE={m['mae']:.4f}  score={s:.3f}")
            print(f"    ({time.time() - t0:.0f}s)")

    return {"results": results, "hyper": chosen, "eval_perts": eval_perts,
            "lines": [cl.name for cl in lines]}


def context_ablation(n_eval: int = 600, hp: Hyper | None = None,
                     seed: int = 0, verbose: bool = True) -> dict:
    """Does adding *contexts* help zero-shot transfer, or just adding cells?

    For every held-out line, fit on each possible subset of the remaining lines
    and score it.  If the claim that virtual cells are limited by contextual
    coverage rather than by scale is right, one source line should be clearly
    worse than three even though the perturbation panel is nearly identical in
    all of them -- the extra lines add no new knockdowns, only new contexts.

    Hyperparameters are held fixed across subset sizes so the comparison
    measures the data, not the tuning.
    """
    from itertools import combinations

    lines, _, symbols = load_all()
    common = shared_perturbations(lines)
    rng = np.random.default_rng(seed)
    sel = np.sort(rng.choice(common, min(n_eval, common.size), replace=False))
    hp = hp or Hyper()

    out: dict[str, dict[int, list[dict]]] = {}
    for i, target in enumerate(lines):
        pool = [s for j, s in enumerate(lines) if j != i]
        truth = de_truth(target, sel)
        base = GlobalMeanBaseline().fit(pool, target.mu, symbols)
        base_m = evaluate(base.predict(sel), target, sel, truth, symbols)

        per_size: dict[int, list[dict]] = {}
        for k in range(1, len(pool) + 1):
            rows = []
            for subset in combinations(pool, k):
                model = ContextTransferModel(hp).fit(list(subset), target.mu, symbols)
                m = evaluate(model.predict(sel), target, sel, truth, symbols)
                m["sources"] = ",".join(s.name for s in subset)
                m["vcc_score"] = M.vcc_score(m, base_m)["avg_score"]
                rows.append(m)
            per_size[k] = rows
            if verbose:
                print(f"  {target.name:7s} n_sources={k}: "
                      f"PDS={np.mean([r['discrimination_score_l1'] for r in rows]):.4f} "
                      f"ovl@100={np.mean([r['overlap_at_100'] for r in rows]):.4f} "
                      f"pearson={np.mean([r['pearson_delta'] for r in rows]):.4f} "
                      f"score={np.mean([r['vcc_score'] for r in rows]):.4f}", flush=True)
        out[target.name] = per_size
    return out
