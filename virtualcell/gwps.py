"""The genome-wide K562 screen as a source atlas, and the benchmark it needs.

Everything else in this project transfers between four cell lines whose
knockdown panels are the *same* ~2,000 common-essential genes.  The Virtual Cell
Challenge 2026 validation panel is not drawn from that set at all: of its 300
target genes, **zero** appear in the four-line atlas.  They are non-essential
regulators -- ACLY, AGO1, AKT2, ADNP, SOCS5 -- which an essential-gene screen
deliberately excludes.

Replogle's genome-wide K562 arm does cover them: 9,866 knockdowns including
272 of the 300.  So a model that is to say anything at all about the official
panel has to be built on that screen instead.  The cost is that it is a single
cell line, which removes the cross-context averaging the four-line benchmark
measured.  Whether what remains -- transfer plus context modulation plus
on-target knockdown -- still beats the challenge baseline is an empirical
question this module answers rather than assumes:

    python -m virtualcell.gwps --benchmark

holds out RPE1, HepG2 and Jurkat one at a time, gives the model *only* the
genome-wide K562 screen and the held-out line's controls, tunes on the other two
held-out lines, and scores with the challenge metrics.  K562 itself is never a
target, because the source screen is the same cell line and scoring it would be
reading back the training data.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth, evaluate
from .data import DATA, CellLine, load_all, load_replogle
from .model import (ContextTransferModel, ControlBaseline, GlobalMeanBaseline,
                    Hyper, NaiveTransferBaseline, SourceBank)

GWPS_H5AD = DATA / "replogle" / "K562_gwps_raw_bulk.h5ad"
CACHE = DATA / "pseudobulk" / "k562_gwps.npz"


def load_gwps(cache: bool = True) -> CellLine:
    """Replogle's genome-wide K562 CRISPRi screen as a :class:`CellLine`.

    Reading and pseudobulking the h5ad takes about a minute, so the result is
    cached as a ``.npz`` next to the other prepared pseudobulk.
    """
    if cache and CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        return CellLine(name="K562-gw", ctrl=z["ctrl"], pert=z["pert"],
                        names=z["names"], ncells=z["ncells"],
                        ctrl_ncells=z["ctrl_ncells"], genes=z["genes"],
                        symbols=z["symbols"])

    line = load_replogle(GWPS_H5AD, "K562-gw")
    if cache:
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(CACHE, ctrl=line.ctrl, pert=line.pert,
                            names=line.names, ncells=line.ncells,
                            ctrl_ncells=line.ctrl_ncells, genes=line.genes,
                            symbols=line.symbols)
    return line


def project(line: CellLine, genes: np.ndarray) -> CellLine:
    """Restrict a cell line to a gene axis it already contains, in that order."""
    idx = {g: i for i, g in enumerate(line.genes)}
    take = np.array([idx[g] for g in genes])
    return replace(line, ctrl=line.ctrl[:, take], pert=line.pert[:, take],
                   genes=genes, symbols=line.symbols[take],
                   _mu=None, _delta=None)


def frame(targets: list[str] | None = None
          ) -> tuple[CellLine, list[CellLine], np.ndarray, np.ndarray]:
    """Genome-wide source and held-out targets on one shared gene axis.

    Returns ``(source, targets, genes, symbols)``.
    """
    lines, _, _ = load_all()
    keep = targets or ["RPE1", "HepG2", "Jurkat"]
    lines = [cl for cl in lines if cl.name in keep]

    source = load_gwps()
    shared = set(source.genes)
    for cl in lines:
        shared &= set(cl.genes)
    genes = np.array(sorted(shared))

    source = project(source, genes)
    return source, [project(cl, genes) for cl in lines], genes, source.symbols


# --------------------------------------------------------------------------
# tuning, with the held-out line excluded from selection
# --------------------------------------------------------------------------

# ``temp`` is dropped: with one source line the context-similarity softmax has
# nothing to choose between and returns weight 1 whatever its value.  So is
# ``rank``/``rank_mix`` order -- kept, but searched after the switches that
# matter more, since coordinate ascent takes the first improvement it finds.
SEARCH: dict[str, list] = {
    "beta":        [0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0],
    "gamma":       [0.0, 0.25, 0.5, 0.75, 1.0],
    "shrink":      [0.0, 0.5, 1.0, 2.0],
    "smooth":      [0.0, 0.15, 0.3, 0.5],
    "use_global":  [0.0, 0.1, 0.2, 0.35],
    "n_neighbors": [5, 15, 30],
    "mod_clip":    [2.0, 3.0, 5.0],
    "rank_mix":    [0.0, 0.25, 0.5, 0.75],
    "rank":        [40, 80, 160, 320],
}
# ``unseen_k`` is absent on purpose: every knockdown scored here is present in
# the genome-wide screen, so the functional-neighbour path is never taken and
# searching it would spend a fifth of the budget measuring nothing.  It still
# matters for the 26 official panel targets the screen misses, which is why the
# submission keeps the four-line benchmark's value for it.


def _score_on(hp: Hyper, source: CellLine, targets: list[CellLine],
              symbols: np.ndarray, cache: dict) -> float:
    """Mean unclipped leaderboard score of ``hp`` over the given target lines."""
    fitted = None
    out = []
    for t in targets:
        sel, truth, base_m = cache[t.name]
        # The consensus is target-independent, so it is built for the first
        # target and re-pointed at the rest.  That is the whole cost of a fit.
        fitted = (fitted.retarget([source], t.mu) if fitted is not None else
                  ContextTransferModel(hp).fit([source], t.mu, symbols,
                                               bank=cache["bank"]))
        m = evaluate(fitted.predict(sel), t, sel, truth, symbols)
        out.append(M.vcc_score(m, base_m, clip=False)["avg_score"])
    return float(np.mean(out))


def _prepare(source: CellLine, targets: list[CellLine], symbols: np.ndarray,
             n_eval: int, seed: int, cache: dict) -> None:
    """Pick the evaluation knockdowns and score the reference model, once.

    With a single source line the context weights are fixed at 1, so the
    :class:`SourceBank` is fully hyperparameter-independent and is built here
    once for the whole search rather than per candidate.
    """
    rng = np.random.default_rng(seed)
    have = set(source.names)
    if "bank" not in cache:
        cache["bank"] = SourceBank.build([source], np.ones(1))
    for t in targets:
        if t.name in cache:
            continue
        avail = np.array(sorted(set(t.names) & have))
        sel = np.sort(rng.choice(avail, n_eval, replace=False)) \
            if avail.size > n_eval else avail
        truth = de_truth(t, sel)
        base = GlobalMeanBaseline().fit([source], t.mu, symbols)
        cache[t.name] = (sel, truth,
                         evaluate(base.predict(sel), t, sel, truth, symbols))


def tune(source: CellLine, targets: list[CellLine], symbols: np.ndarray,
         n_eval: int = 300, passes: int = 2, seed: int = 0,
         verbose: bool = True) -> Hyper:
    """Coordinate ascent over the switches, scored on ``targets`` only."""
    cache: dict = {}
    _prepare(source, targets, symbols, n_eval, seed, cache)
    hp = Hyper(temp=np.inf)
    best = _score_on(hp, source, targets, symbols, cache)
    if verbose:
        print(f"      start avg_score={best:+.4f}", flush=True)

    for p in range(passes):
        improved = False
        for key, values in SEARCH.items():
            for v in values:
                if v == getattr(hp, key):
                    continue
                cand = replace(hp, **{key: v})
                s = _score_on(cand, source, targets, symbols, cache)
                if s > best + 1e-6:
                    best, hp, improved = s, cand, True
            if verbose:
                print(f"      pass{p} {key:<12} -> {getattr(hp, key)!r:>7} "
                      f"(avg_score={best:+.4f})", flush=True)
        if not improved:
            break
    return hp


# --------------------------------------------------------------------------
# the benchmark
# --------------------------------------------------------------------------

def run(n_eval: int = 400, tune_eval: int = 300, seed: int = 0,
        verbose: bool = True) -> dict:
    """Hold out each non-K562 line; the only source is the genome-wide screen."""
    source, targets, genes, symbols = frame()
    if verbose:
        print(f"\n{'=' * 78}\nsingle-source transfer from the genome-wide K562 "
              f"screen\n{source.names.size:,} knockdowns, {genes.size:,} shared "
              f"genes, {len(targets)} held-out lines\n{'=' * 78}")

    rng = np.random.default_rng(seed)
    results: dict[str, dict] = {}
    chosen: dict[str, Hyper] = {}
    for i, target in enumerate(targets):
        t0 = time.time()
        others = [t for j, t in enumerate(targets) if j != i]
        if verbose:
            print(f"\n--- held out: {target.name} "
                  f"(tuned on {', '.join(t.name for t in others)}) ---")

        avail = np.array(sorted(set(target.names) & set(source.names)))
        sel = np.sort(rng.choice(avail, n_eval, replace=False)) \
            if avail.size > n_eval else avail
        truth = de_truth(target, sel)
        if verbose:
            print(f"    {avail.size:,} knockdowns shared with the screen, "
                  f"{sel.size} scored; measured DE genes/knockdown: median "
                  f"{np.median(truth.sig.sum(1)):.0f}", flush=True)

        hp = tune(source, others, symbols, n_eval=tune_eval, seed=seed,
                  verbose=verbose)
        chosen[target.name] = hp

        models = {
            "control (delta=0)": ControlBaseline(),
            "global mean [challenge baseline]": GlobalMeanBaseline(),
            "naive transfer": NaiveTransferBaseline(),
            "ContextTransfer (single source)": ContextTransferModel(hp),
        }
        fold = {}
        for name, model in models.items():
            model.fit([source], target.mu, symbols)
            fold[name] = evaluate(model.predict(sel), target, sel, truth, symbols)
        results[target.name] = fold

        if verbose:
            base = fold["global mean [challenge baseline]"]
            for name, m in fold.items():
                print(f"    {name:<34} PDS={m['discrimination_score_l1']:.3f}  "
                      f"ovl@100={m['overlap_at_100']:.3f}  "
                      f"MAE={m['mae']:.4f}  "
                      f"score={M.vcc_score(m, base)['avg_score']:.3f}  "
                      f"balanced={M.vcc_score(m, base, clip=False)['avg_score']:+.3f}")
            print(f"    ({time.time() - t0:.0f}s)", flush=True)

        Path("results").mkdir(exist_ok=True)
        Path("results/gwps_single_source.json").write_text(json.dumps(
            {"results": results, "hyper": {k: asdict(v) for k, v in chosen.items()}},
            indent=2, default=str))

    return {"results": results, "hyper": chosen}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmark", action="store_true",
                    help="run the leave-one-line-out single-source benchmark")
    ap.add_argument("--n-eval", type=int, default=400)
    ap.add_argument("--tune-eval", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.benchmark:
        run(n_eval=args.n_eval, tune_eval=args.tune_eval, seed=args.seed)
        return

    src, targets, genes, symbols = frame()
    print(src)
    for t in targets:
        print(f"  {t.name}: {len(set(t.names) & set(src.names)):,} knockdowns "
              f"shared with the screen")
    print(f"shared gene axis: {genes.size:,}")


if __name__ == "__main__":
    main()
