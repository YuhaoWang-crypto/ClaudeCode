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
from .data import (DATA, CellLine, load_all, load_replogle,
                   shared_perturbations)
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


JIANG_PRIOR = DATA / "pseudobulk" / "jiang_transferability.npz"


def jiang_transferability(symbols: np.ndarray, min_perturbations: int = 50
                          ) -> np.ndarray | None:
    """The same prior estimated on Jiang's six cell lines instead of four.

    Built by :mod:`virtualcell.prep_jiang`.  Two reasons to prefer it where it
    is available: six lines put the pure-noise floor near 0.17 rather than the
    0.33 that three lines impose, which widens the usable range; and four of the
    six are epithelial, which is what the official contexts B and C are and what
    the four-line atlas contains none of.

    Genes it cannot speak for -- unmeasured, or with too few usable
    perturbations -- take the median, which is neutral.  Returns ``None`` if the
    prior has not been built.
    """
    if not JIANG_PRIOR.exists():
        return None
    z = np.load(JIANG_PRIOR, allow_pickle=True)
    frac, n = z["fraction"], z["n_perturbations"]
    ok = n >= min_perturbations
    lut = {g: float(f) for g, f, k in
           zip(z["genes"].astype(str), frac, ok) if k}
    med = float(np.median(frac[ok]))
    return np.array([lut.get(str(s), med) for s in symbols])


def gene_transferability(symbols: np.ndarray,
                         exclude: str | None = None) -> np.ndarray:
    """Per-gene fraction of response variance that survives a change of context.

    The four-line atlas covers *none* of the official 2026 knockdowns, which is
    what forced the switch to the genome-wide screen.  It can still answer a
    different question, and one the single-source screen cannot: for each
    measured gene, how much of its response is shared across contexts rather
    than specific to one.  Estimating that needs several contexts but not the
    *same* knockdowns as the target panel, so an atlas disjoint from the panel is
    still admissible evidence about it.

    ``var(mean over lines) / var(pooled)`` per gene, the same definition the data
    package exports.  Genes absent from the atlas get the median, which is
    neutral.  Pass ``exclude`` to drop a line -- necessary when the value is used
    while predicting that line, or the prior has seen the answer.

    Read the scale carefully: with ``k`` lines a gene whose response is pure
    noise still averages to ``1/k``, so with three lines the floor is 0.33 and
    the observed median of 0.405 sits only a little above it.  The quantity is
    therefore as much a per-gene *noise level* as a transferability, and the
    dynamic range is narrow.  Whether it earns its place is settled by
    :func:`frontier`, not by the argument for it.
    """
    lines, _, atlas_sym = load_all()
    if exclude is not None:
        lines = [cl for cl in lines if cl.name != exclude]
    common = shared_perturbations(lines)
    stack = np.stack([cl.delta[np.array([{n: i for i, n in enumerate(cl.names)}[p]
                                         for p in common])] for cl in lines])
    shared = stack.mean(axis=0).var(axis=0)
    total = stack.reshape(-1, stack.shape[-1]).var(axis=0)
    frac = shared / np.maximum(total, 1e-12)

    lut: dict[str, float] = {}
    for s, f in zip(atlas_sym.astype(str), frac):
        lut.setdefault(s, float(f))
    med = float(np.median(frac))
    return np.array([lut.get(str(s), med) for s in symbols])


# --------------------------------------------------------------------------
# how good could any model be?
# --------------------------------------------------------------------------

def replicate_ceiling(verbose: bool = True) -> dict:
    """Agreement between the two K562 arms, which bounds what any model can score.

    Replogle ran K562 twice: an essential-gene arm and a genome-wide arm, same
    cell line, same lab, same CRISPRi library, same processing.  2,053 genes are
    targeted in both.  Their disagreement is therefore *measurement*, not
    biology, and it is the ceiling a cross-context model is being asked to beat
    -- without it, a correlation of 0.33 reads as a failure rather than as
    saturation.

    Nothing here involves the model.  It is a property of the data.
    """
    ess = load_replogle(DATA / "replogle" / "K562_essential_raw_bulk.h5ad", "ess")
    gw = load_gwps()

    genes = np.array(sorted(set(ess.genes) & set(gw.genes)))
    ie = {x: i for i, x in enumerate(ess.genes)}
    ig = {x: i for i, x in enumerate(gw.genes)}
    shared = np.array(sorted(set(ess.names) & set(gw.names)))
    re_ = {n: i for i, n in enumerate(ess.names)}
    rg = {n: i for i, n in enumerate(gw.names)}

    A = ess.delta[[re_[n] for n in shared]][:, [ie[x] for x in genes]]
    B = gw.delta[[rg[n] for n in shared]][:, [ig[x] for x in genes]]
    r = M.pearson_delta(A, B)

    out = {
        "n_knockdowns": int(shared.size), "n_genes": int(genes.size),
        "median_r": float(np.median(r)), "mean_r": float(np.mean(r)),
        "q25_r": float(np.percentile(r, 25)), "q75_r": float(np.percentile(r, 75)),
        "l2_essential": float(np.median(np.linalg.norm(A, axis=1))),
        "l2_genomewide": float(np.median(np.linalg.norm(B, axis=1))),
        "cells_essential": float(np.median(ess.ncells[[re_[n] for n in shared]])),
        "cells_genomewide": float(np.median(gw.ncells[[rg[n] for n in shared]])),
        "norm_cv_essential": M.norm_cv(A), "norm_cv_genomewide": M.norm_cv(B),
    }
    if verbose:
        print(f"{out['n_knockdowns']:,} knockdowns targeted in both K562 arms, "
              f"{out['n_genes']:,} shared genes\n")
        print(f"  effect agreement between arms   median r={out['median_r']:.3f} "
              f"(IQR {out['q25_r']:.3f}-{out['q75_r']:.3f})")
        print(f"  median L2 effect                essential {out['l2_essential']:.2f}"
              f", genome-wide {out['l2_genomewide']:.2f}")
        print(f"  median cells per pseudobulk     essential "
              f"{out['cells_essential']:.0f}, genome-wide "
              f"{out['cells_genomewide']:.0f}")
        print(f"  magnitude spread (norm CV)      essential "
              f"{out['norm_cv_essential']:.3f}, genome-wide "
              f"{out['norm_cv_genomewide']:.3f}")
        print("\nA cross-context model scoring Pearson ~0.33 is therefore at the "
              "level\nwhere the same screen repeated in the same cell line agrees "
              "with itself.")
    Path("results").mkdir(exist_ok=True)
    Path("results/replicate_ceiling.json").write_text(json.dumps(out, indent=2))
    return out


# --------------------------------------------------------------------------
# the magnitude frontier
# --------------------------------------------------------------------------

BETAS = (0.0, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0, 1.25, 1.5, 2.0)
# ``none`` is the no-prior arm; the two priors are otherwise identical in role
# but nearly uncorrelated with each other (Spearman 0.08 on the genes both
# cover), so which one is better is an empirical question, not a choice.
ARMS = (("none", 0.0), ("atlas", 0.5), ("atlas", 1.0),
        ("jiang", 0.5), ("jiang", 1.0), ("oracle", 0.5), ("oracle", 1.0))
BENCH_FILE = Path("results/gwps_single_source.json")


def build_prior(name: str, symbols: np.ndarray, exclude: str | None,
                target: CellLine | None = None, perts: np.ndarray | None = None
                ) -> np.ndarray | None:
    """The named per-gene transferability prior on this gene axis.

    ``oracle`` is not a usable prior -- it is built from the held-out truth and
    exists only to bound the route.  Without it, "this prior has saturated" and
    "this prior is weak" are indistinguishable, which is how a bounded gain gets
    mistaken for a ceiling.  It answers: if the per-gene weighting were perfect,
    how much would the whole mechanism be worth?
    """
    if name == "atlas":
        return gene_transferability(symbols, exclude=exclude)
    if name == "jiang":
        p = jiang_transferability(symbols)
        if p is None:
            raise SystemExit("run `python -m virtualcell.prep_jiang` first")
        return p
    if name == "oracle":
        if target is None or perts is None:
            raise ValueError("the oracle prior needs the held-out truth")
        index = {n: i for i, n in enumerate(target.names)}
        d = target.delta[[index[p] for p in perts]]
        # Per gene: how much of its measured response in the held-out line is
        # shared across these knockdowns rather than knockdown-specific --  the
        # same shape of quantity the real priors estimate, but read off the
        # answer.
        v = d.var(axis=0)
        return v / max(float(np.median(v)), 1e-12)
    return None


def _tuned_per_fold(names: list[str]) -> dict[str, Hyper]:
    """Each fold's cross-validated settings, or defaults if the run is missing."""
    if not BENCH_FILE.exists():
        print("  no benchmark output; sweeping on default settings instead")
        return {n: Hyper(temp=np.inf) for n in names}
    blob = json.loads(BENCH_FILE.read_text())["hyper"]
    out = {}
    for n in names:
        raw = blob.get(n)
        if raw is None:
            out[n] = Hyper(temp=np.inf)
            continue
        clean = {}
        for k, v in raw.items():
            d = getattr(Hyper(), k)
            clean[k] = (bool(v) if isinstance(d, bool) else
                        int(float(v)) if isinstance(d, int) else float(v))
        out[n] = Hyper(**clean)
    return out


def frontier(n_eval: int = 400, seed: int = 0, verbose: bool = True,
             arms: tuple[tuple[str, float], ...] = ARMS) -> dict:
    """Trace error against discrimination as predicted effects are scaled.

    With four source contexts the model cleared the challenge baseline on all
    three headline metrics at once.  With one it does not: transferring a K562
    effect into an unrelated line is, on average, *less* accurate than predicting
    no change at all, so any magnitude large enough to win discrimination and
    differential expression costs mean absolute error.

    That trade-off is a property of the data, not a tuning accident, and the
    challenge enforces minimum thresholds on every metric -- so the useful output
    is not one tuned point but the whole curve, and the largest scale that still
    keeps every metric at or better than baseline.

    Neither ``beta`` nor ``gene_w`` enters the consensus -- both act only in
    :meth:`ContextTransferModel.predict` -- so the entire sweep reuses one fit
    per fold.  ``gene_w`` weights each gene by how much of its response survives
    a change of context, measured on the four-line atlas; it is swept alongside
    the scale so the two are compared at each one's own best point rather than
    at a scale chosen for the other.
    """
    source, targets, genes, symbols = frame()
    rng = np.random.default_rng(seed)
    bank = SourceBank.build([source], np.ones(1))

    # Sweep on top of each fold's *tuned* settings, not the defaults.  The
    # MAE-against-scale curve depends on the smoothing and projection in force,
    # so a scale measured at default settings is not a scale that is safe at
    # tuned ones -- and it is the tuned ones a submission ships.  Each fold's
    # values were selected on the other two lines, so this stays leak-free.
    tuned = _tuned_per_fold([t.name for t in targets])

    out: dict[str, dict[str, dict]] = {}
    for target in targets:
        hp0 = tuned[target.name]
        avail = np.array(sorted(set(target.names) & set(source.names)))
        sel = np.sort(rng.choice(avail, n_eval, replace=False)) \
            if avail.size > n_eval else avail
        truth = de_truth(target, sel)
        base_m = evaluate(GlobalMeanBaseline().fit([source], target.mu, symbols)
                          .predict(sel), target, sel, truth, symbols)
        # a prior must not have seen the line it is used to predict; the Jiang
        # screen contains none of these four lines, so only the atlas one needs
        # the exclusion
        priors = {name: build_prior(name, symbols, exclude=target.name,
                                    target=target, perts=sel)
                  for name, _ in arms}

        # one fit per fold: beta and gene_w act only in predict()
        model = ContextTransferModel(hp0).fit(
            [source], target.mu, symbols, bank=bank)

        per_w: dict[str, dict] = {}
        for name, w in arms:
            key = f"{name}:{w:g}"
            model._gene_prior = priors[name]
            rows = {"beta": [], "pds": [], "ovl": [], "mae": [], "score": [],
                    "balanced": [], "mae_vs_base": [], "norm_cv": []}
            for b in BETAS:
                model.hp = replace(hp0, beta=b, gene_w=w)
                pred = model.predict(sel)
                m = evaluate(pred, target, sel, truth, symbols)
                rows["beta"].append(b)
                rows["pds"].append(m["discrimination_score_l1"])
                rows["ovl"].append(m["overlap_at_100"])
                rows["mae"].append(m["mae"])
                rows["mae_vs_base"].append(m["mae"] - base_m["mae"])
                rows["score"].append(M.vcc_score(m, base_m)["avg_score"])
                rows["balanced"].append(
                    M.vcc_score(m, base_m, clip=False)["avg_score"])
                rows["norm_cv"].append(M.norm_cv(pred))
            per_w[key] = rows

        index = {n: i for i, n in enumerate(target.names)}
        truth_cv = M.norm_cv(target.delta[[index[p] for p in sel]])
        out[target.name] = {"arm": per_w, "truth_norm_cv": truth_cv}
        if verbose:
            print(f"\n--- {target.name} ({sel.size} knockdowns; measured "
                  f"magnitude spread {truth_cv:.3f}) ---")
            for key, rows in per_w.items():
                print(f"  prior:gene_w = {key}")
                print(f"    {'beta':>5}{'PDS':>8}{'ovl@100':>10}{'MAE':>9}"
                      f"{'ΔMAE vs base':>14}{'norm CV':>10}{'score':>8}"
                      f"{'balanced':>10}")
                for i, b in enumerate(BETAS):
                    flag = "" if rows["mae_vs_base"][i] <= 0 else "  ✗ MAE"
                    print(f"    {b:>5.2f}{rows['pds'][i]:>8.3f}"
                          f"{rows['ovl'][i]:>10.3f}{rows['mae'][i]:>9.4f}"
                          f"{rows['mae_vs_base'][i]:>+14.4f}"
                          f"{rows['norm_cv'][i]:>10.3f}{rows['score'][i]:>8.3f}"
                          f"{rows['balanced'][i]:>+10.3f}{flag}", flush=True)

    best: dict[str, dict] = {}
    for name, w in arms:
        key = f"{name}:{w:g}"
        worst = [max(out[t]["arm"][key]["mae_vs_base"][i] for t in out)
                 for i in range(len(BETAS))]
        safe = [b for b, x in zip(BETAS, worst) if x <= 0]
        entry = {"safe_beta": max(safe) if safe else None,
                 "prior": name, "gene_w": w}
        if safe:
            i = BETAS.index(max(safe))
            for m in ("pds", "ovl", "score"):
                entry[m] = float(np.mean([out[t]["arm"][key][m][i] for t in out]))
        best[key] = entry
    if verbose:
        print("\nlargest scale keeping MAE at or below baseline on every fold")
        print(f"  {'prior:gene_w':>14}{'beta':>7}{'PDS':>8}{'ovl@100':>10}"
              f"{'score':>8}")
        for key, e in best.items():
            if e["safe_beta"] is None:
                print(f"  {key:>14}{'none':>7}")
                continue
            print(f"  {key:>14}{e['safe_beta']:>7.2f}{e['pds']:>8.3f}"
                  f"{e['ovl']:>10.3f}{e['score']:>8.3f}")

    Path("results").mkdir(exist_ok=True)
    Path("results/gwps_frontier.json").write_text(json.dumps(
        {"betas": list(BETAS), "arms": [list(a) for a in arms], "folds": out,
         "best": best}, indent=2))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmark", action="store_true",
                    help="run the leave-one-line-out single-source benchmark")
    ap.add_argument("--frontier", action="store_true",
                    help="trace error against discrimination over effect scale")
    ap.add_argument("--ceiling", action="store_true",
                    help="agreement between the two K562 arms, the noise ceiling")
    ap.add_argument("--n-eval", type=int, default=400)
    ap.add_argument("--tune-eval", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.ceiling:
        replicate_ceiling()
        return
    if args.frontier:
        frontier(n_eval=args.n_eval, seed=args.seed)
        return
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
