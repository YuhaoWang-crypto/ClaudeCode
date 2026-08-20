"""Score a submission the way the challenge does: DE from the cells you submit.

Everything else in this project scores at pseudobulk, comparing mean effects.
``cell-eval`` does not: it runs a **Wilcoxon rank-sum test on the submitted
cells** against real control cells, per gene per condition, and builds the
differential expression score from those calls.

That difference is not cosmetic and it does not move all models the same way.
A predictor that emits one mean replicated across 400 rows has *zero*
within-group variance, so the rank-sum test sees a perfectly separated
distribution and returns significance for any nonzero shift — a model that
predicts nothing can score well, and a good model can lose by submitting
honestly. This module measures that gap directly on held-out data, which is the
only way to know which side of it a given model sits on.

    python -m virtualcell.vcc_eval --target K562 --n-knockdowns 40

``--source gwps`` swaps the four-line atlas for Replogle's genome-wide K562
screen, which is the only source that covers the official 2026 panel and so the
only one a real submission can be built from.  Paired with ``--match-panel`` --
which picks knockdowns whose measured effect size matches the official panel's
rather than the much stronger essential-gene default -- it is an end-to-end dry
run of the submission pipeline, on data whose ground truth is not withheld.
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy import sparse, stats

from . import metrics as M
from .data import load_all, shared_perturbations
from .model import (ContextTransferModel, ControlBaseline, GlobalMeanBaseline,
                    NaiveTransferBaseline)
from .predict import DEFAULT_HYPER
from .vcc_submit import sample_cells

SC = "/home/user/vcc_data/singlecell"
TAG = {"K562": "k562", "RPE1": "rpe1", "HepG2": "hepg2", "Jurkat": "jurkat"}


def de_from_cells(pert: sparse.csr_matrix, ctrl: sparse.csr_matrix,
                  fdr: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Rank-sum DE of perturbed cells against control cells, as cell-eval does.

    Returns ``(significant, log2fc)``.  Both matrices are normalised per cell to
    counts-per-10k first, matching the scorer.
    """
    def norm(M: sparse.csr_matrix) -> np.ndarray:
        X = np.asarray(M.todense(), dtype=np.float64)
        totals = X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1.0
        return np.log1p(X / totals * 1e4)

    A, B = norm(pert), norm(ctrl)
    # Mann-Whitney U across all genes at once
    stat, p = stats.mannwhitneyu(A, B, axis=0, alternative="two-sided")
    lfc = (A.mean(axis=0) - B.mean(axis=0)) / M.LOG2
    order = np.argsort(p)
    n = p.size
    thresh = fdr * np.arange(1, n + 1) / n
    passing = p[order] <= thresh
    sig = np.zeros(n, dtype=bool)
    if passing.any():
        sig[order[:np.flatnonzero(passing)[-1] + 1]] = True
    return sig, lfc


def _panel_matched(source, candidates: np.ndarray, n: int,
                   rng: np.random.Generator) -> np.ndarray:
    """Pick knockdowns whose source-measured effect matches the official panel's.

    The essential-gene panel this project benchmarks on is unusually strong:
    median L2 effect 4.26 in K562 against 3.16 for the 272 official targets the
    genome-wide screen covers.  Scoring on it would overstate what a submission
    can expect.  Sampling candidates in proportion to how the official panel's
    magnitudes are distributed removes that bias.
    """
    from .vcc2026 import panel as official_panel

    idx = {p: i for i, p in enumerate(source.names)}
    mag = np.linalg.norm(source.delta, axis=1)
    ref = np.array([mag[idx[p]] for p in official_panel() if p in idx])
    if ref.size < 20:
        raise SystemExit(
            f"--match-panel needs a source that measures the official panel; "
            f"this one covers {ref.size} of 300. Use --source gwps.")
    edges = np.quantile(ref, np.linspace(0, 1, 6))
    edges[0], edges[-1] = -np.inf, np.inf

    cand_mag = np.array([mag[idx[c]] for c in candidates])
    out: list[str] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        pool = candidates[(cand_mag >= lo) & (cand_mag < hi)]
        if pool.size:
            out.extend(rng.choice(pool, min(n // 5, pool.size), replace=False))
    return np.sort(np.array(out))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="K562",
                    choices=["K562", "RPE1", "HepG2", "Jurkat"])
    ap.add_argument("--source", default="atlas", choices=["atlas", "gwps"],
                    help="four-line essential atlas, or the genome-wide K562 screen")
    ap.add_argument("--match-panel", action="store_true",
                    help="match the knockdowns' effect sizes to the official panel")
    ap.add_argument("--n-knockdowns", type=int, default=40)
    ap.add_argument("--n-cells", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import anndata as ad

    rng = np.random.default_rng(args.seed)
    if args.source == "gwps":
        from .gwps import frame
        if args.target == "K562":
            raise SystemExit("K562 is the source screen's own line; scoring it "
                             "would read back the training data")
        src, targets, genes, symbols = frame([args.target])
        target, sources = targets[0], [src]
        common = np.array(sorted(set(target.names) & set(src.names)))
    else:
        lines, genes, symbols = load_all()
        target = next(cl for cl in lines if cl.name == args.target)
        sources = [cl for cl in lines if cl.name != args.target]
        common = shared_perturbations(lines)

    if args.match_panel:
        perts = _panel_matched(sources[0], common, args.n_knockdowns, rng)
    else:
        perts = np.sort(rng.choice(common, args.n_knockdowns, replace=False))

    # real cells of the held-out line, on this project's gene axis
    cells = ad.read_h5ad(f"{SC}/{TAG[args.target]}.h5ad")
    lookup: dict[str, int] = {}
    for i, s in enumerate(cells.var["gene_name"].astype(str).values):
        lookup.setdefault(s, i)
    take = np.array([lookup.get(s, -1) for s in symbols])
    have = take >= 0
    cols = take[have]
    X = sparse.csr_matrix(cells.X)
    label = cells.obs["gene"].values

    ctrl_cells = X[label == "non-targeting"][:, cols]
    print(f"held out {args.target}: {ctrl_cells.shape[0]:,} control cells, "
          f"{have.sum():,} genes, {perts.size} knockdowns\n")

    mu = np.log1p(np.asarray(ctrl_cells.sum(0)).ravel()
                  / max(ctrl_cells.sum(), 1.0) * 1e4)

    # measured ground truth, from real perturbed cells vs real controls
    truth_sig, truth_lfc = [], []
    keep = []
    for p in perts:
        rows = np.flatnonzero(label == p)
        if rows.size < 10:
            continue
        s, l = de_from_cells(X[rows][:, cols], ctrl_cells)
        truth_sig.append(s); truth_lfc.append(l); keep.append(p)
    perts = np.array(keep)
    truth_sig = np.array(truth_sig); truth_lfc = np.array(truth_lfc)
    print(f"measured DE (rank-sum on real cells): median "
          f"{np.median(truth_sig.sum(1)):.0f} genes/knockdown\n")

    hp = DEFAULT_HYPER
    if args.source == "gwps":
        # the four-line defaults were tuned with three source contexts to average
        # over; the single-source run has its own tuned settings
        from .vcc2026 import single_source_hyper
        hp = single_source_hyper(exclude=args.target)

    models = {
        "control (delta=0)": ControlBaseline(),
        "global mean [challenge baseline]": GlobalMeanBaseline(),
        "naive transfer": NaiveTransferBaseline(),
        "ContextTransfer (ours)": ContextTransferModel(hp),
    }
    sub_sym = symbols[have]
    sub_sources = []
    from dataclasses import replace as dc_replace
    for s_ in sources:
        sub_sources.append(dc_replace(s_, ctrl=s_.ctrl[:, have],
                                      pert=s_.pert[:, have],
                                      genes=s_.genes[have],
                                      symbols=s_.symbols[have],
                                      _mu=None, _delta=None))

    head = (f"{'model':<34}{'submitted as':<16}{'DE ovl@100':>12}"
            f"{'DE genes':>10}{'direction':>11}")
    print(head); print("-" * len(head))
    rows_out = []
    for name, model in models.items():
        model.fit(sub_sources, mu, sub_sym)
        effects = model.predict(perts)

        for mode in ("cells", "replicated mean"):
            sigs, lfcs = [], []
            for i, p in enumerate(perts):
                if mode == "cells":
                    block = sample_cells(ctrl_cells, effects[i], args.n_cells,
                                         np.random.default_rng(args.seed + i),
                                         control_mean=mu)
                else:
                    # the degenerate case: one predicted mean, replicated
                    target_mean = np.expm1(np.maximum(mu + effects[i], 0.0))
                    row = np.round(target_mean * 100).astype(np.int32)
                    block = sparse.csr_matrix(
                        np.repeat(row[None, :], args.n_cells, axis=0))
                s, l = de_from_cells(block, ctrl_cells)
                sigs.append(s); lfcs.append(l)
            sigs = np.array(sigs); lfcs = np.array(lfcs)
            ovl = float(np.mean(M.overlap_at_k(lfcs, truth_lfc, truth_sig, k=100)))
            direction = float(np.nanmean(
                M.de_direction_match(lfcs, truth_lfc, truth_sig)))
            print(f"{name if mode=='cells' else '':<34}{mode:<16}{ovl:>12.4f}"
                  f"{np.mean(sigs.sum(1)):>10.0f}{direction:>11.3f}")
            rows_out.append((name, mode, ovl))
        print()

    print("The gap between the two rows of each model is the degenerate-variance")
    print("trap: replicated means have no within-group spread, so the rank-sum")
    print("test separates perfectly and calls far more genes significant than the")
    print("prediction earns. Submit cells.")


if __name__ == "__main__":
    main()
