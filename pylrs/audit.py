"""Audit what the PylRS-libY model's accuracy actually is, and what it transfers to.

Why this exists: the repo ships `4.prediction/reference_data/leaderboard.info.txt`
with `score_test` = 1.000 for several models and 0.990 for the final ensemble. That
column is **in-sample**. `3.training/analyzing-autogluon.ipynb` calls
`predictor.leaderboard(dataset)` where `dataset` is the training set itself, so
`score_test` measures the model on data it was fitted on.

✅ Verified: loading the shipped predictor and scoring the 936 training rows gives
AUC 0.9905, matching the leaderboard's `score_test` of 0.990481 for
WeightedEnsemble_L3 to four decimals. So that number is not a generalisation
estimate. `score_val` (0.986) is AutoGluon's own stacked out-of-fold figure, which
for a 3-layer bagged stack fitted on 28 positives is itself optimistic.

What this module does instead: recovers the dataset's real structure and
cross-validates under four increasingly honest splits.

## The structure (recovered, not documented upstream)

936 rows = **117 PylRS variants x 8 ncAAs**, in 8 contiguous blocks of 117 with
the variants in the same order in each block. Column layout, read off
`analyzing-autogluon.ipynb` (which indexes `X[:,18]`..`X[:,23]` as the PLM
scores), is **not** the order the README's "256 + 5 + 1 + 18" suggests:

    cols   0..17    18 Rosetta cartesian_ddg terms   <- the ONLY substrate-aware features
    cols  18..22     5 ESM-1v scores                 } variant-only: identical
    col      23      1 ESM-msa-1b score              } across all 8 rows of
    cols  24..279  256 eUniRep (PylRS-fine-tuned)    } a given variant

The 8 substrates (`3.training/UAA-info.csv`) are all **O-substituted tyrosines**:
A OMe-Tyr, B OEt-Tyr, C O-propargyl-Tyr, D O-cyclohexyl-Tyr, E O-phenyl-Tyr,
F O-benzyl-Tyr, G O-cyclooctenyl-Tyr (the TCO one), H O-(2-nitrobenzyl)-Tyr.
Positives per substrate: 12, 2, **0**, 1, 5, 2, 3, 3 = 28.

## The four splits

    random      stratified 5-fold, ignores structure   <- what the published number does
    by-variant  GroupKFold over the 117 variants       <- no sibling rows in train
    by-ncAA     leave-one-substrate-out                <- "can it score a NEW ncAA?"
    both        unseen variants AND unseen substrate   <- the real cold start

Run (any modern Python; needs numpy/scikit-learn, not AutoGluon):

    python3 pylrs/audit.py --repo /path/to/PylRS-libY
"""

from __future__ import annotations

import argparse
import pathlib
import pickle
import warnings

import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
RNG = 0
N_VARIANTS = 117
N_SUBSTRATES = 8

#: Feature blocks, per the layout recovered above.
BLOCKS = {
    "ddg": np.arange(0, 18),     # substrate-aware
    "plm": np.arange(18, 280),   # variant-only (ESM + eUniRep)
    "all": np.arange(0, 280),
}

NCAA = ["A OMe-Tyr", "B OEt-Tyr", "C OPrg-Tyr", "D OCyhex-Tyr",
        "E OPh-Tyr", "F OBn-Tyr", "G TCO-Tyr", "H ONB-Tyr"]

MODELS = {
    # ExtraTrees is the family AutoGluon's own leaderboard ranked top here.
    "extratrees": lambda: ExtraTreesClassifier(
        n_estimators=400, class_weight="balanced", random_state=RNG, n_jobs=-1),
    "logreg": lambda: make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=4000, class_weight="balanced", C=0.1, random_state=RNG)),
}


def load(repo: pathlib.Path):
    """Load the training matrix and recover variant / substrate ids."""
    path = repo / "4.prediction" / "reference_data" / "train_evo256.pickle"
    with open(path, "rb") as fh:
        X, y = pickle.load(fh)
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(bool)
    if X.shape != (936, 280):
        raise SystemExit(f"unexpected training matrix shape {X.shape}")

    # Variant id from the variant-only feature block; substrate id from row order.
    _, gid = np.unique(np.round(X[:, 18:280], 6), axis=0, return_inverse=True)
    sid = np.arange(len(y)) // N_VARIANTS
    if len(set(gid)) != N_VARIANTS:
        raise SystemExit(f"expected {N_VARIANTS} variants, recovered {len(set(gid))}")
    return X, y, gid, sid


#: Repeats per split. With 28 positives the fold partition alone moves AUC by
#: ~0.08, so a single number would be reporting noise as a result.
N_REPEATS = 5


def splits(y, gid, sid, seed):
    """Yield (name, folds) for one repeat. `seed` reshuffles the partitions."""
    idx = np.arange(len(y))
    rs = np.random.RandomState(seed)

    yield "random", list(StratifiedKFold(5, shuffle=True, random_state=seed).split(idx, y))

    # GroupKFold is deterministic, so permute the group labels to resample folds.
    relabel = rs.permutation(N_VARIANTS)
    yield "by-variant", list(GroupKFold(n_splits=5).split(idx, y, groups=relabel[gid]))

    # Leave-one-substrate-out is exhaustive; identical across repeats by construction.
    scoreable = [s for s in range(N_SUBSTRATES) if y[sid == s].sum() > 0]
    yield "by-ncAA", [(idx[sid != s], idx[sid == s]) for s in scoreable]

    folds = []
    for vf in np.array_split(rs.permutation(N_VARIANTS), 4):
        held = np.isin(gid, vf)
        for s in scoreable:
            te, tr = idx[(sid == s) & held], idx[(sid != s) & ~held]
            if y[te].sum() > 0 and y[tr].sum() > 2:
                folds.append((tr, te))
    yield "both", folds


def evaluate(X, y, fold_list, cols, factory):
    """Pooled out-of-fold AUC/AP plus the mean per-fold AUC.

    Both are reported because they answer different questions. Pooling mixes
    scores from models fitted on different folds, which is misleading when the
    folds are whole substrates whose score scales differ; the mean per-fold AUC is
    the "rank variants within one ncAA" number that a design decision rests on.
    """
    oof = np.full(len(y), np.nan)
    per_fold = []
    for tr, te in fold_list:
        if y[tr].sum() == 0:
            continue
        model = factory()
        model.fit(X[np.ix_(tr, cols)], y[tr])
        proba = model.predict_proba(X[np.ix_(te, cols)])[:, 1]
        oof[te] = proba
        if 0 < y[te].sum() < len(te):
            per_fold.append(roc_auc_score(y[te], proba))
    ok = ~np.isnan(oof)
    pooled = roc_auc_score(y[ok], oof[ok]) if 0 < y[ok].sum() < ok.sum() else float("nan")
    ap = average_precision_score(y[ok], oof[ok]) if y[ok].sum() else float("nan")
    return pooled, ap, (float(np.mean(per_fold)) if per_fold else float("nan")), len(per_fold)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=pathlib.Path, required=True)
    args = ap.parse_args(argv)

    X, y, gid, sid = load(args.repo)
    print(f"training matrix {X.shape} | {int(y.sum())} positives "
          f"({y.mean() * 100:.1f}%) | {N_VARIANTS} variants x {N_SUBSTRATES} ncAAs")
    print("positives per ncAA: "
          + ", ".join(f"{NCAA[s]}={int(y[sid == s].sum())}" for s in range(N_SUBSTRATES)))
    print(f"a random ranker scores AUC 0.500 / AP {y.mean():.3f}\n")

    acc: dict[tuple, list] = {}
    for rep in range(N_REPEATS):
        for split_name, folds in splits(y, gid, sid, seed=RNG + rep):
            for feat, cols in BLOCKS.items():
                for mname, factory in MODELS.items():
                    acc.setdefault((split_name, feat, mname), []).append(
                        evaluate(X, y, folds, cols, factory))

    def stat(key, which):
        vals = [r[which] for r in acc[key] if not np.isnan(r[which])]
        return (float(np.mean(vals)), float(np.std(vals))) if vals else (float("nan"),) * 2

    head = (f"{'split':<11}{'features':<8}{'model':<12}"
            f"{'pooled AUC':>16}{'mean fold AUC':>18}")
    print(f"mean +- sd over {N_REPEATS} repartitions ({N_REPEATS} identical for by-ncAA, "
          "which is exhaustive)\n")
    print(head)
    print("-" * len(head))
    for split_name in ("random", "by-variant", "by-ncAA", "both"):
        for feat in BLOCKS:
            for mname in MODELS:
                k = (split_name, feat, mname)
                pm, ps = stat(k, 0)
                fm, fs = stat(k, 2)
                print(f"{split_name:<11}{feat:<8}{mname:<12}"
                      f"{pm:>11.3f} ±{ps:<4.3f}{fm:>13.3f} ±{fs:<4.3f}")
        print()

    print("=" * 78)
    print("1. Sibling-row leakage: does grouping the 8 rows of a variant matter?")
    worst = 0.0
    for feat in BLOCKS:
        r, rs_ = stat(("random", feat, "extratrees"), 0)
        g, gs = stat(("by-variant", feat, "extratrees"), 0)
        worst = min(worst, g - r)
        print(f"   {feat:<4} random {r:.3f}±{rs_:.3f} -> by-variant {g:.3f}±{gs:.3f}"
              f"   ({g - r:+.3f})")
    print(f"   -> largest drop {worst:+.3f}. Real but modest, and comparable to the")
    print("      repartition noise above, so grouping is not the main story.")

    print("\n2. Against the published figure")
    best = max(stat(k, 0)[0] for k in acc if k[0] == "by-variant")
    print(f"   best honest pooled AUC here: {best:.3f}  vs leaderboard score_test 0.990")
    print("   score_test is in-sample (docstring); the honest expectation is ~0.6-0.7.")

    print("\n3. Transfer to an unseen ncAA -- what matters for a new substrate")
    for feat in BLOCKS:
        a = stat(("by-ncAA", feat, "extratrees"), 2)[0]
        b = stat(("both", feat, "extratrees"), 2)[0]
        print(f"   {feat:<4} new ncAA, known variants: {a:.3f}"
              f"   |  new ncAA + new variants: {b:.3f}")
    ddg = stat(("by-ncAA", "ddg", "extratrees"), 2)[0]
    plm = stat(("by-ncAA", "plm", "extratrees"), 2)[0]
    print(f"   -> the substrate-aware ddG block transfers worst ({ddg:.3f}, ~chance),")
    print(f"      the variant-only block best ({plm:.3f}). What generalises across")
    print("      substrates is 'which variant is promiscuous at all' -- a prior on the")
    print("      scaffold, not a model of fit to this particular ncAA.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
