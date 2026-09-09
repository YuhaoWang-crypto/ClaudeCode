#!/usr/bin/env python
"""Full reproduction run: the paper's three evaluation settings, its baselines,
the trivial controls it omits, and architecture ablations.

Usage
-----
    python scripts/run_experiments.py --quick        # small, ~minutes
    python scripts/run_experiments.py                # default scale
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteintalks import (  # noqa: E402
    SyntheticPerturbationProteome,
    bootstrap_ci,
    classification_metrics,
    make_splits,
    predict,
    train_proteintalks,
    trajectory_metrics,
)
from proteintalks.baselines import (  # noqa: E402
    CellLineMeanBaseline,
    DrugMeanBaseline,
    GlobalPriorBaseline,
    flatten_features,
    sklearn_baselines,
    train_deepsynergy,
)

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(RESULTS, exist_ok=True)


def build_dataset(args):
    sim = SyntheticPerturbationProteome(
        n_proteins=args.n_proteins,
        n_cell_lines=args.n_cell_lines,
        n_drugs=args.n_drugs,
        n_moa=args.n_moa,
        seed=args.seed,
    )
    ds = sim.generate()
    print(f"[data] {ds.summary()}")
    return ds


# --------------------------------------------------------------------------- #
def run_setting1(ds, args):
    print("\n=== Setting 1: random 0.7/0.2/0.1 split over conditions ===")
    tr, va, te = make_splits(ds, setting=1, seed=args.seed)
    print(f"[split] train={len(tr)} val={len(va)} test={len(te)}")
    rows = []

    # -- trivial controls (no proteomics at all) --
    for name, cls in [
        ("Drug-mean (no proteome)", DrugMeanBaseline),
        ("Cell-line-mean (no proteome)", CellLineMeanBaseline),
        ("Constant prior", GlobalPriorBaseline),
    ]:
        b = cls().fit(ds, tr)
        p = b.predict_proba(ds, te)
        m = classification_metrics(ds.label[te], p)
        m["model"] = name
        m["family"] = "trivial control"
        rows.append(m)
        print(f"  {name:32s} AUROC={m['auroc']:.3f} AUPRC={m['auprc']:.3f} ACC={m['accuracy']:.3f}")

    # -- the paper's classical comparators --
    xtr = flatten_features(ds, tr)
    xva = flatten_features(ds, va)
    xte = flatten_features(ds, te)
    for name, clf in sklearn_baselines(args.seed).items():
        t0 = time.time()
        clf.fit(xtr, ds.label[tr].astype(int))
        p = clf.predict_proba(xte)[:, 1]
        m = classification_metrics(ds.label[te], p)
        m.update(model=name, family="paper baseline", seconds=time.time() - t0)
        rows.append(m)
        print(f"  {name:32s} AUROC={m['auroc']:.3f} AUPRC={m['auprc']:.3f} ACC={m['accuracy']:.3f}")

    # -- DeepSynergy --
    t0 = time.time()
    p = train_deepsynergy(
        xtr, ds.label[tr], xva, ds.label[va], xte,
        hidden=args.ds_hidden, epochs=args.ds_epochs, seed=args.seed,
    )
    m = classification_metrics(ds.label[te], p)
    m.update(model="DeepSynergy", family="paper baseline", seconds=time.time() - t0)
    rows.append(m)
    print(f"  {'DeepSynergy':32s} AUROC={m['auroc']:.3f} AUPRC={m['auprc']:.3f} ACC={m['accuracy']:.3f}")

    # -- ProteinTalks and its ablations --
    for backbone, adaptive, label in [
        ("neural_ode", True, "ProteinTalks (neural ODE)"),
        ("mlp", True, "Ablation: MLP dynamics"),
        ("identity", True, "Ablation: no dynamics"),
        ("neural_ode", False, "Ablation: fixed task weights"),
    ]:
        t0 = time.time()
        model, hist = train_proteintalks(
            ds, tr, va,
            dynamics_backbone=backbone,
            adaptive_weights=adaptive,
            epochs=args.epochs,
            lam=args.lam,
            seed=args.seed,
            verbose=args.verbose,
        )
        pf, prob = predict(model, ds, te)
        m = classification_metrics(ds.label[te], prob)
        tm = trajectory_metrics(ds.p_future[te], pf, ds.p0[te])
        lo, hi = bootstrap_ci(ds.label[te], prob, "auroc")
        m.update(
            model=label, family="ProteinTalks",
            seconds=time.time() - t0,
            n_params=hist["n_params"]["total"],
            auroc_ci=[lo, hi], **{f"traj_{k}": v for k, v in tm.items()},
        )
        rows.append(m)
        print(
            f"  {label:32s} AUROC={m['auroc']:.3f} [{lo:.3f},{hi:.3f}] "
            f"AUPRC={m['auprc']:.3f} ACC={m['accuracy']:.3f} "
            f"| traj skill vs no-change={tm['skill_vs_nochange']:.3f} "
            f"delta_r={tm['delta_pearson_r']:.3f} | {hist['n_params']['total']} params"
        )
    return rows


# --------------------------------------------------------------------------- #
def run_setting2(ds, args):
    print("\n=== Setting 2: leave-one-cell-line-out ===")
    rows = []
    lines = sorted(set(ds.cell_line.tolist()))[: args.max_holdout]
    for l in lines:
        tr, va, te = make_splits(ds, setting=2, seed=args.seed, holdout=l)
        model, _ = train_proteintalks(
            ds, tr, va, epochs=args.epochs, lam=args.lam,
            seed=args.seed, verbose=False,
        )
        _, prob = predict(model, ds, te)
        m = classification_metrics(ds.label[te], prob)
        b = DrugMeanBaseline().fit(ds, tr)
        mb = classification_metrics(ds.label[te], b.predict_proba(ds, te))
        m.update(
            model="ProteinTalks", holdout=ds.cell_line_names[l],
            drugmean_auroc=mb["auroc"], drugmean_accuracy=mb["accuracy"],
        )
        rows.append(m)
        print(
            f"  {ds.cell_line_names[l]:8s} AUROC={m['auroc']:.3f} "
            f"AUPRC={m['auprc']:.3f} ACC={m['accuracy']:.3f} "
            f"| drug-mean AUROC={mb['auroc']:.3f}"
        )
    aur = [r["auroc"] for r in rows if not np.isnan(r["auroc"])]
    dm = [r["drugmean_auroc"] for r in rows if not np.isnan(r["drugmean_auroc"])]
    print(f"  --> ProteinTalks mean AUROC={np.mean(aur):.3f}; "
          f"cell lines with AUROC>=0.9: {sum(a >= 0.9 for a in aur)}/{len(aur)}")
    print(f"  --> drug-mean control mean AUROC={np.mean(dm):.3f}")
    return rows


# --------------------------------------------------------------------------- #
def run_setting3(ds, args):
    """Leave-one-drug-out, split by whether the drug's MOA is in the training set."""
    print("\n=== Setting 3: leave-one-drug-out (unseen chemistry) ===")
    rows = []
    drugs = sorted(set(ds.drug.tolist()))
    rng = np.random.default_rng(args.seed)
    drugs = list(rng.permutation(drugs))[: args.max_holdout]
    for k in drugs:
        tr, va, te = make_splits(ds, setting=3, seed=args.seed, holdout=k)
        moa_k = int(ds.moa[ds.drug == k][0])
        moa_seen = bool(np.any(ds.moa[tr] == moa_k))
        model, _ = train_proteintalks(
            ds, tr, va, epochs=args.epochs, lam=args.lam,
            seed=args.seed, verbose=False,
        )
        _, prob = predict(model, ds, te)
        m = classification_metrics(ds.label[te], prob)
        m.update(
            model="ProteinTalks", holdout=ds.drug_names[k],
            moa=moa_k, moa_in_training=moa_seen,
        )
        rows.append(m)
    seen = [r for r in rows if r["moa_in_training"]]
    unseen = [r for r in rows if not r["moa_in_training"]]
    for tag, grp in [("MOA seen in training", seen), ("MOA absent from training", unseen)]:
        if not grp:
            continue
        acc = np.mean([r["accuracy"] for r in grp])
        aur = np.nanmean([r["auroc"] for r in grp])
        print(f"  {tag:28s} n={len(grp):3d}  accuracy={acc:.3f}  AUROC={aur:.3f}")
    return rows


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-proteins", type=int, default=1200)
    ap.add_argument("--n-cell-lines", type=int, default=18)
    ap.add_argument("--n-drugs", type=int, default=63)
    ap.add_argument("--n-moa", type=int, default=20)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--ds-epochs", type=int, default=150)
    ap.add_argument("--ds-hidden", type=int, nargs=2, default=[2048, 2048])
    ap.add_argument("--lam", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-holdout", type=int, default=18)
    ap.add_argument("--settings", type=str, default="1,2,3")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--tag", type=str, default="main")
    args = ap.parse_args()

    if args.quick:
        args.n_proteins, args.n_cell_lines, args.n_drugs = 300, 8, 20
        args.epochs, args.ds_epochs, args.max_holdout = 40, 40, 3
        args.ds_hidden = [256, 256]

    args.ds_hidden = tuple(args.ds_hidden)
    ds = build_dataset(args)

    out = {"config": vars(args)}
    want = set(args.settings.split(","))
    if "1" in want:
        out["setting1"] = run_setting1(ds, args)
    if "2" in want:
        out["setting2"] = run_setting2(ds, args)
    if "3" in want:
        out["setting3"] = run_setting3(ds, args)

    path = os.path.join(RESULTS, f"results_{args.tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[saved] {path}")


if __name__ == "__main__":
    main()
