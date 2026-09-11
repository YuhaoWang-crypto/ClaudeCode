#!/usr/bin/env python
"""Ablation ladder: does each proposed change to ProteinTalks actually help?

Starts from the released architecture and adds one change at a time, so any
improvement is attributable. Every configuration is trained with the *same*
budget, optimiser and seeds; only the flags differ.

Metrics
-------
``skill``      1 - MSE/MSE_nochange on the 6/24/48 h proteome. Positive means the
               model beats "the proteome does not move"; the released
               architecture scores about -105 here, which is the problem the
               residual decoder is meant to fix.
``delta_r``    Pearson correlation between predicted and true *change* from
               baseline. Immune to the static component of the proteome.
``AUROC``      drug-efficacy classification, the paper's headline task.

Run:
    python scripts/optimize_model.py --n-proteins 300 --epochs 80 --seeds 0 1
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteintalks import (  # noqa: E402
    SyntheticPerturbationProteome, classification_metrics, make_splits,
    trajectory_metrics,
)
from proteintalks.baselines import DrugMeanBaseline  # noqa: E402
from proteintalks.improved import ProteinTalksR  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

# The ladder. Each row adds one change to the row above it.
LADDER = [
    ("released architecture", dict(
        delta_decode=False, real_time=False, time_conditioned=False,
        pert_conditioned=False, coupling_rank=0, head_rank=0, substeps=1)),
    ("+ residual decoding", dict(
        delta_decode=True, real_time=False, time_conditioned=False,
        pert_conditioned=False, coupling_rank=0, head_rank=0, substeps=1)),
    ("+ real elapsed time", dict(
        delta_decode=True, real_time=True, time_conditioned=False,
        pert_conditioned=False, coupling_rank=0, head_rank=0, substeps=2)),
    ("+ f(z,t,D) conditioning", dict(
        delta_decode=True, real_time=True, time_conditioned=True,
        pert_conditioned=True, coupling_rank=0, head_rank=0, substeps=2)),
    ("+ protein coupling", dict(
        delta_decode=True, real_time=True, time_conditioned=True,
        pert_conditioned=True, coupling_rank=8, head_rank=0, substeps=2)),
    ("+ low-rank head (full)", dict(
        delta_decode=True, real_time=True, time_conditioned=True,
        pert_conditioned=True, coupling_rank=8, head_rank=64, substeps=2)),
]


def to_t(a, device):
    return torch.tensor(a, device=device)


def train_one(ds, tr, va, cfg, args, seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = ProteinTalksR(
        pro_feats=ds.n_proteins,
        hidden_feats=args.hidden,
        drug_feature_feats=ds.drug_feats.shape[-1],
        dropout=0.0,
        **cfg,
    )
    # The released training configuration (config.py): AdamW, 5e-4, wd 1e-4,
    # batch 64, fixed lambda = 0.8. The SI confirms the fixed weighting is what
    # the reported results used.
    opt = torch.optim.AdamW(model.parameters(), lr=5e-4, weight_decay=1e-4)
    mse, bce = nn.MSELoss(), nn.BCEWithLogitsLoss()
    lam = 0.8

    def pack(idx):
        return (to_t(ds.p0[idx][..., None], "cpu"), to_t(ds.pert[idx][..., None], "cpu"),
                to_t(ds.drug_feats[idx][:, 0][..., None], "cpu"),
                to_t(ds.drug_feats[idx][:, 1][..., None], "cpu"),
                to_t(ds.p_future[idx], "cpu"), to_t(ds.label[idx], "cpu"))

    Xtr, Xva = pack(tr), pack(va)
    best, best_state, bad = np.inf, None, 0
    for ep in range(args.epochs):
        model.train()
        perm = torch.randperm(len(tr))
        for i in range(0, len(perm), args.batch):
            b = perm[i : i + args.batch]
            x, p, a, c, y, lab = (t[b] for t in Xtr)
            opt.zero_grad()
            pf, logit = model(x, p, a, c)
            loss = (1 - lam) * mse(pf, y) + lam * bce(logit, lab)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
        model.eval()
        with torch.no_grad():
            x, p, a, c, y, lab = Xva
            pf, logit = model(x, p, a, c)
            v = (1 - lam) * mse(pf, y).item() + lam * bce(logit, lab).item()
        if v < best - 1e-6:
            best, bad = v, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            bad += 1
            if bad >= args.patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model


@torch.no_grad()
def evaluate(model, ds, idx):
    model.eval()
    x = to_t(ds.p0[idx][..., None], "cpu")
    p = to_t(ds.pert[idx][..., None], "cpu")
    a = to_t(ds.drug_feats[idx][:, 0][..., None], "cpu")
    c = to_t(ds.drug_feats[idx][:, 1][..., None], "cpu")
    pf, logit = model(x, p, a, c)
    prob = torch.sigmoid(logit).numpy()
    tm = trajectory_metrics(ds.p_future[idx], pf.numpy(), ds.p0[idx])
    cm = classification_metrics(ds.label[idx], prob)
    return tm, cm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-proteins", type=int, default=300)
    ap.add_argument("--n-cell-lines", type=int, default=18)
    ap.add_argument("--n-drugs", type=int, default=62)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--patience", type=int, default=25)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--seeds", type=int, nargs="*", default=[0, 1])
    ap.add_argument("--tag", type=str, default="optimize")
    args = ap.parse_args()

    torch.set_num_threads(os.cpu_count() or 4)
    sim = SyntheticPerturbationProteome(
        n_proteins=args.n_proteins, n_cell_lines=args.n_cell_lines,
        n_drugs=args.n_drugs, seed=0)
    ds = sim.generate()
    print(f"[data] {ds.summary()}")

    tr, va, te = make_splits(ds, setting=1, seed=0)
    dm = DrugMeanBaseline().fit(ds, np.concatenate([tr, va]))
    ctrl = classification_metrics(ds.label[te], dm.predict_proba(ds, te))
    print(f"[control] drug-mean AUROC={ctrl['auroc']:.3f}\n")

    print(f"{'configuration':26s} {'skill':>9s} {'delta_r':>9s} {'AUROC':>8s} "
          f"{'params':>9s} {'sec':>6s}")
    print("-" * 76)

    rows = []
    for name, cfg in LADDER:
        t0 = time.time()
        sk, dr, au = [], [], []
        n_par = None
        for seed in args.seeds:
            model = train_one(ds, tr, va, cfg, args, seed)
            n_par = model.n_parameters()
            tm, cm = evaluate(model, ds, te)
            sk.append(tm["skill_vs_nochange"])
            dr.append(tm["delta_pearson_r"])
            au.append(cm["auroc"])
        row = {
            "config": name, "flags": cfg,
            "skill": float(np.mean(sk)), "skill_sd": float(np.std(sk)),
            "delta_r": float(np.nanmean(dr)),
            "auroc": float(np.mean(au)), "auroc_sd": float(np.std(au)),
            "params": n_par, "seconds": time.time() - t0,
            "n_seeds": len(args.seeds),
        }
        rows.append(row)
        print(f"{name:26s} {row['skill']:9.3f} {row['delta_r']:9.3f} "
              f"{row['auroc']:8.3f} {n_par['total']:9,d} {row['seconds']:6.0f}")

    out = {"control_auroc": ctrl["auroc"], "config": vars(args), "ladder": rows}
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"results_{args.tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[saved] {path}")
    print("\nskill = 1 - MSE/MSE_nochange; positive beats 'the proteome does not move'.")


if __name__ == "__main__":
    main()
