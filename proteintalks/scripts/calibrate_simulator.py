#!/usr/bin/env python
"""Check that the simulator's label structure resembles the real one.

Target statistics, measured on the paper's own label matrix (Supplementary
Table S1, sheet C_Efficacy; 1,116 conditions, 18 cell lines x 62 perturbations):

    positive rate            0.334
    drug-mean AUROC          0.909   (leave-out random split)
    cell-line-mean AUROC     0.463

The point is that in the real data, *drug identity carries almost all of the
label information and cell-line identity carries essentially none*. A simulator
that inverts this makes every downstream comparison meaningless, so it is worth
checking explicitly rather than assuming.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks import SyntheticPerturbationProteome, make_splits  # noqa: E402
from proteintalks.baselines import CellLineMeanBaseline, DrugMeanBaseline  # noqa: E402
from proteintalks.evaluate import classification_metrics  # noqa: E402

REAL = {"positive_rate": 0.334, "drug_mean_auroc": 0.909, "cellline_mean_auroc": 0.463}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-proteins", type=int, default=600)
    ap.add_argument("--cl-sd", type=float, nargs="*", default=[0.05, 0.15, 0.3, 0.6])
    ap.add_argument("--n-rep", type=int, default=20)
    args = ap.parse_args()

    print(f"target (real data):  positive rate {REAL['positive_rate']:.3f}  "
          f"drug-mean AUROC {REAL['drug_mean_auroc']:.3f}  "
          f"cell-line-mean AUROC {REAL['cellline_mean_auroc']:.3f}\n")
    print(f"{'cl_resistance_sd':>18s} {'pos rate':>9s} {'drug-mean':>11s} {'cl-mean':>9s}")

    for sd in args.cl_sd:
        sim = SyntheticPerturbationProteome(
            n_proteins=args.n_proteins, n_cell_lines=18, n_drugs=62,
            cl_resistance_sd=sd, seed=0,
        )
        ds = sim.generate()
        dm, cm = [], []
        for r in range(args.n_rep):
            tr, va, te = make_splits(ds, setting=1, seed=r)
            tr = np.concatenate([tr, va])
            dm.append(classification_metrics(
                ds.label[te], DrugMeanBaseline().fit(ds, tr).predict_proba(ds, te))["auroc"])
            cm.append(classification_metrics(
                ds.label[te], CellLineMeanBaseline().fit(ds, tr).predict_proba(ds, te))["auroc"])
        print(f"{sd:18.2f} {ds.label.mean():9.3f} {np.nanmean(dm):11.3f} {np.nanmean(cm):9.3f}")


if __name__ == "__main__":
    main()
