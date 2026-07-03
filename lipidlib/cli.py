"""Command-line featurizer: CSV of SMILES -> feature matrix.

Usage:
    python -m lipidlib.cli INPUT.csv --smiles-col smiles --fp morgan -o out.npy

Writes an ``.npy`` matrix (aligned to input row order) and prints a short report.
With ``--csv`` it instead writes a labelled CSV (feature columns named).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from .featurize import FEATURIZERS, featurize_many, descriptor_names, canonicalize


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="input CSV containing a SMILES column")
    p.add_argument("--smiles-col", default="smiles", help="name of the SMILES column")
    p.add_argument("--fp", default="morgan", choices=list(FEATURIZERS),
                   help="featurizer to use")
    p.add_argument("-o", "--out", default="features.npy", help="output path (.npy or .csv)")
    p.add_argument("--csv", action="store_true", help="write labelled CSV instead of .npy")
    args = p.parse_args(argv)

    df = pd.read_csv(args.input)
    if args.smiles_col not in df.columns:
        print(f"error: column {args.smiles_col!r} not in {list(df.columns)}", file=sys.stderr)
        return 2

    smiles = df[args.smiles_col].astype(str).tolist()
    X = featurize_many(smiles, kind=args.fp)

    n_bad = int(np.isnan(X).any(axis=1).sum())
    print(f"featurized {len(smiles)} molecules with '{args.fp}' -> {X.shape}")
    if n_bad:
        print(f"  warning: {n_bad} row(s) failed to parse (all-NaN)")

    if args.csv or args.out.endswith(".csv"):
        cols = descriptor_names(args.fp)
        out = pd.concat([df.reset_index(drop=True), pd.DataFrame(X, columns=cols)], axis=1)
        out.to_csv(args.out, index=False)
    else:
        np.save(args.out, X)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
