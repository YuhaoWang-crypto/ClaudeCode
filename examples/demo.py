"""End-to-end sanity demo: featurize reference ionizable lipids and show that the
embeddings behave sensibly (e.g. Tanimoto/cosine similarity separates lipids from
the small helper sterol).

Run after building the dataset:
    python scripts/fetch_reference_lipids.py
    python examples/demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lipidlib import featurize_many, canonicalize

DATA = Path(__file__).resolve().parents[1] / "data" / "reference_lipids.csv"


def cosine_matrix(X: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    U = X / norm
    return U @ U.T


def main() -> None:
    df = pd.read_csv(DATA)
    print(f"loaded {len(df)} reference molecules from {DATA.name}\n")

    # canonicalize as a validity check
    df["canonical"] = df["smiles"].map(canonicalize)
    bad = df[df["canonical"].isna()]
    if len(bad):
        print("WARNING unparseable:", list(bad["name"]))

    X = featurize_many(df["smiles"], kind="morgan")
    print(f"Morgan features: {X.shape}")

    S = cosine_matrix(X)
    names = df["name"].tolist()
    print("\nMorgan cosine similarity (rounded):")
    header = "".join(f"{n[:8]:>9}" for n in names)
    print(" " * 14 + header)
    for i, n in enumerate(names):
        row = "".join(f"{S[i, j]:9.2f}" for j in range(len(names)))
        print(f"{n:>14}{row}")


if __name__ == "__main__":
    main()
