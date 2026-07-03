"""lipidlib — AI-guided lipid / targeting-ligand library screening.

Phase 1 focus: turn a library of SMILES into numerical features/embeddings that
downstream models (potency QSAR for Problem A, binding retrieval for Problem B)
can consume.
"""

from .featurize import (
    FEATURIZERS,
    featurize_smiles,
    featurize_many,
    canonicalize,
)

__all__ = [
    "FEATURIZERS",
    "featurize_smiles",
    "featurize_many",
    "canonicalize",
]

__version__ = "0.0.1"
