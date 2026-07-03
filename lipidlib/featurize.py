"""SMILES -> numerical features / embeddings.

Provides several featurizers useful as a Phase-1 baseline for both problems:

* ``morgan``   : Morgan / ECFP-style circular fingerprint (bit vector).
* ``maccs``    : 167-bit MACCS structural keys.
* ``rdkit2d``  : the full RDKit 2D descriptor panel (physchem properties).

These are deterministic, CPU-only, and dependency-light (RDKit + NumPy). Learned
embeddings (Chemprop D-MPNN, MolFormer, Uni-Mol) are added in a later phase and
will expose the same ``featurize_*`` interface so callers do not change.

All functions return ``numpy.ndarray`` (float32) and never raise on a single bad
SMILES: invalid entries yield an all-NaN row so a whole library still vectorizes.
"""

from __future__ import annotations

import logging
from typing import Callable, Iterable, Sequence

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator

# RDKit is noisy on odd SMILES; silence its C++ logger, we handle errors here.
RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


def canonicalize(smiles: str) -> str | None:
    """Return the canonical SMILES, or ``None`` if the string cannot be parsed."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


# --- individual featurizers ------------------------------------------------

def _morgan(mol: Chem.Mol, n_bits: int = 2048, radius: int = 2) -> np.ndarray:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.GetFingerprint(mol)
    arr = np.zeros((n_bits,), dtype=np.float32)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _maccs(mol: Chem.Mol) -> np.ndarray:
    fp = MACCSkeys.GenMACCSKeys(mol)
    arr = np.zeros((167,), dtype=np.float32)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


# Precompute the RDKit 2D descriptor list once.
_DESC_NAMES: list[str] = [name for name, _ in Descriptors.descList]
_DESC_FUNCS: list[Callable[[Chem.Mol], float]] = [fn for _, fn in Descriptors.descList]


def _rdkit2d(mol: Chem.Mol) -> np.ndarray:
    out = np.empty((len(_DESC_FUNCS),), dtype=np.float32)
    for i, fn in enumerate(_DESC_FUNCS):
        try:
            out[i] = fn(mol)
        except Exception:  # a few descriptors throw on exotic structures
            out[i] = np.nan
    return out


# name -> (function, output dimensionality)
FEATURIZERS: dict[str, tuple[Callable[[Chem.Mol], np.ndarray], int]] = {
    "morgan": (_morgan, 2048),
    "maccs": (_maccs, 167),
    "rdkit2d": (_rdkit2d, len(_DESC_NAMES)),
}


def descriptor_names(kind: str) -> list[str]:
    """Human-readable column names for a featurizer (useful for ``rdkit2d``)."""
    if kind == "rdkit2d":
        return list(_DESC_NAMES)
    _, dim = FEATURIZERS[kind]
    return [f"{kind}_{i}" for i in range(dim)]


# --- public API ------------------------------------------------------------

def featurize_smiles(smiles: str, kind: str = "morgan", **kwargs) -> np.ndarray:
    """Featurize a single SMILES. Returns an all-NaN vector on parse failure."""
    if kind not in FEATURIZERS:
        raise ValueError(f"unknown featurizer {kind!r}; choose from {list(FEATURIZERS)}")
    fn, dim = FEATURIZERS[kind]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("could not parse SMILES: %s", smiles)
        return np.full((dim,), np.nan, dtype=np.float32)
    if kwargs:
        return fn(mol, **kwargs)
    return fn(mol)


def featurize_many(
    smiles_list: Sequence[str] | Iterable[str],
    kind: str = "morgan",
    **kwargs,
) -> np.ndarray:
    """Featurize a library. Returns an ``(n, dim)`` float32 matrix.

    Rows for unparseable SMILES are all-NaN so the matrix stays rectangular and
    aligned with the input order.
    """
    smiles_list = list(smiles_list)
    rows = [featurize_smiles(s, kind=kind, **kwargs) for s in smiles_list]
    if not rows:
        _, dim = FEATURIZERS[kind]
        return np.empty((0, dim), dtype=np.float32)
    return np.vstack(rows)
