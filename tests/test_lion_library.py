"""Regression test: our LiON screen-library builder must match the exact column
schema the upstream LNP_ML `predict` command expects.

Run:  python -m pytest tests/ -q     (or just `python tests/test_lion_library.py`)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lipidlib.lion_library import build_library, Formulation, BASELINES, EXTRA_X_COLUMNS

REF = (Path(__file__).resolve().parents[1] / "external" / "LNP_ML" / "data"
       / "libraries" / "test_screen" / "test_screen_extra_x.csv")

SMILES = ["CCN(C)C", "CCCCCCCCN(C)C"]


def test_extra_x_columns_constant_matches_reference():
    if not REF.exists():
        return  # LNP_ML not cloned in this env; skip silently
    ref_cols = list(pd.read_csv(REF, nrows=1).columns)
    assert EXTRA_X_COLUMNS == ref_cols, "extra_x schema drifted from LNP_ML"


def test_build_library_shapes_and_onehots(tmp_path):
    d = build_library(SMILES, tmp_path, name="t", formulation=BASELINES["liver_IV"])
    smi = pd.read_csv(d / "t.csv")
    ex = pd.read_csv(d / "t_extra_x.csv")
    assert len(smi) == len(ex) == len(SMILES)
    assert list(ex.columns) == EXTRA_X_COLUMNS
    # every one-hot group sums to exactly 1 in each row
    for prefix in ["Delivery_target", "Helper_lipid_ID", "Route_of_administration",
                   "Cargo_type", "Model_type", "Batch_or_individual_or_barcoded"]:
        cols = [c for c in ex.columns if c.startswith(prefix)]
        assert (ex[cols].sum(axis=1) == 1).all(), f"{prefix} one-hot invalid"
    # liver_IV context encoded correctly
    assert ex["Delivery_target_liver"].iloc[0] == 1
    assert ex["Route_of_administration_intravenous"].iloc[0] == 1


def test_bad_formulation_rejected():
    import pytest
    with pytest.raises(ValueError):
        Formulation(delivery_target="brain").validate()  # not an allowed target
    with pytest.raises(ValueError):
        Formulation(cationic_mol=10, phospholipid_mol=10,
                    cholesterol_mol=10, peg_mol=10).validate()  # ratios != 100


if __name__ == "__main__":
    test_extra_x_columns_constant_matches_reference()
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_build_library_shapes_and_onehots(Path(td))
    test_bad_formulation_rejected()
    print("all tests passed")
