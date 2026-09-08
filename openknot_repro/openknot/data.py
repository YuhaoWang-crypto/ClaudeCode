"""Loading helpers for the OpenKnot benchmark release (eternagame/OpenKnotAIDesignData v4.5.x)."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

DATA_ROOT = Path(os.environ.get("OPENKNOT_DATA", Path(__file__).resolve().parent.parent / "data"))

BENCH_CSV = DATA_ROOT / "Data" / "OpenKnotBench_data.v4.5.1.csv"
M2R_CSV = DATA_ROOT / "Data" / "OK7a_M2R_data.v4.5.1.csv"
M2_CSV = DATA_ROOT / "Data" / "OK7a_M2_data.v4.5.2.csv"

META_COLUMNS = [
    "id",
    "sequence",
    "reads",
    "signal_to_noise",
    "SN_filter",
    "round",
    "puzzle",
    "method",
    "target_openknot_score",
    "sub_start",
    "sub_end",
    "design_length",
    "design_sequence",
    "target_structure",
    "RNet_structure",
    "RNet_F1",
    "RNet_F1_crossed_pair",
]

# Design methods compared in the paper, in the order used for the figures.
AI_METHODS = [
    "Rosetta",
    "Rosetta-LoRes",
    "3DRNA",
    "MPNN-fixbb",
    "MPNN-RFdiff",
    "codesign-RFdiff",
    "gRNAde",
    "gRNAde-no3d",
    "Struct2SeQ",
    "Struct2SeQ-SHAPE",
]
HUMAN_METHOD = "Eterna"
BASELINE_METHOD = "Starting sequence"


def reactivity_columns(csv_path: Path) -> list[str]:
    header = pd.read_csv(csv_path, nrows=0)
    return [
        c
        for c in header.columns
        if c.startswith("reactivity_") and not c.startswith("reactivity_error_")
    ]


def load_benchmark(
    csv_path: Path = BENCH_CSV,
    with_reactivity: bool = True,
    nrows: int | None = None,
) -> tuple[pd.DataFrame, np.ndarray | None]:
    """Load the benchmark table, optionally with the reactivity matrix.

    Returns ``(meta, reactivity)`` where ``reactivity`` is a float array of
    shape (n_designs, n_positions) with NaN for unprobed positions, or None.
    """
    if not csv_path.exists():
        raise FileNotFoundError(
            f"{csv_path} not found. Run scripts/download_data.sh first."
        )
    rcols = reactivity_columns(csv_path) if with_reactivity else []
    usecols = META_COLUMNS + rcols
    frame = pd.read_csv(csv_path, usecols=usecols, nrows=nrows, low_memory=False)
    meta = frame[META_COLUMNS].copy()
    reactivity = frame[rcols].to_numpy(dtype=float) if with_reactivity else None
    return meta, reactivity


def load_targets(round_name: str) -> pd.DataFrame:
    """Target secondary structures given to designers. round_name: '1and2', '3' or '4'."""
    name = {
        "1and2": "Rounds1and2_targets.csv",
        "3": "Round3_targets.csv",
        "4": "Round4_targets.csv",
    }[round_name]
    return pd.read_csv(DATA_ROOT / "Targets" / name)
