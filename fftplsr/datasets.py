"""The PylRS datasets from the FP4COM paper, vendored as portable CSV/FASTA.

Source: Hu et al., "Machine learning-guided evolution of pyrrolysyl-tRNA
synthetase for improved incorporation efficiency of diverse noncanonical amino
acids", Nature Communications 16 (2025), doi:10.1038/s41467-025-61952-2, and its
code repository https://github.com/zjuhaoran/FPFORCOM (MIT). Converted from
`input/data.xlsx` and the two EnzymeML records in `data/` -- see
`data/PROVENANCE.md`.

Fitness is relative stop-codon-suppression efficiency, parent = 1.0.

Backgrounds
-----------
IFRS       the starting synthetase (454 aa)
Com1-IFRS  the round-2 winner, D2N/T56P/R61K/H62Y/H63Y relative to IFRS

Training rounds (as used in the paper)
--------------------------------------
trainset1  13   IFRS + 12 single mutants                    -> round 1 model
trainset2  38   trainset1 + 25 measured double/triple combos -> round 2 model (picked Com1)
trainset3  96   Com1-IFRS + 95 DL-nominated singles          -> site triage
trainset4  120  Com1-IFRS + 27 singles + 92 doubles          -> round 3 model (picked Com2)

Measured panels (every variant with a measurement, for held-out scoring)
------------------------------------------------------------------------
measured_ifrs_panel  102 variants on the IFRS background
measured_com1_panel  140 variants on the Com1-IFRS background
"""

from __future__ import annotations

import pathlib

import pandas as pd

DATA_DIR = pathlib.Path(__file__).parent / "data"

#: The 12 single mutations recombined in paper rounds 1-2 (IFRS background).
ROUND12_SITES = [
    "D2N", "K3N", "R19H", "H29R", "V31I", "T56P",
    "R61K", "H62Y", "H63Y", "A100E", "T122S", "S193R",
]

__all__ = [
    "DATA_DIR",
    "ROUND12_SITES",
    "load_sequence",
    "load_table",
    "ifrs",
    "com1",
    "trainset",
    "measured_panel",
]


def load_sequence(name: str) -> str:
    """Read a vendored FASTA and return the bare sequence."""
    path = DATA_DIR / f"{name}.fasta"
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def ifrs() -> str:
    """The IFRS parent sequence (454 aa)."""
    return load_sequence("IFRS")


def com1() -> str:
    """The Com1-IFRS parent sequence (454 aa)."""
    return load_sequence("Com1-IFRS")


def load_table(stem: str) -> pd.DataFrame:
    """Load a vendored CSV with `Variants` / `Fitness` columns."""
    df = pd.read_csv(DATA_DIR / f"{stem}.csv")
    return df[df["Variants"].notna()].reset_index(drop=True)


_TRAINSETS = {
    1: ("trainset1_ifrs_singles", "IFRS"),
    2: ("trainset2_ifrs_combos", "IFRS"),
    3: ("trainset3_com1_singles", "Com1-IFRS"),
    4: ("trainset4_com1_combos", "Com1-IFRS"),
}


def trainset(n: int) -> tuple[pd.DataFrame, str]:
    """Return ``(table, parent_sequence)`` for paper training set `n` (1-4)."""
    stem, background = _TRAINSETS[n]
    parent = ifrs() if background == "IFRS" else com1()
    return load_table(stem), parent


def measured_panel(background: str = "IFRS") -> tuple[pd.DataFrame, str]:
    """Every measured variant on `background` ('IFRS' or 'Com1-IFRS'), plus its sequence."""
    if background.upper() == "IFRS":
        return load_table("measured_ifrs_panel"), ifrs()
    if background.lower().startswith("com1"):
        return load_table("measured_com1_panel"), com1()
    raise KeyError(f"unknown background {background!r}")
