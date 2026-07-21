"""
prototypes.py -- crystal-structure prototypes + ionic data for the STRUCTURAL
(substitution) candidate pipeline, plus a small curated reference set for the demo.

GNoME's structural pipeline generates candidates by substituting ions into known
crystal prototypes. We encode a handful of classic prototypes with their site
"roles" (A/B/X cation/anion slots), common oxidation states, and Shannon ionic
radii so we can (a) enumerate chemically-sensible substitutions and (b) apply
cheap physical filters (charge neutrality, perovskite tolerance factor) BEFORE any
expensive scoring -- exactly the "reduce the candidate set" logic of the paper.
"""
from __future__ import annotations
from dataclasses import dataclass, field


# --- common oxidation states (subset relevant to oxide/chalcogenide/halide solids)
OXIDATION_STATES = {
    "Li": [1], "Na": [1], "K": [1], "Ag": [1], "Cu": [1, 2],
    "Mg": [2], "Ca": [2], "Sr": [2], "Ba": [2], "Zn": [2],
    "Al": [3], "Ga": [3], "Sc": [3], "Y": [3], "La": [3], "In": [3],
    "Ti": [4], "Zr": [4], "Sn": [4], "Ge": [4], "Si": [4], "Hf": [4],
    "Nb": [5], "Ta": [5], "V": [5],
    "P": [5], "As": [5], "W": [6], "Mo": [6],
    "O": [-2], "S": [-2], "Se": [-2],
    "F": [-1], "Cl": [-1], "Br": [-1], "I": [-1],
}

# Shannon ionic radii (Angstrom), 6-coordinate unless noted; enough for tolerance factor
SHANNON_RADII = {
    ("Li", 1): 0.76, ("Na", 1): 1.02, ("K", 1): 1.38, ("Ag", 1): 1.15,
    ("Mg", 2): 0.72, ("Ca", 2): 1.00, ("Sr", 2): 1.18, ("Ba", 2): 1.35, ("Zn", 2): 0.74,
    ("Al", 3): 0.535, ("Ga", 3): 0.62, ("Sc", 3): 0.745, ("Y", 3): 0.90,
    ("La", 3): 1.032, ("In", 3): 0.80,
    ("Ti", 4): 0.605, ("Zr", 4): 0.72, ("Sn", 4): 0.69, ("Ge", 4): 0.53,
    ("Hf", 4): 0.71, ("Nb", 5): 0.64, ("Ta", 5): 0.64,
    ("O", -2): 1.40, ("S", -2): 1.84, ("F", -1): 1.33, ("Cl", -1): 1.81,
}


@dataclass(frozen=True)
class Prototype:
    """A structural template with cation/anion roles and per-role site multiplicity."""
    name: str
    roles: tuple            # ordered role names, e.g. ("A","B","X")
    counts: dict            # role -> atoms per formula unit
    role_sign: dict         # role -> +1 cation / -1 anion
    formula_template: str   # display, e.g. "ABX3"
    note: str = ""


PROTOTYPES = {
    "rock_salt":   Prototype("rock_salt", ("A", "X"), {"A": 1, "X": 1},
                             {"A": +1, "X": -1}, "AX", "NaCl-type MX"),
    "perovskite":  Prototype("perovskite", ("A", "B", "X"), {"A": 1, "B": 1, "X": 3},
                             {"A": +1, "B": +1, "X": -1}, "ABX3",
                             "ABO3 perovskite (tolerance-factor filtered)"),
    "spinel":      Prototype("spinel", ("A", "B", "X"), {"A": 1, "B": 2, "X": 4},
                             {"A": +1, "B": +1, "X": -1}, "AB2X4", "AB2O4 spinel"),
    "antiperovskite": Prototype("antiperovskite", ("X", "B", "A"),
                                {"X": 3, "B": 1, "A": 1}, {"X": +1, "B": -1, "A": -1},
                                "X3BA", "Li3OCl-type Li-ion conductor"),
    "olivine":     Prototype("olivine", ("A", "B", "P", "O"),
                             {"A": 1, "B": 1, "P": 1, "O": 4},
                             {"A": +1, "B": +1, "P": +1, "O": -1}, "ABPO4",
                             "LiFePO4-type cathode"),
}


def tolerance_factor(rA: float, rB: float, rX: float) -> float:
    """Goldschmidt tolerance factor t for ABX3 perovskite; 0.8<t<1.0 ~ formable."""
    from math import sqrt
    return (rA + rX) / (sqrt(2.0) * (rB + rX))


# --- curated reference set for the DEMO: alkaline-earth titanate/zirconate oxides.
# References are the BINARY oxides only; the pipeline "discovers" the ternary ABO3
# perovskites as new convex-hull points below the binary tie-lines -- a faithful
# GNoME-style setup (known simple phases in, new ternaries out).
# ⚠️ ILLUSTRATIVE formation energies (eV/atom), order-of-magnitude literature values,
# NOT a substitute for Materials Project / DFT. Real runs replace these via score.py.
DEMO_REFERENCE_ENERGIES = {
    "CaO":  -3.28,
    "SrO":  -3.02,
    "BaO":  -2.82,
    "TiO2": -3.30,
    "ZrO2": -3.74,
}
DEMO_SYSTEM = ("Ca", "Sr", "Ba", "Ti", "Zr", "O")
