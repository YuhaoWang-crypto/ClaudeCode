"""
composition.py -- minimal chemical-composition algebra (no pymatgen dependency).

Parses formulas like "Li3PS4", "LiLaTiO6", "Al2O3" into element->amount maps and
provides the fractional-composition vector used by the convex-hull code.
Pure, deterministic, offline.
"""
from __future__ import annotations
import re
from dataclasses import dataclass

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")

# Pauling electronegativities -- shared constant (score.py imports this too).
# Used to order reduced formulas electropositive-first (SrTiO3, not O3SrTi).
ELECTRONEG = {
    "Li": 0.98, "Na": 0.93, "K": 0.82, "Ag": 1.93,
    "Mg": 1.31, "Ca": 1.00, "Sr": 0.95, "Ba": 0.89, "Zn": 1.65,
    "Al": 1.61, "Ga": 1.81, "Sc": 1.36, "Y": 1.22, "La": 1.10, "In": 1.78,
    "Ti": 1.54, "Zr": 1.33, "Sn": 1.96, "Ge": 2.01, "Si": 1.90, "Hf": 1.30,
    "Nb": 1.60, "Ta": 1.50, "V": 1.63, "P": 2.19, "W": 2.36, "Mo": 2.16,
    "O": 3.44, "S": 2.58, "Se": 2.55, "F": 3.98, "Cl": 3.16, "Br": 2.96,
    "I": 2.66, "N": 3.04,
}


def _formula_order(el: str):
    """Sort key: electropositive (low electronegativity) first, then alphabetical."""
    return (ELECTRONEG.get(el, 2.2), el)


def parse_formula(formula: str) -> dict:
    """'Li3PS4' -> {'Li':3.0,'P':1.0,'S':4.0}. Rejects unparseable strings."""
    formula = formula.strip()
    if not formula:
        raise ValueError("empty formula")
    counts: dict = {}
    consumed = 0
    for m in _TOKEN.finditer(formula):
        el, num = m.group(1), m.group(2)
        counts[el] = counts.get(el, 0.0) + (float(num) if num else 1.0)
        consumed += len(m.group(0))
    if consumed != len(formula):
        raise ValueError(f"could not fully parse formula: {formula!r}")
    return counts


@dataclass(frozen=True)
class Composition:
    """An immutable composition; amounts are per-formula-unit (not normalized)."""
    amounts: tuple  # tuple of (element, amount) sorted by element

    @classmethod
    def from_formula(cls, formula: str) -> "Composition":
        return cls.from_dict(parse_formula(formula))

    @classmethod
    def from_dict(cls, d: dict) -> "Composition":
        items = tuple(sorted((el, float(a)) for el, a in d.items() if a > 0))
        if not items:
            raise ValueError("composition has no positive amounts")
        return cls(items)

    @property
    def elements(self) -> tuple:
        return tuple(el for el, _ in self.amounts)

    def as_dict(self) -> dict:
        return {el: a for el, a in self.amounts}

    @property
    def num_atoms(self) -> float:
        return sum(a for _, a in self.amounts)

    def fractional(self) -> dict:
        n = self.num_atoms
        return {el: a / n for el, a in self.amounts}

    def reduced_formula(self) -> str:
        """Divide by gcd-like common factor for display (integer amounts only)."""
        ordered = sorted(self.amounts, key=lambda t: _formula_order(t[0]))
        amounts = [a for _, a in ordered]
        if all(abs(a - round(a)) < 1e-9 for a in amounts):
            ints = [int(round(a)) for a in amounts]
            from math import gcd
            g = 0
            for x in ints:
                g = gcd(g, x)
            g = g or 1
            parts = []
            for (el, _), x in zip(ordered, ints):
                q = x // g
                parts.append(el if q == 1 else f"{el}{q}")
            return "".join(parts)
        return "".join(f"{el}{a:g}" for el, a in ordered)

    def __str__(self) -> str:
        return self.reduced_formula()
