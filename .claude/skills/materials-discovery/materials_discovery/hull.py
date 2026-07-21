"""
hull.py -- convex-hull / energy-above-hull, the GNoME stability criterion.

A candidate crystal is "stable" (a real discovery) iff it sits ON the convex hull
of formation energy vs composition -- i.e. it does not lower its energy by
decomposing into any mixture of known phases. GNoME's GNN predicts this number
(energy above hull, E_hull); we compute it exactly by linear programming given a
set of reference phases.

Convention (Materials Project standard): energies are FORMATION energy per atom
(eV/atom); the pure elements are the zero reference. E_hull(x) =
  eps(x) - min over phase mixtures matching composition x of sum_i a_i eps_i,
solved as an LP over atom-fraction weights a_i >= 0, sum a_i = 1, sum a_i f_i = x.

Pure/deterministic/offline (numpy + scipy.optimize.linprog).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import linprog
from .composition import Composition


@dataclass(frozen=True)
class Entry:
    composition: Composition
    formation_energy_per_atom: float   # eV/atom, relative to the elements
    label: str = ""


class PhaseDiagram:
    """Convex hull over a set of reference Entry objects."""

    def __init__(self, entries: list):
        self.entries = list(entries)
        # universe of elements
        els = set()
        for e in self.entries:
            els.update(e.composition.elements)
        # ensure every element has a zero-energy reference (the elemental phase)
        have = {e.composition.reduced_formula() for e in self.entries}
        for el in sorted(els):
            if el not in have:
                self.entries.append(Entry(Composition.from_formula(el), 0.0, label=el))
        # rebuild element list (stable order)
        self.elements = sorted({el for e in self.entries
                                for el in e.composition.elements})
        self._eps = np.array([e.formation_energy_per_atom for e in self.entries])
        # fractional-composition matrix F[i, j] = fraction of element j in phase i
        F = np.zeros((len(self.entries), len(self.elements)))
        for i, e in enumerate(self.entries):
            frac = e.composition.fractional()
            for j, el in enumerate(self.elements):
                F[i, j] = frac.get(el, 0.0)
        self._F = F

    def _augmented_system(self, comp: Composition):
        """Return (eps, F, entry_refs, elements) with an elemental phase guaranteed
        for every element in `comp` (so the LP is always feasible). A query element
        absent from all references gets an elemental reference at energy 0."""
        elements = list(self.elements)
        entries = list(self.entries)
        eps = list(self._eps)
        Frows = [self._F[i].tolist() for i in range(len(self.entries))]
        for el in comp.elements:
            if el not in elements:
                elements.append(el)
                for row in Frows:              # pad existing rows with 0 for new column
                    row.append(0.0)
                # add the elemental phase (energy 0, pure element)
                entries.append(Entry(Composition.from_formula(el), 0.0, label=el))
                eps.append(0.0)
                Frows.append([1.0 if e == el else 0.0 for e in elements])
        return np.array(eps), np.array(Frows), entries, elements

    def hull_energy(self, comp: Composition, exclude_label: str | None = None) -> float:
        """Minimum formation energy/atom reachable by decomposing into references."""
        eps, F, entries, elements = self._augmented_system(comp)
        x = np.array([comp.fractional().get(el, 0.0) for el in elements])
        mask = np.ones(len(entries), dtype=bool)
        if exclude_label is not None:
            for i, e in enumerate(entries):
                if e.label == exclude_label or e.composition.reduced_formula() == exclude_label:
                    mask[i] = False
        c = eps[mask]
        Fm = F[mask]
        n = c.shape[0]
        # equality: F^T a = x  (per element)  and  sum a = 1
        A_eq = np.vstack([Fm.T, np.ones((1, n))])
        b_eq = np.concatenate([x, [1.0]])
        res = linprog(c=c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * n, method="highs")
        if not res.success:
            raise RuntimeError(f"hull LP failed for {comp}: {res.message}")
        return float(res.fun)

    def e_above_hull(self, comp: Composition, formation_energy_per_atom: float,
                     exclude_self: bool = True) -> float:
        """E_hull >= 0 always for a hull built from references; < ~0 means the
        candidate DEFINES a new hull point (a predicted-stable discovery)."""
        label = comp.reduced_formula()
        hull = self.hull_energy(comp, exclude_label=label if exclude_self else None)
        return formation_energy_per_atom - hull

    def decomposition(self, comp: Composition, exclude_label: str | None = None) -> list:
        """Return the stable phases (and weights) the composition decomposes into."""
        eps, F, entries, elements = self._augmented_system(comp)
        x = np.array([comp.fractional().get(el, 0.0) for el in elements])
        mask = np.ones(len(entries), dtype=bool)
        if exclude_label is not None:
            for i, e in enumerate(entries):
                if e.label == exclude_label or e.composition.reduced_formula() == exclude_label:
                    mask[i] = False
        idx = np.where(mask)[0]
        c = eps[mask]
        Fm = F[mask]
        n = c.shape[0]
        A_eq = np.vstack([Fm.T, np.ones((1, n))])
        b_eq = np.concatenate([x, [1.0]])
        res = linprog(c=c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * n, method="highs")
        out = []
        for k, w in enumerate(res.x):
            if w > 1e-6:
                out.append((entries[idx[k]].composition.reduced_formula(), round(float(w), 4)))
        return sorted(out, key=lambda t: -t[1])
