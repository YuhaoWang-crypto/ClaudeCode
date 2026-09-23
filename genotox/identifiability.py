"""
Diagnostics: what can the assay's own readout actually constrain?

The package's stated next step has been "fit a core to a measured dose
series".  The honest step *before* that is to ask how many parameters such a
fit could recover at all — because a least-squares routine will happily
return a confident value for a parameter the data cannot see.

Two distinct questions, answered separately:

**Structural identifiability** — is a parameter invisible by construction, no
matter how good the data?  The umu readout is an induction *ratio*, so any
parameter that scales the reporter linearly cancels exactly: transcription
rate, translation rate and the instrument gain all have identically zero
influence.  Worse, a model can carry a continuous symmetry that is not
obvious from the equations; :func:`check_symmetry` tests a candidate
direction by applying it at finite size rather than trusting the
linearisation.

**Practical identifiability** — given a real noise floor, how many
*directions* in parameter space are estimable to a useful precision?  This is
a property of the sensitivity spectrum, not of any single parameter, and it
is almost always far smaller than the parameter count.  Reporting it as a
count of directions rather than a condition number is deliberate: a condition
number spanning many orders of magnitude says a model is sloppy without
saying what could be measured.

Everything here is local (a linearisation about one operating point) and says
nothing about whether the model is *right* — only about what a fit to this
readout could and could not learn.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np


@dataclass
class Sensitivity:
    """Relative sensitivities d log(observable) / d log(parameter)."""

    matrix: np.ndarray          # (n_observations, n_parameters)
    names: tuple                # parameter names, column order
    label: str = ""

    @property
    def spectrum(self) -> np.ndarray:
        return np.linalg.svd(self.matrix, compute_uv=False)

    def column_norms(self) -> dict:
        return {n: float(np.linalg.norm(self.matrix[:, j]))
                for j, n in enumerate(self.names)}

    def structural_nulls(self, tol: float = 1e-6) -> list:
        """Directions the observable cannot see at all.

        ``tol`` is relative to the leading singular value.  A direction below
        it is flat to within integration error, which for a well-converged
        solver means flat, full stop — but :func:`check_symmetry` is what
        turns that suspicion into a demonstration.
        """
        U, sv, Vt = np.linalg.svd(self.matrix, full_matrices=False)
        cut = sv[0] * tol
        out = []
        for k in range(len(sv)):
            if sv[k] <= cut:
                v = Vt[k]
                loads = sorted(zip(self.names, v), key=lambda t: -abs(t[1]))
                out.append({"singular_value": float(sv[k]),
                            "loadings": [(n, float(c)) for n, c in loads
                                         if abs(c) > 0.05]})
        return out

    def practical_rank(self, noise: float = 0.05,
                       factor: float = 1.65) -> int:
        """Number of directions estimable to within ``factor`` at ``noise``.

        With relative sensitivities and relative measurement noise ``noise``,
        the standard-error along the direction with singular value ``s`` is
        ``noise / s`` in log-parameter units.  A direction is counted when
        that is smaller than ``log(factor)`` — i.e. the parameter combination
        is pinned to better than a factor of ``factor``.
        """
        sv = self.spectrum
        return int((sv > noise / np.log(factor)).sum())


def relative_sensitivity(make, observe, params: dict,
                         rel_step: float = 0.02, label: str = "") -> Sensitivity:
    """Central-difference sensitivity of log(observable) to log(parameter).

    ``make(**overrides)`` builds the object under test, ``observe(obj)``
    returns a positive observation vector.  Working in logs on both sides
    makes the columns dimensionless and therefore comparable, which a raw
    Fisher matrix over parameters with different units is not.
    """
    cols, names = [], []
    for name, value in params.items():
        if value == 0:
            continue                       # a log step off zero is undefined
        up = np.log(observe(make(**{name: value * (1 + rel_step)})))
        dn = np.log(observe(make(**{name: value * (1 - rel_step)})))
        cols.append((up - dn) / (2 * rel_step))
        names.append(name)
    return Sensitivity(np.stack(cols, axis=1), tuple(names), label)


def check_symmetry(make, observe, direction: dict,
                   scales=(0.5, 1.5, 4.0)) -> dict:
    """Apply a candidate null direction at finite size and measure the drift.

    A linearisation can only ever say "flat to first order". Scaling the
    parameters by a large factor and finding the observable unmoved is the
    difference between a small singular value and an exact symmetry — and an
    exact symmetry means one parameter can be fixed by convention rather than
    fitted.
    """
    y0 = np.asarray(observe(make()), float)
    worst = 0.0
    rows = []
    for c in scales:
        over = {n: v0 * (c ** e) for n, (v0, e) in direction.items()}
        y = np.asarray(observe(make(**over)), float)
        drift = float(np.max(np.abs(y - y0) / np.abs(y0)))
        rows.append({"scale": c, "max_rel_drift": drift})
        worst = max(worst, drift)
    return {"per_scale": rows, "max_rel_drift": worst,
            "exact": worst < 1e-6}


def summarise(sens: Sensitivity, noise: float = 0.05,
              factor: float = 1.65) -> dict:
    sv = sens.spectrum
    nulls = sens.structural_nulls()
    return {
        "label": sens.label,
        "n_observations": int(sens.matrix.shape[0]),
        "n_parameters": int(sens.matrix.shape[1]),
        "spectrum": sv,
        "practical_rank": sens.practical_rank(noise, factor),
        "n_structural_nulls": len(nulls),
        "structural_nulls": nulls,
        "column_norms": sens.column_norms(),
        "noise": noise, "factor": factor,
    }


def print_summary(s: dict) -> None:
    print(f"\n{s['label']}")
    print(f"  observations {s['n_observations']}   parameters "
          f"{s['n_parameters']}")
    sv = s["spectrum"]
    shown = np.array2string(sv[:8], precision=3, suppress_small=True)
    print(f"  singular values (first 8): {shown}")
    print(f"  estimable directions at {100*s['noise']:.0f}% noise, to within "
          f"a factor of {s['factor']}: "
          f"{s['practical_rank']} of {s['n_parameters']}")
    dead = [n for n, v in s["column_norms"].items() if v < 1e-9]
    if dead:
        print(f"  zero influence on the readout: {', '.join(dead)}")
    for nul in s["structural_nulls"]:
        loads = "  ".join(f"{n}:{c:+.2f}" for n, c in nul["loadings"])
        print(f"  null direction (sv={nul['singular_value']:.2e}): {loads}")
