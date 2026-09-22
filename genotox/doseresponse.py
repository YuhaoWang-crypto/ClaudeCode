"""
Dose-response analysis and the assay's decision rule.

Deliberately separate from the biology: the induction ratio, the 1.5
threshold and the growth-factor gate are *operational definitions* from the
ISO 13829 / Oda-style umu protocol, not biological constants.  Keeping them
here means the core can be refitted without touching the call logic, and the
call logic can be changed without anyone suspecting the biology moved.
"""
from __future__ import annotations

import numpy as np

#: Induction-ratio threshold for a positive call (operational, ISO 13829).
IR_THRESHOLD = 1.5
#: Growth factor below which a well is not interpretable (operational).
GROWTH_GATE = 0.5


def dose_series(assay, compound, doses, s9=False) -> dict:
    """Run a solvent control plus a dose series; return per-well rows.

    The induction ratio is formed from growth-normalised specific activity
    (``units``), never from raw A410, and every row carries its own growth
    factor so that a failed validity gate travels with the number it
    invalidates.
    """
    ctrl = assay.control(compound, s9=s9)
    u0 = ctrl["units"]
    n0 = ctrl["A600"]

    rows = []
    for d in doses:
        w = assay.well(compound, float(d), s9=s9)
        ir = w["units"] / u0 if u0 > 0 else np.nan
        gf = w["A600"] / n0 if n0 > 0 else np.nan
        rows.append({
            "dose_uM": float(d),
            "IR": ir,
            "growth_factor": gf,
            "valid": bool(gf >= GROWTH_GATE),
            "units": w["units"],
            "A410": w["A410"],
            "lesions_end": w["lesions_end"],
            "viability_end": w["viability_end"],
        })
    return {"compound": compound.name, "s9": s9, "control_units": u0,
            "rows": rows}


def ec_ir(series: dict, threshold: float = IR_THRESHOLD):
    """Lowest dose at which IR crosses ``threshold``, among valid wells.

    Linear interpolation between the bracketing doses.  Returns ``None`` when
    the threshold is never crossed inside the valid window — which is a
    meaningfully different outcome from "not crossed at all", so callers
    should also look at :func:`call_result`.
    """
    rows = [r for r in series["rows"] if r["valid"] and np.isfinite(r["IR"])]
    for a, b in zip(rows, rows[1:]):
        if a["IR"] < threshold <= b["IR"]:
            span = b["IR"] - a["IR"]
            if span <= 0:
                return b["dose_uM"]
            f = (threshold - a["IR"]) / span
            return a["dose_uM"] + f * (b["dose_uM"] - a["dose_uM"])
    return None


def call_result(series: dict, threshold: float = IR_THRESHOLD) -> dict:
    """The assay's verdict, with the reason attached.

    Three outcomes, not two.  "Inconclusive" is a real result here: if every
    well that reached the threshold also failed the growth gate, the honest
    report is that the test could not decide, and collapsing that into
    "negative" is how cytotoxic compounds get mis-scored in practice.
    """
    rows = series["rows"]
    valid = [r for r in rows if r["valid"]]
    max_ir_valid = max((r["IR"] for r in valid), default=float("nan"))
    max_ir_any = max((r["IR"] for r in rows), default=float("nan"))

    if valid and max_ir_valid >= threshold:
        verdict, why = "POSITIVE", "IR >= threshold within the valid window"
    elif not valid:
        verdict, why = "INCONCLUSIVE", "no well passed the growth gate"
    elif max_ir_any >= threshold:
        verdict, why = ("INCONCLUSIVE",
                        "threshold only reached in growth-gated wells")
    else:
        verdict, why = "NEGATIVE", "IR below threshold across the valid window"

    return {
        "verdict": verdict,
        "reason": why,
        "max_IR_valid": max_ir_valid,
        "max_IR_any": max_ir_any,
        "ec_ir": ec_ir(series, threshold),
        "n_valid": len(valid),
        "n_wells": len(rows),
    }


def log_doses(lo: float, hi: float, n: int = 13):
    return np.concatenate([[0.0], np.geomspace(lo, hi, n)])
