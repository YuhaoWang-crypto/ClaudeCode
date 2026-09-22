"""
Dose-response analysis and the assay's decision rule.

Deliberately separate from the biology: induction thresholds and validity
gates are *operational definitions* from a protocol, not biological
constants.  Keeping them here means a core can be refitted without touching
the call logic, and the call logic can be changed without anyone suspecting
the biology moved.

Each protocol supplies its own numbers; the machinery below is shared, and
reads only the generic ``signal`` / ``biomass`` contract that the readout
layer guarantees.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Protocol:
    """Operational decision rule for one assay format."""

    name: str
    ir_threshold: float
    growth_gate: float
    gate_label: str = "growth factor"
    source: str = ""


#: umu test, ISO 13829-style: IR >= 1.5, wells with growth factor < 0.5
#: discarded.
UMU = Protocol("umu", 1.5, 0.50, "growth factor",
               "ISO 13829 / Oda-style operational criteria")

#: Mammalian GADD45a-GFP reporter line.  A lower induction threshold with a
#: stricter cytotoxicity limit is the usual trade for the higher background
#: variability of a fluorescent reporter in cycling cells.  Numbers here are
#: representative of that format rather than any one vendor's criteria.
GADD45A_GFP = Protocol("GADD45a-GFP", 1.3, 0.80, "relative cell density",
                       "representative reporter-line criteria (H)")

# Backwards-compatible aliases for callers written against the umu-only API.
IR_THRESHOLD = UMU.ir_threshold
GROWTH_GATE = UMU.growth_gate


def dose_series(assay, compound, doses, s9=False, protocol: Protocol = UMU) -> dict:
    """Run a solvent control plus a dose series; return per-well rows.

    The induction ratio is formed from the biomass-normalised ``signal``,
    never from a raw instrument reading, and every row carries its own growth
    factor so that a failed validity gate travels with the number it
    invalidates.
    """
    ctrl = assay.control(compound, s9=s9)
    s0, b0 = ctrl["signal"], ctrl["biomass"]

    rows = []
    for d in doses:
        w = assay.well(compound, float(d), s9=s9)
        ir = w["signal"] / s0 if s0 > 0 else np.nan
        gf = w["biomass"] / b0 if b0 > 0 else np.nan
        row = {
            "dose_uM": float(d),
            "IR": ir,
            "growth_factor": gf,
            "valid": bool(gf >= protocol.growth_gate),
            "signal": w["signal"],
            "lesions_end": w.get("lesions_end", float("nan")),
            "viability_end": w.get("viability_end", float("nan")),
        }
        for k in ("p53_peak", "p53_auc", "p53_pulses", "A410",
                  "fluorescence_total"):
            if k in w:
                row[k] = w[k]
        rows.append(row)
    return {"compound": compound.name, "s9": s9, "protocol": protocol,
            "control_signal": s0, "rows": rows}


def ec_ir(series: dict, threshold: float | None = None):
    """Lowest dose at which IR crosses ``threshold``, among valid wells.

    Linear interpolation between the bracketing doses.  Returns ``None`` when
    the threshold is never crossed inside the valid window — which is a
    meaningfully different outcome from "not crossed at all", so callers
    should also look at :func:`call_result`.
    """
    if threshold is None:
        threshold = series["protocol"].ir_threshold
    rows = [r for r in series["rows"] if r["valid"] and np.isfinite(r["IR"])]
    for a, b in zip(rows, rows[1:]):
        if a["IR"] < threshold <= b["IR"]:
            span = b["IR"] - a["IR"]
            if span <= 0:
                return b["dose_uM"]
            f = (threshold - a["IR"]) / span
            return a["dose_uM"] + f * (b["dose_uM"] - a["dose_uM"])
    return None


def call_result(series: dict, threshold: float | None = None) -> dict:
    """The assay's verdict, with the reason attached.

    Three outcomes, not two.  "Inconclusive" is a real result here: if every
    well that reached the threshold also failed the validity gate, the honest
    report is that the test could not decide, and collapsing that into
    "negative" is how cytotoxic compounds get mis-scored in practice.
    """
    proto = series["protocol"]
    if threshold is None:
        threshold = proto.ir_threshold
    rows = series["rows"]
    valid = [r for r in rows if r["valid"]]
    max_ir_valid = max((r["IR"] for r in valid), default=float("nan"))
    max_ir_any = max((r["IR"] for r in rows), default=float("nan"))

    if valid and max_ir_valid >= threshold:
        verdict, why = "POSITIVE", "IR >= threshold within the valid window"
    elif not valid:
        verdict, why = "INCONCLUSIVE", f"no well passed the {proto.gate_label} gate"
    elif max_ir_any >= threshold:
        verdict, why = ("INCONCLUSIVE",
                        "threshold only reached in gated wells")
    else:
        verdict, why = "NEGATIVE", "IR below threshold across the valid window"

    return {
        "verdict": verdict,
        "reason": why,
        "protocol": proto.name,
        "max_IR_valid": max_ir_valid,
        "max_IR_any": max_ir_any,
        "ec_ir": ec_ir(series, threshold),
        "n_valid": len(valid),
        "n_wells": len(rows),
    }


def log_doses(lo: float, hi: float, n: int = 13):
    return np.concatenate([[0.0], np.geomspace(lo, hi, n)])
