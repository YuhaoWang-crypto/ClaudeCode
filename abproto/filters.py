"""Physics + sequence-liability filtering (plan step 5, and the hard
constraints of step 3).

Sequence liabilities are detected with the standard motif rules used in
developability triage. The physics term here is a C-alpha clash proxy; a
production system swaps in a real all-atom interface energy.
"""

from __future__ import annotations

import re

from .geometry import self_clash_count
from .types import CDRS, Design

# --- sequence liability motifs (scanned within CDRs only) -----------------
# N-linked glycosylation sequon: N-X-[S/T], X != P
_NGLYC = re.compile(r"N[^P][ST]")
# Asn deamidation hotspots.
_DEAMIDATION = re.compile(r"N[GSTN]")
# Asp isomerization hotspots.
_ISOMERIZATION = re.compile(r"D[GSTD]")
# Met oxidation / free-thiol handled separately.
_UNPAIRED_CYS = "C"


def cdr_liabilities(design: Design) -> dict[str, int]:
    """Count liability motifs across all six CDRs."""
    counts = {"nglyc": 0, "deamidation": 0, "isomerization": 0, "cys": 0}
    for cdr in CDRS:
        if cdr not in design.cdr_spans:
            continue
        seq = design.cdr_sequence(cdr)
        counts["nglyc"] += len(_NGLYC.findall(seq))
        counts["deamidation"] += len(_DEAMIDATION.findall(seq))
        counts["isomerization"] += len(_ISOMERIZATION.findall(seq))
        counts["cys"] += seq.count(_UNPAIRED_CYS)
    return counts


def clash_count(design: Design) -> int:
    """Steric-clash proxy over the antibody C-alpha trace."""
    return self_clash_count(design.antibody_ca(), contact=3.0)


# --- thresholds -----------------------------------------------------------
DEFAULT_THRESHOLDS = {
    "max_recovery_rmsd": 2.0,     # angstrom (self-consistency gate)
    "max_nglyc": 0,               # hard: no glycosylation sequons in CDRs
    "max_deamidation": 2,
    "max_isomerization": 2,
    "max_cys": 0,                 # no unpaired cysteines in CDRs
    "max_clashes": 5,
}


def passes(recovery_rmsd: float, liabilities: dict[str, int], clashes: int,
           thresholds: dict | None = None) -> bool:
    t = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    return (
        recovery_rmsd <= t["max_recovery_rmsd"]
        and liabilities["nglyc"] <= t["max_nglyc"]
        and liabilities["deamidation"] <= t["max_deamidation"]
        and liabilities["isomerization"] <= t["max_isomerization"]
        and liabilities["cys"] <= t["max_cys"]
        and clashes <= t["max_clashes"]
    )
