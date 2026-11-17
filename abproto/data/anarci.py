r"""Exact CDR / framework annotation via ANARCI (optimisation #7).

ANARCI (Dunbar & Deane, 2016) assigns IMGT numbering to a VH/VL domain, which
gives exact CDR boundaries. We use it when importable; otherwise we fall back to
an anchor-based locator that finds CDR3 reliably (between the conserved pre-CDR3
cysteine and the J-motif WGxG/FGxG) and CDR1/2 by approximate IMGT offsets.

    pip install anarci        # also needs HMMER (hmmscan) + germline HMMs

CDR boundaries let the splitter dedup and measure novelty on the CDRs (what the
paper clusters on) rather than whole chains.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

try:
    import anarci as _anarci_mod
    _HAS_ANARCI = True
except Exception:          # pragma: no cover - import is environment dependent
    _HAS_ANARCI = False

# IMGT CDR position ranges (inclusive) used to slice ANARCI numbering.
_IMGT_CDR = {"1": (27, 38), "2": (56, 65), "3": (105, 117)}

# Fallback anchors.
_J_HEAVY = re.compile(r"WG.G")
_J_LIGHT = re.compile(r"FG.G")


@dataclass
class ChainAnnotation:
    chain_type: str                 # "H" or "L"
    cdr_spans: dict[str, tuple[int, int]]   # {"1":(lo,hi), "2":..., "3":...} 0-based
    numbered: bool                  # True if from ANARCI, False if heuristic


def annotate_chain(seq: str, prefer: str = "H") -> Optional[ChainAnnotation]:
    """Return CDR spans for one antibody V-domain sequence."""
    if _HAS_ANARCI:
        a = _annotate_anarci(seq)
        if a is not None:
            return a
    return _annotate_heuristic(seq, prefer)


def _annotate_anarci(seq: str) -> Optional[ChainAnnotation]:
    try:
        numbering, chain_type = _anarci_mod.number(seq)
    except Exception:
        return None
    if not numbering:
        return None
    ct = "H" if chain_type == "H" else "L"
    # numbering: list of ((imgt_pos, ins), aa); map to sequence indices.
    spans: dict[str, tuple[int, int]] = {}
    seq_idx = []
    for (pos, _ins), aa in numbering:
        if aa != "-":
            seq_idx.append((pos, len(seq_idx)))
    for cdr, (lo, hi) in _IMGT_CDR.items():
        idxs = [si for pos, si in seq_idx if lo <= pos <= hi]
        if idxs:
            spans[cdr] = (min(idxs), max(idxs) + 1)
    return ChainAnnotation(ct, spans, numbered=True)


def _annotate_heuristic(seq: str, prefer: str) -> Optional[ChainAnnotation]:
    """Anchor-based fallback: CDR3 between the last pre-J cysteine and WGxG/FGxG;
    CDR1/2 by approximate offsets. Approximate — install ANARCI for exact spans."""
    jmotif = _J_HEAVY if prefer == "H" else _J_LIGHT
    m = jmotif.search(seq)
    spans: dict[str, tuple[int, int]] = {}
    if m:
        j = m.start()
        cys = seq.rfind("C", 0, j)
        if cys != -1 and j - cys > 2:
            spans["3"] = (cys + 1, j)
    # CDR1/2 ~ IMGT-canonical windows relative to the domain start.
    if len(seq) > 70:
        spans["1"] = (25, min(34, len(seq)))
        spans["2"] = (50, min(58, len(seq)))
    if not spans:
        return None
    return ChainAnnotation(prefer, spans, numbered=False)


def cdr_sequences(seq: str, ann: ChainAnnotation) -> dict[str, str]:
    return {c: seq[lo:hi] for c, (lo, hi) in ann.cdr_spans.items()}


def concat_cdrs(heavy: str, light: str) -> str:
    """Concatenated H+L CDRs (for CDR-level dedup / novelty). Empty if none found."""
    out = []
    for s, pref in ((heavy, "H"), (light, "L")):
        if not s:
            continue
        ann = annotate_chain(s, prefer=pref)
        if ann:
            out.extend(cdr_sequences(s, ann).values())
    return "".join(out)


def has_anarci() -> bool:
    return _HAS_ANARCI
