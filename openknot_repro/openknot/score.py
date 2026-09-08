"""OpenKnot score: Eterna Classic Score (ECS) + Crossed Pair Quality (CPQ).

Faithful re-implementation of the reference pipeline
(eternagame/OpenKnotScorePipeline, src/openknotscore/pipeline/scoring.py),
with the base-pair utilities that the reference imports from arnie inlined so
this module has no dependency beyond the standard library.

Reference semantics that are easy to get wrong and are preserved here:

* ECS unpaired threshold is 0.25 * 0.5 + 0.75 * 0.0 = 0.125; ECS paired
  threshold is 0.5.  The dataset README states 0.25 for paired residues, but
  0.25 is the CPQ threshold, not the ECS one -- the code is authoritative.
* Positions with missing (NaN) reactivity are counted as *incorrect* in ECS:
  the reference skips them in the loop but keeps them in the denominator
  ``len(correct_hit)``.
* CPQ enumerates crossed residues from the *unfiltered* structure, then gives
  a residue full weight only if its pair survives restriction to the scored
  region, and half weight otherwise.  NaN crossed residues count towards the
  denominator ``max_count`` but never towards the numerator.
* OKS = 0.5 * ECS + 0.5 * CPQ_quality (the second element of the reference's
  CPQ return value, not the first).
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

OPEN2CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
CLOSE2OPEN = {v: k for k, v in OPEN2CLOSE.items()}

ECS_PAIRED_THRESHOLD = 0.5
ECS_UNPAIRED_THRESHOLD = 0.25 * ECS_PAIRED_THRESHOLD
CPQ_THRESHOLD = 0.25


def _is_missing(value) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def bp_list(structure: str) -> list[tuple[int, int]]:
    """Dot-bracket to (i, j) base-pair list, pseudoknots included.

    Matches arnie's ``convert_dotbracket_to_bp_list(allow_pseudoknots=True)``:
    one LIFO stack per bracket type, the result sorted by opening index.  That
    ordering matters downstream: ``get_helices`` walks the list and only stacks
    a pair onto the previous one, so an unsorted list would shatter every helix.

    Unbalanced closing brackets are dropped rather than raising: a handful of
    released target structures carry one.
    """
    stacks: dict[str, list[int]] = {k: [] for k in OPEN2CLOSE}
    pairs: list[tuple[int, int]] = []
    for i, char in enumerate(structure):
        if char in OPEN2CLOSE:
            stacks[char].append(i)
        elif char in CLOSE2OPEN:
            stack = stacks[CLOSE2OPEN[char]]
            if stack:
                pairs.append((stack.pop(), i))
    return sorted(pairs)


def get_helices(structure: str, allowed_bulge_len: int = 0) -> list[list[tuple[int, int]]]:
    """Group base pairs into helices, port of arnie's ``get_helices``.

    The grouping walks the base-pair list (sorted by opening index) and starts
    a new helix whenever a pair does not stack onto the previous one.
    """
    pairs = bp_list(structure)
    helices: list[list[tuple[int, int]]] = []
    current: list[tuple[int, int]] = []
    for pair in pairs:
        if not current:
            current = [pair]
            continue
        left = range(current[-1][0] + 1, current[-1][0] + allowed_bulge_len + 2)
        right = range(current[-1][1] - allowed_bulge_len - 1, current[-1][1])
        if pair[0] in left and pair[1] in right:
            current.append(pair)
        else:
            helices.append(current)
            current = [pair]
    helices.append(current)
    return [h for h in helices if h]


def crossing_residues(pairs: Iterable[tuple[int, int]]) -> set[int]:
    """Residues belonging to at least one crossed (pseudoknotted) pair.

    A pair (i, j) is crossed when some other pair (m, n) satisfies
    i < m < j < n or m < i < n < j.
    """
    pairs = list(pairs)
    crossed: set[int] = set()
    for a in pairs:
        for b in pairs:
            if (a[0] < b[0] < a[1] < b[1]) or (b[0] < a[0] < b[1] < a[1]):
                crossed.update((a[0], a[1], b[0], b[1]))
    return crossed


def filter_singlet_pairs(
    structure: str, min_len_helix: int = 2, allowed_bulge_len: int = 0
) -> list[tuple[int, int]]:
    """Base pairs surviving arnie's ``post_process_struct``: helices shorter than
    ``min_len_helix`` are dropped.  Returned sorted."""
    kept: list[tuple[int, int]] = []
    for helix in get_helices(structure, allowed_bulge_len):
        if len(helix) >= min_len_helix:
            kept.extend(helix)
    return sorted(kept)


def eterna_classic_score(
    structure: str,
    data: Sequence[float],
    score_start_idx: int,
    score_end_idx: int,
    tol: float = 0.0,
) -> float:
    """Percentage of scored residues whose reactivity matches the structure.

    ``score_start_idx``/``score_end_idx`` are inclusive 0-based indices into
    ``structure`` and ``data``, which must be the same length.

    ``tol`` widens both comparisons by an absolute amount.  It exists because
    the released reactivities are rounded to three decimals, so a residue whose
    true value sat just below a threshold can land exactly on it; the reference
    scored the unrounded values.  Use 0.0 for the literal reference rule.
    """
    if len(structure) != len(data):
        raise ValueError(f"structure/data length mismatch: {len(structure)} != {len(data)}")
    total = score_end_idx - score_start_idx + 1
    if total <= 0:
        raise ValueError("empty scoring region")
    hits = 0
    for i in range(score_start_idx, score_end_idx + 1):
        value = data[i]
        if _is_missing(value):
            continue  # counted as a miss: it stays in the denominator
        if structure[i] == ".":
            if value > ECS_UNPAIRED_THRESHOLD - tol:
                hits += 1
        elif value < ECS_PAIRED_THRESHOLD + tol:
            hits += 1
    return 100.0 * hits / total


def crossed_pair_scores(
    structure: str,
    data: Sequence[float],
    score_start_idx: int,
    score_end_idx: int,
    tol: float = 0.0,
    filter_singlets: bool = True,
) -> tuple[float, float]:
    """Return (crossed_pair_score, crossed_pair_quality_score).

    Only the second is used by the OpenKnot score.

    ``filter_singlets`` drops helices of length 1 before crossed pairs are
    enumerated.  It defaults to True because that is what reproduces the
    released ``target_openknot_score`` column; the reference function's own
    default is False, and the released scores were evidently produced with the
    flag turned on for this half of the score only (the ECS half is scored on
    the unfiltered structure).
    """
    if len(structure) != len(data):
        raise ValueError(f"structure/data length mismatch: {len(structure)} != {len(data)}")
    pairs = filter_singlet_pairs(structure) if filter_singlets else bp_list(structure)
    crossed = crossing_residues(pairs)
    in_region = [
        bp
        for bp in pairs
        if bp[0] >= score_start_idx and bp[1] >= score_start_idx
        and bp[0] <= score_end_idx and bp[1] <= score_end_idx
    ]
    crossed_in_region = crossing_residues(in_region)

    numerator = 0.0
    max_count = 0
    for i in sorted(crossed):
        if i < score_start_idx or i > score_end_idx:
            continue
        max_count += 1
        value = data[i]
        if _is_missing(value):
            continue
        if value < CPQ_THRESHOLD + tol:
            numerator += 1.0 if i in crossed_in_region else 0.5

    region_length = 1 + score_end_idx - score_start_idx
    max_crossed = 0.7 * max(region_length - 20, 20)
    cp_score = 100.0 * min(numerator / max_crossed, 1.0)
    cpq_score = 100.0 * numerator / max_count if max_count > 0 else 0.0
    return cp_score, cpq_score


def openknot_score(
    structure: str,
    data: Sequence[float],
    score_start_idx: int,
    score_end_idx: int,
    tol: float = 0.0,
    filter_singlets: bool = True,
) -> float:
    """OKS = 0.5 * ECS + 0.5 * crossed-pair-quality."""
    ecs = eterna_classic_score(structure, data, score_start_idx, score_end_idx, tol)
    _, cpq = crossed_pair_scores(
        structure, data, score_start_idx, score_end_idx, tol, filter_singlets
    )
    return 0.5 * ecs + 0.5 * cpq


def embed_target(target_structure: str, sub_start: int, sub_end: int, full_length: int) -> str:
    """Place a design-length target structure into a full-length dot-bracket string.

    ``sub_start``/``sub_end`` are the 1-based inclusive bounds of the design
    within the padded construct, as released in the benchmark CSV.  Flanking
    pad positions are unpaired, which is what the constructs are built for.
    """
    start = sub_start - 1
    end = sub_end - 1
    if len(target_structure) != end - start + 1:
        raise ValueError(
            f"target length {len(target_structure)} does not match sub_start/sub_end span "
            f"{end - start + 1}"
        )
    return "." * start + target_structure + "." * (full_length - end - 1)
