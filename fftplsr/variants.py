"""Mutation-string parsing, sequence construction and combinatorial space enumeration.

Variant naming follows the convention used in the FP4COM paper data:

    "D2N"              single substitution, 1-based position
    "D2N/K3N"          multiple substitutions, '/'-separated
    "IFRS", "Com1"     the reference (parent) sequence, i.e. no substitution

Unlike the original implementation, a substitution whose stated wild-type
residue disagrees with the parent sequence raises instead of printing a warning
and silently continuing -- a silently-ignored mutation produces a variant whose
label and sequence disagree, which poisons the training set.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass

AA20 = "ACDEFGHIKLMNPQRSTVWY"

_SUB_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


class VariantError(ValueError):
    """Raised when a variant string is malformed or inconsistent with its parent."""


@dataclass(frozen=True)
class Substitution:
    """A single substitution: `wt` at 1-based `pos` replaced by `mut`."""

    wt: str
    pos: int
    mut: str

    def __str__(self) -> str:
        return f"{self.wt}{self.pos}{self.mut}"

    @property
    def index(self) -> int:
        """0-based sequence index."""
        return self.pos - 1


def parse_variant(name: str) -> tuple[Substitution, ...]:
    """Parse a variant label into substitutions; the parent label yields ``()``.

    A label is treated as the parent whenever none of its '/'-separated tokens
    looks like a substitution (e.g. "IFRS", "Com1-IFRS", "WT").
    """
    tokens = [t.strip() for t in name.strip().split("/") if t.strip()]
    if not tokens:
        raise VariantError(f"empty variant label: {name!r}")

    matches = [_SUB_RE.match(t) for t in tokens]
    if all(m is None for m in matches):
        return ()
    if any(m is None for m in matches):
        bad = [t for t, m in zip(tokens, matches) if m is None]
        raise VariantError(f"variant {name!r} mixes substitutions with unparsable tokens {bad}")

    subs = tuple(Substitution(m.group(1), int(m.group(2)), m.group(3)) for m in matches)
    positions = [s.pos for s in subs]
    if len(set(positions)) != len(positions):
        raise VariantError(f"variant {name!r} substitutes the same position twice")
    return subs


def apply_variant(parent: str, name: str) -> str:
    """Build the variant sequence from `parent`, validating each stated wild-type residue."""
    seq = list(parent)
    for sub in parse_variant(name):
        if not 1 <= sub.pos <= len(parent):
            raise VariantError(
                f"{sub} in variant {name!r} is outside the parent sequence (length {len(parent)})"
            )
        if seq[sub.index] != sub.wt:
            raise VariantError(
                f"{sub} in variant {name!r} states wild-type {sub.wt} but the parent has "
                f"{seq[sub.index]} at position {sub.pos}"
            )
        seq[sub.index] = sub.mut
    return "".join(seq)


def build_sequences(parent: str, names) -> dict[str, str]:
    """Map each variant label to its full sequence."""
    return {name: apply_variant(parent, name) for name in names}


def enumerate_combinations(
    parent: str,
    singles,
    min_order: int = 2,
    max_order: int | None = None,
) -> list[str]:
    """Enumerate every combination of `singles` of order in ``[min_order, max_order]``.

    Combinations touching the same position twice are skipped (two substitutions
    at one site cannot co-occur). Returns labels in the canonical position-sorted
    order used by `combination_label`.
    """
    subs = []
    for name in singles:
        parsed = parse_variant(name)
        if len(parsed) != 1:
            raise VariantError(f"{name!r} is not a single substitution")
        apply_variant(parent, name)  # validate against the parent
        subs.append(parsed[0])

    hi = len(subs) if max_order is None else min(max_order, len(subs))
    out = []
    for order in range(max(min_order, 1), hi + 1):
        for combo in itertools.combinations(subs, order):
            if len({s.pos for s in combo}) != order:
                continue
            out.append(combination_label(combo))
    return out


def combination_label(subs) -> str:
    """Canonical label for a set of substitutions: position-sorted, '/'-joined."""
    return "/".join(str(s) for s in sorted(subs, key=lambda s: s.pos))


def saturation_scan(parent: str, positions=None, length: int | None = None) -> list[str]:
    """All 19 single substitutions at each requested position.

    Give either explicit 1-based `positions`, or `length` to scan the first
    `length` residues (the paper scanned the 240-residue tRNA-binding domain).
    The parent residue itself is excluded, so a scan of *n* positions yields
    19*n labels rather than the original implementation's 20*n (which emitted
    self-substitutions such as "A86A").
    """
    if positions is None:
        if length is None:
            raise ValueError("pass either positions or length")
        positions = range(1, length + 1)
    out = []
    for pos in positions:
        wt = parent[pos - 1]
        for aa in AA20:
            if aa != wt:
                out.append(f"{wt}{pos}{aa}")
    return out


def mutated_positions(names) -> list[int]:
    """Sorted union of the 1-based positions touched by any of `names`."""
    positions: set[int] = set()
    for name in names:
        positions.update(s.pos for s in parse_variant(name))
    return sorted(positions)
