"""Developability and specificity gates.

Two things sink affinity-matured TCRs that a binding score cannot see:

1. **Sequence liabilities** introduced into a CDR loop. A new N-glycosylation
   sequon or a free cysteine in a paratope loop can cost you the molecule
   regardless of its affinity.
2. **Where the new affinity came from.** Affinity gained at MHC-contacting
   positions raises recognition of the MHC surface itself, which is
   peptide-independent. That is the mechanism behind the best-known clinical
   failures of affinity-enhanced TCRs, so it is gated here rather than scored.

These gates are rule-based and conservative. They flag risk; they do not
measure it.
"""

from __future__ import annotations

import re

from .constants import HYDROPATHY
from .scan import Mutation
from .structure import Complex

# Liabilities that disqualify a candidate outright.
HARD_LIABILITIES = ("new N-glycosylation sequon", "new free cysteine")

# Mean-hydropathy ceiling for a CDR loop, above which aggregation and
# non-specific binding risk rises.
LOOP_HYDROPATHY_CEILING = 0.8


def mutate_sequence(sequence: str, index: int, mutant: str) -> str:
    return sequence[:index] + mutant + sequence[index + 1 :]


def _loop_span(complex_obj: Complex, mutation: Mutation) -> tuple[int, int] | None:
    if not mutation.loop:
        return None
    span = complex_obj.loops.get(mutation.loop)
    if span is None:
        return None
    _, start, end = span
    return start, end


def _residue_index(complex_obj: Complex, mutation: Mutation) -> int | None:
    residue = complex_obj.residue_by_label(mutation.chain, mutation.number)
    return residue.index if residue else None


def assess_liabilities(complex_obj: Complex, mutation: Mutation) -> tuple[str, ...]:
    """Find sequence liabilities that the mutation creates but the wild type lacks."""
    index = _residue_index(complex_obj, mutation)
    if index is None:
        return ()
    wild = complex_obj.sequences[mutation.chain]
    mutant = mutate_sequence(wild, index, mutation.mutant)

    # Compare only a local window, so pre-existing motifs elsewhere in the
    # chain are not blamed on this mutation.
    lo, hi = max(0, index - 3), min(len(wild), index + 4)
    wild_window, mutant_window = wild[lo:hi], mutant[lo:hi]

    found: list[str] = []

    def created(pattern: str) -> bool:
        return len(re.findall(pattern, mutant_window)) > len(re.findall(pattern, wild_window))

    if created(r"N[^P][ST]"):
        found.append("new N-glycosylation sequon")
    if mutation.mutant == "C":
        found.append("new free cysteine")
    if created(r"NG"):
        found.append("new Asn deamidation motif (NG)")
    if created(r"D[GP]"):
        found.append("new Asp isomerisation motif (DG/DP)")
    if mutation.mutant == "M":
        found.append("new methionine (oxidation-prone)")

    span = _loop_span(complex_obj, mutation)
    if span:
        start, end = span
        loop_mutant = mutant[start:end]
        mean_hydropathy = sum(HYDROPATHY[a] for a in loop_mutant if a in HYDROPATHY) / max(
            1, len(loop_mutant)
        )
        loop_wild = wild[start:end]
        mean_wild = sum(HYDROPATHY[a] for a in loop_wild if a in HYDROPATHY) / max(
            1, len(loop_wild)
        )
        if mean_hydropathy > LOOP_HYDROPATHY_CEILING and mean_hydropathy > mean_wild:
            found.append(
                f"{mutation.loop} mean hydropathy {mean_hydropathy:+.2f} "
                f"exceeds the aggregation-risk ceiling"
            )
    return tuple(found)


def apply_gates(
    complex_obj: Complex,
    mutations: list[Mutation],
    z_threshold: float = 1.0,
) -> list[Mutation]:
    """Attach liabilities and a final verdict to every scored mutation.

    Verdicts
    --------
    ``advance``            clean, peptide-focused and above the score threshold
    ``watch-specificity``  worth testing, but gains lean on MHC contacts
    ``reject-specificity`` affinity would come almost entirely from MHC contact
    ``reject-liability``   introduces a disqualifying sequence liability
    ``below-threshold``    no predicted gain
    """
    for mutation in mutations:
        mutation.liabilities = assess_liabilities(complex_obj, mutation)
        hard = [item for item in mutation.liabilities if item in HARD_LIABILITIES]

        # Peptide-focused candidates are judged against the other
        # peptide-focused candidates, not against the MHC-contacting
        # positions that dominate the interface.
        effective_z = (
            mutation.peptide_z if mutation.peptide_z is not None else mutation.z_score
        )

        if hard:
            mutation.verdict = "reject-liability"
        elif mutation.specificity_flag == "mhc-dominated":
            mutation.verdict = "reject-specificity"
        elif effective_z < z_threshold:
            mutation.verdict = "below-threshold"
        elif mutation.specificity_flag == "mhc-leaning":
            mutation.verdict = "watch-specificity"
        else:
            mutation.verdict = "advance"
    return mutations


def summarise_gates(mutations: list[Mutation]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for mutation in mutations:
        counts[mutation.verdict] = counts.get(mutation.verdict, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))
