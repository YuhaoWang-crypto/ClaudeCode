"""Retrospective benchmarking of the triage score.

A scoring function with no benchmark is a rumour. This module runs the only
retrospective check that can be built from public structures alone: take a
wild-type TCR-pMHC complex and a deposited structure of an engineered variant
of the same receptor, recover the substitutions that distinguish them, and ask
where the triage score ranked them in a blind scan of the wild type.

The reference case is the DMF5 TCR against the MART-1 decamer ELAGIGILTV on
HLA-A2:

* ``3QDG`` - wild-type DMF5 bound to the decamer
* ``4L3E`` - the computationally designed high-affinity variant

The designed variant carries two substitutions and is reported to bind roughly
400-fold more tightly than wild type while retaining specificity for the
cognate ligand (Pierce et al., PLoS Comput Biol 2014,
doi:10.1371/journal.pcbi.1003478). That design came from Rosetta-style physics
plus experimental iteration, so it is a fair target for any cheaper method to
be measured against.

Two mutations is not a benchmark. It is a smoke test, and the honest result is
reported whether or not it flatters the score. For a real benchmark, point
``benchmark_from_csv`` at the ATLAS database, which curates wild-type and
mutant TCR-pMHC affinities alongside structures.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .liability import apply_gates
from .scan import Mutation, ranking_key as _ranking_key, scan_complex
from .structure import Complex, load_complex

# Known wild-type / engineered-variant structure pairs.
REFERENCE_PAIRS = {
    "DMF5": {
        "wild_type": "3QDG",
        "variant": "4L3E",
        "description": (
            "DMF5 TCR vs MART-1 decamer ELAGIGILTV / HLA-A2; designed variant "
            "reported at ~400-fold improved affinity with retained specificity"
        ),
        "reference": "Pierce et al., PLoS Comput Biol 2014, 10.1371/journal.pcbi.1003478",
    }
}


@dataclass
class RankedKnown:
    """Where a known mutation landed in a blind scan of the wild type."""

    name: str
    loop: str | None
    rank: int
    total: int
    percentile: float
    position_rank: int
    position_total: int
    triage_score: float
    z_score: float
    predicted_ddg: float | None
    verdict: str
    contacts_peptide: int
    contacts_mhc: int
    mhc_contact_fraction: float
    notes: tuple[str, ...]
    scanned: bool = True

    def as_row(self) -> dict:
        row = self.__dict__.copy()
        row["notes"] = "; ".join(self.notes)
        return row


def diff_structures(
    wild_type: Complex, variant: Complex
) -> list[tuple[str, str, int, str, str, str | None]]:
    """Recover substitutions between two structures of the same receptor.

    Residues are matched by author numbering within the corresponding chain
    role, which is what makes the comparison meaningful across two separate
    depositions whose chain letters may differ.

    Returns tuples of ``(role, chain_id_in_wild_type, number, wt_aa, var_aa, loop)``.
    """
    differences = []
    for role in ("tcr_alpha", "tcr_beta", "peptide", "mhc"):
        wt_chain = wild_type.chains.get(role)
        var_chain = variant.chains.get(role)
        if not wt_chain or not var_chain:
            continue
        wt_map = {
            r.number: r for r in wild_type.residues if r.chain == wt_chain
        }
        var_map = {r.number: r for r in variant.residues if r.chain == var_chain}
        for number in sorted(set(wt_map) & set(var_map)):
            if wt_map[number].aa != var_map[number].aa:
                differences.append(
                    (
                        role,
                        wt_chain,
                        number,
                        wt_map[number].aa,
                        var_map[number].aa,
                        wt_map[number].loop or var_map[number].loop,
                    )
                )
    return differences


def locate_in_scan(
    mutations: list[Mutation], chain: str, number: int, mutant_aa: str
) -> RankedKnown | None:
    """Find a specific substitution in a completed scan and report its rank."""
    ranked = sorted(mutations, key=_ranking_key)
    total = len(ranked)
    match = None
    rank = -1
    for index, mutation in enumerate(ranked, start=1):
        if (
            mutation.chain == chain
            and mutation.number == number
            and mutation.mutant == mutant_aa
        ):
            match, rank = mutation, index
            break
    if match is None:
        return None

    same_position = sorted(
        [m for m in mutations if m.chain == chain and m.number == number],
        key=_ranking_key,
    )
    position_rank = next(
        i for i, m in enumerate(same_position, start=1) if m.mutant == mutant_aa
    )
    return RankedKnown(
        name=match.name,
        loop=match.loop,
        rank=rank,
        total=total,
        percentile=100.0 * rank / total,
        position_rank=position_rank,
        position_total=len(same_position),
        triage_score=match.triage_score,
        z_score=match.z_score,
        predicted_ddg=match.predicted_ddg,
        verdict=match.verdict,
        contacts_peptide=match.contacts_peptide,
        contacts_mhc=match.contacts_mhc,
        mhc_contact_fraction=match.mhc_contact_fraction,
        notes=match.notes,
    )


def retrospective_check(
    wild_type_id: str,
    variant_id: str,
    cache_dir: Path,
    min_contacts: int = 3,
    z_threshold: float = 1.0,
    ddg_model=None,
) -> dict:
    """Score a wild-type complex blind, then look up the known mutations.

    When ``ddg_model`` is supplied the ranking uses its predicted ddG instead
    of the hand-weighted score. For the result to mean anything, that model
    must have been fitted without the structures being tested here; see
    ``benchmark.fit_excluding``.
    """
    wild_type = load_complex(wild_type_id, cache_dir)
    variant = load_complex(variant_id, cache_dir)

    differences = diff_structures(wild_type, variant)
    tcr_differences = [d for d in differences if d[0] in ("tcr_alpha", "tcr_beta")]
    other_differences = [d for d in differences if d[0] not in ("tcr_alpha", "tcr_beta")]

    mutations = apply_gates(
        wild_type,
        scan_complex(wild_type, min_contacts=min_contacts, ddg_model=ddg_model),
        z_threshold=z_threshold,
    )

    results: list[RankedKnown] = []
    missing: list[str] = []
    for _role, chain, number, wt_aa, var_aa, loop in tcr_differences:
        found = locate_in_scan(mutations, chain, number, var_aa)
        if found is None:
            missing.append(f"{chain}:{wt_aa}{number}{var_aa}")
            continue
        results.append(found)

    percentiles = [r.percentile for r in results]
    median_percentile = (
        sorted(percentiles)[len(percentiles) // 2] if percentiles else None
    )
    advanced = [r for r in results if r.verdict in ("advance", "watch-specificity")]

    return {
        "wild_type": wild_type_id.upper(),
        "variant": variant_id.upper(),
        "peptide": wild_type.peptide,
        "mutations_scanned": len(mutations),
        "known_mutations": [r.as_row() for r in results],
        "known_mutations_not_scanned": missing,
        "non_tcr_differences": [
            f"{role}:{wt}{number}{var}" for role, _c, number, wt, var, _l in other_differences
        ],
        "median_percentile": median_percentile,
        "would_have_advanced": [r.name for r in advanced],
        "verdict": _summarise(results, median_percentile),
    }


def _summarise(results: list[RankedKnown], median_percentile: float | None) -> str:
    if not results:
        return (
            "No known mutations were recoverable from this structure pair, so "
            "nothing was tested."
        )
    advanced = [r for r in results if r.verdict in ("advance", "watch-specificity")]
    n = len(results)
    if median_percentile is None:
        return "inconclusive"

    # Being ranked first among the 19 substitutions available at a position is
    # informative even when the absolute threshold was not met, because that is
    # the choice a designer actually faces once a position is selected.
    top_at_position = [r for r in results if r.position_rank <= 3]
    position_note = (
        f" {len(top_at_position)} of {n} ranked in the top 3 of the 19 "
        f"substitutions available at their own position"
        f" ({', '.join(f'{r.name} {r.position_rank}/{r.position_total}' for r in top_at_position)})."
        if top_at_position
        else ""
    )

    if median_percentile <= 10 and len(advanced) == n:
        return (
            f"Enriched: all {n} experimentally validated mutations landed in the "
            f"top {median_percentile:.0f}% of the blind scan and all would have "
            f"been shortlisted.{position_note} With n this small that is "
            "encouraging, not proof."
        )
    if len(advanced) == 0:
        verdict = (
            f"FAILED: none of the {n} experimentally validated mutations would "
            f"have been shortlisted; median rank was the "
            f"{median_percentile:.0f}th percentile.{position_note}"
        )
        if top_at_position:
            verdict += (
                " So there is some signal in which substitution to choose at a "
                "given position, but not enough to clear the threshold for "
                "calling a mutation an improvement at all."
            )
        else:
            verdict += (
                " That is indistinguishable from picking at random, and the "
                "score should not be trusted to enrich for real affinity gains."
            )
        return verdict
    return (
        f"Mixed: {len(advanced)} of {n} validated mutations would have been "
        f"shortlisted, median percentile {median_percentile:.0f}.{position_note} "
        "Too weak to rely on and too small a sample to characterise."
    )


def benchmark_from_csv(
    pdb_id: str,
    csv_path: Path,
    cache_dir: Path,
    min_contacts: int = 3,
    z_threshold: float = 1.0,
) -> dict:
    """Benchmark the score against a user-supplied set of measured mutations.

    The CSV needs columns ``chain``, ``number``, ``mutant`` and ``ddg``, where
    ``ddg`` is the measured change in binding free energy in kcal/mol using the
    convention that a negative value means tighter binding. ATLAS is the obvious
    source for TCR-pMHC mutant affinities.

    Reports Spearman correlation between the triage score and measured ddG. Any
    positive correlation here would be new information, since no TCR-specific
    ddG benchmark of this score exists.
    """
    complex_obj = load_complex(pdb_id, cache_dir)
    mutations = apply_gates(
        complex_obj,
        scan_complex(complex_obj, min_contacts=min_contacts, include_wild_type=False),
        z_threshold=z_threshold,
    )
    lookup = {(m.chain, m.number, m.mutant): m for m in mutations}

    paired: list[tuple[float, float, str]] = []
    unmatched: list[str] = []
    with csv_path.open() as handle:
        for row in csv.DictReader(handle):
            key = (row["chain"].strip(), int(row["number"]), row["mutant"].strip())
            measured = float(row["ddg"])
            mutation = lookup.get(key)
            if mutation is None:
                unmatched.append(f"{key[0]}:{key[1]}{key[2]}")
                continue
            paired.append((mutation.triage_score, measured, mutation.name))

    spearman = _spearman([p[0] for p in paired], [p[1] for p in paired])
    return {
        "pdb_id": complex_obj.pdb_id,
        "measured_mutations_matched": len(paired),
        "measured_mutations_unmatched": unmatched,
        # A higher triage score is meant to mean tighter binding, and tighter
        # binding is a more negative ddG, so a working score gives a negative
        # correlation here.
        "spearman_triage_vs_ddg": spearman,
        "interpretation": (
            "The score is meant to increase as binding tightens, and ddG becomes "
            "more negative as binding tightens, so a useful score gives a "
            "negative Spearman value. Values near zero mean no signal."
        ),
        "pairs": [{"mutation": name, "triage": t, "ddg": d} for t, d, name in paired],
    }


def _spearman(a: list[float], b: list[float]) -> float | None:
    if len(a) < 3:
        return None
    ranks_a, ranks_b = _rank(a), _rank(b)
    n = len(a)
    mean_a = sum(ranks_a) / n
    mean_b = sum(ranks_b) / n
    cov = sum((x - mean_a) * (y - mean_b) for x, y in zip(ranks_a, ranks_b))
    var_a = sum((x - mean_a) ** 2 for x in ranks_a) ** 0.5
    var_b = sum((y - mean_b) ** 2 for y in ranks_b) ** 0.5
    if var_a == 0 or var_b == 0:
        return None
    return cov / (var_a * var_b)


def _rank(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        stop = index
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[index]]:
            stop += 1
        average = (index + stop) / 2 + 1
        for position in range(index, stop + 1):
            ranks[order[position]] = average
        index = stop + 1
    return ranks


def format_report(result: dict) -> str:
    """Render a retrospective check as markdown."""
    lines = [
        f"# Retrospective check: {result['wild_type']} -> {result['variant']}",
        "",
        f"Peptide: `{result['peptide']}`",
        "",
        f"Scanned **{result['mutations_scanned']} single mutants** of "
        f"{result['wild_type']} with no knowledge of the engineered variant, then "
        "looked up where the known substitutions ranked.",
        "",
        "| mutation | loop | rank | percentile | rank at position | predicted ddG | "
        "verdict | pep contacts | MHC contacts |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in result["known_mutations"]:
        predicted = row.get("predicted_ddg")
        predicted_text = "-" if predicted is None else f"{predicted:+.2f}"
        lines.append(
            f"| `{row['name']}` | {row['loop'] or '-'} | "
            f"{row['rank']}/{row['total']} | top {row['percentile']:.0f}% | "
            f"{row['position_rank']}/{row['position_total']} | {predicted_text} | "
            f"{row['verdict']} | "
            f"{row['contacts_peptide']} | {row['contacts_mhc']} |"
        )
    lines += [
        "",
        f"**Outcome:** {result['verdict']}",
        "",
    ]
    if result["known_mutations_not_scanned"]:
        lines += [
            "Not reachable by the scan (the position fell below the contact "
            "threshold or sits outside the paratope): "
            + ", ".join(f"`{m}`" for m in result["known_mutations_not_scanned"]),
            "",
        ]
    if result["non_tcr_differences"]:
        lines += [
            "Differences outside the TCR chains, ignored here: "
            + ", ".join(f"`{m}`" for m in result["non_tcr_differences"]),
            "",
        ]
    return "\n".join(lines) + "\n"
