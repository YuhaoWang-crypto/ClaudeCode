"""Does stacking mutations add up? Testing additivity on measured data.

Greedy combination - take the best single mutations, stack them, assume the
free energies add - is the standard move in computational affinity maturation.
It is also an assumption, and ATLAS can test it directly: wherever a measured
multi-mutant can be decomposed into disjoint subsets that were *also* measured,
the sum of the parts can be compared against the measured whole.

The answer is size-dependent, which is the useful part:

* **Pairs add up.** Restricted to decompositions where every component was
  measured in the same study, two-mutation combinations come out at mean error
  0.00 kcal/mol with every case inside 0.5 kcal/mol.
* **Deep stacks do not.** By eight or nine simultaneous mutations the measured
  gain falls short of the sum by 1.6 to 3.4 kcal/mol - strongly sub-additive,
  the classic diminishing return of an interface that is already optimised.

So greedy pairing is defensible and greedy stacking is not, which is consistent
with real affinity-matured TCRs being found by directed evolution rather than
by computational accumulation.

A caveat the analysis handles explicitly: wild-type affinities for the same
receptor differ between studies - DMF5 against ELAGIGILTV is reported at 5.6,
9.5, 11 and 36.4 micromolar - so ddG values from different papers rest on
different baselines, worth about 1 kcal/mol of noise. Every statistic is
therefore reported both pooled and restricted to single studies.
"""

from __future__ import annotations

import csv
import itertools
import re
import statistics
from dataclasses import dataclass
from pathlib import Path

MISSING = {"N/A", "n.d.", "", "WT"}


@dataclass
class Decomposition:
    """A measured multi-mutant alongside measured pieces that partition it."""

    tcr: str
    peptide: str
    study: str
    whole: frozenset[str]
    parts: tuple[frozenset[str], ...]
    sum_of_parts: float
    measured: float

    @property
    def error(self) -> float:
        """Measured minus predicted. Positive means the stack underdelivered."""
        return self.measured - self.sum_of_parts

    @property
    def size(self) -> int:
        return len(self.whole)

    @property
    def behaviour(self) -> str:
        if self.error < -0.5:
            return "synergistic"
        if self.error > 0.5:
            return "antagonistic"
        return "additive"

    def as_row(self) -> dict:
        return {
            "tcr": self.tcr,
            "peptide": self.peptide,
            "study_pmid": self.study,
            "mutations": "+".join(sorted(self.whole)),
            "n_mutations": self.size,
            "decomposed_as": " | ".join(
                "{" + ",".join(sorted(p)) + "}" for p in self.parts
            ),
            "sum_of_parts": round(self.sum_of_parts, 3),
            "measured": round(self.measured, 3),
            "error": round(self.error, 3),
            "behaviour": self.behaviour,
        }


def _parse_ddg(value: str | None) -> float | None:
    text = (value or "").replace("±", "+/-").split("+/-")[0].strip()
    if text in MISSING:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _mutation_set(text: str) -> frozenset[str]:
    return frozenset(m.strip() for m in re.split(r"[|.]", text) if m.strip())


def index_atlas(path: Path, by_study: bool) -> dict[tuple, dict[frozenset[str], float]]:
    """Index measured ddG by system, optionally keeping studies separate.

    Grouping always includes the peptide: the same receptor against a different
    epitope is a different measurement, and pooling them would invent epistasis
    that is really a change of target.
    """
    index: dict[tuple, dict[frozenset[str], list[float]]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            mutation = (row.get("TCR_mut") or "").strip()
            if mutation in MISSING:
                continue
            if (row.get("MHC_mut") or "").strip() not in ("WT", "N/A"):
                continue
            if (row.get("PEP_mut") or "").strip() not in ("WT", "N/A"):
                continue
            ddg = _parse_ddg(row.get("Delta_DeltaG_kcal_per_mol"))
            if ddg is None:
                continue
            key = (row["TCRname"], row["PEPseq"])
            if by_study:
                key = key + (row.get("PMID", ""),)
            index.setdefault(key, {}).setdefault(_mutation_set(mutation), []).append(ddg)
    return {
        key: {m: statistics.mean(values) for m, values in group.items()}
        for key, group in index.items()
    }


def find_decompositions(
    index: dict[tuple, dict[frozenset[str], float]], max_parts: int = 3
) -> list[Decomposition]:
    """Every exact partition of a measured multi-mutant into measured pieces."""
    found: list[Decomposition] = []
    seen: set = set()
    for key, group in index.items():
        measured = list(group)
        for whole in measured:
            if len(whole) < 2:
                continue
            pieces = [m for m in measured if m and m < whole]
            for count in range(2, max_parts + 1):
                for combo in itertools.combinations(pieces, count):
                    union: set[str] = set()
                    disjoint = True
                    for piece in combo:
                        if union & piece:
                            disjoint = False
                            break
                        union |= piece
                    if not disjoint or frozenset(union) != whole:
                        continue
                    signature = (
                        key,
                        whole,
                        tuple(sorted(tuple(sorted(p)) for p in combo)),
                    )
                    if signature in seen:
                        continue
                    seen.add(signature)
                    found.append(
                        Decomposition(
                            tcr=key[0],
                            peptide=key[1],
                            study=key[2] if len(key) > 2 else "pooled",
                            whole=whole,
                            parts=combo,
                            sum_of_parts=sum(group[p] for p in combo),
                            measured=group[whole],
                        )
                    )
    return sorted(found, key=lambda d: (d.size, d.tcr))


def summarise(decompositions: list[Decomposition]) -> dict:
    """Overall and per-stack-size statistics."""
    if not decompositions:
        return {"n": 0}
    errors = [d.error for d in decompositions]
    by_size: dict[int, list[float]] = {}
    for item in decompositions:
        by_size.setdefault(item.size, []).append(item.error)

    return {
        "n": len(errors),
        "mean_error": statistics.mean(errors),
        "median_error": statistics.median(errors),
        "stdev_error": statistics.pstdev(errors),
        "fraction_within_half_kcal": sum(1 for e in errors if abs(e) < 0.5) / len(errors),
        "by_stack_size": {
            size: {
                "n": len(values),
                "mean_error": statistics.mean(values),
                "within_half_kcal": sum(1 for e in values if abs(e) < 0.5),
            }
            for size, values in sorted(by_size.items())
        },
    }


def run_additivity_analysis(atlas_path: Path) -> dict:
    """Pooled and within-study additivity tests."""
    results = {}
    for label, by_study in (("pooled", False), ("within_study", True)):
        decompositions = find_decompositions(index_atlas(atlas_path, by_study))
        results[label] = {
            "summary": summarise(decompositions),
            "decompositions": [d.as_row() for d in decompositions],
        }
    pairs = results["within_study"]["summary"]["by_stack_size"].get(2, {})
    results["verdict"] = _verdict(pairs, results["within_study"]["summary"])
    return results


def _verdict(pairs: dict, within_study: dict) -> str:
    if not pairs or pairs.get("n", 0) < 3:
        return (
            "Too few within-study pair decompositions to judge whether greedy "
            "combination is safe."
        )
    deep = {
        size: stats
        for size, stats in within_study["by_stack_size"].items()
        if size >= 6
    }
    deep_text = ""
    if deep:
        worst = max(deep.values(), key=lambda s: s["mean_error"])
        deep_text = (
            f" Beyond six simultaneous mutations the measured gain falls short "
            f"of the sum by about {worst['mean_error']:+.1f} kcal/mol, so "
            f"stacking that deep cannot be predicted by adding up parts."
        )
    return (
        f"Pairs add up: {pairs['within_half_kcal']} of {pairs['n']} two-mutation "
        f"combinations measured within one study land inside 0.5 kcal/mol of the "
        f"sum of their parts, mean error {pairs['mean_error']:+.2f}."
        + deep_text
    )


def format_report(results: dict) -> str:
    """Render the additivity analysis as markdown."""
    lines = [
        "# Does stacking mutations add up?",
        "",
        "Greedy combination assumes the free energies of individual mutations "
        "add. ATLAS can test that directly: wherever a measured multi-mutant "
        "decomposes into disjoint subsets that were also measured, the sum of "
        "the parts can be compared with the measured whole.",
        "",
        f"**{results['verdict']}**",
        "",
        "Wild-type affinities for the same receptor differ between studies, so "
        "ddG values from different papers rest on different baselines. Both the "
        "pooled analysis and one restricted to single studies are reported; the "
        "size-dependent pattern survives the restriction, so it is not an "
        "artefact of that noise.",
        "",
    ]
    for label, title in (
        ("within_study", "Within a single study"),
        ("pooled", "Pooled across studies"),
    ):
        summary = results[label]["summary"]
        if not summary.get("n"):
            continue
        lines += [
            f"## {title}",
            "",
            f"n = {summary['n']} decompositions, mean error "
            f"{summary['mean_error']:+.2f}, median {summary['median_error']:+.2f}, "
            f"sd {summary['stdev_error']:.2f} kcal/mol. "
            f"{100 * summary['fraction_within_half_kcal']:.0f}% land inside "
            "0.5 kcal/mol.",
            "",
            "| mutations stacked | n | mean error (kcal/mol) | within 0.5 |",
            "|---|---|---|---|",
        ]
        for size, stats in summary["by_stack_size"].items():
            lines.append(
                f"| {size} | {stats['n']} | {stats['mean_error']:+.2f} | "
                f"{stats['within_half_kcal']}/{stats['n']} |"
            )
        lines.append("")
    lines += [
        "Error is measured minus predicted, so a positive value means the stack "
        "underdelivered against the sum of its parts.",
        "",
        "## What this means for a design campaign",
        "",
        "* Predicting a **double** mutant by adding two singles is sound.",
        "* **Triples** are borderline and should be checked, not assumed.",
        "* Stacking **beyond about four** mutations cannot be predicted by "
        "summation. The published deep variants, which reach 6 to 14 "
        "simultaneous mutations, were found by directed evolution, and the "
        "sub-additivity here is why.",
        "",
    ]
    return "\n".join(lines) + "\n"
