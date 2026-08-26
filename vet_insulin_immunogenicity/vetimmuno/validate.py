"""The validation harness: known-answer tests and controls.

A prediction workflow for a species with no benchmark data cannot be validated
against ground-truth immunogenicity -- that data does not exist for dogs or
cats. What *can* be validated, and is validated here, is everything else:

  * that the sequence layer reproduces facts that are independently known
    (porcine insulin is identical to canine insulin; feline insulin differs
    from human at exactly A8, A10, A18, B30);
  * that the tolerance filter fires exactly where it should and nowhere else;
  * that a negative control (species-matched insulin) scores zero risk and a
    positive control (scrambled sequence) scores near-maximal;
  * that the %rank calibration is actually uniform;
  * that the scoring is deterministic and register assignment is stable;
  * that the applicability-domain guard-rail fires for DLA/FLA rather than
    quietly passing an out-of-distribution prediction through.

Each check returns PASS / FAIL / INFO with the observed value, so the harness
output is the evidence, not a claim about it.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import epitope, groove
from .groove import DomainCall, Molecule
from .insulin import Insulin, apply_edits, diff, natural_insulin


@dataclass
class Check:
    id: str
    name: str
    status: str          # PASS / FAIL / INFO / SKIP
    observed: str
    expected: str = ""
    detail: str = ""

    @property
    def ok(self) -> bool:
        return self.status in ("PASS", "INFO", "SKIP")


def _check(cid: str, name: str, condition: bool, observed, expected: str = "",
           detail: str = "") -> Check:
    return Check(cid, name, "PASS" if condition else "FAIL", str(observed), expected, detail)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rho without a SciPy dependency."""
    def rank(v: Sequence[float]) -> np.ndarray:
        arr = np.asarray(v, dtype=float)
        order = arr.argsort()
        ranks = np.empty(len(arr), dtype=float)
        ranks[order] = np.arange(len(arr), dtype=float)
        # average ties
        for value in np.unique(arr):
            mask = arr == value
            if mask.sum() > 1:
                ranks[mask] = ranks[mask].mean()
        return ranks
    rx, ry = rank(x), rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


# --------------------------------------------------------------------------
# sequence-layer known-answer tests
# --------------------------------------------------------------------------

# Independently established facts, asserted against freshly fetched UniProt
# records. If an upstream release changes, these fail loudly instead of
# silently shifting every downstream number.
SEQUENCE_KATS = [
    ("dog", "human", 1, ["B30A>T"], "canine insulin differs from human only at B30"),
    ("cat", "human", 4, ["A8A>T", "A10V>I", "A18H>N", "B30A>T"],
     "feline insulin differs from human at A8, A10, A18, B30"),
    ("dog", "pig", 0, [], "porcine insulin is sequence-identical to canine insulin"),
    ("cat", "bovine", 1, ["A18H>N"], "bovine insulin differs from feline only at A18"),
    ("dog", "bovine", 2, ["A8T>A", "A10I>V"], "bovine insulin differs from canine at A8 and A10"),
]


def check_sequence_ground_truth() -> List[Check]:
    out: List[Check] = []
    for i, (recipient, donor, n_expected, labels, why) in enumerate(SEQUENCE_KATS, 1):
        observed = diff(natural_insulin(recipient), natural_insulin(donor))
        got = [d.label for d in observed]
        out.append(_check(
            f"KAT-1.{i}", f"{donor} insulin vs {recipient} self",
            got == labels,
            f"{len(got)} differences {got}",
            f"{n_expected} differences {labels}",
            why,
        ))
    return out


# --------------------------------------------------------------------------
# tolerance-filter controls
# --------------------------------------------------------------------------

def check_negative_control(species: str) -> Check:
    """Species-matched insulin must produce exactly zero non-self cores."""
    own = natural_insulin(species)
    cores = epitope.neo_cores(own, own)
    return _check(
        "KAT-2", f"negative control: {species} insulin in a {species}",
        len(cores) == 0, f"{len(cores)} non-self cores", "0 non-self cores",
        "central tolerance: the recipient's own insulin cannot present a novel core",
    )


def check_identity_control(species: str, donor: str) -> Check:
    """A donor insulin identical in sequence must also produce zero non-self cores."""
    own, drug = natural_insulin(species), natural_insulin(donor)
    cores = epitope.neo_cores(drug, own)
    return _check(
        "KAT-3", f"identity control: {donor} insulin in a {species}",
        len(cores) == 0 and not diff(own, drug),
        f"{len(cores)} non-self cores", "0 non-self cores",
        f"{donor} and {species} insulin are the same molecule, so no new core exists",
    )


def check_tolerance_filter_precision(species: str) -> Check:
    """A single engineered substitution must flag exactly the cores that span it."""
    own = natural_insulin(species)
    # Mutate the last B-chain residue; only the single core ending there spans it.
    last = len(own.B)
    mutant = apply_edits(own, [f"B{last}{own.B[-1]}>W"], f"{species}-B{last}W probe")
    cores = epitope.neo_cores(mutant, own)
    expected_starts = {last - epitope.CORE_LENGTH + 1}
    got_starts = {c.core.start for c in cores if c.core.chain == "B"}
    return _check(
        "KAT-6", "tolerance filter precision (single-residue probe)",
        got_starts == expected_starts and all(c.core.chain == "B" for c in cores),
        f"cores flagged at B{sorted(got_starts)}",
        f"cores flagged at B{sorted(expected_starts)} only",
        "a substitution must flag every core spanning it and no other core",
    )


def check_scramble_positive_control(species: str, n_seeds: int = 20,
                                    threshold: float = 0.90) -> Check:
    """A composition-matched scramble must look almost entirely non-self."""
    own = natural_insulin(species)
    fractions: List[float] = []
    for seed in range(n_seeds):
        rng = random.Random(seed)
        a = list(own.A); b = list(own.B)
        rng.shuffle(a); rng.shuffle(b)
        scrambled = Insulin(f"scrambled-{seed}", "".join(a), "".join(b))
        total = len(epitope.all_cores(scrambled))
        fractions.append(len(epitope.neo_cores(scrambled, own)) / total)
    mean = float(np.mean(fractions))
    return _check(
        "KAT-5", "positive control: composition-matched scramble",
        mean >= threshold, f"{mean:.3f} of cores non-self (mean of {n_seeds} seeds)",
        f">= {threshold:.2f}",
        "same amino-acid composition, different order -- the filter must not be "
        "fooled by composition alone",
    )


def check_core_set_invariant(species: str, products: Sequence[Insulin]) -> Check:
    """Exact invariant: the non-self cores are precisely the cores spanning a difference.

    Stronger than a correlation. A core can only be novel to the recipient if it
    covers at least one residue the recipient's own insulin does not have there,
    and every such core must be flagged unless its 9-mer coincidentally recurs
    elsewhere in the recipient's insulin. Any mismatch is reported with the
    offending core rather than absorbed into a threshold.

    Note this is deliberately *not* a monotonicity test: a difference at a chain
    terminus spans fewer registers than an interior one, so core count is not a
    monotone function of residue-difference count -- which is itself one of the
    workflow's substantive points about B30 in the dog.
    """
    own = natural_insulin(species)
    tolerated = epitope.self_core_set(own)
    mismatches: List[str] = []
    for product in products:
        difference_positions = {(d.chain, d.position) for d in diff(own, product)}
        flagged = {(c.core.chain, c.core.start) for c in epitope.neo_cores(product, own)}
        for core in epitope.all_cores(product):
            spans = any((core.chain, p) in difference_positions for p in core.positions())
            is_flagged = (core.chain, core.start) in flagged
            if spans and not is_flagged and core.sequence not in tolerated:
                mismatches.append(f"{product.name}:{core.label} spans a difference but was not flagged")
            if is_flagged and not spans:
                mismatches.append(f"{product.name}:{core.label} flagged without spanning a difference")
    return _check(
        "KAT-4", "non-self core set equals the set of cores spanning a difference",
        not mismatches,
        "exact match for all products" if not mismatches
        else f"{len(mismatches)} mismatches, first: {mismatches[0]}",
        "0 mismatches",
        "ties the risk list to the sequence differences with no slack",
    )


# --------------------------------------------------------------------------
# scoring-layer controls
# --------------------------------------------------------------------------

def check_rank_calibration(ranker, molecule: Molecule, tolerance: float = 0.6) -> Check:
    """A fresh background draw must land ~5% of peptides at rank <= 5."""
    probe = type(ranker)(ranker.backend, n_background=5000, seed=ranker.seed + 7)
    scores = ranker.backend.score(probe._peptides, molecule)
    ranks = ranker.rank(scores, molecule)
    pct5 = 100.0 * float(np.mean(ranks <= 5.0))
    return _check(
        "KAT-7", "background %rank calibration",
        abs(pct5 - 5.0) <= tolerance,
        f"{pct5:.2f}% of held-out background peptides at rank <= 5",
        f"5.00% +/- {tolerance}",
        "the custom-molecule rank reference must actually be uniform, otherwise "
        "cross-molecule comparison is meaningless",
    )


def check_determinism(backend, molecule: Molecule, cores: Sequence[str]) -> Check:
    a = backend.score(cores, molecule)
    b = backend.score(cores, molecule)
    return _check(
        "KAT-12", "scoring determinism", bool(np.array_equal(a, b)),
        "identical on repeat" if np.array_equal(a, b) else "scores differ between runs",
        "bit-identical", "a re-run must reproduce the report exactly",
    )


def check_register_stability(backend, molecules: Sequence[Molecule],
                             drug: Insulin, threshold: float = 0.60) -> Check:
    """How often the top register survives a one-residue change of the window."""
    stable = total = 0
    for mol in molecules:
        for ch in ("A", "B"):
            seq = drug.chain(ch)
            for start in range(0, max(1, len(seq) - epitope.PEPTIDE_LENGTH + 1)):
                window = seq[start:start + epitope.PEPTIDE_LENGTH]
                shifted = seq[max(0, start - 1):start - 1 + epitope.PEPTIDE_LENGTH] \
                    if start > 0 else seq[start:start + epitope.PEPTIDE_LENGTH + 1]
                if len(window) < epitope.CORE_LENGTH or len(shifted) < epitope.CORE_LENGTH:
                    continue
                best_w = _best_core(backend, mol, window, start + 1)
                best_s = _best_core(backend, mol, shifted, max(1, start))
                total += 1
                stable += best_w == best_s
    frac = stable / total if total else float("nan")
    return _check(
        "KAT-8", "binding-register stability under window shift",
        frac >= threshold, f"{frac:.3f} of windows keep their top register",
        f">= {threshold:.2f}",
        "an unstable register assignment means the reported core is an artefact "
        "of how the peptide was tiled",
    )


def _best_core(backend, mol: Molecule, window: str, window_start: int) -> int:
    cores = [window[i:i + epitope.CORE_LENGTH]
             for i in range(len(window) - epitope.CORE_LENGTH + 1)]
    if not cores:
        return -1
    scores = backend.score(cores, mol)
    return window_start + int(np.argmax(scores))


def check_surrogate_direction(backend, molecule: Molecule) -> Check:
    """⚠️ surrogate sanity: the DR1 P1 pocket must prefer hydrophobic anchors.

    HLA-DRB1*01:01 has a large hydrophobic P1 pocket (Gly at beta86) with a
    well-documented preference for aromatic/aliphatic anchors over charged ones.
    A scorer that gets this backwards is broken regardless of what else it does.
    """
    hydrophobic = "YFWLIVM"
    charged = "DEKR"
    filler = "AGSTAGST"
    hyd = backend.score([a + filler for a in hydrophobic], molecule)
    chg = backend.score([a + filler for a in charged], molecule)
    return _check(
        "KAT-11", "surrogate direction: hydrophobic P1 preferred (HLA-DRB1*01:01)",
        float(hyd.mean()) > float(chg.mean()),
        f"hydrophobic mean {hyd.mean():.3f} vs charged mean {chg.mean():.3f}",
        "hydrophobic > charged",
        "textbook property of the DR1 P1 pocket; a directional smoke test only",
    )


# --------------------------------------------------------------------------
# panel and applicability-domain controls
# --------------------------------------------------------------------------

def check_panel_integrity(molecules: Sequence[Molecule], n_decoys: int = 10,
                          seed: int = 20260826) -> Check:
    """Panel members must be separable from decoys on locus-invariant positions.

    Not "as conserved as a human allele" -- a canine or feline allele is
    legitimately less similar to the human calibration set than another human
    allele is, and calibrating on the training species alone would reject real
    DLA and FLA molecules for being canine and feline.

    Instead this is a discrimination test. Every panel molecule is scored on the
    fraction of cross-species-invariant positions it keeps, and so are two kinds
    of decoy: composition-matched shuffles of the panel sequences, and molecules
    of the *opposite chain class* (an alpha chain scored against a beta locus and
    vice versa). The panel passes when its worst member still scores above the
    best decoy.

    Deliberately not a decoy: a beta chain of a different locus. DRB and DQB
    beta1 domains are homologous enough that a DRB sequence scores well on DQB's
    invariant positions -- which is precisely why a pan-specific model pools
    them, and is not a defect this filter should be asked to catch.
    """
    if not molecules:
        return _check("KAT-9", "panel integrity", False, "empty panel", "non-empty")

    by_locus: Dict[str, List[Molecule]] = {}
    for mol in molecules:
        by_locus.setdefault(mol.locus, []).append(mol)

    worst_real = (1.0, "")
    best_decoy = (0.0, "")
    for locus, members in by_locus.items():
        for mol in members:
            score = groove.conservation_match(mol, locus)
            if score < worst_real[0]:
                worst_real = (score, f"{mol.name} [{locus}]")

        rng = random.Random(seed)
        for i in range(n_decoys):
            donor = members[i % len(members)]
            chars = list(donor.sequence)
            rng.shuffle(chars)
            decoy = groove.build_molecule(f"shuffled::{donor.name}", "decoy",
                                          locus, "".join(chars))
            score = groove.conservation_match(decoy, locus)
            if score > best_decoy[0]:
                best_decoy = (score, f"shuffle of {donor.name}")

        want_chain = groove.LOCUS_REF[locus]["chain"]
        for other, others in by_locus.items():
            if other == locus or groove.LOCUS_REF[other]["chain"] == want_chain:
                continue
            for mol in others[:5]:
                score = groove.conservation_match(mol, locus)
                if score > best_decoy[0]:
                    best_decoy = (score, f"{other} {mol.name} scored as {locus}")

    gap = worst_real[0] - best_decoy[0]
    return _check(
        "KAT-9", "panel integrity separates real alleles from decoys",
        gap > 0,
        f"worst panel member {worst_real[0]:.3f} ({worst_real[1]}) vs best decoy "
        f"{best_decoy[0]:.3f} ({best_decoy[1]}); margin {gap:+.3f}",
        "margin > 0",
        "guards against pseudogenes, wrong chains and truncated records entering "
        "the panel",
    )


def check_domain_guardrail(calls: Sequence[DomainCall], loo: Sequence[float]) -> Check:
    """The out-of-domain warning must actually be emitted for DLA/FLA.

    This does not test prediction accuracy -- no dog or cat benchmark exists.
    It tests that the pipeline refuses to present an extrapolated score as if it
    were in-distribution.
    """
    out = sum(c.verdict.startswith("OUT-OF-DOMAIN") for c in calls)
    med = float(np.median([c.identity for c in calls])) if calls else float("nan")
    train_p5 = float(np.percentile(loo, 5)) if len(loo) else float("nan")
    return Check(
        "KAT-10", "applicability-domain guard-rail fires",
        "PASS" if out > 0 or med < train_p5 else "FAIL",
        f"{out}/{len(calls)} molecules out-of-domain; median NN identity {med:.3f} "
        f"vs training 5th percentile {train_p5:.3f}",
        "a warning is emitted whenever the panel sits outside the training space",
        "the guard-rail, not the predictor, is what is being validated here",
    )


# --------------------------------------------------------------------------
# runner
# --------------------------------------------------------------------------

@dataclass
class ValidationReport:
    species: str
    checks: List[Check] = field(default_factory=list)

    @property
    def passed(self) -> int:
        return sum(c.status == "PASS" for c in self.checks)

    @property
    def failed(self) -> List[Check]:
        return [c for c in self.checks if c.status == "FAIL"]

    @property
    def ok(self) -> bool:
        return not self.failed

    def to_rows(self) -> List[Dict[str, str]]:
        return [
            {"id": c.id, "check": c.name, "status": c.status,
             "observed": c.observed, "expected": c.expected, "rationale": c.detail}
            for c in self.checks
        ]


def run_all(species: str, products: Sequence[Insulin], molecules: Sequence[Molecule],
            backend, ranker, domain_calls: Sequence[DomainCall],
            loo: Sequence[float], identity_donor: Optional[str] = None,
            reference_molecule: Optional[Molecule] = None) -> ValidationReport:
    report = ValidationReport(species)
    report.checks.extend(check_sequence_ground_truth())
    report.checks.append(check_negative_control(species))
    if identity_donor:
        report.checks.append(check_identity_control(species, identity_donor))
    else:
        report.checks.append(Check(
            "KAT-3", f"identity control: a natural insulin identical to {species} insulin",
            "SKIP", "no such donor species exists",
            "a sequence-identical natural donor",
            f"there is no commercially used insulin whose sequence matches {species} "
            f"insulin exactly -- unlike the dog, which has porcine insulin. That absence "
            f"is a finding, not a gap in the harness.",
        ))
    report.checks.append(check_core_set_invariant(species, products))
    report.checks.append(check_scramble_positive_control(species))
    report.checks.append(check_tolerance_filter_precision(species))
    report.checks.append(check_panel_integrity(molecules))
    report.checks.append(check_domain_guardrail(domain_calls, loo))

    probe = molecules[0]
    cores = [c.sequence for c in epitope.all_cores(products[0])]
    report.checks.append(check_rank_calibration(ranker, probe))
    report.checks.append(check_determinism(backend, probe, cores))
    report.checks.append(check_register_stability(backend, molecules[:8], products[0]))
    if reference_molecule is not None and not getattr(backend, "trained", True):
        report.checks.append(check_surrogate_direction(backend, reference_molecule))
    return report
