"""Pre-scoring filters — the four gates every candidate clears before co-folding.

The prompt requires each candidate to be assessed, *before* any co-folding
spend, for

1. **novelty** — MMseqs2 vs UniRef90 plus the known-binder corpus plus every
   chain of the target's own campaign reference structure and positive-control
   complex.  REJECT at ``>60 %`` identity over ``>50 %`` coverage, OR at
   ``>=30 %`` gapped local identity over ``>=40`` aligned residues, OR
   ``TM-score >= 0.5`` to any target or control chain (so target-mimic
   protomers are caught here, not downstream);
2. **liability** — Cys parity, homopolymer runs, surface hydrophobic patches;
3. **monomer-foldability** — binder alone at a per-target mean-pLDDT threshold
   frozen at the validation gate (default 70 on the 0-100 scale / 0.7 on 0-1);
4. **structural-plausibility** — backbone geometry, steric clashes, core
   packing, at thresholds chosen at kickoff and frozen.

Every rejected ``design_id`` is recorded, and a gate counts as run only when its
rejects are traceably absent downstream — :func:`assert_rejects_absent`.

The homology search itself is an external tool (MMseqs2); it is injected as a
:class:`HomologySearcher`.  The *decision rules* on top of its hits, and the
gapped-local-identity check (a real Smith-Waterman, not a stub), live here so
the gate is reproducible and testable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Protocol, Sequence

__all__ = [
    "GateThresholds",
    "GateResult",
    "HomologyHit",
    "HomologySearcher",
    "liability_check",
    "novelty_check",
    "monomer_foldability_check",
    "structural_plausibility_check",
    "smith_waterman_identity",
    "run_prescoring_gates",
    "assert_rejects_absent",
]

HYDROPHOBIC = set("AVILMFWYC")


@dataclass(frozen=True)
class GateThresholds:
    """Chosen at kickoff and FROZEN; recorded on every design row."""

    # novelty
    max_global_identity: float = 0.60          # > 60 % identity ...
    max_global_identity_coverage: float = 0.50  # ... over > 50 % coverage
    max_local_identity: float = 0.30           # >= 30 % gapped local identity ...
    min_local_aligned: int = 40                # ... over >= 40 aligned residues
    max_tm_score_to_target: float = 0.50       # TM-score >= 0.5 to any chain
    # liability
    max_homopolymer_run: int = 4               # a run longer than this is flagged
    require_even_cys: bool = True
    max_hydrophobic_window_fraction: float = 0.75
    hydrophobic_window: int = 9
    # monomer foldability
    min_monomer_plddt: float = 70.0            # 0-100 scale
    # structural plausibility
    max_clashes: int = 0
    min_core_packing: float = 0.40
    max_ca_bond_deviation_a: float = 0.20
    # miniprotein length window
    min_binder_len: int = 50
    max_binder_len: int = 120
    permitted_min_len: int = 35                # where epitope geometry motivates
    permitted_max_len: int = 160

    def as_row(self) -> dict:
        return {f"thr_{k}": v for k, v in self.__dict__.items()}


@dataclass(frozen=True)
class GateResult:
    gate: str
    verdict: str  # "PASS" or "REJECT"
    reasons: tuple[str, ...] = ()
    metrics: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.verdict == "PASS"


@dataclass(frozen=True)
class HomologyHit:
    """One MMseqs2 (or equivalent) hit."""

    subject_id: str
    database: str          # "uniref90" | "binder_corpus" | "target_chain" | "control_chain"
    identity: float        # fraction, 0-1
    coverage: float        # fraction of the query covered, 0-1
    aligned_residues: int = 0


class HomologySearcher(Protocol):
    """The injected MMseqs2 wrapper."""

    def search(self, sequence: str) -> Sequence[HomologyHit]: ...


# --------------------------------------------------------------------------- #
# 2. liability
# --------------------------------------------------------------------------- #


def _longest_run(sequence: str) -> tuple[str, int]:
    best_aa, best = "", 0
    cur_aa, cur = "", 0
    for aa in sequence:
        if aa == cur_aa:
            cur += 1
        else:
            cur_aa, cur = aa, 1
        if cur > best:
            best_aa, best = cur_aa, cur
    return best_aa, best


def _max_hydrophobic_window_fraction(sequence: str, window: int) -> float:
    if len(sequence) < window:
        if not sequence:
            return 0.0
        return sum(aa in HYDROPHOBIC for aa in sequence) / len(sequence)
    flags = [1 if aa in HYDROPHOBIC else 0 for aa in sequence]
    run = sum(flags[:window])
    best = run
    for i in range(window, len(flags)):
        run += flags[i] - flags[i - window]
        best = max(best, run)
    return best / window


def liability_check(
    sequence: str, thresholds: GateThresholds = GateThresholds()
) -> GateResult:
    """Cys parity, homopolymer runs, surface hydrophobic patches."""
    seq = sequence.strip().upper()
    reasons: list[str] = []

    n_cys = seq.count("C")
    if thresholds.require_even_cys and n_cys % 2 == 1:
        reasons.append(f"unpaired cysteine (n_cys={n_cys})")

    run_aa, run_len = _longest_run(seq)
    if run_len > thresholds.max_homopolymer_run:
        reasons.append(f"homopolymer run {run_aa}x{run_len}")

    hyd = _max_hydrophobic_window_fraction(seq, thresholds.hydrophobic_window)
    if hyd > thresholds.max_hydrophobic_window_fraction:
        reasons.append(
            f"hydrophobic patch: {hyd:.2f} of a {thresholds.hydrophobic_window}-mer"
        )

    return GateResult(
        gate="liability",
        verdict="REJECT" if reasons else "PASS",
        reasons=tuple(reasons),
        metrics={
            "n_cys": n_cys,
            "longest_run_aa": run_aa,
            "longest_run_len": run_len,
            "max_hydrophobic_window_fraction": hyd,
        },
    )


# --------------------------------------------------------------------------- #
# 1. novelty
# --------------------------------------------------------------------------- #


def smith_waterman_identity(
    query: str,
    subject: str,
    match: int = 2,
    mismatch: int = -1,
    gap_open: int = -11,
    gap_extend: int = -1,
) -> tuple[float, int]:
    """Gapped local alignment; returns ``(identity_fraction, aligned_residues)``.

    Affine-gap Smith-Waterman with traceback.  Used for the prompt's second
    novelty rule (``>=30 %`` gapped local identity over ``>=40`` aligned
    residues) and, specifically, to catch Ubiquitin, which "often emerges with
    short terminal extensions, so detect by local alignment rather than exact
    match".
    """
    q, s = query.strip().upper(), subject.strip().upper()
    n, m = len(q), len(s)
    if n == 0 or m == 0:
        return 0.0, 0

    NEG = float("-inf")
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    E = [[NEG] * (m + 1) for _ in range(n + 1)]  # gap in query (move along s)
    F = [[NEG] * (m + 1) for _ in range(n + 1)]  # gap in subject (move along q)
    ptr = [[0] * (m + 1) for _ in range(n + 1)]  # 0 stop, 1 diag, 2 left(E), 3 up(F)

    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            E[i][j] = max(H[i][j - 1] + gap_open, E[i][j - 1] + gap_extend)
            F[i][j] = max(H[i - 1][j] + gap_open, F[i - 1][j] + gap_extend)
            diag = H[i - 1][j - 1] + (match if q[i - 1] == s[j - 1] else mismatch)
            cell = max(0.0, diag, E[i][j], F[i][j])
            H[i][j] = cell
            if cell == 0.0:
                ptr[i][j] = 0
            elif cell == diag:
                ptr[i][j] = 1
            elif cell == E[i][j]:
                ptr[i][j] = 2
            else:
                ptr[i][j] = 3
            if cell > best:
                best, bi, bj = cell, i, j

    if best == 0.0:
        return 0.0, 0

    identities = 0
    aligned = 0
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i][j] != 0:
        d = ptr[i][j]
        if d == 1:
            aligned += 1
            if q[i - 1] == s[j - 1]:
                identities += 1
            i, j = i - 1, j - 1
        elif d == 2:
            aligned += 1
            j -= 1
        else:
            aligned += 1
            i -= 1
    return (identities / aligned if aligned else 0.0), aligned


def novelty_check(
    sequence: str,
    hits: Iterable[HomologyHit],
    structural_tm_scores: dict[str, float] | None = None,
    control_sequences: dict[str, str] | None = None,
    thresholds: GateThresholds = GateThresholds(),
) -> GateResult:
    """Apply all three novelty rules.

    ``structural_tm_scores`` maps a target/control chain id to the TM-score of
    the binder against it; ``control_sequences`` are chains to run the gapped
    local-identity rule against directly (target chains, control chains, and the
    Ubiquitin sequences P0CG47/P0CG48 when staged).
    """
    reasons: list[str] = []
    metrics: dict = {}
    worst_identity = 0.0
    worst_local = (0.0, 0)

    for hit in hits:
        worst_identity = max(worst_identity, hit.identity)
        # rule 1: > 60 % identity over > 50 % coverage
        if (
            hit.identity > thresholds.max_global_identity
            and hit.coverage > thresholds.max_global_identity_coverage
        ):
            reasons.append(
                f"{hit.database}:{hit.subject_id} identity {hit.identity:.2f} "
                f"over coverage {hit.coverage:.2f}"
            )
        # rule 2 as reported by the search, when it gives aligned counts
        if (
            hit.aligned_residues >= thresholds.min_local_aligned
            and hit.identity >= thresholds.max_local_identity
        ):
            reasons.append(
                f"{hit.database}:{hit.subject_id} local identity "
                f"{hit.identity:.2f} over {hit.aligned_residues} aligned residues"
            )

    # rule 2 computed directly against the chains we hold
    for name, subject in (control_sequences or {}).items():
        ident, aligned = smith_waterman_identity(sequence, subject)
        if (ident, aligned) > worst_local:
            worst_local = (ident, aligned)
        if (
            aligned >= thresholds.min_local_aligned
            and ident >= thresholds.max_local_identity
        ):
            reasons.append(
                f"local identity {ident:.2f} over {aligned} aligned residues to {name}"
            )

    # rule 3: TM-score >= 0.5 to any target or control chain
    for chain, tm in (structural_tm_scores or {}).items():
        if tm >= thresholds.max_tm_score_to_target:
            reasons.append(f"TM-score {tm:.2f} to {chain} (target-mimic)")

    metrics.update(
        max_identity=worst_identity,
        max_local_identity=worst_local[0],
        max_local_aligned=worst_local[1],
        max_tm_score=max((structural_tm_scores or {}).values(), default=0.0),
    )
    return GateResult(
        gate="novelty",
        verdict="REJECT" if reasons else "PASS",
        reasons=tuple(reasons),
        metrics=metrics,
    )


# --------------------------------------------------------------------------- #
# 3 & 4. foldability and plausibility
# --------------------------------------------------------------------------- #


def monomer_foldability_check(
    mean_plddt: float, thresholds: GateThresholds = GateThresholds()
) -> GateResult:
    """Binder alone at the per-target mean-pLDDT threshold frozen at the gate.

    ESMFold2 emits pLDDT on 0-1; anything at or below 1.0 is rescaled to 0-100
    so the frozen threshold means the same thing on either scale.
    """
    plddt = mean_plddt * 100.0 if mean_plddt <= 1.0 else mean_plddt
    ok = plddt >= thresholds.min_monomer_plddt
    return GateResult(
        gate="monomer_foldability",
        verdict="PASS" if ok else "REJECT",
        reasons=() if ok else (
            f"mean pLDDT {plddt:.1f} < {thresholds.min_monomer_plddt:.1f}",
        ),
        metrics={"monomer_plddt": plddt},
    )


def structural_plausibility_check(
    n_clashes: int,
    core_packing: float,
    max_bond_deviation_a: float,
    binder_len: int,
    thresholds: GateThresholds = GateThresholds(),
    length_rationale: str | None = None,
) -> GateResult:
    """Backbone geometry, steric clashes, core packing, and the length window."""
    reasons: list[str] = []
    if n_clashes > thresholds.max_clashes:
        reasons.append(f"{n_clashes} steric clashes")
    if core_packing < thresholds.min_core_packing:
        reasons.append(f"core packing {core_packing:.2f} < {thresholds.min_core_packing}")
    if max_bond_deviation_a > thresholds.max_ca_bond_deviation_a:
        reasons.append(f"backbone bond deviation {max_bond_deviation_a:.2f} A")

    if not (thresholds.min_binder_len <= binder_len <= thresholds.max_binder_len):
        if thresholds.permitted_min_len <= binder_len <= thresholds.permitted_max_len:
            if not length_rationale:
                reasons.append(
                    f"binder_len {binder_len} outside 50-120 without a recorded "
                    "rationale in the design sheet target metadata"
                )
        else:
            reasons.append(f"binder_len {binder_len} outside the permitted 35-160")

    return GateResult(
        gate="structural_plausibility",
        verdict="REJECT" if reasons else "PASS",
        reasons=tuple(reasons),
        metrics={
            "n_clashes": n_clashes,
            "core_packing": core_packing,
            "max_bond_deviation_a": max_bond_deviation_a,
            "binder_len": binder_len,
        },
    )


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #


def run_prescoring_gates(
    design_id: str,
    sequence: str,
    *,
    hits: Iterable[HomologyHit] = (),
    structural_tm_scores: dict[str, float] | None = None,
    control_sequences: dict[str, str] | None = None,
    mean_plddt: float | None = None,
    n_clashes: int = 0,
    core_packing: float = 1.0,
    max_bond_deviation_a: float = 0.0,
    length_rationale: str | None = None,
    thresholds: GateThresholds = GateThresholds(),
) -> dict:
    """Run all four gates and return a traceable verdict record.

    The record is what gets written to ``novelty_verdict_path`` and what the
    sheet writer recomputes at write time.
    """
    results = [
        novelty_check(sequence, hits, structural_tm_scores, control_sequences,
                      thresholds),
        liability_check(sequence, thresholds),
    ]
    if mean_plddt is not None:
        results.append(monomer_foldability_check(mean_plddt, thresholds))
    results.append(
        structural_plausibility_check(
            n_clashes, core_packing, max_bond_deviation_a, len(sequence.strip()),
            thresholds, length_rationale,
        )
    )

    rejected_by = [r.gate for r in results if not r.passed]
    return {
        "design_id": design_id,
        "verdict": "REJECT" if rejected_by else "PASS",
        "rejected_by": rejected_by,
        "gates": {
            r.gate: {"verdict": r.verdict, "reasons": list(r.reasons),
                     "metrics": r.metrics}
            for r in results
        },
        "thresholds": thresholds.as_row(),
    }


def assert_rejects_absent(
    rejected_design_ids: Iterable[str], downstream_design_ids: Iterable[str]
) -> None:
    """A gate counts as run only when its rejects are traceably absent downstream."""
    leaked = sorted(set(rejected_design_ids) & set(downstream_design_ids))
    if leaked:
        raise AssertionError(
            f"{len(leaked)} gate-rejected designs reached a downstream scoring "
            f"pool: {leaked[:10]}{'...' if len(leaked) > 10 else ''}"
        )
