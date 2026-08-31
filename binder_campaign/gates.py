"""The per-target validation gate: ``/state/gates/{target}.json``.

Production scoring on a target is blocked until this file exists with
``status == "PASS"``.  ``submit_gate()`` reads the file when ``target=<TARGET>``
is passed for a scoring job and returns False if it is absent or not PASS — a
verbal claim, a Slack post, or an in-memory variable does not satisfy the gate.

The file records, for each ranking instrument:

(a) **target fold recapitulation** against a named CA-RMSD threshold (any
    core-scoping disclosed); and
(b) **positive-control separation** — a known literature binder at full native
    stoichiometry scoring clearly above negative controls.

For GDF-8 the file MUST additionally record a GDF-11 arm.  If no literature
control exists after a comprehensive search, condition (b) is dropped and the
gate passes on (a) alone (the ``no_literature_control`` path, which
``instrument_realization.csv`` then reports as ``control_separation_value=NA``).
If controls exist but (b) fails on every instrument, self-authorisation is
barred: the three prescribed remedies must be attempted first.

The gate is also where the per-target instrument mask, the scoring construct
(chains, crop, cofactors) and the frozen thresholds are pinned; they are
identical at every seed tier thereafter.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

from .scoring import DEFAULT_ARMS, DEFAULT_POSE_DOCKQ_THRESHOLD, InstrumentMask

__all__ = [
    "ScoringConstruct",
    "InstrumentGateRow",
    "ValidationGate",
    "write_gate",
    "read_gate",
    "gate_is_pass",
]

#: The three remedies that must be tried before a failing separation check can
#: be waived.  Self-authorisation past (b) is barred.
SEPARATION_REMEDIES = (
    "rerun_control_with_control_ligand_msa",
    "alternative_cofolder",
    "af_unmasked_template_injection",
)


@dataclass(frozen=True)
class ScoringConstruct:
    """Frozen at gate time; identical at every seed tier.

    Asserted in the output structure file at the validation gate AND on every
    production scoring row: the scoring job parses its own output structure and
    compares asserted chain count and cofactor atoms against this record.  On
    mismatch the score is written NaN with ``construct_status=CONSTRUCT_FAIL``,
    never a numeric.
    """

    chains: tuple[str, ...]
    residue_range: str
    cofactors: tuple[str, ...] = ()
    n_target_chains: int = 1
    native_oligomer_n: int = 1
    crop: str | None = None
    reference_structure: str | None = None  # PDB id, from a lookup, never memory

    def matches(self, n_chains_folded: int, cofactor_atoms: int) -> bool:
        """Both counts must agree: cofactor atoms and target chain count."""
        if n_chains_folded != self.n_target_chains + 1:  # +1 for the binder
            return False
        if self.cofactors and cofactor_atoms <= 0:
            return False
        if not self.cofactors and cofactor_atoms > 0:
            return False
        return True


@dataclass(frozen=True)
class InstrumentGateRow:
    """One ranking instrument's gate evidence."""

    arm: str
    # (a) fold recapitulation
    ca_rmsd: float
    ca_rmsd_threshold: float
    core_scoping: str | None = None
    # (b) positive-control separation
    control_name: str | None = None
    control_is_antibody: bool = False
    control_score: float | None = None
    negative_control_scores: tuple[float, ...] = ()
    # construct realisation
    cofactors_present: tuple[str, ...] = ()
    n_target_chains_folded: int = 1
    remedies_attempted: tuple[str, ...] = ()

    @property
    def fold_pass(self) -> bool:
        """(a): recapitulates the target fold within the named CA-RMSD threshold."""
        return self.ca_rmsd <= self.ca_rmsd_threshold

    @property
    def separation_value(self) -> float | None:
        """control_score minus the best negative control; None when no control."""
        if self.control_score is None or not self.negative_control_scores:
            return None
        return self.control_score - max(self.negative_control_scores)

    @property
    def separation_pass(self) -> bool | None:
        """(b): ``None`` when there is no literature control to test."""
        sep = self.separation_value
        if sep is None:
            return None
        return sep > 0.0


@dataclass
class ValidationGate:
    """``/state/gates/{target}.json``."""

    target: str
    instruments: list[InstrumentGateRow]
    construct: ScoringConstruct
    mask: InstrumentMask = field(default_factory=InstrumentMask)
    #: frozen thresholds, identical at every seed tier from here on
    pose_dockq_threshold: float = DEFAULT_POSE_DOCKQ_THRESHOLD
    monomer_plddt_threshold: float = 70.0
    rank_zscore_formula: str = "weighted_z_4x_ipsae_1x_scdockq"
    final_score_formula: str = "raw_mean_of_mask_terms"
    #: set when a comprehensive search found no literature control
    no_literature_control: bool = False
    #: GDF-8 specifically must additionally record a GDF-11 arm
    counter_target: str | None = None
    counter_target_instruments: list[InstrumentGateRow] = field(default_factory=list)
    alternatives_tried: list[str] = field(default_factory=list)
    set_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    # ---- status ------------------------------------------------------------ #
    def blockers(self) -> list[str]:
        problems: list[str] = []
        if not self.instruments:
            problems.append("no ranking instrument recorded")

        for row in self.instruments:
            if not row.fold_pass:
                problems.append(
                    f"{row.arm}: fold recapitulation CA-RMSD {row.ca_rmsd:.2f} > "
                    f"{row.ca_rmsd_threshold:.2f}"
                )

        if not self.no_literature_control:
            seps = [r.separation_pass for r in self.instruments]
            tested = [s for s in seps if s is not None]
            if tested and not any(tested):
                # controls exist but (b) fails on every instrument
                missing = [
                    r for r in SEPARATION_REMEDIES
                    if not any(r in i.remedies_attempted for i in self.instruments)
                ]
                if missing:
                    problems.append(
                        "positive-control separation fails on every instrument and "
                        f"self-authorisation is barred; remedies not yet tried: {missing}"
                    )
                else:
                    problems.append(
                        "positive-control separation fails on every instrument after "
                        "all three remedies"
                    )
            # an antibody-only control that fails separation is weak evidence and
            # the gate passes on fold check (a) alone
            if tested and not any(tested):
                if all(r.control_is_antibody for r in self.instruments
                       if r.separation_pass is not None):
                    problems = [p for p in problems
                                if "positive-control separation" not in p]

        if self.counter_target and not self.counter_target_instruments:
            problems.append(
                f"{self.target} is counter-screened against {self.counter_target}: "
                "the gate file must record that arm"
            )
        return problems

    @property
    def status(self) -> str:
        return "PASS" if not self.blockers() else "FAIL"

    def to_json(self) -> dict:
        return {
            "target": self.target,
            "status": self.status,
            "set_at": self.set_at,
            "blockers": self.blockers(),
            "instruments": [asdict(r) for r in self.instruments],
            "counter_target": self.counter_target,
            "counter_target_instruments": [
                asdict(r) for r in self.counter_target_instruments
            ],
            "construct": asdict(self.construct),
            "frozen_mask": {
                "name": self.mask.describe(),
                "arms": list(self.mask.arms),
                "counter_screened": self.mask.counter_screened,
            },
            "frozen_thresholds": {
                "pose_dockq": self.pose_dockq_threshold,
                "monomer_plddt": self.monomer_plddt_threshold,
            },
            "rank_zscore_formula": self.rank_zscore_formula,
            "final_score_formula": self.final_score_formula,
            "no_literature_control": self.no_literature_control,
            "alternatives_tried": list(self.alternatives_tried),
        }


def write_gate(state_root: str, gate: ValidationGate) -> str:
    path = os.path.join(state_root, "gates", f"{gate.target}.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(gate.to_json(), fh, indent=2, sort_keys=True)
    return path


def read_gate(state_root: str, target: str) -> dict | None:
    path = os.path.join(state_root, "gates", f"{target}.json")
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def gate_is_pass(state_root: str, target: str) -> bool:
    gate = read_gate(state_root, target)
    return bool(gate and str(gate.get("status", "")).upper() == "PASS")


def default_instrument_rows(
    ca_rmsd: float, threshold: float, arms: tuple[str, ...] = DEFAULT_ARMS
) -> list[InstrumentGateRow]:
    return [
        InstrumentGateRow(arm=a, ca_rmsd=ca_rmsd, ca_rmsd_threshold=threshold)
        for a in arms
    ]
