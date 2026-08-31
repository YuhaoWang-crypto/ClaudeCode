"""The campaign scoring instrument.

Implements the prompt's scoring section as executable code:

* per-arm ``ipSAE_min`` = min over both alignment directions (Dunbrack 2025),
  then **max over seeds**;
* ``sc_DockQ_{arm}`` taken on that arm's ``argmax_seed(ipSAE_min)`` seed, not on
  a separately chosen max-DockQ seed;
* ``pose_dockq`` = min over the arms that ran, ``pose_PASS`` at the per-target
  threshold (default 0.23);
* ``final_score`` = raw mean of the six terms (three ``ipSAE_min``, three
  ``sc_DockQ``), the headline value; under a disclosed ``REDUCED_MASK`` the mean
  runs over the realized terms only;
* ``rank_zscore`` = per-target weighted z-score average of the same terms, each
  ``ipSAE_min`` z-term weighted 4 and each ``sc_DockQ`` z-term weighted 1, used
  for ranking only;
* for counter-screened targets both formulas extend from six to nine terms with
  ``selectivity_delta_{arm} = ipSAE_min(on-target) - ipSAE_min(off-target)``,
  the z-term carrying the ipSAE weight of 4.

z-scores are transductive: mu and sigma come from the scored pool, so the raw
term vector is carried on every row and any batch can be re-standardised.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np

__all__ = [
    "DEFAULT_ARMS",
    "DEFAULT_MASK_NAME",
    "DEFAULT_POSE_DOCKQ_THRESHOLD",
    "IPSAE_WEIGHT",
    "DOCKQ_WEIGHT",
    "SeedRecord",
    "ArmAggregate",
    "InstrumentMask",
    "DesignScore",
    "aggregate_arm",
    "score_pool",
]

#: The campaign scoring default (the instrument): three arms, one sample per seed.
DEFAULT_ARMS: tuple[str, ...] = ("ef2full", "ef2fast", "ptxv2")
DEFAULT_MASK_NAME = "default_3arm"

#: pose_PASS threshold unless the validation gate freezes a stricter value.
DEFAULT_POSE_DOCKQ_THRESHOLD = 0.23

IPSAE_WEIGHT = 4.0
DOCKQ_WEIGHT = 1.0


# --------------------------------------------------------------------------- #
# per-seed records and per-arm aggregation
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SeedRecord:
    """One design x scoring arm x seed row of ``per_seed_metrics``.

    ``ipsae_ab`` / ``ipsae_ba`` are the two alignment directions; ipSAE for the
    seed is their minimum.  ``dockq`` is the self-consistency DockQ of this
    seed's predicted complex against the designed complex, under the best
    symmetric chain relabelling.
    """

    seed: int
    ipsae_ab: float
    ipsae_ba: float
    dockq: float
    structure_path: str | None = None
    iptm: float | None = None
    lis: float | None = None

    @property
    def ipsae_min(self) -> float:
        """min over both alignment directions, per Dunbrack 2025."""
        return min(self.ipsae_ab, self.ipsae_ba)


@dataclass(frozen=True)
class ArmAggregate:
    """One design x arm summary: what the sheet's per-arm columns carry."""

    arm: str
    ipsae_min: float
    sc_dockq: float
    n_seeds: int
    seeds: tuple[int, ...]
    argmax_seed_ipsae: int
    argmax_seed_dockq: int
    seedbest_structure_path: str | None = None

    @property
    def seeds_are_distinct(self) -> bool:
        return len(set(self.seeds)) == len(self.seeds)


def aggregate_arm(arm: str, records: Sequence[SeedRecord]) -> ArmAggregate:
    """Collapse a design's per-seed rows for one arm into its sheet values.

    ``ipSAE_min`` is the max over seeds of the per-seed two-direction minimum;
    ``sc_DockQ`` is read off the *same* seed (``argmax_seed(ipSAE_min)``), while
    ``argmax_seed(DockQ)`` is recorded separately so seed concordance stays
    auditable.
    """
    if not records:
        raise ValueError(f"arm {arm!r}: no seed records to aggregate")

    seeds = tuple(r.seed for r in records)
    if len(set(seeds)) != len(seeds):
        raise ValueError(
            f"arm {arm!r}: seeds must be DISTINCT integers, got {seeds}"
        )

    best_ipsae = max(records, key=lambda r: r.ipsae_min)
    best_dockq = max(records, key=lambda r: r.dockq)
    return ArmAggregate(
        arm=arm,
        ipsae_min=best_ipsae.ipsae_min,
        sc_dockq=best_ipsae.dockq,  # same seed as the ipSAE that became the score
        n_seeds=len(records),
        seeds=seeds,
        argmax_seed_ipsae=best_ipsae.seed,
        argmax_seed_dockq=best_dockq.seed,
        seedbest_structure_path=best_ipsae.structure_path,
    )


# --------------------------------------------------------------------------- #
# frozen instrument mask
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class InstrumentMask:
    """The per-target instrument mask, FROZEN at validation-gate time.

    ``arms`` are the arms whose terms enter ``final_score`` / ``rank_zscore``.
    ``substitutions`` records which arm each substitute stood in for.  A mask
    with fewer than three arms must carry a ``reduced`` name so a reader can see
    from ``score_instrument`` which arms produced the live ranking score.
    """

    name: str = DEFAULT_MASK_NAME
    arms: tuple[str, ...] = DEFAULT_ARMS
    counter_screened: bool = False
    pose_dockq_threshold: float = DEFAULT_POSE_DOCKQ_THRESHOLD
    substitutions: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.arms:
            raise ValueError("instrument mask must name at least one arm")
        if len(set(self.arms)) != len(self.arms):
            raise ValueError(f"duplicate arms in mask: {self.arms}")

    @property
    def reduced(self) -> bool:
        return len(self.arms) < len(DEFAULT_ARMS)

    @property
    def n_terms(self) -> int:
        """Six terms on the default mask; nine when counter-screened."""
        per_arm = 3 if self.counter_screened else 2
        return per_arm * len(self.arms)

    def describe(self) -> str:
        """The string written into the sheet's ``score_instrument`` column."""
        if self.name != DEFAULT_MASK_NAME:
            return self.name
        return "_".join(self.arms) if self.reduced else DEFAULT_MASK_NAME


# --------------------------------------------------------------------------- #
# per-design scores
# --------------------------------------------------------------------------- #


@dataclass
class DesignScore:
    """Everything the instrument computes for one design on one target."""

    design_id: str
    target: str
    on_target: dict[str, ArmAggregate]
    off_target: dict[str, ArmAggregate] = field(default_factory=dict)

    # ---- per-arm term accessors -------------------------------------------
    def ipsae(self, arm: str) -> float:
        agg = self.on_target.get(arm)
        return agg.ipsae_min if agg else float("nan")

    def sc_dockq(self, arm: str) -> float:
        agg = self.on_target.get(arm)
        return agg.sc_dockq if agg else float("nan")

    def selectivity_delta(self, arm: str) -> float:
        """``ipSAE_min(on-target) - ipSAE_min(off-target paralog)`` for ``arm``.

        Both sides must have been computed under the same seed count and
        aggregation; an instrument-mismatched subtraction is invalid and raises.
        """
        on = self.on_target.get(arm)
        off = self.off_target.get(arm)
        if on is None or off is None:
            return float("nan")
        if on.n_seeds != off.n_seeds:
            raise ValueError(
                f"{self.design_id}/{arm}: selectivity_delta needs matched seed "
                f"counts, got on={on.n_seeds} off={off.n_seeds}"
            )
        return on.ipsae_min - off.ipsae_min

    def n_seeds(self, arm: str) -> int:
        agg = self.on_target.get(arm)
        return agg.n_seeds if agg else 0

    def min_n_seeds(self, mask: InstrumentMask) -> int:
        return min((self.n_seeds(a) for a in mask.arms), default=0)

    # ---- composite values --------------------------------------------------
    def realized_terms(self, mask: InstrumentMask) -> dict[str, float]:
        """The raw terms of the frozen mask that actually have a value.

        Under a disclosed REDUCED_MASK the mean runs over these only.  A NaN
        from an arm that ran but failed is *not* silently dropped here: arms are
        dropped by re-freezing the mask, never per row.
        """
        terms: dict[str, float] = {}
        for arm in mask.arms:
            terms[f"ipsae_{arm}"] = self.ipsae(arm)
            terms[f"sc_DockQ_{arm}"] = self.sc_dockq(arm)
            if mask.counter_screened:
                terms[f"selectivity_delta_{arm}"] = self.selectivity_delta(arm)
        return terms

    def final_score(self, mask: InstrumentMask) -> float:
        """Headline value: the raw mean of the mask's terms."""
        vals = [v for v in self.realized_terms(mask).values() if not math.isnan(v)]
        if not vals:
            return float("nan")
        return float(np.mean(vals))

    def pose_dockq(self, mask: InstrumentMask) -> float:
        """min over the arms that ran; NaN is a missing measurement, not a flag."""
        vals = [self.sc_dockq(a) for a in mask.arms if a in self.on_target]
        if not vals or any(math.isnan(v) for v in vals):
            return float("nan")
        return float(min(vals))

    def pose_pass(self, mask: InstrumentMask) -> bool | None:
        pd = self.pose_dockq(mask)
        if math.isnan(pd):
            return None
        return bool(pd >= mask.pose_dockq_threshold)


# --------------------------------------------------------------------------- #
# transductive z-scoring over a pool
# --------------------------------------------------------------------------- #


def _zscore(values: np.ndarray) -> np.ndarray:
    """Population z-score, NaN-safe; a zero-variance term contributes 0."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros_like(values)
    mu = float(finite.mean())
    sigma = float(finite.std(ddof=0))
    if sigma == 0.0:
        return np.where(np.isfinite(values), 0.0, np.nan)
    return (values - mu) / sigma


def term_weights(mask: InstrumentMask) -> dict[str, float]:
    """4 for every ipSAE term (and every selectivity term), 1 for every DockQ."""
    w: dict[str, float] = {}
    for arm in mask.arms:
        w[f"ipsae_{arm}"] = IPSAE_WEIGHT
        w[f"sc_DockQ_{arm}"] = DOCKQ_WEIGHT
        if mask.counter_screened:
            w[f"selectivity_delta_{arm}"] = IPSAE_WEIGHT
    return w


def score_pool(
    designs: Iterable[DesignScore],
    mask: InstrumentMask,
) -> list[dict]:
    """Compute ``final_score``, ``rank_zscore``, pose columns for a whole pool.

    z-scoring is done here, over exactly this pool, because mu and sigma depend
    on it: the returned rows carry the raw term vector so any other batch can
    re-standardise.  Ranking uses ``rank_zscore``; reporting uses ``final_score``.
    """
    designs = list(designs)
    if not designs:
        return []

    weights = term_weights(mask)
    term_names = list(weights)

    raw = {
        name: np.array([d.realized_terms(mask)[name] for d in designs], dtype=float)
        for name in term_names
    }
    z = {name: _zscore(raw[name]) for name in term_names}

    rows: list[dict] = []
    for i, d in enumerate(designs):
        num = 0.0
        den = 0.0
        for name in term_names:
            zi = z[name][i]
            if math.isnan(zi):
                continue
            num += weights[name] * zi
            den += weights[name]
        rank_z = num / den if den else float("nan")

        row: dict = {
            "design_id": d.design_id,
            "target": d.target,
            "final_score": d.final_score(mask),
            "rank_zscore": rank_z,
            "score_instrument": mask.describe(),
            "pose_dockq": d.pose_dockq(mask),
            "pose_PASS": d.pose_pass(mask),
            "n_seeds": d.min_n_seeds(mask),
        }
        for name in term_names:
            row[name] = raw[name][i]
        for arm in mask.arms:
            agg = d.on_target.get(arm)
            row[f"n_seeds_{arm}"] = agg.n_seeds if agg else 0
            row[f"predicted_structure_path_{arm}"] = (
                agg.seedbest_structure_path if agg else None
            )
            if mask.counter_screened:
                off = d.off_target.get(arm)
                row[f"ipsae_offtarget_{arm}"] = off.ipsae_min if off else float("nan")
        rows.append(row)
    return rows
