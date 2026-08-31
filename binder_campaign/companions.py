"""Design-sheet companion files.

* ``per_seed_metrics.parquet`` — one row per ``design_id x scoring_arm x seed``
  with the raw metric values, so any sheet score can be reproduced without
  re-compute.  Companion ``design_id``s cover 100 % of the ranked sheet's
  ``design_id``s (exact, case-sensitive match), and for counter-screened targets
  the off-target arm at the same seed count.  No ``(design_id, scoring_arm)``
  group is uniformly null on the ranking metric.
* ``instrument_realization.csv`` — one row per ``(target, ranking_arm)``.  Every
  numeric is **derived** from the companion and ``/state/gates/{target}.json``
  at write time, never hand-authored, and the writer asserts row-for-row
  agreement with a recount over the companion.

Close-out is transactional: freeze the ranked sheet, derive every companion from
*that* snapshot, and assert the mandatory gate columns are non-null.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from .scoring import ArmAggregate, DesignScore, InstrumentMask, SeedRecord

__all__ = [
    "build_per_seed_metrics",
    "assert_companion_coverage",
    "recompute_final_scores_from_companion",
    "build_instrument_realization",
    "InstrumentRealizationRow",
]


# --------------------------------------------------------------------------- #
# per_seed_metrics
# --------------------------------------------------------------------------- #


def build_per_seed_metrics(
    designs: Iterable[DesignScore], mask: InstrumentMask
) -> pd.DataFrame:
    """One row per design x arm x seed, on-target and (when counter-screened)
    off-target at the same seed count."""
    rows: list[dict] = []
    for d in designs:
        for side, book in (("on_target", d.on_target), ("off_target", d.off_target)):
            for arm, agg in book.items():
                rows.append(_agg_row(d, side, arm, agg))
    return pd.DataFrame(rows)


def _agg_row(d: DesignScore, side: str, arm: str, agg: ArmAggregate) -> dict:
    return {
        "design_id": d.design_id,
        "target": d.target,
        "side": side,
        "scoring_arm": arm,
        "n_seeds": agg.n_seeds,
        "seeds": ",".join(str(s) for s in agg.seeds),
        "ipsae_min": agg.ipsae_min,
        "sc_dockq": agg.sc_dockq,
        "argmax_seed_ipsae": agg.argmax_seed_ipsae,
        "argmax_seed_dockq": agg.argmax_seed_dockq,
        "structure_path": agg.seedbest_structure_path,
    }


def build_per_seed_metrics_from_seeds(
    seed_book: Mapping[tuple[str, str, str], Sequence[SeedRecord]],
    targets: Mapping[str, str],
) -> pd.DataFrame:
    """Expanded form: one row per *individual seed*, keyed
    ``(design_id, side, arm) -> [SeedRecord, ...]``."""
    rows: list[dict] = []
    for (design_id, side, arm), records in seed_book.items():
        for r in records:
            rows.append({
                "design_id": design_id,
                "target": targets.get(design_id, ""),
                "side": side,
                "scoring_arm": arm,
                "seed": r.seed,
                "ipsae_ab": r.ipsae_ab,
                "ipsae_ba": r.ipsae_ba,
                "ipsae_min": r.ipsae_min,
                "dockq": r.dockq,
                "iptm": r.iptm,
                "lis": r.lis,
                "structure_path": r.structure_path,
            })
    return pd.DataFrame(rows)


def assert_companion_coverage(
    sheet: pd.DataFrame, companion: pd.DataFrame, mask: InstrumentMask
) -> None:
    """100 % coverage, exact and case-sensitive; no uniformly-null arm group."""
    sheet_ids = set(sheet["design_id"].astype(str))
    comp_ids = set(companion["design_id"].astype(str))
    missing = sheet_ids - comp_ids
    if missing:
        raise AssertionError(
            f"companion covers {len(comp_ids & sheet_ids)}/{len(sheet_ids)} ranked "
            f"design_ids; missing {sorted(missing)[:10]} — lower coverage is a "
            "deliverable defect"
        )

    on = companion[companion["side"] == "on_target"]
    for (design_id, arm), grp in on.groupby(["design_id", "scoring_arm"]):
        if design_id not in sheet_ids:
            continue
        if grp["ipsae_min"].isna().all():
            raise AssertionError(
                f"({design_id}, {arm}) is uniformly null on the ranking metric"
            )

    if mask.counter_screened:
        off = companion[companion["side"] == "off_target"]
        for design_id in sheet_ids:
            on_seeds = on[on["design_id"] == design_id]["n_seeds"]
            off_seeds = off[off["design_id"] == design_id]["n_seeds"]
            if off_seeds.empty:
                raise AssertionError(
                    f"{design_id}: counter-screened target has no off-target arm"
                )
            if not on_seeds.empty and set(on_seeds) != set(off_seeds):
                raise AssertionError(
                    f"{design_id}: off-target arm is not seed-matched to on-target"
                )


def recompute_final_scores_from_companion(
    sheet: pd.DataFrame,
    companion: pd.DataFrame,
    mask: InstrumentMask,
    n_sample: int = 30,
    tolerance: float = 1e-4,
    rng: random.Random | None = None,
) -> pd.DataFrame:
    """Recompute ``final_score`` for a sample of rows from the companion alone.

    Recomputed values must agree with the sheet to within ``1e-4`` or the
    companion is regenerated.  Returns the comparison frame so the check is
    reportable, and raises on disagreement.
    """
    rng = rng or random.Random(0)
    out: list[dict] = []
    for target, tgt_rows in sheet.groupby("target"):
        ids = list(tgt_rows["design_id"].astype(str))
        sample = ids if len(ids) <= n_sample else rng.sample(ids, n_sample)
        for design_id in sample:
            terms = _companion_terms(companion, design_id, mask)
            live = [t for t in terms if not math.isnan(t)]
            recomputed = sum(live) / len(live) if live else float("nan")
            carried = float(
                tgt_rows.loc[tgt_rows["design_id"] == design_id, "final_score"].iloc[0]
            )
            out.append({
                "target": target, "design_id": design_id,
                "sheet_final_score": carried,
                "companion_final_score": recomputed,
                "abs_diff": abs(carried - recomputed),
            })
    frame = pd.DataFrame(out)
    bad = frame[frame["abs_diff"] > tolerance]
    if not bad.empty:
        raise AssertionError(
            f"{len(bad)} sampled rows do not reproduce from the companion within "
            f"{tolerance}: {bad.head().to_dict('records')}"
        )
    return frame


def _companion_terms(
    companion: pd.DataFrame, design_id: str, mask: InstrumentMask
) -> list[float]:
    rows = companion[
        (companion["design_id"] == design_id) & (companion["side"] == "on_target")
    ]
    off = companion[
        (companion["design_id"] == design_id) & (companion["side"] == "off_target")
    ]
    terms: list[float] = []
    for arm in mask.arms:
        r = rows[rows["scoring_arm"] == arm]
        ipsae = float(r["ipsae_min"].iloc[0]) if len(r) else float("nan")
        dockq = float(r["sc_dockq"].iloc[0]) if len(r) else float("nan")
        terms += [ipsae, dockq]
        if mask.counter_screened:
            o = off[off["scoring_arm"] == arm]
            off_ipsae = float(o["ipsae_min"].iloc[0]) if len(o) else float("nan")
            terms.append(ipsae - off_ipsae)
    return terms


# --------------------------------------------------------------------------- #
# instrument_realization
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class InstrumentRealizationRow:
    target: str
    arm_name: str
    gate_status: str
    n_seeds_run: int
    n_ranked_with_score: int
    used_in_final_score: bool
    control_separation_value: Any
    drop_reason: str
    cofactors_present: str
    n_target_chains_folded: int
    native_oligomer_n: int


def build_instrument_realization(
    sheet: pd.DataFrame,
    companion: pd.DataFrame,
    gates: Mapping[str, Mapping[str, Any]],
    masks: Mapping[str, InstrumentMask],
) -> pd.DataFrame:
    """Derive the table from the companion and the gate files; then recount.

    Asserts: ``control_separation_value`` non-null on every ``gate_status=PASS``
    row except rows whose gate file records the no-literature-control path
    (those write ``NA`` and name that path in ``drop_reason``); ``drop_reason``
    non-empty wherever ``used_in_final_score`` is false.
    """
    rows: list[InstrumentRealizationRow] = []

    for target, gate in gates.items():
        mask = masks.get(target, InstrumentMask())
        sheet_t = sheet[sheet["target"] == target]
        comp_t = companion[
            (companion["target"] == target) & (companion["side"] == "on_target")
        ]
        construct = gate.get("construct", {})
        no_control = bool(gate.get("no_literature_control"))

        by_arm = {r["arm"]: r for r in gate.get("instruments", [])}
        all_arms = sorted(set(by_arm) | set(mask.arms) |
                          set(comp_t["scoring_arm"].unique()))

        for arm in all_arms:
            grp = comp_t[comp_t["scoring_arm"] == arm]
            n_seeds_run = int(grp["n_seeds"].max()) if len(grp) else 0
            scored_ids = set(grp.loc[grp["ipsae_min"].notna(), "design_id"])
            n_with_score = len(set(sheet_t["design_id"].astype(str)) & scored_ids)

            grow = by_arm.get(arm)
            used = arm in mask.arms
            gate_status = "PASS" if (grow and _fold_pass(grow)) else "FAIL"
            if grow is None:
                gate_status = "NOT_RUN"

            sep: Any = "NA"
            if grow and not no_control:
                sep = _separation(grow)

            drop_reason = ""
            if not used:
                drop_reason = (
                    f"arm not in the frozen mask {mask.describe()!r}"
                    if gate_status == "PASS"
                    else f"gate_status={gate_status}"
                )
            if no_control:
                drop_reason = (drop_reason + "; " if drop_reason else "") + \
                    "no_literature_control path (separation check (b) dropped)"

            if gate_status == "PASS" and (sep is None) and not no_control:
                raise AssertionError(
                    f"{target}/{arm}: gate_status=PASS but control_separation_value "
                    "is null and the gate does not record the no-literature-control "
                    "path"
                )
            if not used and not drop_reason:
                raise AssertionError(
                    f"{target}/{arm}: used_in_final_score=false with empty drop_reason"
                )

            rows.append(InstrumentRealizationRow(
                target=target,
                arm_name=arm,
                gate_status=gate_status,
                n_seeds_run=n_seeds_run,
                n_ranked_with_score=n_with_score,
                used_in_final_score=used,
                control_separation_value=sep,
                drop_reason=drop_reason,
                cofactors_present=";".join(
                    (grow or {}).get("cofactors_present", ())
                    or construct.get("cofactors", ())
                ),
                n_target_chains_folded=int(
                    (grow or {}).get("n_target_chains_folded",
                                     construct.get("n_target_chains", 1))
                ),
                native_oligomer_n=int(construct.get("native_oligomer_n", 1)),
            ))

    frame = pd.DataFrame([r.__dict__ for r in rows])
    _assert_recount(frame, sheet, companion)
    return frame


def _fold_pass(grow: Mapping[str, Any]) -> bool:
    return float(grow.get("ca_rmsd", 1e9)) <= float(grow.get("ca_rmsd_threshold", 0))


def _separation(grow: Mapping[str, Any]) -> Any:
    negs = grow.get("negative_control_scores") or []
    ctrl = grow.get("control_score")
    if ctrl is None or not negs:
        return None
    return float(ctrl) - max(float(n) for n in negs)


def _assert_recount(
    frame: pd.DataFrame, sheet: pd.DataFrame, companion: pd.DataFrame
) -> None:
    """Row-for-row agreement with an independent recount over the companion."""
    on = companion[companion["side"] == "on_target"]
    for _, row in frame.iterrows():
        grp = on[(on["target"] == row["target"]) &
                 (on["scoring_arm"] == row["arm_name"])]
        scored = set(grp.loc[grp["ipsae_min"].notna(), "design_id"])
        ranked = set(
            sheet.loc[sheet["target"] == row["target"], "design_id"].astype(str)
        )
        want = len(ranked & scored)
        if int(row["n_ranked_with_score"]) != want:
            raise AssertionError(
                f"{row['target']}/{row['arm_name']}: n_ranked_with_score "
                f"{row['n_ranked_with_score']} != recount {want}"
            )
