"""The design-sheet writer: selection caps, rank assignment, relaxation ladder.

The prompt makes the sheet writer the mechanical enforcement point for most of
the campaign's invariants, so this module is where they live:

**Selection caps**, applied in selection order:

===  =======================================================================
(a)  reject exact-sequence duplicates
(b)  reject pairs within Levenshtein distance five
(c)  cap any single ``root_backbone_id`` at 5 % of rows (rounded up)
(d)  cap any single TM-score-0.90 single-linkage cluster at 10 % (rounded up)
(e)  cap any single ``structure_method`` at 50 % (max 15 of 30) and require at
     least three distinct ``structure_method`` values
(f)  cap any single ``seq_method`` at two-thirds, backfilling from the next-best
     alternate ``seq_method``; relaxed only as a last resort
===  =======================================================================

These are ceilings, not quotas.

**Rank** is assigned by ``ORDER BY (n_seeds >= 5 on all arms) DESC,
pose_PASS DESC, rank_zscore DESC`` — the only rank-assignment path.  Rows with
``n_seeds < 5`` after a failed top-up stay *ranked*, below the 5-seed tier, with
their true seed count disclosed; the unranked section is reserved for rows with
missing or zero scores and for rows missing a required provenance column.

**Relaxation ladder**, one step at a time, recorded on every affected row:
(i) diversity caps — never past per-root 25 %, per-method 50 %, or fewer than 3
distinct structure methods; (ii) liability flags; and, only after regeneration
and the full ladder have failed, ``NOVELTY_LAST_RESORT``.  ``pose_dockq``
non-null and ``final_score`` non-null are never relaxed, and padding a short
sheet with duplicates is forbidden — ship the actual N.

**Gate recompute**: every gate is recomputed from the row's own sequence and
predicted structure at write time and must match the carried value to within
1e-4; a mismatch halts the writer with the row id.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping, Sequence

from .filters import GateThresholds, liability_check, novelty_check
from .lcp import lcp_score
from .schema import MethodVocab, SheetSchema, default_method_vocab
from .scoring import InstrumentMask

__all__ = [
    "ROWS_PER_TARGET",
    "SelectionCaps",
    "SheetResult",
    "levenshtein_within",
    "select_and_rank",
    "fold_diversity_report",
]

ROWS_PER_TARGET = 30
TOLERANCE = 1e-4

#: The ladder, in the prompt's fixed order.
RELAXATION_LADDER = ("DIVERSITY_CAPS", "LIABILITY_FLAGS", "NOVELTY_LAST_RESORT")


@dataclass(frozen=True)
class SelectionCaps:
    """Caps (c)-(f) plus the Levenshtein floor, and their relaxation bounds."""

    max_root_backbone_frac: float = 0.05
    max_tm90_cluster_frac: float = 0.10
    max_structure_method_frac: float = 0.50
    max_seq_method_frac: float = 2.0 / 3.0
    min_distinct_structure_methods: int = 3
    min_levenshtein: int = 5  # reject pairs *within* distance five

    # relaxation bounds - the ladder may never go past these
    relaxed_max_root_backbone_frac: float = 0.25
    relaxed_max_tm90_cluster_frac: float = 0.25
    relaxed_max_structure_method_frac: float = 0.50

    def cap(self, frac: float, n_rows: int) -> int:
        """A cap of ``frac`` of ``n_rows``, rounded up."""
        return max(1, math.ceil(frac * n_rows))

    def relaxed(self) -> "SelectionCaps":
        """Ladder step (i): loosen the diversity caps to their hard bounds."""
        return SelectionCaps(
            max_root_backbone_frac=self.relaxed_max_root_backbone_frac,
            max_tm90_cluster_frac=self.relaxed_max_tm90_cluster_frac,
            max_structure_method_frac=self.relaxed_max_structure_method_frac,
            max_seq_method_frac=1.0,  # (f) relaxed only as a last resort
            min_distinct_structure_methods=self.min_distinct_structure_methods,
            min_levenshtein=self.min_levenshtein,
        )


@dataclass
class SheetResult:
    target: str
    ranked: list[dict]
    unranked: list[dict]
    deviations: list[dict] = field(default_factory=list)
    relaxation_steps: list[str] = field(default_factory=list)
    diagnostics: dict = field(default_factory=dict)

    @property
    def n_ranked(self) -> int:
        return len(self.ranked)


# --------------------------------------------------------------------------- #
# (b) bounded Levenshtein
# --------------------------------------------------------------------------- #


def levenshtein_within(a: str, b: str, limit: int) -> bool:
    """True when ``distance(a, b) <= limit``.  Banded DP with early exit."""
    if abs(len(a) - len(b)) > limit:
        return False
    if a == b:
        return True
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        cur = [i] + [0] * len(b)
        lo = max(1, i - limit)
        hi = min(len(b), i + limit)
        for j in range(1, len(b) + 1):
            if j < lo or j > hi:
                cur[j] = limit + 1
                continue
            cost = 0 if ca == b[j - 1] else 1
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
        if min(cur[lo : hi + 1] or [limit + 1]) > limit:
            return False
        prev = cur
    return prev[len(b)] <= limit


# --------------------------------------------------------------------------- #
# eligibility and ordering
# --------------------------------------------------------------------------- #


def _is_null(v: Any) -> bool:
    return v is None or v == "" or (isinstance(v, float) and v != v)


def _eligibility(row: Mapping[str, Any], schema: SheetSchema,
                 mask: InstrumentMask) -> str | None:
    """``None`` when rank-eligible, else the reason it goes to the unranked set."""
    # missing or zero scores -> unranked (the only thing the unranked set is for,
    # besides a missing required provenance column)
    for col in ("final_score", "rank_zscore"):
        if _is_null(row.get(col)):
            return f"missing {col}"
    if row.get("final_score") == 0 and row.get("rank_zscore") == 0:
        return "zero score"

    # the pose check must have RUN: a NaN pose is a missing measurement, not a flag
    if _is_null(row.get("pose_dockq")):
        return "pose check not run (pose_dockq null)"

    # every pre-registered gate must carry a recorded PASS
    for gate in ("novelty_verdict", "liability_verdict",
                 "monomer_foldability_verdict", "structural_plausibility_verdict"):
        v = row.get(gate)
        if v is None:
            return f"gate result missing: {gate}"
        if v != "PASS" and not row.get("relaxation_step"):
            return f"gate {gate} = {v}"

    # diversity gates need these columns non-null regardless of score
    for col in schema.diversity_required:
        if _is_null(row.get(col)):
            return f"missing {col}"

    return None


def _sort_key(row: Mapping[str, Any], mask: InstrumentMask) -> tuple:
    """``(n_seeds>=5 on all arms) DESC, pose_PASS DESC, rank_zscore DESC``."""
    full_seeds = all(int(row.get(f"n_seeds_{a}", 0) or 0) >= 5 for a in mask.arms)
    pose = bool(row.get("pose_PASS"))
    z = float(row.get("rank_zscore") or float("-inf"))
    return (not full_seeds, not pose, -z)


# --------------------------------------------------------------------------- #
# greedy capped selection
# --------------------------------------------------------------------------- #


def _select(
    ordered: Sequence[Mapping[str, Any]],
    caps: SelectionCaps,
    n_rows: int,
    allow_liability_flagged: bool,
) -> tuple[list[dict], dict[str, int]]:
    """Walk the ordered pool once, admitting rows that keep every cap satisfied."""
    root_cap = caps.cap(caps.max_root_backbone_frac, n_rows)
    clus_cap = caps.cap(caps.max_tm90_cluster_frac, n_rows)
    smeth_cap = caps.cap(caps.max_structure_method_frac, n_rows)
    qmeth_cap = caps.cap(caps.max_seq_method_frac, n_rows)

    chosen: list[dict] = []
    seqs: set[str] = set()
    root_n: dict[str, int] = {}
    clus_n: dict[str, int] = {}
    smeth_n: dict[str, int] = {}
    qmeth_n: dict[str, int] = {}
    blocked = {"dup": 0, "lev": 0, "root": 0, "cluster": 0,
               "structure_method": 0, "seq_method": 0, "liability": 0}

    for row in ordered:
        if len(chosen) >= n_rows:
            break
        seq = row["sequence"]
        # (a) exact-sequence duplicates
        if seq in seqs:
            blocked["dup"] += 1
            continue
        # (b) pairs within Levenshtein distance five
        if any(levenshtein_within(seq, c["sequence"], caps.min_levenshtein)
               for c in chosen):
            blocked["lev"] += 1
            continue
        if not allow_liability_flagged and row.get("liability_verdict") != "PASS":
            blocked["liability"] += 1
            continue
        # (c)
        root = row["root_backbone_id"]
        if root_n.get(root, 0) >= root_cap:
            blocked["root"] += 1
            continue
        # (d)
        clus = row["tm90_cluster_id"]
        if clus_n.get(clus, 0) >= clus_cap:
            blocked["cluster"] += 1
            continue
        # (e)
        sm = row["structure_method"]
        if smeth_n.get(sm, 0) >= smeth_cap:
            blocked["structure_method"] += 1
            continue
        # (f) - backfill from the next-best alternate seq_method
        qm = row.get("seq_method")
        if qm and qmeth_n.get(qm, 0) >= qmeth_cap:
            blocked["seq_method"] += 1
            continue

        chosen.append(dict(row))
        seqs.add(seq)
        root_n[root] = root_n.get(root, 0) + 1
        clus_n[clus] = clus_n.get(clus, 0) + 1
        smeth_n[sm] = smeth_n.get(sm, 0) + 1
        if qm:
            qmeth_n[qm] = qmeth_n.get(qm, 0) + 1

    return chosen, blocked


# --------------------------------------------------------------------------- #
# gate recompute (fail-loud)
# --------------------------------------------------------------------------- #


def recompute_row(
    row: Mapping[str, Any],
    mask: InstrumentMask,
    thresholds: GateThresholds,
) -> dict[str, Any]:
    """Recompute at write time what the row claims to carry.

    ``final_score`` is recomputed from the frozen mask's raw terms;
    ``pose_dockq`` from the per-arm ``sc_DockQ`` values; liability and LCP from
    the sequence.  A carried value that disagrees by more than 1e-4 halts the
    writer with the row id.
    """
    terms: list[float] = []
    for arm in mask.arms:
        terms.append(float(row.get(f"ipsae_{arm}", float("nan"))))
        terms.append(float(row.get(f"sc_DockQ_{arm}", float("nan"))))
        if mask.counter_screened:
            terms.append(float(row.get(f"selectivity_delta_{arm}", float("nan"))))
    live = [t for t in terms if not math.isnan(t)]
    final = sum(live) / len(live) if live else float("nan")

    dockqs = [
        float(row.get(f"sc_DockQ_{a}", float("nan")))
        for a in mask.arms
        if not _is_null(row.get(f"sc_DockQ_{a}"))
    ]
    pose = min(dockqs) if dockqs and not any(math.isnan(d) for d in dockqs) else float("nan")

    liab = liability_check(row["sequence"], thresholds)
    return {
        "final_score": final,
        "pose_dockq": pose,
        "pose_PASS": (pose >= mask.pose_dockq_threshold) if not math.isnan(pose) else None,
        "liability_verdict": liab.verdict,
        "lcp_score": lcp_score(row["sequence"]),
        "binder_len": len(row["sequence"].strip()),
    }


def _assert_recompute(row: Mapping[str, Any], recomputed: Mapping[str, Any]) -> None:
    for col, want in recomputed.items():
        have = row.get(col)
        if isinstance(want, float):
            if math.isnan(want) and (_is_null(have) or
                                     (isinstance(have, float) and math.isnan(have))):
                continue
            if _is_null(have) or abs(float(have) - want) > TOLERANCE:
                raise ValueError(
                    f"sheet writer HALT on row {row.get('design_id')!r}: "
                    f"carried {col}={have!r} disagrees with recomputed {want!r}"
                )
        else:
            if have != want:
                raise ValueError(
                    f"sheet writer HALT on row {row.get('design_id')!r}: "
                    f"carried {col}={have!r} disagrees with recomputed {want!r}"
                )


# --------------------------------------------------------------------------- #
# the writer
# --------------------------------------------------------------------------- #


def select_and_rank(
    rows: Iterable[Mapping[str, Any]],
    *,
    target: str,
    mask: InstrumentMask,
    schema: SheetSchema,
    vocab: MethodVocab | None = None,
    caps: SelectionCaps = SelectionCaps(),
    thresholds: GateThresholds = GateThresholds(),
    n_rows: int = ROWS_PER_TARGET,
    regeneration_tried: bool = True,
    recompute: Callable[[Mapping[str, Any]], Mapping[str, Any]] | None = None,
) -> SheetResult:
    """Produce one target's ranked sheet.

    ``regeneration_tried`` gates the relaxation ladder: the prompt requires
    upstream regeneration (more designs, more methods) to have been tried before
    any cap is loosened, and relaxation is never the first response to an
    under-diverse pool or to low scores.
    """
    vocab = vocab or default_method_vocab()
    pool = [dict(r) for r in rows]
    deviations: list[dict] = []
    unranked: list[dict] = []

    # --- canonicalise provenance keys, reject unknown tokens ---------------- #
    canonical: list[dict] = []
    for row in pool:
        try:
            row["structure_method"] = vocab.canonical_structure_method(
                row["structure_method"]
            )
            if row.get("seq_method"):
                row["seq_method"] = vocab.canonical_seq_method(row["seq_method"])
        except (KeyError, ValueError) as exc:
            unranked.append({**row, "unranked_reason": str(exc)})
            continue
        canonical.append(row)

    # --- eligibility --------------------------------------------------------- #
    # Routing comes BEFORE the recompute: a *missing* measurement (null score,
    # un-run pose check, missing provenance column) belongs in the unranked
    # section, whereas a *present but wrong* value is corruption and halts the
    # writer.  Recomputing first would turn every missing value into a halt.
    eligible: list[dict] = []
    for row in canonical:
        reason = _eligibility(row, schema, mask)
        if reason:
            unranked.append({**row, "unranked_reason": reason})
        else:
            eligible.append(row)

    # --- counter-screen fail-loud ------------------------------------------- #
    # "A null off-target or selectivity value on any ranked row of a
    # counter-screened target is a sheet-writer fail-loud, not a documented
    # gap" - so it is an error on a rank-eligible row, not a quiet demotion.
    if mask.counter_screened:
        for row in eligible:
            for arm in mask.arms:
                if _is_null(row.get(f"ipsae_offtarget_{arm}")) or _is_null(
                    row.get(f"selectivity_delta_{arm}")
                ):
                    raise ValueError(
                        f"{target} is counter-screened: rank-eligible row "
                        f"{row.get('design_id')!r} has no off-target score on {arm}; "
                        "all 30 ranked rows must have run the counter-screen "
                        "before the sheet is canonical"
                    )

    # --- fail-loud recompute of every gate ---------------------------------- #
    recompute = recompute or (lambda r: recompute_row(r, mask, thresholds))
    for row in eligible:
        _assert_recompute(row, recompute(row))

    ordered = sorted(eligible, key=lambda r: _sort_key(r, mask))

    # --- selection, then the relaxation ladder ------------------------------ #
    relaxation_steps: list[str] = []
    chosen, blocked = _select(ordered, caps, n_rows, allow_liability_flagged=False)

    if len(chosen) < n_rows and regeneration_tried:
        # (i) diversity caps
        relaxed = caps.relaxed()
        cand, blocked = _select(ordered, relaxed, n_rows, allow_liability_flagged=False)
        if len(cand) > len(chosen):
            relaxation_steps.append("DIVERSITY_CAPS")
            for r in cand:
                if r["design_id"] not in {c["design_id"] for c in chosen}:
                    r["relaxation_step"] = "DIVERSITY_CAPS"
            chosen = cand
            deviations.append({
                "kind": "relaxation", "target": target,
                "what": "diversity caps relaxed to per-root 25% / per-method 50%",
                "action": "recorded on each admitted row",
            })

    if len(chosen) < n_rows and regeneration_tried:
        # (ii) liability flags
        pool2 = ordered + [
            r for r in unranked
            if r.get("unranked_reason", "").startswith("gate liability_verdict")
        ]
        pool2.sort(key=lambda r: _sort_key(r, mask))
        cand, blocked = _select(pool2, caps.relaxed(), n_rows,
                                allow_liability_flagged=True)
        if len(cand) > len(chosen):
            relaxation_steps.append("LIABILITY_FLAGS")
            prev_ids = {c["design_id"] for c in chosen}
            for r in cand:
                if r["design_id"] not in prev_ids:
                    r["relaxation_step"] = "LIABILITY_FLAGS"
            chosen = cand
            deviations.append({
                "kind": "relaxation", "target": target,
                "what": "liability flags relaxed",
                "action": "recorded on each admitted row",
            })

    # Novelty is relaxed only as a last resort, and only by the caller passing a
    # pool that already contains novelty-flagged rows with the absolute floors
    # still clear (no identical rows, no target-mimic, no natural-sequence copy).
    if len(chosen) < n_rows:
        deviations.append({
            "kind": "short_sheet", "target": target,
            "what": f"only {len(chosen)} of {n_rows} rank-eligible rows after the "
                    f"relaxation ladder (blocked: {blocked})",
            "action": "shipping the actual N; padding with duplicates is forbidden",
        })

    # --- (e) at least three distinct structure methods ---------------------- #
    n_methods = len({r["structure_method"] for r in chosen})
    if chosen and n_methods < caps.min_distinct_structure_methods:
        deviations.append({
            "kind": "diversity_floor", "target": target,
            "what": f"only {n_methods} distinct structure_methods in the "
                    f"selected set (floor {caps.min_distinct_structure_methods})",
            "action": "fix upstream: more designs, more methods",
        })

    # --- rank assignment ----------------------------------------------------- #
    chosen.sort(key=lambda r: _sort_key(r, mask))
    for i, row in enumerate(chosen, start=1):
        row["rank"] = i
        row["target"] = target
        row["n_seeds"] = min(
            int(row.get(f"n_seeds_{a}", 0) or 0) for a in mask.arms
        )
        row["score_instrument"] = mask.describe()
        row.setdefault("relaxation_step", "")
        if row["n_seeds"] < 5:
            deviations.append({
                "kind": "seed_shortfall", "target": target,
                "what": f"{row['design_id']} ships at n_seeds={row['n_seeds']}",
                "action": "seed top-up wave attempted and failed; disclosed on the row",
            })

    # --- closing checks ------------------------------------------------------ #
    for row in chosen:
        schema.validate_row(
            {k: row.get(k) for k in schema.columns}, ranked=True
        )
        if _is_null(row.get("pose_dockq")):
            raise ValueError(
                f"{row['design_id']}: pose_dockq is mandatory_nonnull and is "
                "never relaxed"
            )
    seqs = [r["sequence"] for r in chosen]
    if len(set(seqs)) != len(seqs):
        raise ValueError(f"{target}: duplicate sequences in the shipped set")

    diagnostics = {
        "n_pool": len(pool),
        "n_eligible": len(eligible),
        "n_selected": len(chosen),
        "blocked_by": blocked,
        "n_distinct_structure_methods": n_methods,
        **fold_diversity_report(chosen),
    }
    return SheetResult(
        target=target,
        ranked=chosen,
        unranked=unranked,
        deviations=deviations,
        relaxation_steps=relaxation_steps,
        diagnostics=diagnostics,
    )


def fold_diversity_report(rows: Sequence[Mapping[str, Any]]) -> dict:
    """Fold diversity is a *reported* target, not a ranking gate.

    At least 10 % of each order sheet should be non-all-alpha (a design counts
    as not-all-alpha under DSSP if it has at least one beta strand of >= 3
    consecutive E/B residues, or its helical fraction is below 70 %).  If the
    ranked pool cannot supply it without displacing materially higher-scoring
    designs, ship fewer and state the count and the reason.
    """
    n = len(rows)
    n_non_alpha = sum(1 for r in rows if r.get("fold_class") == "not_all_alpha")
    frac = n_non_alpha / n if n else 0.0
    return {
        "n_non_all_alpha": n_non_alpha,
        "fraction_non_all_alpha": frac,
        "fold_diversity_target_met": frac >= 0.10,
    }
