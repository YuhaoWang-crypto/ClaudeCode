"""The design driver: measured singles in, ranked combinatorial pick-list out.

This is the part you reuse on a new enzyme. One call does what the paper's
notebook spread over four hand-edited task blocks:

    select AAindex descriptors -> fit PLSR -> score the whole combinatorial
    space -> compare against baselines -> emit a diversity-filtered pick-list

The pick-list, not the ranking, is the deliverable. Two things it does that
ranking by predicted fitness alone does not:

* **Diversity filter.** PLS top-N lists are typically near-duplicates of one
  mutation set. Capping the Jaccard overlap between picks spends a fixed assay
  budget on distinguishable hypotheses instead of ten spellings of one.
* **Honest headline.** Every report carries the held-out or cross-validated R2
  *and* the best baseline's, because a pick-list from a model that cannot beat
  log-additive is a pick-list you could have written by hand.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .baselines import ALL_BASELINES, LogAdditiveModel
from .encode import MutationEncoder
from .model import IndexSelection, fit_predict, r2, screen_indices
from .variants import enumerate_combinations, mutated_positions, parse_variant

__all__ = ["DesignReport", "design_round", "improved_sites", "evaluate_against_baselines"]


@dataclass
class DesignReport:
    """Everything one design round produced."""

    selection: IndexSelection
    ranking: pd.DataFrame
    picks: pd.DataFrame
    baselines: pd.DataFrame
    n_train: int
    n_space: int
    parent_fitness: float
    notes: list[str] = field(default_factory=list)

    @property
    def headline(self) -> str:
        best_base = self.baselines["cvR2"].max() if len(self.baselines) else float("nan")
        return (
            f"FFT-PLSR [{self.selection.label}] n_components={self.selection.n_components}: "
            f"cvR2={self.selection.cv_r2:.3f} on {self.n_train} variants; "
            f"best baseline cvR2={best_base:.3f}; "
            f"{self.n_space} variants scored, {len(self.picks)} picked"
        )

    def summary(self) -> str:
        lines = [self.headline, ""]
        lines.append("Baselines (leave-one-out on the same training set):")
        lines.append(self.baselines.to_string(index=False))
        lines.append("")
        lines.append(f"Pick-list (top {len(self.picks)}, diversity-filtered):")
        lines.append(self.picks.to_string(index=False))
        if self.notes:
            lines.append("")
            lines.extend(f"note: {n}" for n in self.notes)
        return "\n".join(lines)


def improved_sites(
    measured: dict[str, float],
    parent_fitness: float = 1.0,
    threshold: float = 1.10,
    per_position: int | None = None,
) -> list[str]:
    """Pick which single mutations are worth recombining.

    Keeps every single whose fitness is at least `threshold` x parent (the paper
    used +10%). Prune *positions* this way, not substitutions within a position.

    `per_position` additionally keeps only the N best substitutions at each
    position. **It defaults to off, and you should usually leave it off.** A
    single mutation's solo effect is a poor guide to its value in combination,
    which is the whole premise of fitting an epistasis-aware model in the first
    place. Measured on this dataset: the best variant in the released panel,
    Com2-IFRS = N7Y/H63L/K67N/V74W at 2.75x parent, is built from singles ranking
    3rd of 3, 1st of 7, 4th of 5 and 2nd of 2 at their positions -- three of the
    four are in the bottom half. `per_position=3` removes K67N and makes that
    variant unreachable before the model is even fitted.

    Keeping everything above `threshold` here reproduces the 11,520-variant space
    the paper enumerated (11,492 combinations of order >= 2, plus 27 singles and
    the parent).
    """
    singles: dict[int, list[tuple[float, str]]] = {}
    for label, value in measured.items():
        subs = parse_variant(label)
        if len(subs) != 1:
            continue
        if value >= threshold * parent_fitness:
            singles.setdefault(subs[0].pos, []).append((float(value), label))

    out = []
    for pos in sorted(singles):
        ranked = sorted(singles[pos], reverse=True)
        out.extend(label for _, label in ranked[: per_position or len(ranked)])
    return out


def _jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 1.0


def diversify(ranking: pd.DataFrame, pick: int, max_jaccard: float = 0.6) -> pd.DataFrame:
    """Greedily take the top-ranked variants whose mutation sets are not near-duplicates."""
    chosen: list[int] = []
    chosen_sets: list[set[str]] = []
    for row in ranking.itertuples():
        muts = {str(s) for s in parse_variant(row.variant)}
        if any(_jaccard(muts, prev) > max_jaccard for prev in chosen_sets):
            continue
        chosen.append(row.Index)
        chosen_sets.append(muts)
        if len(chosen) == pick:
            break
    return ranking.loc[chosen].reset_index(drop=True)


def evaluate_against_baselines(labels, y, cv=None) -> pd.DataFrame:
    """Leave-one-out (or `cv`-fold) R2/MSE for every baseline, on the given training set."""
    from sklearn.model_selection import KFold, LeaveOneOut

    labels = list(labels)
    y = np.asarray(y, dtype=float)
    splitter = LeaveOneOut() if cv is None else KFold(n_splits=cv, shuffle=True, random_state=0)
    rows = []
    for cls in ALL_BASELINES:
        oof = np.full(len(y), np.nan)
        for train, test in splitter.split(np.zeros((len(y), 1))):
            try:
                mdl = cls().fit([labels[i] for i in train], y[train])
                oof[test] = mdl.predict([labels[i] for i in test])
            except Exception:
                pass
        ok = ~np.isnan(oof)
        rows.append(
            {
                "model": cls.name,
                "cvR2": r2(y[ok], oof[ok]) if ok.sum() > 1 else float("nan"),
                "cvMSE": float(((y[ok] - oof[ok]) ** 2).mean()) if ok.any() else float("nan"),
                "n_scored": int(ok.sum()),
            }
        )
    return pd.DataFrame(rows).sort_values("cvR2", ascending=False).reset_index(drop=True)


def design_round(
    parent: str,
    measured: dict[str, float],
    sites=None,
    *,
    space=None,
    n_rounds: int = 1,
    cv=None,
    pick: int = 8,
    max_order: int | None = None,
    max_jaccard: float = 0.6,
    candidates=None,
    n_jobs: int = -1,
    verbose: bool = True,
) -> DesignReport:
    """Run one FFT-PLSR design round.

    Parameters
    ----------
    parent:
        Background sequence all labels in `measured` are expressed against.
    measured:
        ``{variant label: relative fitness}``. The parent itself (label without
        substitutions, e.g. "WT") should be included at its reference value.
    sites:
        Single mutations to recombine. Defaults to every single mutation present
        in `measured`; pass `improved_sites(...)` to keep only the beneficial ones.
    space:
        Explicit variant labels to score, overriding combinatorial enumeration --
        use this to score a pre-built list or a saturation scan.
    n_rounds:
        How many AAindex entries to select greedily (paper used 1 or 3).
    cv:
        Fold count; ``None`` means leave-one-out, as the paper used for small sets.
    pick / max_jaccard:
        Size of the pick-list and the diversity cap applied when building it.
    max_order:
        Cap on how many mutations may be combined. ``None`` means no cap.
    """
    labels = list(measured)
    y = np.asarray([measured[k] for k in labels], dtype=float)
    notes: list[str] = []

    parent_fitness = next((v for k, v in measured.items() if not parse_variant(k)), 1.0)

    if sites is None:
        sites = [k for k in labels if len(parse_variant(k)) == 1]
    sites = list(sites)
    if not sites:
        raise ValueError("no single mutations to recombine -- pass `sites` explicitly")

    if space is None:
        space = enumerate_combinations(parent, sites, min_order=2, max_order=max_order)
    space = list(space)
    # Score the singles too, so the ranking shows what recombination adds over them.
    query = list(dict.fromkeys([*sites, *space]))

    positions = sorted(set(mutated_positions(labels)) | set(mutated_positions(query)))
    encoder = MutationEncoder(parent, positions=positions)

    if verbose:
        print(
            f"FFT-PLSR design round: {len(labels)} measured variants, {len(sites)} sites, "
            f"{len(space)} combinations, {len(positions)} mutable positions"
        )

    selection = screen_indices(
        lambda codes: encoder.encode(labels, codes),
        y,
        n_rounds=n_rounds,
        cv=cv,
        candidates=candidates,
        n_jobs=n_jobs,
        verbose=verbose,
    )

    X_train = encoder.encode(labels, selection.indices)
    X_query = encoder.encode(query, selection.indices)
    pred = fit_predict(X_train, y, X_query, selection.n_components)

    additive = LogAdditiveModel().fit(labels, y).predict(query)
    ranking = (
        pd.DataFrame(
            {
                "variant": query,
                "predicted": pred,
                "log_additive": additive,
                "order": [max(len(parse_variant(q)), 1) for q in query],
                "measured": [measured.get(q, np.nan) for q in query],
            }
        )
        .sort_values("predicted", ascending=False)
        .reset_index(drop=True)
    )
    ranking["epistasis"] = ranking["predicted"] - ranking["log_additive"]

    untested = ranking[ranking["measured"].isna()].reset_index(drop=True)
    if len(untested) < len(ranking):
        notes.append(
            f"{len(ranking) - len(untested)} of {len(ranking)} scored variants already have "
            "measurements and were excluded from the pick-list"
        )
    picks = diversify(untested, pick=pick, max_jaccard=max_jaccard)

    base = evaluate_against_baselines(labels, y, cv=cv)
    if selection.cv_r2 <= base["cvR2"].max():
        notes.append(
            "FFT-PLSR did not beat the best baseline on cross-validation -- treat the "
            "ranking as a tie-breaker, not evidence, and prefer the cheaper model"
        )

    return DesignReport(
        selection=selection,
        ranking=ranking,
        picks=picks,
        baselines=base,
        n_train=len(labels),
        n_space=len(query),
        parent_fitness=parent_fitness,
        notes=notes,
    )
