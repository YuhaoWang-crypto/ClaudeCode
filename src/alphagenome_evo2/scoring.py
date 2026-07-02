"""Cell-type-specificity fitness.

Turns AlphaGenome read-outs into a single scalar the loop maximises. The design
goal from the spec is an element that is **active in the positive context** and
**silent in every negative context**, without cryptic splicing / off-target
activity. Concretely::

    fitness = w_expr   * expr(positive)
            + w_acc    * accessibility(positive)
            + w_tf     * tf_binding(positive)
            - w_off    * max_over_negatives( expr(negative) )
            - w_splice * max_over_contexts( splice_anomaly )

Using the *max* over negative contexts (rather than the mean) enforces silence in
the single worst offending line, which is what specificity actually requires.
"""

from __future__ import annotations

from typing import Dict, Tuple

from .config import FitnessWeights, PipelineConfig
from .types import Candidate


def score_candidate(
    candidate: Candidate,
    config: PipelineConfig,
) -> Tuple[float, Dict[str, float]]:
    """Return ``(fitness, term_breakdown)`` for a candidate with predictions."""
    w: FitnessWeights = config.weights
    pos = candidate.prediction_for(config.positive_context)
    if pos is None:
        raise ValueError(
            f"candidate {candidate.id} lacks a prediction for the positive "
            f"context {config.positive_context.name!r}; score it with the oracle first"
        )

    on_expr = w.on_target_expression * pos.expression
    on_acc = w.on_target_accessibility * pos.accessibility
    on_tf = w.tf_binding_bonus * pos.tf_binding

    # Worst-case leakage into any negative context.
    off_expr = 0.0
    for ctx in config.negative_contexts:
        neg = candidate.prediction_for(ctx)
        if neg is not None:
            off_expr = max(off_expr, neg.expression)
    off_term = w.off_target_expression * off_expr

    # Worst-case splicing anomaly across all contexts.
    splice = pos.splice_anomaly
    for ctx in config.negative_contexts:
        neg = candidate.prediction_for(ctx)
        if neg is not None:
            splice = max(splice, neg.splice_anomaly)
    splice_term = w.splice_penalty * splice

    fitness = on_expr + on_acc + on_tf - off_term - splice_term
    terms = {
        "on_target_expression": on_expr,
        "on_target_accessibility": on_acc,
        "on_target_tf_binding": on_tf,
        "off_target_penalty": -off_term,
        "splice_penalty": -splice_term,
        "specificity_gap": pos.expression - off_expr,
    }
    return fitness, terms


def apply_scores(candidates, config: PipelineConfig):
    """Return candidates with ``fitness``/``fitness_terms`` populated."""
    out = []
    for c in candidates:
        fitness, terms = score_candidate(c, config)
        out.append(c.scored(fitness, terms))
    return out
