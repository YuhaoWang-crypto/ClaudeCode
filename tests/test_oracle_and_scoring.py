from alphagenome_evo2 import (
    CellContext,
    MockOracle,
    PipelineConfig,
    Sequence,
    score_candidate,
)
from alphagenome_evo2.oracle.mock import _motifs_for_context


def _config():
    return PipelineConfig(
        positive_context=CellContext("A", "ONT:A"),
        negative_contexts=[CellContext("B", "ONT:B")],
    )


def test_mock_oracle_is_deterministic():
    oracle = MockOracle()
    ctx = [CellContext("A", "ONT:A")]
    s = Sequence("s", "ACGTACGTAC" * 20)
    p1 = oracle.predict(s, ctx)["A"]
    p2 = oracle.predict(s, ctx)["A"]
    assert p1.expression == p2.expression
    assert p1.accessibility == p2.accessibility


def test_positive_motif_rich_sequence_scores_higher_than_random():
    config = _config()
    oracle = MockOracle()
    pos_motifs = _motifs_for_context(config.positive_context)

    # Sequence packed with positive-context motifs.
    rich = Sequence("rich", ("".join(pos_motifs) + "AAA") * 15)
    # A low-signal poly-A sequence.
    poor = Sequence("poor", "A" * 200)

    rich_cand = oracle.predict_batch([rich], config.all_contexts())[0]
    poor_cand = oracle.predict_batch([poor], config.all_contexts())[0]

    rich_fit, _ = score_candidate(rich_cand, config)
    poor_fit, _ = score_candidate(poor_cand, config)
    assert rich_fit > poor_fit


def test_scoring_penalises_off_target_activity():
    config = _config()
    oracle = MockOracle()
    pos_motifs = _motifs_for_context(config.positive_context)
    neg_motifs = _motifs_for_context(config.negative_contexts[0])

    specific = Sequence("spec", ("".join(pos_motifs)) * 15)
    leaky = Sequence("leaky", ("".join(pos_motifs + neg_motifs)) * 10)

    spec_cand = oracle.predict_batch([specific], config.all_contexts())[0]
    leaky_cand = oracle.predict_batch([leaky], config.all_contexts())[0]

    spec_fit, spec_terms = score_candidate(spec_cand, config)
    leaky_fit, leaky_terms = score_candidate(leaky_cand, config)

    # The leaky sequence has a smaller specificity gap.
    assert spec_terms["specificity_gap"] >= leaky_terms["specificity_gap"]
