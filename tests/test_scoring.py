"""The three-arm scoring instrument."""

import math

import pytest

from binder_campaign.scoring import (
    DEFAULT_ARMS,
    DesignScore,
    InstrumentMask,
    SeedRecord,
    aggregate_arm,
    score_pool,
    term_weights,
)


def _recs(vals):
    """(seed, ipsae_ab, ipsae_ba, dockq) tuples -> SeedRecords."""
    return [SeedRecord(seed=s, ipsae_ab=a, ipsae_ba=b, dockq=d) for s, a, b, d in vals]


def test_ipsae_min_is_the_minimum_over_both_alignment_directions():
    r = SeedRecord(seed=1, ipsae_ab=0.61, ipsae_ba=0.48, dockq=0.5)
    assert r.ipsae_min == 0.48


def test_arm_score_is_max_over_seeds_and_dockq_comes_from_that_same_seed():
    agg = aggregate_arm("ef2full", _recs([
        (11, 0.30, 0.30, 0.90),   # best DockQ ...
        (12, 0.55, 0.52, 0.40),   # ... but this seed has the best ipSAE
        (13, 0.41, 0.40, 0.60),
    ]))
    assert agg.ipsae_min == pytest.approx(0.52)
    assert agg.sc_dockq == pytest.approx(0.40)  # NOT 0.90
    assert agg.argmax_seed_ipsae == 12
    assert agg.argmax_seed_dockq == 11  # recorded so seed concordance is auditable


def test_seeds_must_be_distinct():
    with pytest.raises(ValueError, match="DISTINCT"):
        aggregate_arm("ef2fast", _recs([(5, 0.4, 0.4, 0.4), (5, 0.5, 0.5, 0.5)]))


def _design(ipsae, dockq, off=None, n_seeds=5):
    on = {
        arm: aggregate_arm(arm, _recs([
            (100 * i + k, ipsae[arm], ipsae[arm], dockq[arm]) for k in range(n_seeds)
        ]))
        for i, arm in enumerate(DEFAULT_ARMS)
    }
    off_book = {}
    if off:
        off_book = {
            arm: aggregate_arm(arm, _recs([
                (100 * i + k, off[arm], off[arm], 0.3) for k in range(n_seeds)
            ]))
            for i, arm in enumerate(DEFAULT_ARMS)
        }
    return DesignScore("d1", "PD-L1", on, off_book)


def test_final_score_is_the_raw_mean_of_the_six_terms():
    d = _design(
        {"ef2full": 0.60, "ef2fast": 0.50, "ptxv2": 0.40},
        {"ef2full": 0.80, "ef2fast": 0.70, "ptxv2": 0.60},
    )
    mask = InstrumentMask()
    assert mask.n_terms == 6
    assert d.final_score(mask) == pytest.approx((0.6 + 0.5 + 0.4 + 0.8 + 0.7 + 0.6) / 6)


def test_pose_dockq_is_the_min_over_arms_and_pass_uses_the_frozen_threshold():
    d = _design(
        {"ef2full": 0.6, "ef2fast": 0.6, "ptxv2": 0.6},
        {"ef2full": 0.80, "ef2fast": 0.22, "ptxv2": 0.60},
    )
    mask = InstrumentMask()
    assert d.pose_dockq(mask) == pytest.approx(0.22)
    assert d.pose_pass(mask) is False  # below the 0.23 default
    assert d.pose_pass(InstrumentMask(pose_dockq_threshold=0.20)) is True


def test_reduced_mask_means_over_realized_terms_only_and_names_them():
    mask = InstrumentMask(name="default_3arm", arms=("ef2full", "ptxv2"))
    assert mask.reduced
    assert mask.describe() == "ef2full_ptxv2"
    d = _design(
        {"ef2full": 0.60, "ef2fast": 0.50, "ptxv2": 0.40},
        {"ef2full": 0.80, "ef2fast": 0.70, "ptxv2": 0.60},
    )
    assert d.final_score(mask) == pytest.approx((0.6 + 0.4 + 0.8 + 0.6) / 4)


def test_counter_screened_mask_has_nine_terms_and_ipsae_weight_on_selectivity():
    mask = InstrumentMask(counter_screened=True)
    assert mask.n_terms == 9
    w = term_weights(mask)
    assert w["ipsae_ef2full"] == 4.0
    assert w["sc_DockQ_ef2full"] == 1.0
    assert w["selectivity_delta_ef2full"] == 4.0


def test_selectivity_delta_is_on_minus_off():
    d = _design(
        {"ef2full": 0.60, "ef2fast": 0.50, "ptxv2": 0.40},
        {"ef2full": 0.80, "ef2fast": 0.70, "ptxv2": 0.60},
        off={"ef2full": 0.25, "ef2fast": 0.20, "ptxv2": 0.15},
    )
    assert d.selectivity_delta("ef2full") == pytest.approx(0.35)


def test_selectivity_delta_refuses_an_instrument_mismatched_subtraction():
    on = {"ef2full": aggregate_arm("ef2full", _recs([(i, 0.6, 0.6, 0.5) for i in range(5)]))}
    off = {"ef2full": aggregate_arm("ef2full", _recs([(i, 0.3, 0.3, 0.4) for i in range(3)]))}
    d = DesignScore("d", "GDF-8", on, off)
    with pytest.raises(ValueError, match="matched seed counts"):
        d.selectivity_delta("ef2full")


def test_rank_zscore_is_the_4_to_1_weighted_z_average_over_the_pool():
    hi = _design({a: 0.9 for a in DEFAULT_ARMS}, {a: 0.1 for a in DEFAULT_ARMS})
    hi.design_id = "hi"
    lo = _design({a: 0.1 for a in DEFAULT_ARMS}, {a: 0.9 for a in DEFAULT_ARMS})
    lo.design_id = "lo"
    mask = InstrumentMask()
    rows = {r["design_id"]: r for r in score_pool([hi, lo], mask)}

    # both designs have identical final_score (mean of the same six numbers) ...
    assert rows["hi"]["final_score"] == pytest.approx(rows["lo"]["final_score"])
    # ... but ipSAE is weighted 4:1 over DockQ, so ranking separates them
    assert rows["hi"]["rank_zscore"] > rows["lo"]["rank_zscore"]
    # two-point pool: z = +/-1, so the average is +/-(4-1)/5 = +/-0.6
    assert rows["hi"]["rank_zscore"] == pytest.approx(0.6)
    assert rows["lo"]["rank_zscore"] == pytest.approx(-0.6)


def test_rank_zscore_is_transductive_and_raw_terms_are_carried():
    designs = [
        _design({a: v for a in DEFAULT_ARMS}, {a: 0.5 for a in DEFAULT_ARMS})
        for v in (0.2, 0.4, 0.6)
    ]
    for i, d in enumerate(designs):
        d.design_id = f"d{i}"
    mask = InstrumentMask()
    rows = score_pool(designs, mask)
    # raw values survive so any other batch can re-standardise
    assert [r["ipsae_ef2full"] for r in rows] == pytest.approx([0.2, 0.4, 0.6])
    # z-scores are pool-relative: dropping the top design shifts them
    rows2 = score_pool(designs[:2], mask)
    assert rows2[0]["rank_zscore"] != rows[0]["rank_zscore"]


def test_zero_variance_term_contributes_zero_not_nan():
    designs = [
        _design({a: 0.5 for a in DEFAULT_ARMS}, {a: 0.5 for a in DEFAULT_ARMS})
        for _ in range(3)
    ]
    for i, d in enumerate(designs):
        d.design_id = f"d{i}"
    rows = score_pool(designs, InstrumentMask())
    assert all(r["rank_zscore"] == pytest.approx(0.0) for r in rows)
    assert not any(math.isnan(r["rank_zscore"]) for r in rows)


def test_mask_rejects_duplicate_arms():
    with pytest.raises(ValueError, match="duplicate arms"):
        InstrumentMask(arms=("ef2full", "ef2full"))
