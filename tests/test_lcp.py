"""LCP restraint — Figure 1 of the prompt bundle."""

import math

import pytest

from binder_campaign.lcp import (
    DEFAULT_THRESHOLD_ENTROPY,
    DEFAULT_WINDOW,
    lcp_per_position_penalty,
    lcp_report,
    lcp_score,
    rank_by_lcp,
    window_entropies,
)

# Human ubiquitin (P0CG47), a real natural sequence: no 30-mer window falls
# below the 5th-percentile PDB entropy, so LCP leaves it alone.
NATIVE_LIKE = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)

# An idealised E/K/A/L-rich designed helical bundle -- exactly the low-complexity
# regime the restraint exists to push sequence design away from.
IDEALISED_BUNDLE = (
    "SEEELKKLAEEAKRLAEEIKRQGYSAEEVREAIRLARENGNEELAKRVLEEAKKLGVDPEEIARRA"
    "LEIYKKTGDEKLAEEIRRLAEEHGLTPEQVREAVRLAKEQ"
)


def test_figure1_constants():
    assert DEFAULT_WINDOW == 30
    assert DEFAULT_THRESHOLD_ENTROPY == pytest.approx(2.32)


def test_homopolymer_is_penalised_and_native_sequence_is_not():
    assert lcp_score("A" * 60) > 0
    assert lcp_score(NATIVE_LIKE) == 0.0


def test_restraint_orders_native_below_idealised_bundle_below_homopolymer():
    """The whole point: it pushes design away from low-complexity sequence."""
    assert lcp_score(NATIVE_LIKE) < lcp_score(IDEALISED_BUNDLE) < lcp_score("A" * 60)


def test_penalty_is_the_squared_perplexity_gap_on_sub_threshold_windows():
    """C3 = L/(L-w+1) * sum (e^Shat - e^Si)^2 over windows with Si < Shat."""
    poly = "A" * 30  # exactly one window, entropy 0, perplexity 1
    expected = (30 / 1) * (math.exp(DEFAULT_THRESHOLD_ENTROPY) - 1.0) ** 2
    assert lcp_score(poly) == pytest.approx(expected)


def test_indicator_switches_the_penalty_off_above_threshold():
    """A window at or above the threshold entropy contributes nothing."""
    # 30 residues, 20 distinct amino acids -> entropy well above 2.32 nats
    seq = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    ents = window_entropies(seq)
    assert len(ents) == 1
    assert ents[0] > DEFAULT_THRESHOLD_ENTROPY
    assert lcp_score(seq) == 0.0


def test_window_count_is_L_minus_w_plus_1():
    seq = NATIVE_LIKE[:80]
    assert len(window_entropies(seq)) == len(seq) - DEFAULT_WINDOW + 1


def test_per_position_penalty_sums_to_the_total():
    seq = "A" * 20 + NATIVE_LIKE[:60]
    per_pos = lcp_per_position_penalty(seq)
    assert len(per_pos) == len(seq)
    assert sum(per_pos) == pytest.approx(lcp_score(seq), rel=1e-9)


def test_per_position_penalty_localises_to_the_low_complexity_stretch():
    seq = "K" * 30 + NATIVE_LIKE
    per_pos = lcp_per_position_penalty(seq)
    assert sum(per_pos[:30]) > sum(per_pos[60:])


def test_report_counts_penalised_windows():
    rep = lcp_report("A" * 60)
    assert rep.n_windows == 31
    assert rep.n_windows_penalised == 31
    assert rep.fraction_penalised == 1.0
    assert rep.min_window_entropy == 0.0


def test_short_sequence_still_scores():
    assert lcp_score("AAAA") > 0
    assert lcp_score("") == 0.0


def test_rank_by_lcp_puts_the_cleanest_sequence_first():
    ranked = rank_by_lcp([("A" * 60), NATIVE_LIKE, "K" * 30 + NATIVE_LIKE])
    assert ranked[0][0] == NATIVE_LIKE
    assert ranked[0][1] == 0.0
    assert ranked[-1][0] == "A" * 60
