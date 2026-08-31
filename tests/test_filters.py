"""The four pre-scoring gates."""

import pytest

from binder_campaign.filters import (
    GateThresholds,
    HomologyHit,
    assert_rejects_absent,
    liability_check,
    monomer_foldability_check,
    novelty_check,
    run_prescoring_gates,
    smith_waterman_identity,
    structural_plausibility_check,
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)
CLEAN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGA"
)


# --- liability --------------------------------------------------------------- #

def test_odd_cysteine_count_is_an_unpaired_cysteine():
    assert liability_check("ACDEFGHIK").verdict == "REJECT"
    assert "unpaired cysteine" in liability_check("ACDEFGHIK").reasons[0]
    assert liability_check("ACDEFGHIKC").passed


def test_homopolymer_runs_longer_than_four_are_flagged():
    assert liability_check("DEKR" + "K" * 5 + "DEKR").verdict == "REJECT"
    assert liability_check("DEKR" + "K" * 4 + "DEKR").passed


def test_surface_hydrophobic_patch_is_flagged():
    res = liability_check("DEKRDEKR" + "LIVFLIVFL" + "DEKRDEKR")
    assert res.verdict == "REJECT"
    assert any("hydrophobic patch" in r for r in res.reasons)


def test_liability_metrics_are_recorded_even_on_a_pass():
    res = liability_check(UBIQUITIN)
    assert res.passed
    assert res.metrics["n_cys"] == 0
    assert res.metrics["longest_run_len"] >= 1


# --- novelty ----------------------------------------------------------------- #

def test_rule_one_rejects_high_identity_over_high_coverage():
    hit = HomologyHit("UniRef90_X", "uniref90", identity=0.75, coverage=0.80)
    assert novelty_check("SEQ", [hit]).verdict == "REJECT"
    # 60% identity but only 40% coverage does not trip rule 1
    near = HomologyHit("UniRef90_Y", "uniref90", identity=0.75, coverage=0.40)
    assert novelty_check("SEQ", [near]).passed


def test_rule_three_catches_a_target_mimic_protomer_by_tm_score():
    res = novelty_check("SEQ", [], structural_tm_scores={"target_A": 0.62})
    assert res.verdict == "REJECT"
    assert "target-mimic" in res.reasons[0]
    assert novelty_check("SEQ", [], structural_tm_scores={"target_A": 0.49}).passed


def test_ubiquitin_with_terminal_extensions_is_caught_by_local_alignment():
    """The prompt calls this out by name: detect by local alignment, not exact match."""
    disguised = "GSHMGSGS" + UBIQUITIN + "GGSGGSLEHHHHHH"
    assert disguised != UBIQUITIN  # an exact-match check would let it through
    ident, aligned = smith_waterman_identity(disguised, UBIQUITIN)
    assert aligned >= 40 and ident >= 0.30
    res = novelty_check(disguised, [], control_sequences={"P0CG47": UBIQUITIN})
    assert res.verdict == "REJECT"


def test_an_unrelated_sequence_passes_the_local_alignment_rule():
    unrelated = "DEKRDEKRQNSTQNSTDEKRDEKRQNSTQNSTDEKRDEKRQNSTQNSTDEKRDEKRQNST"
    res = novelty_check(unrelated, [], control_sequences={"P0CG47": UBIQUITIN})
    assert res.passed


def test_smith_waterman_identity_is_one_for_a_self_alignment():
    ident, aligned = smith_waterman_identity(UBIQUITIN, UBIQUITIN)
    assert ident == pytest.approx(1.0)
    assert aligned == len(UBIQUITIN)


# --- foldability and plausibility --------------------------------------------- #

def test_plddt_is_normalised_so_the_frozen_threshold_means_one_thing():
    """ESMFold2 emits 0-1; the frozen threshold is on the 0-100 scale."""
    assert monomer_foldability_check(0.85).passed          # 85 on 0-100
    assert monomer_foldability_check(85.0).passed
    assert not monomer_foldability_check(0.65).passed
    assert not monomer_foldability_check(65.0).passed


def test_binder_length_window_needs_a_recorded_rationale_outside_50_to_120():
    ok = structural_plausibility_check(0, 0.8, 0.05, binder_len=80)
    assert ok.passed
    no_rationale = structural_plausibility_check(0, 0.8, 0.05, binder_len=40)
    assert not no_rationale.passed
    with_rationale = structural_plausibility_check(
        0, 0.8, 0.05, binder_len=40, length_rationale="flat epitope, short helix hairpin"
    )
    assert with_rationale.passed
    assert not structural_plausibility_check(0, 0.8, 0.05, binder_len=200).passed


def test_clashes_and_poor_packing_are_rejected():
    assert not structural_plausibility_check(3, 0.8, 0.05, 80).passed
    assert not structural_plausibility_check(0, 0.2, 0.05, 80).passed


# --- orchestration ------------------------------------------------------------ #

def test_run_prescoring_gates_records_every_gate_and_the_frozen_thresholds():
    rec = run_prescoring_gates("d1", CLEAN, mean_plddt=88.0)
    assert rec["verdict"] == "PASS"
    assert set(rec["gates"]) == {
        "novelty", "liability", "monomer_foldability", "structural_plausibility"
    }
    assert rec["thresholds"]["thr_min_monomer_plddt"] == 70.0


def test_a_failing_gate_names_itself_in_rejected_by():
    rec = run_prescoring_gates("d2", "DEKR" + "K" * 8 + "DEKR", mean_plddt=90.0)
    assert rec["verdict"] == "REJECT"
    assert "liability" in rec["rejected_by"]


def test_a_gate_counts_as_run_only_when_its_rejects_are_absent_downstream():
    assert_rejects_absent(["d1", "d2"], ["d3", "d4"])  # fine
    with pytest.raises(AssertionError, match="reached a downstream scoring pool"):
        assert_rejects_absent(["d1", "d2"], ["d2", "d3"])


def test_thresholds_are_serialisable_onto_every_design_row():
    row = GateThresholds().as_row()
    assert row["thr_max_global_identity"] == 0.60
    assert row["thr_max_tm_score_to_target"] == 0.50
