"""Selection caps, rank assignment, relaxation ladder, write-time gate recompute."""

import random

import pytest

from binder_campaign.filters import GateThresholds, liability_check
from binder_campaign.lcp import lcp_score
from binder_campaign.schema import default_method_vocab, default_sheet_schema
from binder_campaign.scoring import DEFAULT_ARMS, InstrumentMask
from binder_campaign.sheet_writer import (
    SelectionCaps,
    fold_diversity_report,
    levenshtein_within,
    recompute_row,
    select_and_rank,
)

SCHEMA = default_sheet_schema()
VOCAB = default_method_vocab()
MASK = InstrumentMask()
AA = "ACDEFGHIKLMNPQRSTVWY"
METHODS = ("rfdiffusion", "boltzgen", "genie3")


def _seq(i: int, n: int = 80) -> str:
    """Distinct, well-separated, liability-clean sequences.

    Seeded per-index rather than a cyclic pattern (rotations of a repeating
    motif are only a few edits apart and would trip cap (b) on every pair), and
    resampled until the real liability gate passes, so these fixtures exercise
    the *selection* rules rather than the liability filter.
    """
    rng = random.Random(1000 + i)
    for _ in range(200):
        seq = "".join(rng.choice(AA) for _ in range(n))
        if liability_check(seq).passed:
            return seq
    raise RuntimeError("could not sample a liability-clean fixture sequence")


def make_row(i, *, ipsae=0.5, dockq=0.5, root=None, cluster=None,
             method="rfdiffusion", seq_method="solublempnn", n_seeds=5,
             sequence=None, fold_class="all_alpha", liability="PASS"):
    seq = sequence if sequence is not None else _seq(i)
    row = {
        "design_id": f"d{i:04d}",
        "target": "PD-L1",
        "sequence": seq,
        "binder_len": len(seq),
        "rank": None,
        "rank_zscore": ipsae,          # ordering proxy
        "final_score": (ipsae + dockq) / 2,
        "score_instrument": "default_3arm",
        "pose_PASS": dockq >= MASK.pose_dockq_threshold,
        "pose_dockq": dockq,
        "structure_method": method,
        "seq_method": seq_method,
        "opt_round": 0,
        "root_backbone_id": root or f"bb{i:04d}",
        "parent_design_id": "",
        "n_seeds": n_seeds,
        "novelty_verdict_path": f"novelty/d{i}.json",
        "tm90_cluster_id": cluster or f"c{i:04d}",
        "fold_class": fold_class,
        "designed_structure_path": f"ordered/d{i}/designed.pdb",
        "monomer_plddt": 85.0,
        "lcp_score": lcp_score(seq),
        "esmc_ll": -1.0,
        "relaxation_step": "",
        "construct_status": "OK",
        "binder_binder_clashes_NN": 0,
        "novelty_verdict": "PASS",
        "liability_verdict": liability,
        "monomer_foldability_verdict": "PASS",
        "structural_plausibility_verdict": "PASS",
    }
    for arm in DEFAULT_ARMS:
        row[f"ipsae_{arm}"] = ipsae
        row[f"sc_DockQ_{arm}"] = dockq
        row[f"n_seeds_{arm}"] = n_seeds
        row[f"predicted_structure_path_{arm}"] = f"ordered/d{i}/predicted_{arm}.cif"
        row[f"selectivity_delta_{arm}"] = float("nan")
        row[f"ipsae_offtarget_{arm}"] = float("nan")
        row[f"ipsae_NN_{arm}"] = float("nan")
        row[f"sc_DockQ_NN_{arm}"] = float("nan")
    # keep the carried values consistent with the write-time recompute
    row.update(recompute_row(row, MASK, GateThresholds()))
    return row


def run(rows, **kw):
    kw.setdefault("target", "PD-L1")
    kw.setdefault("mask", MASK)
    kw.setdefault("schema", SCHEMA)
    kw.setdefault("vocab", VOCAB)
    return select_and_rank(rows, **kw)


# --- (b) bounded Levenshtein -------------------------------------------------- #

def test_levenshtein_within_is_exact_at_the_boundary():
    assert levenshtein_within("ABCDE", "ABCDE", 5)
    assert levenshtein_within("ABCDEFGHIJ", "ABCDEFGHIX", 5)
    assert not levenshtein_within(_seq(1), _seq(2), 5)
    assert not levenshtein_within("A" * 10, "A" * 20, 5)  # length gap alone


# --- caps --------------------------------------------------------------------- #

def test_exact_sequence_duplicates_are_rejected():
    dup = _seq(1)
    rows = [make_row(i, sequence=dup, ipsae=0.9 - 0.01 * i) for i in range(5)]
    rows += [make_row(100 + i, ipsae=0.5) for i in range(40)]
    res = run(rows)
    assert sum(1 for r in res.ranked if r["sequence"] == dup) == 1
    assert res.diagnostics["blocked_by"]["dup"] == 4


def test_pairs_within_levenshtein_five_are_rejected():
    base = _seq(1)
    near = base[:-3] + "AAA"  # distance 3
    rows = [make_row(1, sequence=base, ipsae=0.9),
            make_row(2, sequence=near, ipsae=0.85)]
    rows += [make_row(100 + i, ipsae=0.5) for i in range(40)]
    res = run(rows)
    ids = {r["design_id"] for r in res.ranked}
    assert "d0001" in ids and "d0002" not in ids
    assert res.diagnostics["blocked_by"]["lev"] >= 1


def test_root_backbone_capped_at_five_percent_rounded_up():
    """5 % of 30 rows = 1.5 -> 2 rows per root_backbone_id."""
    rows = [make_row(i, root="shared_bb", ipsae=0.9 - 0.001 * i) for i in range(20)]
    rows += [make_row(100 + i, ipsae=0.4) for i in range(40)]
    res = run(rows)
    assert sum(1 for r in res.ranked if r["root_backbone_id"] == "shared_bb") == 2


def test_tm90_cluster_capped_at_ten_percent_rounded_up():
    rows = [make_row(i, cluster="shared_c", ipsae=0.9 - 0.001 * i) for i in range(20)]
    rows += [make_row(100 + i, ipsae=0.4) for i in range(40)]
    res = run(rows)
    assert sum(1 for r in res.ranked if r["tm90_cluster_id"] == "shared_c") == 3


def test_structure_method_capped_at_fifty_percent_max_15_of_30():
    rows = [make_row(i, method="rfdiffusion", ipsae=0.9 - 0.001 * i)
            for i in range(40)]
    rows += [make_row(100 + i, method="boltzgen", ipsae=0.4) for i in range(10)]
    rows += [make_row(200 + i, method="genie3", ipsae=0.3) for i in range(10)]
    res = run(rows)
    assert sum(1 for r in res.ranked if r["structure_method"] == "rfdiffusion") == 15
    assert len(res.ranked) == 30


def test_seq_method_capped_at_two_thirds_with_backfill():
    """2/3 of 30 = 20 rows per seq_method; the rest backfills from alternates."""
    rows = [make_row(i, seq_method="solublempnn", ipsae=0.9 - 0.001 * i,
                     method=METHODS[i % 3])
            for i in range(40)]
    rows += [make_row(100 + i, seq_method="solublecaliby", ipsae=0.4,
                      method=METHODS[i % 3])
             for i in range(20)]
    res = run(rows)
    n_mpnn = sum(1 for r in res.ranked if r["seq_method"] == "solublempnn")
    assert n_mpnn == 20
    assert len(res.ranked) == 30


def test_fewer_than_three_structure_methods_is_a_logged_deviation():
    rows = [make_row(i, method="rfdiffusion" if i % 2 else "boltzgen", ipsae=0.5)
            for i in range(40)]
    res = run(rows)
    assert any(d["kind"] == "diversity_floor" for d in res.deviations)


# --- rank ordering ------------------------------------------------------------ #

def test_rank_order_is_full_seeds_then_pose_pass_then_rank_zscore():
    rows = [
        make_row(1, ipsae=0.50, dockq=0.50, n_seeds=5),   # full seeds, pose PASS
        make_row(2, ipsae=0.99, dockq=0.10, n_seeds=5),   # best z but pose FAIL
        make_row(3, ipsae=0.95, dockq=0.90, n_seeds=3),   # below-protocol seeds
        make_row(4, ipsae=0.80, dockq=0.60, n_seeds=5),   # full seeds, pose PASS
    ]
    res = run(rows)
    assert [r["design_id"] for r in res.ranked] == \
        ["d0004", "d0001", "d0002", "d0003"]


def test_below_protocol_seed_rows_stay_ranked_with_their_count_disclosed():
    rows = [make_row(i, ipsae=0.5, n_seeds=3) for i in range(5)]
    res = run(rows)
    assert len(res.ranked) == 5
    assert all(r["n_seeds"] == 3 for r in res.ranked)
    assert any(d["kind"] == "seed_shortfall" for d in res.deviations)


def test_ranks_are_contiguous_from_one():
    rows = [make_row(i, ipsae=0.9 - 0.01 * i, method=METHODS[i % 3])
            for i in range(40)]
    res = run(rows)
    assert [r["rank"] for r in res.ranked] == list(range(1, 31))


# --- unranked routing ---------------------------------------------------------- #

@pytest.mark.parametrize("mutate,expect", [
    ({"final_score": None}, "missing final_score"),
    ({"pose_dockq": None}, "pose check not run"),
    ({"root_backbone_id": None}, "missing root_backbone_id"),
    ({"tm90_cluster_id": None}, "missing tm90_cluster_id"),
    ({"fold_class": None}, "missing fold_class"),
    ({"novelty_verdict": "REJECT"}, "gate novelty_verdict"),
])
def test_rows_missing_a_required_result_go_to_the_unranked_set(mutate, expect):
    bad = make_row(1)
    bad.update(mutate)
    rows = [bad] + [make_row(10 + i) for i in range(3)]
    res = run(rows)
    assert bad["design_id"] not in {r["design_id"] for r in res.ranked}
    reason = next(r["unranked_reason"] for r in res.unranked
                  if r["design_id"] == bad["design_id"])
    assert expect in reason


def test_a_zero_scored_row_is_never_ranked():
    zero = make_row(1, ipsae=0.0, dockq=0.0)
    zero["rank_zscore"] = 0
    res = run([zero] + [make_row(10 + i) for i in range(3)])
    assert zero["design_id"] not in {r["design_id"] for r in res.ranked}


# --- write-time gate recompute -------------------------------------------------- #

def test_a_carried_value_that_disagrees_with_the_recompute_halts_the_writer():
    tampered = make_row(1)
    tampered["final_score"] = 0.99  # not what the six terms say
    with pytest.raises(ValueError, match="sheet writer HALT on row 'd0001'"):
        run([tampered])


def test_recompute_derives_pose_dockq_as_the_min_over_arms():
    row = make_row(1)
    row["sc_DockQ_ef2fast"] = 0.11
    out = recompute_row(row, MASK, GateThresholds())
    assert out["pose_dockq"] == pytest.approx(0.11)
    assert out["pose_PASS"] is False


def test_provenance_tokens_are_canonicalised_and_unknown_ones_rejected():
    aliased = make_row(1, method="RFdiffusion")
    res = run([aliased] + [make_row(10 + i) for i in range(3)])
    assert res.ranked[0]["structure_method"] in {"rfdiffusion"} or any(
        r["structure_method"] == "rfdiffusion" for r in res.ranked
    )

    unknown = make_row(2, method="my_secret_model")
    res2 = run([unknown] + [make_row(20 + i) for i in range(3)])
    assert unknown["design_id"] not in {r["design_id"] for r in res2.ranked}
    assert any("not in the frozen enum" in r["unranked_reason"]
               for r in res2.unranked)


# --- counter-screened targets ---------------------------------------------------- #

def test_a_null_offtarget_score_on_a_counter_screened_target_is_fail_loud():
    mask = InstrumentMask(counter_screened=True)
    row = make_row(1)
    with pytest.raises(ValueError, match="counter-screened"):
        run([row], target="GDF-8", mask=mask)


# --- relaxation ladder ------------------------------------------------------------ #

def test_diversity_caps_relax_first_and_the_step_is_recorded_on_each_row():
    # everything shares one root: strict cap admits 2, relaxed cap admits 8
    rows = [make_row(i, root="only_bb", ipsae=0.9 - 0.001 * i,
                     method=METHODS[i % 3])
            for i in range(40)]
    res = run(rows)
    assert "DIVERSITY_CAPS" in res.relaxation_steps
    assert len(res.ranked) == 8  # ceil(0.25 * 30)
    assert any(r["relaxation_step"] == "DIVERSITY_CAPS" for r in res.ranked)
    assert any(d["kind"] == "relaxation" for d in res.deviations)


def test_relaxation_never_goes_past_the_hard_bounds():
    rows = [make_row(i, root="only_bb", cluster="only_c", ipsae=0.5,
                     method=METHODS[i % 3])
            for i in range(40)]
    res = run(rows)
    per_root = sum(1 for r in res.ranked if r["root_backbone_id"] == "only_bb")
    assert per_root <= 8  # never past 25 % of 30


def test_a_short_sheet_ships_the_actual_n_and_is_never_padded():
    rows = [make_row(i, ipsae=0.5) for i in range(7)]
    res = run(rows)
    assert len(res.ranked) == 7
    assert len({r["sequence"] for r in res.ranked}) == 7
    assert any(d["kind"] == "short_sheet" for d in res.deviations)


def test_relaxation_is_not_attempted_before_regeneration_has_been_tried():
    rows = [make_row(i, root="only_bb", ipsae=0.5) for i in range(40)]
    res = run(rows, regeneration_tried=False)
    assert res.relaxation_steps == []
    assert len(res.ranked) == 2  # the strict 5 % cap


# --- reported (not gated) fold diversity ------------------------------------------- #

def test_fold_diversity_is_reported_not_enforced():
    rows = [make_row(i, fold_class="all_alpha", ipsae=0.9 - 0.001 * i,
                     method=METHODS[i % 3])
            for i in range(40)]
    res = run(rows)
    assert len(res.ranked) == 30  # still ships
    assert res.diagnostics["fold_diversity_target_met"] is False
    assert res.diagnostics["n_non_all_alpha"] == 0


def test_fold_diversity_report_uses_the_ten_percent_target():
    rows = [{"fold_class": "not_all_alpha"}] * 3 + [{"fold_class": "all_alpha"}] * 27
    rep = fold_diversity_report(rows)
    assert rep["fraction_non_all_alpha"] == pytest.approx(0.10)
    assert rep["fold_diversity_target_met"] is True


def test_selection_caps_round_up():
    caps = SelectionCaps()
    assert caps.cap(0.05, 30) == 2
    assert caps.cap(0.10, 30) == 3
    assert caps.cap(0.50, 30) == 15
