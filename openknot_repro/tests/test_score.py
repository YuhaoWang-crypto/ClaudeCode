"""Tests for the OpenKnot score implementation.

Run:  python -m pytest openknot_repro/tests -q
The last test needs the benchmark CSV and is skipped when it is absent.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openknot.score import (
    bp_list,
    crossed_pair_scores,
    crossing_residues,
    embed_target,
    eterna_classic_score,
    filter_singlet_pairs,
    get_helices,
    openknot_score,
)


def test_bp_list_handles_all_bracket_types():
    assert bp_list("((..))") == [(0, 5), (1, 4)]
    assert bp_list("([{<..>}])") == [(0, 9), (1, 8), (2, 7), (3, 6)]


def test_bp_list_is_sorted_by_opening_index():
    # helix grouping depends on this ordering
    pairs = bp_list("(((...)))")
    assert pairs == sorted(pairs)


def test_bp_list_drops_unbalanced_closing_bracket():
    assert bp_list("((..)))") == [(0, 5), (1, 4)]


def test_crossing_residues_finds_a_pseudoknot():
    # stem A crosses stem B
    assert crossing_residues(bp_list("((..[[..))..]]")) == {0, 1, 4, 5, 8, 9, 12, 13}


def test_crossing_residues_empty_for_nested_structure():
    assert crossing_residues(bp_list("((((....))))")) == set()


def test_get_helices_groups_a_stacked_run():
    helices = get_helices("(((...)))")
    assert len(helices) == 1 and len(helices[0]) == 3


def test_filter_singlet_pairs_removes_lone_pairs():
    # one 3-bp helix plus an isolated pair
    structure = "(((...)))..(.)"
    kept = filter_singlet_pairs(structure)
    assert (11, 13) not in kept
    assert len(kept) == 3


def test_ecs_thresholds():
    structure = "(()."
    data = [0.49, 0.10, 0.51, 0.13]
    # paired hits below 0.5:      0.49 hit, 0.10 hit, 0.51 miss
    # unpaired hits above 0.125:  0.13 hit
    assert eterna_classic_score(structure, data, 0, 3) == pytest.approx(75.0)


def test_ecs_threshold_comparisons_are_strict():
    # a residue sitting exactly on its threshold is a miss, both ways
    assert eterna_classic_score("(", [0.5], 0, 0) == 0.0
    assert eterna_classic_score(".", [0.125], 0, 0) == 0.0


def test_ecs_counts_missing_data_as_a_miss():
    structure = "...."
    data = [0.5, 0.5, float("nan"), 0.5]
    # NaN stays in the denominator
    assert eterna_classic_score(structure, data, 0, 3) == pytest.approx(75.0)


def test_ecs_scoring_region_is_respected():
    structure = "....."
    data = [0.0, 0.5, 0.5, 0.5, 0.0]  # ends would miss
    assert eterna_classic_score(structure, data, 1, 3) == pytest.approx(100.0)


def test_cpq_quality_counts_only_crossed_residues():
    structure = "((..[[..))..]]"
    protected = [0.1] * len(structure)
    _, cpq = crossed_pair_scores(structure, protected, 0, len(structure) - 1)
    assert cpq == pytest.approx(100.0)
    reactive = [0.9] * len(structure)
    _, cpq = crossed_pair_scores(structure, reactive, 0, len(structure) - 1)
    assert cpq == pytest.approx(0.0)


def test_cpq_is_zero_without_crossed_pairs():
    structure = "((((....))))"
    _, cpq = crossed_pair_scores(structure, [0.1] * len(structure), 0, len(structure) - 1)
    assert cpq == 0.0


def test_openknot_score_is_the_mean_of_its_halves():
    structure = "((..[[..))..]]"
    data = [0.1, 0.9, 0.3, 0.3, 0.1, 0.1, 0.3, 0.3, 0.1, 0.1, 0.3, 0.3, 0.1, 0.1]
    ecs = eterna_classic_score(structure, data, 0, 13)
    cpq = crossed_pair_scores(structure, data, 0, 13)[1]
    assert openknot_score(structure, data, 0, 13) == pytest.approx(0.5 * ecs + 0.5 * cpq)


def test_embed_target_places_the_design_inside_its_pads():
    assert embed_target("((..))", 3, 8, 10) == "..((..)).."


def test_embed_target_rejects_a_span_mismatch():
    with pytest.raises(ValueError):
        embed_target("((..))", 3, 9, 12)


@pytest.mark.parametrize("n_rows", [2000])
def test_matches_released_scores_on_the_benchmark(n_rows):
    """Regression against the released target_openknot_score column."""
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")
    from openknot.data import BENCH_CSV, load_benchmark

    if not BENCH_CSV.exists():
        pytest.skip("benchmark CSV not downloaded")

    meta, reactivity = load_benchmark(nrows=n_rows)
    exact = 0
    scored = 0
    for i, row in enumerate(meta.itertuples(index=False)):
        if not isinstance(row.target_structure, str):
            continue
        start, end = int(row.sub_start) - 1, int(row.sub_end) - 1
        length = len(row.sequence)
        try:
            full = embed_target(row.target_structure, start + 1, end + 1, length)
        except ValueError:
            continue
        value = openknot_score(full, reactivity[i][:length], start, end)
        scored += 1
        if math.isclose(value, row.target_openknot_score, abs_tol=1e-6):
            exact += 1
    assert scored > 0
    # Residuals are single residues whose released reactivity was rounded onto a
    # threshold; the release-wide exact rate is 93.6%.
    assert exact / scored > 0.9
