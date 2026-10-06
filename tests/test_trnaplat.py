"""Tests for the trnaplat demo layer.

Run with pytest, or directly: ``python3 tests/test_trnaplat.py``.

The load-bearing ones are:

* `test_com1_transplants_onto_wild_type_mm_pylrs` -- Demo 1's whole premise is
  that the activity mutations are defined on the AzK target with no remapping.
  If that stops being true, the demo is built on nothing.
* `test_check_transplant_refuses_wrong_target` -- it refuses rather than
  silently mis-numbering, the same guard as `pylrs/literature.verify_numbering`.
* `test_crosscheck_catches_disagreement` -- the one cross-check the hand-curated
  chemotype table has (pPRF and libY-C are the same molecule) actually fires.
* `test_subst_class_is_id_like` -- pins the finding that made the categorical
  gate testable: a column with one level per molecule gates nothing.
* `test_chemotype_block_is_constant_within_ncaa` -- the structural reason the
  chemotype block scores exactly 0.5 under leave-one-ncAA-out.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "pylrs"))

from trnaplat import demo1_libY_control as control  # noqa: E402
from trnaplat import demo1_multisite as d1  # noqa: E402
from trnaplat import demo3_evidence as d3  # noqa: E402
from trnaplat import ncaa_chemotype as nc  # noqa: E402


# ------------------------------------------------------------------ Demo 1

def test_com1_transplants_onto_wild_type_mm_pylrs():
    table = d1.check_transplant(d1.MM_PYLRS_WT, d1.COM1)
    assert len(table) == 7
    assert table["defined"].all()


def test_six_of_seven_activity_sites_are_outside_the_pocket():
    table = d1.check_transplant(d1.MM_PYLRS_WT, d1.COM1)
    outside = (~table["domain"].str.startswith("catalytic")).sum()
    assert outside == 6, "the ncAA-agnostic argument rests on this count"


def test_ifrs_differs_from_wild_type_only_in_the_catalytic_domain():
    from fftplsr import datasets
    ifrs = datasets.ifrs()
    diff = [i + 1 for i in range(len(ifrs)) if ifrs[i] != d1.MM_PYLRS_WT[i]]
    assert diff == [346, 348]
    assert all(d1.domain_of(p).startswith("catalytic") for p in diff)


def test_check_transplant_refuses_wrong_target():
    scrambled = "A" * len(d1.MM_PYLRS_WT)
    with pytest.raises(d1.TransplantError) as err:
        d1.check_transplant(scrambled, d1.COM1, target_name="scrambled")
    assert "not defined" in str(err.value)


def test_check_transplant_refuses_position_past_the_end():
    with pytest.raises(d1.TransplantError):
        d1.check_transplant(d1.MM_PYLRS_WT[:100], d1.COM1)


def test_domain_boundaries_are_contiguous_and_cover_the_protein():
    assert d1.domain_of(1).startswith("N-terminal")
    assert d1.domain_of(149).startswith("N-terminal")
    assert d1.domain_of(150) == "linker"
    assert d1.domain_of(185).startswith("catalytic")
    assert d1.domain_of(454).startswith("catalytic")
    assert d1.domain_of(455) == "out of range"


def test_multiplicative_null_is_p_to_the_n():
    frame = d1.multiplicative_null(0.5, 4)
    assert frame["relative_yield"].tolist() == [1.0, 0.5, 0.25, 0.125, 0.0625]


def test_compounding_gain_is_the_ratio_to_the_n():
    frame = d1.compounding_table(0.5, 0.75, 3)
    np.testing.assert_allclose(frame["fold_gain"].to_numpy(),
                               [1.0, 1.5, 2.25, 3.375])


def test_fit_recovers_p_and_zero_curvature_from_clean_multiplicative_data():
    n = np.arange(5)
    y = 0.6 ** n
    fit = d1.fit_site_efficiency(n, y)
    assert fit["p_per_site"] == pytest.approx(0.6, abs=1e-9)
    assert fit["curvature_b"] == pytest.approx(0.0, abs=1e-9)
    assert fit["mechanism"].startswith("sites independent")


def test_fit_recovers_curvature_and_calls_the_shared_resource():
    n = np.arange(1, 6)
    y = np.exp(np.log(0.6) * n - 0.09 * n ** 2)
    fit = d1.fit_site_efficiency(n, y)
    assert fit["curvature_b"] == pytest.approx(-0.09, abs=1e-6)
    assert fit["p_per_site_quadratic"] == pytest.approx(0.6, abs=1e-6)
    assert fit["mechanism"].startswith("shared resource")


def test_fit_needs_three_points_to_separate_slope_from_curvature():
    with pytest.raises(ValueError, match="at least 3"):
        d1.fit_site_efficiency([1, 2], [0.5, 0.25])


def test_power_table_reports_a_false_positive_row():
    power = d1.power_table(0.55, [1, 2, 3], curvatures=(0.0, -0.08),
                           replicates=(3,), n_sim=200)
    fp = power[power["true_b"] == 0.0].iloc[0]
    assert fp["kind"] == "FALSE POSITIVE"
    assert 0.0 <= fp["calls_tRNA_limited"] < 0.5


def test_construct_series_has_a_denominator_and_both_controls():
    plate = d1.construct_series([0, 1, 2], ["arm A"], replicates=2)
    assert (plate["n_sites"] == 0).any(), "no n=0 arm means no denominator"
    roles = set(plate["role"])
    assert any(r.startswith("control: truncation") for r in roles)
    assert any(r.startswith("control: background") for r in roles)


def test_ingest_refuses_an_arm_without_an_n0_denominator(tmp_path):
    frame = pd.DataFrame({"arm": ["a"] * 3, "n_sites": [1, 2, 3],
                          "yield": [100.0, 50.0, 25.0]})
    path = tmp_path / "y.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="no n=0"):
        d1.ingest(path)


def test_ingest_refuses_an_unmatched_baseline_arm(tmp_path):
    frame = pd.DataFrame({"arm": ["a"] * 4, "n_sites": [0, 1, 2, 3],
                          "yield": [100.0, 60.0, 36.0, 21.6]})
    path = tmp_path / "y.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="matched none"):
        d1.ingest(path, baseline_arm="does-not-exist")


def test_ingest_uses_file_order_not_alphabetical_for_the_default_baseline(tmp_path):
    rows = []
    for arm in ["zebra", "alpha"]:          # deliberately not alphabetical
        for n, y in [(0, 100.0), (1, 60.0), (2, 36.0), (3, 21.6)]:
            rows.append({"arm": arm, "n_sites": n, "yield": y})
    path = tmp_path / "y.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    fits = d1.ingest(path)
    assert fits[fits["is_baseline"]].iloc[0]["arm"] == "zebra"


# -------------------------------------------------------- chemotype table

def test_crosscheck_passes_on_the_shipped_table():
    nc.verify_crosscheck()          # must not raise


def test_crosscheck_catches_disagreement():
    frame = nc.full_table()
    frame.loc[frame["ncAA"] == "pPRF", "n_heavy"] = 99
    with pytest.raises(ValueError, match="same molecule"):
        nc.verify_crosscheck(frame)


def test_full_table_rejects_duplicates(monkeypatch):
    monkeypatch.setattr(nc, "AZK", [("pAzF", "Lys", "x", "y", 1, 0, 0, 0, 0, 0)])
    with pytest.raises(ValueError, match="duplicate"):
        nc.full_table()


def test_distance_matrix_is_symmetric_with_zero_diagonal():
    dist = nc.distance_matrix()
    np.testing.assert_allclose(dist.to_numpy(), dist.to_numpy().T, atol=1e-12)
    np.testing.assert_allclose(np.diag(dist.to_numpy()), 0.0, atol=1e-12)


def test_nearest_rejects_an_unknown_ncaa():
    with pytest.raises(KeyError, match="not in the descriptor table"):
        nc.nearest("not-a-real-ncAA")


def test_product_and_chemistry_groups_only_name_known_ncaas():
    known = set(nc.table()["ncAA"])
    for groups in (nc.PRODUCT_GROUPS, nc.CHEMISTRY_GROUPS):
        for name, spec in groups.items():
            unknown = set(spec["members"]) - known
            assert not unknown, f"{name} names unknown ncAA(s) {unknown}"


# ------------------------------------------------------------------ Demo 3

def test_subst_class_is_id_like_and_excluded_from_the_gate():
    usable, id_like = d3.gateable_columns()
    assert id_like == ["subst_class"]
    assert usable == ["parent", "attach"]


def test_unseen_levels_counts_only_gateable_columns():
    others = [n for n in nc.table()["ncAA"] if n != "2-NPA"]
    # 2-NPA is the only ortho-substituted entry, so `attach` is unseen when it
    # is held out, while `parent` (Phe) is not.
    assert d3.unseen_levels("2-NPA", others) == 1
    # pAzF shares both gateable levels with several others.
    assert d3.unseen_levels("pAzF", [n for n in nc.table()["ncAA"]
                                     if n != "pAzF"]) == 0


def test_chemotype_block_is_constant_within_one_ncaa():
    frame = pd.DataFrame({"ncAA": ["pAzF"] * 4,
                          "mutations": ["Y32A", "Y32G", "D158T", "L162Q"]})
    block = d3.chemotype_block(frame)
    assert np.ptp(block, axis=0).max() == 0.0


def test_chemotype_block_rejects_an_undescribed_ncaa():
    frame = pd.DataFrame({"ncAA": ["mystery-AA"], "mutations": ["Y32A"]})
    with pytest.raises(KeyError, match="no chemotype descriptors"):
        d3.chemotype_block(frame)


def test_interaction_block_is_the_flattened_outer_product():
    onehot = np.array([[1.0, 0.0], [0.0, 1.0]])
    chem = np.array([[2.0, 3.0], [4.0, 5.0]])
    out = d3.interaction_block(onehot, chem)
    assert out.shape == (2, 4)
    np.testing.assert_allclose(out[0], [2.0, 3.0, 0.0, 0.0])
    np.testing.assert_allclose(out[1], [0.0, 0.0, 4.0, 5.0])


def test_auc_standard_error_shrinks_with_more_negatives():
    wide = d3.auc_standard_error(0.8, 7, 10)
    tight = d3.auc_standard_error(0.8, 7, 100)
    assert tight < wide


def test_evidence_tiers_are_monotone_and_total():
    limits = [limit for _, limit, _ in d3.EVIDENCE_TIERS]
    assert limits == sorted(limits)
    assert limits[-1] == float("inf"), "the tiers must cover every distance"
    assert d3.assign_tier(0.0).startswith("A")
    assert d3.assign_tier(1.5).startswith("B")
    assert d3.assign_tier(2.5).startswith("C")
    assert d3.assign_tier(99.0).startswith("D")


def test_onehot_block_marks_wild_type_where_unmutated():
    wt = "A" * 200
    frame = pd.DataFrame({"mutations": ["A32G"]})
    block = d3.onehot_block(frame, [32, 158], wt)
    assert block[0, d3.AA20.index("G")] == 1.0            # position 32 -> G
    assert block[0, len(d3.AA20) + d3.AA20.index("A")] == 1.0   # 158 stays A


# ------------------------------------------------- libY out-of-domain control

def test_azk_is_out_of_libY_domain_on_every_categorical_column():
    frame = nc.full_table()
    liby = [n for n in frame["ncAA"] if n.startswith("libY-")]
    coverage = control.categorical_coverage(frame, "AzK", liby)
    assert len(coverage) == 3
    assert not coverage["in_domain"].any(), "all three levels must be unseen"


def test_libY_training_substrates_are_one_chemotype_cell():
    frame = nc.full_table()
    liby = frame[frame["ncAA"].str.startswith("libY-")]
    assert set(liby["parent"]) == {"Tyr"}
    assert set(liby["attach"]) == {"para_O"}


def test_only_one_feature_block_is_substrate_aware():
    aware = [b for b in control.FEATURE_BLOCKS if b[3]]
    assert len(aware) == 1
    assert aware[0][1] == 18
    assert sum(b[1] for b in control.FEATURE_BLOCKS) == 280


def test_distance_fails_to_separate_azk_from_libY_training_spread():
    """Pins the ❌ finding: the metric gate does not detect a parent change."""
    frame = nc.full_table()
    dist = nc.distance_matrix(frame)
    liby = [n for n in frame["ncAA"] if n.startswith("libY-")]
    import itertools
    within = [dist.loc[a, b] for a, b in itertools.combinations(liby, 2)]
    azk_gap = dist.loc["AzK", liby].min()
    assert azk_gap < max(within), (
        "if this ever fails, Euclidean distance has started detecting the "
        "out-of-domain case and demo1_libY_control's conclusion needs revisiting"
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
