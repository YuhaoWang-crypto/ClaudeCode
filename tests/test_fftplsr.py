"""Tests for the fftplsr package.

Run with pytest, or directly: ``python3 tests/test_fftplsr.py``.

The load-bearing ones are:

* `test_fast_encoder_matches_reference` -- the rank-k spectrum update is exact, so
  the fast path used everywhere is not a separate approximate implementation.
* `test_pls_path_matches_per_k_refits` -- slicing one PLS fit reproduces independent
  refits, which is what makes the 566-descriptor screen affordable.
* `test_dc_bin_is_roundoff` / `test_drop_dc_makes_features_precision_independent` --
  the M2 finding, pinned so a future change cannot quietly reintroduce it.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneOut

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from fftplsr import datasets, variants  # noqa: E402
from fftplsr.baselines import LogAdditiveModel, MeanPredictor, OneHotRidge  # noqa: E402
from fftplsr.encode import MutationEncoder, encode_sequences, spectrum  # noqa: E402
from fftplsr.model import cv_score, pls_path_predict  # noqa: E402
from fftplsr.variants import VariantError  # noqa: E402

CODES = ["OOBM850103", "QIAN880114", "RADA880104"]


# --------------------------------------------------------------------- variants


def test_parse_single_and_multi():
    assert variants.parse_variant("D2N") == (variants.Substitution("D", 2, "N"),)
    assert len(variants.parse_variant("D2N/K3N/H63Y")) == 3
    assert variants.parse_variant("IFRS") == ()
    assert variants.parse_variant("Com1-IFRS") == ()


def test_apply_variant_validates_wild_type():
    parent = datasets.ifrs()
    assert variants.apply_variant(parent, "D2N")[1] == "N"
    assert variants.apply_variant(parent, "IFRS") == parent
    # The parent has D at position 2, not A -- the original code printed and carried on.
    with pytest.raises(VariantError, match="states wild-type"):
        variants.apply_variant(parent, "A2N")
    with pytest.raises(VariantError, match="outside the parent"):
        variants.apply_variant(parent, "D9999N")


def test_duplicate_position_rejected():
    with pytest.raises(VariantError, match="same position twice"):
        variants.parse_variant("H63Y/H63L")


def test_enumerate_combinations_counts():
    parent = datasets.ifrs()
    sites = datasets.ROUND12_SITES
    combos = variants.enumerate_combinations(parent, sites, min_order=2)
    # 2^12 subsets, minus the empty set and the 12 singletons.
    assert len(combos) == 2 ** len(sites) - 1 - len(sites) == 4083
    assert all(len(variants.parse_variant(c)) >= 2 for c in combos)
    assert len(set(combos)) == len(combos)


def test_enumerate_skips_same_position_clashes():
    parent = datasets.com1()
    # Two substitutions at 63 can never co-occur.
    combos = variants.enumerate_combinations(parent, ["H63L", "H63A", "K67N"], min_order=2)
    assert "H63A/K67N" in combos and "H63L/K67N" in combos
    assert not any(c.count("63") == 2 for c in combos)


def test_saturation_scan_excludes_self_substitution():
    parent = datasets.com1()
    scan = variants.saturation_scan(parent, positions=[63, 67])
    assert len(scan) == 19 * 2  # not 20 * 2
    assert f"{parent[62]}63{parent[62]}" not in scan


# ----------------------------------------------------------------------- encode


def _round1():
    parent = datasets.ifrs()
    table, _ = datasets.trainset(1)
    labels = table["Variants"].tolist()
    y = table["Fitness"].to_numpy(dtype=float)
    encoder = MutationEncoder(
        parent, positions=variants.mutated_positions(labels + datasets.ROUND12_SITES)
    )
    return parent, labels, y, encoder


def test_fast_encoder_matches_reference():
    parent, labels, _, encoder = _round1()
    space = variants.enumerate_combinations(parent, datasets.ROUND12_SITES, min_order=2)
    query = labels + space[:150]
    fast = encoder.encode(query, CODES)
    ref = encode_sequences([variants.apply_variant(parent, q) for q in query], CODES)
    assert fast.shape == ref.shape
    assert np.abs(fast - ref).max() < 1e-12


def test_feature_width():
    _, labels, _, encoder = _round1()
    # 454 residues -> 227 half-spectrum bins, minus the dropped DC bin, per descriptor.
    assert encoder.encode(labels, CODES).shape == (len(labels), 3 * 226)
    assert encoder.encode(labels, CODES, drop_dc=False).shape == (len(labels), 3 * 227)


def test_dc_bin_is_roundoff():
    _, labels, _, encoder = _round1()
    kept = encoder.encode(labels, ["OOBM850103"], drop_dc=False)
    # Mean-centering zeroes the DC term, so bin 0 is round-off, not signal.
    assert np.abs(kept[:, 0]).max() < 1e-12
    dropped = encoder.encode(labels, ["OOBM850103"], drop_dc=True)
    assert np.abs(dropped - kept[:, 1:]).max() == 0.0


def test_drop_dc_makes_scoring_precision_independent():
    """The M2 finding, pinned: bin 0 turns FFT round-off into a fitted predictor.

    Re-drawing bin 0 anywhere inside its round-off range is what a different FFT
    backend, precision or summation order does. With the bin kept, that moves the
    cross-validated error enough to change which descriptor wins; with it dropped,
    the error is invariant.
    """
    _, labels, y, encoder = _round1()
    rng = np.random.default_rng(0)

    kept = encoder.encode(labels, ["OOBM850103"], drop_dc=False)
    scores = []
    for _ in range(30):
        X = kept.copy()
        X[:, 0] = rng.uniform(0.0, 6e-16, size=len(y))
        scores.append(cv_score(X, y, cv=None).cv_mse)
    assert np.ptp(scores) > 0.1, f"expected DC-bin instability, spread was {np.ptp(scores)}"

    dropped = encoder.encode(labels, ["OOBM850103"], drop_dc=True)
    base = cv_score(dropped, y, cv=None).cv_mse
    for _ in range(10):
        X = dropped.copy()
        X[:, 0] += rng.uniform(-6e-16, 6e-16, size=len(y))
        assert abs(cv_score(X, y, cv=None).cv_mse - base) < 1e-9


def test_spectrum_handles_constant_descriptor():
    assert np.count_nonzero(spectrum(np.ones(20))) == 0


def test_candidate_pool_excludes_incomplete_indices():
    """Entries missing residue values must leave the pool loudly, not silently.

    `aaindex` 1.0.5 supplies all 566; 1.3.2 returns None for 13 of them, so the
    counts are asserted relationally rather than as fixed numbers.
    """
    from fftplsr.encode import all_indices, available_indices, incomplete_indices

    assert len(all_indices()) == 566
    assert set(available_indices()) == set(all_indices()) - set(incomplete_indices())
    assert set(available_indices(complete_only=False)) == set(all_indices())
    assert len(available_indices()) + len(incomplete_indices()) == 566


def test_incomplete_index_raises_a_named_error():
    from fftplsr.encode import IncompleteIndexError, _lookup_table, incomplete_indices

    _, labels, _, encoder = _round1()
    for code in incomplete_indices():
        with pytest.raises(IncompleteIndexError, match="no value for"):
            _lookup_table(code)
        with pytest.raises(IncompleteIndexError):
            encoder.encode(labels, [code])
    # Whatever the package version reports as usable must actually encode.
    from fftplsr.encode import available_indices

    for code in available_indices()[:25]:
        assert np.isfinite(encoder.encode(labels, [code])).all()


# ------------------------------------------------------------------------ model


def test_pls_path_matches_per_k_refits():
    _, labels, y, encoder = _round1()
    X = encoder.encode(labels, ["RADA880104"])
    train, test = np.arange(0, 10), np.arange(10, len(y))
    path = pls_path_predict(X[train], y[train], X[test], 8)
    for k in range(1, 9):
        ref = PLSRegression(n_components=k, scale=True).fit(X[train], y[train]).predict(X[test])
        assert np.abs(path[k - 1] - ref.ravel()).max() < 1e-8, f"mismatch at k={k}"


def test_cv_score_picks_the_cv_minimum():
    _, labels, y, encoder = _round1()
    X = encoder.encode(labels, ["RADA880104"])
    best = cv_score(X, y, cv=None)
    oof = np.empty((10, len(y)))
    for train, test in LeaveOneOut().split(X):
        oof[:, test] = pls_path_predict(X[train], y[train], X[test], 10)
    mse = ((y[None, :] - oof) ** 2).mean(axis=1)
    assert best.n_components == int(np.argmin(mse[1:])) + 2  # the grid starts at 2
    assert best.cv_mse == pytest.approx(mse[best.n_components - 1])


def test_cv_score_reports_out_of_fold_not_in_sample():
    """The original implementation reported in-sample R2; these must differ here."""
    _, labels, y, encoder = _round1()
    res = cv_score(encoder.encode(labels, ["RADA880104"]), y, cv=None)
    assert res.train_r2 > res.cv_r2
    assert res.train_mse < res.cv_mse


# -------------------------------------------------------------------- baselines


def test_log_additive_multiplies_measured_singles():
    measured = {"WT": 1.0, "D2N": 2.0, "H62Y": 3.0}
    mdl = LogAdditiveModel().fit(list(measured), np.array(list(measured.values())))
    assert mdl.predict(["D2N/H62Y"])[0] == pytest.approx(6.0)
    assert mdl.predict(["D2N"])[0] == pytest.approx(2.0)
    assert mdl.predict(["WT"])[0] == pytest.approx(1.0)


def test_log_additive_flags_unmeasured_singles():
    mdl = LogAdditiveModel().fit(["WT", "D2N"], np.array([1.0, 2.0]))
    mdl.predict(["D2N/K3N"])
    assert mdl.n_missing == 1


def test_mean_and_onehot_baselines_run():
    _, labels, y, _ = _round1()
    assert MeanPredictor().fit(labels, y).predict(labels) == pytest.approx(y.mean())
    pred = OneHotRidge().fit(labels, y).predict(labels)
    assert pred.shape == y.shape and np.isfinite(pred).all()


# --------------------------------------------------------------------- datasets


def test_all_dataset_labels_are_consistent_with_their_parent():
    """Every shipped label must apply cleanly to its background sequence."""
    for n in (1, 2, 3, 4):
        table, parent = datasets.trainset(n)
        for label in table["Variants"]:
            variants.apply_variant(parent, label)  # raises on any inconsistency
    for background in ("IFRS", "Com1-IFRS"):
        table, parent = datasets.measured_panel(background)
        for label in table["Variants"]:
            variants.apply_variant(parent, label)


def test_dataset_shapes_and_parent_reference():
    assert len(datasets.ifrs()) == len(datasets.com1()) == 454
    assert [len(datasets.trainset(n)[0]) for n in (1, 2, 3, 4)] == [13, 38, 96, 120]
    assert len(datasets.measured_panel("IFRS")[0]) == 102
    assert len(datasets.measured_panel("Com1-IFRS")[0]) == 140
    for n in (1, 2, 3, 4):
        table, _ = datasets.trainset(n)
        parent_rows = [
            f for f, v in zip(table["Fitness"], table["Variants"]) if not variants.parse_variant(v)
        ]
        assert parent_rows == [1.0], f"trainset {n} parent should be the 1.0 reference"


def test_com1_is_seven_mutations_from_ifrs():
    """Com1-IFRS carries 7 substitutions; pinned so the rank metric targets the right label."""
    from fftplsr.m1_reproduce import COM1

    assert COM1 == "D2N/V31I/T56P/R61K/H62Y/T122S/S193R"
    assert len(variants.parse_variant(COM1)) == 7
    assert variants.apply_variant(datasets.ifrs(), COM1) == datasets.com1()
    assert set(variants.mutated_positions([COM1])) <= set(
        variants.mutated_positions(datasets.ROUND12_SITES)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
