"""Tests for the properties this project's conclusions actually depend on.

These are not coverage tests.  Each one pins a claim made in ``docs/`` to a
computation, so that a refactor which quietly breaks the claim fails here rather
than in a paragraph nobody re-runs.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from mri_cd8_til.experiments.simulate import (
    SimulationConfig,
    effect_size_for_auc,
    simulate_cohort,
    theoretical_feature_auc,
)
from mri_cd8_til.labels.schema import Immunophenotype, LabelTask, to_binary
from mri_cd8_til.labels.spatial import (
    PhenotypeThresholds,
    SpatialTILMetrics,
    compute_spatial_metrics,
)
from mri_cd8_til.labels.transcriptome import SIGNATURES, score_signatures
from mri_cd8_til.modeling.combat import ComBat, measure_batch_effect
from mri_cd8_til.modeling.metrics import (
    auc_interval_from_reported,
    calibration,
    decision_curve,
)
from mri_cd8_til.modeling.pipeline import (
    UnivariateAUCSelector,
    build_pipeline,
    columnwise_auc,
    run_reference_protocol,
)
from mri_cd8_til.modeling.stability import icc_2_1


# --------------------------------------------------------------- fast AUC


def test_columnwise_auc_matches_sklearn_including_ties():
    rng = np.random.default_rng(0)
    X = np.round(rng.standard_normal((80, 40)), 1)  # rounding forces ties
    y = (rng.random(80) < 0.4).astype(int)
    expected = np.array([roc_auc_score(y, X[:, j]) for j in range(X.shape[1])])
    assert np.allclose(columnwise_auc(X, y), expected)


def test_columnwise_auc_handles_constant_columns():
    X = np.zeros((20, 3))
    X[:, 1] = np.arange(20)
    y = np.r_[np.zeros(10), np.ones(10)].astype(int)
    auc = columnwise_auc(X, y)
    assert auc[0] == 0.5 and auc[2] == 0.5
    assert auc[1] > 0.9


# ------------------------------------------------------------ the core claim


def test_selection_before_splitting_inflates_auc_on_pure_noise():
    """The headline claim: with no signal, the published protocol still scores high."""
    config = SimulationConfig(n_patients=182, n_features=2000, n_informative=0, effect_size=0.0)
    leaky, clean = [], []
    for seed in range(12):
        cohort = simulate_cohort(SimulationConfig(**{**config.__dict__, "random_state": seed}))
        train, valid = cohort.train_validation_split(random_state=seed)
        leaky.append(
            run_reference_protocol(
                cohort.X, cohort.y, train, valid, select_on_all_data=True
            ).validation_auc
        )
        clean.append(
            run_reference_protocol(
                cohort.X, cohort.y, train, valid, select_on_all_data=False
            ).validation_auc
        )
    assert np.median(clean) == pytest.approx(0.5, abs=0.12)
    assert np.median(leaky) > 0.65
    assert np.median(leaky) - np.median(clean) > 0.15


def test_pipeline_selection_is_refitted_per_fold():
    """The leakage-free pipeline must not retain a selection across folds."""
    rng = np.random.default_rng(1)
    # More features than the selector keeps, otherwise "select everything"
    # trivially produces identical supports and the test proves nothing.
    X, y = rng.standard_normal((60, 300)), (rng.random(60) < 0.5).astype(int)
    pipeline = build_pipeline()
    first = pipeline.fit(X[:40], y[:40]).named_steps["select"].get_support(indices=True)
    second = pipeline.fit(X[20:], y[20:]).named_steps["select"].get_support(indices=True)
    assert not np.array_equal(first, second)


def test_univariate_selector_keeps_anticorrelated_features():
    """Selection must rank on |AUC - 0.5|, not on AUC."""
    rng = np.random.default_rng(2)
    X = rng.standard_normal((100, 5))
    y = (rng.random(100) < 0.5).astype(int)
    X[:, 3] -= 3.0 * y  # strongly anti-correlated
    support = UnivariateAUCSelector(k=1).fit(X, y).get_support(indices=True)
    assert support.tolist() == [3]


# ------------------------------------------------------------------ ComBat


def test_combat_removes_batch_variance():
    rng = np.random.default_rng(3)
    X = rng.standard_normal((200, 50))
    batch = rng.integers(0, 4, 200)
    for b in range(4):
        X[batch == b] = X[batch == b] * (1.0 + 0.4 * b) + 2.0 * b

    before = measure_batch_effect(X, batch)
    after = measure_batch_effect(ComBat().fit_transform(X, batch=batch), batch)
    assert before.median_eta_squared > 0.3
    assert after.median_eta_squared < 0.02


def test_combat_refuses_asymmetric_covariate_use():
    """Fitting with covariates and transforming without them would leak."""
    rng = np.random.default_rng(4)
    X = rng.standard_normal((60, 10))
    batch = rng.integers(0, 3, 60)
    covariates = rng.standard_normal((60, 1))
    harmoniser = ComBat().fit(X, batch=batch, covariates=covariates)
    with pytest.raises(ValueError, match="leak"):
        harmoniser.transform(X, batch=batch)


def test_combat_unseen_batch_policy_error_raises():
    rng = np.random.default_rng(5)
    X = rng.standard_normal((60, 10))
    batch = np.repeat([0, 1], 30)
    harmoniser = ComBat(unseen_batch="error").fit(X, batch=batch)
    with pytest.raises(ValueError, match="not seen during fit"):
        harmoniser.transform(X[:10], batch=np.full(10, 99))


def test_combat_handles_unseen_batch_under_reference_policy():
    rng = np.random.default_rng(6)
    X = rng.standard_normal((90, 10))
    batch = np.repeat([0, 1, 2], 30)
    harmoniser = ComBat(unseen_batch="reference").fit(X[:60], batch=batch[:60])
    out = harmoniser.transform(X[60:], batch=batch[60:])
    assert out.shape == (30, 10) and np.isfinite(out).all()


# -------------------------------------------------------------- metrics


def test_hanley_mcneil_interval_widens_as_n_shrinks():
    wide = auc_interval_from_reported(0.985, 15, 30)
    narrow = auc_interval_from_reported(0.985, 150, 300)
    assert wide.width > narrow.width
    assert wide.lower < 0.95  # the published number cannot exclude a much weaker model


def test_calibration_detects_overconfident_predictions():
    rng = np.random.default_rng(7)
    y = (rng.random(400) < 0.4).astype(int)
    honest = np.where(y == 1, rng.uniform(0.5, 0.8, 400), rng.uniform(0.2, 0.5, 400))
    overconfident = np.clip((honest - 0.5) * 6 + 0.5, 0.001, 0.999)
    assert calibration(y, honest).slope > calibration(y, overconfident).slope


def test_decision_curve_bounds():
    rng = np.random.default_rng(8)
    y = (rng.random(200) < 0.3).astype(int)
    p = np.clip(0.3 + 0.4 * y + rng.normal(0, 0.1, 200), 0.01, 0.99)
    curve = decision_curve(y, p)
    assert curve["model"].shape == curve["threshold"].shape
    assert np.all(curve["treat_none"] == 0)


# -------------------------------------------------------------- stability


def test_icc_is_one_for_perfectly_repeated_measurements():
    rng = np.random.default_rng(9)
    subjects = rng.standard_normal((30, 1, 8))
    measurements = np.repeat(subjects, 3, axis=1)
    assert np.allclose(icc_2_1(measurements), 1.0, atol=1e-6)


def test_icc_is_low_for_pure_noise():
    rng = np.random.default_rng(10)
    assert icc_2_1(rng.standard_normal((30, 3, 8))).mean() < 0.3


# ---------------------------------------------------------------- labels


def test_inflamed_case_is_excluded_from_the_peritumoural_task():
    assert to_binary(Immunophenotype.INFLAMED, LabelTask.PERI) is None
    assert to_binary(Immunophenotype.DESERT, LabelTask.PERI) == 0
    assert to_binary(Immunophenotype.EXCLUDED, LabelTask.PERI) == 1
    assert to_binary(Immunophenotype.INFLAMED, LabelTask.WHOLE) == 1


def test_peritumoural_task_requires_spatial_labels():
    assert LabelTask.PERI.requires_spatial_label
    assert not LabelTask.WHOLE.requires_spatial_label


def test_spatial_metrics_separate_excluded_from_inflamed():
    """A rim-only infiltrate must not be scored as inflamed."""
    coords = np.array([(x, y) for x in range(30) for y in range(30)], dtype=float) * 50.0
    x_um, y_um = coords[:, 0], coords[:, 1]
    centre = np.hypot(x_um - 725.0, y_um - 725.0)
    tumour = (centre < 600.0).astype(float)

    excluded = ((centre > 400.0) & (centre < 800.0)).astype(float)
    inflamed = tumour.copy()

    m_excluded = compute_spatial_metrics("EX", x_um, y_um, excluded, tumour, margin_width_um=200.0)
    m_inflamed = compute_spatial_metrics("IN", x_um, y_um, inflamed, tumour, margin_width_um=200.0)

    assert m_excluded.core_density < m_inflamed.core_density
    assert m_excluded.margin_density > m_excluded.core_density
    assert m_excluded.core_to_margin_ratio < m_inflamed.core_to_margin_ratio


def test_phenotype_thresholds_classify_all_three_classes():
    # Three cases per phenotype, with the fitted quantiles set to the group
    # proportions so the cut-points land between groups rather than inside one.
    groups = {
        Immunophenotype.DESERT: [(0.01, 0.02), (0.02, 0.03), (0.03, 0.04)],
        Immunophenotype.EXCLUDED: [(0.04, 0.60), (0.05, 0.65), (0.06, 0.70)],
        Immunophenotype.INFLAMED: [(0.50, 0.65), (0.55, 0.70), (0.60, 0.75)],
    }
    metrics, truth = [], []
    for phenotype, cases in groups.items():
        for i, (core, margin) in enumerate(cases):
            metrics.append(
                SpatialTILMetrics(
                    f"{phenotype.name}{i}", core_density=core, margin_density=margin,
                    stroma_density=0.0, n_core_patches=200, n_margin_patches=50,
                    n_tumour_patches=300,
                )
            )
            truth.append(phenotype)

    thresholds = PhenotypeThresholds.fit(metrics, core_quantile=2 / 3, margin_quantile=0.5)
    assigned = [thresholds.classify(m) for m in metrics]
    assert set(assigned) == set(Immunophenotype)
    # Every inflamed case must be recovered: a missed inflamed case would be
    # scored as desert, which is the most consequential possible label error.
    for phenotype, got in zip(truth, assigned):
        if phenotype is Immunophenotype.INFLAMED:
            assert got is Immunophenotype.INFLAMED


def test_signature_scoring_tracks_injected_immune_signal():
    genes = list(dict.fromkeys(sum((list(v) for v in SIGNATURES.values()), []))) + ["ACTB", "GAPDH"]
    rng = np.random.default_rng(11)
    expression = rng.normal(8.0, 1.0, (len(genes), 60))
    hot = np.zeros(60, dtype=bool)
    hot[:20] = True
    for gene in SIGNATURES["cd8_core"]:
        expression[genes.index(gene), hot] += 3.0

    scores = score_signatures(expression, genes, [f"s{i}" for i in range(60)])
    assert scores["cd8_core"][hot].mean() > scores["cd8_core"][~hot].mean() + 1.0


# ------------------------------------------------------------- simulator


def test_effect_size_round_trips_through_theoretical_auc():
    for target in (0.55, 0.62, 0.70, 0.80):
        assert theoretical_feature_auc(effect_size_for_auc(target)) == pytest.approx(target, abs=1e-9)


def test_simulated_informative_features_reach_their_target_auc():
    target = 0.70
    cohort = simulate_cohort(
        SimulationConfig(
            n_patients=4000, n_features=400, n_informative=8,
            effect_size=effect_size_for_auc(target), random_state=12,
        )
    )
    observed = columnwise_auc(cohort.X[:, cohort.informative_idx], cohort.y)
    assert observed.mean() == pytest.approx(target, abs=0.03)


def test_simulated_blocks_are_correlated():
    cohort = simulate_cohort(
        SimulationConfig(
            n_patients=800, n_features=200, block_size=40,
            within_block_correlation=0.75, n_informative=0, effect_size=0.0, random_state=13,
        )
    )
    corr = np.corrcoef(cohort.X, rowvar=False)
    within = corr[:40, :40][np.triu_indices(40, k=1)].mean()
    between = corr[:40, 40:80].mean()
    assert within == pytest.approx(0.75, abs=0.06)
    assert abs(between) < 0.1


def test_center_effects_are_detectable_before_harmonisation():
    cohort = simulate_cohort(
        SimulationConfig(
            n_patients=600, n_features=200, n_centers=5,
            center_shift_sd=0.8, center_scale_sd=0.3, random_state=14,
        )
    )
    assert measure_batch_effect(cohort.X, cohort.center).median_eta_squared > 0.05


# ------------------------------------------------------------ TIL geometry


def _phantom(kind: str, n: int = 240, radius: int = 100):
    """Idealised desert / excluded / inflamed maps on one tissue disc."""
    yy, xx = np.mgrid[:n, :n]
    distance = np.hypot(yy - n / 2, xx - n / 2)
    tissue = distance < radius
    rng = np.random.default_rng(0)
    if kind == "desert":
        til = tissue & (rng.random((n, n)) < 0.01)
    elif kind == "excluded":
        til = tissue & (distance > radius * 0.80) & (distance < radius * 0.98)
    else:
        til = tissue & (rng.random((n, n)) < 0.16)
    grid = np.where(tissue, 40, 0).astype(np.uint8)
    grid[til] = 200
    return grid


def test_penetration_separates_excluded_from_inflamed():
    """The descriptor that stands in for a tumour boundary must actually work."""
    from mri_cd8_til.labels.geometry import describe

    excluded = describe("excluded", _phantom("excluded"))
    inflamed = describe("inflamed", _phantom("inflamed"))

    # Excluded has MORE infiltrate here, so amount alone cannot separate them.
    assert excluded.til_fraction > inflamed.til_fraction
    # Penetration must, and by a wide margin.
    assert excluded.unreached_fraction > 0.3
    assert inflamed.unreached_fraction < 0.05
    assert excluded.penetration_p90 > 3 * inflamed.penetration_p90


def test_band_index_flags_banded_infiltrate():
    from mri_cd8_til.labels.geometry import describe

    assert describe("excluded", _phantom("excluded")).band_index > 5.0


def test_geometry_descriptors_are_scale_free():
    """Doubling the specimen must not change the normalised descriptors."""
    from mri_cd8_til.labels.geometry import describe

    small = describe("small", _phantom("excluded", n=240, radius=100))
    large = describe("large", _phantom("excluded", n=480, radius=200))
    assert small.til_fraction == pytest.approx(large.til_fraction, abs=0.03)
    assert small.unreached_fraction == pytest.approx(large.unreached_fraction, abs=0.08)
    assert small.penetration_p90 == pytest.approx(large.penetration_p90, abs=0.06)


def test_geometry_handles_a_map_with_no_infiltrate():
    from mri_cd8_til.labels.geometry import describe

    grid = np.zeros((120, 120), dtype=np.uint8)
    grid[20:100, 20:100] = 40  # tissue, no TIL anywhere
    geometry = describe("empty", grid)
    assert geometry.til_fraction == 0.0
    assert geometry.unreached_fraction == 1.0
    assert np.isfinite(geometry.penetration_p90)



# ------------------------------------------------------------------ I-SPY2


def test_multiframe_geometry_reads_orientation_from_shared_groups():
    """The bug that silently misregistered every I-SPY2 mask.

    PlaneOrientationSequence lives in SharedFunctionalGroupsSequence for these
    masks. Reading only the per-frame groups leaves SimpleITK's identity
    default, which flips x and y against the (-1, -1, 1) of the intensity
    volumes and moves the mask off the target grid entirely -- without raising.
    """
    pydicom = pytest.importorskip("pydicom")
    from mri_cd8_til.experiments.ispy2 import _image_from_multiframe

    n_frames, size = 4, 8

    measures = pydicom.Dataset()
    measures.PixelSpacing = [0.5, 0.5]
    measures.SliceThickness = 2.0
    orientation = pydicom.Dataset()
    orientation.ImageOrientationPatient = [-1, 0, 0, 0, -1, 0]
    shared_item = pydicom.Dataset()
    shared_item.PixelMeasuresSequence = pydicom.Sequence([measures])
    shared_item.PlaneOrientationSequence = pydicom.Sequence([orientation])

    # Frames stored back-to-front, to exercise the sort.
    per_frame = []
    for z in reversed(range(n_frames)):
        position = pydicom.Dataset()
        position.ImagePositionPatient = [0.0, 0.0, float(z) * 2.0]
        item = pydicom.Dataset()
        item.PlanePositionSequence = pydicom.Sequence([position])
        per_frame.append(item)

    class _MultiFrame:
        """Minimal stand-in: pydicom's pixel_array is read-only."""

        def __init__(self):
            self.pixel_array = np.arange(
                n_frames * size * size, dtype=np.uint8
            ).reshape(n_frames, size, size)
            self._groups = {
                "SharedFunctionalGroupsSequence": pydicom.Sequence([shared_item]),
                "PerFrameFunctionalGroupsSequence": pydicom.Sequence(per_frame),
            }

        def get(self, key, default=None):
            return self._groups.get(key, default)

    dataset = _MultiFrame()

    image = _image_from_multiframe(dataset)
    direction = np.array(image.GetDirection()).reshape(3, 3)
    # Columns are the axis cosines: x and y must be flipped, not identity.
    assert direction[0, 0] == pytest.approx(-1.0)
    assert direction[1, 1] == pytest.approx(-1.0)
    assert direction[2, 2] == pytest.approx(1.0)
    # Origin comes from the lowest frame after sorting, not the stored first.
    assert image.GetOrigin()[2] == pytest.approx(0.0)
    assert image.GetSpacing()[2] == pytest.approx(2.0)


def test_ispy2_analysis_mask_is_not_a_tumour_segmentation():
    """Pins the negative finding that blocked the I-SPY2 imaging cohort.

    The published mask is a breast analysis mask: every bit plane is either
    breast-scale and connected, or tiny and shattered into thousands of
    components. Neither is a tumour. This test encodes what a tumour-shaped
    plane would have to look like, so the claim in docs/ISPY2.md stays checkable.
    """
    from scipy.ndimage import label

    rng = np.random.default_rng(0)
    shape = (16, 64, 64)

    # bit 5: breast-scale, one component.
    breast = np.zeros(shape, dtype=bool)
    breast[2:14, 8:56, 8:56] = True
    # bit 1: scattered threshold flags.
    scattered = rng.random(shape) < 0.002

    def profile(mask):
        labelled, n = label(mask)
        sizes = np.bincount(labelled.ravel())[1:]
        return n, (sizes.max() / mask.sum() if mask.any() and len(sizes) else 0.0)

    breast_components, breast_share = profile(breast)
    scattered_components, scattered_share = profile(scattered)

    assert breast_components == 1 and breast_share == pytest.approx(1.0)
    assert scattered_components > 50 and scattered_share < 0.05
    # A tumour-shaped plane would be neither: few components, a dominant one,
    # and a small fraction of the breast. Nothing in the real mask matches.
    assert not (scattered_components < 10 and 0.5 < scattered_share < 1.0)
