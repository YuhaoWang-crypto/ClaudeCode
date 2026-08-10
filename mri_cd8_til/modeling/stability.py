"""Feature-stability filtering by intraclass correlation.

The reference study drops features with ICC < 0.8.  That is standard and
sensible, but the ICC there is computed between two human observers'
segmentations -- which measures only one of the several things that move a
radiomic feature.  A feature can be perfectly reproducible between two
radiologists at one centre and still be useless across centres because it moves
with the scanner.

This module therefore separates two different questions that both get called
"stability":

``ICC(2,1)`` **over perturbations**
    Re-extract features under small, label-preserving perturbations of the
    input -- mask erosion/dilation, in-plane rotation, added noise -- and keep
    features whose value survives.  Cheap: needs no second observer.

**Batch dominance**
    A separate filter, in ``modeling.combat.measure_batch_effect``: drop
    features whose variance is mostly explained by centre.  A feature can pass
    the ICC filter with flying colours and fail this one badly.

Both filters are unsupervised -- they never touch the outcome -- so applying
them to the full dataset does not leak labels.  That is worth stating precisely,
because it is the one preprocessing step in this pipeline that legitimately may
be done outside the cross-validation loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


def icc_2_1(measurements: np.ndarray) -> np.ndarray:
    """ICC(2,1) -- two-way random effects, absolute agreement, single measure.

    Parameters
    ----------
    measurements:
        ``(n_subjects, n_repeats, n_features)`` array of repeated measurements.

    Returns
    -------
    ``(n_features,)`` array of ICC values, clipped to ``[0, 1]``.
    """
    values = np.asarray(measurements, dtype=float)
    if values.ndim != 3:
        raise ValueError("measurements must be (n_subjects, n_repeats, n_features)")
    n, k, _ = values.shape
    if n < 2 or k < 2:
        raise ValueError("ICC needs at least 2 subjects and 2 repeats")

    subject_mean = values.mean(axis=1)
    repeat_mean = values.mean(axis=0)
    grand_mean = values.mean(axis=(0, 1))

    ms_rows = k * ((subject_mean - grand_mean) ** 2).sum(axis=0) / (n - 1)
    ms_cols = n * ((repeat_mean - grand_mean) ** 2).sum(axis=0) / (k - 1)
    residual = values - subject_mean[:, None, :] - repeat_mean[None, :, :] + grand_mean
    ms_error = (residual**2).sum(axis=(0, 1)) / ((n - 1) * (k - 1))

    denominator = ms_rows + (k - 1) * ms_error + k * (ms_cols - ms_error) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        icc = (ms_rows - ms_error) / denominator
    return np.clip(np.nan_to_num(icc, nan=0.0, posinf=0.0, neginf=0.0), 0.0, 1.0)


@dataclass
class PerturbationPlan:
    """Label-preserving perturbations used to probe feature stability.

    Chosen to imitate the things that actually differ between two acquisitions
    of the same tumour rather than to be maximally destructive.
    """

    mask_dilate_mm: tuple[float, ...] = (-1.0, 0.0, 1.0)
    rotation_degrees: tuple[float, ...] = (-5.0, 0.0, 5.0)
    noise_sigma_fraction: tuple[float, ...] = (0.0, 0.02)
    random_seed: int = 0

    @property
    def n_repeats(self) -> int:
        return (
            len(self.mask_dilate_mm) * len(self.rotation_degrees) * len(self.noise_sigma_fraction)
        )


class StabilitySelector(BaseEstimator, TransformerMixin):
    """Keep only features whose ICC over perturbations clears a threshold.

    ``icc_values`` are supplied by the caller because computing them requires
    re-running feature extraction, which is far too expensive to hide inside a
    ``fit``.  ``scripts/`` computes them once and caches them.
    """

    def __init__(self, icc_values: np.ndarray, threshold: float = 0.8) -> None:
        self.icc_values = np.asarray(icc_values, dtype=float)
        self.threshold = threshold

    def fit(self, X: np.ndarray, y: np.ndarray | None = None) -> "StabilitySelector":
        X = np.asarray(X)
        if self.icc_values.shape[0] != X.shape[1]:
            raise ValueError(
                f"icc_values has {self.icc_values.shape[0]} entries but X has {X.shape[1]} features"
            )
        self.support_ = self.icc_values >= self.threshold
        if not self.support_.any():
            raise ValueError(
                f"no feature reached ICC >= {self.threshold}; "
                "either the perturbations are too aggressive or the features are noise"
            )
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        return np.asarray(X)[:, self.support_]

    def get_support(self, indices: bool = False) -> np.ndarray:
        return np.where(self.support_)[0] if indices else self.support_


@dataclass
class StabilityReport:
    n_features: int
    n_stable: int
    threshold: float
    median_icc: float

    @property
    def fraction_stable(self) -> float:
        return self.n_stable / self.n_features if self.n_features else 0.0

    def describe(self) -> str:
        return (
            f"{self.n_stable}/{self.n_features} features with ICC >= {self.threshold:g} "
            f"({self.fraction_stable:.1%}); median ICC {self.median_icc:.3f}"
        )


def summarize(icc_values: np.ndarray, threshold: float = 0.8) -> StabilityReport:
    icc_values = np.asarray(icc_values, dtype=float)
    return StabilityReport(
        n_features=int(icc_values.size),
        n_stable=int((icc_values >= threshold).sum()),
        threshold=threshold,
        median_icc=float(np.median(icc_values)),
    )
