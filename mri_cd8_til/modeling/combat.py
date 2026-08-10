"""ComBat batch harmonisation, written to survive cross-validation honestly.

ComBat (Johnson, Li & Rabinovic, Biostatistics 2007) models each feature as

    y_ijg = alpha_g + X_ij*beta_g + gamma_ig + delta_ig * eps_ijg

for feature ``g`` in batch ``i``, and shrinks the batch effects ``gamma``,
``delta`` toward a common prior by empirical Bayes.  In radiomics the batch is
the scanner or the centre.

Two implementation choices here matter for multi-centre validation:

**Covariate preservation.**  If the immune phenotype is unevenly distributed
across centres -- and it always is -- then removing "centre" without protecting
the biological variable removes some of the signal with it.  Passing the outcome
as a covariate protects it, but the outcome is unavailable for the test fold, so
using it at fit time and not at transform time is a subtle form of leakage.
This implementation therefore defaults to covariate-free harmonisation and, if
covariates are supplied, requires them for transform as well, so the asymmetry
cannot happen silently.

**Unseen batches.**  Leave-one-centre-out validation guarantees the test batch
was never seen during fit.  Two defensible answers, both provided:

``reference``
    M-ComBat: map every batch onto a designated reference batch's location and
    scale.  An unseen batch is standardised using its *own* features only,
    which uses no labels and is therefore not leakage.
``error``
    Refuse.  Appropriate when the deployment scenario is a single new scan,
    where no batch-level statistics can be estimated at all -- which is the
    honest statement of the deployment limitation of ComBat-based radiomics.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


@dataclass
class _BatchParams:
    gamma_star: np.ndarray
    delta_star: np.ndarray
    n_samples: int = 0


class ComBat(BaseEstimator, TransformerMixin):
    """Empirical-Bayes ComBat as a scikit-learn transformer.

    Parameters
    ----------
    unseen_batch:
        ``"reference"`` standardises an unseen batch from its own statistics;
        ``"error"`` raises.
    parametric:
        Parametric empirical Bayes (inverse-gamma/normal priors).  The
        non-parametric variant is materially slower and rarely changes anything
        at radiomics feature counts, so it is not implemented.
    """

    def __init__(
        self,
        unseen_batch: str = "reference",
        parametric: bool = True,
        eps: float = 1e-8,
    ) -> None:
        self.unseen_batch = unseen_batch
        self.parametric = parametric
        self.eps = eps

    # --------------------------------------------------------------- sklearn

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
        batch: np.ndarray | None = None,
        covariates: np.ndarray | None = None,
    ) -> "ComBat":
        X = np.asarray(X, dtype=float)
        if batch is None:
            raise ValueError("ComBat.fit requires `batch`")
        batch = np.asarray(batch)
        if X.shape[0] != batch.shape[0]:
            raise ValueError("X and batch must have the same number of rows")

        self.batches_ = np.unique(batch)
        self.n_features_in_ = X.shape[1]
        self.uses_covariates_ = covariates is not None

        design = self._design(batch, covariates)
        # Grand mean and pooled variance, with batch effects projected out.
        coefficients, _, _, _ = np.linalg.lstsq(design, X, rcond=None)
        n_batch = len(self.batches_)
        batch_sizes = np.array([(batch == b).sum() for b in self.batches_], dtype=float)
        grand_mean = (batch_sizes / batch_sizes.sum()) @ coefficients[:n_batch]

        residual = X - design @ coefficients
        pooled_var = (residual**2).sum(axis=0) / X.shape[0]
        pooled_var = np.maximum(pooled_var, self.eps)

        self.grand_mean_ = grand_mean
        if self.uses_covariates_:
            self.covariate_coef_ = coefficients[n_batch:]
            covariate_effect = np.asarray(covariates, dtype=float) @ self.covariate_coef_
        else:
            self.covariate_coef_ = None
            covariate_effect = 0.0
        self.pooled_sd_ = np.sqrt(pooled_var)

        standardised = (X - grand_mean - covariate_effect) / self.pooled_sd_

        self.params_: dict[object, _BatchParams] = {}
        for b in self.batches_:
            rows = batch == b
            self.params_[b] = self._empirical_bayes(standardised[rows])
        return self

    def transform(
        self,
        X: np.ndarray,
        batch: np.ndarray | None = None,
        covariates: np.ndarray | None = None,
    ) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        if batch is None:
            raise ValueError("ComBat.transform requires `batch`")
        batch = np.asarray(batch)

        if self.uses_covariates_ and covariates is None:
            raise ValueError(
                "this ComBat was fitted with covariates, so transform needs them too; "
                "supplying them only at fit time would leak the outcome into training"
            )

        covariate_effect = (
            np.asarray(covariates, dtype=float) @ self.covariate_coef_
            if self.uses_covariates_
            else 0.0
        )
        standardised = (X - self.grand_mean_ - covariate_effect) / self.pooled_sd_

        out = np.empty_like(standardised)
        for b in np.unique(batch):
            rows = batch == b
            params = self.params_.get(b)
            if params is None:
                params = self._handle_unseen(standardised[rows], b)
            out[rows] = (standardised[rows] - params.gamma_star) / params.delta_star

        return out * self.pooled_sd_ + self.grand_mean_ + covariate_effect

    def fit_transform(  # type: ignore[override]
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
        batch: np.ndarray | None = None,
        covariates: np.ndarray | None = None,
    ) -> np.ndarray:
        return self.fit(X, y, batch=batch, covariates=covariates).transform(
            X, batch=batch, covariates=covariates
        )

    # --------------------------------------------------------------- internals

    def _design(self, batch: np.ndarray, covariates: np.ndarray | None) -> np.ndarray:
        indicators = np.column_stack([(batch == b).astype(float) for b in self.batches_])
        if covariates is None:
            return indicators
        return np.column_stack([indicators, np.asarray(covariates, dtype=float)])

    def _empirical_bayes(self, standardised_batch: np.ndarray) -> _BatchParams:
        n = standardised_batch.shape[0]
        gamma_hat = standardised_batch.mean(axis=0)
        delta_hat = np.maximum(standardised_batch.var(axis=0, ddof=1) if n > 1 else np.ones(
            standardised_batch.shape[1]
        ), self.eps)

        if not self.parametric or n < 3:
            # Too few samples for the priors to mean anything; use the raw
            # estimates rather than pretending to shrink.
            return _BatchParams(gamma_hat, np.sqrt(delta_hat), n)

        gamma_bar = gamma_hat.mean()
        tau2 = gamma_hat.var(ddof=1)
        a_prior, b_prior = _inverse_gamma_prior(delta_hat)

        gamma_star, delta_star = gamma_hat.copy(), delta_hat.copy()
        for _ in range(100):
            new_gamma = (n * tau2 * gamma_hat + delta_star * gamma_bar) / (n * tau2 + delta_star)
            sum_sq = ((standardised_batch - new_gamma) ** 2).sum(axis=0)
            new_delta = (0.5 * sum_sq + b_prior) / (n / 2.0 + a_prior - 1.0)
            new_delta = np.maximum(new_delta, self.eps)
            if np.max(np.abs(new_gamma - gamma_star)) < 1e-6 and np.max(
                np.abs(new_delta - delta_star)
            ) < 1e-6:
                gamma_star, delta_star = new_gamma, new_delta
                break
            gamma_star, delta_star = new_gamma, new_delta

        return _BatchParams(gamma_star, np.sqrt(delta_star), n)

    def _handle_unseen(self, standardised_batch: np.ndarray, batch_label: object) -> _BatchParams:
        if self.unseen_batch == "error":
            raise ValueError(
                f"batch {batch_label!r} was not seen during fit. Under "
                "unseen_batch='error' this is refused: no batch statistics can be "
                "estimated for a single unseen scan, which is a real deployment "
                "limitation of ComBat-harmonised radiomics, not a bug."
            )
        if self.unseen_batch != "reference":
            raise ValueError(f"unknown unseen_batch policy {self.unseen_batch!r}")
        if standardised_batch.shape[0] < 2:
            # One sample: centre it, but do not invent a scale.
            return _BatchParams(
                standardised_batch.mean(axis=0), np.ones(standardised_batch.shape[1]), 1
            )
        return self._empirical_bayes(standardised_batch)


def _inverse_gamma_prior(delta_hat: np.ndarray) -> tuple[float, float]:
    """Method-of-moments inverse-gamma hyperparameters, as in the original paper."""
    mean = float(delta_hat.mean())
    variance = float(delta_hat.var(ddof=1))
    if variance <= 0:
        return 2.0, mean
    a = (2.0 * variance + mean**2) / variance
    b = (mean * variance + mean**3) / variance
    return a, b


@dataclass
class BatchEffectReport:
    """How much of a feature set is explained by batch rather than biology."""

    n_features: int
    median_eta_squared: float
    fraction_features_batch_dominated: float
    per_feature_eta_squared: np.ndarray = field(repr=False, default_factory=lambda: np.array([]))


def measure_batch_effect(
    X: np.ndarray, batch: np.ndarray, dominance_threshold: float = 0.10
) -> BatchEffectReport:
    """One-way ANOVA effect size (eta^2) of batch on each feature.

    Run this before and after harmonisation.  If the median eta^2 before
    correction is large, any pooled-data model that skips harmonisation is
    partly a scanner classifier, and its cross-centre performance will not
    survive.
    """
    X = np.asarray(X, dtype=float)
    batch = np.asarray(batch)
    groups = np.unique(batch)
    grand_mean = X.mean(axis=0)

    ss_between = np.zeros(X.shape[1])
    for g in groups:
        rows = batch == g
        ss_between += rows.sum() * (X[rows].mean(axis=0) - grand_mean) ** 2
    ss_total = ((X - grand_mean) ** 2).sum(axis=0)
    ss_total = np.where(ss_total <= 0, np.nan, ss_total)

    eta2 = ss_between / ss_total
    return BatchEffectReport(
        n_features=X.shape[1],
        median_eta_squared=float(np.nanmedian(eta2)),
        fraction_features_batch_dominated=float(np.nanmean(eta2 > dominance_threshold)),
        per_feature_eta_squared=eta2,
    )
