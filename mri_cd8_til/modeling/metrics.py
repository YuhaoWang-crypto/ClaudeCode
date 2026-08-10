"""Discrimination, calibration and clinical-utility metrics -- with intervals.

An AUC reported without an interval on a validation set of 45 patients is not a
measurement, it is a point.  The single most useful thing this module does is
:func:`auc_confidence_interval`: on the reference study's validation split, an
AUC of 0.985 carries a confidence interval wide enough to be equally consistent
with a considerably more ordinary model, and that fact is visible from the
published numbers alone, without re-analysing any data.

Calibration and net benefit are here because a model intended to select patients
for immunotherapy is a decision rule, and a decision rule with excellent ranking
and terrible calibration will still make bad decisions at any fixed threshold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from sklearn.metrics import roc_auc_score


@dataclass(frozen=True)
class AUCInterval:
    auc: float
    lower: float
    upper: float
    n_positive: int
    n_negative: int
    method: str

    @property
    def width(self) -> float:
        return self.upper - self.lower

    def describe(self) -> str:
        return (
            f"AUC {self.auc:.3f} (95% CI {self.lower:.3f}-{self.upper:.3f}, "
            f"width {self.width:.3f}; n+={self.n_positive}, n-={self.n_negative}, {self.method})"
        )


def auc_confidence_interval(
    y_true: np.ndarray, y_score: np.ndarray, alpha: float = 0.05
) -> AUCInterval:
    """Hanley-McNeil confidence interval for the AUC.

    Uses the standard exponential approximation for the variance:

        Q1 = A / (2 - A),   Q2 = 2A^2 / (1 + A)

    which needs only the AUC and the two class counts.  That makes it applicable
    to *published* results, not just to data in hand -- so a reader can put an
    interval on someone else's headline number.
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n_pos = int((y_true == 1).sum())
    n_neg = int((y_true == 0).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUC needs both classes present")
    auc = float(roc_auc_score(y_true, y_score))
    return _hanley_mcneil(auc, n_pos, n_neg, alpha)


def auc_interval_from_reported(
    auc: float, n_positive: int, n_negative: int, alpha: float = 0.05
) -> AUCInterval:
    """Interval for an AUC taken from a paper, given only the class counts."""
    return _hanley_mcneil(auc, n_positive, n_negative, alpha)


def _hanley_mcneil(auc: float, n_pos: int, n_neg: int, alpha: float) -> AUCInterval:
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc**2 / (1.0 + auc)
    variance = (
        auc * (1 - auc) + (n_pos - 1) * (q1 - auc**2) + (n_neg - 1) * (q2 - auc**2)
    ) / (n_pos * n_neg)
    se = math.sqrt(max(variance, 0.0))
    z = _normal_quantile(1 - alpha / 2)
    return AUCInterval(
        auc=auc,
        lower=max(0.0, auc - z * se),
        upper=min(1.0, auc + z * se),
        n_positive=n_pos,
        n_negative=n_neg,
        method="Hanley-McNeil",
    )


def _normal_quantile(p: float) -> float:
    """Inverse standard normal CDF via the error function."""
    return math.sqrt(2.0) * _erfinv(2.0 * p - 1.0)


def _erfinv(x: float) -> float:
    # Winitzki's approximation, refined by two Newton steps: accurate to ~1e-12
    # over the range we use, and avoids a SciPy import in a hot path.
    a = 0.147
    ln1mx2 = math.log(1 - x * x)
    term = 2 / (math.pi * a) + ln1mx2 / 2
    y = math.copysign(math.sqrt(math.sqrt(term * term - ln1mx2 / a) - term), x)
    for _ in range(2):
        err = math.erf(y) - x
        y -= err / (2 / math.sqrt(math.pi) * math.exp(-y * y))
    return y


@dataclass(frozen=True)
class CalibrationReport:
    brier: float
    intercept: float
    slope: float
    n: int

    @property
    def is_well_calibrated(self) -> bool:
        """Slope near 1 and intercept near 0 (the usual rough reading)."""
        return abs(self.slope - 1.0) < 0.2 and abs(self.intercept) < 0.2

    def describe(self) -> str:
        verdict = "ok" if self.is_well_calibrated else "MISCALIBRATED"
        return (
            f"Brier {self.brier:.4f}; calibration intercept {self.intercept:+.3f}, "
            f"slope {self.slope:.3f} [{verdict}]"
        )


def calibration(y_true: np.ndarray, y_prob: np.ndarray) -> CalibrationReport:
    """Brier score plus the calibration intercept and slope.

    Slope is the coefficient of the linear predictor in a logistic regression of
    the outcome on ``logit(p)``.  Slope < 1 is the classic signature of an
    overfitted model: predictions are too extreme in both directions.
    """
    y_true = np.asarray(y_true).astype(float)
    p = np.clip(np.asarray(y_prob, dtype=float), 1e-6, 1 - 1e-6)
    logit = np.log(p / (1 - p))

    slope, intercept = _logistic_fit(logit, y_true)
    return CalibrationReport(
        brier=float(np.mean((p - y_true) ** 2)),
        intercept=float(intercept),
        slope=float(slope),
        n=int(y_true.size),
    )


def _logistic_fit(x: np.ndarray, y: np.ndarray, iterations: int = 100) -> tuple[float, float]:
    """Two-parameter logistic regression by Newton-Raphson (slope, intercept)."""
    design = np.column_stack([x, np.ones_like(x)])
    beta = np.zeros(2)
    for _ in range(iterations):
        eta = design @ beta
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-9, None)
        gradient = design.T @ (y - mu)
        hessian = design.T @ (design * w[:, None])
        try:
            step = np.linalg.solve(hessian + 1e-9 * np.eye(2), gradient)
        except np.linalg.LinAlgError:  # pragma: no cover
            break
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    return float(beta[0]), float(beta[1])


def net_benefit(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> float:
    """Net benefit at one threshold probability (decision-curve analysis).

    ``NB = TP/n - FP/n * pt/(1-pt)``.  Compare against treat-all and treat-none:
    a model whose net benefit never exceeds both is clinically useless however
    good its AUC.
    """
    y_true = np.asarray(y_true).astype(int)
    p = np.asarray(y_prob, dtype=float)
    n = y_true.size
    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must be in (0, 1)")
    predicted = p >= threshold
    tp = int(((predicted == 1) & (y_true == 1)).sum())
    fp = int(((predicted == 1) & (y_true == 0)).sum())
    return tp / n - (fp / n) * (threshold / (1 - threshold))


def decision_curve(
    y_true: np.ndarray, y_prob: np.ndarray, thresholds: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Net benefit of the model, treat-all and treat-none across thresholds."""
    y_true = np.asarray(y_true).astype(int)
    thresholds = np.asarray(
        thresholds if thresholds is not None else np.linspace(0.01, 0.99, 99), dtype=float
    )
    prevalence = float(y_true.mean())
    return {
        "threshold": thresholds,
        "model": np.array([net_benefit(y_true, y_prob, t) for t in thresholds]),
        "treat_all": np.array([prevalence - (1 - prevalence) * t / (1 - t) for t in thresholds]),
        "treat_none": np.zeros_like(thresholds),
    }


def bootstrap_auc_ci(
    y_true: np.ndarray,
    y_score: np.ndarray,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 0,
) -> AUCInterval:
    """Stratified bootstrap percentile interval for the AUC.

    Slower than Hanley-McNeil but makes no distributional assumption; use it
    when reporting your own results and the analytic interval when putting an
    interval on someone else's.
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    pos = np.where(y_true == 1)[0]
    neg = np.where(y_true == 0)[0]
    if pos.size == 0 or neg.size == 0:
        raise ValueError("AUC needs both classes present")

    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate(
            [rng.choice(pos, pos.size, replace=True), rng.choice(neg, neg.size, replace=True)]
        )
        draws[i] = roc_auc_score(y_true[idx], y_score[idx])

    return AUCInterval(
        auc=float(roc_auc_score(y_true, y_score)),
        lower=float(np.quantile(draws, alpha / 2)),
        upper=float(np.quantile(draws, 1 - alpha / 2)),
        n_positive=int(pos.size),
        n_negative=int(neg.size),
        method=f"stratified bootstrap ({n_boot} resamples)",
    )
