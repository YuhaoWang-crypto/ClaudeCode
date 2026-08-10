"""Validation schemes, ordered by how hard they are to fool.

============================  =========================================
Scheme                        What it estimates
============================  =========================================
single random split           performance on *this* split; nothing more
repeated stratified CV        internal performance, averaged
nested CV                     internal performance with tuning accounted
leave-one-centre-out          performance at a centre you have never seen
permutation null              what the pipeline scores on shuffled labels
============================  =========================================

The reference study reports the first.  For a model intended to travel between
hospitals, leave-one-centre-out is the estimate that answers the actual
question, and it is routinely 0.10-0.20 AUC below the internal estimate --
which is a property of the problem, not a failure of the model.

The permutation null deserves particular attention with p >> n.  If a pipeline
scores 0.75 on shuffled labels, then its 0.85 on real labels means considerably
less than it appears to, and the difference between the two is the only part
that is real.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.base import BaseEstimator, clone
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score


@dataclass
class ValidationResult:
    scheme: str
    scores: np.ndarray = field(repr=False)
    fold_labels: tuple[str, ...] = ()

    @property
    def mean(self) -> float:
        return float(np.nanmean(self.scores))

    @property
    def std(self) -> float:
        return float(np.nanstd(self.scores, ddof=1)) if self.scores.size > 1 else 0.0

    def percentile_interval(self, alpha: float = 0.05) -> tuple[float, float]:
        return (
            float(np.nanquantile(self.scores, alpha / 2)),
            float(np.nanquantile(self.scores, 1 - alpha / 2)),
        )

    def describe(self) -> str:
        low, high = self.percentile_interval()
        return (
            f"{self.scheme}: AUC {self.mean:.3f} +/- {self.std:.3f} "
            f"[{low:.3f}, {high:.3f}] over {self.scores.size} folds"
        )


def repeated_cv(
    estimator: BaseEstimator,
    X: np.ndarray,
    y: np.ndarray,
    n_splits: int = 5,
    n_repeats: int = 20,
    random_state: int = 0,
) -> ValidationResult:
    """Repeated stratified k-fold, scoring AUC.

    Repeats matter: a single 5-fold CV on n=182 has enough split-to-split
    variance that two runs with different seeds can differ by 0.05 AUC. Reporting
    one run's number invites reporting the best one.
    """
    X, y = np.asarray(X, dtype=float), np.asarray(y).astype(int)
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    scores = [
        _fit_score(estimator, X, y, train, test) for train, test in cv.split(X, y)
    ]
    return ValidationResult(scheme=f"repeated {n_splits}-fold CV x{n_repeats}", scores=np.array(scores))


def nested_cv(
    estimator: BaseEstimator,
    X: np.ndarray,
    y: np.ndarray,
    param_grid: dict[str, list],
    outer_splits: int = 5,
    inner_splits: int = 5,
    random_state: int = 0,
) -> ValidationResult:
    """Nested CV: hyperparameters chosen inside each outer training fold.

    Choosing the regularisation strength on the same data used to report
    performance is a smaller leak than selecting features that way, but it is
    the same kind of leak and it compounds with it.
    """
    from sklearn.model_selection import GridSearchCV

    X, y = np.asarray(X, dtype=float), np.asarray(y).astype(int)
    outer = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=random_state)
    scores = []
    for train, test in outer.split(X, y):
        inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=random_state)
        search = GridSearchCV(clone(estimator), param_grid, scoring="roc_auc", cv=inner, n_jobs=1)
        search.fit(X[train], y[train])
        scores.append(_score(search.best_estimator_, X[test], y[test]))
    return ValidationResult(scheme=f"nested CV ({outer_splits}x{inner_splits})", scores=np.array(scores))


def leave_one_center_out(
    estimator: BaseEstimator,
    X: np.ndarray,
    y: np.ndarray,
    center: np.ndarray,
    min_test_size: int = 10,
) -> ValidationResult:
    """Hold out one centre at a time -- the estimate that answers the question.

    Centres with too few patients, or with only one class present, are skipped
    and reported rather than silently contributing a degenerate AUC.
    """
    X, y = np.asarray(X, dtype=float), np.asarray(y).astype(int)
    center = np.asarray(center)

    scores, labels = [], []
    for held_out in np.unique(center):
        test = center == held_out
        train = ~test
        if test.sum() < min_test_size or len(np.unique(y[test])) < 2 or len(np.unique(y[train])) < 2:
            continue
        scores.append(_fit_score(estimator, X, y, np.where(train)[0], np.where(test)[0]))
        labels.append(f"{held_out} (n={int(test.sum())})")

    if not scores:
        raise ValueError(
            "no centre had enough patients of both classes to be held out; "
            "leave-one-centre-out is not estimable on this cohort"
        )
    return ValidationResult(
        scheme="leave-one-centre-out", scores=np.array(scores), fold_labels=tuple(labels)
    )


def permutation_null(
    estimator: BaseEstimator,
    X: np.ndarray,
    y: np.ndarray,
    n_permutations: int = 200,
    n_splits: int = 5,
    random_state: int = 0,
) -> ValidationResult:
    """Distribution of CV AUC when the labels carry no information.

    Anything the pipeline scores here it can score on noise. The honest measure
    of a model is its real score minus this distribution, not its real score.
    """
    X, y = np.asarray(X, dtype=float), np.asarray(y).astype(int)
    rng = np.random.default_rng(random_state)
    scores = np.empty(n_permutations)
    for i in range(n_permutations):
        permuted = rng.permutation(y)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=int(rng.integers(1 << 30)))
        fold_scores = [
            _fit_score(estimator, X, permuted, train, test) for train, test in cv.split(X, permuted)
        ]
        scores[i] = float(np.nanmean(fold_scores))
    return ValidationResult(scheme=f"permutation null ({n_permutations} shuffles)", scores=scores)


def permutation_p_value(observed: float, null_result: ValidationResult) -> float:
    """One-sided permutation p-value with the standard +1 correction."""
    null = null_result.scores
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def _fit_score(
    estimator: BaseEstimator, X: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray
) -> float:
    model = clone(estimator)
    try:
        model.fit(X[train], y[train])
    except ValueError:  # pragma: no cover - degenerate fold
        return float("nan")
    return _score(model, X[test], y[test])


def _score(model: BaseEstimator, X: np.ndarray, y: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    if hasattr(model, "decision_function"):
        scores = model.decision_function(X)
    else:
        scores = model.predict_proba(X)[:, 1]
    return float(roc_auc_score(y, scores))
