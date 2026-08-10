"""Model pipelines: the reference protocol, and a leakage-free counterpart.

Both are built here side by side on purpose.  The interesting quantity is not
how well the leakage-free pipeline performs -- it is the *gap* between the two
on identical data, which is exactly the amount of performance the protocol
manufactures.  Keeping both implementations in one module makes that comparison
a two-line experiment instead of an argument.

The reference protocol, as published:

1. drop features with ICC < 0.8;
2. keep the **50 features with the highest AUC**;
3. fit LASSO with 5-fold CV on the training split;
4. report AUC on the held-out validation split.

Step 2 is where it goes wrong, and only if it is run on all 182 patients rather
than on the 137 training patients.  The paper does not state which it did.  The
distinction is invisible in the write-up and worth roughly the entire headline
result -- see ``mri_cd8_til.experiments.optimism``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class UnivariateAUCSelector(BaseEstimator, TransformerMixin):
    """Keep the ``k`` features with the highest univariate AUC.

    This is the reference study's step 2.  It is a *supervised* selector: it
    reads ``y``.  Placed inside a scikit-learn ``Pipeline`` and driven by
    cross-validation it is perfectly legitimate, because each fold re-selects
    using only that fold's training labels.  Run once on the pooled data before
    splitting, it is not.
    """

    def __init__(self, k: int = 50) -> None:
        self.k = k

    def fit(self, X: np.ndarray, y: np.ndarray) -> "UnivariateAUCSelector":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y).astype(int)
        if len(np.unique(y)) < 2:
            raise ValueError("univariate AUC selection needs both classes")

        scores = columnwise_auc(X, y)
        # Distance from 0.5: a feature that is strongly anti-correlated is just
        # as informative as one strongly correlated, and the linear model can
        # flip its sign.
        strength = np.abs(scores - 0.5)
        k = min(self.k, X.shape[1])
        self.support_ = np.zeros(X.shape[1], dtype=bool)
        self.support_[np.argsort(strength)[::-1][:k]] = True
        self.scores_ = scores
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        return np.asarray(X)[:, self.support_]

    def get_support(self, indices: bool = False) -> np.ndarray:
        return np.where(self.support_)[0] if indices else self.support_


def columnwise_auc(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """AUC of every column of ``X`` against ``y``, computed in one pass.

    Uses the Mann-Whitney identity ``AUC = (R+ - n+(n+ + 1)/2) / (n+ n-)`` on
    mid-ranks, which handles ties correctly and matches ``roc_auc_score``.

    This is not a micro-optimisation: the selection step is run thousands of
    times across the optimism experiment, and a per-column Python loop over
    3,332 features would dominate the entire runtime.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(int)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return np.full(X.shape[1], 0.5)

    ranks = _midrank(X)
    positive_rank_sum = ranks[y == 1].sum(axis=0)
    auc = (positive_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)

    constant = ~np.isfinite(X).all(axis=0) | (X.max(axis=0) == X.min(axis=0))
    auc[constant] = 0.5
    return auc


def _midrank(X: np.ndarray) -> np.ndarray:
    """Column-wise ranks starting at 1, averaging ties."""
    n = X.shape[0]
    order = np.argsort(X, axis=0, kind="mergesort")
    ranks = np.empty_like(X, dtype=float)
    positions = np.arange(1, n + 1, dtype=float)[:, None]
    np.put_along_axis(ranks, order, np.broadcast_to(positions, X.shape).copy(), axis=0)

    sorted_X = np.take_along_axis(X, order, axis=0)
    sorted_ranks = np.take_along_axis(ranks, order, axis=0)
    for j in range(X.shape[1]):
        column = sorted_X[:, j]
        i = 0
        while i < n:
            k = i + 1
            while k < n and column[k] == column[i]:
                k += 1
            if k - i > 1:
                sorted_ranks[i:k, j] = sorted_ranks[i:k, j].mean()
            i = k
    np.put_along_axis(ranks, order, sorted_ranks, axis=0)
    return ranks


def _safe_auc(y: np.ndarray, x: np.ndarray) -> float:
    if not np.isfinite(x).all() or np.allclose(x, x[0]):
        return 0.5
    try:
        return float(roc_auc_score(y, x))
    except ValueError:  # pragma: no cover
        return 0.5


def make_lasso_logistic(
    C: float = 1.0, max_iter: int = 10_000, random_state: int = 0
) -> LogisticRegression:
    """L1-penalised logistic regression, across scikit-learn API generations.

    scikit-learn 1.8 deprecated ``penalty='l1'`` in favour of ``l1_ratio=1``.
    Both spellings are accepted by different versions in the wild, so the
    supported one is chosen by inspecting the signature rather than by pinning a
    version -- a pin here would silently break every downstream reproduction.
    """
    import inspect

    kwargs: dict[str, object] = {
        "C": C,
        "solver": "liblinear",
        "max_iter": max_iter,
        "random_state": random_state,
    }
    parameters = inspect.signature(LogisticRegression.__init__).parameters
    if "l1_ratio" in parameters:
        kwargs["l1_ratio"] = 1.0
    else:  # pragma: no cover - scikit-learn < 1.8
        kwargs["penalty"] = "l1"
    return LogisticRegression(**kwargs)  # type: ignore[arg-type]


@dataclass
class PipelineConfig:
    """Knobs shared by both pipelines, so the comparison stays fair."""

    n_select: int = 50
    lasso_cs: int = 20
    lasso_cv: int = 5
    max_iter: int = 10_000
    random_state: int = 0
    #: ``None`` skips variance filtering.
    variance_threshold: float | None = 1e-10


def build_pipeline(config: PipelineConfig | None = None) -> Pipeline:
    """The leakage-free pipeline: every supervised step lives inside the fold.

    Standardisation, univariate selection and the LASSO fit are all pipeline
    steps, so a ``cross_val_score`` over this object re-runs all three on each
    training fold and never lets the test fold influence any of them.
    """
    config = config or PipelineConfig()
    return Pipeline(
        [
            ("scale", StandardScaler()),
            ("select", UnivariateAUCSelector(k=config.n_select)),
            ("lasso", make_lasso_logistic(1.0, config.max_iter, config.random_state)),
        ]
    )


@dataclass
class LeakyProtocolResult:
    """Outcome of running the reference protocol as published."""

    train_auc: float
    validation_auc: float
    n_selected: int
    selected_indices: np.ndarray

    def describe(self) -> str:
        return (
            f"train AUC {self.train_auc:.3f} / validation AUC {self.validation_auc:.3f} "
            f"({self.n_selected} features selected)"
        )


def run_reference_protocol(
    X: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    valid_idx: np.ndarray,
    config: PipelineConfig | None = None,
    select_on_all_data: bool = True,
) -> LeakyProtocolResult:
    """Run the published protocol, with the leaky step switchable.

    ``select_on_all_data=True`` reproduces the protocol as a reader would most
    naturally implement it from the paper: rank features by AUC over the whole
    cohort, then split.  ``False`` restricts the ranking to the training split,
    which is the only defensible version.  Everything else is held identical, so
    the difference between the two runs is attributable to that one choice.
    """
    config = config or PipelineConfig()
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(int)

    scaler = StandardScaler().fit(X[train_idx])
    Xs = scaler.transform(X)

    selector = UnivariateAUCSelector(k=config.n_select)
    if select_on_all_data:
        selector.fit(Xs, y)  # <-- sees the validation labels
    else:
        selector.fit(Xs[train_idx], y[train_idx])

    Xsel = selector.transform(Xs)
    model = make_lasso_logistic(1.0, config.max_iter, config.random_state).fit(
        Xsel[train_idx], y[train_idx]
    )

    return LeakyProtocolResult(
        train_auc=_safe_auc(y[train_idx], model.decision_function(Xsel[train_idx])),
        validation_auc=_safe_auc(y[valid_idx], model.decision_function(Xsel[valid_idx])),
        n_selected=int(selector.support_.sum()),
        selected_indices=selector.get_support(indices=True),
    )
