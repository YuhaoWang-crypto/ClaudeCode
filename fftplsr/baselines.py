"""Baselines the FFT-PLSR model has to beat to be worth its complexity.

The paper reports FFT-PLSR's accuracy but not what a simpler model would have
scored on the same split. That comparison matters here more than usual: in a
combinatorial-recombination task only a handful of positions ever change, so the
protein spectrum is a deterministic function of a short binary mutation pattern.
Anything the 227-bin x 3-index spectrum can express is a function of those same
few bits -- so the honest question is whether the Fourier detour buys predictive
accuracy over regressing directly on the bits.

All models share one interface: ``fit(labels, y)`` / ``predict(labels)``, where
`labels` are variant strings like ``"D2N/H63Y"``.

* `MeanPredictor`      -- predict the training mean. Floor; R2 = 0 by construction.
* `LogAdditiveModel`   -- multiply the measured single-mutant fold-changes. The
                          textbook no-epistasis model of mutation stacking, and
                          the model a protein engineer uses without a computer.
* `OneHotRidge`        -- ridge on mutation indicators (additive, but *fitted*,
                          so it can use multi-mutant data and shrink noise).
* `OneHotPLS`          -- PLS on the same indicators: isolates how much of
                          FFT-PLSR's performance is PLS rather than the FFT.
* `PairwiseRidge`      -- indicators plus all pairwise products: the cheapest
                          model that can represent two-way epistasis at all.
"""

from __future__ import annotations

import itertools

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import RidgeCV

from .variants import parse_variant

__all__ = [
    "MeanPredictor",
    "LogAdditiveModel",
    "OneHotRidge",
    "OneHotPLS",
    "PairwiseRidge",
    "ALL_BASELINES",
]

ALPHAS = np.logspace(-3, 4, 40)


class MeanPredictor:
    name = "mean"

    def fit(self, labels, y):
        self._mu = float(np.mean(y))
        return self

    def predict(self, labels):
        return np.full(len(list(labels)), self._mu)


class LogAdditiveModel:
    """Fold-changes multiply: ``f(A/B) = f(A) * f(B)``, fitted only on singles.

    Fitness is a *relative* activity (parent = 1.0), so multiplying fold-changes
    -- adding in log space -- is the natural epistasis-free null. Singles absent
    from the training data are treated as neutral (fold-change 1.0) and counted
    in `n_missing`, so you can tell whether a prediction is informed or a guess.
    """

    name = "log-additive"

    def __init__(self, floor: float = 1e-3):
        self.floor = floor

    def fit(self, labels, y):
        self._effect: dict[str, float] = {}
        self._parent = 1.0
        for label, value in zip(labels, y):
            subs = parse_variant(label)
            if not subs:
                self._parent = float(value)
            elif len(subs) == 1:
                self._effect[str(subs[0])] = max(float(value), self.floor)
        self.n_missing = 0
        return self

    def predict(self, labels):
        out = []
        for label in labels:
            subs = parse_variant(label)
            value = self._parent
            for sub in subs:
                eff = self._effect.get(str(sub))
                if eff is None:
                    self.n_missing += 1
                    continue
                value *= eff / self._parent
            out.append(value)
        return np.asarray(out, dtype=float)


class _IndicatorModel:
    """Shared machinery: build a binary design matrix over the observed substitutions."""

    def _columns(self, labels):
        seen: list[str] = []
        for label in labels:
            for sub in parse_variant(label):
                if str(sub) not in seen:
                    seen.append(str(sub))
        return seen

    def _design(self, labels):
        rows = []
        for label in labels:
            present = {str(s) for s in parse_variant(label)}
            rows.append([1.0 if c in present else 0.0 for c in self._cols])
        return np.asarray(rows, dtype=float).reshape(len(rows), len(self._cols))

    def _expand(self, X):
        return X


class OneHotRidge(_IndicatorModel):
    name = "onehot-ridge"

    def fit(self, labels, y):
        labels = list(labels)
        self._cols = self._columns(labels)
        X = self._expand(self._design(labels))
        n_folds = min(5, max(2, len(y)))
        self._model = RidgeCV(alphas=ALPHAS, cv=n_folds).fit(X, np.asarray(y, dtype=float))
        return self

    def predict(self, labels):
        return self._model.predict(self._expand(self._design(list(labels))))


class PairwiseRidge(OneHotRidge):
    name = "pairwise-ridge"

    def _expand(self, X):
        n = X.shape[1]
        if n < 2:
            return X
        pairs = [X[:, i] * X[:, j] for i, j in itertools.combinations(range(n), 2)]
        return np.column_stack([X] + pairs)


class OneHotPLS(_IndicatorModel):
    name = "onehot-pls"

    def __init__(self, components=range(2, 11)):
        self.components = components

    def fit(self, labels, y):
        from .model import cv_score

        labels = list(labels)
        self._cols = self._columns(labels)
        X = self._design(labels)
        y = np.asarray(y, dtype=float)
        best = cv_score(X, y, cv=None, components=self.components)
        self.n_components = best.n_components
        self._model = PLSRegression(n_components=best.n_components, scale=True).fit(X, y)
        return self

    def predict(self, labels):
        return self._model.predict(self._design(list(labels))).ravel()


ALL_BASELINES = [MeanPredictor, LogAdditiveModel, OneHotRidge, OneHotPLS, PairwiseRidge]
