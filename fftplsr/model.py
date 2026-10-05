"""PLS regression on protein spectra, with greedy AAindex selection and honest scoring.

Three deliberate departures from the reference FP4COM implementation:

1. **Scores are cross-validated, not in-sample.** The original `regscore` fit on
   all of X and then scored `predict(X)` against `y`, so its reported MSE/R2 were
   training-set numbers. Here `R2` and `MSE` come from out-of-fold predictions;
   the in-sample values are still available but labelled `train_r2`/`train_mse`.

2. **Selection bias is measurable.** Choosing the best of 566 AAindex entries *by*
   cross-validated error makes that error optimistic. `nested_cv_r2` repeats the
   whole selection inside an outer loop to put a number on the gap.

3. **n_components is capped** at what each fold can support, instead of relying on
   a ValueError fallback that silently scored a candidate as -100.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import KFold, LeaveOneOut

from .encode import available_indices

__all__ = [
    "CVResult",
    "cv_score",
    "screen_indices",
    "IndexSelection",
    "fit_predict",
    "nested_cv_r2",
]

DEFAULT_COMPONENTS = range(2, 11)


def _splitter(cv, n_samples):
    if cv is None or (isinstance(cv, int) and cv >= n_samples):
        return LeaveOneOut()
    if isinstance(cv, int):
        return KFold(n_splits=cv, shuffle=True, random_state=0)
    return cv


def r2(y_true, y_pred) -> float:
    """Coefficient of determination, 1 - SS_res/SS_tot (the paper's definition)."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_tot = float(((y_true - y_true.mean()) ** 2).sum())
    if ss_tot == 0.0:
        return float("nan")
    return 1.0 - float(((y_true - y_pred) ** 2).sum()) / ss_tot


@dataclass
class CVResult:
    """Out-of-fold performance of one (feature matrix, n_components) choice."""

    n_components: int
    cv_mse: float
    cv_r2: float
    train_mse: float
    train_r2: float
    oof: np.ndarray = field(repr=False)

    def as_row(self) -> dict:
        return {
            "n_components": self.n_components,
            "cvMSE": self.cv_mse,
            "cvR2": self.cv_r2,
            "trainMSE": self.train_mse,
            "trainR2": self.train_r2,
        }


def pls_path_predict(X_train, y_train, X_test, max_components: int) -> np.ndarray:
    """Predictions for **every** component count 1..`max_components` from a single fit.

    PLS builds its components sequentially, so a model fitted with K components
    already contains the 1..K-component models: the k-component regression matrix
    is ``x_rotations_[:, :k] @ y_loadings_[:, :k].T``. Slicing that instead of
    refitting per k makes a 566-descriptor screen ~9x cheaper, which is what turns
    the exhaustive screen from minutes into seconds.

    Standardization is done here rather than via ``scale=True`` so that the
    rotations are expressed in a fixed space; this reproduces ``scale=True``
    exactly (`test_model.py` asserts it against per-k refits).

    Returns an array of shape ``(max_components, n_test)``.
    """
    X_train = np.asarray(X_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float).ravel()
    X_test = np.asarray(X_test, dtype=float)

    x_mean = X_train.mean(axis=0)
    x_std = X_train.std(axis=0, ddof=1)
    x_std[x_std == 0.0] = 1.0
    y_mean = y_train.mean()
    y_std = y_train.std(ddof=1)
    if y_std == 0.0:
        y_std = 1.0

    Xs = (X_train - x_mean) / x_std
    ys = (y_train - y_mean) / y_std
    pls = PLSRegression(n_components=max_components, scale=False).fit(Xs, ys)

    Zs = (X_test - x_mean) / x_std
    out = np.empty((max_components, len(X_test)))
    for k in range(1, max_components + 1):
        beta = pls.x_rotations_[:, :k] @ pls.y_loadings_[:, :k].T
        out[k - 1] = (Zs @ beta).ravel() * y_std + y_mean
    return out


def cv_score(X, y, cv=None, components=DEFAULT_COMPONENTS) -> CVResult:
    """Pick `n_components` by out-of-fold MSE and report both OOF and in-sample scores."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    splitter = _splitter(cv, len(y))
    splits = list(splitter.split(X))
    # A PLS component needs at least one sample and one feature to spare.
    max_comp = min(X.shape[1], min(len(tr) for tr, _ in splits) - 1)
    usable = [c for c in components if 1 <= c <= max_comp] or [1]
    top = max(usable)

    oof = np.empty((top, len(y)))
    for train, test in splits:
        oof[:, test] = pls_path_predict(X[train], y[train], X[test], top)
    mse = ((y[None, :] - oof) ** 2).mean(axis=1)

    k = min(usable, key=lambda c: mse[c - 1])
    in_sample = pls_path_predict(X, y, X, k)[k - 1]
    return CVResult(
        n_components=k,
        cv_mse=float(mse[k - 1]),
        cv_r2=r2(y, oof[k - 1]),
        train_mse=float(((y - in_sample) ** 2).mean()),
        train_r2=r2(y, in_sample),
        oof=oof[k - 1],
    )


@dataclass
class IndexSelection:
    """Outcome of greedy forward AAindex selection."""

    indices: list[str]
    n_components: int
    cv_mse: float
    cv_r2: float
    rounds: list[pd.DataFrame] = field(default_factory=list, repr=False)

    @property
    def label(self) -> str:
        return "_".join(self.indices)


def screen_indices(
    encode,
    y,
    n_rounds: int = 1,
    cv=None,
    candidates=None,
    components=DEFAULT_COMPONENTS,
    n_jobs: int = -1,
    verbose: bool = True,
) -> IndexSelection:
    """Greedy forward selection of AAindex entries, scored by out-of-fold MSE.

    `encode` is a callable taking a list of accessions and returning the feature
    matrix. This mirrors the paper: screen one entry, fix it, screen the next
    against it, and so on for `n_rounds`.
    """
    from joblib import Parallel, delayed

    y = np.asarray(y, dtype=float).ravel()
    pool = list(candidates) if candidates is not None else list(available_indices())
    chosen: list[str] = []
    rounds: list[pd.DataFrame] = []
    best: CVResult | None = None

    for rnd in range(n_rounds):
        remaining = [c for c in pool if c not in chosen]
        if not remaining:
            break

        def score_one(code):
            try:
                return code, cv_score(encode(chosen + [code]), y, cv=cv, components=components)
            except Exception as exc:
                return code, exc

        # Threads, not processes: the per-candidate work is numpy/BLAS (which releases
        # the GIL) and the shared encoder caches the parent spectrum, so re-pickling
        # it into 566 worker processes would cost more than the fits themselves.
        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(score_one)(code) for code in remaining
        )
        scored = {code: res for code, res in results if isinstance(res, CVResult)}
        failed = {code: res for code, res in results if not isinstance(res, CVResult)}
        rows = {"_".join(chosen + [code]): res.as_row() for code, res in scored.items()}
        if not rows:
            first = next(iter(failed.values()), None)
            raise RuntimeError(
                f"every one of {len(remaining)} candidate AAindex entries failed in round "
                f"{rnd + 1}; first error: {type(first).__name__}: {first}"
            )
        if failed and verbose:
            # Don't let a candidate vanish from the pool without saying so.
            kinds = sorted({type(e).__name__ for e in failed.values()})
            print(
                f"  round {rnd + 1}: skipped {len(failed)} of {len(remaining)} candidates "
                f"({', '.join(kinds)}); {len(rows)} scored"
            )
        table = pd.DataFrame(rows).T.sort_values("cvMSE")
        table["n_components"] = table["n_components"].astype(int)
        rounds.append(table)

        winner = dict(results)[table.index[0].split("_")[-1]]
        if best is not None and winner.cv_mse >= best.cv_mse:
            if verbose:
                print(
                    f"  round {rnd + 1}: no candidate improved cvMSE "
                    f"({winner.cv_mse:.4f} >= {best.cv_mse:.4f}); stopping with {chosen}"
                )
            break
        chosen.append(table.index[0].split("_")[-1])
        best = winner
        if verbose:
            print(
                f"  round {rnd + 1}: {'_'.join(chosen)}  "
                f"cvMSE={best.cv_mse:.4f} cvR2={best.cv_r2:.3f} "
                f"n_components={best.n_components}"
            )

    if best is None:
        raise RuntimeError("index selection produced no usable model")
    return IndexSelection(
        indices=chosen,
        n_components=best.n_components,
        cv_mse=best.cv_mse,
        cv_r2=best.cv_r2,
        rounds=rounds,
    )


def fit_predict(X_train, y_train, X_predict, n_components: int) -> np.ndarray:
    """Fit PLSR on the training block and predict the (usually much larger) query block."""
    pls = PLSRegression(n_components=n_components, scale=True)
    pls.fit(np.asarray(X_train, dtype=float), np.asarray(y_train, dtype=float).ravel())
    return pls.predict(np.asarray(X_predict, dtype=float)).ravel()


def nested_cv_r2(encode, y, n_rounds: int = 1, inner_cv=None, outer_cv=None,
                 candidates=None, components=DEFAULT_COMPONENTS, n_jobs: int = -1):
    """Honest generalisation estimate: AAindex selection repeated inside an outer loop.

    Returns ``(r2, mse, oof, picked)`` where `picked` lists the accessions chosen
    in each outer fold -- how stable that list is tells you whether the winning
    descriptor is a real signal or a coin flip among 566 near-ties.
    """
    y = np.asarray(y, dtype=float).ravel()
    splits = list(_splitter(outer_cv, len(y)).split(np.zeros((len(y), 1))))
    oof = np.empty(len(y))
    picked = []
    for i, (train, test) in enumerate(splits, 1):
        sel = screen_indices(
            lambda codes: encode(codes)[train],
            y[train],
            n_rounds=n_rounds,
            cv=inner_cv,
            candidates=candidates,
            components=components,
            n_jobs=n_jobs,
            verbose=False,
        )
        X = encode(sel.indices)
        oof[test] = fit_predict(X[train], y[train], X[test], sel.n_components)
        picked.append(sel.indices)
        print(f"  outer fold {i}/{len(splits)}: picked {sel.label}")
    return r2(y, oof), float(((y - oof) ** 2).mean()), oof, picked
