"""Evaluation metrics matching the paper's reported quantities.

Classification (drug efficacy / synergy): AUROC, AUPRC, accuracy at 0.5.
Proteome trajectory: MSE, per-condition and per-protein Pearson r, and -- most
importantly -- MSE relative to a *no-change* control, because a min-max
normalised proteome barely moves between 6 h and 48 h and an unnormalised MSE
can look excellent while carrying no dynamical information.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    roc_auc_score,
)

__all__ = ["classification_metrics", "trajectory_metrics", "bootstrap_ci"]


def classification_metrics(y_true, y_prob, threshold: float = 0.5) -> dict:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob, dtype=float)
    out = {"n": int(len(y_true)), "positive_rate": float(y_true.mean())}
    if len(np.unique(y_true)) < 2:
        out.update(auroc=float("nan"), auprc=float("nan"))
    else:
        out["auroc"] = float(roc_auc_score(y_true, y_prob))
        out["auprc"] = float(average_precision_score(y_true, y_prob))
    out["accuracy"] = float(accuracy_score(y_true, (y_prob >= threshold).astype(int)))
    return out


def trajectory_metrics(y_true, y_pred, p0=None) -> dict:
    """y_true / y_pred: (n, T, N).  ``p0``: (n, N) for the no-change control."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mse = float(np.mean((y_true - y_pred) ** 2))
    out = {"mse": mse, "rmse": float(np.sqrt(mse))}

    n, T, N = y_true.shape
    # Per-condition correlation across proteins, averaged over conditions/times.
    rs = []
    for i in range(n):
        for t in range(T):
            a, b = y_true[i, t], y_pred[i, t]
            if a.std() > 1e-9 and b.std() > 1e-9:
                rs.append(stats.pearsonr(a, b)[0])
    out["pearson_r_per_condition"] = float(np.mean(rs)) if rs else float("nan")

    # The control that actually matters: how much better than "nothing changed"?
    if p0 is not None:
        p0 = np.asarray(p0, dtype=float)
        base = np.repeat(p0[:, None, :], T, axis=1)
        mse_base = float(np.mean((y_true - base) ** 2))
        out["mse_nochange_baseline"] = mse_base
        out["skill_vs_nochange"] = float(1.0 - mse / mse_base) if mse_base > 0 else float("nan")

        # Delta-space correlation: correlate the *change* from baseline, which
        # removes the trivially predictable static component of the proteome.
        d_true = (y_true - base).reshape(n * T, N)
        d_pred = (y_pred - base).reshape(n * T, N)
        rs = [
            stats.pearsonr(a, b)[0]
            for a, b in zip(d_true, d_pred)
            if a.std() > 1e-9 and b.std() > 1e-9
        ]
        out["delta_pearson_r"] = float(np.mean(rs)) if rs else float("nan")
    return out


def bootstrap_ci(y_true, y_prob, metric: str = "auroc", n_boot: int = 1000, seed: int = 0):
    """Percentile bootstrap 95% CI for a classification metric."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob, dtype=float)
    fn = {
        "auroc": roc_auc_score,
        "auprc": average_precision_score,
        "accuracy": lambda a, b: accuracy_score(a, (b >= 0.5).astype(int)),
    }[metric]
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        vals.append(fn(y_true[idx], y_prob[idx]))
    if not vals:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
