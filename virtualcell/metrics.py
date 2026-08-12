"""Virtual Cell Challenge metrics, ported from Arc's ``cell-eval`` (v0.8.2).

The three headline metrics in the challenge commentary are the perturbation
discrimination score, the differential expression score, and mean absolute
error.  Definitions here follow ``cell_eval.metrics`` so that numbers computed
against public data are comparable in kind to leaderboard numbers:

``discrimination_score``
    For each perturbation, rank every *real* perturbation effect by its L1
    distance to the *predicted* effect for that perturbation, then report
    ``1 - rank / n_perturbations``.  The knocked-down gene itself is dropped
    from the feature set, so a model cannot win by only nailing on-target
    knockdown.  1.0 is perfect, 0.5 is chance.

``mae``
    Mean absolute error on absolute expression (not the delta).

``overlap_at_k`` / ``precision_at_k``
    Fraction of the top-k differentially expressed genes recovered, where genes
    are ranked by |log2 fold change| after filtering to significant ones.

``vcc_score``
    The leaderboard aggregation: every metric is renormalised against a
    baseline model, ``1 - user/base`` when zero is best and
    ``(user - base) / (1 - base)`` when one is best, clipped below at zero.

One deviation is documented rather than hidden.  The challenge computes
significance from submitted single cells; this pipeline works on pseudobulk, so
significance for the *measured* data comes from an empirical null built out of
the non-targeting control replicates (below), and predictions are ranked by
predicted |log2 fold change| without a significance filter of their own.  That
is the pseudobulk analogue of ``overlap_at_k`` and it is applied identically to
every model and baseline being compared.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.spatial.distance import cdist

LOG2 = np.log(2.0)


# --------------------------------------------------------------------------
# expression-space metrics
# --------------------------------------------------------------------------

def discrimination_score(pred_delta: np.ndarray, real_delta: np.ndarray,
                         pert_names: np.ndarray, gene_symbols: np.ndarray,
                         metric: str = "l1") -> np.ndarray:
    """Per-perturbation discrimination score; 1.0 best, 0.5 chance.

    Port of ``cell_eval.metrics._anndata.discrimination_score`` with
    ``exclude_target_gene=True``.
    """
    if metric != "l1":
        raise ValueError("only the challenge default (l1) is implemented")

    n_pert = pert_names.size
    pred = np.ascontiguousarray(pred_delta, dtype=np.float32)
    real = np.ascontiguousarray(real_delta, dtype=np.float32)

    # dist[i, j] = L1 between predicted effect i and measured effect j
    dist = cdist(pred, real, "cityblock")

    # Excluding the target gene removes exactly one column from each row's sum,
    # so subtract its contribution instead of rebuilding the matrix per row.
    sym_to_col = {s: i for i, s in enumerate(gene_symbols)}
    for i, p in enumerate(pert_names):
        drop = sym_to_col.get(p)
        if drop is not None:
            dist[i] -= np.abs(real[:, drop] - pred[i, drop])

    order = np.argsort(dist, axis=1, kind="stable")
    rank = np.argmax(order == np.arange(n_pert)[:, None], axis=1)
    return 1.0 - rank / n_pert


def mae(pred_expr: np.ndarray, real_expr: np.ndarray) -> np.ndarray:
    """Mean absolute error per perturbation, on absolute expression."""
    return np.abs(pred_expr - real_expr).mean(axis=1)


def mae_delta(pred_delta: np.ndarray, real_delta: np.ndarray) -> np.ndarray:
    return np.abs(pred_delta - real_delta).mean(axis=1)


def mse(pred_expr: np.ndarray, real_expr: np.ndarray) -> np.ndarray:
    return ((pred_expr - real_expr) ** 2).mean(axis=1)


def pearson_delta(pred_delta: np.ndarray, real_delta: np.ndarray) -> np.ndarray:
    """Pearson correlation between predicted and measured effect, per perturbation."""
    out = np.zeros(pred_delta.shape[0])
    for i in range(pred_delta.shape[0]):
        a, b = pred_delta[i], real_delta[i]
        if a.std() < 1e-12 or b.std() < 1e-12:
            out[i] = 0.0
        else:
            out[i] = float(np.corrcoef(a, b)[0, 1])
    return np.nan_to_num(out)


# --------------------------------------------------------------------------
# differential expression
# --------------------------------------------------------------------------

def control_variance_model(ctrl: np.ndarray, ctrl_ncells: np.ndarray
                           ) -> tuple[np.ndarray, np.ndarray]:
    """Split control-replicate variance into batch and per-cell components.

    A pseudobulk profile built from ``n`` cells has variance
    ``batch + sampling / n``: the batch term is gem-group-to-gem-group drift that
    does not shrink with more cells, the sampling term does.  Control replicates
    span a range of cell counts, so regressing their squared deviations on
    ``1/n`` separates the two.  Without this split, perturbations profiled from
    ~50 cells would be judged against a null built from control replicates of
    ~1000 cells and almost every gene would look significant.

    Returns the per-gene ``(batch, sampling)`` variance components.
    """
    dev = ctrl - ctrl.mean(axis=0, keepdims=True)
    n_rep = ctrl.shape[0]
    # unbiased-ish rescaling for having estimated the mean from the same rows
    y = dev ** 2 * (n_rep / max(n_rep - 1, 1))
    x = 1.0 / np.maximum(ctrl_ncells.astype(float), 1.0)

    xc = x - x.mean()
    denom = float((xc ** 2).sum())
    if denom < 1e-12:                      # all replicates the same size
        batch = y.mean(axis=0)
        return batch, np.zeros_like(batch)

    slope = (xc[:, None] * y).sum(axis=0) / denom       # sampling component
    intercept = y.mean(axis=0) - slope * x.mean()       # batch component
    return np.maximum(intercept, 0.0), np.maximum(slope, 0.0)


def _expression_trend(var: np.ndarray, mu: np.ndarray, n_bins: int = 25) -> np.ndarray:
    """Median variance among genes of comparable expression.

    Variance in count data tracks expression level, so the right thing to
    moderate a noisy per-gene estimate towards is the trend at that expression,
    not a global constant.  A global constant would inflate the variance of
    lowly expressed genes by orders of magnitude and silence them entirely.
    """
    order = np.argsort(mu)
    trend = np.empty_like(var)
    for chunk in np.array_split(order, n_bins):
        if chunk.size:
            trend[chunk] = np.median(var[chunk])
    return trend


def de_significance(delta: np.ndarray, ncells: np.ndarray,
                    batch_var: np.ndarray, samp_var: np.ndarray,
                    mu: np.ndarray, fdr: float = 0.05,
                    moderation: float = 0.3) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided z-test of each perturbation effect against the control null.

    Returns ``(significant_mask, log2_fold_change)``.
    """
    # Moderate each variance component towards its expression-matched trend, in
    # the spirit of limma's eBayes, so a gene that happens to look quiet across
    # control replicates cannot produce a spuriously enormous z-score.
    batch_var = ((1 - moderation) * batch_var
                 + moderation * _expression_trend(batch_var, mu))
    samp_var = ((1 - moderation) * samp_var
                + moderation * _expression_trend(samp_var, mu))

    var = batch_var[None, :] + samp_var[None, :] / np.maximum(ncells, 1)[:, None]
    sd = np.sqrt(np.maximum(var, 1e-12))

    z = delta / sd
    p = 2.0 * stats.norm.sf(np.abs(z))
    sig = np.zeros_like(p, dtype=bool)
    for i in range(p.shape[0]):
        sig[i] = _bh(p[i], fdr)
    return sig, delta / LOG2      # natural-log delta -> log2 fold change


def _bh(pvals: np.ndarray, alpha: float) -> np.ndarray:
    """Benjamini-Hochberg rejection mask."""
    n = pvals.size
    order = np.argsort(pvals)
    thresh = alpha * np.arange(1, n + 1) / n
    passing = pvals[order] <= thresh
    out = np.zeros(n, dtype=bool)
    if passing.any():
        cutoff = np.flatnonzero(passing)[-1]
        out[order[:cutoff + 1]] = True
    return out


def overlap_at_k(pred_lfc: np.ndarray, real_lfc: np.ndarray,
                 real_sig: np.ndarray, k: int | None = 100) -> np.ndarray:
    """Top-k differential expression overlap, per perturbation.

    Real genes are those passing FDR, ranked by |log2FC|; predicted genes are
    ranked by predicted |log2FC|.  Score is ``|intersection| / k_eff`` with
    ``k_eff = min(k, n_real_significant)``, matching ``cell-eval``'s
    ``overlap_at_k``.
    """
    n_pert = pred_lfc.shape[0]
    out = np.zeros(n_pert)
    for i in range(n_pert):
        real_idx = np.flatnonzero(real_sig[i])
        if real_idx.size == 0:
            out[i] = 0.0
            continue
        real_idx = real_idx[np.argsort(-np.abs(real_lfc[i, real_idx]))]
        k_eff = real_idx.size if k is None else min(k, real_idx.size)
        pred_idx = np.argsort(-np.abs(pred_lfc[i]))[:k_eff]
        out[i] = np.intersect1d(real_idx[:k_eff], pred_idx).size / k_eff
    return out


def de_direction_match(pred_lfc: np.ndarray, real_lfc: np.ndarray,
                       real_sig: np.ndarray) -> np.ndarray:
    """Fraction of truly-significant genes whose predicted sign is correct."""
    out = np.zeros(pred_lfc.shape[0])
    for i in range(pred_lfc.shape[0]):
        idx = np.flatnonzero(real_sig[i])
        if idx.size == 0:
            out[i] = np.nan
            continue
        out[i] = float(np.mean(np.sign(pred_lfc[i, idx]) == np.sign(real_lfc[i, idx])))
    return out


def de_spearman_lfc(pred_lfc: np.ndarray, real_lfc: np.ndarray,
                    real_sig: np.ndarray) -> np.ndarray:
    """Spearman correlation of log2FC restricted to truly-significant genes."""
    out = np.zeros(pred_lfc.shape[0])
    for i in range(pred_lfc.shape[0]):
        idx = np.flatnonzero(real_sig[i])
        if idx.size < 3:
            out[i] = np.nan
            continue
        out[i] = float(stats.spearmanr(pred_lfc[i, idx], real_lfc[i, idx]).statistic)
    return np.nan_to_num(out)


# --------------------------------------------------------------------------
# leaderboard aggregation
# --------------------------------------------------------------------------

BEST_ZERO = {"mae", "mae_delta", "mse"}
BEST_ONE = {"discrimination_score_l1", "pearson_delta", "overlap_at_100",
            "overlap_at_50", "de_direction_match", "de_spearman_lfc"}

#: The three metrics the challenge commentary names, one per row of its Figure
#: 1C: differential expression score, perturbation discrimination score, and
#: mean absolute error.  Scoring on these keeps MAE at a third of the total.
#: Averaging the whole suite instead would put four DE-derived metrics against
#: one error metric and quietly make error almost free to give up.
HEADLINE = ("discrimination_score_l1", "overlap_at_100", "mae")


def vcc_score(user: dict[str, float], base: dict[str, float],
              clip: bool = True,
              metrics: tuple[str, ...] | None = HEADLINE) -> dict[str, float]:
    """Renormalise a model's metrics against a baseline, as ``cell_eval._score``.

    ``clip=True`` reproduces the leaderboard exactly: a metric worse than the
    reference model contributes zero rather than a negative number.  That floor
    is right for *reporting* but wrong for *tuning* -- once a metric is already
    worse than baseline, making it arbitrarily worse costs nothing, so a search
    will happily trade all of MAE away for a little discrimination.  The
    challenge guards against that by enforcing minimum thresholds on every
    metric ("discouraging models that perform well on one metric at the expense
    of the others"); since the thresholds are not published, model selection
    here uses ``clip=False``, where overshooting is genuinely expensive.

    ``metrics`` defaults to :data:`HEADLINE`, the three the commentary names.
    Pass ``None`` to average every metric present.
    """
    out: dict[str, float] = {}
    for name, value in user.items():
        if metrics is not None and name not in metrics:
            continue
        b = base[name]
        if name in BEST_ZERO:
            s = 1.0 - value / b if b != 0 else 0.0
        elif name in BEST_ONE:
            s = (value - b) / (1.0 - b) if b != 1 else 0.0
        else:
            continue
        s = float(np.nan_to_num(s))
        out[name] = float(np.clip(s, 0.0, None)) if clip else s
    out["avg_score"] = float(np.mean(list(out.values()))) if out else 0.0
    return out
