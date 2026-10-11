"""
discriminate — honest scoring of candidate readouts.

auc()      Mann-Whitney AUC, orientation-free (reports max(AUC, 1-AUC) and the
           direction), with percentile bootstrap CI.
cv_auc()   cross-validated logistic regression AUC for a multi-feature panel
           (so a panel never gets credit for in-sample overfitting).
"""
import numpy as np
from scipy.stats import mannwhitneyu


def _auc(a, b):
    u = mannwhitneyu(a, b, alternative="two-sided").statistic
    return u / (len(a) * len(b))


def auc(pos, neg, n_boot=1000, rng=0):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    a = _auc(pos, neg)
    sign = 1 if a >= 0.5 else -1
    rng = np.random.default_rng(rng)
    boots = []
    for _ in range(n_boot):
        bp = rng.choice(pos, len(pos))
        bn = rng.choice(neg, len(neg))
        boots.append(_auc(bp, bn))
    boots = np.array(boots)
    if sign < 0:
        a, boots = 1 - a, 1 - boots
    p = mannwhitneyu(pos, neg, alternative="two-sided").pvalue
    return {"auc": a, "ci": tuple(np.percentile(boots, [2.5, 97.5])),
            "direction": "higher in positives" if sign > 0
            else "lower in positives", "p": p}


def cv_auc(X, y, folds=5, rng=0, C=1.0):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=2000))
    skf = StratifiedKFold(folds, shuffle=True, random_state=rng)
    p = cross_val_predict(clf, X, y, cv=skf, method="predict_proba")[:, 1]
    return roc_auc_score(y, p)
