"""Independent re-evaluation harness for the VC-suite QSAR models.

Tests three things the suite's own report does not:
  1. split severity  -- random vs Murcko-scaffold vs sphere-exclusion cluster vs temporal
  2. q-RASAR gain    -- does the similarity-feature boost survive a stricter split?
  3. applicability   -- performance as a continuous function of nearest-neighbour similarity
"""
import os, math, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 42
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# ---------------------------------------------------------------- features
def fps(smiles):
    X = np.zeros((len(smiles), 2048), dtype=np.float32)
    ok = np.zeros(len(smiles), dtype=bool)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        arr = GEN.GetFingerprintAsNumPy(m)
        X[i] = arr
        ok[i] = True
    return X, ok


def tanimoto_block(A, B, blk=2000):
    """Exact Tanimoto between binary float32 fingerprint matrices."""
    na, nb = A.sum(1), B.sum(1)
    out = np.zeros((A.shape[0], B.shape[0]), dtype=np.float32)
    for i in range(0, A.shape[0], blk):
        j = min(i + blk, A.shape[0])
        inter = A[i:j] @ B.T
        out[i:j] = inter / (na[i:j, None] + nb[None, :] - inter + 1e-9)
    return out


# ---------------------------------------------------------------- splits
def sphere_exclusion(X, sim_cut=0.4, seed=SEED):
    """Leader clustering: no two cluster leaders are more similar than sim_cut."""
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    order = np.argsort(-X.sum(1))  # densest fingerprints first, deterministic
    labels = np.full(n, -1)
    leaders = []
    for idx in order:
        if labels[idx] != -1:
            continue
        if leaders:
            s = tanimoto_block(X[idx:idx + 1], X[leaders])[0]
            if s.max() >= sim_cut:
                labels[idx] = int(np.argmax(s))
                continue
        labels[idx] = len(leaders)
        leaders.append(idx)
    # assign remaining
    L = X[leaders]
    un = np.where(labels == -1)[0]
    if len(un):
        s = tanimoto_block(X[un], L)
        labels[un] = s.argmax(1)
    return labels


def group_folds(groups, k=5, seed=SEED):
    """Balanced group k-fold: assign whole groups to folds, largest first."""
    rng = np.random.default_rng(seed)
    uniq, counts = np.unique(groups, return_counts=True)
    order = np.argsort(-counts)
    sizes = np.zeros(k)
    gmap = {}
    for gi in order:
        f = int(np.argmin(sizes))
        gmap[uniq[gi]] = f
        sizes[f] += counts[gi]
    fold = np.array([gmap[g] for g in groups])
    return [(np.where(fold != f)[0], np.where(fold == f)[0]) for f in range(k)]


def random_folds(n, k=5, seed=SEED):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return [(np.setdiff1d(idx, p), p) for p in np.array_split(idx, k)]


def temporal_split(years, frac=0.3):
    y = pd.Series(years).fillna(9999).values.astype(float)
    cut = np.quantile(y[y < 9999], 1 - frac)
    tr = np.where(y <= cut)[0]
    te = np.where(y > cut)[0]
    return [(tr, te)] if len(te) > 30 and len(tr) > 50 else []


# ---------------------------------------------------------------- q-RASAR
def _topk_stable(T, k, margin=64, block=1024):
    """Top-k reference indices per row, ties broken by ascending reference index
    (the report's stated (-t, idx) contract). Row-blocked so memory stays O(block*n_ref)
    instead of materialising an n_query x n_ref index matrix."""
    n, m = T.shape
    out = np.empty((n, k), dtype=np.int64)
    cut = min(margin, m)
    for i in range(0, n, block):
        j = min(i + block, n)
        B = T[i:j]
        if m <= cut:
            order = np.argsort(-B, axis=1, kind="stable")[:, :k]
        else:
            cand = np.argpartition(-B, cut - 1, axis=1)[:, :cut]
            cv = np.take_along_axis(B, cand, axis=1)
            # order the candidate window by (-sim, reference index); argpartition
            # returns ties in arbitrary order, so index must be an explicit key
            loc = np.lexsort((cand, -cv), axis=1)
            cand = np.take_along_axis(cand, loc, axis=1)
            cv = np.take_along_axis(cv, loc, axis=1)
            # a tie spanning the candidate cut could hide a lower-index equal neighbour
            bad = cv[:, k - 1] == cv[:, cut - 1]
            order = cand[:, :k]
            if bad.any():
                rows = np.where(bad)[0]
                order = order.copy()
                order[rows] = np.argsort(-B[rows], axis=1, kind="stable")[:, :k]
        # within each row the candidate indices must still ascend on ties
        out[i:j] = order
    return out


def qrasar_feats(T, ytr, k=5, loo=False):
    """T: (n_query, n_ref) Tanimoto. Returns 6 similarity features (q-RASAR style).
    Modifies T in place when loo=True (caller owns a scratch matrix)."""
    if loo:
        np.fill_diagonal(T, -1.0)
    kk = max(k, 5)
    idx = _topk_stable(T, kk)
    S = np.take_along_axis(T, idx, axis=1)
    L = ytr[idx]
    w = np.exp(-(1 - S) / 0.1)
    w = w / (w.sum(1, keepdims=True) + 1e-12)
    return np.column_stack([
        S[:, 0], S[:, 1], S[:, :3].mean(1), S[:, :5].mean(1),
        (w * L).sum(1), L[:, :5].std(1)
    ]).astype(np.float32)


# ---------------------------------------------------------------- metrics
def ef_at(y_true, y_pred, top_frac=0.05):
    """Enrichment: of the top predicted `top_frac`, how many are in the true top `top_frac`."""
    n = len(y_true)
    k = max(1, int(round(n * top_frac)))
    true_top = set(np.argsort(-y_true)[:k])
    pred_top = np.argsort(-y_pred)[:k]
    hit = sum(1 for i in pred_top if i in true_top)
    return (hit / k) / top_frac


def metrics(y, p, tag=""):
    r = spearmanr(y, p).statistic
    rmse = float(np.sqrt(np.mean((y - p) ** 2)))
    ss = 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
    return dict(tag=tag, n=len(y), spearman=float(r), rmse=rmse, r2=float(ss),
                ef5=float(ef_at(y, p, 0.05)), ef1=float(ef_at(y, p, 0.01)),
                within1log=float(np.mean(np.abs(y - p) <= 1.0)),
                within2fold=float(np.mean(np.abs(y - p) <= math.log10(2))))


# ---------------------------------------------------------------- model
def make_model(n, kind="auto"):
    """N_TREES env var lets the large endpoints run at the 100-tree size the suite's
    pages actually embed (the report measures r>=0.999 against its 300-tree parents)."""
    if kind == "ridge" or (kind == "auto" and n < 400):
        return Ridge(alpha=10.0)
    return RandomForestRegressor(n_estimators=int(os.environ.get("N_TREES", 300)),
                                 n_jobs=4, random_state=SEED, min_samples_leaf=1)


def run_cv(X, y, folds, use_qrasar=False, kind="auto"):
    oof = np.full(len(y), np.nan)
    nnsim = np.full(len(y), np.nan)
    for tr, te in folds:
        T_te = tanimoto_block(X[te], X[tr])
        nnsim[te] = T_te.max(1)
        if use_qrasar:
            T_tr = tanimoto_block(X[tr], X[tr])
            Ftr = qrasar_feats(T_tr, y[tr], loo=True)
            Fte = qrasar_feats(T_te, y[tr], loo=False)
            Xtr = np.hstack([X[tr], Ftr]); Xte = np.hstack([X[te], Fte])
        else:
            Xtr, Xte = X[tr], X[te]
        m = make_model(len(tr), kind)
        m.fit(Xtr, y[tr])
        oof[te] = m.predict(Xte)
    ok = ~np.isnan(oof)
    return oof, nnsim, ok
