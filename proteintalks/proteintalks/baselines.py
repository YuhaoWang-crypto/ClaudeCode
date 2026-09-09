"""Baselines for the drug-efficacy task.

Two families:

1. The paper's own comparators (Methods, "Implementation of relevant methods
   for benchmarking"): Bootstrap/Bagging, random forest, logistic regression,
   SGD, KNN, and DeepSynergy.  Hyperparameters are copied from the paper.

2. Trivial baselines that the paper does not report.  These matter because the
   2025 benchmarking literature on perturbation models (Ahlmann-Eltze, Huber &
   Anders, Nat Methods 2025; Kedzierska et al., Genome Biol 2025) repeatedly
   found that deep perturbation models fail to beat constant or linear
   predictors once the evaluation split is examined.  For a drug-efficacy label
   defined per (cell line, drug), the sharpest such control is:

       DrugMeanBaseline      predict the training-set efficacy rate of that drug
       CellLineMeanBaseline  predict the training-set efficacy rate of that line

   These use NO proteomics at all.  If a proteome-driven model does not clearly
   beat them, the reported AUROC is measuring drug identity, not cell biology.
"""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier

__all__ = [
    "flatten_features",
    "sklearn_baselines",
    "DeepSynergy",
    "train_deepsynergy",
    "DrugMeanBaseline",
    "CellLineMeanBaseline",
    "GlobalPriorBaseline",
]


def flatten_features(ds, idx, use_future: bool = True) -> np.ndarray:
    """Concatenate proteome (all timepoints) with drug descriptors.

    This mirrors the paper's statement that for DeepSynergy it "concatenat[ed]
    data from multiple time points and combin[ed] it with drug data".
    """
    parts = [ds.p0[idx]]
    if use_future:
        parts.append(ds.p_future[idx].reshape(len(idx), -1))
    parts.append(ds.drug_feats[idx].reshape(len(idx), -1))
    return np.concatenate(parts, axis=1).astype(np.float32)


# --------------------------------------------------------------------------- #
# Classical ML, hyperparameters as specified in the paper's Methods
# --------------------------------------------------------------------------- #
def sklearn_baselines(seed: int = 0) -> dict:
    return {
        "Bootstrap (Bagging+DT)": BaggingClassifier(
            estimator=DecisionTreeClassifier(random_state=seed),
            n_estimators=100,
            bootstrap=True,
            random_state=seed,
        ),
        "Random forest": RandomForestClassifier(random_state=seed),
        "Logistic regression": LogisticRegression(max_iter=1000),
        "SGD": SGDClassifier(max_iter=1000, loss="log_loss", random_state=seed),
        "KNN": KNeighborsClassifier(),
    }


# --------------------------------------------------------------------------- #
# DeepSynergy
# --------------------------------------------------------------------------- #
class DeepSynergy(nn.Module):
    """Preuer et al. (2018) architecture: 8182-8182-1 fully connected, ReLU.

    The paper says it "modified the input data dimensions of the DeepSynergy
    code to fit our dataset" and kept other parameters as the original.
    """

    def __init__(self, in_dim: int, hidden=(8182, 8182), dropout: float = 0.5):
        super().__init__()
        layers, prev = [], in_dim
        layers += [nn.Dropout(0.2)]        # input dropout, as in the original
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers += [nn.Linear(prev, 1)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_deepsynergy(
    x_tr, y_tr, x_va, y_va, x_te,
    hidden=(2048, 2048), epochs: int = 200, lr: float = 1e-5,
    batch_size: int = 64, device: str = "cpu", seed: int = 0, verbose: bool = False,
):
    """Train DeepSynergy with early stopping on the validation loss.

    ``hidden`` defaults to a narrower net than the original 8182-8182 purely for
    CPU tractability; pass ``hidden=(8182, 8182)`` to match the paper exactly.
    """
    torch.manual_seed(seed)
    model = DeepSynergy(x_tr.shape[1], hidden=hidden).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    lossf = nn.BCEWithLogitsLoss()

    xt = torch.tensor(x_tr, device=device)
    yt = torch.tensor(y_tr, device=device)
    xv = torch.tensor(x_va, device=device)
    yv = torch.tensor(y_va, device=device)
    xe = torch.tensor(x_te, device=device)

    best, best_state, patience, bad = np.inf, None, 25, 0
    for ep in range(epochs):
        model.train()
        perm = torch.randperm(len(xt), device=device)
        for i in range(0, len(xt), batch_size):
            b = perm[i : i + batch_size]
            opt.zero_grad()
            loss = lossf(model(xt[b]), yt[b])
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = lossf(model(xv), yv).item()
        if vl < best - 1e-5:
            best, bad = vl, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= patience:
                break
        if verbose and ep % 20 == 0:
            print(f"  DeepSynergy ep{ep} val_loss={vl:.4f}")

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(xe)).cpu().numpy()


# --------------------------------------------------------------------------- #
# Trivial baselines (not in the paper)
# --------------------------------------------------------------------------- #
class _GroupRateBaseline:
    """Predict the training-set positive rate of a grouping variable."""

    attr = "drug"

    def fit(self, ds, train_idx):
        key = getattr(ds, self.attr)[train_idx]
        y = ds.label[train_idx]
        self.global_rate = float(y.mean())
        self.rates = {}
        for g in np.unique(key):
            m = key == g
            self.rates[int(g)] = float(y[m].mean())
        return self

    def predict_proba(self, ds, idx):
        key = getattr(ds, self.attr)[idx]
        return np.array(
            [self.rates.get(int(g), self.global_rate) for g in key], dtype=float
        )


class DrugMeanBaseline(_GroupRateBaseline):
    """Efficacy rate of this drug across training cell lines. Uses no proteome."""

    attr = "drug"


class CellLineMeanBaseline(_GroupRateBaseline):
    """Efficacy rate of this cell line across training drugs. Uses no proteome."""

    attr = "cell_line"


class GlobalPriorBaseline:
    """Constant predictor: the training-set positive rate. AUROC is 0.5 by construction."""

    def fit(self, ds, train_idx):
        self.rate = float(ds.label[train_idx].mean())
        return self

    def predict_proba(self, ds, idx):
        return np.full(len(idx), self.rate, dtype=float)
