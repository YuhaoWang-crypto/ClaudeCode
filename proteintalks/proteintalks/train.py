"""Training loop for ProteinTalks (Eq. 6, 11, 12 and the Eq. 13-18 weighting)."""

from __future__ import annotations

import copy
import time

import numpy as np
import torch
import torch.nn as nn

from .model import ProteinTalks
from .multitask import GradientConflictWeighter

__all__ = ["train_proteintalks", "predict"]


def _tensors(ds, idx, device):
    return (
        torch.tensor(ds.p0[idx], device=device),
        torch.tensor(ds.pert[idx], device=device),
        torch.tensor(ds.drug_feats[idx], device=device),
        torch.tensor(ds.p_future[idx], device=device),
        torch.tensor(ds.label[idx], device=device),
    )


def train_proteintalks(
    ds,
    train_idx,
    val_idx,
    *,
    lam: float = 0.8,
    epochs: int = 300,
    lr: float = 1e-3,
    batch_size: int = 32,
    weight_decay: float = 1e-5,
    dynamics_backbone: str = "neural_ode",
    adaptive_weights: bool = True,
    patience: int = 40,
    device: str = "cpu",
    seed: int = 0,
    verbose: bool = True,
    hidden_dim: int = 128,
    embed_dim: int = 32,
):
    """Train the two-module model with the paper's multi-task objective.

    Returns ``(model, history)``.  Early stopping monitors the validation value
    of the *combined* objective, so neither task can be silently sacrificed.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = ProteinTalks(
        n_proteins=ds.n_proteins,
        n_drug_features=ds.drug_feats.shape[-1],
        embed_dim=embed_dim,
        hidden_dim=hidden_dim,
        dynamics_backbone=dynamics_backbone,
    ).to(device)

    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    mse = nn.MSELoss()
    bce = nn.BCEWithLogitsLoss()
    weighter = GradientConflictWeighter(lam=lam)

    # Gradient cosine similarity (Eq. 13) is measured on the shared trunk, i.e.
    # the parameters both losses actually touch. Loss2 reaches the dynamics
    # module through the predicted proteome, so the whole dynamics net is shared.
    shared_params = [p for p in model.dynamics.parameters() if p.requires_grad]

    xtr = _tensors(ds, train_idx, device)
    xva = _tensors(ds, val_idx, device)

    best_val, best_state, bad = np.inf, None, 0
    history = {"train": [], "val": [], "w1": [], "w2": [], "cos": []}
    t0 = time.time()

    for ep in range(epochs):
        model.train()
        perm = torch.randperm(len(train_idx), device=device)
        ep_losses = []
        for i in range(0, len(perm), batch_size):
            b = perm[i : i + batch_size]
            p0, d, df, pf, y = (t[b] for t in xtr)

            opt.zero_grad()
            pred_pf, logits = model(p0, d, df)
            loss1 = mse(pred_pf, pf)      # Eq. 6
            loss2 = bce(logits, y)        # Eq. 11

            if adaptive_weights and shared_params:
                # Eq. 13-18: adapt the weights from the gradient conflict.
                weighter.update(loss1, loss2, shared_params)
            loss = weighter.combine(loss1, loss2)   # Eq. 12
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            ep_losses.append(
                (loss.item(), loss1.item(), loss2.item())
            )

        model.eval()
        with torch.no_grad():
            p0, d, df, pf, y = xva
            pred_pf, logits = model(p0, d, df)
            v1 = mse(pred_pf, pf).item()
            v2 = bce(logits, y).item()
            vtot = weighter.w1 * v1 + weighter.w2 * v2

        history["train"].append(float(np.mean([l[0] for l in ep_losses])))
        history["val"].append(vtot)
        history["w1"].append(weighter.w1)
        history["w2"].append(weighter.w2)
        if weighter.history:
            history["cos"].append(weighter.history[-1]["cos"])

        if vtot < best_val - 1e-6:
            best_val, bad = vtot, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            bad += 1
            if bad >= patience:
                if verbose:
                    print(f"  early stop at epoch {ep} (best val {best_val:.5f})")
                break

        if verbose and ep % 25 == 0:
            print(
                f"  ep{ep:4d} train={history['train'][-1]:.5f} "
                f"val={vtot:.5f} (mse={v1:.5f}, bce={v2:.5f}) "
                f"w=({weighter.w1:.3f},{weighter.w2:.3f})"
            )

    if best_state is not None:
        model.load_state_dict(best_state)
    history["seconds"] = time.time() - t0
    history["n_params"] = model.n_parameters()
    history["best_val"] = best_val
    return model, history


@torch.no_grad()
def predict(model, ds, idx, device: str = "cpu", batch_size: int = 64):
    """Return ``(predicted_proteome, efficacy_probability)`` for ``idx``."""
    model.eval()
    probs, pfs = [], []
    for i in range(0, len(idx), batch_size):
        b = idx[i : i + batch_size]
        p0 = torch.tensor(ds.p0[b], device=device)
        d = torch.tensor(ds.pert[b], device=device)
        df = torch.tensor(ds.drug_feats[b], device=device)
        pf, logits = model(p0, d, df)
        pfs.append(pf.cpu().numpy())
        probs.append(torch.sigmoid(logits).cpu().numpy())
    return np.concatenate(pfs), np.concatenate(probs)
