"""SHAP-based protein prioritisation (paper: "SHAP value calculation/filtration").

The paper computes SHAP values for Loss2 (the efficacy head), using resistant /
non-synergistic conditions as the background set, then keeps the 30 proteins
with the most positive and the 30 with the most negative mean SHAP value.

We reproduce that procedure with ``shap.GradientExplainer`` over the proteome
input, and add one control the paper does not report: a **label-permutation
null**. Re-running the identical selection on shuffled efficacy labels shows how
large a mean |SHAP| the top-60 cut produces by chance. Without that null, a
ranked protein list is not evidence of anything -- any model with any accuracy
produces a top-60 list.
"""

from __future__ import annotations

import numpy as np
import torch

__all__ = ["proteome_shap", "top_proteins", "permutation_null"]


class _EfficacyWrapper(torch.nn.Module):
    """Exposes efficacy logits as a function of the baseline proteome alone."""

    def __init__(self, model, pert, drug_feats):
        super().__init__()
        self.model = model
        self.register_buffer("pert", pert)
        self.register_buffer("drug_feats", drug_feats)

    def forward(self, p0):
        n = p0.shape[0]
        pert = self.pert[:n] if self.pert.shape[0] >= n else self.pert.expand(n, -1)
        df = (
            self.drug_feats[:n]
            if self.drug_feats.shape[0] >= n
            else self.drug_feats.expand(n, -1, -1)
        )
        _, logits = self.model(p0, pert, df)
        return logits.unsqueeze(-1)


def proteome_shap(model, ds, target_idx, background_idx, n_background: int = 50,
                  device: str = "cpu") -> np.ndarray:
    """Mean SHAP value per protein for the conditions in ``target_idx``.

    ``background_idx`` should be the resistant / non-synergistic conditions, as
    in the paper.
    """
    import shap

    model.eval()
    bg = np.asarray(background_idx)[:n_background]
    wrapper = _EfficacyWrapper(
        model,
        torch.tensor(ds.pert[target_idx], device=device),
        torch.tensor(ds.drug_feats[target_idx], device=device),
    ).to(device)

    bg_t = torch.tensor(ds.p0[bg], device=device)
    tg_t = torch.tensor(ds.p0[target_idx], device=device)
    explainer = shap.GradientExplainer(wrapper, bg_t)
    vals = explainer.shap_values(tg_t)
    if isinstance(vals, list):
        vals = vals[0]
    vals = np.asarray(vals)
    if vals.ndim == 3:          # (n, n_proteins, 1)
        vals = vals[..., 0]
    return vals.mean(axis=0)


def top_proteins(mean_shap: np.ndarray, protein_names, k: int = 30) -> dict:
    """Paper's filtration: top-k positive and top-k negative mean SHAP values."""
    order = np.argsort(mean_shap)
    neg = order[:k]
    pos = order[-k:][::-1]
    return {
        "positive": [(protein_names[i], float(mean_shap[i])) for i in pos],
        "negative": [(protein_names[i], float(mean_shap[i])) for i in neg],
    }


def permutation_null(model, ds, target_idx, background_idx, n_permutations: int = 20,
                     k: int = 30, seed: int = 0, device: str = "cpu") -> dict:
    """How large is the top-60 mean |SHAP| under shuffled labels?

    Returns the observed statistic, the null distribution, and an empirical
    p-value. This is the control that separates "the model ranked some proteins"
    from "the ranking carries label information".
    """
    rng = np.random.default_rng(seed)
    observed = proteome_shap(model, ds, target_idx, background_idx, device=device)
    obs_stat = float(np.abs(np.sort(observed)[-k:]).mean()
                     + np.abs(np.sort(observed)[:k]).mean()) / 2.0

    null = []
    all_idx = np.concatenate([np.asarray(target_idx), np.asarray(background_idx)])
    for _ in range(n_permutations):
        perm = rng.permutation(all_idx)
        t = perm[: len(target_idx)]
        b = perm[len(target_idx):]
        v = proteome_shap(model, ds, t, b, device=device)
        null.append(
            float(np.abs(np.sort(v)[-k:]).mean() + np.abs(np.sort(v)[:k]).mean()) / 2.0
        )
    null = np.asarray(null)
    p = float((np.sum(null >= obs_stat) + 1) / (len(null) + 1))
    return {
        "observed": obs_stat,
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "p_value": p,
        "null": null.tolist(),
    }
