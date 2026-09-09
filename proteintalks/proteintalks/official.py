"""Exact port of the released reference implementation (guomics-lab/PTV-1, MIT).

The reference code (``ProteinTalks/model.py``, class ``ppODE``) differs from the
paper's Equations 1–11 in several ways that change what the model can represent.
This module reproduces the *code*, so both can be run side by side against the
released checkpoint. ``proteintalks/model.py`` reproduces the *equations*.

Differences from the published equations, all verified against the source
-------------------------------------------------------------------------
1. **`C1`/`C2` use `kernel_size=1`.** They mix channels, never proteins. Combined
   with `linear_input` (per-protein), an ODE field that is `nn.Linear` over the
   channel dimension, and a per-protein `layer_final`, **module 1 is entirely
   protein-independent**: each protein's 6/24/48 h trajectory is a function of
   its own baseline abundance and its own perturbation entry, through a shared
   map. No protein-protein coupling exists anywhere in the dynamics module.
   Cross-protein mixing appears only in `drugsens_conv1`, in the phenotype head.
   Kernel 1 is the defensible choice given that protein row order is arbitrary,
   but it means the "protein network dynamics" the model is said to learn is not
   represented in the dynamics module.

2. **The ODE integrates over `torch.linspace(0, 3, 4)`** — integer ticks
   0, 1, 2, 3 — not over 0, 6, 24, 48 hours. The solver never sees that the
   24→48 h gap is four times the 0→6 h gap, so the "continuous time" component
   is a uniform three-step recurrence, not a model of real elapsed time.

3. **`C3`/`C4` output 32 channels, not 128** as Eq. 7–8 state
   (`mid_feats_drugsens = mid_feats_drugs = 32`).

4. **`drugsens_conv1` is `Conv1d(n_proteins → 32, kernel_size=4)`**: proteins are
   *channels* and the four timepoints are the spatial axis. So the proteome head
   is a linear map over all 5,585 proteins with a length-4 kernel over time.
   Similarly `drugs_conv2 = Conv1d(935 → 32, kernel_size=2)` over the drug pair.

5. **Default `hidden_size` is 64**, not the 128 of Eq. 1's "convolutional network
   C1 … to 128 dimensions"; the released config default is 64.

6. **Default `dropout_rate` is 0.0**, not the 0.1 of Eq. 3.

7. **LayerNorm and GroupNorm are used throughout** and appear nowhere in the paper.

8. **The second ODE layer has no activation** (Eq. 3 specifies Softplus), and
   `conv1` is followed by ReLU (the paper implies Softplus).

9. **SWAG** (Stochastic Weight Averaging-Gaussian) is implemented and the
   published Supplementary Table S5C reports `ppODE (SWA)` and `ppODE (Non-SWA)`
   separately, but SWAG is not described in the paper's Methods.

10. **Task weights are recomputed from `[1-λ, λ]` at every step**, not
    accumulated across steps as Eq. 16–18 imply, and per-tensor gradient clipping
    is applied *before* the cosine similarity is measured.

Training hyperparameters from the released `config.py`, none of which appear in
the paper: AdamW, lr 5e-4, weight decay 1e-4, batch size 64, 1000 epochs,
patience 500, 10 warm-up epochs, grad-norm clip 1.0, seed 1995.
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

__all__ = ["ppODE", "FullyConnectedLayer", "OFFICIAL_DEFAULTS", "rk4_classical"]

OFFICIAL_DEFAULTS = {
    "hidden_size": 64,
    "dropout_rate": 0.0,
    "optimizer": "adamw",
    "learning_rate": 5e-4,
    "weight_decay": 1e-4,
    "batch_size": 64,
    "total_epoch": 1000,
    "patience": 500,
    "warmup_epochs": 10,
    "clip_grad_norm": 1.0,
    "lambda_pheno": 0.8,
    "random_seed": 1995,
    "solver": "rk4",
}


class FullyConnectedLayer(nn.Module):
    """Reference `FullyConnectedLayer`: dropout -> linear -> LayerNorm -> act."""

    def __init__(self, in_feats, out_feats, activation, dropout, bias=True):
        super().__init__()
        self.fc = nn.Linear(in_feats, out_feats, bias=bias)
        self.norm = nn.LayerNorm(out_feats)
        self.activation = activation
        self.dropout = nn.Dropout(p=dropout) if dropout else None
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1.0 / math.sqrt(self.fc.weight.size(1))
        self.fc.weight.data.uniform_(-stdv, stdv)
        if self.fc.bias is not None:
            self.fc.bias.data.uniform_(-stdv, stdv)

    def forward(self, h):
        if self.dropout:
            h = self.dropout(h)
        h = self.fc(h)
        h = self.norm(h)
        if self.activation:
            h = self.activation(h)
        return h


def rk4_classical(func, y0, t):
    """Fixed-step classical Runge-Kutta 4, one step per interval of ``t``.

    This is deliberately *not* ``torchdiffeq.odeint(..., method="rk4")``.
    torchdiffeq's ``rk4`` implements the 3/8-rule tableau, whereas torchdyn's
    ``solver="rk4"`` -- what the reference implementation uses -- is the classical
    tableau. Both are fourth-order, so they agree as the step size shrinks, but
    this model integrates with h = 1 over only three steps, and there the two
    tableaux differ by ~0.2% in the trajectory and ~2.6% in the decoded proteome.

    "rk4", on its own, therefore does not pin down this model. Reproducing the
    reference requires the classical tableau specifically.
    """
    ys = [y0]
    y = y0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = func(t[i], y)
        k2 = func(t[i] + h / 2, y + h / 2 * k1)
        k3 = func(t[i] + h / 2, y + h / 2 * k2)
        k4 = func(t[i] + h, y + h * k3)
        y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        ys.append(y)
    return torch.stack(ys)


class _ODEFunc(nn.Module):
    """The reference vector field.

    The reference uses `torchdyn.core.NeuralODE(func, solver='rk4')`, whose
    `func` is autonomous (no explicit t argument). We accept and ignore t.
    """

    def __init__(self, hidden_feats, dropout):
        super().__init__()
        self.net = nn.Sequential(
            FullyConnectedLayer(hidden_feats, hidden_feats, nn.Softplus(), dropout),
            FullyConnectedLayer(hidden_feats, hidden_feats, None, dropout),
        )

    def forward(self, t, h):
        return self.net(h)


class ppODE(nn.Module):
    """Reference architecture, with torchdiffeq standing in for torchdyn.

    Parameters mirror the reference signature: ``pro_feats`` is the number of
    proteins (5,585 in the released checkpoint) and ``drug_feature_feats`` is 935.
    """

    def __init__(
        self,
        pro_feats: int,
        hidden_feats: int = 64,
        out_feats: int = 1,
        drug_feature_feats: int = 935,
        dropout: float = 0.0,
        time_tick_num: int = 4,
    ):
        super().__init__()
        self.mid_feats = 32
        self.mid_feats_drugsens = 32
        self.mid_feats_drugs = 32
        self.mid_feats_drugsens_drugs = 32
        self.pro_feats = pro_feats
        self.time_tick_num = time_tick_num

        self.linear_input = FullyConnectedLayer(2, self.mid_feats, nn.Softplus(), dropout)
        self.conv1 = nn.Conv1d(self.mid_feats, hidden_feats, kernel_size=1)
        self.conv1_norm = nn.GroupNorm(1, hidden_feats)
        self.conv2 = nn.Conv1d(hidden_feats, self.mid_feats, kernel_size=1)
        self.convdrug1 = nn.Conv1d(2, hidden_feats, kernel_size=2)  # unused, as in ref

        self.ode_func = _ODEFunc(hidden_feats, dropout)
        self.layer_final = nn.Linear(self.mid_feats, out_feats)

        self.drugsens_conv1 = nn.Conv1d(pro_feats, self.mid_feats_drugsens, kernel_size=4)
        self.drugs_conv2 = nn.Conv1d(drug_feature_feats, self.mid_feats_drugs, kernel_size=2)
        self.pheno_fc1 = nn.Linear(
            self.mid_feats_drugsens + self.mid_feats_drugs, self.mid_feats_drugsens_drugs
        )
        self.pheno_fc2 = nn.Linear(self.mid_feats_drugsens_drugs, 1)

    def forward(self, x, pert, fp_phA, fp_phB):
        """x, pert: (B, n_proteins, 1); fp_phA/B: (B, 935, 1)."""
        emb = self.linear_input(torch.cat([x, pert], dim=-1))    # (B, N, 32)
        emb = torch.transpose(emb, 1, 2)                          # (B, 32, N)
        emb = F.relu(self.conv1_norm(self.conv1(emb)))            # (B, H, N)
        emb = torch.transpose(emb, 1, 2)                          # (B, N, H)

        t = torch.linspace(0, self.time_tick_num - 1, self.time_tick_num,
                           device=x.device, dtype=x.dtype)
        traj = rk4_classical(self.ode_func, emb, t)[1:]             # (3, B, N, H)

        T, B, N, H = traj.shape
        z = self.conv2(traj.reshape(T * B, N, H).transpose(1, 2)).transpose(1, 2)
        y = self.layer_final(z).reshape(T, B, N, -1)               # (3, B, N, 1)
        y = y.permute(1, 0, 2, 3).squeeze(-1)                      # (B, 3, N)

        # [P0, P6, P24, P48] as 4 positions along the conv's spatial axis.
        xy = torch.cat([x, y.transpose(1, 2)], dim=-1)             # (B, N, 4)
        xy = F.relu(self.drugsens_conv1(xy)).squeeze(-1)           # (B, 32)

        fp = torch.cat([fp_phA, fp_phB], dim=-1)                   # (B, 935, 2)
        fp = F.relu(self.drugs_conv2(fp)).squeeze(-1)              # (B, 32)

        h = F.relu(self.pheno_fc1(torch.cat([xy, fp], dim=1)))
        pheno = torch.sigmoid(self.pheno_fc2(h)).squeeze(1)
        return y, pheno

    def n_parameters(self) -> dict[str, int]:
        head = sum(
            p.numel()
            for m in [self.drugsens_conv1, self.drugs_conv2, self.pheno_fc1, self.pheno_fc2]
            for p in m.parameters()
        )
        return {"total": sum(p.numel() for p in self.parameters()), "phenotype_head": head}
