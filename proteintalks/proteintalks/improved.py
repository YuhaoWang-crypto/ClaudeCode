"""ProteinTalks-R: the released architecture with five targeted changes.

Every change below is motivated by something this reproduction measured, or by
something the published Supplementary Information says the model *should* do but
the released code does not. Each is an independent flag so it can be ablated.

Motivation, change by change
----------------------------

**A. Residual (delta) decoding.** `delta_decode=True`

Measured problem: on our simulator the released architecture's proteome
prediction scores **-105 against a no-change baseline**, i.e. its trajectory MSE
is 106x worse than simply asserting the proteome does not move. The decoder has
to reconstruct the entire absolute abundance vector from a 32-d per-protein code
at every timepoint, and most of that vector is the unchanged baseline.

Change: predict the *change* and add it to the observed baseline,
``P(t) = P(0) + g(z(t))``. The model now starts at the no-change baseline by
construction and spends its capacity on the perturbation response, which is the
quantity of interest. Skill can no longer be worse than the baseline unless the
learned delta is actively harmful.

**B. Real elapsed time.** `real_time=True`

The released code integrates over ``torch.linspace(0, 3, 4)`` — ticks 0, 1, 2, 3
— so the solver is never told that the 24->48 h gap is four times the 0->6 h gap.
The SI's own justification for the four-timepoint design is a Taylor expansion
``x(t) = x(0) + tv + (t^2/2)a + O(t^3)`` over the *real* 6/24/48 h grid, where
those terms differ by orders of magnitude between the three points.

Change: integrate over ``[0, 6, 24, 48] / 48``. Rescaling to O(1) keeps the
fixed-step solve stable; the *ratios* between intervals are what matters and
they are now correct. Sub-steps per interval are configurable because the 24->48
leap is half the whole horizon.

**C. Non-autonomous, perturbation-conditioned vector field.**
`time_conditioned=True`, `pert_conditioned=True`

The SI writes the dynamics as ``f_theta(z_ij(t), t, D_j)`` and the integral
relation as ``z(t_k) = z(0) + int f_theta(z(s), s, D_j) ds``. The released field
is ``nn.Sequential`` over the latent width alone: its first layer has
``in_features = hidden``, it receives neither ``s`` nor ``D_j``, and torchdyn
warns on construction that the callable lacks a time argument. The perturbation
reaches the dynamics only through the initial state.

Change: implement what the SI describes. The field concatenates a Fourier time
embedding and a per-protein perturbation embedding onto the state.

**D. Low-rank protein coupling.** `coupling_rank=r`

Every convolution in the released dynamics module uses ``kernel_size=1``, so each
protein's trajectory is a function of its own baseline and its own perturbation
entry through a shared map. There is no protein-protein coupling anywhere in
module 1. A model described as learning "protein network dynamics" cannot
represent a network.

Naively widening the kernel is wrong: protein row order is arbitrary, so a
kernel-3 convolution mixes arbitrary neighbours. Instead we use a **low-rank
global coupling**: project the protein axis onto ``r`` learned factors, mix, and
project back. ``U`` and ``V`` are dense over proteins, which is meaningful
because protein identity is fixed across samples (row i is always the same
protein), and costs ``2*N*r`` parameters rather than ``N^2``.

**E. Low-rank phenotype head.** `head_rank=r`

`drugsens_conv1` is ``Conv1d(5585 -> 32, kernel_size=4)``: 714,880 parameters,
**90.5% of the whole model**, mapping every protein into a 32-d feature. On
~1,100 labelled conditions that is a very wide layer.

Change: factor it as ``5585 -> r -> 32``, cutting that layer to ``N*r + r*32*4``.
At r=64 the model drops from 790,306 to roughly 150,000 parameters.
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .official import FullyConnectedLayer, rk4_classical

__all__ = ["ProteinTalksR", "DEFAULT_HOURS"]

DEFAULT_HOURS = (0.0, 6.0, 24.0, 48.0)


# --------------------------------------------------------------------------- #
class _TimeEmbedding(nn.Module):
    """Fourier features of scalar t, so the field can be non-autonomous."""

    def __init__(self, dim: int = 8):
        super().__init__()
        self.dim = dim
        self.register_buffer(
            "freqs", torch.tensor([2.0**k for k in range(dim // 2)]), persistent=False
        )

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        t = t.reshape(1)
        ang = t * self.freqs * math.pi
        return torch.cat([torch.sin(ang), torch.cos(ang)], dim=-1)  # (dim,)


class _LowRankCoupling(nn.Module):
    """Cross-protein mixing through r learned factors.

    h: (B, N, H) -> (B, N, H). Costs 2*N*r parameters instead of N^2.
    """

    def __init__(self, n_proteins: int, rank: int, hidden: int):
        super().__init__()
        self.down = nn.Parameter(torch.randn(n_proteins, rank) / math.sqrt(n_proteins))
        self.up = nn.Parameter(torch.randn(n_proteins, rank) / math.sqrt(rank))
        self.mix = nn.Linear(hidden, hidden)
        self.scale = nn.Parameter(torch.zeros(1))  # start as a no-op

    def forward(self, h: torch.Tensor) -> torch.Tensor:
        # (B, N, H) -> factors (B, r, H) -> back to (B, N, H)
        f = torch.einsum("bnh,nr->brh", h, self.down) / h.shape[1] ** 0.5
        f = self.mix(f)
        return self.scale * torch.einsum("brh,nr->bnh", f, self.up)


class _ConditionedField(nn.Module):
    """The SI's f(z, t, D): state, time and perturbation all enter the field."""

    def __init__(
        self,
        hidden: int,
        dropout: float,
        n_proteins: int,
        time_conditioned: bool,
        pert_conditioned: bool,
        coupling_rank: int,
        time_dim: int = 8,
        pert_dim: int = 8,
    ):
        super().__init__()
        self.time_conditioned = time_conditioned
        self.pert_conditioned = pert_conditioned

        in_dim = hidden
        if time_conditioned:
            self.t_embed = _TimeEmbedding(time_dim)
            in_dim += time_dim
        if pert_conditioned:
            self.pert_proj = nn.Linear(1, pert_dim)
            in_dim += pert_dim

        self.net = nn.Sequential(
            FullyConnectedLayer(in_dim, hidden, nn.Softplus(), dropout),
            FullyConnectedLayer(hidden, hidden, None, dropout),
        )
        self.coupling = (
            _LowRankCoupling(n_proteins, coupling_rank, hidden)
            if coupling_rank > 0
            else None
        )
        self._pert = None  # set per batch by the enclosing module

    def set_perturbation(self, pert: torch.Tensor | None) -> None:
        """pert: (B, N, 1), held fixed across the solve (it is time-invariant)."""
        self._pert = pert

    def forward(self, t: torch.Tensor, h: torch.Tensor) -> torch.Tensor:
        parts = [h]
        if self.time_conditioned:
            te = self.t_embed(t.detach() if torch.is_tensor(t) else torch.tensor(t))
            parts.append(te.expand(h.shape[0], h.shape[1], -1))
        if self.pert_conditioned and self._pert is not None:
            parts.append(self.pert_proj(self._pert))
        z = torch.cat(parts, dim=-1) if len(parts) > 1 else h
        out = self.net(z)
        if self.coupling is not None:
            out = out + self.coupling(h)
        return out


# --------------------------------------------------------------------------- #
class ProteinTalksR(nn.Module):
    """The released ppODE with the five changes above, each independently gated.

    With every flag off and ``coupling_rank=0``, ``head_rank=0``, this is the
    released architecture (see ``official.ppODE``), so the ablation ladder starts
    from the published model rather than from a different one.
    """

    def __init__(
        self,
        pro_feats: int,
        hidden_feats: int = 64,
        drug_feature_feats: int = 935,
        dropout: float = 0.0,
        hours: tuple[float, ...] = DEFAULT_HOURS,
        *,
        delta_decode: bool = True,
        real_time: bool = True,
        time_conditioned: bool = True,
        pert_conditioned: bool = True,
        coupling_rank: int = 8,
        head_rank: int = 64,
        substeps: int = 2,
    ):
        super().__init__()
        self.pro_feats = pro_feats
        self.mid = 32
        self.delta_decode = delta_decode
        self.real_time = real_time
        self.hours = tuple(hours)
        self.substeps = max(1, substeps)

        # ---- module 1: encoder, field, decoder (as released, plus the flags)
        self.linear_input = FullyConnectedLayer(2, self.mid, nn.Softplus(), dropout)
        self.conv1 = nn.Conv1d(self.mid, hidden_feats, kernel_size=1)
        self.conv1_norm = nn.GroupNorm(1, hidden_feats)
        self.conv2 = nn.Conv1d(hidden_feats, self.mid, kernel_size=1)
        self.field = _ConditionedField(
            hidden_feats, dropout, pro_feats,
            time_conditioned, pert_conditioned, coupling_rank,
        )
        self.layer_final = nn.Linear(self.mid, 1)
        if delta_decode:
            # Start as an exact no-change predictor; the model must earn a delta.
            nn.init.zeros_(self.layer_final.weight)
            nn.init.zeros_(self.layer_final.bias)

        # ---- module 2: phenotype head
        n_t = len(self.hours)
        if head_rank > 0:
            self.prot_down = nn.Linear(pro_feats, head_rank, bias=False)
            self.drugsens_conv1 = nn.Conv1d(head_rank, self.mid, kernel_size=n_t)
        else:
            self.prot_down = None
            self.drugsens_conv1 = nn.Conv1d(pro_feats, self.mid, kernel_size=n_t)
        self.drugs_conv2 = nn.Conv1d(drug_feature_feats, self.mid, kernel_size=2)
        self.pheno_fc1 = nn.Linear(self.mid * 2, self.mid)
        self.pheno_fc2 = nn.Linear(self.mid, 1)

    # -- time grid ---------------------------------------------------------- #
    def _time_grid(self, device, dtype):
        if self.real_time:
            # Real hours, rescaled to O(1). The ratios between intervals are what
            # the released integer-tick grid gets wrong; rescaling preserves them.
            base = torch.tensor(self.hours, device=device, dtype=dtype)
            base = base / max(self.hours)
        else:
            base = torch.linspace(0, len(self.hours) - 1, len(self.hours),
                                  device=device, dtype=dtype)
        if self.substeps == 1:
            return base, None
        # Refine each interval, then keep only the observation indices.
        pieces, keep = [base[:1]], [0]
        for i in range(len(base) - 1):
            step = (base[i + 1] - base[i]) / self.substeps
            pieces.append(base[i] + step * torch.arange(
                1, self.substeps + 1, device=device, dtype=dtype))
            keep.append(keep[-1] + self.substeps)
        return torch.cat(pieces), keep

    # -- forward ------------------------------------------------------------ #
    def forward(self, x, pert, fp_phA, fp_phB):
        """x, pert: (B, N, 1); fp_phA/B: (B, 935, 1)."""
        emb = self.linear_input(torch.cat([x, pert], dim=-1))      # (B, N, 32)
        emb = torch.transpose(emb, 1, 2)
        emb = F.relu(self.conv1_norm(self.conv1(emb)))
        emb = torch.transpose(emb, 1, 2)                            # (B, N, H)

        self.field.set_perturbation(pert)
        t, keep = self._time_grid(x.device, x.dtype)
        traj = rk4_classical(self.field, emb, t)
        if keep is not None:
            traj = traj[keep]
        traj = traj[1:]                                             # (T, B, N, H)
        self.field.set_perturbation(None)

        T, B, N, H = traj.shape
        z = self.conv2(traj.reshape(T * B, N, H).transpose(1, 2)).transpose(1, 2)
        out = self.layer_final(z).reshape(T, B, N).permute(1, 0, 2)  # (B, T, N)
        if self.delta_decode:
            out = out + x.squeeze(-1).unsqueeze(1)                   # P(t) = P0 + delta

        xy = torch.cat([x, out.transpose(1, 2)], dim=-1)             # (B, N, T+1)
        if self.prot_down is not None:
            xy = self.prot_down(xy.transpose(1, 2)).transpose(1, 2)  # (B, r, T+1)
        xy = F.relu(self.drugsens_conv1(xy)).squeeze(-1)

        fp = F.relu(self.drugs_conv2(torch.cat([fp_phA, fp_phB], dim=-1))).squeeze(-1)
        h = F.relu(self.pheno_fc1(torch.cat([xy, fp], dim=1)))
        return out, self.pheno_fc2(h).squeeze(1)   # logits; sigmoid in the loss

    def n_parameters(self) -> dict[str, int]:
        dyn = ["linear_input", "conv1", "conv1_norm", "conv2", "field", "layer_final"]
        head = ["prot_down", "drugsens_conv1", "drugs_conv2", "pheno_fc1", "pheno_fc2"]
        g = lambda names: sum(  # noqa: E731
            p.numel() for k, p in self.named_parameters()
            if any(k.startswith(n) for n in names)
        )
        return {"dynamics": g(dyn), "head": g(head),
                "total": sum(p.numel() for p in self.parameters())}
