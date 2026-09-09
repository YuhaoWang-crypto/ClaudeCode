"""Faithful re-implementation of the ProteinTalks architecture.

Reference
---------
Sun R., Qian L., Li Y., et al.
"A perturbation proteomics-based foundation model for virtual cell construction"
bioRxiv 2025.02.07.637070 (published in Nature as "An operational perturbation
proteomics-based virtual cell model", doi:10.1038/s41586-026-11001-9).

Every layer below is annotated with the numbered equation from the paper's
"Model architecture and parameters" section that it implements.  Where the text
is dimensionally ambiguous, the chosen interpretation is flagged with a
``# INTERPRETATION`` comment and justified in ``docs/ARCHITECTURE_NOTES.md``.

IMPORTANT
---------
This module implements the paper's **equations**.  The authors' released code
(``guomics-lab/PTV-1``) differs from those equations in ten identified places,
several of which change what the model can represent -- most importantly, its
dynamics module has no protein-protein coupling at all, and its ODE integrates
over uniform integer ticks rather than over 0/6/24/48 hours.

``proteintalks/official.py`` is an exact port of that released code, verified
numerically equivalent to the reference implementation (max difference 1.8e-7)
using the released trained checkpoint.  **Use ``official.py`` if you want the
published model; use this module if you want the published equations.**  See
``docs/REPRODUCTION_STATUS.md`` section 4 for the full list of differences.

Module 1 (proteome dynamics, Eq. 1-6)
    [P0, D]  -> L1 (2 -> 32 per protein)          Eq. 1
             -> C1 (32 -> 128 channels)
             -> f = L3(Dropout(L2(.))) with Softplus, dropout 0.1   Eq. 2-3
             -> ODESolve(f, rk4) over t = 0, 6, 24, 48 h            Eq. 4
             -> C2 (128 -> 32 channels)
             -> L4 (32 -> 1 per protein)                            Eq. 5
             -> MSE against measured proteomes                      Eq. 6

Module 2 (drug efficacy / synergy, Eq. 7-11)
    [P0, P6, P24, P48] -> C3 -> 128-d vector                        Eq. 7
    [D1, D2] (935 x 2) -> C4 -> 128-d vector                        Eq. 8
    concat -> L5 (ReLU, -> 32)                                      Eq. 9
           -> L6 (sigmoid, -> 1)                                    Eq. 10
           -> BCE against efficacy / synergy label                  Eq. 11

Total loss (Eq. 12):  Loss = (1 - lambda) * Loss1 + lambda * Loss2, lambda = 0.8.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from torchdiffeq import odeint

__all__ = ["ProteinTalks", "ODEFunc", "ProteomeDynamicsModule", "EfficacyModule"]


# --------------------------------------------------------------------------- #
# Module 1: proteome dynamics
# --------------------------------------------------------------------------- #
class ODEFunc(nn.Module):
    """The vector field f(t, h) parameterised by the two-layer net of Eq. 2-3.

    The paper states: "we have engineered a two-layer linear network, utilizing
    the Softplus as the activation function and implementing a dropout ratio of
    0.1", and "The next step involves parameterizing this two-layer linear
    network with Neural Ordinary Differential Equations".

    So L2 and L3 *are* the ODE vector field, acting on the 128-d per-protein
    channel embedding.  Softplus is used rather than ReLU because a smooth
    vector field is required for a stable fixed-step RK4 solve.
    """

    def __init__(self, hidden_dim: int = 128, dropout: float = 0.1):
        super().__init__()
        self.l2 = nn.Linear(hidden_dim, hidden_dim)  # Eq. 2
        self.l3 = nn.Linear(hidden_dim, hidden_dim)  # Eq. 3
        self.act = nn.Softplus()
        self.dropout = nn.Dropout(dropout)
        # nfe = number of function evaluations, useful for reporting solver cost
        self.nfe = 0

    def forward(self, t: torch.Tensor, h: torch.Tensor) -> torch.Tensor:  # noqa: D401
        """h: (B, n_proteins, hidden_dim) -> dh/dt of the same shape."""
        self.nfe += 1
        z = self.act(self.l2(h))            # Eq. 2
        z = self.act(self.l3(self.dropout(z)))  # Eq. 3
        return z


class ProteomeDynamicsModule(nn.Module):
    """Eq. 1-5: predict P(6h), P(24h), P(48h) from P(0) and the perturbation."""

    def __init__(
        self,
        n_proteins: int,
        embed_dim: int = 32,
        hidden_dim: int = 128,
        dropout: float = 0.1,
        timepoints: tuple[float, ...] = (0.0, 6.0, 24.0, 48.0),
        solver: str = "rk4",
        time_scale: float = 48.0,
    ):
        super().__init__()
        self.n_proteins = n_proteins
        self.timepoints = timepoints
        self.solver = solver
        # INTERPRETATION: the raw times (0..48 h) are rescaled to O(1) before the
        # RK4 solve. Integrating a learned field over 48 unscaled time units with
        # a fixed-step solver diverges in fp32; rescaling is equivalent up to a
        # constant absorbed into the vector field.
        self.time_scale = time_scale

        # Eq. 1: L1 lifts the 2 input features (expression, perturbation) of each
        # protein from 2 -> 32.  Applied identically to all proteins, i.e. a
        # weight-shared per-protein linear map.
        self.l1 = nn.Linear(2, embed_dim)

        # C1: convolutional network raising 32 -> 128 channels along the protein
        # axis.  Conv1d with kernel 3 lets neighbouring proteins in the (fixed,
        # arbitrary) protein ordering mix; this is what the paper describes.
        self.c1 = nn.Sequential(
            nn.Conv1d(embed_dim, hidden_dim, kernel_size=3, padding=1),
            nn.Softplus(),
        )

        # Eq. 2-4: the neural ODE vector field.
        self.ode_func = ODEFunc(hidden_dim, dropout)

        # C2: convolutional network reducing 128 -> 32 channels.
        self.c2 = nn.Sequential(
            nn.Conv1d(hidden_dim, embed_dim, kernel_size=3, padding=1),
            nn.Softplus(),
        )

        # Eq. 5: L4 decodes the 32-d per-protein representation back to a scalar
        # abundance.
        self.l4 = nn.Linear(embed_dim, 1)

    def forward(self, p0: torch.Tensor, d: torch.Tensor) -> torch.Tensor:
        """p0, d: (B, n_proteins) -> (B, n_future_timepoints, n_proteins)."""
        x = torch.stack([p0, d], dim=-1)          # (B, N, 2)
        h = self.l1(x)                            # Eq. 1 -> (B, N, 32)
        h = self.c1(h.transpose(1, 2)).transpose(1, 2)   # -> (B, N, 128)

        t = torch.tensor(
            self.timepoints, dtype=p0.dtype, device=p0.device
        ) / self.time_scale
        self.ode_func.nfe = 0
        # Eq. 4: ODESolve with the fixed-step 4th-order Runge-Kutta solver.
        traj = odeint(self.ode_func, h, t, method=self.solver)  # (T, B, N, 128)
        traj = traj[1:]                            # drop t=0, keep 6/24/48 h

        T, B, N, H = traj.shape
        z = traj.reshape(T * B, N, H).transpose(1, 2)
        z = self.c2(z).transpose(1, 2)             # -> (T*B, N, 32)
        out = self.l4(z).squeeze(-1)               # Eq. 5 -> (T*B, N)
        return out.reshape(T, B, N).permute(1, 0, 2)   # (B, T, N)


# --------------------------------------------------------------------------- #
# Module 2: drug efficacy / synergy
# --------------------------------------------------------------------------- #
class EfficacyModule(nn.Module):
    """Eq. 7-10: predict drug efficacy / combination synergy."""

    def __init__(
        self,
        n_proteins: int,
        n_drug_features: int = 935,
        hidden_dim: int = 128,
        mlp_dim: int = 32,
        n_timepoints: int = 3,
    ):
        super().__init__()
        in_ch = n_timepoints + 1  # [P0, P6h, P24h, P48h]

        # Eq. 7: C3 elevates the stacked proteome tensor to 128 dimensions.
        # INTERPRETATION: the paper gives no pooling step, but L5/L6 must consume
        # a fixed-length vector to emit a scalar.  A strided conv stack followed
        # by global average pooling over the protein axis is the standard way to
        # read "convolutional network -> 128 dimensions" and keeps the parameter
        # count independent of n_proteins.
        self.c3 = nn.Sequential(
            nn.Conv1d(in_ch, 32, kernel_size=7, stride=4, padding=3),
            nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=7, stride=4, padding=3),
            nn.ReLU(),
            nn.Conv1d(64, hidden_dim, kernel_size=7, stride=4, padding=3),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        # Eq. 8: C4 expands the (935 x 2) drug feature tensor to 128 dimensions.
        self.c4 = nn.Sequential(
            nn.Conv1d(2, 32, kernel_size=7, stride=4, padding=3),
            nn.ReLU(),
            nn.Conv1d(32, hidden_dim, kernel_size=7, stride=4, padding=3),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        # Eq. 9-10.
        self.l5 = nn.Sequential(nn.Linear(2 * hidden_dim, mlp_dim), nn.ReLU())
        self.l6 = nn.Linear(mlp_dim, 1)

    def forward(
        self, p0: torch.Tensor, p_future: torch.Tensor, drug_feats: torch.Tensor
    ) -> torch.Tensor:
        """p0: (B, N); p_future: (B, T, N); drug_feats: (B, 2, 935) -> logits (B,)."""
        prot = torch.cat([p0.unsqueeze(1), p_future], dim=1)  # (B, T+1, N)
        z_prot = self.c3(prot)                                # Eq. 7 -> (B, 128)
        z_drug = self.c4(drug_feats)                          # Eq. 8 -> (B, 128)
        z = self.l5(torch.cat([z_prot, z_drug], dim=-1))      # Eq. 9 -> (B, 32)
        return self.l6(z).squeeze(-1)   # Eq. 10 logits; sigmoid applied in loss


# --------------------------------------------------------------------------- #
# Full model
# --------------------------------------------------------------------------- #
class ProteinTalks(nn.Module):
    """The two-module ProteinTalks foundation model.

    Parameters
    ----------
    dynamics_backbone
        ``"neural_ode"``   - the published model (Eq. 4).
        ``"mlp"``          - ablation: replace the ODE solve with a plain
                             residual MLP applied once per timepoint. Isolates
                             how much the continuous-time formulation buys.
        ``"identity"``     - ablation: no dynamics at all, P(t) = P(0).
    """

    def __init__(
        self,
        n_proteins: int,
        n_drug_features: int = 935,
        embed_dim: int = 32,
        hidden_dim: int = 128,
        mlp_dim: int = 32,
        dropout: float = 0.1,
        timepoints: tuple[float, ...] = (0.0, 6.0, 24.0, 48.0),
        dynamics_backbone: str = "neural_ode",
    ):
        super().__init__()
        self.n_proteins = n_proteins
        self.n_future = len(timepoints) - 1
        self.dynamics_backbone = dynamics_backbone

        if dynamics_backbone == "neural_ode":
            self.dynamics = ProteomeDynamicsModule(
                n_proteins, embed_dim, hidden_dim, dropout, timepoints
            )
        elif dynamics_backbone == "mlp":
            self.dynamics = _MLPDynamics(
                n_proteins, embed_dim, hidden_dim, dropout, self.n_future
            )
        elif dynamics_backbone == "identity":
            self.dynamics = _IdentityDynamics(self.n_future)
        else:
            raise ValueError(f"unknown dynamics_backbone: {dynamics_backbone}")

        self.efficacy = EfficacyModule(
            n_proteins, n_drug_features, hidden_dim, mlp_dim, self.n_future
        )

    def forward(self, p0, d, drug_feats):
        p_future = self.dynamics(p0, d)                 # (B, T, N)
        logits = self.efficacy(p0, p_future, drug_feats)  # (B,)
        return p_future, logits

    def n_parameters(self) -> dict[str, int]:
        return {
            "dynamics": sum(p.numel() for p in self.dynamics.parameters()),
            "efficacy": sum(p.numel() for p in self.efficacy.parameters()),
            "total": sum(p.numel() for p in self.parameters()),
        }


# --------------------------------------------------------------------------- #
# Ablation backbones
# --------------------------------------------------------------------------- #
class _MLPDynamics(nn.Module):
    """Discrete-time ablation: same encoder/decoder, no ODE solve."""

    def __init__(self, n_proteins, embed_dim, hidden_dim, dropout, n_future):
        super().__init__()
        self.n_future = n_future
        self.l1 = nn.Linear(2, embed_dim)
        self.c1 = nn.Sequential(
            nn.Conv1d(embed_dim, hidden_dim, 3, padding=1), nn.Softplus()
        )
        self.blocks = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.Softplus(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.Softplus(),
                )
                for _ in range(n_future)
            ]
        )
        self.c2 = nn.Sequential(
            nn.Conv1d(hidden_dim, embed_dim, 3, padding=1), nn.Softplus()
        )
        self.l4 = nn.Linear(embed_dim, 1)

    def forward(self, p0, d):
        x = torch.stack([p0, d], dim=-1)
        h = self.l1(x)
        h = self.c1(h.transpose(1, 2)).transpose(1, 2)
        outs = []
        for blk in self.blocks:
            h = h + blk(h)  # residual step, the Euler analogue of the ODE
            z = self.c2(h.transpose(1, 2)).transpose(1, 2)
            outs.append(self.l4(z).squeeze(-1))
        return torch.stack(outs, dim=1)


class _IdentityDynamics(nn.Module):
    """Null ablation: the proteome never moves."""

    def __init__(self, n_future):
        super().__init__()
        self.n_future = n_future
        self._dummy = nn.Parameter(torch.zeros(1))

    def forward(self, p0, d):
        return p0.unsqueeze(1).repeat(1, self.n_future, 1) + 0.0 * self._dummy
