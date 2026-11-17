r"""Trainable shared-backbone network for joint coordinate + sequence co-design.

One transformer backbone, two heads:
  - velocity head  -> v_theta(x_t, t, c)         (continuous flow matching)
  - sequence head  -> p_theta(x_1 | x_t, c)       (discrete CTMC denoiser)

The designed residues (CDRs) and the conditioning antigen residues are embedded
as a single token set; self-attention lets the CDR tokens read the antigen, so
the heads are conditioned on the target. Heads read out only the designed
positions.

This is the last piece before training: the numpy samplers in ``interpolants.py``
and ``ctmc.py`` drive this network through the ``velocity_field`` / ``denoiser``
adapters, so the same orchestration now runs on a learned model.

NOTE on equivariance: this reference backbone featurises raw coordinates and is
*not* SE(3)-equivariant. For production, replace the encoder with frame-based
Invariant Point Attention (IPA) / a structure module; the training objective,
masking path, heads and adapters here are unchanged by that swap.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable, Optional

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    _HAS_TORCH = True
except ImportError:  # pragma: no cover
    _HAS_TORCH = False

from .ctmc import AA, MASK
from .interpolants import structure_first_time

N_AA = len(AA)
VOCAB = N_AA + 1            # + MASK


if _HAS_TORCH:

    def sinusoidal_time_embedding(t: "torch.Tensor", dim: int) -> "torch.Tensor":
        """t: (B,) in [0,1] -> (B, dim) sinusoidal features."""
        half = dim // 2
        freqs = torch.exp(-math.log(10000.0) * torch.arange(half, device=t.device) / half)
        ang = t[:, None] * freqs[None, :] * 1000.0
        emb = torch.cat([torch.sin(ang), torch.cos(ang)], dim=-1)
        if dim % 2:
            emb = F.pad(emb, (0, 1))
        return emb

    class SharedBackboneNet(nn.Module):
        def __init__(self, d_model: int = 128, n_layers: int = 4, n_heads: int = 4,
                     dropout: float = 0.0):
            super().__init__()
            self.d_model = d_model
            self.tok_emb = nn.Embedding(VOCAB, d_model)
            self.chain_emb = nn.Embedding(2, d_model)       # 0 design, 1 antigen
            self.coord_proj = nn.Linear(3, d_model)
            self.time_mlp = nn.Sequential(
                nn.Linear(d_model, d_model), nn.SiLU(), nn.Linear(d_model, d_model))
            layer = nn.TransformerEncoderLayer(
                d_model, n_heads, dim_feedforward=4 * d_model, dropout=dropout,
                batch_first=True, activation="gelu")
            self.encoder = nn.TransformerEncoder(layer, n_layers)
            self.norm = nn.LayerNorm(d_model)
            self.velocity_head = nn.Linear(d_model, 3)
            self.sequence_head = nn.Linear(d_model, N_AA)

        def forward(self, design_coords, design_tokens, antigen_coords, antigen_tokens, t,
                    design_pad=None, antigen_pad=None):
            B, Ld, _ = design_coords.shape
            La = antigen_coords.shape[1]
            src_key_padding = None
            if design_pad is not None or antigen_pad is not None:
                dp = design_pad if design_pad is not None else torch.ones(B, Ld, dtype=torch.bool, device=t.device)
                ap = antigen_pad if antigen_pad is not None else torch.ones(B, La, dtype=torch.bool, device=t.device)
                src_key_padding = ~torch.cat([dp, ap], dim=1)   # True = pad (ignored)

            d_feat = (self.tok_emb(design_tokens) + self.coord_proj(design_coords)
                      + self.chain_emb(torch.zeros(B, Ld, dtype=torch.long, device=t.device)))
            a_feat = (self.tok_emb(antigen_tokens) + self.coord_proj(antigen_coords)
                      + self.chain_emb(torch.ones(B, La, dtype=torch.long, device=t.device)))

            temb = self.time_mlp(sinusoidal_time_embedding(t, self.d_model))[:, None, :]
            x = torch.cat([d_feat, a_feat], dim=1) + temb
            x = self.norm(self.encoder(x, src_key_padding_mask=src_key_padding))
            d_out = x[:, :Ld]
            return self.velocity_head(d_out), self.sequence_head(d_out)

    class IPABackbone(nn.Module):
        """SE(3)-equivariant backbone: frames from CA (Gram-Schmidt) + IPA blocks.

        Velocity output v_i = R_i @ u_i is EQUIVARIANT; sequence logits are
        INVARIANT (both verified in tests). The velocity's rotation-only form
        means it transforms as a displacement (no translation component)."""

        def __init__(self, d_model: int = 128, n_layers: int = 4, n_heads: int = 4):
            super().__init__()
            from .ipa import IPABlock, frames_from_ca, rotate
            self._frames = frames_from_ca
            self._rotate = rotate
            self.d_model = d_model
            self.tok_emb = nn.Embedding(VOCAB, d_model)
            self.chain_emb = nn.Embedding(2, d_model)
            self.time_mlp = nn.Sequential(
                nn.Linear(d_model, d_model), nn.SiLU(), nn.Linear(d_model, d_model))
            self.blocks = nn.ModuleList([IPABlock(d_model, n_heads=n_heads) for _ in range(n_layers)])
            self.vel_local = nn.Linear(d_model, 3)       # local displacement u_i
            self.sequence_head = nn.Linear(d_model, N_AA)

        def forward(self, design_coords, design_tokens, antigen_coords, antigen_tokens, t,
                    design_pad=None, antigen_pad=None):
            B, Ld, _ = design_coords.shape
            La = antigen_coords.shape[1]
            # Per-block frames (a residue's frame depends only on its own chain,
            # so padding / the other block never corrupts it).
            R_d, t_d = self._frames(design_coords)
            R_a, t_a = self._frames(antigen_coords)
            R = torch.cat([R_d, R_a], dim=1)
            tr = torch.cat([t_d, t_a], dim=1)

            tok = torch.cat([design_tokens, antigen_tokens], dim=1)
            chain = torch.cat([torch.zeros(B, Ld, dtype=torch.long, device=t.device),
                               torch.ones(B, La, dtype=torch.long, device=t.device)], dim=1)
            temb = self.time_mlp(sinusoidal_time_embedding(t, self.d_model))[:, None, :]
            s = self.tok_emb(tok) + self.chain_emb(chain) + temb

            kpm = None
            if design_pad is not None or antigen_pad is not None:
                dp = design_pad if design_pad is not None else torch.ones(B, Ld, dtype=torch.bool, device=t.device)
                ap = antigen_pad if antigen_pad is not None else torch.ones(B, La, dtype=torch.bool, device=t.device)
                kpm = torch.cat([dp, ap], dim=1)
            for blk in self.blocks:
                s = blk(s, R, tr, kpm)

            s_d = s[:, :Ld]
            u = self.vel_local(s_d)                      # (B,Ld,3) invariant local vector
            v = self._rotate(R_d, u)                     # (B,Ld,3) equivariant
            return v, self.sequence_head(s_d)


@dataclass
class CoDesignModel:
    """Wraps the network with the training objective and sampler adapters."""

    d_model: int = 128
    n_layers: int = 4
    lam: float = 1.0           # weight on the sequence CTMC loss
    struct_first_delta: float = 0.2
    lr: float = 1e-3
    device: str = "cpu"
    backbone: str = "mlp"      # "mlp" (transformer) | "ipa" (SE(3)-equivariant)

    def __post_init__(self):
        if not _HAS_TORCH:
            raise RuntimeError("PyTorch required: pip install torch")
        if self.backbone == "ipa":
            self.net = IPABackbone(self.d_model, self.n_layers).to(self.device)
        elif self.backbone == "mlp":
            self.net = SharedBackboneNet(self.d_model, self.n_layers).to(self.device)
        else:
            raise ValueError(f"unknown backbone '{self.backbone}'")
        self.opt = torch.optim.Adam(self.net.parameters(), lr=self.lr)

    # -- training --------------------------------------------------------
    def training_step(self, batch: dict) -> dict:
        r"""One optimisation step on L = L_flow + lam * L_ctmc.

        batch: design_coords (B,Ld,3), design_tokens (B,Ld) in [0,20),
               antigen_coords (B,La,3), antigen_tokens (B,La).
        """
        self.net.train()
        x1 = batch["design_coords"].to(self.device)
        seq1 = batch["design_tokens"].to(self.device)
        ag_c = batch["antigen_coords"].to(self.device)
        ag_t = batch["antigen_tokens"].to(self.device)
        B, Ld, _ = x1.shape
        dev = self.device

        # generate_mask: positions being designed (CDRs). Default: all positions.
        gen = batch.get("generate_mask")
        gen = (torch.ones(B, Ld, dtype=torch.bool, device=dev) if gen is None
               else gen.to(dev).bool())
        design_pad = batch.get("design_pad")
        design_pad = (torch.ones(B, Ld, dtype=torch.bool, device=dev) if design_pad is None
                      else design_pad.to(dev).bool())
        antigen_pad = batch.get("antigen_pad")
        antigen_pad = None if antigen_pad is None else antigen_pad.to(dev).bool()
        gen = gen & design_pad                              # never design padding

        t = torch.rand(B, device=dev)
        x0 = torch.randn_like(x1)
        tb = t[:, None, None]
        # Noise ONLY the generated (CDR) positions; framework/antigen stay clean.
        xt = torch.where(gen[..., None], tb * x1 + (1 - tb) * x0, x1)
        target_v = x1 - x0

        # Sequence: mask a fraction of generated positions at the lagged time.
        t_seq = torch.tensor([structure_first_time(float(ti), self.struct_first_delta)
                              for ti in t], device=dev)
        masked = gen & (torch.rand(B, Ld, device=dev) >= t_seq[:, None])
        seq_in = torch.where(masked, torch.full_like(seq1, MASK), seq1)

        v_pred, logits = self.net(xt, seq_in, ag_c, ag_t, t,
                                  design_pad=design_pad, antigen_pad=antigen_pad)

        if gen.any():
            d = (v_pred - target_v)[gen]
            l_flow = (d * d).mean()
        else:
            l_flow = torch.zeros((), device=dev)
        l_ctmc = (F.cross_entropy(logits[masked], seq1[masked]) if masked.any()
                  else torch.zeros((), device=dev))
        loss = l_flow + self.lam * l_ctmc

        self.opt.zero_grad(); loss.backward(); self.opt.step()
        return {"loss": float(loss.detach()), "l_flow": float(l_flow.detach()),
                "l_ctmc": float(l_ctmc.detach())}

    def fit(self, batch: dict, steps: int = 300, log_every: int = 0) -> list[dict]:
        hist = []
        for s in range(steps):
            m = self.training_step(batch)
            hist.append(m)
            if log_every and s % log_every == 0:
                print(f"step {s}: loss={m['loss']:.4f} flow={m['l_flow']:.4f} ctmc={m['l_ctmc']:.4f}")
        return hist

    # -- adapters for the numpy samplers --------------------------------
    def velocity_field(self, cond: dict) -> Callable[[np.ndarray, float, Any], np.ndarray]:
        """Return a VelocityField closure usable by interpolants.sample_ode."""
        ag_c = torch.tensor(cond["antigen_coords"], dtype=torch.float32, device=self.device)[None]
        ag_t = torch.tensor(cond["antigen_tokens"], dtype=torch.long, device=self.device)[None]
        seq = torch.tensor(cond.get("design_tokens_fixed",
                                    [MASK] * cond["design_len"]),
                           dtype=torch.long, device=self.device)[None]

        def vfield(x_np, t, _cond=None):
            self.net.eval()
            with torch.no_grad():
                x = torch.tensor(x_np, dtype=torch.float32, device=self.device)[None]
                tt = torch.tensor([t], dtype=torch.float32, device=self.device)
                v, _ = self.net(x, seq, ag_c, ag_t, tt)
            return v[0].cpu().numpy()
        return vfield

    def denoiser(self, cond: dict) -> Callable[[np.ndarray, float, Any], np.ndarray]:
        """Return a Denoiser closure usable by ctmc.ConstrainedCTMC."""
        ag_c = torch.tensor(cond["antigen_coords"], dtype=torch.float32, device=self.device)[None]
        ag_t = torch.tensor(cond["antigen_tokens"], dtype=torch.long, device=self.device)[None]
        coords = torch.tensor(cond["design_coords"], dtype=torch.float32, device=self.device)[None]

        def denoise(tokens_np, t, _cond=None):
            self.net.eval()
            with torch.no_grad():
                seq = torch.tensor(tokens_np, dtype=torch.long, device=self.device)[None]
                tt = torch.tensor([t], dtype=torch.float32, device=self.device)
                _, logits = self.net(coords, seq, ag_c, ag_t, tt)
                p = torch.softmax(logits[0], dim=-1)
            return p.cpu().numpy()
        return denoise
