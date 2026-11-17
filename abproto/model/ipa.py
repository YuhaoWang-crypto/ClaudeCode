r"""SE(3)-equivariant backbone: Invariant Point Attention (IPA).

Replaces the plain transformer in net.py. Reference: Jumper et al., AlphaFold2
(2021), Algorithm 22. The construction:

  - Each residue carries a frame T_i = (R_i, x_i), R_i in SO(3), x_i in R^3.
    Here frames are built *equivariantly* from the C-alpha trace by Gram-Schmidt
    on neighbouring CA differences, so the network input stays CA-only.
  - Attention mixes scalar features AND points. Query/key/value points are
    predicted in each residue's LOCAL frame and mapped to the global frame by
    T_i. Pairwise point distances are computed in the global frame, so they are
    INVARIANT to a global rigid motion applied to all frames. Output points are
    mapped back to the local frame, giving INVARIANT node features.
  - A vector readout v_i = R_i @ u_i (u_i a local displacement) is EQUIVARIANT:
    a global rotation R of the complex rotates every v_i by R. That is the
    flow-matching velocity for the coordinates.

Guarantees (verified in tests):
  node features / sequence logits : invariant under global SE(3)
  coordinate velocity             : equivariant under global SE(3)
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


# --------------------------------------------------------------------------
# Equivariant frame construction from a C-alpha trace.
# --------------------------------------------------------------------------
def frames_from_ca(ca: torch.Tensor, eps: float = 1e-6):
    r"""Build per-residue frames (R, t) from CA coordinates ``ca`` (B, L, 3).

    Interior residue i: e1 along (CA_{i+1} - CA_i), e2 orthogonalised from
    (CA_{i-1} - CA_i), e3 = e1 x e2 (Gram-Schmidt). The two end residues copy
    their neighbour's rotation (copying an equivariant frame is equivariant).

    Equivariant under a proper rigid motion ca -> ca @ Rg^T + tg (Rg in SO(3)):
    R_i -> Rg R_i, t_i -> Rg t_i + tg. (Reflections flip e3, as for any frame
    built with a cross product — expected; we only claim SO(3)-equivariance.)
    Requires L >= 3; CA-only frames degenerate for collinear triples, so a
    production model uses N-CA-C frames instead.
    Returns R (B, L, 3, 3), t (B, L, 3).
    """
    B, L, _ = ca.shape
    if L < 3:
        raise ValueError("frames_from_ca needs L >= 3 (use N-CA-C frames otherwise)")

    fwd = ca[:, 2:] - ca[:, 1:-1]          # CA_{i+1}-CA_i for interior i=1..L-2
    bwd = ca[:, :-2] - ca[:, 1:-1]         # CA_{i-1}-CA_i for interior i=1..L-2
    e1 = F.normalize(fwd, dim=-1, eps=eps)
    proj = (bwd * e1).sum(-1, keepdim=True) * e1
    e2 = F.normalize(bwd - proj, dim=-1, eps=eps)
    e3 = torch.cross(e1, e2, dim=-1)
    R_int = torch.stack([e1, e2, e3], dim=-1)      # (B, L-2, 3, 3), columns = basis

    # End residues copy the adjacent interior rotation.
    R = torch.cat([R_int[:, :1], R_int, R_int[:, -1:]], dim=1)
    return R, ca


def rotate(R: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    """Apply R (B,L,3,3) to points p (B,L,...,3) along the last axis."""
    return torch.einsum("blij,bl...j->bl...i", R, p)


def to_global(R, t, p_local):
    """Local points p_local (B,L,...,3) -> global via T_i."""
    extra = p_local.dim() - 3
    t_exp = t[(slice(None), slice(None)) + (None,) * extra + (slice(None),)]
    return rotate(R, p_local) + t_exp


def to_local(R, t, p_global):
    """Global points -> local frame of residue i (inverse of to_global)."""
    extra = p_global.dim() - 3
    t_exp = t[(slice(None), slice(None)) + (None,) * extra + (slice(None),)]
    Rt = R.transpose(-1, -2)
    return rotate(Rt, p_global - t_exp)


# --------------------------------------------------------------------------
# Invariant Point Attention
# --------------------------------------------------------------------------
class InvariantPointAttention(nn.Module):
    def __init__(self, c_s: int, n_heads: int = 4, c_hidden: int = 16,
                 n_qk_points: int = 4, n_v_points: int = 6, c_pair: int = 16):
        super().__init__()
        self.H, self.c = n_heads, c_hidden
        self.Pqk, self.Pv = n_qk_points, n_v_points

        self.to_q = nn.Linear(c_s, n_heads * c_hidden)
        self.to_k = nn.Linear(c_s, n_heads * c_hidden)
        self.to_v = nn.Linear(c_s, n_heads * c_hidden)
        self.to_qp = nn.Linear(c_s, n_heads * n_qk_points * 3)
        self.to_kp = nn.Linear(c_s, n_heads * n_qk_points * 3)
        self.to_vp = nn.Linear(c_s, n_heads * n_v_points * 3)

        # Pair bias from invariant CA-distance RBFs.
        self.n_rbf = c_pair
        self.pair_bias = nn.Linear(c_pair, n_heads)

        self.head_weight = nn.Parameter(torch.zeros(n_heads))     # softplus -> w_h
        out_dim = n_heads * c_hidden + n_heads * n_v_points * 4
        self.to_out = nn.Linear(out_dim, c_s)

    def _rbf(self, t):
        # pairwise CA distances -> gaussian RBF features (B,L,L,n_rbf), invariant
        d = torch.cdist(t, t)                                     # (B,L,L)
        centers = torch.linspace(0.0, 20.0, self.n_rbf, device=t.device)
        return torch.exp(-((d[..., None] - centers) ** 2) / 2.0)

    def forward(self, s, R, t):
        B, L, _ = s.shape
        H, c = self.H, self.c

        q = self.to_q(s).view(B, L, H, c).permute(0, 2, 1, 3)    # (B,H,L,c)
        k = self.to_k(s).view(B, L, H, c).permute(0, 2, 1, 3)
        v = self.to_v(s).view(B, L, H, c).permute(0, 2, 1, 3)
        a_scalar = torch.einsum("bhic,bhjc->bhij", q, k) / (c ** 0.5)

        qp = self.to_qp(s).view(B, L, H, self.Pqk, 3)
        kp = self.to_kp(s).view(B, L, H, self.Pqk, 3)
        vp = self.to_vp(s).view(B, L, H, self.Pv, 3)
        qpg = to_global(R, t, qp).permute(0, 2, 1, 3, 4)         # (B,H,L,P,3)
        kpg = to_global(R, t, kp).permute(0, 2, 1, 3, 4)
        vpg = to_global(R, t, vp).permute(0, 2, 1, 3, 4)         # (B,H,L,Pv,3)

        diff = qpg[:, :, :, None] - kpg[:, :, None]             # (B,H,Li,Lj,P,3)
        sq = (diff ** 2).sum(dim=(-1, -2))                      # (B,H,Li,Lj)
        w = F.softplus(self.head_weight).view(1, H, 1, 1)
        a_point = -0.5 * w * sq

        b_pair = self.pair_bias(self._rbf(t)).permute(0, 3, 1, 2)  # (B,H,L,L)

        logits = (a_scalar + a_point + b_pair) * (1.0 / (3.0 ** 0.5))
        att = torch.softmax(logits, dim=-1)                     # (B,H,Li,Lj)

        o_scalar = torch.einsum("bhij,bhjc->bhic", att, v)      # (B,H,L,c)
        o_scalar = o_scalar.permute(0, 2, 1, 3).reshape(B, L, H * c)

        o_pt_g = torch.einsum("bhij,bhjpx->bhipx", att, vpg)    # (B,H,L,Pv,3)
        o_pt_g = o_pt_g.permute(0, 2, 1, 3, 4)                  # (B,L,H,Pv,3)
        o_pt_local = to_local(R, t, o_pt_g.reshape(B, L, H * self.Pv, 3))
        o_pt_local = o_pt_local.reshape(B, L, H, self.Pv, 3)
        o_norm = o_pt_local.norm(dim=-1)                        # (B,L,H,Pv) invariant

        out = torch.cat([
            o_scalar,
            o_pt_local.reshape(B, L, H * self.Pv * 3),
            o_norm.reshape(B, L, H * self.Pv),
        ], dim=-1)
        return self.to_out(out)                                 # invariant node update


class IPABlock(nn.Module):
    def __init__(self, c_s: int, **kw):
        super().__init__()
        self.ipa = InvariantPointAttention(c_s, **kw)
        self.norm1 = nn.LayerNorm(c_s)
        self.norm2 = nn.LayerNorm(c_s)
        self.ff = nn.Sequential(nn.Linear(c_s, 2 * c_s), nn.GELU(), nn.Linear(2 * c_s, c_s))

    def forward(self, s, R, t):
        s = self.norm1(s + self.ipa(s, R, t))
        s = self.norm2(s + self.ff(s))
        return s
