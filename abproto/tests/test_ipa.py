"""SE(3)-equivariant IPA backbone: the defining symmetry tests + training.

Skipped if torch is absent.
"""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from abproto.model.ipa import frames_from_backbone, frames_from_ca, to_global, to_local
from abproto.model.net import CoDesignModel, IPABackbone, IPABackboneBB, N_AA


def _random_rotation(seed=0):
    g = torch.Generator().manual_seed(seed)
    a = torch.randn(3, 3, generator=g)
    q, r = torch.linalg.qr(a)
    q = q * torch.sign(torch.diag(r))            # make it proper
    if torch.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def test_frames_are_equivariant():
    ca = torch.randn(2, 10, 3)
    Rg = _random_rotation(1)
    tg = torch.randn(3)
    R, t = frames_from_ca(ca)
    R2, t2 = frames_from_ca(ca @ Rg.T + tg)
    # R_i -> Rg R_i and t_i -> Rg t_i + tg
    assert torch.allclose(R2, Rg @ R, atol=1e-4)
    assert torch.allclose(t2, t @ Rg.T + tg, atol=1e-4)


def test_to_global_local_roundtrip():
    ca = torch.randn(1, 6, 3)
    R, t = frames_from_ca(ca)
    p = torch.randn(1, 6, 4, 3)
    assert torch.allclose(to_local(R, t, to_global(R, t, p)), p, atol=1e-5)


def _inputs(B=2, Ld=6, La=8, seed=0):
    g = torch.Generator().manual_seed(seed)
    return dict(
        design_coords=torch.randn(B, Ld, 3, generator=g),
        design_tokens=torch.randint(0, N_AA, (B, Ld), generator=g),
        antigen_coords=torch.randn(B, La, 3, generator=g),
        antigen_tokens=torch.randint(0, N_AA, (B, La), generator=g),
        t=torch.rand(B, generator=g),
    )


def test_ipa_velocity_equivariant_and_logits_invariant():
    torch.manual_seed(0)
    net = IPABackbone(d_model=48, n_layers=2).eval()
    x = _inputs()
    Rg = _random_rotation(2)
    tg = torch.randn(3)

    with torch.no_grad():
        v0, logits0 = net(**x)
        xr = dict(x)
        xr["design_coords"] = x["design_coords"] @ Rg.T + tg
        xr["antigen_coords"] = x["antigen_coords"] @ Rg.T + tg
        v1, logits1 = net(**xr)

    # velocity (a displacement) rotates with Rg; translation tg does not affect it
    assert torch.allclose(v1, v0 @ Rg.T, atol=1e-4), (v1 - v0 @ Rg.T).abs().max().item()
    # sequence logits are invariant
    assert torch.allclose(logits1, logits0, atol=1e-4)


def test_backbone_frames_equivariant_and_orthonormal():
    rng = torch.Generator().manual_seed(3)
    N = torch.randn(2, 7, 3, generator=rng)
    CA = torch.randn(2, 7, 3, generator=rng)
    C = torch.randn(2, 7, 3, generator=rng)
    R, t = frames_from_backbone(N, CA, C)
    # proper orthonormal frames
    eye = torch.eye(3)
    assert torch.allclose(R.transpose(-1, -2) @ R, eye.expand_as(R), atol=1e-4)
    assert torch.allclose(torch.det(R), torch.ones(2, 7), atol=1e-4)
    assert torch.allclose(t, CA)
    # equivariance under proper rigid motion
    Rg = _random_rotation(5); tg = torch.randn(3)
    R2, t2 = frames_from_backbone(N @ Rg.T + tg, CA @ Rg.T + tg, C @ Rg.T + tg)
    assert torch.allclose(R2, Rg @ R, atol=1e-4)
    assert torch.allclose(t2, CA @ Rg.T + tg, atol=1e-4)


def test_backbone_frames_handle_collinear_ca():
    # Three residues whose CAs are collinear (x-axis). CA-only frames degenerate
    # here (e2 -> 0); real N-CA-C frames stay orthonormal.
    CA = torch.tensor([[[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]]])     # (1,3,3)
    N = CA + torch.tensor([0.0, 1.46, 0.0])
    C = CA + torch.tensor([0.55, -1.42, 0.3])
    R, _ = frames_from_backbone(N, CA, C)
    eye = torch.eye(3)
    assert torch.allclose(R.transpose(-1, -2) @ R, eye.expand_as(R), atol=1e-4)
    # CA-only frames on the same collinear CAs are NOT orthonormal (degenerate)
    Rca, _ = frames_from_ca(CA)
    err = (Rca.transpose(-1, -2) @ Rca - eye).abs().max()
    assert err > 0.1


def test_bb_model_velocity_equivariant():
    torch.manual_seed(0)
    net = IPABackboneBB(d_model=48, n_layers=2).eval()
    g = torch.Generator().manual_seed(1)
    design = torch.randn(2, 6, 3, 3, generator=g)
    antigen = torch.randn(2, 8, 3, 3, generator=g)
    dtok = torch.randint(0, N_AA, (2, 6), generator=g)
    atok = torch.randint(0, N_AA, (2, 8), generator=g)
    t = torch.rand(2, generator=g)
    Rg = _random_rotation(7); tg = torch.randn(3)

    with torch.no_grad():
        v0, lg0 = net(design, dtok, antigen, atok, t)
        v1, lg1 = net(design @ Rg.T + tg, dtok, antigen @ Rg.T + tg, atok, t)
    # each atom's velocity (a displacement) rotates with Rg
    assert torch.allclose(v1, v0 @ Rg.T, atol=1e-4)
    assert torch.allclose(lg1, lg0, atol=1e-4)


def test_ipa_model_trains():
    torch.manual_seed(0)
    model = CoDesignModel(d_model=48, n_layers=2, lr=2e-3, backbone="ipa")
    g = torch.Generator().manual_seed(1)
    B, Ld, La = 12, 6, 8
    batch = dict(
        design_coords=torch.randn(1, Ld, 3, generator=g).repeat(B, 1, 1),
        design_tokens=torch.randint(0, N_AA, (1, Ld), generator=g).repeat(B, 1),
        antigen_coords=torch.randn(1, La, 3, generator=g).repeat(B, 1, 1),
        antigen_tokens=torch.randint(0, N_AA, (1, La), generator=g).repeat(B, 1),
    )
    hist = model.fit(batch, steps=200)
    early = np.mean([h["loss"] for h in hist[:20]])
    late = np.mean([h["loss"] for h in hist[-20:]])
    assert late < 0.7 * early, f"{early:.3f} -> {late:.3f}"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all ipa tests passed")
