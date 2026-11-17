"""SE(3)-equivariant IPA backbone: the defining symmetry tests + training.

Skipped if torch is absent.
"""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from abproto.model.ipa import frames_from_ca, to_global, to_local
from abproto.model.net import CoDesignModel, IPABackbone, N_AA


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
