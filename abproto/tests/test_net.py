"""Torch network: trains (loss decreases) and drives the numpy samplers.

Skipped automatically if torch is not installed.
"""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from abproto.model.ctmc import AA, ConstrainedCTMC, MASK, tokens_to_seq
from abproto.model.interpolants import sample_ode
from abproto.model.net import CoDesignModel, N_AA


def _toy_batch(B=16, Ld=8, La=12, seed=0):
    g = torch.Generator().manual_seed(seed)
    # A single fixed "data" complex, repeated: the net should memorise it.
    design_coords = torch.randn(1, Ld, 3, generator=g).repeat(B, 1, 1)
    design_tokens = torch.randint(0, N_AA, (1, Ld), generator=g).repeat(B, 1)
    antigen_coords = torch.randn(1, La, 3, generator=g).repeat(B, 1, 1)
    antigen_tokens = torch.randint(0, N_AA, (1, La), generator=g).repeat(B, 1)
    return {
        "design_coords": design_coords, "design_tokens": design_tokens,
        "antigen_coords": antigen_coords, "antigen_tokens": antigen_tokens,
    }


def test_training_reduces_loss():
    torch.manual_seed(0)
    model = CoDesignModel(d_model=64, n_layers=2, lr=2e-3)
    batch = _toy_batch()
    hist = model.fit(batch, steps=400)
    early = np.mean([h["loss"] for h in hist[:20]])
    late = np.mean([h["loss"] for h in hist[-20:]])
    # Note: L_flow has an irreducible floor (the velocity target x1-x0 is
    # stochastic given x_t), so we check a clear downward trend, not loss -> 0.
    assert late < 0.65 * early, f"loss did not drop enough: {early:.3f} -> {late:.3f}"
    # sequence head should learn the fixed tokens (CTMC loss well below random CE)
    assert np.mean([h["l_ctmc"] for h in hist[-20:]]) < np.log(N_AA) * 0.7


def test_adapters_drive_numpy_samplers():
    torch.manual_seed(0)
    model = CoDesignModel(d_model=64, n_layers=2)
    Ld, La = 8, 12
    cond = {
        "antigen_coords": np.random.randn(La, 3).astype("float32"),
        "antigen_tokens": np.random.randint(0, N_AA, La).tolist(),
        "design_coords": np.random.randn(Ld, 3).astype("float32"),
        "design_len": Ld,
    }
    # continuous: the velocity field integrates to a finite structure
    vfield = model.velocity_field(cond)
    x0 = np.random.randn(Ld, 3).astype("float32")
    xT = sample_ode(vfield, x0, steps=20)
    assert xT.shape == (Ld, 3) and np.isfinite(xT).all()

    # discrete: the denoiser drives the constrained CTMC; liabilities stay zero
    denoise = model.denoiser(cond)
    sampler = ConstrainedCTMC(denoiser=denoise, steps=40)
    region = np.ones(Ld, dtype=bool)
    toks = sampler.sample(Ld, region, seed=0)
    seq = tokens_to_seq(toks)
    assert MASK not in toks
    assert "C" not in seq                     # hard constraint holds with a learned denoiser


if __name__ == "__main__":
    test_training_reduces_loss(); print("PASS training")
    test_adapters_drive_numpy_samplers(); print("PASS adapters")
