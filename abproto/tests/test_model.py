"""Tests for the continuous interpolant and the constrained CTMC."""

from __future__ import annotations

import numpy as np

from abproto.model.ctmc import (
    AA, MASK, ConstrainedCTMC, MotifConstraint, seq_to_tokens, tokens_to_seq,
    uniform_denoiser,
)
from abproto.model.interpolants import (
    LinearInterpolant, flow_matching_loss, sample_ode, structure_first_time,
)


# ---- continuous interpolant ---------------------------------------------
def test_interpolant_endpoints():
    interp = LinearInterpolant(latent_scale=0.0)
    x0 = np.random.default_rng(0).normal(size=(10, 3))
    x1 = np.random.default_rng(1).normal(size=(10, 3))
    assert np.allclose(interp.interpolate(x0, x1, 0.0), x0)
    assert np.allclose(interp.interpolate(x0, x1, 1.0), x1)
    assert np.allclose(interp.target_velocity(x0, x1, 0.5), x1 - x0)


def test_flow_loss_zero_when_velocity_matches_target():
    interp = LinearInterpolant(latent_scale=0.0)
    x0 = np.random.default_rng(2).normal(size=(8, 3))
    x1 = np.random.default_rng(3).normal(size=(8, 3))

    def oracle(xt, t, cond):
        return x1 - x0            # exact conditional velocity

    loss = flow_matching_loss(oracle, interp, x0, x1, t=0.4)
    assert loss < 1e-20


def test_ode_transports_noise_to_data():
    # Integrating the conditional velocity from x0 must land exactly on x1.
    x0 = np.random.default_rng(4).normal(size=(12, 3))
    x1 = np.random.default_rng(5).normal(size=(12, 3))

    def oracle(xt, t, cond):
        return x1 - x0

    out = sample_ode(oracle, x0, steps=50)
    assert np.allclose(out, x1, atol=1e-6)


def test_guidance_pulls_toward_target():
    target = np.ones((5, 3)) * 3.0
    x0 = np.zeros((5, 3))

    def zero_v(xt, t, cond):
        return np.zeros_like(xt)

    def grad_U(xt, t, cond):
        return xt - target        # U = 0.5||x-target||^2

    out = sample_ode(zero_v, x0, steps=200, guidance=grad_U, guidance_scale=1.0)
    assert np.linalg.norm(out - target) < np.linalg.norm(x0 - target)


def test_structure_first_time_lags():
    for t in np.linspace(0, 1, 11):
        ts = structure_first_time(float(t), delta=0.2)
        assert ts <= t + 1e-9
    assert structure_first_time(0.1, 0.2) == 0.0
    assert structure_first_time(1.0, 0.2) == 1.0


# ---- constrained CTMC ----------------------------------------------------
def _nglyc_count(seq: str) -> int:
    import re
    return len(re.findall(r"N[^P][ST]", seq))


def _biased_denoiser(tokens, t, cond=None):
    """Posterior biased toward N,S,T,G,D,C so motifs are LIKELY without constraints."""
    L = len(tokens)
    p = np.full((L, len(AA)), 0.2 / len(AA))
    for a in "NSTGDC":
        p[:, AA.index(a)] += 0.8 / 6
    p /= p.sum(axis=1, keepdims=True)
    return p


def test_constrained_ctmc_yields_zero_liabilities():
    L = 40
    region = np.ones(L, dtype=bool)
    sampler = ConstrainedCTMC(denoiser=_biased_denoiser, steps=60)
    for seed in range(25):
        toks = sampler.sample(L, region, seed=seed)
        assert MASK not in toks                      # all positions decided
        seq = tokens_to_seq(toks)
        assert _nglyc_count(seq) == 0                # exact: never a sequon
        assert "C" not in seq                        # no cysteine in region


def test_unconstrained_ctmc_does_create_motifs():
    # With the rules removed, the same biased denoiser DOES produce liabilities,
    # proving the constraint is what removes them (not the denoiser).
    L = 40
    region = np.ones(L, dtype=bool)
    sampler = ConstrainedCTMC(denoiser=_biased_denoiser,
                              constraint=MotifConstraint(rules=[]), steps=60)
    total = 0
    for seed in range(25):
        seq = tokens_to_seq(sampler.sample(L, region, seed=seed))
        total += _nglyc_count(seq) + seq.count("C")
    assert total > 0


def test_framework_positions_preserved():
    L = 20
    region = np.zeros(L, dtype=bool)
    region[5:10] = True                              # only 5 design positions
    init = seq_to_tokens("ACDEFGHIKLMNPQRSTVWY")      # a fixed "framework"
    sampler = ConstrainedCTMC(denoiser=uniform_denoiser, steps=40)
    toks = sampler.sample(L, region, seed=0, init_tokens=init)
    # Positions outside the region are untouched.
    for i in range(L):
        if not region[i]:
            assert toks[i] == init[i]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all model tests passed")
