"""FR-zero-contact guidance: gradient correctness + steering effect."""

from __future__ import annotations

import numpy as np

from abproto.model.guidance import (
    FRContactGuidance, count_fr_contacts, fr_contact_grad, fr_contact_potential,
    framework_mask_from_spans,
)
from abproto.model.interpolants import sample_ode


def test_grad_matches_finite_difference():
    rng = np.random.default_rng(0)
    L, A = 8, 6
    x = rng.normal(scale=4.0, size=(L, 3))
    antigen = rng.normal(scale=4.0, size=(A, 3))
    fr_mask = np.array([True, True, False, True, False, True, True, False])
    g = fr_contact_grad(x, antigen, fr_mask, contact=8.0)

    h = 1e-5
    num = np.zeros_like(x)
    for i in range(L):
        for c in range(3):
            xp = x.copy(); xp[i, c] += h
            xm = x.copy(); xm[i, c] -= h
            num[i, c] = (fr_contact_potential(xp, antigen, fr_mask, 8.0)
                         - fr_contact_potential(xm, antigen, fr_mask, 8.0)) / (2 * h)
    assert np.allclose(g, num, atol=1e-4), np.abs(g - num).max()
    # CDR rows carry no gradient
    assert np.allclose(g[~fr_mask], 0.0)


def test_framework_mask_from_spans():
    m = framework_mask_from_spans(10, {"H1": (2, 4), "H3": (7, 9)})
    assert m.tolist() == [True, True, False, False, True, True, True, False, False, True]


def test_guidance_clears_framework_contacts():
    # antigen cluster at origin; all design residues start on top of it.
    antigen = np.zeros((5, 3)) + np.array([0.0, 0.0, 0.0])
    antigen[:, 0] = np.linspace(-2, 2, 5)
    L = 8
    fr_mask = np.array([True, True, True, True, True, False, False, False])  # 5 FR, 3 CDR
    x0 = np.random.default_rng(1).normal(scale=1.0, size=(L, 3))             # near antigen

    before = count_fr_contacts(x0, antigen, fr_mask, contact=8.0)
    guide = FRContactGuidance(antigen, fr_mask, contact=8.0)

    def zero_v(x, t, cond):
        return np.zeros_like(x)

    out = sample_ode(zero_v, x0, steps=400, guidance=guide, guidance_scale=5.0)
    after = count_fr_contacts(out, antigen, fr_mask, contact=8.0)

    assert before > 0
    assert after < before
    assert guide.potential(out) < 0.05 * guide.potential(x0)   # framework pushed to the 8A rim
    # CDR residues (no gradient) are untouched by pure guidance
    assert np.allclose(out[~fr_mask], x0[~fr_mask], atol=1e-6)


def test_scale_zero_is_inert():
    antigen = np.random.default_rng(2).normal(size=(5, 3))
    L = 6
    fr_mask = np.ones(L, dtype=bool)
    x0 = np.random.default_rng(3).normal(size=(L, 3))

    def zero_v(x, t, cond):
        return np.zeros_like(x)

    out = sample_ode(zero_v, x0, steps=50,
                     guidance=FRContactGuidance(antigen, fr_mask), guidance_scale=0.0)
    assert np.allclose(out, x0)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all guidance tests passed")
