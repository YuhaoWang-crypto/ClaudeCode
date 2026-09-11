#!/usr/bin/env python
"""Prove the ablation ladder starts from the released architecture.

`scripts/optimize_model.py` claims its first rung is the published model. That
claim is only worth anything if it is checked: an ablation ladder whose baseline
is subtly *not* the thing being improved on attributes its gains to the wrong
change.

This loads the same weights into `official.ppODE` and into
`improved.ProteinTalksR` with every flag off, and compares outputs. It should
print a difference of exactly zero — not merely small — because with the flags
off the two forward passes are the same sequence of operations.
"""
from __future__ import annotations

import os
import sys

import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks.improved import ProteinTalksR  # noqa: E402
from proteintalks.official import ppODE  # noqa: E402

RELEASED_FLAGS = dict(
    delta_decode=False, real_time=False, time_conditioned=False,
    pert_conditioned=False, coupling_rank=0, head_rank=0, substeps=1,
)


def main(n_proteins: int = 120, hidden: int = 64, batch: int = 3) -> int:
    torch.manual_seed(0)
    ref = ppODE(pro_feats=n_proteins, hidden_feats=hidden)
    new = ProteinTalksR(pro_feats=n_proteins, hidden_feats=hidden, **RELEASED_FLAGS)

    # The only naming difference is where the ODE field lives.
    mapped = {}
    for k, v in ref.state_dict().items():
        if k.startswith("ode_func.net."):
            mapped["field.net." + k[len("ode_func.net."):]] = v
        elif k.startswith("convdrug1"):
            continue  # declared but unused in the reference forward pass
        else:
            mapped[k] = v
    missing, unexpected = new.load_state_dict(mapped, strict=False)
    assert not missing, f"unmapped parameters in the improved model: {missing}"
    assert not unexpected, f"weights with nowhere to go: {unexpected}"

    ref.eval()
    new.eval()
    x = torch.rand(batch, n_proteins, 1)
    p = torch.rand(batch, n_proteins, 1)
    a = torch.rand(batch, 935, 1)
    b = torch.rand(batch, 935, 1)
    with torch.no_grad():
        y_ref, ph_ref = ref(x, p, a, b)
        y_new, logit = new(x, p, a, b)
        ph_new = torch.sigmoid(logit)

    d_prot = (y_ref - y_new).abs().max().item()
    d_phen = (ph_ref - ph_new).abs().max().item()
    print(f"proteome  max |released - ladder rung 0| = {d_prot:.3e}")
    print(f"phenotype max |released - ladder rung 0| = {d_phen:.3e}")
    ok = d_prot == 0.0 and d_phen == 0.0
    print("\nRESULT:", "ladder rung 0 is the released architecture"
          if ok else "rung 0 DIFFERS from the released architecture")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
