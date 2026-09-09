#!/usr/bin/env python
"""Validate the port in ``proteintalks/official.py`` against the released weights.

Loads ``best_checkpoint.pth`` from guomics-lab/PTV-1, matches every tensor in the
released state dict to a parameter of our port, and reports any shape mismatch.
A clean run means the port is structurally identical to the trained model.
"""

from __future__ import annotations

import argparse
import os
import sys

import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks.official import ppODE  # noqa: E402

# Reference name -> our name. The only renaming is the ODE field, because the
# reference wraps it in torchdyn's NeuralODE (neuralDE.vf.vf.vf.N) while we call
# torchdiffeq directly (ode_func.net.N).
RENAME = [
    ("neuralDE.vf.vf.vf.", "ode_func.net."),
]


def canon(name: str) -> str:
    for a, b in RENAME:
        if name.startswith(a):
            return b + name[len(a):]
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--proteins", type=int, default=5585)
    ap.add_argument("--hidden", type=int, default=64)
    args = ap.parse_args()

    ckpt = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    sd = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt))
    sd = {k: v for k, v in sd.items() if hasattr(v, "shape")}
    print(f"[checkpoint] {len(sd)} tensors, "
          f"{sum(v.numel() for v in sd.values()):,} elements")
    print(f"[checkpoint] keys present: {sorted(ckpt.keys())[:8]}")

    model = ppODE(pro_feats=args.proteins, hidden_feats=args.hidden)
    ours = dict(model.state_dict())
    print(f"[port]       {len(ours)} tensors, {model.n_parameters()['total']:,} parameters")
    print(f"[port]       phenotype head: {model.n_parameters()['phenotype_head']:,} "
          f"({100*model.n_parameters()['phenotype_head']/model.n_parameters()['total']:.1f}% "
          "of all parameters)")

    matched, mismatched, missing, extra = [], [], [], []
    for k, v in sd.items():
        ck = canon(k)
        if ck in ours:
            (matched if tuple(ours[ck].shape) == tuple(v.shape) else mismatched).append(
                (k, ck, tuple(v.shape), tuple(ours[ck].shape))
            )
        else:
            missing.append((k, tuple(v.shape)))
    seen = {canon(k) for k in sd}
    extra = [(k, tuple(v.shape)) for k, v in ours.items() if k not in seen]

    print(f"\nmatched shapes : {len(matched)}")
    print(f"shape mismatch : {len(mismatched)}")
    print(f"in ckpt only   : {len(missing)}")
    print(f"in port only   : {len(extra)}")

    for tag, items in [("MISMATCH", mismatched), ("CKPT-ONLY", missing), ("PORT-ONLY", extra)]:
        if items:
            print(f"\n[{tag}]")
            for it in items[:20]:
                print("  ", it)

    print("\n[key tensors]")
    for k in ["drugsens_conv1.weight", "drugs_conv2.weight", "pheno_fc1.weight",
              "linear_input.fc.weight", "conv1.weight", "layer_final.weight"]:
        for src, name in [(sd, "ckpt"), (ours, "port")]:
            hit = [kk for kk in src if canon(kk) == k]
            if hit:
                print(f"  {name:5s} {k:26s} {tuple(src[hit[0]].shape)}")

    ok = not mismatched and not missing and not extra
    print("\nRESULT:", "port is structurally identical to the released model"
          if ok else "port differs from the released model (see above)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
