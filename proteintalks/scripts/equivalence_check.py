#!/usr/bin/env python
"""Numerical equivalence: our port vs the reference implementation itself.

Instantiates the reference ``ppODE`` from a checkout of guomics-lab/PTV-1 and our
``proteintalks.official.ppODE``, loads the *released trained weights* into both,
and compares outputs on identical inputs. The reference integrates with
torchdyn; we use torchdiffeq with the same rk4 solver over the same time grid.

Also prints where the model's parameters actually live, which is the single most
informative number about what this architecture is.
"""

from __future__ import annotations

import argparse
import os
import sys

import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks.official import ppODE as OurPPODE  # noqa: E402


def load_reference(repo_dir):
    sys.path.insert(0, os.path.join(repo_dir, "ProteinTalks"))
    import importlib
    mod = importlib.import_module("model")
    return mod.ppODE


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="checkout of guomics-lab/PTV-1")
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--proteins", type=int, default=5585)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--batch", type=int, default=3)
    args = ap.parse_args()

    torch.manual_seed(0)
    RefPPODE = load_reference(args.repo)

    ref = RefPPODE(node_feats=1, pert_feats=1, hidden_feats=args.hidden, out_feats=1,
                   pro_feats=args.proteins, drug_feature_feats=935, dropout=0.0)
    ours = OurPPODE(pro_feats=args.proteins, hidden_feats=args.hidden, dropout=0.0)

    ckpt = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    sd = ckpt["model_state_dict"]
    ref.load_state_dict(sd, strict=True)
    ours.load_state_dict(
        {k.replace("neuralDE.vf.vf.vf.", "ode_func.net."): v for k, v in sd.items()},
        strict=True,
    )
    ref.eval()
    ours.eval()

    B, N = args.batch, args.proteins
    x = torch.rand(B, N, 1)
    pert = torch.rand(B, N, 1)
    fa = torch.rand(B, 935, 1)
    fb = torch.rand(B, 935, 1)

    with torch.no_grad():
        y_ref, ph_ref, _ = ref(x, pert, fa, fb, "6_24_48")
        y_our, ph_our = ours(x, pert, fa, fb)

    y_ref = y_ref.squeeze(-1)
    if y_ref.shape != y_our.shape:
        y_ref = y_ref.reshape(y_our.shape)

    dy = (y_ref - y_our).abs().max().item()
    dp = (ph_ref - ph_our).abs().max().item()
    print(f"proteome  max |ref - port| = {dy:.3e}   (shape {tuple(y_our.shape)})")
    print(f"phenotype max |ref - port| = {dp:.3e}")
    print("phenotype ref :", [round(v, 6) for v in ph_ref.tolist()])
    print("phenotype port:", [round(v, 6) for v in ph_our.tolist()])
    ok = dy < 1e-4 and dp < 1e-5
    print("\nRESULT:", "numerically equivalent to the reference implementation"
          if ok else "outputs differ")

    print("\n--- where the parameters are ---")
    groups = {
        "module 1: dynamics (neural ODE + encoder/decoder)": [
            "linear_input", "conv1", "conv1_norm", "conv2", "ode_func", "layer_final"],
        "module 2: phenotype head": [
            "drugsens_conv1", "drugs_conv2", "pheno_fc1", "pheno_fc2"],
        "unused in forward": ["convdrug1"],
    }
    total = sum(p.numel() for p in ours.parameters())
    for name, prefixes in groups.items():
        n = sum(p.numel() for k, p in ours.named_parameters()
                if any(k.startswith(pre) for pre in prefixes))
        print(f"  {name:52s} {n:>9,}  ({100*n/total:5.1f}%)")
    print(f"  {'TOTAL':52s} {total:>9,}")
    n_ds = ours.drugsens_conv1.weight.numel()
    print(f"\n  of which drugsens_conv1 alone (32 x {args.proteins} x 4): "
          f"{n_ds:,} ({100*n_ds/total:.1f}%)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
