r"""Training loop for the IPA co-design model on real antibody-antigen structures.

Demo (downloads a few real PDB complexes, trains on CPU):
    python -m abproto.train.train --steps 150 --batch-size 4 --out /tmp/abproto_ckpt.pt

Full run: pass a large --pdb-ids list (or wire in the Stage 1 split), use a GPU
(--device cuda), raise --d-model/--n-layers/--steps. The data, batching, loss,
validation and checkpointing are identical; only scale changes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ..data.seeds import SEED_COMPLEXES
from ..data.structures import load_structures
from .dataset import StructureDataset, collate


def _minibatches(ds: StructureDataset, batch_size: int, rng: np.random.Generator):
    idx = rng.permutation(len(ds))
    for i in range(0, len(idx), batch_size):
        yield [ds[j] for j in idx[i : i + batch_size]]


def evaluate(model, ds: StructureDataset, batch_size: int) -> float:
    import torch
    if len(ds) == 0:
        return float("nan")
    model.net.eval()
    losses = []
    with torch.no_grad():
        rng = np.random.default_rng(0)
        for mb in _minibatches(ds, batch_size, rng):
            # reuse training_step's loss without stepping: call net + losses here
            batch = collate(mb)
            losses.append(_loss_only(model, batch))
    return float(np.mean(losses)) if losses else float("nan")


def _loss_only(model, batch) -> float:
    import torch
    import torch.nn.functional as F
    from ..model.interpolants import structure_first_time
    from ..model.ctmc import MASK

    dev = model.device
    x1 = batch["design_coords"].to(dev); seq1 = batch["design_tokens"].to(dev)
    ag_c = batch["antigen_coords"].to(dev); ag_t = batch["antigen_tokens"].to(dev)
    gen = (batch["generate_mask"].to(dev) & batch["design_pad"].to(dev))
    B, Ld, _ = x1.shape
    t = torch.full((B,), 0.5, device=dev)
    x0 = torch.randn_like(x1)
    xt = torch.where(gen[..., None], 0.5 * x1 + 0.5 * x0, x1)
    tseq = structure_first_time(0.5, model.struct_first_delta)
    masked = gen & (torch.rand(B, Ld, device=dev) >= tseq)
    seq_in = torch.where(masked, torch.full_like(seq1, MASK), seq1)
    v, logits = model.net(xt, seq_in, ag_c, ag_t, t,
                          design_pad=batch["design_pad"].to(dev),
                          antigen_pad=batch["antigen_pad"].to(dev))
    lf = ((v - (x1 - x0))[gen] ** 2).mean() if gen.any() else torch.zeros((), device=dev)
    lc = (F.cross_entropy(logits[masked], seq1[masked]) if masked.any()
          else torch.zeros((), device=dev))
    return float(lf + model.lam * lc)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Train IPA co-design model on real structures")
    p.add_argument("--pdb-ids", nargs="*", default=None,
                   help="PDB ids to train on (default: the seed Ab-Ag complexes)")
    p.add_argument("--cache", default="structures_cache")
    p.add_argument("--steps", type=int, default=150)
    p.add_argument("--batch-size", type=int, default=4)
    p.add_argument("--d-model", type=int, default=64)
    p.add_argument("--n-layers", type=int, default=3)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--val-frac", type=float, default=0.2)
    p.add_argument("--device", default="cpu")
    p.add_argument("--out", type=Path, default=Path("abproto_ckpt.pt"))
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)

    import torch
    from ..model.net import CoDesignModel

    ids = args.pdb_ids or SEED_COMPLEXES
    print(f"downloading + parsing {len(ids)} structures ...")
    structs, skipped = load_structures(ids, args.cache)
    print(f"  usable complexes: {len(structs)}  (skipped {len(skipped)})")
    ds = StructureDataset.from_structures(structs)
    print(f"  training examples: {len(ds)}")
    if len(ds) < 2:
        print("not enough examples; pass more --pdb-ids")
        return 1

    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(ds))
    n_val = max(1, int(len(ds) * args.val_frac))
    val = StructureDataset([ds[i] for i in perm[:n_val]])
    train = StructureDataset([ds[i] for i in perm[n_val:]])
    print(f"  train={len(train)}  val={len(val)}")

    torch.manual_seed(args.seed)
    model = CoDesignModel(d_model=args.d_model, n_layers=args.n_layers, lr=args.lr,
                          backbone="ipa", device=args.device)

    step = 0
    hist = []
    while step < args.steps:
        for mb in _minibatches(train, args.batch_size, rng):
            m = model.training_step(collate(mb))
            hist.append(m["loss"]); step += 1
            if step % 25 == 0 or step == 1:
                print(f"step {step:4d}  loss={m['loss']:.4f}  flow={m['l_flow']:.4f}  ctmc={m['l_ctmc']:.4f}")
            if step >= args.steps:
                break

    early = float(np.mean(hist[:10])); late = float(np.mean(hist[-10:]))
    val_loss = evaluate(model, val, args.batch_size)
    print(f"\ntrain loss {early:.3f} -> {late:.3f}   val loss {val_loss:.3f}")

    torch.save({"state_dict": model.net.state_dict(),
                "config": {"d_model": args.d_model, "n_layers": args.n_layers,
                           "backbone": "ipa"}},
               args.out)
    print(f"saved checkpoint -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
