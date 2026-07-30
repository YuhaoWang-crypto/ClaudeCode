#!/usr/bin/env python
"""
Raygun demo (CPU): template-based protein redesign.

Takes an existing protein (mCherry, 236 aa) as a template, encodes it into
Raygun's fixed-length latent with ESM-2 650M embeddings, then decodes back to
NEW sequences at chosen target lengths -- demonstrating length-controlled
insertion/deletion + substitution, the core Raygun operation.
"""
import time, torch
from raygun.pretrained import raygun_2_2mil_800M
from esm.pretrained import esm2_t33_650M_UR50D

DEVICE = "cpu"
torch.manual_seed(0)

# mCherry (red fluorescent protein) -- the design template
TEMPLATE = ("MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQF"
            "MYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPV"
            "MQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHN"
            "EDYTIVEQYERAEGRHSTGGMDELYK")

def pct_identity(a, b):
    """Ungapped identity over the overlapping prefix (quick, position-wise)."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(x == y for x, y in zip(a[:n], b[:n]))
    return 100.0 * same / n

def main():
    print(f"Template: mCherry, length={len(TEMPLATE)} aa\n")

    t0 = time.time()
    print("Loading ESM-2 650M ...", flush=True)
    esm, alphabet = esm2_t33_650M_UR50D()
    esm = esm.to(DEVICE).eval()
    bc = alphabet.get_batch_converter()

    print("Loading Raygun (2.2M) ...", flush=True)
    ray = raygun_2_2mil_800M().to(DEVICE).eval()
    print(f"  models ready ({time.time()-t0:.1f}s)\n", flush=True)

    # ESM-2 layer-33 per-residue embedding of the template (drop <cls>/<eos>)
    _, _, toks = bc([("mcherry", TEMPLATE)])
    with torch.no_grad():
        emb = esm(toks.to(DEVICE), repr_layers=[33])["representations"][33][:, 1:-1, :]
    print(f"ESM embedding shape: {tuple(emb.shape)}  (batch, residues, dim)\n")

    # Generate variants at several target lengths (shrink / keep / expand)
    targets = [200, 236, 270]
    noise = 0.1
    print(f"Generating variants  (noise={noise})")
    print("-" * 72)
    for tl in targets:
        with torch.no_grad():
            out = ray(emb,
                      target_lengths=torch.tensor([tl], dtype=int),
                      noise=noise,
                      return_logits_and_seqs=True)
        seq = out["generated-sequences"][0]
        delta = tl - len(TEMPLATE)
        tag = "same length" if delta == 0 else (f"+{delta} (insertions)" if delta > 0 else f"{delta} (deletions)")
        print(f"target_len={tl:>3}  got={len(seq):>3}  Δ={tag}")
        print(f"  identity-to-template: {pct_identity(seq, TEMPLATE):5.1f}%")
        print(f"  {seq}")
        print("-" * 72)

    print(f"\nDone in {time.time()-t0:.1f}s total (CPU).")

if __name__ == "__main__":
    main()
