#!/usr/bin/env python
"""
Raygun 8.8M (species/function-aware) demo on CPU.

Template: mCherry (236 aa). We encode it with ESM-2 650M, then decode variants
at several target lengths / noise levels, and align each variant back to the
template with a gapped Needleman-Wunsch (BLOSUM62) to quantify the real
insertions / deletions / substitutions. Results are dumped to JSON for the
interactive viewer.
"""
import time, json, torch
from raygun.pretrained import raygun_8_8mil_800M
from esm.pretrained import esm2_t33_650M_UR50D
from Bio.Align import PairwiseAligner, substitution_matrices

DEVICE = "cpu"
torch.manual_seed(7)

TEMPLATE = ("MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQF"
            "MYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPV"
            "MQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHN"
            "EDYTIVEQYERAEGRHSTGGMDELYK")

# (label, target_length, noise)
JOBS = [
    ("v1_same_lownoise",  236, 0.10),
    ("v2_same_highnoise", 236, 0.30),
    ("v3_shortened",      210, 0.15),
    ("v4_extended",       260, 0.15),
]

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = "global"


def analyze(ref, var):
    aln = aligner.align(ref, var)[0]
    a, b = str(aln[0]), str(aln[1])              # gapped strings, equal length
    cols = len(a)
    matches = mism = gaps = pos = 0              # pos = positive BLOSUM (similar)
    bl = aligner.substitution_matrix
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matches += 1
            pos += 1
        else:
            mism += 1
            try:
                if bl[x, y] > 0:
                    pos += 1
            except Exception:
                pass
    aligned = matches + mism
    return {
        "ref_aln": a, "var_aln": b,
        "columns": cols,
        "matches": matches, "mismatches": mism, "gap_cols": gaps,
        "pct_identity": round(100 * matches / aligned, 1) if aligned else 0.0,
        "pct_similarity": round(100 * pos / aligned, 1) if aligned else 0.0,
        "score": round(float(aln.score), 1),
    }


def main():
    t0 = time.time()
    print("Loading ESM-2 650M ...", flush=True)
    esm, alphabet = esm2_t33_650M_UR50D()
    esm = esm.to(DEVICE).eval()
    bc = alphabet.get_batch_converter()

    print("Loading Raygun 8.8M (species/function-aware) ...", flush=True)
    ray = raygun_8_8mil_800M().to(DEVICE).eval()
    print(f"  ready ({time.time()-t0:.1f}s)\n", flush=True)

    _, _, toks = bc([("mcherry", TEMPLATE)])
    with torch.no_grad():
        emb = esm(toks.to(DEVICE), repr_layers=[33])["representations"][33][:, 1:-1, :]

    results = {"template": {"name": "mCherry", "seq": TEMPLATE, "len": len(TEMPLATE)},
               "model": "raygun_8_8mil_800M", "variants": []}

    for label, tl, noise in JOBS:
        with torch.no_grad():
            out = ray(emb, target_lengths=torch.tensor([tl], dtype=int),
                      noise=noise, return_logits_and_seqs=True)
        seq = out["generated-sequences"][0]
        stats = analyze(TEMPLATE, seq)
        rec = {"label": label, "target_len": tl, "noise": noise,
               "seq": seq, "len": len(seq), **stats}
        results["variants"].append(rec)
        print(f"[{label}] target={tl} noise={noise} -> len={len(seq)} "
              f"id={stats['pct_identity']}% sim={stats['pct_similarity']}% "
              f"gaps={stats['gap_cols']} subs={stats['mismatches']}", flush=True)

    with open("/home/user/ClaudeCode/scratchpad/raygun_8_8M_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved JSON. Done in {time.time()-t0:.1f}s (CPU).")


if __name__ == "__main__":
    main()
