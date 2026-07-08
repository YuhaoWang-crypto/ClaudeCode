---
name: evo2-modal
description: Deploy the Evo2 genomic foundation model (Arc Institute, evo2_7b) on Modal H100 for zero-shot variant-effect scoring via sequence log-likelihood — the right tool for CODING loss-of-function (missense/nonsense/frameshift) and long-context genomic-LM likelihoods (e.g. siRNA off-target / seed scoring). Use alongside alphagenome-modal (which covers regulatory/splice). STATUS - deployment recipe is developed and the model loads + runs on H100, but final scoring validation is gated on a one-time flash-attn/Transformer-Engine build; treat as WIP until a validated log-likelihood number is produced.
---

# Evo2-7B on Modal — genomic language-model variant scoring (WIP)

Evo2 is a DNA foundation model; zero-shot variant effect = difference in sequence log-likelihood between ref and alt. Best for **coding LoF** and **genomic-LM likelihoods** (siRNA off-target/seed), complementing AlphaGenome's regulatory/splice predictions.

## Status (be honest with the user)
- ✅ Weights staged: 13.8 GB `arcinstitute/evo2_7b` committed to Modal Volume `evo2-weights`.
- ✅ Model loads and runs on H100 (`from evo2 import Evo2; Evo2('evo2_7b')` → "Model loaded"; forward pass executes through StripedHyena + Transformer-Engine).
- 🟡 Final scoring not yet validated — blocked on the flash-attn ↔ Transformer-Engine ↔ cuBLAS version chain. The definitive fix is baked into `scripts/evo2_modal.py`.

## The dependency trap (why this is hard)
Evo2 needs `transformer_engine` + `flash-attn` compiled against a matching torch/cuBLAS. Pitfalls hit and their fixes:
1. cuda-devel + pip torch (old ABI) + prebuilt flash-attn 2.8 wheel → `undefined symbol …__cxx11…` (ABI mismatch).
2. NGC `nvcr.io/nvidia/pytorch:24.12-py3` ships flash-attn 2.4.2 → missing `softcap` arg vortex passes.
3. Prebuilt flash-attn 2.8 wheels → `c10::StorageImpl::throw_data_ptr_access_error` undefined on NGC's alpha torch.
4. cu126 torch + TE 2.5 → `cublasLt` "unsupported value" (TE targets newer cuBLAS).

## Recommended recipe (in scripts/evo2_modal.py)
**Definitive (guaranteed match):** base `nvcr.io/nvidia/pytorch:24.12-py3` (ships torch + cuBLAS + TE + Apex matched) + install evo2/vtx `--no-deps` + **compile flash-attn 2.8.0.post2 from source inside the image** (`--no-deps --no-binary flash-attn`, arch sm_90, MAX_JOBS=8). One-time ~1–2h build, cached after; leaves NGC's torch/TE untouched so the cuBLAS GEMM works. GPU function `timeout` ≥ 2h for the first triton autotune warmup.
**Faster alt (untested to completion):** cu128 torch (cuBLAS 12.8 matches TE 2.5) + cu12torch2.7 flash-attn wheel — builds/deploys but a validation run hung at triton warmup; needs `PYTHONUNBUFFERED=1` + a long timeout.

## How to run
```bash
modal run scripts/evo2_modal.py::download_weights   # CPU; weights already staged → no-op
modal deploy scripts/evo2_modal.py                  # triggers the one-time flash-attn source build
modal run scripts/evo2_modal.py::variant_effect     # H100 → {ref_logL, alt_logL, delta_alt_minus_ref}
```
Requires: `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET`, `HF_TOKEN`, H100 (evo2_7b uses FP8/Transformer-Engine — H100-class only). Cost: builds are CPU (~free); one short H100 validation run is a few $.

## To finish it
Run the definitive recipe (accept the ~1–2h one-time flash-attn compile), then confirm `variant_effect` returns a real delta. Until then, use `alphagenome-modal` for the regulatory/splice half of any variant analysis.
