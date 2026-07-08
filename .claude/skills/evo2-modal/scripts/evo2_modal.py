"""
Evo2-7B genomic foundation model on Modal.

Supports an obesity drug program's gene-layer work:
  - siRNA off-target likelihoods
  - zero-shot variant-effect scoring of INHBE / ACVR1C / GPR75 alleles

Design
------
- Persistent modal.Volume "evo2-weights" mounted at /root/.cache/huggingface so
  the HF weight download (arcinstitute/evo2_7b, ~14 GB bf16) persists across runs.
- Image: NVIDIA CUDA 12.4.1 devel base + Python 3.11, torch 2.6 (cu124),
  flash-attn built for sm_90 (H100 only, to keep the build tractable), then the
  Arc Institute `evo2` package (light install path -> 7B only, no transformer_engine).
- download_weights(): CPU function, snapshot_download into the Volume. Cheap.
- variant_effect(ref, alt): H100 function, loads Evo2('evo2_7b'), scores both
  sequences (mean sequence log-likelihood) and returns the alt-ref delta.

Invoke
------
  modal run scratchpad/evo2_modal.py::download_weights
  modal run scratchpad/evo2_modal.py::variant_effect        # uses built-in demo pair
  modal run scratchpad/evo2_modal.py::score --seq ACGT...   # single-seq likelihood

Zero-shot variant-effect convention: score = LL(alt) - LL(ref).
More-negative delta => the alternate allele is less likely under the genomic
prior (candidate deleterious / disruptive); ~0 => tolerated.
"""

import modal

APP_NAME = "evo2-7b"
HF_CACHE = "/root/.cache/huggingface"
MODEL_REPO = "arcinstitute/evo2_7b"
MODEL_NAME = "evo2_7b"

app = modal.App(APP_NAME)

# Persistent volume for HF weight cache (survives across runs).
weights_volume = modal.Volume.from_name("evo2-weights", create_if_missing=True)

# HF token as a Modal secret, taken from the local environment at deploy time.
hf_secret = modal.Secret.from_local_environ(["HF_TOKEN"])

# ---------------------------------------------------------------------------
# Image  (ACTIVE recipe: NGC base + flash-attn compiled from source)
# ---------------------------------------------------------------------------
# This is the GUARANTEED recipe. The NGC PyTorch container ships torch + cuBLAS +
# Transformer Engine + Apex all built together and matched. In earlier testing this
# NGC stack LOADED evo2_7b on H100 ("Model loaded.") and ran TE's GEMM in the forward
# pass — the ONLY failure was flash-attn: NGC bundles flash_attn 2.4.2, whose
# fwd() lacks the `softcap` arg vortex 1.1.0 passes. So we ONLY replace flash-attn,
# compiling 2.8.0.post2 FROM SOURCE against NGC's exact torch (symbols + C++ ABI +
# softcap all correct by construction). No prebuilt flash-attn 2.8 wheel works on
# NGC's torch 2.6.0a0 alpha (they need a stable-torch c10 symbol it lacks).
#
# COST: a one-time ~1-2h flash-attn source compile (cached thereafter). TE/torch/
# cuBLAS are left exactly as NGC ships them, so there is no TE build and no cuBLAS
# version mismatch.
#
# (An alternative all-prebuilt-wheels recipe on stable torch 2.7 is documented at the
#  bottom of this file as `evo2_image_stable_wheels`. It builds+imports+loads but was
#  not validated end-to-end: the cu126 variant hit a TE<->cuBLAS-12.6 GEMM mismatch,
#  and the cu128 variant deployed but the H100 run hung at first-call triton/JIT of
#  vortex's hyena kernels and was stopped before producing a score.)
evo2_image = (
    modal.Image.from_registry("nvcr.io/nvidia/pytorch:24.12-py3", add_python=None)
    .env(
        {
            "HF_HOME": HF_CACHE,
            "HF_HUB_ENABLE_HF_TRANSFER": "0",
            "TORCH_CUDA_ARCH_LIST": "9.0",  # H100 only -> smaller flash-attn build
            "MAX_JOBS": "8",
            "FLASH_ATTENTION_FORCE_BUILD": "TRUE",
        }
    )
    # Pure-python deps only (--no-deps on evo2/vtx keeps NGC's torch/flash-attn/TE).
    .pip_install(
        "huggingface_hub", "hf_transfer", "biopython", "einops==0.8.1", "rich",
        "ninja", "packaging",
    )
    .pip_install("vtx==1.1.0", extra_options="--no-deps")
    .pip_install("evo2", extra_options="--no-deps")
    # Compile flash-attn 2.8.0.post2 against NGC's exact torch (--no-binary forces the
    # source build; --no-deps keeps torch untouched). Fixes the vortex `softcap` fwd().
    .pip_install(
        "flash-attn==2.8.0.post2",
        extra_options="--no-deps --no-build-isolation --no-binary flash-attn",
    )
)


# ---------------------------------------------------------------------------
# Weight download (CPU, cheap)
# ---------------------------------------------------------------------------
@app.function(
    image=evo2_image,
    volumes={HF_CACHE: weights_volume},
    secrets=[hf_secret],
    timeout=60 * 60,
)
def download_weights():
    """Snapshot-download the Evo2-7B repo into the persistent Volume."""
    import os
    from huggingface_hub import snapshot_download

    print(f"Downloading {MODEL_REPO} into {HF_CACHE} ...")
    path = snapshot_download(
        repo_id=MODEL_REPO,
        token=os.environ.get("HF_TOKEN"),
        # allow all files; Evo2 loader picks the checkpoint it needs
    )
    print("snapshot_download complete:", path)

    # List what we got.
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            try:
                total += os.path.getsize(fp)
            except OSError:
                pass
    print(f"Cached files under {path}, total ~{total / 1e9:.1f} GB")

    weights_volume.commit()
    print("Volume committed.")
    return {"path": path, "gb": round(total / 1e9, 1)}


# ---------------------------------------------------------------------------
# Scoring helpers (GPU)
# ---------------------------------------------------------------------------
def _load_model():
    from evo2 import Evo2

    print(f"Loading Evo2('{MODEL_NAME}') ...")
    model = Evo2(MODEL_NAME)
    print("Model loaded.")
    return model


@app.function(
    image=evo2_image,
    gpu="H100",
    volumes={HF_CACHE: weights_volume},
    secrets=[hf_secret],
    timeout=60 * 120,  # generous: vortex triton hyena kernels JIT/autotune on 1st call
)
def score(seqs: list[str] | None = None):
    """Return per-sequence mean log-likelihoods under Evo2-7B."""
    if seqs is None:
        # tiny default so `modal run ...::score` works with no args
        seqs = ["ACGT" * 25]
    if isinstance(seqs, str):
        seqs = [seqs]

    model = _load_model()
    scores = model.score_sequences(seqs)
    scores = [float(s) for s in scores]
    for i, (s, sc) in enumerate(zip(seqs, scores)):
        print(f"seq[{i}] len={len(s)} logL_mean={sc:.6f}")
    return scores


@app.function(
    image=evo2_image,
    gpu="H100",
    volumes={HF_CACHE: weights_volume},
    secrets=[hf_secret],
    timeout=60 * 120,  # generous: vortex triton hyena kernels JIT/autotune on 1st call
)
def variant_effect(ref_seq: str | None = None, alt_seq: str | None = None):
    """
    Zero-shot variant-effect score = LL(alt) - LL(ref).

    If ref/alt are not supplied, a short human demo pair is used (a ~249 bp
    window with a single-nucleotide substitution) so the GPU path can be
    validated end to end in a couple of minutes.
    """
    if ref_seq is None or alt_seq is None:
        # ~249 bp window of human sequence; alt differs by one SNV (pos 124: G->A).
        ref_seq = (
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT"
            "CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC"
            "GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT"
            "ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA"
            "ACAATTGAAT"
        )
        # single-nucleotide substitution at index 124
        alt_list = list(ref_seq)
        alt_list[124] = "A" if alt_list[124] != "A" else "T"
        alt_seq = "".join(alt_list)
        print(f"Using demo pair; ref len={len(ref_seq)}, single SNV at idx 124.")

    model = _load_model()
    ref_ll, alt_ll = (float(x) for x in model.score_sequences([ref_seq, alt_seq]))
    delta = alt_ll - ref_ll

    result = {
        "ref_logL_mean": ref_ll,
        "alt_logL_mean": alt_ll,
        "delta_alt_minus_ref": delta,
        "ref_len": len(ref_seq),
        "alt_len": len(alt_seq),
    }
    print("VARIANT-EFFECT RESULT:", result)
    return result


@app.function(image=evo2_image, timeout=600)
def debug_import():
    """CPU-only: verify the full evo2 import chain (flash_attn / TE / vortex).
    No GPU cost. Model construction still needs a GPU, so we only import here."""
    import torch
    print("torch.__version__ =", torch.__version__)
    print("torch _GLIBCXX_USE_CXX11_ABI =", torch._C._GLIBCXX_USE_CXX11_ABI)
    import flash_attn_2_cuda  # noqa  -> the exact import that was failing
    print("flash_attn_2_cuda import: OK")
    import transformer_engine.pytorch  # noqa  -> TE torch binding
    print("transformer_engine.pytorch import: OK")
    from evo2 import Evo2  # noqa  -> full vortex import chain
    print("from evo2 import Evo2: OK")
    return "import-ok"


@app.local_entrypoint()
def main():
    """`modal run scratchpad/evo2_modal.py` -> download then validate."""
    print("Step 1: download weights (CPU)")
    print(download_weights.remote())
    print("Step 2: variant-effect validation (H100)")
    print(variant_effect.remote())


# ===========================================================================
# ALTERNATIVE RECIPE (all-prebuilt-wheels on STABLE torch 2.7) — NOT validated
# ===========================================================================
# Avoids the ~1-2h flash-attn source compile of the active recipe: every
# torch-adjacent piece is a prebuilt wheel on stable torch 2.7, EXCEPT the small
# transformer_engine_torch binding (sdist -> compiles with clang, ~40 min once).
# Status: builds + deploys cleanly; on H100 flash_attn / transformer_engine.pytorch
# import and Evo2('evo2_7b') loads. Two variants were tried but NEITHER produced a
# score end-to-end:
#   * cu126 (cuBLAS 12.6): TE GEMM raised "cuBLAS Error: an unsupported value or
#     parameter" — TE 2.5's prebuilt core targets a newer cuBLASLt.
#   * cu128 (cuBLAS 12.8, shown below — the intended fix): deployed, but the first
#     H100 run HUNG at first-call JIT/autotune of vortex's triton hyena kernels
#     (never reached the GEMM) and was stopped under the H100 budget. If you retry,
#     give variant_effect a much larger `timeout` (e.g. 2-3h) for the one-time
#     triton warmup, or pre-warm on a tiny sequence.
#
# evo2_image_stable_wheels = (
#     modal.Image.from_registry(
#         "nvidia/cuda:12.8.1-cudnn-devel-ubuntu22.04", add_python="3.12")
#     .apt_install("git", "build-essential", "ninja-build", "cmake")
#     .env({"HF_HOME": HF_CACHE, "HF_HUB_ENABLE_HF_TRANSFER": "0",
#           "CUDA_HOME": "/usr/local/cuda", "TORCH_CUDA_ARCH_LIST": "9.0",
#           "MAX_JOBS": "8", "NVTE_FRAMEWORK": "pytorch"})
#     .pip_install("torch==2.7.0", index_url="https://download.pytorch.org/whl/cu128")
#     .pip_install("packaging", "ninja", "cmake", "pybind11", "wheel", "setuptools",
#                  "huggingface_hub", "hf_transfer", "biopython", "einops==0.8.1", "rich")
#     .pip_install(
#         "flash-attn @ https://github.com/Dao-AILab/flash-attention/releases/download/"
#         "v2.8.0.post2/flash_attn-2.8.0.post2+cu12torch2.7cxx11abiTRUE-cp312-cp312-"
#         "linux_x86_64.whl", extra_options="--no-deps --no-build-isolation")
#     .apt_install("clang")  # transformer_engine_torch build invokes clang++
#     .pip_install("transformer_engine[pytorch]==2.5.0", extra_options="--no-build-isolation")
#     .pip_install("vtx==1.1.0", extra_options="--no-deps")
#     .pip_install("evo2", extra_options="--no-deps")
# )
