"""
evo2_modal.py — Modal adapter that serves Evo2 (evo2_7b) as a scoring service.

WHAT IT DOES
------------
Deploys the Arc Institute Evo2 genome foundation model on a Modal H100 and
exposes ``score(sequences)`` which returns, per sequence, the **mean
log-likelihood** the model assigns to it. Higher (less negative) = the sequence
is more "natural" / consistent with the functional constraints Evo2 learned
across the tree of life. This is the quantity both dashboards flagged as the
production replacement for the divergence-time *proxy* on the fidelity axis:

    G / D  sequence-constraint layer  <-  Evo2 log-likelihood   (this file)

WHY MODAL
---------
Evo2-7B needs an H100 (StripedHyena2 + FP8). This container has no GPU, but the
session ships Modal credentials, so we run the heavy part on Modal's cloud and
pull back only the scalar scores.

RUN
---
Local container must first do (once):
    pip install modal python-socks
    export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt   # trust the agent-proxy CA

Smoke test (deploys image on first call, then scores two toy sequences):
    modal run psmos/evo2_modal.py

From Python (see psmos/evo2_client.py):
    from psmos.evo2_client import score_sequences
    score_sequences(["ACGT...", "ACGT..."])

COST NOTE
---------
The image build (Evo2 + vortex + transformer-engine + flash-attn) is slow the
first time (~15-30 min) and cached afterwards. Scoring a handful of pathway
genes for a handful of species is a few H100-minutes. Keep batches small.
"""

import os
import modal

MODEL_NAME = os.environ.get("EVO2_MODEL", "evo2_7b")
GPU = os.environ.get("EVO2_GPU", "H100")

app = modal.App("psmos-evo2")

# ---------------------------------------------------------------------------
# Image: NVIDIA NGC PyTorch container ships torch + transformer-engine +
# flash-attn PREBUILT, so we skip the slow/fragile source compilation entirely.
# evo2's real dependency tree is light (vtx/vortex + einops + biopython), so on
# top of NGC we only `pip install evo2`. The build is CPU-only (no GPU-seconds);
# GPU is used solely by the scoring class below.
# ---------------------------------------------------------------------------
evo2_image = (
    # 25.04-py3 is the exact base the Arc Institute evo2 Dockerfile validates
    # against; its flash-attn/transformer-engine match vortex's kernels (25.01
    # gave a fwd() signature mismatch during scoring).
    modal.Image.from_registry("nvcr.io/nvidia/pytorch:25.04-py3")
    .apt_install("git")
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
    .pip_install("hf_transfer", "huggingface_hub", "biopython")
    .run_commands("pip install --no-build-isolation evo2 || pip install evo2")
)

# Cache HF weights across runs so we download evo2_7b (~15 GB) only once.
weights_vol = modal.Volume.from_name("psmos-evo2-weights", create_if_missing=True)
WEIGHTS_DIR = "/weights"


@app.cls(
    image=evo2_image,
    gpu=GPU,
    volumes={WEIGHTS_DIR: weights_vol},
    secrets=[
        modal.Secret.from_dict(
            {
                "HF_HOME": WEIGHTS_DIR,
                "HF_TOKEN": os.environ.get("HF_TOKEN", ""),
            }
        )
    ],
    timeout=60 * 60,
    scaledown_window=300,
)
class Evo2Scorer:
    @modal.enter()
    def load(self):
        import torch
        from evo2 import Evo2

        os.environ.setdefault("HF_HOME", WEIGHTS_DIR)
        # NGC's PyTorch >=2.6 defaults torch.load(weights_only=True), which
        # rejects Evo2's checkpoint because it pickles a transformer_engine
        # global (_OverrideLinearPrecision). The Arc Institute checkpoint is
        # trusted, so force weights_only=False for the model load.
        _orig_load = torch.load

        def _trusting_load(*a, **k):
            k["weights_only"] = False
            return _orig_load(*a, **k)

        torch.load = _trusting_load
        try:
            self.torch = torch
            self.model = Evo2(MODEL_NAME)
        finally:
            torch.load = _orig_load
        # warm CUDA
        _ = torch.zeros(1, device="cuda")

    @modal.method()
    def score(self, sequences: list[str]) -> list[float]:
        """Return per-sequence mean log-likelihood (natural log)."""
        # Evo2 exposes score_sequences returning a list of mean log-likelihoods.
        scores = self.model.score_sequences(sequences)
        out = []
        for s in scores:
            try:
                out.append(float(s))
            except Exception:
                out.append(float(getattr(s, "item", lambda: s)()))
        return out

    @modal.method()
    def info(self) -> dict:
        return {"model": MODEL_NAME, "gpu": GPU, "cuda": self.torch.cuda.is_available()}


@app.local_entrypoint()
def main():
    """Smoke test: score two toy sequences and print scores + model info."""
    scorer = Evo2Scorer()
    print("model info:", scorer.info.remote())
    toy = [
        "ATGGCGCCCGCGGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGA",   # ~GC-rich ORF-ish
        "ATATATATATATATATATATATATATATATATATATATATATAT",   # low-complexity
    ]
    scores = scorer.score.remote(toy)
    for seq, sc in zip(toy, scores):
        print(f"len={len(seq):4d}  meanLL={sc:.4f}  {seq[:24]}...")
