"""
BioNeMo on Modal Cloud GPU
- Case 1: ESM-2 tokenizer + 批量 embedding（真实模型）
- Case 2: TransformerEngine FP8 矩阵乘法加速演示
- Case 3: ESM-2 650M via NVIDIA NIM API（网络畅通）
"""
import modal

app = modal.App("bionemo-demo")

# 带 PyTorch + transformers 的镜像
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch",
        "transformers",
        "sentencepiece",
        "requests",
        "numpy",
    )
)

# ── Case 1: ESM-2 真实模型 embedding ──────────────────────────────────────
@app.function(gpu="T4", image=image, timeout=300)
def esm2_embedding():
    import torch
    from transformers import AutoTokenizer, AutoModel

    print("=== Case 1: ESM-2 650M Real Model Embedding ===")
    model_name = "facebook/esm2_t33_650M_UR50D"
    print(f"Loading {model_name}...")

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to("cuda").eval()

    seqs = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "ACDEFGHIKLMNPQRSTVWY",
        "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNG",
    ]

    inputs = tokenizer(seqs, return_tensors="pt", padding=True, truncation=True, max_length=512)
    inputs = {k: v.to("cuda") for k, v in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)

    # mean pooling（忽略 padding）
    mask = inputs["attention_mask"].unsqueeze(-1).float()
    embeddings = (outputs.last_hidden_state * mask).sum(1) / mask.sum(1)

    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Embedding shape: {embeddings.shape}")

    # 余弦相似度
    from torch.nn.functional import cosine_similarity
    sim_01 = cosine_similarity(embeddings[0:1], embeddings[1:2]).item()
    sim_02 = cosine_similarity(embeddings[0:1], embeddings[2:3]).item()
    print(f"Similarity seq0 vs seq1: {sim_01:.4f}")
    print(f"Similarity seq0 vs seq2: {sim_02:.4f}")
    return embeddings.cpu().tolist()


# ── Case 2: TransformerEngine FP8 加速 ────────────────────────────────────
te_image = (
    modal.Image.from_registry(
        "nvcr.io/nvidia/pytorch:25.03-py3",
        add_python="3.11",
    )
    .pip_install("transformer-engine[pytorch]")
)

@app.function(gpu="A10G", image=te_image, timeout=300)
def transformer_engine_demo():
    import torch
    import transformer_engine.pytorch as te
    from transformer_engine.common.recipe import Format, DelayedScaling

    print("=== Case 2: TransformerEngine FP8 Acceleration ===")
    print(f"GPU: {torch.cuda.get_device_name(0)}")

    hidden_size = 1024
    batch, seq = 4, 512

    layer = te.TransformerLayer(
        hidden_size=hidden_size,
        ffn_hidden_size=hidden_size * 4,
        num_attention_heads=16,
        layer_type="encoder",
    ).cuda()

    x = torch.randn(seq, batch, hidden_size, device="cuda")
    mask = torch.zeros(batch, 1, seq, seq, device="cuda", dtype=torch.bool)

    fp8_recipe = DelayedScaling(fp8_format=Format.HYBRID, amax_history_len=16, amax_compute_algo="max")

    # BF16 基线
    import time
    with torch.no_grad():
        t0 = time.perf_counter()
        for _ in range(10):
            y_bf16 = layer(x, attention_mask=mask)
        torch.cuda.synchronize()
        bf16_ms = (time.perf_counter() - t0) * 100

    # FP8 加速
    with torch.no_grad(), te.fp8_autocast(enabled=True, fp8_recipe=fp8_recipe):
        t0 = time.perf_counter()
        for _ in range(10):
            y_fp8 = layer(x, attention_mask=mask)
        torch.cuda.synchronize()
        fp8_ms = (time.perf_counter() - t0) * 100

    speedup = bf16_ms / fp8_ms
    print(f"Input shape: {list(x.shape)}  (seq={seq}, batch={batch}, hidden={hidden_size})")
    print(f"BF16 latency: {bf16_ms:.2f} ms")
    print(f"FP8  latency: {fp8_ms:.2f} ms")
    print(f"Speedup: {speedup:.2f}x")
    return {"bf16_ms": bf16_ms, "fp8_ms": fp8_ms, "speedup": speedup}


# ── Case 3: NVIDIA NIM API（Modal 网络畅通）──────────────────────────────
@app.function(image=image, timeout=120, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def nim_esm2_embedding():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    print("=== Case 3: NVIDIA NIM ESM-2 650M (hosted API) ===")

    url = "https://health.api.nvidia.com/v1/biology/meta/esm2-650m/embeddings"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    payload = {"sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK", "pooling": "mean"}

    r = requests.post(url, headers=headers, json=payload, timeout=60)
    r.raise_for_status()
    data = r.json()
    emb = data.get("embedding", [])
    print(f"Embedding dim: {len(emb)}")
    print(f"First 5 values: {emb[:5]}")
    return emb


# ── 入口 ──────────────────────────────────────────────────────────────────
@app.local_entrypoint()
def main():
    print("\n" + "="*60)
    print("BioNeMo on Modal Cloud GPU")
    print("="*60)

    print("\n[1/2] ESM-2 650M real model embedding (T4)...")
    emb = esm2_embedding.remote()
    print(f"Done. Embedding dim={len(emb[0])}")

    print("\n[2/2] TransformerEngine FP8 acceleration (A10G)...")
    te_result = transformer_engine_demo.remote()
    print(f"FP8 speedup: {te_result['speedup']:.2f}x")

    print("\nAll done!")
