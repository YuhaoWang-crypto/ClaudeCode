"""ESM-2 650M real model embedding on Modal T4 GPU"""
import modal

app = modal.App("esm2-embedding")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch", "transformers", "sentencepiece")
)

@app.function(gpu="T4", image=image, timeout=600)
def esm2_embedding():
    import torch
    from transformers import AutoTokenizer, AutoModel
    from torch.nn.functional import cosine_similarity

    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    model_name = "facebook/esm2_t33_650M_UR50D"
    print(f"\nLoading {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).cuda().eval()
    params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {params/1e6:.0f}M")

    seqs = [
        ("spike-RBD",   "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"),
        ("20 std AA",   "ACDEFGHIKLMNPQRSTVWY"),
        ("short peptide","MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNG"),
    ]
    names = [n for n, _ in seqs]
    sequences = [s for _, s in seqs]

    inputs = tokenizer(sequences, return_tensors="pt", padding=True, truncation=True, max_length=512)
    inputs = {k: v.cuda() for k, v in inputs.items()}

    with torch.no_grad():
        out = model(**inputs)

    mask = inputs["attention_mask"].unsqueeze(-1).float()
    embeddings = (out.last_hidden_state * mask).sum(1) / mask.sum(1)

    print(f"\nEmbedding shape per sequence: {embeddings.shape[1]}D")
    for i, name in enumerate(names):
        print(f"  [{i}] {name}: len={len(sequences[i])}")

    print("\nCosine similarities:")
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            sim = cosine_similarity(embeddings[i:i+1], embeddings[j:j+1]).item()
            print(f"  {names[i]} vs {names[j]}: {sim:.4f}")

    return {
        "dim": embeddings.shape[1],
        "gpu": torch.cuda.get_device_name(0),
        "params_M": params // 1_000_000,
    }

@app.local_entrypoint()
def main():
    print("Running ESM-2 650M on Modal cloud GPU...")
    result = esm2_embedding.remote()
    print(f"\n✅ Done! Embedding dim={result['dim']}, GPU={result['gpu']}, Params={result['params_M']}M")
