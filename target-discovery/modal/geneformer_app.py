"""
Modal app: get Geneformer genuinely running on a GPU.

Stage 1 (this file, `smoke`): build an image with the geneformer package,
clone the pretrained checkpoint + gene/token dictionaries into a cached Volume,
import geneformer, and confirm a checkpoint loads on the GPU.

Run:  modal run modal_geneformer.py::smoke
"""
import modal

app = modal.App("gf-ipf")

# geneformer + its scientific deps. Pin nothing exotic; let pip resolve.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "git-lfs")
    # transformers pinned: Geneformer does `from transformers import
    # SpecialTokensMixin`, which newer transformers dropped from the top-level
    # namespace. 4.46.3 still exports it and matches Geneformer-V2's era.
    .pip_install(
        "torch", "transformers==4.46.3", "datasets", "numpy", "pandas", "scipy",
        "scikit-learn", "loompy", "anndata", "scanpy", "pyarrow",
        "tdigest", "peft", "optuna", "ray",
    )
    # Full clone (pip's default --filter=blob:none partial clone fails against
    # HF's git promisor remote). LFS smudge skipped so the clone stays small;
    # then pull only the gene/token dictionaries so the installed package has
    # real .pkl files (not LFS pointers), and pip install from the local dir.
    .run_commands(
        "GIT_LFS_SKIP_SMUDGE=1 git clone "
        "https://huggingface.co/ctheodoris/Geneformer.git /opt/Geneformer",
        "cd /opt/Geneformer && git lfs pull "
        '--include="geneformer/*.pkl" '
        '--include="geneformer/gene_dictionaries_30m/*"',
        "pip install /opt/Geneformer",
    )
)

# Persisted cache for the (large, LFS) model weights so we clone once.
model_vol = modal.Volume.from_name("gf-models", create_if_missing=True)
MODEL_DIR = "/models"
# Current default V2 checkpoint; its gc104M dictionaries ship in the package.
CKPT = "Geneformer-V2-104M"


@app.function(image=image, volumes={MODEL_DIR: model_vol}, timeout=1800)
def fetch_checkpoint():
    """Clone only the chosen checkpoint dir + dictionaries into the Volume."""
    import os
    from huggingface_hub import snapshot_download

    dest = f"{MODEL_DIR}/Geneformer"
    weights = f"{dest}/{CKPT}/model.safetensors"

    def real_weights() -> bool:
        # A real safetensors file is 100s of MB; an LFS pointer is ~130 bytes.
        return os.path.exists(weights) and os.path.getsize(weights) > 10_000_000

    if real_weights():
        print(f"[cache] {CKPT} weights already present "
              f"({os.path.getsize(weights)/1e6:.0f} MB)")
        return sorted(os.listdir(f"{dest}/{CKPT}"))

    os.makedirs(dest, exist_ok=True)
    # HF CDN download (reliable) of only the chosen checkpoint + gc dictionaries.
    snapshot_download(
        repo_id="ctheodoris/Geneformer",
        allow_patterns=[f"{CKPT}/*", "geneformer/*.pkl"],
        local_dir=dest,
    )
    model_vol.commit()
    assert real_weights(), f"download did not produce real weights at {weights}"
    print(f"[done] {CKPT} weights: {os.path.getsize(weights)/1e6:.0f} MB")
    return sorted(os.listdir(f"{dest}/{CKPT}"))


@app.function(image=image, timeout=120)
def show_reqs():
    """Print Geneformer's declared deps + installed transformers version."""
    import subprocess
    for f in ("/opt/Geneformer/setup.py", "/opt/Geneformer/requirements.txt"):
        print(f"===== {f} =====")
        try:
            print(open(f).read())
        except FileNotFoundError:
            print("(missing)")
    print("===== installed =====")
    print(subprocess.run(["pip", "show", "transformers"], capture_output=True, text=True).stdout)


@app.function(image=image, volumes={MODEL_DIR: model_vol}, timeout=300)
def list_repo():
    """List top-level dirs in the cloned repo to find real checkpoint names."""
    import os
    dest = f"{MODEL_DIR}/Geneformer"
    if not os.path.exists(dest):
        return {"error": "repo not cloned yet"}
    top = sorted(os.listdir(dest))
    # candidate model dirs = those containing a config.json
    models = []
    for d in top:
        p = f"{dest}/{d}"
        if os.path.isdir(p) and "config.json" in os.listdir(p):
            models.append(d)
    print("TOP_LEVEL:", top)
    print("MODEL_DIRS:", models)
    # also show geneformer package dictionaries present
    gf = f"{dest}/geneformer"
    if os.path.exists(gf):
        print("GENEFORMER_PKGDIR:", sorted(x for x in os.listdir(gf) if x.endswith(('.pkl', '.txt')) or 'dict' in x))
    return {"top_level": top, "model_dirs": models}


@app.function(image=image, volumes={MODEL_DIR: model_vol}, gpu="T4", timeout=900)
def smoke():
    """Load the checkpoint on the GPU and run one forward pass."""
    import os
    import torch
    from transformers import AutoModel

    ckpt_path = f"{MODEL_DIR}/Geneformer/{CKPT}"
    files = sorted(os.listdir(ckpt_path))
    model = AutoModel.from_pretrained(ckpt_path).to("cuda").eval()
    n_params = sum(p.numel() for p in model.parameters())
    # Minimal forward pass with dummy token ids in-range.
    vocab = model.get_input_embeddings().num_embeddings
    ids = torch.randint(0, vocab, (1, 16), device="cuda")
    with torch.no_grad():
        out = model(ids)
    import geneformer
    return {
        "geneformer_version": getattr(geneformer, "__version__", "unknown"),
        "ckpt_files": files,
        "n_params_millions": round(n_params / 1e6, 1),
        "vocab_size": vocab,
        "hidden_state_shape": tuple(out.last_hidden_state.shape),
        "cuda_device": torch.cuda.get_device_name(0),
    }


@app.local_entrypoint()
def main():
    print("Fetching checkpoint into Volume...")
    print("CKPT FILES:", fetch_checkpoint.remote())
    print("Running GPU smoke...")
    res = smoke.remote()
    for k, v in res.items():
        print(f"  {k}: {v}")
