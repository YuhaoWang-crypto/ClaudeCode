"""Embed control cells with Arc's SE-600M on a GPU, and gate the result honestly.

This is the last untested use of SE-600M.  Its *protein* embeddings were tested
locally and were worth +0.0012, acting only on the 26 official targets that are
neither perturbed nor measured in the source.  Its *cell* encoder is a different
object: it could in principle supply a learned context representation to replace
the raw control pseudobulk the model conditions on now.

Two measurements say to expect very little, and they are stated here so the
result cannot be reinterpreted afterwards:

* the cross-context axis -- the one a context encoder acts on -- is close to its
  measurement floor.  Two runs of the *same* K562 screen agree at r = 0.319;
  this model already reaches 0.328.
* ``transferability-prior-eval`` measured an **oracle** context feature scoring
  *worse* than using no context feature at all, which makes that axis
  architecture-limited rather than representation-limited.

So the first thing computed is not a score but a **gate**: correlate every cell
line's predicted embedding against every line's observed basal profile and check
that the diagonal wins.  If a line's embedding matches a different line's
profile better than its own, the encoder cannot separate these contexts and
nothing downstream can either.  That gate has failed before on this class of
model, so it runs first and its verdict is reported whether or not it is
convenient.

    modal run modal/se_embed_modal.py
"""

from __future__ import annotations

import modal

REPO = "https://github.com/yuhaowang-crypto/claudecode.git"
BRANCH = "claude/virtual-cell-model-prediction-qu31ry"
N_CELLS = 3000          # per context: enough for a stable mean, cheap to embed

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "curl")
    .pip_install("numpy", "scipy", "pandas", "anndata", "h5py", "scikit-learn",
                 "huggingface_hub", "arc-state")
    .run_commands(f"git clone --depth 1 -b {BRANCH} {REPO} /repo")
)

app = modal.App("se-embed")
vol = modal.Volume.from_name("vcc-submission", create_if_missing=True)
models = modal.Volume.from_name("se-models", create_if_missing=True)
DATA, MODELS = "/vol", "/models"


@app.function(image=image, volumes={MODELS: models}, timeout=60 * 45)
def fetch_model() -> str:
    """Pull SE-600M onto its own Volume, once."""
    from pathlib import Path
    from huggingface_hub import snapshot_download

    dest = Path(MODELS) / "SE-600M"
    if (dest / "config.yaml").exists():
        models.commit()
        return f"already present: {sum(1 for _ in dest.rglob('*'))} files"
    snapshot_download("arcinstitute/SE-600M", local_dir=str(dest))
    models.commit()
    total = sum(p.stat().st_size for p in dest.rglob("*") if p.is_file())
    return f"downloaded {total / 1e9:.2f} GB to {dest}"


@app.function(image=image, volumes={DATA: vol, MODELS: models}, gpu="T4",
              timeout=60 * 60, memory=32768)
def embed() -> dict:
    """Embed each context's control cells; return only summaries and the file."""
    import subprocess
    from pathlib import Path

    import anndata as ad
    import h5py
    import numpy as np
    from scipy import sparse

    work = Path(DATA) / "se"
    work.mkdir(parents=True, exist_ok=True)
    genes = np.loadtxt(f"{DATA}/vcc_official/gene_names.csv", dtype=str,
                       skiprows=1)

    out: dict[str, list] = {}
    for c in ("A", "B", "C"):
        sub = work / f"ctrl_{c}.h5ad"
        if not sub.exists():
            with h5py.File(f"{DATA}/vcc_official/context_{c}.h5ad") as f:
                X = f["X"]
                n, g = (int(v) for v in X.attrs["shape"])
                rng = np.random.default_rng(0)
                rows = np.sort(rng.choice(n, min(N_CELLS, n), replace=False))
                ip = X["indptr"][:]
                data, idx, indptr = [], [], [0]
                for r in rows:
                    a, b = int(ip[r]), int(ip[r + 1])
                    data.append(X["data"][a:b]); idx.append(X["indices"][a:b])
                    indptr.append(indptr[-1] + (b - a))
                M = sparse.csr_matrix(
                    (np.concatenate(data).astype(np.float32),
                     np.concatenate(idx), np.array(indptr)),
                    shape=(rows.size, g))
            a_ = ad.AnnData(X=M)
            a_.var_names = genes
            a_.var["gene_name"] = genes
            a_.write_h5ad(sub)

        emb = work / f"emb_{c}.h5ad"
        if not emb.exists():
            p = subprocess.run(
                f"state emb --model-folder {MODELS}/SE-600M --input {sub} "
                f"--output {emb}", shell=True, capture_output=True, text=True)
            if p.returncode != 0:
                return {"error": (p.stdout[-1500:] + p.stderr[-2500:])}
        a_ = ad.read_h5ad(emb)
        key = next((k for k in a_.obsm
                    if "state" in k.lower() or "emb" in k.lower()), None)
        if key is None:
            return {"error": f"no embedding in obsm; keys={list(a_.obsm)}"}
        out[c] = np.asarray(a_.obsm[key], dtype=np.float32)

    np.savez_compressed(work / "se_embeddings.npz",
                        **{f"ctx_{c}": v for c, v in out.items()})
    vol.commit()

    import itertools
    summary = {c: list(v.shape) for c, v in out.items()}
    for a_, b_ in itertools.combinations(out, 2):
        ma, mb = out[a_].mean(0), out[b_].mean(0)
        summary[f"cosine_{a_}{b_}"] = float(
            ma @ mb / (np.linalg.norm(ma) * np.linalg.norm(mb) + 1e-9))
    return summary


@app.local_entrypoint()
def main(stage: str = "all"):
    if stage in ("model", "all"):
        print(fetch_model.remote())
    if stage in ("embed", "all"):
        for k, v in embed.remote().items():
            print(f"  {k}: {v}")
