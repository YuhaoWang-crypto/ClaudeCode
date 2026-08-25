"""Embed control cells with Arc's SE-600M, so the context encoder can be tested.

Runtime: **T4 GPU** is enough (the checkpoint is 2.7 GB). This is the one use of
SE-600M this project could not reach on CPU.

The GPU part is only the embedding. The analysis stays in the repo, because the
question is not "what do the embeddings look like" but two specific things the
``transferability-prior-eval`` skill names:

1. **Track resolution.** Correlate every cell line's *predicted* embedding
   against every line's *observed* one and check the diagonal wins. If a line's
   embedding matches another line's profile better than its own, the encoder
   cannot tell these contexts apart and nothing built on it will either. That
   gate has failed before, so it runs before anything else.
2. **Whether a learned context vector beats the raw control pseudobulk** the
   model already uses, on the cross-context axis.

Expect little. The cross-context axis is close to its measurement floor -- two
runs of the same K562 screen agree at r = 0.319 against this model's 0.328 --
and the skill measured an *oracle* context feature scoring worse than no context
feature at all. This is worth testing and worth not expecting much from.

Output: one ``se_embeddings.npz``, a few MB, which is small enough to attach or
push back. Nothing large moves.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORK = Path("/content/se")
OFFICIAL = Path("/content/vcc_official")
MODEL = WORK / "SE-600M"
CONTEXTS = ("A", "B", "C")
N_CELLS = 3000          # per context; enough for a stable mean, cheap to embed


def sh(cmd: str) -> int:
    print(f"$ {re.sub(r'vcc_pat_[A-Za-z0-9]+', '***', cmd)}", flush=True)
    return subprocess.call(cmd, shell=True)


def token() -> str:
    try:
        from google.colab import userdata            # type: ignore
        t = userdata.get("VCC_TOKEN")
        if t:
            return t.strip()
    except Exception:
        pass
    t = os.environ.get("VCC_TOKEN", "").strip()
    if not t:
        sys.exit("Set VCC_TOKEN in Colab's Secrets panel (key icon, left bar).")
    return t


def main() -> None:
    import numpy as np

    WORK.mkdir(parents=True, exist_ok=True)
    OFFICIAL.mkdir(parents=True, exist_ok=True)

    try:
        import torch
        print(f"GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NONE'}")
    except Exception:
        print("torch not importable yet")

    sh("pip -q install arc-state anndata h5py scipy vcc-cli huggingface_hub")

    # 1. the checkpoint
    if not (MODEL / "config.yaml").exists():
        sh(f"huggingface-cli download arcinstitute/SE-600M --local-dir {MODEL}")

    # 2. the official controls, via your own key
    if not (OFFICIAL / "gene_names.csv").exists():
        p = subprocess.run("vcc login --token-stdin", shell=True,
                           input=token(), text=True, capture_output=True)
        if p.returncode != 0:
            sys.exit(re.sub(r"vcc_pat_[A-Za-z0-9]+", "***", p.stderr)[-800:])
        sh(f"cd {OFFICIAL} && vcc datasets download controls -o controls.zip "
           f"&& unzip -oq controls.zip")

    # 3. subsample each context's controls so the embedding job is small
    import anndata as ad
    import h5py
    from scipy import sparse

    genes = np.loadtxt(OFFICIAL / "gene_names.csv", dtype=str, skiprows=1)
    for c in CONTEXTS:
        sub = WORK / f"ctrl_{c}.h5ad"
        if sub.exists():
            continue
        with h5py.File(OFFICIAL / f"context_{c}.h5ad") as f:
            X = f["X"]
            n, g = (int(v) for v in X.attrs["shape"])
            rng = np.random.default_rng(0)
            rows = np.sort(rng.choice(n, min(N_CELLS, n), replace=False))
            ip = X["indptr"][:]
            data, indices, indptr = [], [], [0]
            for r in rows:
                a, b = int(ip[r]), int(ip[r + 1])
                data.append(X["data"][a:b]); indices.append(X["indices"][a:b])
                indptr.append(indptr[-1] + (b - a))
            M = sparse.csr_matrix(
                (np.concatenate(data).astype(np.float32),
                 np.concatenate(indices), np.array(indptr)),
                shape=(rows.size, g))
        a = ad.AnnData(X=M)
        a.var_names = genes
        a.var["gene_name"] = genes
        a.obs["context"] = c
        a.write_h5ad(sub)
        print(f"  wrote {sub.name}: {rows.size} cells x {g} genes", flush=True)

    # 4. embed
    out = {}
    for c in CONTEXTS:
        emb = WORK / f"emb_{c}.h5ad"
        if not emb.exists():
            rc = sh(f"state emb --model-folder {MODEL} "
                    f"--input {WORK}/ctrl_{c}.h5ad --output {emb}")
            if rc != 0:
                sys.exit(f"state emb failed on context {c}")
        a = ad.read_h5ad(emb)
        key = next((k for k in a.obsm if "state" in k.lower() or "emb" in k.lower()),
                   None)
        if key is None:
            sys.exit(f"no embedding in {emb}; obsm keys: {list(a.obsm)}")
        out[c] = np.asarray(a.obsm[key], dtype=np.float32)
        print(f"  context {c}: {out[c].shape} from obsm[{key!r}]")

    dest = WORK / "se_embeddings.npz"
    np.savez_compressed(dest, **{f"ctx_{c}": v for c, v in out.items()})
    size = dest.stat().st_size / 1e6

    print("\n=== PASTE BACK ===")
    for c, v in out.items():
        m = v.mean(0)
        print(f"context {c}: {v.shape[0]} cells x {v.shape[1]}d, "
              f"|mean|={np.linalg.norm(m):.3f}")
    import itertools
    for a_, b_ in itertools.combinations(CONTEXTS, 2):
        ma, mb = out[a_].mean(0), out[b_].mean(0)
        r = float(ma @ mb / (np.linalg.norm(ma) * np.linalg.norm(mb) + 1e-9))
        print(f"cosine({a_},{b_}) = {r:+.4f}")
    print(f"file: se_embeddings.npz, {size:.1f} MB")
    print("=== END ===")
    print(f"\nDownload {dest} and attach it, or:")
    print(f"  from google.colab import files; files.download('{dest}')")


if __name__ == "__main__":
    main()
