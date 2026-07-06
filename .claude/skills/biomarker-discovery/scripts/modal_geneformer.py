"""
modal_geneformer.py — 轻量 GPU Geneformer 真实嵌入 + de novo 细胞态发现（在 Modal 上跑）。

做的是真实 Geneformer：加载官方预训练权重 + 官方 token/gene-median 字典，
按 Geneformer 的 rank-value 编码把细胞转成 token 序列，在 GPU 上前向取细胞嵌入，
再做 de novo 聚类 + 炎症富集。数据用合成但**真实 Ensembl ID**(从 Geneformer 词表采样)，
保证 tokenizer 兼容、无需下载大 scRNA 数据 —— 真正"轻量"。

运行：
  MODAL_TOKEN_ID/SECRET 已在环境中。
  python3 -m modal run modal_geneformer.py

已验证：真实 V1-10M 预训练权重 + gc30M 官方词表 + rank-value 编码 + GPU 前向端到端跑通
(见 geneformer-gpu-run-results.md)。注意：合成随机基因数据无真实共表达结构，Geneformer
预训练嵌入不会将其分离(silhouette≈随机)——这是 Geneformer 编码"真实生物先验"的本质。
真正的 de novo 发现需真实 scRNA(GSE134809/SCP259)：把 build_synthetic 段换成读 .h5ad 的
rank-value 编码(counts∩词表)，其余流程不变。
"""
import modal

app = modal.App("gf-ibd-denovo")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch", "transformers", "huggingface_hub", "safetensors",
        "numpy", "scipy", "scikit-learn", "pandas",
    )
)
cache = modal.Volume.from_name("gf-hf-cache", create_if_missing=True)
REPO = "ctheodoris/Geneformer"


@app.function(image=image, gpu="T4", volumes={"/cache": cache}, timeout=1800)
def run(n_cells: int = 1200, top_k: int = 1024, seed: int = 0):
    import os, glob, pickle, json
    import numpy as np
    import torch
    from huggingface_hub import HfApi, hf_hub_download
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    os.environ["HF_HOME"] = "/cache/hf"
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[env] torch {torch.__version__} | device={dev} | "
          f"gpu={torch.cuda.get_device_name(0) if dev=='cuda' else 'none'}")

    # ---- 1) 下载通用预训练模型 + 匹配词表的官方字典 ----
    # 钉死 V1-10M（通用预训练、最轻量）+ V1 的 gc30M 字典（词表必须与模型匹配！）
    files = set(HfApi().list_repo_files(REPO))
    model_dir = "Geneformer-V1-10M"
    tok_f = "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    med_f = "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl"
    assert model_dir + "/config.json" in files and tok_f in files and med_f in files, "repo 结构变化，请核对"
    print(f"[repo] model={model_dir} (general pretrained) token={tok_f}")

    dl = lambda f: hf_hub_download(REPO, f, cache_dir="/cache/hf")
    tok_path, med_path = dl(tok_f), dl(med_f)
    for f in [model_dir + "/config.json", model_dir + "/model.safetensors"]:
        dl(f)
    mdir = os.path.dirname(dl(model_dir + "/config.json"))
    print(f"[model] dir={mdir}")

    with open(tok_path, "rb") as fh: token_dict = pickle.load(fh)
    with open(med_path, "rb") as fh: median_dict = pickle.load(fh)
    ens = [g for g in token_dict if isinstance(g, str) and g.startswith("ENSG") and g in median_dict]
    print(f"[dict] token_dict={len(token_dict)} usable ENSG={len(ens)}")

    # ---- 2) 合成 IBD 样细胞(真实 Ensembl ID)，植入炎症程序 ----
    rng = np.random.default_rng(seed)
    n_genes = 2000                         # 远大于 top_k，使 rank 截断真正生效
    panel = list(rng.choice(ens, size=n_genes, replace=False))
    infl_idx = list(range(60))             # 前 60 个为"炎症程序"
    med = np.array([median_dict[g] for g in panel], dtype=np.float64)
    med[med <= 0] = 1.0
    tokens = np.array([token_dict[g] for g in panel])

    counts = rng.poisson(0.2, size=(n_cells, n_genes)).astype(np.float64)
    cond = np.array(["inflamed" if i % 2 else "healthy" for i in range(n_cells)])
    for c in range(n_cells):
        if cond[c] == "inflamed":          # 炎症细胞：强烈激活炎症程序 → 挤进 top-k
            for i in rng.choice(infl_idx, size=35, replace=False):
                counts[c, i] += rng.poisson(25)
        else:                              # 健康细胞：激活其它随机程序
            for i in rng.choice(range(60, n_genes), size=35, replace=False):
                counts[c, i] += rng.poisson(25)

    # ---- 3) Geneformer rank-value 编码 ----
    cls = token_dict.get("<cls>")
    def encode(row):
        norm = row / max(row.sum(), 1.0) * 1e4
        scaled = norm / med
        nz = np.where(scaled > 0)[0]
        order = nz[np.argsort(scaled[nz])[::-1]][:top_k]
        ids = tokens[order].tolist()
        return ([cls] + ids) if cls is not None else ids

    seqs = [encode(counts[c]) for c in range(n_cells)]
    maxlen = min(max(len(s) for s in seqs), top_k + 1)
    pad = token_dict.get("<pad>", 0)

    # ---- 4) 加载真实模型，GPU 前向取细胞嵌入 ----
    from transformers import AutoModel
    model = AutoModel.from_pretrained(mdir).to(dev).eval()
    print(f"[model] loaded {model.__class__.__name__}")

    embs = []
    bs = 32
    with torch.no_grad():
        for i in range(0, n_cells, bs):
            batch = seqs[i:i + bs]
            L = min(max(len(s) for s in batch), maxlen)
            ids = torch.full((len(batch), L), pad, dtype=torch.long)
            mask = torch.zeros((len(batch), L), dtype=torch.long)
            for j, s in enumerate(batch):
                s = s[:L]; ids[j, :len(s)] = torch.tensor(s); mask[j, :len(s)] = 1
            out = model(input_ids=ids.to(dev), attention_mask=mask.to(dev))
            h = out.last_hidden_state
            m = mask.unsqueeze(-1).to(dev).float()
            cell = (h * m).sum(1) / m.sum(1).clamp(min=1)   # mean-pool
            embs.append(cell.cpu().numpy())
    Z = np.concatenate(embs, 0)
    print(f"[embed] Geneformer cell embeddings: {Z.shape}")

    # ---- 5) de novo 聚类 + 炎症富集 ----
    best = None
    for k in range(4, 9):
        lab = KMeans(k, random_state=seed, n_init=10).fit_predict(Z)
        sil = silhouette_score(Z, lab, sample_size=min(1000, len(lab)), random_state=seed)
        if best is None or sil > best[0]: best = (sil, k, lab)
    sil, k, labels = best

    # 每个细胞的炎症 token 占比（用真实编码里炎症基因 token 的存在度）
    infl_tokens = set(tokens[infl_idx].tolist())
    infl_frac = np.array([len(infl_tokens & set(s)) / max(len(s), 1) for s in seqs])
    summary = {}
    for cl in sorted(set(labels)):
        msk = labels == cl
        summary[int(cl)] = {
            "n": int(msk.sum()),
            "mean_infl_token_frac": round(float(infl_frac[msk].mean()), 4),
            "inflamed_frac": round(float((cond[msk] == "inflamed").mean()), 3),
        }
    top = max(summary, key=lambda c: summary[c]["mean_infl_token_frac"])
    report = {
        "backend": "geneformer_gpu", "device": dev,
        "model_class": model.__class__.__name__, "model_dir": os.path.basename(mdir),
        "n_cells": n_cells, "emb_dim": int(Z.shape[1]),
        "k_states": int(k), "silhouette": round(float(sil), 3),
        "inflammatory_state_cluster": int(top),
        "state_summary": summary,
    }
    print("[result] " + json.dumps(report, indent=2))
    return report


@app.local_entrypoint()
def main(n_cells: int = 1200, top_k: int = 1024):
    import json
    rep = run.remote(n_cells=n_cells, top_k=top_k)
    print("\n==== RETURNED REPORT ====")
    print(json.dumps(rep, indent=2, ensure_ascii=False))
    with open("geneformer_gpu_report.json", "w") as f:
        json.dump(rep, f, indent=2, ensure_ascii=False)
    print("saved geneformer_gpu_report.json")
