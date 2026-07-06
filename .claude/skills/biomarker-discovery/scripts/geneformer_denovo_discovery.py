#!/usr/bin/env python3
"""
geneformer_denovo_discovery.py
IBD 黏膜 scRNA 的 de novo 细胞态发现管线（用途 C 诊断/分型 & A/B 应答态挖掘）。

设计：一条完整、可运行的管线，两个嵌入后端
  --backend pca         轻量、CPU、本环境可跑（numpy/scipy/pandas/scikit-learn）
  --backend geneformer  真实 Geneformer 嵌入（需 torch + 预训练权重 + GPU；代码已写好，见 geneformer_embed()）

流程：load → QC/normalize → embed(backend) → 邻接图聚类(de novo 态)
      → 每个态 one-vs-rest marker → 炎症模块打分 → 与条件(inflamed/healthy 或 responder/NR)关联
      → in-silico 扰动钩子（geneformer 后端）

数据
  --synthetic           生成 IBD 样合成数据(含植入的 IAF 炎症态)，用于自检/演示端到端
  --h5ad PATH           真实数据(如 GSE134809 Martin2019 / SCP259 Smillie2019)，需要 anndata
真实数据获取见 scRNA_datasets.md。

署名：炎症模块基因参考 Smillie 2019 Cell (IL13RA2/IL11/TNFRSF11B IAF, anti-TNF抵抗,
doi:10.1016/j.cell.2019.06.029) 与 Martin 2019 Cell GIMATS (doi:10.1016/j.cell.2019.08.008)。
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))

# 炎症相关模块（IAF + GIMATS + 髓系）——用于给 de novo 态打"炎症"分并锚定
INFLAMMATION_MODULE = [
    "IL13RA2", "IL11", "TNFRSF11B",            # Smillie inflammatory fibroblasts (anti-TNF抵抗)
    "OSM", "OSMR",                             # West/Kokkotis anti-TNF抵抗基质轴
    "TREM1", "S100A8", "S100A9", "IL1B", "CXCL13",  # 炎性髓系/GIMATS
    "TIMP1", "TNF", "CXCL9", "CXCL10",
]

# ---------------- 合成 IBD 样数据(自检用) ----------------
def make_synthetic(n_cells=3000, seed=0):
    rng = np.random.default_rng(seed)
    # 细胞类型及其标记程序
    programs = {
        "Epithelial":        ["EPCAM", "KRT8", "KRT18", "BEST4"],
        "T_cell":            ["CD3D", "CD3E", "CD8A", "IL17A"],
        "Plasma_B":          ["MS4A1", "MZB1", "IGHG1", "CD79A"],
        "Myeloid":           ["LYZ", "CD14", "S100A8", "S100A9", "TREM1", "IL1B"],
        "Fibroblast":        ["COL1A1", "COL3A1", "PDGFRB", "TIMP1"],
        "Inflam_Fibroblast": ["IL13RA2", "IL11", "TNFRSF11B", "OSM", "OSMR", "CXCL13", "CXCL9", "CXCL10", "TNF"],
    }
    genes = sorted({g for gs in programs.values() for g in gs}
                   | {f"NOISE{i}" for i in range(40)})
    gidx = {g: i for i, g in enumerate(genes)}
    # 样本：healthy(低IAF) vs inflamed(高IAF)
    donors = [(f"HC{i}", "healthy") for i in range(6)] + [(f"UC{i}", "inflamed") for i in range(6)]
    base_frac = {"Epithelial": .40, "T_cell": .22, "Plasma_B": .14, "Myeloid": .12, "Fibroblast": .10, "Inflam_Fibroblast": .02}
    X = np.zeros((n_cells, len(genes)), dtype=np.float32)
    ct, cond, don = [], [], []
    for c in range(n_cells):
        d, dc = donors[rng.integers(len(donors))]
        frac = dict(base_frac)
        if dc == "inflamed":                      # 炎症样本里 IAF/髓系扩张
            frac["Inflam_Fibroblast"] = .14; frac["Myeloid"] = .18; frac["Epithelial"] = .30
        types = list(frac); p = np.array([frac[t] for t in types]); p = p / p.sum()
        t = types[rng.choice(len(types), p=p)]
        # 该类型的标记高表达 + 全局噪声
        expr = rng.poisson(0.3, size=len(genes)).astype(np.float32)
        for g in programs[t]:
            expr[gidx[g]] += rng.poisson(8)
        X[c] = expr; ct.append(t); cond.append(dc); don.append(d)
    obs = pd.DataFrame({"cell_type_true": ct, "condition": cond, "donor": don})
    return X, genes, obs

# ---------------- QC + 归一化 ----------------
def normalize(X):
    lib = X.sum(1, keepdims=True); lib[lib == 0] = 1
    Xn = X / lib * 1e4          # CPM-like
    return np.log1p(Xn)

def select_hvg(Xln, n_top=2000):
    v = Xln.var(0)
    return np.argsort(v)[::-1][:min(n_top, Xln.shape[1])]

# ---------------- 嵌入后端 ----------------
def embed_pca(Xln, n_comp=30, seed=0):
    Z = StandardScaler().fit_transform(Xln)
    n_comp = min(n_comp, Z.shape[1] - 1, Z.shape[0] - 1)
    return PCA(n_components=n_comp, random_state=seed).fit_transform(Z)

def geneformer_embed(h5ad_path, model_dir=None, out_dir=None):
    """真实 Geneformer 嵌入。需 GPU + `pip install geneformer` + 预训练权重。
    在本无 GPU 环境不会执行；此处是可直接迁移到 GPU 机的代码路径。"""
    from geneformer import TranscriptomeTokenizer, EmbExtractor  # noqa
    out_dir = out_dir or os.path.join(HERE, "gf_out")
    os.makedirs(out_dir, exist_ok=True)
    # 1) tokenize（需 .loom/.h5ad，含 ensembl_id 与 n_counts）
    tk = TranscriptomeTokenizer({"condition": "condition", "donor": "donor"}, nproc=4)
    tk.tokenize_data(os.path.dirname(h5ad_path), out_dir, "gf", file_format="h5ad")
    # 2) 提取细胞嵌入
    emb = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cell",
                       max_ncells=None, emb_layer=-1, forward_batch_size=16, nproc=4)
    embs = emb.extract_embs(model_dir, os.path.join(out_dir, "gf.dataset"), out_dir, "gf_emb")
    return embs  # DataFrame: 细胞 x 嵌入维

# ---------------- de novo 聚类 ----------------
def cluster_denovo(Z, kmin=4, kmax=10, seed=0):
    best = None
    for k in range(kmin, kmax + 1):
        lab = KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(Z)
        if len(set(lab)) < 2:
            continue
        sil = silhouette_score(Z, lab, sample_size=min(2000, len(lab)), random_state=seed)
        if best is None or sil > best[0]:
            best = (sil, k, lab)
    return best[2], best[1], best[0]

# ---------------- marker: one-vs-rest Welch t ----------------
def markers(Xln, genes, labels, topn=8):
    out = {}
    genes = np.array(genes)
    for cl in sorted(set(labels)):
        m = labels == cl
        t, _ = stats.ttest_ind(Xln[m], Xln[~m], axis=0, equal_var=False)
        t = np.nan_to_num(t)
        top = np.argsort(t)[::-1][:topn]
        out[int(cl)] = [(str(genes[i]), round(float(t[i]), 2)) for i in top]
    return out

# ---------------- 炎症模块打分 + 条件关联 ----------------
def module_score(Xln, genes, module):
    idx = [i for i, g in enumerate(genes) if g in set(module)]
    if not idx:
        return np.zeros(Xln.shape[0])
    z = StandardScaler().fit_transform(Xln[:, idx])
    return z.mean(1)

def condition_enrichment(labels, obs, pos="inflamed"):
    df = obs.copy(); df["cluster"] = labels
    res = {}
    for cl in sorted(set(labels)):
        sub = df[df.cluster == cl]
        frac_pos = (sub.condition == pos).mean()
        overall_pos = (df.condition == pos).mean()
        res[int(cl)] = {"n": int(len(sub)),
                        "frac_" + pos: round(float(frac_pos), 3),
                        "enrichment_vs_overall": round(float(frac_pos / max(overall_pos, 1e-6)), 2)}
    return res

def insilico_perturbation_hook():
    return ("in-silico 扰动(真实版需 geneformer)：对候选态删除/过表达某基因 token，"
            "用模型嵌入度量该态是否被推离'炎症态'——用于区分驱动基因 vs 旁观者。"
            "本环境用 PCA 后端仅提供占位；GPU 上换 geneformer 后端执行。")

# ---------------- 主流程 ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["pca", "geneformer"], default="pca")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--h5ad", default=None)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    print(f"== Geneformer de novo discovery | backend={a.backend} ==")
    if a.synthetic or not a.h5ad:
        print("[data] 合成 IBD 样数据(自检/演示)")
        X, genes, obs = make_synthetic(seed=a.seed)
    else:
        import anndata as ad
        A = ad.read_h5ad(a.h5ad)
        X = A.X.toarray() if hasattr(A.X, "toarray") else np.asarray(A.X)
        genes = list(A.var_names); obs = A.obs.reset_index(drop=True)
        print(f"[data] {a.h5ad}: {X.shape[0]} cells x {X.shape[1]} genes")

    Xln = normalize(X)
    hvg = select_hvg(Xln)
    print(f"[qc] {X.shape[0]} cells, HVG={len(hvg)}")

    if a.backend == "geneformer":
        if not a.h5ad:
            sys.exit("geneformer 后端需 --h5ad 真实数据 + GPU + 权重")
        Z = geneformer_embed(a.h5ad).values
    else:
        Z = embed_pca(Xln[:, hvg], seed=a.seed)
    print(f"[embed] {a.backend} -> {Z.shape[1]} dims")

    labels, k, sil = cluster_denovo(Z, seed=a.seed)
    print(f"[cluster] de novo 态数 k={k} (silhouette={sil:.3f})")

    mk = markers(Xln, genes, labels)
    infl = module_score(Xln, genes, INFLAMMATION_MODULE)
    enr = condition_enrichment(labels, obs)

    # 每个态的平均炎症分
    infl_by_cluster = {int(cl): round(float(infl[labels == cl].mean()), 3) for cl in sorted(set(labels))}
    top_infl = max(infl_by_cluster, key=infl_by_cluster.get)

    print("\n=== de novo 细胞态汇总 ===")
    for cl in sorted(set(labels)):
        star = "  <== 炎症相关态(候选诊断/应答 biomarker)" if cl == top_infl else ""
        print(f"  态{cl}: n={enr[cl]['n']:4d}  炎症分={infl_by_cluster[cl]:+.2f}  "
              f"inflamed富集={enr[cl]['enrichment_vs_overall']}x  "
              f"top marker={[g for g,_ in mk[cl][:5]]}{star}")

    print(f"\n[结论] 炎症相关态 = 态{top_infl}："
          f"炎症模块最高 + inflamed 富集 {enr[top_infl]['enrichment_vs_overall']}x")
    print("[in-silico 扰动] " + insilico_perturbation_hook())

    report = {
        "backend": a.backend, "n_cells": int(X.shape[0]), "k_states": int(k),
        "silhouette": round(float(sil), 3),
        "inflammation_module": INFLAMMATION_MODULE,
        "inflammatory_state_cluster": int(top_infl),
        "state_summary": {int(cl): {"n": enr[cl]["n"],
                                    "inflammation_score": infl_by_cluster[cl],
                                    "inflamed_enrichment": enr[cl]["enrichment_vs_overall"],
                                    "top_markers": [g for g, _ in mk[cl][:8]]}
                          for cl in sorted(set(labels))},
    }
    out = os.path.join(HERE, "denovo_discovery_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    mrows = [{"cluster": cl, "rank": r + 1, "gene": g, "t_stat": t}
             for cl in sorted(set(labels)) for r, (g, t) in enumerate(mk[cl])]
    pd.DataFrame(mrows).to_csv(os.path.join(HERE, "denovo_markers.csv"), index=False)
    print(f"\n✓ 写出 {os.path.relpath(out, HERE)} 与 denovo_markers.csv")


if __name__ == "__main__":
    main()
