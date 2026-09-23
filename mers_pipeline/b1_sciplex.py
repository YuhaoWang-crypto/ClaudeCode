"""B1 —— sci-Plex 化学扰动数据适配。

数据：GSE139944 / GSM4150378 sci-Plex3（A549 / MCF7 / K562 × 188 化合物 × 4 剂量，24 h）

适配到"零样本跨上下文扰动预测"的数据模式：
    细胞系 (cell_type)      -> 上下文 (context)
    化合物 × 剂量            -> 扰动标签 (perturbation)
    溶剂对照 (vehicle=TRUE)  -> 对照 ("non-targeting" 的等价物)

实现要点：UMI 计数矩阵是 2.9 GB 的 `gene cell count` 三元组文本，不做单细胞分析，
只需 pseudobulk，因此**流式累加**，不把矩阵读进内存。

产出：results/trackB/sciplex_pseudobulk.npz（分组 × 基因 的计数矩阵 + 分组元数据）
      results/trackB/b1_sciplex.json
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd

from .config import RES_B, SCIPLEX_DIR, log

PDATA = SCIPLEX_DIR / "GSM4150378_sciPlex3_pData.txt.gz"
GENES = SCIPLEX_DIR / "GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt.gz"
CELLS = SCIPLEX_DIR / "GSM4150378_sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz"
MATRIX = SCIPLEX_DIR / "GSM4150378_UMI.count.matrix.gz"

# 细胞层 QC（沿用 sci-Plex 原文的哈希判定口径）
MIN_UMI = 500
MAX_HASH_QVAL = 0.01
MIN_HASH_RATIO = 5.0
MIN_CELLS_PER_GROUP = 25      # pseudobulk 每组的最少细胞数


def load_pdata() -> pd.DataFrame:
    cols = ["cell", "n.umi", "qval_W", "top_to_second_best_ratio_W", "cell_type",
            "replicate", "dose", "treatment", "vehicle", "product_name", "target",
            "pathway", "catalog_number"]
    log("读取 pData…")
    p = pd.read_csv(PDATA, sep=" ", quotechar='"', usecols=lambda c: c in cols,
                    low_memory=False)
    log(f"  {len(p)} 个细胞，字段 {list(p.columns)}")
    return p


def apply_qc(p: pd.DataFrame) -> pd.DataFrame:
    n0 = len(p)
    p = p[p["n.umi"] >= MIN_UMI]
    n1 = len(p)
    p = p[p["qval_W"].astype(float) < MAX_HASH_QVAL]
    n2 = len(p)
    p = p[p["top_to_second_best_ratio_W"].astype(float) >= MIN_HASH_RATIO]
    n3 = len(p)
    p = p[p["cell_type"].isin(["A549", "MCF7", "K562"])]
    log(f"  QC: {n0} -> UMI≥{MIN_UMI} {n1} -> hash qval<{MAX_HASH_QVAL} {n2} "
        f"-> hash ratio≥{MIN_HASH_RATIO} {n3} -> 三细胞系 {len(p)}")
    return p


def build_groups(p: pd.DataFrame) -> pd.DataFrame:
    """定义 pseudobulk 分组：(细胞系, 扰动)。溶剂对照单独成组。"""
    veh = p["vehicle"].astype(str).str.upper().isin(["TRUE", "T"])
    pert = np.where(
        veh, "vehicle",
        p["product_name"].astype(str) + "__" + p["dose"].astype(str),
    )
    p = p.assign(perturbation=pert, is_vehicle=veh)
    p["group"] = p["cell_type"].astype(str) + "|" + p["perturbation"]
    return p


def stream_pseudobulk(cell_to_group: dict[str, int], n_groups: int,
                      n_genes: int) -> tuple[np.ndarray, np.ndarray]:
    """流式扫描三元组文件，累加每个分组的基因计数。

    文件格式：每行 `gene_idx cell_idx count`，1-based，
    gene_idx 对应 gene.annotations 行号，cell_idx 对应 cell.annotations 行号。
    """
    log("建立 cell_idx -> group 映射…")
    cell_order = pd.read_csv(CELLS, sep="\t", header=None, usecols=[0]).iloc[:, 0]
    grp_of_cell = np.full(len(cell_order) + 1, -1, dtype=np.int32)   # 1-based
    for i, bc in enumerate(cell_order, start=1):
        g = cell_to_group.get(bc)
        if g is not None:
            grp_of_cell[i] = g
    n_keep = int((grp_of_cell >= 0).sum())
    log(f"  矩阵中 {len(cell_order)} 个细胞，其中 {n_keep} 个通过 QC 并归入分组")

    pb = np.zeros((n_groups, n_genes), dtype=np.int64)
    n_lines = 0
    chunk = 20_000_000
    log("流式累加 pseudobulk（不载入完整矩阵）…")
    reader = pd.read_csv(MATRIX, sep="\t", header=None, chunksize=chunk,
                         names=["gene", "cell", "count"],
                         dtype={"gene": np.int32, "cell": np.int32, "count": np.int32})
    for k, ch in enumerate(reader, 1):
        g = grp_of_cell[ch["cell"].to_numpy()]
        keep = g >= 0
        if keep.any():
            np.add.at(pb, (g[keep], ch["gene"].to_numpy()[keep] - 1), ch["count"].to_numpy()[keep])
        n_lines += len(ch)
        log(f"  已处理 {n_lines/1e6:.0f}M 行（第 {k} 块）")
    log(f"  完成，共 {n_lines/1e6:.1f}M 个非零元素")
    n_cells = np.zeros(n_groups, dtype=np.int64)
    return pb, n_cells


def main() -> dict:
    genes = pd.read_csv(GENES, sep=" ", quotechar='"')
    n_genes = len(genes)
    log(f"基因注释 {n_genes} 行")

    p = load_pdata()
    p = apply_qc(p)
    p = build_groups(p)

    counts = p.groupby("group").size()
    keep_groups = counts[counts >= MIN_CELLS_PER_GROUP].index
    p = p[p["group"].isin(keep_groups)]
    log(f"分组: {len(counts)} 个 -> 细胞数≥{MIN_CELLS_PER_GROUP} 的 {len(keep_groups)} 个")

    group_list = sorted(p["group"].unique())
    gidx = {g: i for i, g in enumerate(group_list)}
    cell_to_group = dict(zip(p["cell"], p["group"].map(gidx)))

    meta = (p.groupby("group")
            .agg(cell_type=("cell_type", "first"),
                 perturbation=("perturbation", "first"),
                 is_vehicle=("is_vehicle", "first"),
                 product_name=("product_name", "first"),
                 dose=("dose", "first"),
                 target=("target", "first"),
                 pathway=("pathway", "first"),
                 n_cells=("cell", "size"))
            .reindex(group_list).reset_index())

    pb, _ = stream_pseudobulk(cell_to_group, len(group_list), n_genes)
    meta["total_umi"] = pb.sum(axis=1)

    np.savez_compressed(
        RES_B / "sciplex_pseudobulk.npz",
        counts=pb.astype(np.int64),
        groups=np.array(group_list, dtype=object),
        gene_ids=genes.iloc[:, 0].to_numpy().astype(str),
        gene_names=genes.iloc[:, 1].to_numpy().astype(str),
    )
    meta.to_csv(RES_B / "sciplex_group_meta.csv", index=False)

    summary = {
        "dataset": "GSE139944 / GSM4150378 sci-Plex3 (Srivatsan et al. 2020, Science)",
        "cell_lines": sorted(meta["cell_type"].unique().tolist()),
        "n_groups": int(len(group_list)),
        "n_genes": int(n_genes),
        "n_cells_after_qc": int(len(p)),
        "qc": {"min_umi": MIN_UMI, "max_hash_qval": MAX_HASH_QVAL,
               "min_hash_ratio": MIN_HASH_RATIO, "min_cells_per_group": MIN_CELLS_PER_GROUP},
        "per_cell_line": {
            ct: {
                "n_groups": int((meta["cell_type"] == ct).sum()),
                "n_vehicle_groups": int(((meta["cell_type"] == ct) & meta["is_vehicle"]).sum()),
                "n_cells": int(meta.loc[meta["cell_type"] == ct, "n_cells"].sum()),
                "median_cells_per_group": float(meta.loc[meta["cell_type"] == ct, "n_cells"].median()),
            }
            for ct in sorted(meta["cell_type"].unique())
        },
        "n_unique_compounds": int(meta.loc[~meta["is_vehicle"], "product_name"].nunique()),
        "mapping": {
            "context": "cell_type (A549/MCF7/K562)",
            "perturbation": "product_name + dose",
            "control": "vehicle=TRUE 的溶剂对照组（对应遗传扰动数据里的 non-targeting）",
        },
    }
    (RES_B / "b1_sciplex.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    log(f"B1 完成 -> {RES_B/'b1_sciplex.json'}")
    for ct, v in summary["per_cell_line"].items():
        log(f"  {ct}: {v['n_groups']} 组 / {v['n_cells']} 细胞 / 中位 {v['median_cells_per_group']:.0f} 细胞每组")
    return summary


if __name__ == "__main__":
    main()
