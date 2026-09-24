"""B4 —— Tahoe-100M 子集提取。

为什么不下载原始数据
--------------------
Tahoe-100M 的 `metadata/pseudobulk_differential_expression/` 已经是**预计算的
pseudobulk 差异表达**：每行 = (基因, 药物, 浓度, 细胞系) 的 log2FoldChange
（DESeq2，相对同板 DMSO 对照）。这恰好就是 B2 基准需要的 delta，
所以完全不必碰 100M 细胞的原始计数矩阵（3388 个分片）。

数据布局（由本模块实测确认，不靠文档假设）
------------------------------------------
  * 1026 个分片，**每个分片恰好 1 个细胞系** × ~64 个 (药物, 浓度) 条件
    × 62710 个基因；每个条件占约 63 个连续 row group。
  * 因此可以只读需要的 row group + 只读需要的列：单个条件约 0.4 秒，
    而整份 DE 数据是 86.5 GB —— 全量下载既不必要也放不下。
  * 基因顺序**跨条件不一致**，必须按基因名对齐（实测确认）。
  * 每个条件约 55% 的 log2FC 是 NaN（DESeq2 对低表达/离群基因的正常行为）。

产出：data/tahoe/shard_index.json       分片 -> 细胞系 + 条件 -> row group
      results/trackB/tahoe_delta.npz    条件 × 基因 的 log2FC 矩阵
      results/trackB/tahoe_meta.csv     条件元数据
      results/trackB/b4_tahoe.json
"""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .config import DATA, RES_B, log

TAHOE_DIR = DATA / "tahoe"
TAHOE_DIR.mkdir(parents=True, exist_ok=True)

BASE = ("https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/"
        "metadata/pseudobulk_differential_expression")
N_SHARDS = 1026
INDEX_PATH = TAHOE_DIR / "shard_index.json"


def _fs():
    import fsspec
    return fsspec.filesystem("https")


def shard_url(i: int) -> str:
    return f"{BASE}/train-{i:05d}-of-{N_SHARDS:05d}.parquet"


# ------------------------------------------------------------------ 索引
def index_one(i: int) -> dict | None:
    """只读 parquet footer，拿到该分片的细胞系与 (药物,浓度) -> row group 映射。"""
    try:
        with _fs().open(shard_url(i)) as f:
            md = pq.ParquetFile(f).metadata
            names = [md.schema.column(k).name for k in range(md.num_columns)]
            di, ci = names.index("drug"), names.index("concentration")
            cli = names.index("Cell_Name_Vevo")
            pli = names.index("plate")
            # 条件键必须带细胞系：少数边界分片同时含两个细胞系，
            # 只按 (药物,浓度) 建键会把两个细胞系的 row group 混在一起。
            conds: dict[str, list[int]] = {}
            cells, plates = set(), set()
            for r in range(md.num_row_groups):
                rg = md.row_group(r)
                drug = rg.column(di).statistics.min
                conc = round(float(rg.column(ci).statistics.min), 5)
                cell = rg.column(cli).statistics.min
                conds.setdefault(f"{cell}\t{drug}\t{conc}", []).append(r)
                cells.add(cell)
                plates.add(str(rg.column(pli).statistics.min))
            return {"shard": i, "cell_lines": sorted(cells), "plates": sorted(plates),
                    "n_row_groups": md.num_row_groups, "conditions": conds}
    except Exception as exc:  # noqa: BLE001
        log(f"  分片 {i} 索引失败: {exc.__class__.__name__}")
        return None


def build_index(n_workers: int = 16) -> list[dict]:
    if INDEX_PATH.exists():
        log(f"使用缓存索引 {INDEX_PATH.name}")
        return json.loads(INDEX_PATH.read_text())
    log(f"建立分片索引（{N_SHARDS} 个 footer，{n_workers} 线程并发）…")
    out: list[dict] = []
    with ThreadPoolExecutor(n_workers) as ex:
        for k, r in enumerate(ex.map(index_one, range(N_SHARDS)), 1):
            if r:
                out.append(r)
            if k % 100 == 0:
                log(f"  {k}/{N_SHARDS}")
    out.sort(key=lambda d: d["shard"])
    INDEX_PATH.write_text(json.dumps(out))
    log(f"索引完成：{len(out)} 个分片")
    return out


def index_summary(index: list[dict]) -> pd.DataFrame:
    """展开成 (细胞系, 药物, 浓度, 分片) 的长表。"""
    rows = []
    n_multi = 0
    for e in index:
        if len(e["cell_lines"]) != 1:
            n_multi += 1
        for key in e["conditions"]:
            cl, drug, conc = key.split("\t")
            rows.append({"cell_line": cl, "drug": drug, "conc": float(conc),
                         "shard": e["shard"], "key": key})
    if n_multi:
        log(f"  {n_multi} 个分片含多个细胞系（条件键已带细胞系，不受影响）")
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ 提取
def read_condition(pf: pq.ParquetFile, row_groups: list[int]) -> pd.Series | None:
    """读一个条件的 log2FC，按基因名索引。"""
    tb = pf.read_row_groups(row_groups, columns=["gene_name", "log2FoldChange"])
    s = pd.Series(np.asarray(tb.column("log2FoldChange"), dtype=np.float32),
                  index=tb.column("gene_name").to_pylist())
    return s[~s.index.duplicated()]


def read_basemean(pf: pq.ParquetFile, row_groups: list[int]) -> pd.Series:
    """读一个条件的 baseMean，作为基础表达水平的代理。"""
    tb = pf.read_row_groups(row_groups, columns=["gene_name", "baseMean"])
    s = pd.Series(np.asarray(tb.column("baseMean"), dtype=np.float32),
                  index=tb.column("gene_name").to_pylist())
    return s[~s.index.duplicated()]


def gene_list(index: list[dict]) -> list[str]:
    """从第一个分片的第一个条件取基因名全集（各条件顺序不同，故统一到这个列表）。"""
    cache = TAHOE_DIR / "genes.json"
    if cache.exists():
        return json.loads(cache.read_text())
    e = index[0]
    key, rgs = next(iter(e["conditions"].items()))
    with _fs().open(shard_url(e["shard"])) as f:
        s = read_condition(pq.ParquetFile(f), rgs)
    genes = sorted(s.index.astype(str))
    cache.write_text(json.dumps(genes))
    log(f"基因全集 {len(genes)} 个")
    return genes


def gene_universe_from_sample(min_coverage: float = 0.90) -> list[str] | None:
    """用已提取的子集定出基因全集，以控制全量提取时的内存。

    为什么需要：全部 379 个药物 = 56827 个条件；若保留全部 62710 个基因，
    float32 矩阵就是 14.25 GB，超过本机内存。而这些基因里有很大一部分
    在绝大多数条件下都是 NaN（DESeq2 无法检验），保留它们没有信息价值。

    取值口径：在已有的 7500 条件样本（50 个随机药物 × 50 个细胞系 × 3 个浓度）中
    非 NaN 比例 ≥ min_coverage 的基因。这个样本足够大且药物是随机抽取的，
    因此可代表全量。**后续分析中的 ≥98% 过滤在这个全集之内进行**。
    """
    p = RES_B / "tahoe_delta.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    d, g = z["delta"], z["genes"].astype(str)
    frac = (~np.isnan(d)).mean(axis=0)
    keep = sorted(g[frac >= min_coverage].tolist())
    log(f"基因全集由已有 {d.shape[0]} 条件样本定出："
        f"{len(g)} -> {len(keep)}（非 NaN 比例 ≥ {min_coverage:.0%}）")
    log(f"  预计矩阵大小：56827 × {len(keep)} × 4B = "
        f"{56827*len(keep)*4/1e9:.2f} GB")
    return keep


def _open_shard(shard: int, whole: bool):
    """打开分片。

    whole=True  整片下载到内存后本地解析；
    whole=False 通过 HTTP range 按需读取 row group。

    两者的取舍完全取决于**要读这个分片里多大比例的条件**：
      * 只要少数条件时，range 读只传所需数据，胜出；
      * 要几乎全部条件时，逐条件的请求延迟累加起来远超整片传输时间 ——
        实测每片 65 个条件：整片 3.4 秒 vs 逐条件 26 秒（7.7 倍）。
    """
    if not whole:
        return pq.ParquetFile(_fs().open(shard_url(shard))), None
    import io

    import requests
    r = requests.get(shard_url(shard), timeout=600)
    r.raise_for_status()
    buf = io.BytesIO(r.content)
    return pq.ParquetFile(buf), buf


def _extract_shard(args) -> tuple[dict, np.ndarray | None, str | None]:
    """一个分片内的全部所需条件。放在线程里跑（I/O 等待为主）。"""
    shard, entry, keys, genes, want_basemean, whole = args
    out: dict[str, np.ndarray] = {}
    bm = None
    cell = None
    buf = None
    try:
        pf, buf = _open_shard(int(shard), whole)
        # 即使整片已在内存，仍按 row group 逐条件解析：
        # 一次性 to_pandas 整片（约 400 万行含字符串列）会让每线程占用约 1 GB。
        for key in keys:
            rgs = entry["conditions"].get(key)
            if not rgs:
                continue
            out[key] = read_condition(pf, rgs).reindex(genes).to_numpy(dtype=np.float32)
        if want_basemean and keys:
            rgs = entry["conditions"].get(keys[0])
            if rgs:
                bm = read_basemean(pf, rgs).reindex(genes).to_numpy(dtype=np.float32)
                cell = keys[0].split("\t")[0]
    except Exception as exc:  # noqa: BLE001
        log(f"  分片 {shard} 提取失败: {exc.__class__.__name__}: {exc}")
    finally:
        if buf is not None:
            buf.close()
    return out, bm, cell


def extract(index: list[dict], want: pd.DataFrame, genes: list[str],
            n_workers: int = 8) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """按 (细胞系, 药物, 浓度) 提取 log2FC，分片间并行。

    返回 (delta 矩阵, 条件元数据, 各细胞系的 baseMean 代理谱)。

    baseMean 作为"基础表达水平"的代理：它是 DESeq2 在该比较中所有样本的
    归一化计数均值，而对照细胞数远多于处理细胞数（本数据典型为 4862 vs 1378），
    因此以对照为主。用它做 basal gating 与对照相似度加权是近似，
    但不泄漏留出上下文的**扰动方向**信息。
    """
    by_shard = {e["shard"]: e for e in index}
    jobs = []
    for shard, grp in want.groupby("shard"):
        e = by_shard[int(shard)]
        keys = grp["key"].tolist()
        # 该分片里要读的条件占比决定用整片下载还是 range 读（见 _open_shard）
        whole = len(keys) / max(len(e["conditions"]), 1) >= 0.35
        jobs.append((int(shard), e, keys, genes, True, whole))
    n_whole = sum(1 for j in jobs if j[5])
    log(f"提取 {len(want)} 个条件，跨 {len(jobs)} 个分片，{n_workers} 线程并行"
        f"（整片下载 {n_whole} 片 / range 读 {len(jobs)-n_whole} 片）…")

    store: dict[str, np.ndarray] = {}
    basal: dict[str, list[np.ndarray]] = {}
    with ThreadPoolExecutor(n_workers) as ex:
        for k, (got, bm, cell) in enumerate(ex.map(_extract_shard, jobs), 1):
            store.update(got)
            if bm is not None and cell:
                basal.setdefault(cell, []).append(bm)
            if k % 50 == 0:
                log(f"  {k}/{len(jobs)} 个分片，已得 {len(store)} 个条件")

    mat = np.full((len(want), len(genes)), np.nan, dtype=np.float32)
    got_mask = np.zeros(len(want), dtype=bool)
    for i, key in enumerate(want["key"].tolist()):
        v = store.get(key)
        if v is not None:
            mat[i] = v
            got_mask[i] = True

    basal_avg = {k: np.nanmean(np.vstack(v), axis=0) for k, v in basal.items()}
    log(f"提取完成：{int(got_mask.sum())}/{len(want)} 个条件，"
        f"{len(basal_avg)} 个细胞系的基础表达代理")
    return mat, want.assign(extracted=got_mask), basal_avg


# ------------------------------------------------------------------ 主流程
# 与 ChEMBL MERS-CoV 数据结构相同或仅立体异构不同的 Tahoe 药物（由 b5 的比对得出，
# 这里硬编码以保证抽样时一定纳入）
MERS_OVERLAP_DRUGS = [
    "Benztropine (mesylate)", "Carbidopa (monohydrate)", "Gemcitabine",
    "Loperamide (hydrochloride)", "Lopinavir", "Ralimetinib dimesylate",
]


def select_subset(cov: pd.DataFrame, n_drugs: int, seed: int = 42) -> pd.DataFrame:
    """选取药物子集：强制纳入 MERS 交集药物 + 其余**随机**抽取。

    随机抽取而非"挑响应强的"，是为了避免用结果反过来选样本 ——
    响应强弱的筛选交给下游的可重复性门槛（与 sci-Plex 分析保持同一口径）。

    n_drugs <= 0 或 >= 可用药物数时，取**全部**药物（此时不存在抽样问题）。
    """
    per = cov.groupby("drug").agg(n_cell=("cell_line", "nunique"),
                                  n_conc=("conc", "nunique"))
    full = per[(per["n_cell"] == cov["cell_line"].nunique()) & (per["n_conc"] == 3)].index
    forced = [d for d in MERS_OVERLAP_DRUGS if d in set(full)]

    if n_drugs <= 0 or n_drugs >= len(full):
        drugs = sorted(full)
        log(f"药物集：全部 {len(drugs)} 个（无抽样）")
        log(f"  其中 MERS 交集药物 {len(forced)} 个：{forced}")
        sub = cov[cov["drug"].isin(drugs)].copy()
        return sub.drop_duplicates(subset=["cell_line", "drug", "conc"]).reset_index(drop=True)

    rest = sorted(set(full) - set(forced))
    rng = np.random.default_rng(seed)
    n_extra = max(0, n_drugs - len(forced))
    picked = list(rng.choice(rest, size=min(n_extra, len(rest)), replace=False))
    drugs = sorted(set(forced) | set(picked))
    log(f"药物子集：{len(drugs)} 个（强制纳入 MERS 交集 {len(forced)} 个 + 随机 {len(picked)} 个）")
    log(f"  MERS 交集药物：{forced}")
    sub = cov[cov["drug"].isin(drugs)].copy()
    # 同一 (细胞系,药物,浓度) 若出现在多个分片（跨板重复），只保留第一个
    sub = sub.drop_duplicates(subset=["cell_line", "drug", "conc"]).reset_index(drop=True)
    return sub


def main(n_drugs: int = 50, n_workers: int = 8, out_name: str = "tahoe_delta",
         gene_coverage: float | None = None) -> dict:
    index = build_index()
    cov = index_summary(index)
    cov.to_csv(TAHOE_DIR / "coverage.csv", index=False)

    n_cl_all = cov["cell_line"].nunique()
    n_dr_all = cov["drug"].nunique()
    log(f"覆盖：{n_cl_all} 个细胞系 × {n_dr_all} 个药物，共 {len(cov)} 个 (细胞系,药物,浓度) 条件")
    log(f"浓度取值：{sorted(cov['conc'].unique())}")

    want = select_subset(cov, n_drugs)
    log(f"待提取条件：{len(want)}（{want['cell_line'].nunique()} 细胞系 × "
        f"{want['drug'].nunique()} 药物 × {want['conc'].nunique()} 浓度）")

    genes = None
    if gene_coverage is not None:
        genes = gene_universe_from_sample(gene_coverage)
    if genes is None:
        genes = gene_list(index)
    mat, meta, basal = extract(index, want, genes, n_workers)

    cells = sorted(basal)
    np.savez_compressed(
        RES_B / f"{out_name}.npz",
        delta=mat,
        genes=np.array(genes, dtype=object),
        basal=np.vstack([basal[c] for c in cells]) if cells else np.zeros((0, len(genes))),
        basal_cells=np.array(cells, dtype=object),
    )
    meta.to_csv(RES_B / f"{out_name}_meta.csv", index=False)

    frac_nan = float(np.isnan(mat).mean())
    summary = {
        "dataset": "Tahoe-100M (tahoebio/Tahoe-100M, CC0-1.0) pseudobulk_differential_expression",
        "why_not_raw": (
            "DE 表已是相对同板 DMSO 对照的 pseudobulk log2FC，即基准所需的 delta；"
            "全量 DE 为 86.5 GB，按 row group 选择性读取后只取其中极小一部分，"
            "100M 细胞的原始计数矩阵（3388 个分片）完全不需要下载。"),
        "layout_verified": {
            "n_shards": N_SHARDS,
            "n_cell_lines": int(n_cl_all),
            "n_drugs": int(n_dr_all),
            "n_conditions_total": int(len(cov)),
            "concentrations_uM": sorted(float(c) for c in cov["conc"].unique()),
            "n_genes": len(genes),
            "note": "49 个分片跨两个细胞系；条件键含细胞系，已正确处理",
        },
        "subset": {
            "n_drugs": int(want["drug"].nunique()),
            "n_cell_lines": int(want["cell_line"].nunique()),
            "n_conditions_requested": int(len(want)),
            "n_conditions_extracted": int(meta["extracted"].sum()),
            "selection": "强制纳入 MERS 交集药物 + 其余按固定种子随机抽取（不按响应强弱挑选）",
            "mers_overlap_drugs_included": [d for d in MERS_OVERLAP_DRUGS
                                            if d in set(want["drug"])],
        },
        "gene_universe": {
            "n_genes": len(genes),
            "selection": ("由已有子集样本按非 NaN 覆盖率定出（受内存约束：全部 62710 个基因"
                          f"× {len(want)} 个条件的 float32 矩阵为 "
                          f"{len(want)*62710*4/1e9:.1f} GB，超过本机内存）"
                          if gene_coverage is not None else "全部基因"),
            "min_coverage_in_sample": gene_coverage,
        },
        "data_quality": {
            "fraction_nan_log2fc": round(frac_nan, 4),
            "note": "DESeq2 对低表达/离群基因返回 NaN，属正常行为；下游按非 NaN 基因求交集",
        },
        "basal_proxy": {
            "field": "baseMean",
            "n_cell_lines": len(cells),
            "caveat": "baseMean 含处理细胞，但对照细胞数远多于处理细胞（典型 4862 vs 1378），"
                      "故以对照为主；用作 basal gating 与对照相似度的近似，"
                      "不泄漏留出上下文的扰动方向信息",
        },
    }
    (RES_B / f"b4_{out_name.replace('tahoe_delta','tahoe')}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False))
    log(f"NaN 比例 {frac_nan:.1%}")
    log(f"B4 完成 -> {RES_B/'b4_tahoe.json'}")
    return summary


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        # 全量：379 个药物 × 50 个细胞系 × 3 个浓度 = 56827 个条件
        main(n_drugs=0, n_workers=12, out_name="tahoe_delta_full", gene_coverage=0.90)
    else:
        main(int(sys.argv[1]) if len(sys.argv) > 1 else 50)
