# 开源 IBD scRNA 数据集（喂 Geneformer de novo 发现）

> 供 `geneformer_denovo_discovery.py` 使用。据 PubMed + 公开仓库整理（2026-07）。下载大文件走本环境代理即可。

## 旗舰数据集

| 数据集 | 疾病/组织 | 规模 | 获取 | 用途 | 关键论文 (DOI) |
|---|---|---|---|---|---|
| **GSE134809** | 回肠 CD（炎症 vs 未累及） | 22 患者 | GEO：`GSE134809` | A 应答（GIMATS）、C 分型 | Martin 2019 *Cell* [DOI](https://doi.org/10.1016/j.cell.2019.08.008) |
| **SCP259** | 结肠 UC（炎症 vs 健康） | 36.6 万细胞 / 51 亚群 / 30 人 | Broad Single Cell Portal：`SCP259` | C 分型（IAF/BEST4）、B 残留态 | Smillie 2019 *Cell* [DOI](https://doi.org/10.1016/j.cell.2019.06.029) |

- **GSE134809**（Martin 2019）：含 GIMATS 抵抗模块；带 anti-TNF 结局标注，适合**应答分层**（用途 A）与**分子分型**（C）。
- **SCP259**（Smillie 2019）：含 IL13RA2+ IL11+ 炎症成纤维细胞（IAF，anti-TNF 抵抗）、BEST4+ 肠上皮、CD8+IL-17+ T——**炎症分型**的黄金参照。

## 获取方式
- GEO：`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134809`（下 supplementary 的 matrix/barcodes/features 或作者提供的 .h5ad/.loom）。
- Single Cell Portal（SCP259）：`https://singlecell.broadinstitute.org/single_cell/study/SCP259`（需登录下载 expression + metadata）。
- 其他可扩展：Kong 2023 *Immunity* 全肠道 IBD 图谱、Elmentaite 2021 *Nature* 肠道细胞图谱（含儿科 CD）——rsID/accession 使用前请核对。

## 运行（真实数据 + Geneformer 后端）
```bash
# 1) 依赖（GPU 机）
pip install geneformer scanpy anndata torch --extra-index-url https://download.pytorch.org/whl/cu121
# 2) 备好 .h5ad：var 需含 ensembl_id，obs 需含 n_counts / condition / donor
# 3) 跑真实 Geneformer 嵌入 + de novo 发现
python3 geneformer_denovo_discovery.py --h5ad data/GSE134809.h5ad --backend geneformer
```

## 本环境快速自检（无 GPU / 无数据也能跑）
```bash
python3 geneformer_denovo_discovery.py --synthetic --backend pca
# 已验证：恢复炎症相关成纤维细胞态(IL13RA2/CXCL13/TNFRSF11B)，inflamed 富集 1.88x
```

## 预处理要点
- Geneformer 需 **Ensembl 基因 ID** 与每细胞 `n_counts`；用其 `TranscriptomeTokenizer` 分词。
- 严格按**供体/中心**分层切分训练/评估，避免批次效应假信号（IBD 多中心数据最大陷阱）。
- 应答分层需并入临床结局（内镜/组织学缓解）作为标签。

### 署名
数据集论文引用来自 **PubMed**（According to PubMed），DOI 链接见上表。
