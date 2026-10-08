# ArchMap 参考图谱数据库

把 [ArchMap → References → Atlases](https://www.archmap.bio/#/references/atlases) 上全部
**17 个单细胞参考图谱**的数据下载、核对，并整合成一个 SQLite 数据库 `archmap.sqlite`（36 MB；仓库里存的是压缩版 `archmap.sqlite.gz`，
`gunzip -k archmap.sqlite.gz` 解压即可），
同时导出 Excel 可直接打开的 CSV（`exports/`）。

## 数据从哪来

ArchMap 前端是 React 单页应用，背后是一个 REST API（`archmap_api.py`）：

| 端点 | 内容 |
|---|---|
| `GET /temp_auth` | 匿名访客临时 JWT（网站本身就这样拿） |
| `GET /atlases` | 17 个图谱的元数据（物种、细胞数、DOI、batch/cell-type key…） |
| `GET /models` | scVI / scANVI / scPoli 三个映射模型 |
| `GET /scvi-atlases` | 25 个 scvi-hub (HuggingFace) 图谱、50 个预训练模型 |
| `POST /file_download/atlas_files` | 每个图谱的下载文件（7 天有效的签名链接） |

每个图谱提供 `data.h5ad`（模型特征上的参考数据）、`data_only_count.h5ad`（原始计数）和训练好的模型。
**全部文件约 103 GB**（PBMC 原始计数一个就 50 GB），所以数据库构建**不下载表达矩阵**：
`h5ad_remote.py` 用 HTTP Range 请求只读取每个 `data.h5ad` 的 `obs`（细胞注释）和 `var`（基因），
17 个图谱共约 3 分钟。需要完整矩阵时用 `download_atlas.py`。

## 数据库内容

| 表 / 视图 | 行数 | 内容 |
|---|---|---|
| `atlas` | 17 | 门户元数据 + 实际文件的细胞数/基因数/编码、实际使用的 cell-type/batch 列 |
| `model`, `atlas_model` | 3 / 18 | 映射模型及图谱兼容性 |
| `atlas_file` | 67 | 所有可下载文件及大小 |
| `obs_column` | 648 | 每个图谱的每个细胞注释列（类型、类别数、缺失、数值统计） |
| `obs_category_count` | 26,838 | 每个分类列每个取值的细胞数（细胞类型、供体、组织、疾病、测序平台…） |
| `cell_type` | 1,330 | 每个图谱的细胞类型及细胞数、占比 |
| `celltype_batch_count` | 49,458 | 细胞类型 × 批次（样本/供体）细胞数 |
| `celltype_ontology_map` | 781 | 图谱标签 → Cell Ontology (CL) ID（CELLxGENE 格式的图谱） |
| `gene` | 102,787 | 每个图谱的模型特征基因（symbol 与 Ensembl 跨图谱互相补全） |
| `scvi_hub_model` | 50 | scvi-hub 预训练模型及 HuggingFace 链接 |
| `qc_issue` | 14 | 门户元数据与实际文件不一致之处 |
| `v_atlas_summary` | | 一图谱一行的总览 |
| `v_gene_presence` | | 每个基因出现在几个图谱中（按物种） |
| `v_ontology_across_atlases` | | 跨图谱共有的 Cell Ontology 细胞类型 |
| `v_celltype_ontology` | | 每个图谱标签的主要 CL 对应 |

## 发现的数据问题（`qc_issue`）

- **Plaque (refined cell type level)** 的 `data.h5ad` 与 **Breast** 图谱的文件完全相同
  （大小、形状、基因都一致，212 万个细胞），它自己的数据并不在这个文件里。跨图谱视图已排除它。
- **Glioblastoma** 的文件只有门户标称细胞数的 29.8%（338,564 vs 1,135,677）；
  HEOCA 86%、HNOCA / HNOCA Extended 约 98%；Plaque 反而多 40%。
- 门户填写的 cell-type key 在文件里不存在：HEOCA（大小写不同 `Cell_type`，已自动匹配）、
  Retina（改用 `celltype`）、HNOCA（`snapseed_pca_rss_level_123` 不存在，改用 `annot_level_2`）。
  替换都记录在 `atlas.cell_type_key_used` 和 `qc_issue` 中。
- 注意：HLCA 与 HLCA retrained 是同一批细胞配不同模型，跨图谱求和时会重复计数。

## 用法

```bash
pip install h5py fsspec aiohttp numpy
python build_db.py                    # 重新构建（约 3 分钟）
python build_db.py --no-remote        # 只要元数据和文件清单
python download_atlas.py --list
python download_atlas.py HLCA Heart --kinds reference_data model --out /data/archmap
```

查询示例：

```sql
-- 各图谱概览
SELECT * FROM v_atlas_summary;

-- 哪些图谱有 T 细胞，各多少
SELECT a.name, c.cell_type, c.n_cells FROM cell_type c JOIN atlas a USING (atlas_id)
WHERE c.cell_type LIKE '%T cell%' OR c.cell_type LIKE 'CD4%' OR c.cell_type LIKE 'CD8%'
ORDER BY c.n_cells DESC;

-- 某基因被哪些图谱的模型使用
SELECT * FROM v_gene_presence WHERE symbol = 'CD3E';

-- HLCA 每个细胞类型来自多少个数据集
SELECT cell_type, COUNT(*) AS n_batches, SUM(n_cells) FROM celltype_batch_count
WHERE atlas_id = '628668716f930d8b7f44d575' GROUP BY cell_type ORDER BY 3 DESC;
```

`raw/` 是 API 原始返回的快照（已去除签名下载链接）。数据版权归各原始图谱作者，引用请使用 `atlas.doi`。
