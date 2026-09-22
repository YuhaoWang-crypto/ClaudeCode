# FDA IVD marker database

按 marker（检测对象）组织的 FDA 体外诊断产品数据库：目录层来自上传的 `FDA_Marker_Catalog_CN.xlsx`（1,206 个 marker 标签、16,642 条 marker→申报号关联），元数据层来自 openFDA，证据层来自 FDA 510(k) 决策摘要的规则解析，人工层来自 `docs/FDA_510k_IVD_cutoff_atlas.md` 的 25 个靶点。

## 文件

| 路径 | 内容 |
|---|---|
| `fda_ivd_markers.sqlite` | 数据库本体（SQLite，WAL 关闭后单文件） |
| `export/*.csv` | 每张表一份 CSV |
| `export/bundle.json` | 浏览页面用的压缩 JSON |
| `../docs/FDA_510k_marker_db.html` | 单文件浏览页面（筛选 / 搜索 / marker 详情 / 决策摘要字段 / 人工层） |
| `source/FDA_Marker_Catalog_CN.xlsx` | 用户提供的目录工作簿（v1，2026-09-22 快照） |

## 表结构（四层，`source` / `confidence` 字段标明来源）

| 表 | 层 | 说明 |
|---|---|---|
| `marker` | 目录 | 1,206 个标签：FDA 专业分类、层级、CLIA analyte ID、关联计数 |
| `marker_submission` | 目录 | marker ↔ 申报号（510(k) / De Novo / PMA），含限定词 |
| `denovo_catalog` `pma_catalog` `cdx` `nat_catalog` `category_overview` `method_notes` | 目录 | 工作簿的补充表 |
| `submission` | 元数据 | 每个 K/DEN/P 号一行；openFDA 补齐 decision_date、applicant、device_name、product_code、regulation；`in_catalog=0` 表示目录未收录、按产品代码从 openFDA 补入（目录的 510(k) 关联止于 2014 年） |
| `product_code` | 元数据 | openFDA classification：名称、21 CFR、class、panel |
| `document` | 证据 | 已抓取的决策摘要 / 510(k) summary（URL、字符数、文本路径） |
| `extraction` | 证据 | 决策摘要按标题拆出的字段（`confidence='section'`）与正则推导字段（`confidence='regex'`）：measurand、intended_use、indications、specimen_types、assay_cutoff、clinical_cutoff、cutoff_numbers、reference_range、standards、clsi_codes、precision、detection_limit、traceability、method_comparison、clinical_studies、clinical_sensitivity/specificity、sens_pct/spec_pct、sample_n、predicate、conclusion 等 |
| `curated_target` `curated_section` | 人工 | 25 个靶点的 10 节整理、cutoff 原型、代表性 cutoff、映射到的 marker_id 与产品代码 |

## 流水线

```bash
pip install openpyxl pypdf fonttools markdown
python3 db/build_db.py db/source/FDA_Marker_Catalog_CN.xlsx          # 1. 目录 → SQLite
python3 db/openfda_enrich.py db/fda_ivd_markers.sqlite               # 2. openFDA 元数据 + 按产品代码补 2003 年后全部 510(k)
python3 db/fetch_ds.py db/fda_ivd_markers.sqlite txt 2003 0.8 DESC   # 3. 限速抓取决策摘要（可多 worker：ASC / DESC / 偏移）
python3 db/parse_all.py db/fda_ivd_markers.sqlite txt                # 4. 规则解析 → extraction
python3 db/load_curated.py db/fda_ivd_markers.sqlite                 # 5. 载入 25 靶点人工层
python3 db/export.py db/fda_ivd_markers.sqlite db/export             # 6. CSV + JSON
python3 db/build_browser.py db/export/bundle.json docs/FDA_510k_marker_db.html   # 7. 浏览页面
```

## 边界与注意

- 决策摘要（`cdrh_docs/reviews/K*.pdf`）约从 2003 年起才有，且并非每个 K 号都有；1975–2002 年的申报只有目录层和 openFDA 元数据。
- `extraction` 是**规则解析**：按 OIVD 决策摘要的两种标题模板（2003–2018 的 A–R 模板、2018 后的 I–IX 模板）切出原文段落，推导字段（样本类型、CLSI 编号、cutoff 数值、灵敏度/特异度百分比、n）用正则从段落里抓，可能漏抓或抓到相邻数字，用前请回到原文段落核对。它不是人工整理，与 25 靶点人工层的可靠性不同。
- PMA 与 De Novo 只有目录层与元数据层（De Novo 决策摘要已抓取，但模板不同，字段解析覆盖有限）；PMA 的 SSED 未抓取。
- marker ↔ 申报号的关联来自目录（CLIA 文档号），目录说明里已写明这是"发现层"而非注册结论；同一 K 号可能挂在多个 marker 下（多参数系统）。浏览页按产品代码补入 2015 年后申报时，只使用该 marker 目录关联里出现 ≥3 次或占比 ≥25% 的产品代码，避免把多参数系统的无关代码带进来。
- openFDA 快照 2026-09-14；决策摘要抓取 2026-09-22。
