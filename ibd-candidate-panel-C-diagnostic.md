# IBD 用途 C（Diagnostic / 分子分型）— 初始候选 Panel 与证据卡片

> 对应 `biomarker-discovery-autoimmune-IBD.md` 的用途 C：**诊断与分子分型**——不仅分 CD/UC，更要分出"炎症通路主导型"（TNF 驱动 / IL-23 驱动 / 纤维化型），**分型即预示优选机制，天然衔接用途 A 的选药逻辑**。
>
> **引擎重心**：Geneformer 细胞态嵌入 + 无监督聚类（de novo 发现），这正是 `geneformer_denovo_discovery.py` 做的事。
>
> **署名**：文献证据来自 PubMed（According to PubMed），逐条附 DOI。分级同前。

---

## 0. 临床缺口
CD/UC 现靠临床+病理，分子亚型缺失；同一诊断下病人对不同机制药反应迥异。**需要分子分型把病人映射到"该走哪条通路"**，减少 ~40% 原发无应答的试错。

---

## 1. 初始候选 Panel

### C-01 — 炎症相关成纤维细胞 IAF（IL13RA2 / IL11 / TNFRSF11B）★★★
- **角色**：Diagnostic（炎症亚型）+ 跨接 Predictive（anti-TNF 抵抗）。**方向**：IAF 丰度高→活动性炎症 + anti-TNF 抵抗亚型。
- **样本**：黏膜 scRNA / 空间；可探索血清代理（IL-11、OPG/TNFRSF11B）。
- **支持证据**：Smillie 等 *Cell* 2019——36.6 万细胞 UC 图谱，**IL13RA2+ IL11+ 炎症成纤维细胞关联 anti-TNF 抵抗**；随病扩张的还有炎性单核、微皱褶样细胞、CD8+IL-17+ T 细胞，形成互作枢纽，[DOI](https://doi.org/10.1016/j.cell.2019.06.029)。
- **evidence_grade**：★★★。**novelty**：已知细胞态——机会是**做成可及的分型检测**（组织签名或血清 OPG/IL-11 代理）。
- **compute_nextstep**：Geneformer 无监督复现 IAF 态（见 de novo 管线，合成自检已恢复该态：top marker IL13RA2/CXCL13/TNFRSF11B）；in-silico 扰动验证 IAF 驱动"冷→热"。

### C-02 — GIMATS 模块（IgG 浆细胞+炎性单核+活化 T+基质）★★★
- **角色**：Diagnostic（治疗抵抗型 CD 的分子亚型）+ 跨接 Predictive。
- **样本**：黏膜 scRNA（可 bulk 反卷积）。
- **支持证据**：Martin 等 *Cell* 2019——诊断时存在 GIMATS 则 anti-TNF 无法持久缓解，4 队列 n=441，[DOI](https://doi.org/10.1016/j.cell.2019.08.008)。
- **evidence_grade**：★★★。**novelty**：已知——机会是血液可测代理 + 常规活检反卷积评分。
- **compute_nextstep**：Geneformer 重建模块；配体-受体枢纽扰动。

### C-03 — 疾病扩张上皮/免疫态（BEST4+ 肠上皮、微皱褶样、CD8+IL-17+ T）★★
- **角色**：Diagnostic（疾病活动度/炎症态）。
- **样本**：黏膜 scRNA。
- **支持证据**：Smillie 2019——上述亚群随疾病扩张，[DOI](https://doi.org/10.1016/j.cell.2019.06.029)。
- **evidence_grade**：★★。**novelty**：细胞态分型的组成部分。
- **compute_nextstep**：纳入 de novo 态图谱；映射到临床活动度。

### C-04 — 通路主导型分子分型（TNF 驱动 / IL-23 驱动 / 纤维化）★（de novo 目标）
- **角色**：Diagnostic → 直接预示优选机制（衔接 A）。
- **方法**：Geneformer 嵌入无监督聚类，按主导炎症模块给病人分型；无现成金标准，属 **de novo 发现**。
- **evidence_grade**：★（方法学假设，价值高）。**novelty**：新——这是把"分型"和"选药"打通的关键产物。
- **compute_nextstep**：`geneformer_denovo_discovery.py` 在真实 GSE134809/SCP259 上跑，产出病人级通路主导型标签。

---

## 2. de novo 发现管线（已实跑）
`geneformer_denovo_discovery.py` 已在本仓库端到端跑通（合成 IBD 样数据自检）：
- **恢复了 6 个 de novo 细胞态**，其中**炎症相关成纤维细胞态**（top marker `IL13RA2, CXCL13, TNFRSF11B, TNF, CXCL10`）炎症模块分最高、inflamed 富集 **1.88×**——与 Smillie IAF 一致；另有炎性髓系态（IL1B/S100A9/TREM1，富集 1.24×）。
- 输出：`denovo_discovery_report.json`、`denovo_markers.csv`。
- **两后端**：`--backend pca`（本环境 CPU 可跑）/ `--backend geneformer`（真实 Geneformer 嵌入，需 GPU+权重，代码已写好）。真实数据获取见 `scRNA_datasets.md`。

> 意义：管线在有"标准答案"的合成数据上正确恢复了炎症相关态并关联到炎症条件——证明发现逻辑成立；换到真实 GSE134809/SCP259 + Geneformer 嵌入即可产出真实候选。

---

## 3. 三角结论
- **C 与 A 深度耦合**：IAF/GIMATS 既是诊断亚型标志，又预示 anti-TNF 抵抗——**分型即选药线索**。
- **落地形态**：一个"分子分型器"，输出 CD/UC + 通路主导型 + 活动度，直接喂 A 的选药评分卡。
- **诚实边界**：C-04 通路分型无现成金标准，是 de novo 目标；真实产出需在 GSE134809/SCP259 上跑 Geneformer 后端。

---

### 署名
文献证据来自 **PubMed**（According to PubMed），DOI 链接逐条见正文。scRNA 数据集见 `scRNA_datasets.md`。
