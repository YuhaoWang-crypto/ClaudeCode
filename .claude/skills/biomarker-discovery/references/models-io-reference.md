# 模型/工具参考：用途 · 输入 · 输出 + 新适应症迁移手册

> 汇总本项目采用的每个模型/工具：它做什么、要什么输入、给什么输出；末尾给出把整套流程搬到**任意新适应症**的通用操作步骤。

---

## 1. 基础模型层（生成候选）

### 1.1 Geneformer（单细胞转录组基础模型）★核心
| 项 | 内容 |
|---|---|
| **用途** | 单细胞上下文嵌入；**de novo 细胞态发现**、疾病/应答细胞态分型、**in-silico 扰动**定位驱动基因（非旁观者）、网络中心性分析 |
| **输入** | scRNA-seq 计数矩阵（细胞×基因，**Ensembl 基因 ID** + 每细胞 `n_counts`）；经 rank-value 编码（按 归一化表达/基因中位数 降序取 token，截断至 2048）。可选：细胞标签（condition / responder）做有监督探针 |
| **输出** | 每细胞嵌入向量（V1-10M=256维）；聚类/细胞态标签；in-silico 扰动的"状态位移分"；（微调后）分类预测 |
| **本项目用法** | `geneformer_denovo_discovery.py`（管线）、`modal_geneformer.py`（GPU 真跑，V1-10M + gc30M 词表）。用途 C 分型、A/B 应答态 |
| **落地要点** | 词表必须与模型匹配（V1↔gc30M，V2↔gc104M）；按供体/中心分层防批次效应；**需真实 scRNA 才有真实信号**（合成随机基因无共表达结构，见 `geneformer-gpu-run-results.md`） |

### 1.2 ESM-2 / AMPLIFY（蛋白语言模型）
| 项 | 内容 |
|---|---|
| **用途** | 变异效应预测（VEP，零样本）、蛋白功能/家族嵌入、表位/抗原性辅助 |
| **输入** | 蛋白氨基酸序列（野生型 + 变异型），或仅序列 |
| **输出** | 每残基/每序列嵌入；**变异有害性分**（wt vs mut 的掩码 LM 伪对数似然差 Δpll） |
| **本项目用法** | 对候选基因的**编码错义变异**打分（IL23R R381Q、TYK2 P1104A/A928V/I684S 等，见 `proto_variants_il23r_tyk2.tsv`） |

---

## 2. PROTO（核酸/调控层平台）★
> PROTO = 生物学"高级编程语言"，把多个核酸模型编排成"设计+打分"流程。在 biomarker 挖掘里角色是**候选基因的调控/剪接变异功能验证**（判断"基因型层的通路激活度"）。

| 子模型 | 用途 | 输入 | 输出 |
|---|---|---|---|
| **Evo2** | 基因组（DNA）FM；非编码/调控变异效应、序列打分/生成 | DNA 序列窗口（ref/alt，变异居中） | 序列似然/调控影响分 |
| **AlphaGenome / Enformer / Borzoi** | 序列→功能：表达、染色质可及性预测 | DNA 序列窗口（约以变异/位点为中心） | 预测调控轨迹（表达/可及性…）；ref vs alt 差值=变异效应 |
| **SpliceAI / Pangolin** | 剪接改变预测 | DNA 序列 + 变异位置 | 剪接改变概率（供/受体位点增益/丢失 Δ分） |
| **ViennaRNA** | RNA 二级结构 | RNA 序列 | 最小自由能结构、碱基配对概率 |
| **PROTO 整体** | 设计/打分核酸序列（含启动子、调控元件） | 序列（DNA/RNA）、变异、约束 | 每变异功能影响分（表达/染色质/剪接）；优化后的序列 |

**本项目用法**：非编码/内含子/UTR 变异（IL23R rs7517847/rs10889677 等）→ Evo2+AlphaGenome 打表达/染色质，SpliceAI 打剪接；见 `ibd-il23-jak-gap-mining.md` §3。平台文档：`https://proto.evodesign.org`。

---

## 3. 机制验证与配对层

| 工具 | 用途 | 输入 | 输出 |
|---|---|---|---|
| **Boltz** | 结构 + 结合亲和力共折叠；变异后果/耐药机制背书 | 蛋白序列(+配体/伙伴)，可选变异 | 3D 结构、结合亲和力/置信度；ref vs mut 结构改变 |
| **EDEN**（Basecamp） | 免疫原性预测、抗菌肽生成 | 肽/蛋白序列 | 免疫原性分；生成的 AMP 序列 |
| **Inductive Bio** | 小分子性质（ADMET 类）预测 | 分子（SMILES） | 预测性质（吸收/代谢/毒性等） |
| **ChEMBL** | 靶点/生物活性/药物/机制数据库；可成药性、CDx 配对 | 化合物/靶点/适应症查询 | 靶点、活性（IC50/Ki）、机制、药物、ADMET |

---

## 4. 证据与临床数据层（接地、新颖性、竞争基线）

| 工具 | 用途 | 输入 | 输出 |
|---|---|---|---|
| **PubMed** | 已知 biomarker、机制、新颖性核查、证据三角 | 检索词 / PMID | 文章元数据 + 摘要 + **DOI**（须署名引用） |
| **bioRxiv/medRxiv** | 最新预印本、趋势 | 日期/领域/DOI | 预印本元数据、全文 |
| **ClinicalTrials.gov** | 竞争基线、入组标准、终点、验证队列 | 适应症/干预/申办方 | 试验（设计/终点/入组/NCT） |

---

## 5. 编排与算力

| 组件 | 用途 | 输入 | 输出 |
|---|---|---|---|
| **Claude** | 编排全流程、证据三角、候选排序、生成证据卡片/评分卡/代码 | 上述所有 | 结构化候选 panel、`candidates.yaml`、评分卡、报告、脚本 |
| **Modal** | 无服务器 GPU（跑 Geneformer 等） | Python app + 数据 | 运行结果（嵌入/报告） |

---

## 6. 数据如何在模型间流动（一张图）

```
       临床问题(BEST类别) + 对比组
                 │
   ┌─────────────┼───────────────┬───────────────┐
   ▼             ▼               ▼               ▼
基因组变异     scRNA计数        蛋白序列         已知证据
   │             │               │               │
ESM-2 VEP     Geneformer       ESM-2/Boltz    PubMed/Trials
(Δpll)       (细胞态/嵌入)     (功能/结构)     (基线/新颖性)
   │             │               │               │
   └── PROTO ────┘               │               │
   (Evo2/AlphaGenome/SpliceAI    │               │
    非编码变异功能分)             │               │
                 └───────┬───────┴───────┬───────┘
                         ▼               ▼
                   多模态交叉印证 →  Claude 证据三角
                         │
                         ▼
              组合评分卡 (每机制一个应答概率 + 疾病负荷趋势)
```

---

## 7. 迁移到新适应症：通用操作手册

整套流程与代码**与适应症无关**，换 6 样东西即可复用。

### 7.1 六步流程（同 IBD）
1. **定义问题**：选适应症 + 1–2 个 BEST 类别（Risk/Diagnostic/Predictive/Response/Monitoring/Prognostic）+ **明确对比组**（responder vs NR、progressor vs stable、disease vs healthy）。
2. **证据接地**：PubMed/bioRxiv/Trials/ChEMBL → 已知 biomarker、竞争基线、新颖性缺口。
3. **候选生成**：Geneformer（scRNA 细胞态）+ ESM-2 VEP（编码变异）+ PROTO（调控变异）。
4. **机制过滤**：Boltz/EDEN/ChEMBL。
5. **临床筛**：新颖性 / 可及性 / 临床对齐三门槛。
6. **组合评分卡** + 前瞻验证。

### 7.2 换掉的 6 样东西（其余代码不动）
| # | 换什么 | 从哪找 | 对应本项目文件 |
|---|---|---|---|
| 1 | **适应症 + 对比组** | 临床问题 | 各 panel 的 meta |
| 2 | **scRNA/组学数据集** | GEO、Broad Single Cell Portal、**CELLxGENE**、HCA | `scRNA_datasets.md` |
| 3 | **机制/药物类别**（predictive/CDx 用） | 该病在研药 → Trials/ChEMBL | `candidates.yaml` mechanism_classes |
| 4 | **疾病特异基因模块**（"炎症模块"→本病模块） | 文献 + de novo 发现 | `geneformer_denovo_discovery.py` 的 `INFLAMMATION_MODULE` |
| 5 | **候选变异清单**（该病 GWAS 基因） | GWAS Catalog / Open Targets | `proto_variants_*.tsv` |
| 6 | **验证队列** | 该适应症的 Trials（NCT） | 各 panel validation_cohort |

### 7.3 复用的基础设施（零改动或小改）
- `geneformer_denovo_discovery.py`（改 `INFLAMMATION_MODULE` + 喂新 h5ad）
- `modal_geneformer.py`（改数据源；模型/词表不变）
- `candidates.yaml` schema + `score_candidates.py`（结构通用）
- `report.html` 版式；PubMed/Trials/ChEMBL 检索流程

### 7.4 两个迁移示例（示意映射）
| 环节 | 例A：NSCLC 免疫治疗 | 例B：MASH/NASH（肝） |
|---|---|---|
| BEST 重心 | Predictive/CDx（超越 PD-L1/TMB） | Diagnostic + Prognostic（纤维化进展） |
| 对比组 | ICI 应答 vs 原发耐药 | 进展性纤维化 vs 稳定 |
| scRNA 数据 | 肿瘤/免疫 scRNA（GEO/CELLxGENE） | 肝 scRNA 图谱 |
| 疾病模块 | T 细胞耗竭/抗原呈递 | 星状细胞活化/纤维化程序 |
| 变异清单 | HLA/JAK 通路 | PNPLA3/TM6SF2/HSD17B13 |
| 机制类别 | anti-PD1 / anti-CTLA4 / LAG3 | FGF21 / THRβ / GLP1 等 |
| 验证队列 | ICI 试验 NCT | MASH 试验 NCT |

> 核心不变：**按 BEST 定义问题 → 三条组学腿并行生成候选 → Geneformer 扰动 + Boltz/EDEN 机制过滤 → Trials/PubMed 新颖性与临床筛 → 多模态交叉印证组合成评分卡**。多模态交叉印证是降假阳性的通用关键。

---

### 署名
文献/数据引用来自 PubMed 与 ClinicalTrials.gov（见各配套文件的 DOI/NCT）。模型：`ctheodoris/Geneformer`、ESM-2、PROTO（proto.evodesign.org，集成 Evo2/AlphaGenome 等）、Boltz、EDEN、ChEMBL。
