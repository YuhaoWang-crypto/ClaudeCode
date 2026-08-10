# 公开数据清单：MRI × CD8/TIL 标签究竟能配对到多少例

本文所有数字都是**实测**的，不是文献转述。测量脚本是 `scripts/refresh_inventory.py`，
它直接打 TCIA/NBIA、GDC 和 NCBI E-utilities 三个公开 API，把 `mri_cd8_til/cohorts.py`
里记录的每一个数字重新算一遍并报告漂移。最近一次实测：**2026-08-10**，结果见
`results/cohort_inventory.json`。

```
=== TCIA collection counts ===
  TCGA-BRCA                    claimed=137    measured=137    ok
  ISPY1                        claimed=222    measured=222    ok
  ACRIN-6698                   claimed=385    measured=385    ok
  Duke-Breast-Cancer-MRI       claimed=922    measured=922    ok
  Advanced-MRI-Breast-Lesions  claimed=632    measured=632    ok
  Breast-MRI-NACT-Pilot        claimed=64     measured=64     ok
  QIN-BREAST                   claimed=41     measured=41     ok

=== TCGA-BRCA pairing (tier A) ===
  MRI & RNA-seq & diag WSI     claimed=136    measured=136    ok

=== I-SPY2 pairing (tier B) ===
  ISPY2 MRI & expression       claimed=717    measured=717    ok
```

---

## 一句话结论

**乳腺 MRI 不缺 —— 公开可得 2,700–3,100 例。缺的是标签。**

配对了 CD8 IHC 空间表型的公开 MRI 数据，实际上是 **0 例**。所以多中心方案必须换标签，
而换标签这件事本身有代价，代价必须写清楚。

---

## 三层数据

标签质量决定了一个队列能干什么，所以按"标签层级"而不是按样本量组织。

### Tier A：有空间分辨的 TIL 标签 —— TCGA-BRCA，n = 136

这是整个方案的关键，因为它是**唯一**能同时拿到影像和空间免疫表型的公开队列。

| 项目 | 实测值 |
|---|---|
| TCIA 中有 MR 的患者 | **137** |
| 其中同时有 RNA-seq **和**诊断级 FFPE 全切片 | **136** |
| GDC 中 TCGA-BRCA 总病例 | 1,098（1,095 有 RNA-seq，1,062 有诊断切片） |
| 扫描仪型号 | GE SIGNA EXCITE(939 series) / SIEMENS Sonata(234) / GE SIGNA HDx(132) / Philips Achieva(120) |
| 送检机构 | MSKCC、Mayo Clinic、UPMC、Roswell Park 等 TCGA tissue source sites |
| ID 对齐 | TCIA `PatientID` == GDC `submitter_id`，例如 `TCGA-AO-A03M`，**直接可连** |
| 获取门槛 | 完全开放，无需 DUA |

**标签来源**：Saltz 等 (Cell Reports 2018) 的 TIL 图谱，TCIA 编号 `TIL-WSI-TCGA`，
DOI `10.7937/K9/TCIA.2018.Y75F9W1`，CC BY 3.0，覆盖 13 种 TCGA 肿瘤共 4,759 例，
以 patch 级 CSV 发布（全量 73.4 GB，可按 BRCA 子集下载）。这是 patch 级空间图谱，
所以能算"肿瘤中心密度 / 浸润前沿密度"，进而落到沙漠/排斥/炎症三分类 ——
实现见 `mri_cd8_til/labels/spatial.py`。

**必须写进论文的三条局限**：

1. H&E 推断的 TIL 是**全部淋巴细胞**，不是 CD8 特异的。要做 CD8 特异，只能再叠
   RNA-seq 的 CD8A/CD8B 与去卷积分数做锚定（`labels/transcriptome.py`），
   或者接受"TIL 空间表型"这个稍弱的终点。
2. n=136 是**验证队列的量级，不是训练队列的量级**。它的正确用途是当原文缺失的外部验证集。
3. 采集协议本身高度异质（1.5T/3T 混杂、四家厂商）。这既是它的价值（真实的域偏移），
   也意味着不做 harmonization 就直接用会得到一个"扫描仪分类器"。

### Tier B：有定量但无空间信息的免疫标签 —— I-SPY2，n = 717

| 项目 | 实测值 |
|---|---|
| TCIA `ISPY2` 有 MRI 的患者 | **719** |
| GEO `GSE194040` 治疗前表达谱样本 | **987**（Agilent 44K，两平台已用 ComBat 合并） |
| **交集** | **717** |
| ID 对齐 | TCIA `ISPY2-100899` ↔ GEO `!Sample_title = ISPY2_100899`，**1:1 直接可连** |
| 中心数 | I-SPY2 TRIAL 约 20 家美国中心 |
| 附带变量 | HR / HER2 / MammaPrint / **pCR** / 治疗臂（含 pembrolizumab 臂 n=69） |

这是样本量上的解法：一次性把 n 从 182 抬到 717，且天然多中心，还带一个免疫检查点治疗臂。

**但有一条硬约束，不能绕过**：转录组标签**只能支持"炎症型 vs 其余"这一个对比**。
沙漠型和排斥型的肿瘤核心都没有浸润，两者的差别**纯粹是空间排布**；一份 core biopsy
的 bulk 表达谱没有任何机制可以表达这个差别。

也就是说：

| 原文模型 | 对比 | Tier B 能否验证 |
|---|---|---|
| RM-whole | 炎症型 vs（沙漠+排斥） | ✅ 可以 |
| RM-peri | 沙漠型 vs 排斥型 | ❌ **不行** |

用转录组标签去报一个 RM-peri 的 AUC，等于在测量一个和模型声称测量的东西不同的量 ——
这正是本项目要避免的那一类错误。代码层面在 `labels/schema.py` 里把这条约束写死了。

另一条局限：表达谱来自治疗前空芯针穿刺，MRI 拍的是整个肿瘤。**取样错配是真实且不小的噪声源**。

### Tier C：只有影像 —— 用于分割、协调与域偏移压力测试

| 队列 | 患者数 | 用途 |
|---|---|---|
| Duke-Breast-Cancer-MRI | **922** | 单机构、6 种以上扫描仪 —— 测"同一家医院内扫描仪批次效应"的最佳数据 |
| Advanced-MRI-Breast-Lesions | **632** | 单一机型（GE Signa HDxt），做同质对照组 |
| ACRIN-6698 | **385** | I-SPY2 的 DWI 子研究，ADC 图 |
| ISPY1 / ACRIN 6657 | **222** | 9 家 ACRIN 中心 |
| Breast-MRI-NACT-Pilot | **64** | 恰好与原文新辅助队列同样大小，天然替身 |
| QIN-BREAST | 41 | 小样本多模态 |
| MAMA-MIA | 1,506 | **专家分割**再标注（源自上面四个队列，**不可与它们相加**）+ 预训练 nnU-Net |

MAMA-MIA（*Scientific Data* 2025，`github.com/LidiaGarrucho/MAMA-MIA`）的价值不是提供标签，
而是**消除分割瓶颈**：1,506 例专家勾画 + 预训练模型 + 49 个协调过的临床变量。
原文的 VOI 是一位放射科医师手工勾的，多中心复现最大的人力成本就在这里。

**公开乳腺 MRI 总量：2,737 – 3,122 例。**

给区间而不是给一个确数，是因为 ACRIN-6698 的 385 例是 I-SPY2 的 DWI 子研究 ——
大概率就是同一批人，但两个 collection 在 TCIA 上分别去标识化，ID 对不上（实测交集为 0），
所以既不能证明重复也不能证明不重复。下界排除它，上界包含它
（`total_public_mri()` 与 `total_public_mri(count_overlapping=True)`）。
MAMA-MIA 的 1,506 例在两种口径下都不计入 —— 它是再标注，不是新病人。

---

## 两条查证过的**断链**（省得别人再踩一遍）

这两条是实测发现的负面结果，比正面结果更省人时间：

1. **ISPY1 影像 ↔ GSE22226 表达谱：连不上。**
   TCIA 用 `ISPY1_1001`（4 位），GSE22226 用 6 位研究号（如 `226342`）。
   按数字朴素匹配交集为 **0**。需要 ACRIN 的 ID 对照表才能打通。
   打通后可再增加约 150 例多中心 tier-B 样本。

2. **ACRIN-6698 (DWI) ↔ ISPY2 (DCE)：连不上。**
   `ACRIN-6698-102212` 与 `ISPY2-100899` 同为 6 位，但交集实测为 **0** ——
   两个 collection 在 TCIA 上是**分别去标识化**的。
   打通后可给 717 例加上 ADC，而 ADC 对淋巴细胞浸润是有真实生物学信号的。

---

## 汇总

| | Tier A (TCGA-BRCA) | Tier B (I-SPY2) | 原文 |
|---|---|---|---|
| n | 136 | 717 | 182 (+64) |
| 中心数 | ≥4 机构 / 4 厂商 | ~20 中心 | **1** |
| 标签空间分辨 | ✅ | ❌ | ✅ |
| 标签 CD8 特异 | ⚠️ 需 RNA-seq 锚定 | ✅ | ✅ |
| 支持 RM-whole | ✅ | ✅ | ✅ |
| 支持 RM-peri | ✅ | ❌ | ✅ |
| 外部验证 | ✅ **这就是外部验证集** | ✅ | ❌ |
| 获取成本 | 开放 | 开放 | 需原作者 |
