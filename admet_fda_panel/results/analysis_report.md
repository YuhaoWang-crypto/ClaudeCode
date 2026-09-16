# ADMET 与类药性分析报告 — 30 个 FDA 已批准药物面板

生成日期: 2026-09-16  ·  化合物数: 30  ·  预测终点: 41 (ADMET-AI) + 17 理化/规则描述符 + 5 个结构警示目录

## 1. 方法概要

- **标准化**: RDKit `MolStandardize` (Cleanup → 最大片段 → 去电荷) → 规范 SMILES
- **理化描述符**: MW, cLogP (Crippen), TPSA, HBD/HBA, 可旋转键, 芳环数, sp³ 碳比例, 摩尔折射率, QED
- **类药性规则**: Lipinski 五规则、Veber、Ghose、Egan
- **结构警示**: RDKit `FilterCatalog` — PAINS A/B/C、Brenk、NIH
- **ADMET 预测**: ADMET-AI (Swanson et al. 2024), Chemprop-RDKit 集成模型, 在 TDC ADMET Benchmark Group 41 个数据集上训练
- **参考背景**: ADMET-AI 内置 DrugBank 已批准药物参考集 (2,579 个药物) 的百分位数

## 2. 理化性质与类药性

| 性质 | 中位数 | 均值 ± SD | 范围 | 参考阈值 |
|---|---|---|---|---|
| 分子量 (Da) | 364.5 | 363.8 ± 123.7 | 129.2 – 645.3 | ≤ 500 |
| cLogP | 3.47 | 3.50 ± 1.93 | -1.03 – 6.94 | ≤ 5 |
| 拓扑极性表面积 (Å²) | 64.1 | 61.8 ± 26.8 | 12.0 – 113.4 | ≤ 140 |
| 氢键供体 | 1.0 | 1.1 ± 1.1 | 0.0 – 4.0 | ≤ 5 |
| 氢键受体 | 4.0 | 4.0 ± 1.8 | 1.0 – 7.0 | ≤ 10 |
| 可旋转键 | 5.5 | 5.4 ± 3.4 | 0.0 – 13.0 | ≤ 10 |
| sp³ 碳比例 | 0.38 | 0.36 ± 0.14 | 0.11 – 0.76 | — |
| QED 类药性 | 0.62 | 0.61 ± 0.21 | 0.16 – 0.89 | > 0.5 为佳 |

| 规则 | 通过数 | 通过率 | 未通过化合物 |
|---|---|---|---|
| Lipinski 五规则 (≤1 项违反) | 28/30 | 93% | Atorvastatin, Amiodarone |
| Veber 规则 (RotB≤10, TPSA≤140) | 27/30 | 90% | Atorvastatin, Amiodarone, Verapamil |
| Ghose 过滤 | 17/30 | 57% | Aspirin, Ibuprofen, Acetaminophen, Caffeine, Atorvastatin, Metformin |
| Egan 卵 (TPSA/logP) | 24/30 | 80% | Caffeine, Atorvastatin, Metformin, Amiodarone, Terfenadine, Tamoxifen |

## 3. 关键 ADMET 终点

| ADMET 终点 | 关注类别 | 中位预测值 | > 0.5 的化合物 | 最高预测化合物 |
|---|---|---|---|---|
| hERG 钾通道阻断 | 心脏 QT 延长风险 | 0.77 | 21/30 (70%) | Imatinib (0.98) |
| 药物性肝损伤 (DILI) | 肝毒性 | 0.64 | 18/30 (60%) | Caffeine (0.93) |
| Ames 致突变 | 遗传毒性 | 0.18 | 6/30 (20%) | Omeprazole (0.74) |
| 临床试验毒性失败 | 综合临床毒性 | 0.19 | 5/30 (17%) | Imatinib (0.87) |
| 致癌性 | 长期毒性 | 0.13 | 0/30 (0%) | Atorvastatin (0.49) |
| CYP1A2 抑制 | 药物相互作用 | 0.20 | 8/30 (27%) | Diazepam (0.93) |
| CYP2C9 抑制 | 药物相互作用 | 0.11 | 9/30 (30%) | Ketoconazole (0.87) |
| CYP2C19 抑制 | 药物相互作用 | 0.24 | 11/30 (37%) | Loratadine (0.95) |
| CYP2D6 抑制 | 药物相互作用 | 0.32 | 10/30 (33%) | Haloperidol (0.88) |
| CYP3A4 抑制 | 药物相互作用 | 0.32 | 11/30 (37%) | Ketoconazole (0.99) |
| P-糖蛋白抑制 | 外排转运/DDI | 0.53 | 18/30 (60%) | Ketoconazole (0.98) |
| 血脑屏障穿透 | 中枢暴露 | 0.81 | 24/30 (80%) | Risperidone (0.99) |
| 人肠道吸收 | 口服吸收 | 1.00 | 29/30 (97%) | Diazepam (1.00) |
| 口服生物利用度 | 口服吸收 | 0.82 | 29/30 (97%) | Caffeine (0.97) |

## 4. 连续型 (回归) 终点

| 回归终点 | 中位数 | 四分位距 | 最低 | 最高 |
|---|---|---|---|---|
| 水溶解度 log S (mol/L) | -4.26 | -4.56 – -2.99 | Tamoxifen (-6.59) | Metformin (-0.55) |
| 亲脂性 logD7.4 | 2.56 | 1.18 – 3.22 | Metformin (-1.91) | Tamoxifen (4.17) |
| Caco-2 渗透性 log cm/s | -4.72 | -4.93 – -4.59 | Metformin (-6.72) | Caffeine (-4.40) |
| 血浆蛋白结合率 (%) | 90.8 | 76.3 – 95.2 | Metformin (-3.6) | Tamoxifen (105.5) |
| 稳态分布容积 log L/kg | 3.47 | -1.36 – 8.67 | Rosiglitazone (-5.01) | Amiodarone (43.29) |
| 半衰期 (h) | 28.4 | 4.1 – 49.3 | Aspirin (-19.7) | Amiodarone (374.4) |
| 肝细胞清除率 (µL/min/10⁶) | 32.0 | 24.3 – 63.9 | Metformin (-7.3) | Simvastatin (129.1) |
| 微粒体清除率 (mL/min/g) | 22.89 | 3.83 – 56.20 | Ciprofloxacin (-13.12) | Verapamil (110.65) |
| 急性毒性 LD50 (-log mol/kg) | 2.75 | 2.32 – 3.08 | Acetaminophen (1.82) | Amiodarone (4.17) |

## 5. 高风险化合物

| 药物 | 药理类别 | hERG | DILI | Ames | 结构警示 | 触发标记 |
|---|---|---|---|---|---|---|
| Gefitinib | EGFR kinase inhibitor | 0.96 | 0.89 | 0.61 | 1 | hERG, DILI, Ames, ClinTox |
| Imatinib | BCR-ABL kinase inhibitor | 0.98 | 0.91 | 0.18 | 0 | hERG, DILI, ClinTox |
| Cisapride | Prokinetic (withdrawn) | 0.98 | 0.66 | 0.46 | 2 | hERG, DILI, ClinTox |
| Ketoconazole | Azole antifungal | 0.97 | 0.91 | 0.50 | 1 | hERG, DILI, Ames |
| Sildenafil | PDE5 inhibitor | 0.77 | 0.93 | 0.53 | 0 | hERG, DILI, Ames |
| Omeprazole | Proton-pump inhibitor | 0.76 | 0.92 | 0.74 | 1 | hERG, DILI, Ames |
| Amiodarone | Class III antiarrhythmic | 0.95 | 0.59 | 0.39 | 1 | hERG, DILI |
| Risperidone | Atypical antipsychotic | 0.93 | 0.37 | 0.39 | 0 | hERG, ClinTox |
| Quinidine | Class Ia antiarrhythmic | 0.87 | 0.30 | 0.40 | 1 | hERG, ClinTox |
| Atorvastatin | Statin / lipid-lowering | 0.77 | 0.90 | 0.09 | 0 | hERG, DILI |
| Simvastatin | Statin / lipid-lowering | 0.73 | 0.23 | 0.54 | 0 | hERG, Ames |
| Loratadine | Second-generation antihistamine | 0.68 | 0.75 | 0.08 | 0 | hERG, DILI |

## 6. DrugBank 已批准药物参考百分位数

| 药物 | hERG 阻断 | 肝损伤 | Ames 致突变 | CYP3A4 抑制 | 血脑屏障穿透 | 水溶解度 |
|---|---|---|---|---|---|---|
| Ketoconazole | 0.97 (95%) | 0.91 (85%) | 0.50 (79%) | 0.99 (99%) | 0.59 (40%) | -4.32 (26%) |
| Amiodarone | 0.95 (93%) | 0.59 (60%) | 0.39 (69%) | 0.25 (71%) | 0.95 (79%) | -5.59 (7%) |
| Terfenadine | 0.97 (96%) | 0.06 (12%) | 0.01 (2%) | 0.05 (50%) | 0.87 (65%) | -4.53 (22%) |
| Cisapride | 0.98 (97%) | 0.66 (63%) | 0.46 (76%) | 0.83 (92%) | 0.92 (73%) | -4.27 (27%) |
| Imatinib | 0.98 (98%) | 0.91 (84%) | 0.18 (44%) | 0.79 (91%) | 0.51 (34%) | -4.53 (22%) |
| Atorvastatin | 0.77 (76%) | 0.90 (83%) | 0.09 (24%) | 0.33 (76%) | 0.24 (15%) | -5.00 (14%) |
| Metformin | 0.04 (22%) | 0.00 (2%) | 0.11 (30%) | 0.00 (3%) | 0.40 (25%) | -0.55 (89%) |
| Aspirin | 0.02 (13%) | 0.67 (64%) | 0.08 (23%) | 0.00 (14%) | 0.66 (45%) | -1.62 (76%) |
| Caffeine | 0.05 (24%) | 0.93 (87%) | 0.11 (30%) | 0.00 (27%) | 0.98 (91%) | -0.91 (85%) |
| Ciprofloxacin | 0.09 (30%) | 0.78 (71%) | 0.53 (81%) | 0.00 (24%) | 0.47 (31%) | -3.13 (50%) |

## 7. Demo 验证 — 已知阳性/阴性对照药物

| 终点 (阳性对照药物) | 阳性对照预测值 | 阴性对照均值 | 命中 | 结论 |
|---|---|---|---|---|
| hERG 阻断 / TdP 风险 | Terfenadine 0.97, Cisapride 0.98, Amiodarone 0.95, Haloperidol 0.95, Quinidine 0.87 | 0.07 | 5/5 | 全部命中 |
| CYP3A4 抑制 | Ketoconazole 0.99, Verapamil 0.91 | 0.01 | 2/2 | 全部命中 |
| CYP2D6 抑制 | Quinidine 0.60, Fluoxetine 0.66, Sertraline 0.84 | 0.02 | 3/3 | 全部命中 |
| CYP1A2 抑制 | Ciprofloxacin 0.01 | 0.02 | 0/1 | 未命中 |
| P-糖蛋白抑制 | Verapamil 0.87, Ketoconazole 0.98, Quinidine 0.35, Amiodarone 0.90 | 0.01 | 3/4 | 部分命中 |
| 血脑屏障穿透 | Diazepam 0.98, Caffeine 0.98, Haloperidol 0.98, Fluoxetine 0.98, Risperidone 0.99, Propranolol 0.95 | 0.37 | 6/6 | 全部命中 |
| 芳香化酶 (CYP19A1) | Anastrozole 0.27, Ketoconazole 0.68 | 0.00 | 1/2 | 部分命中 |
| 雌激素受体 | Tamoxifen 0.18 | 0.06 | 0/1 | 未命中 |
| PPAR-γ | Rosiglitazone 0.30 | 0.00 | 0/1 | 未命中 |
| 雄激素受体 | Bicalutamide 0.06 | 0.03 | 0/1 | 未命中 |

DMPK / hERG / BBB 类终点共 21 个阳性对照, 命中 19 个;
Tox21 核受体类终点共 5 个阳性对照, 仅命中 1 个。

## 8. 主要局限

- ADMET-AI 的核受体 (NR-*) 与应激反应 (SR-*) 终点来自 Tox21 qHTS 激动模式定量高通量筛选, 数据高度不平衡; 本面板的验证显示它们**无法**重现已批准药物的既定药理活性 (他莫昔芬/ER、比卡鲁胺/AR、罗格列酮/PPARγ、阿那曲唑/芳香化酶均低于 0.5)。
- hERG 模型在本面板上敏感度高但特异性低 (30 个药物中 21 个 > 0.5), 在全为上市药物的集合上会高估风险。
- 个别回归终点超出物理范围 (他莫昔芬血浆蛋白结合率预测 105.5%), 提示适用域边界。
- 所有预测为计算筛选辅助, 不能替代 GLP 体外/体内安全药理学试验。

## 9. Safety Panel 对比与计算覆盖度分析

### 表 8　Safety-44 / Safety-77 / SafetyScreen87 面板对比

| 对比维度 | Safety-44 (Bowes 2012) | Safety-77 (Brennan / IQ DruSafe 2024) | SafetyScreen87 (Eurofins) |
|---|---|---|---|
| 提出背景 | AstraZeneca、GSK、Novartis、Pfizer 四家公司首次公开各自的体外药理谱系策略，凝练出一个「最小推荐面板」 | IQ 联盟 DruSafe 工作组对 **18 家公司**的问卷调查，总结当代实践并提出扩展的推荐靶点集 | CRO 商业产品：在 Bowes-44 基础上扩充，早于 Safety-77，与之相互独立 |
| 靶点总数 | 44 | 77 | 87 |
| 靶点构成 | GPCR 24、离子通道 8、酶（非激酶）6、转运体 3、核受体 2、激酶 1 | 激酶 20（26%）、非激酶酶类约 12%；**完整清单未公开获取**，其余家族构成本报告不作推断 | Bowes-44 全部 44 个 + 额外 43 个（以 GPCR 与离子通道为主）；GPCR 约 35–44、离子通道约 18–19、核受体约 6、激酶 3 |
| 激酶占比 | 1 / 44 = **2.3%** | 20 / 77 = **26%**（第二大靶点类别） | 3 / 87 = 3.4% |
| 核受体 | 仅 **雄激素受体 AR** 与**糖皮质激素受体 GR** 两个 | 文献描述「增强核受体等代表性不足家族的覆盖」，具体成员未获证实 | AR、GR + **ERα、PPARγ、孕激素受体 PR、RARα** |
| 心脏离子通道粒度 | 通用型结合测定：hERG（膜制备）、L-型钙通道（DHP 位点）、钠通道（site 2，脑制备）、通用 Kv | 扩展离子通道与 GABA-A 功能评估（CRO 实现层面） | 分子层面明确：**Nav1.5**、Cav1.2 三个位点（维拉帕米/DHP/地尔硫䓬）、Cav2.2、Kv1.1、KATP 等 |
| 靶器官覆盖 | 以心血管与中枢神经系统为主 | 明确扩展至心血管与中枢以外的器官系统 | 以 GPCR/离子通道扩充为主，器官覆盖随之加宽 |
| 检测模式 | 以放射性配体**结合**测定为主（单浓度，常用 10 µM，n=2），酶类为酶活测定；原文建议结合测定与功能测定并行以提高灵敏度 | 推荐功能性测定与剂量-反应；CRO 实现中激酶在**生理 1 mM ATP** 下检测 | 以放射性配体结合为主，区分拮抗剂/激动剂放射配体；酶类与激酶为酶活测定 |
| 剂量-反应 | 默认单一浓度筛查，命中后再做浓度梯度确证 | 强调命中后的浓度-反应确证与安全裕度（margin）计算 | 标准产品为单浓度（10 µM），IC50/Ki 需另行订购 |
| 法规对齐 | 非法规强制。ICH S7A 仅将「提示潜在不良作用的配体结合或酶测定数据」列为设计核心组合试验时**应考虑**的因素之一 | 论文摘要指出监管机构「**正日益要求**」提供已知不良反应关联靶点的活性数据 | 厂商宣称「ICH S7A/B 对齐」，实指所筛靶点对应 S7A 关注的器官系统，而非满足任何 S7A/S7B 要求 |
| 代表 CRO 产品 | Eurofins SafetyScreen44™ (P270)、Reaction Biology InVEST44 | Reaction Biology InVEST77 | Eurofins SafetyScreen87 (P342 / PP223) |

### 表 9　Safety panel 靶点与 ADMET-AI 终点的覆盖度映射

| 靶点家族 | 面板靶点 | ADMET-AI 终点 | 在 44 / 87 中 | 本面板阳性药物示例 | 判定 | 说明 |
|---|---|---|---|---|---|---|
| 离子通道 | hERG (Kv11.1) | hERG | 是 / 是 | Imatinib 0.98, Cisapride 0.98, Terfenadine 0.97 等 21 个（对照 5/5） | 可映射且已验证 | **唯一干净映射。**面板为膜制备结合测定，TDC 标签多源于膜片钳 IC50，靶点相同但测定原理不同 |
| 核受体 | 雄激素受体 (AR) | NR-AR | 是 / 是 | 无 (最高: Quinidine 0.13)（对照 0/1） | 仅名义映射，验证失败 | 靶点名义对应，但面板为放射配体**结合**测定，Tox21 NR-AR 为**激动模式**转录报告测定；本报告第 4.2 节验证显示比卡鲁胺预测值仅 0.06，**不可用** |
| 核受体 | 糖皮质激素受体 (GR) | — | 是 / 是 | — | 需体外实验 | ADMET-AI / Tox21 终点集中**无对应终点** |
| 核受体 | 雌激素受体 ERα、PPARγ、PR、RARα | NR-ER | 否 / 是 | 无 (最高: Warfarin 0.37)（对照 0/1） | 仅名义映射，验证失败 | 仅存在于 SafetyScreen87，不在 Bowes-44 中。ER 与 PPARγ 有名义对应的 Tox21 终点，但第 4.2 节验证显示他莫昔芬 0.18、罗格列酮 0.30，**均不可用** |
| GPCR | 5-HT2B（瓣膜病）、5-HT1A/1B/2A、多巴胺 D1/D2S、肾上腺素 α1A/α2A/β1/β2、毒蕈碱 M1/M2/M3、组胺 H1/H2、阿片 µ/δ/κ 等共 24 个 | — | 是 / 是 | — | 需体外实验 | **全部无计算终点。**这是覆盖缺口最大的家族（占 Bowes-44 的 54.5%）。5-HT2B 的风险来自**激动**作用，方向性判断本身即超出当前预测能力 |
| 离子通道 | Cav1.2（L-型钙通道）、Nav1.5、GABA-A(BZD)、NMDA、nAChR α4β2、5-HT3、Kv 等其余 7 个 | — | 是 / 是 | — | 需体外实验 | 无计算终点。CiPA 框架下 Cav1.2 与 Nav1.5 对整合致心律失常风险判断至关重要，而 Bowes-44 本身仅以通用型结合测定覆盖，SafetyScreen87 才细化到分子亚型 |
| 酶（非激酶） | COX-1、COX-2、乙酰胆碱酯酶、MAO-A、PDE3A、PDE4D2 | — | 是 / 是 | — | 需体外实验 | 无计算终点。COX-1（消化道出血）、MAO-A（高血压危象）、PDE3A（心衰死亡率）、PDE4D2（恶心呕吐）均为有明确临床后果的靶点 |
| 转运体 | DAT、NET、SERT | — | 是 / 是 | — | 需体外实验 | 无计算终点。ADMET-AI 的 P-糖蛋白终点**不属于**本面板（见下） |
| 激酶 | Lck（Bowes-44）；Safety-77 扩展至 20 个激酶 | — | 是 / 是 | — | 需体外实验 | 无计算终点。这是 Safety-77 相对 Bowes-44 最大的扩充方向，而计算覆盖度为零 |

### 表 10　不属于次级药理面板的 ADMET 终点

| ADMET-AI 终点 | 本面板阳性药物示例 | 归属说明 |
|---|---|---|
| CYP1A2 / 2C9 / 2C19 / 2D6 / 3A4 抑制 | Ketoconazole 0.99, Anastrozole 0.93, Verapamil 0.91 等 11 个 | **不属于**三个面板中的任何一个。CYP 抑制属于 DMPK / 药物相互作用 (DDI) 学科，依 FDA/EMA 的 DDI 指导原则执行，在 CRO 的业务线上也与安全药理学分属不同部门 |
| P-糖蛋白 (ABCB1) 抑制 | Ketoconazole 0.98, Terfenadine 0.95, Tamoxifen 0.94 等 18 个 | **不属于**三个面板。Bowes-44 的 3 个转运体仅为 DAT/NET/SERT；SafetyScreen87 增加的是 Na⁺/K⁺-ATP 酶、ENT1、GAT-1，均非 ABC 外排转运体。P-gp 属于 DMPK 吸收分布与 DDI 终点 |
| 芳香化酶 (CYP19A1) | Ketoconazole 0.68 | **不在**任何一个面板中。属于 Tox21 内分泌干扰筛查终点 |
| 芳香烃受体 (AhR) | Omeprazole 0.83, Gefitinib 0.82, Imatinib 0.69 | **不在**任何一个面板中。AhR 为 bHLH-PAS 类转录因子，并非经典核受体，属于 Tox21 / CALUX 毒理学终点 |
| BSEP (ABCB11) 抑制 | ADMET-AI 无此终点 | **不在**任何一个面板中，且 ADMET-AI 亦**无对应终点**。胆汁盐外排泵抑制与胆汁淤积型肝损伤相关，需专门的肝转运体测定（BSEP 转染的倒置膜囊泡或三明治培养肝细胞） |

**覆盖度结论**：Bowes-44 的 44 个靶点中，名义可映射 2 个（4.5%），经阳性对照验证可用的仅 1 个（hERG，2.3%）。以 SafetyScreen87 为基准，名义可映射 4/87 ≈ 4.6%，经验证可用仍为 hERG 一个。


### 参考文献

1.　Bowes J, Brown AJ, Hamon J, Jarolimek W, Sridhar A, Waldron G, Whitebread S. Reducing safety-related drug attrition: the use of in vitro pharmacological profiling. Nature Reviews Drug Discovery 2012; 11(12): 909–922. DOI: 10.1038/nrd3845.

2.　Brennan RJ, Jenkinson S, Brown AJ, Delaunois A, Dumotier B, Pannirselvam M, Rao M, Ribeiro LR, Schmidt F, Sibony A, Timsit Y, Sales VT, Armstrong D, Lagrutta A, Mittlestadt SW, Naven R, Peri R, Roberts S, Vergis JM, Valentin J-P. The state of the art in secondary pharmacology and its impact on the safety of new medicines. Nature Reviews Drug Discovery 2024; 23(7): 525–545. DOI: 10.1038/s41573-024-00942-3.

3.　Maciag M, Karamyan VT. Enzymes in secondary pharmacology screening panels: is there room for improvement? Nature Reviews Drug Discovery 2025; 24(6): 480–481. DOI: 10.1038/s41573-025-01173-w.

4.　Papoian T, Chiu H-J, Elayan I, et al. Secondary pharmacology data to assess potential off-target activity of new drugs: a regulatory perspective. Nature Reviews Drug Discovery 2015; 14: 294. DOI: 10.1038/nrd3845-c1.

5.　Eurofins Pharma Discovery Services. SafetyScreen44™ Panel product flyer. Ref. P270, Lit. No. EPDSFL420JUNE16, © Eurofins Cerep S.A. 2016. （本报告表 9 中 Bowes-44 的逐条靶点依据此实现版本）

6.　Eurofins Discovery. SafetyScreen87 Panel 产品目录页 (P342 / PP223)。

7.　ICH Harmonised Tripartite Guideline S7A. Safety Pharmacology Studies for Human Pharmaceuticals. ICH Step 4, 8 November 2000; EMA CPMP/ICH/539/00.

8.　ICH Harmonised Tripartite Guideline S7B. The Non-Clinical Evaluation of the Potential for Delayed Ventricular Repolarization (QT Interval Prolongation) by Human Pharmaceuticals. ICH Step 4, May 2005; EMA CHMP/ICH/423/02.

9.　ICH E14/S7B Implementation Working Group. Clinical and Nonclinical Evaluation of QT/QTc Interval Prolongation and Proarrhythmic Potential — Questions and Answers. Step 4, 21 February 2022. (Federal Register 87 FR 52716, 29 August 2022)

10.　Kenna JG, et al. Can Bile Salt Export Pump Inhibition Testing in Drug Discovery and Development Reduce Liver Injury Risk? An International Transporter Consortium Perspective. Clinical Pharmacology &amp; Therapeutics 2018. DOI: 10.1002/cpt.1222.

11.　Swanson K, Walther P, Leitz J, Mukherjee S, Wu JC, Shivnaraine RV, Zou J. ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. Bioinformatics 2024; 40(7): btae416. DOI: 10.1093/bioinformatics/btae416.

12.　Huang K, Fu T, Gao W, et al. Therapeutics Data Commons: machine learning datasets and tasks for drug discovery and development. NeurIPS Datasets and Benchmarks 2021.

13.　Veith H, Southall N, Huang R, et al. Comprehensive characterization of cytochrome P450 isozyme selectivity across chemical libraries. Nature Biotechnology 2009; 27(11): 1050–1055.

14.　Broccatelli F, Carosati E, Neri A, et al. A novel approach for predicting P-glycoprotein (ABCB1) inhibition using molecular interaction fields. Journal of Medicinal Chemistry 2011; 54(6): 1740–1751.

15.　Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the chemical beauty of drugs. Nature Chemistry 2012; 4: 90–98.

16.　Baell JB, Holloway GA. New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries. Journal of Medicinal Chemistry 2010; 53(7): 2719–2740.

17.　Comprehensive In Vitro Proarrhythmia Assay (CiPA) initiative. cipaproject.org.
