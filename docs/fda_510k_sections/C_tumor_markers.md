# 靶点组 C：肿瘤标志物（21 CFR 866.6010 / 866.6050）——FDA 510(k) 申报信息整理

> 数据来源说明：本文所有申报数据均取自 FDA 公开的 510(k) Decision Summary（`cdrh_docs/reviews/<K>.pdf`）及 510(k) Summary（`cdrh_docs/pdfYY/<K>.pdf`）。凡文件未记载的信息一律标注「文件未载明」；来自背景知识（临床指南/生理学）的内容单独标注为「背景（非申报文件）」。openFDA 累计清关数为 `fda_fetch.py list <CODE>` 首行 `total=` 值（检索日期 2026-09-22）。
>
> **本组共同特点（从各文件预期用途原文归纳）**：
> 1. 预期用途一律限定为 "aid in the management / monitoring" 已确诊癌症患者的 **serial testing**，明确 "not intended for screening or diagnosis"（例外：ROMA 是术前风险分层，但同样声明 "not a screening or stand-alone diagnostic assay"；PSA 的"检测辅助"用途属 PMA）。
> 2. 各文件的 "Assay Cut-Off" 栏几乎均为 "Not applicable"/"See clinical cut-off"，判定逻辑不是单一 cutoff，而是 (a) 健康参考人群上限（ULN，通常 95th/97.5th 百分位）+ (b) **serial change**（% change 或 Reference Change Value, RCV）。
> 3. 临床验证设计高度一致：回顾性/前瞻性 serial serum sets（每人 ≥3 次采血），医生根据影像/体检/病理将每次随访归类为 NED / Stable / Responding / Progression，将 Progression vs. No-Progression 与「显著升高」做 2×2 表，报告 sensitivity / specificity / PPV / NPV / total concordance，并与 predicate 做 PPA/NPA。
> 4. 标签均含 "values obtained with different assay methods cannot be used interchangeably…additional sequential testing to confirm baseline values" 警示。

---

## 癌胚抗原（CEA, Carcinoembryonic Antigen）

### 1. 法规定位
- Product code **DHX**（System, Test, Carcinoembryonic Antigen）；21 CFR **866.6010** Tumor-associated antigen immunological test system；**Class II**；Panel: Immunology (82 / IM)。
- openFDA 该代码累计清关数：**total = 28**（含早期对照品/柱等附属品；近 20 年 IVD 主流平台清关约 10 件）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K223921 | 2023-09-22 | Beckman Coulter / Access CEA（DxI 9000, Lumi-Phos PRO 底物） | 顺磁微粒 CLIA，双单抗一步夹心 | 血清 | aid in the management of cancer patients in whom changing CEA concentrations have been observed | 决策摘要 + 510(k) summary |
| K231517 | 2023-08-23 | Ortho Clinical Diagnostics (QuidelOrtho) / VITROS Immunodiagnostic Products CEA Reagent Pack（VITROS 5600） | 微孔板免疫化学发光（HRP–luminol），biotin 抗体预结合孔以消除 biotin 干扰 | 血清、EDTA/heparin 血浆 | aid in the prognosis and management of cancer patients in whom changing concentrations of CEA are observed | 决策摘要 + 510(k) summary |
| K200215 | 2020-04-13 | Siemens / ADVIA Centaur CEA（XP/XPT） | 直接化学发光（吖啶酯），兔多抗 + 鼠单抗 | 血清、EDTA/Li-heparin 血浆（本次新增血浆） | aid in the management of cancer patients in whom changing concentrations of CEA are observed | 决策摘要 + 510(k) summary |
| K071603 | 2008-06-25 | Dade Behring (Siemens) / Dimension Vista CEA Flex（LOCI） | 均相 LOCI 化学发光夹心 | 血清、Na/Li-heparin 血浆 | aid in the management of cancer patients in whom changing CEA concentrations have been observed | 决策摘要（无 510(k) summary） |

### 3. 预期用途与声明类型
- 声明类型：**监测/管理（monitoring / management）**，Ortho 版本另含 "prognosis"。四份文件均为 **Rx only**（处方用）、中心实验室全自动平台，非 POC。
- 标签限制（K231517 510(k) summary 原文）："not recommended as a screening procedure for cancer detection"；"CEA levels…regardless of level, should not be interpreted as absolute evidence of the presence or absence of malignant disease"。
- 三份近年文件（K223921、K231517、K200215）均为已清关产品的修改（平台迁移/底物更换、biotin 干扰改造、新增血浆基质），临床研究与 cutoff 均 "Refer to" 原始清关（K981985、K041322、K981478）；只有 K071603 是新器械并载有完整临床监测研究。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 参考区间/期望值 | 建立方式 | serial change 定义 |
|---|---|---|---|
| K071603 Vista CEA | 非吸烟者 0.0–3.0 ng/mL；吸烟者 0.0–5.0 ng/mL；健康人群 96.4% < 5.0 ng/mL | 非参数法，健康成人 n=347（非吸烟 198，96.0% ≤3.0；吸烟 149，96.6% ≤5.0）+ 结直肠癌 74 例分布 | **RCV = 36.2%**（predicate Access CEA 36.7%）；公式 RCV = 2^1/2 × Z × (CV_A² + CV_I²)^1/2，CV_I（个体内生物学变异）取文献 12.7%（Ricos 1999），CV_A 取本产品总不精密度；决策摘要明确 "Clinical cut-off: Not applicable for serial monitoring assay that looks for a significant rise" |
| K231517 VITROS CEA | 标签期望值：非吸烟者 n=149 中 91.9% 在 0–3.0、6.0% 在 3.0–5.0；吸烟者 n=101 中 67.3%、22.8%、8.9%（5–10）、1.0%（>10）；总 n=250；另列结直肠癌 n=114（72.8% >10）、乳腺 69、肺 56、卵巢 51、GI 40、良性 GI 42、肝硬化 65、肺 50、肝炎 31 | predicate（K041322）期望值基于 768 份健康+疾病样本；本次按 **CLSI EP28-A3c** 用 68 非吸烟 + 72 吸烟健康者验证（分别 89.7%/72.2% 在 0–3.0；无人 >10） | 文件未载明具体 % change；Clinical cut-off "Refer to K041322" |
| K223921 Access CEA | "Refer to K981985" | — | 决策摘要："No cutoff for CEA monitoring has been recommended" |
| K200215 ADVIA Centaur CEA | "established in K981478" | — | 文件未载明 |

- 分层逻辑：非吸烟 3 ng/mL / 吸烟 5 ng/mL 的分层来自健康人群分布（K071603 用非参数百分位；K231517 直接以分布区间列表给出），并非 ROC 最优点。
- 灰区：Ortho/Vista 期望值表均把 >3.0–5.0、>5.0–10.0 ng/mL 作为分布区间列出，但文件未赋予其"灰区"判定含义。

### 5. 生物学 / 生理学依据
- 申报文件内：吸烟者 CEA 分布右移（K071603、K231517 数据）；非恶性疾病（肝硬化 12.3% >10 ng/mL、肝炎、良性 GI/肺）亦可升高（K231517 表）；交叉反应物 NCA（non-specific cross-reacting antigen）/NCA-2 在 500/100 ng/mL 下无显著交叉（K071603、K231517、K200215）。K071603 引用 NACB tumor marker 实践指南作为参考文件。
- 背景（非申报文件）：CEA 属 CEACAM 家族糖蛋白，主要用于结直肠癌术后复发监测；ASCO/NACB 建议每 3 个月检测；生物学变异约 10–13%（与 K071603 引用值一致）。

### 6. 样本类型与样本要求
- K223921：仅血清；矩阵比对 "Not applicable"。
- K231517：血清 + Li-heparin/K2-EDTA 血浆；矩阵等效 n=40 配对，Deming 斜率 0.998（Li-hep）/0.995（EDTA），r 0.999/0.998；注明本次配方改进（BSA 0.5%→3%、加 Tween/EDTA）改善了此前 EDTA 偏差限制；不推荐浑浊样本。
- K200215：新增 K2-EDTA（n=64 配对，Deming 斜率 0.95，r 1.00）和 Li-heparin（n=46，斜率 0.99）；并做抗凝剂过量滴定（3×/5× 标称浓度，回收 96–105%）。
- K071603：血清 + Li/Na-heparin 血浆，n=54 配对，斜率 1.00/0.99，r 0.997/0.998；样本稳定性 4°C 7 天、−20/−70°C 30 天、冻融（偏差 <7%）。
- 溯源：Ortho 溯源至 WHO 1st IRP 72/225（510(k) summary 处写 73/601，两处不一致，如实记录）；Vista 用 predicate Access CEA 赋值（"No information provided on traceability to any reference standard"）；Siemens Centaur 为内部标准品。

### 7. 分析性能验证所依据的标准
| K 号 | 文件列出的标准 | 对应项目 |
|---|---|---|
| K223921 | CLSI EP05-A3；EP06 2nd Ed；EP09c 3rd Ed；EP17-A2 | 精密度（7 水平 2×3×20 天，n≥120，within-lab CV 2.5–5.2%）；线性 0.1–1042.7 ng/mL；方法比对 n=153 Passing-Bablok；LoB 0.09 / LoD 0.1 / LoQ 0.2 ng/mL（20% CV）；AMI 0.2–1,000 ng/mL |
| K231517 | EP05-A3；EP06-Ed2；EP07-A3；EP09c；EP17-A2；EP25-A；EP28-A3c；EP35；EP37；（510(k) summary 另提 EP21、EP34 hook） | 精密度 3 lot n=240，总 CV 3.1–3.6%；线性 0.22–499.6；干扰（biotin 0.351 mg/dL、HAMA 800 µg/L、RF 900 IU/mL、Hb 1000 mg/dL 等）；hook 至 80,000 ng/mL；稳定性 52 周/8 周在机；参考区间验证；矩阵等效；LoB 0.08/LoD 0.31/LoQ 0.31 |
| K200215 | EP05-A3；EP06-A；EP07 3rd；EP17-A2（summary 提 EP09-A3） | LoB 0.5 / LoD 1.0 / LoQ 2.0 ng/mL（AMI 下限从 0.5 提高到 2.0）；抗凝剂干扰 |
| K071603 | NACB 指南；NCCLS H3-A5（采血）；EP5-A2；EP9-A2；EP7-A2；EP17-A；FDA bundling/leftover specimen 指南 | 精密度（within-lab ≤3.6%）；线性 1.3–1207.8；干扰 + HAMA（3 份 HAMA 样本，偏差 < −8.9%）；LoB 0.12/LoD 0.2；hook 225,000 ng/mL |

### 8. 临床验证设计与结果
- **K071603（唯一含完整监测临床研究）**：75 套回顾性结直肠癌 serial serum sets（≥3 次采血），排除 1 例，74 例（36.1–86.2 岁，均 63 岁；I–IV 期）；医生依据体检/影像（CT/MRI/X 线/结肠镜/超声）/手术病理判定 Active-Progressive / Responding / Stable / NED。逐访视配对 n=217。以 RCV >36.2% 升高判"阳性"：Progression 58 对中 32 阳性，No-Progression 159 对中 26 阳性 → **Sensitivity 55.2% (41.5–68.3)、Specificity 83.6% (77.0–89.0)、Overall agreement 76.0% (69.8–81.6)**；predicate（Access CEA，RCV 36.7%）Sens 55.2%、Spec 79.9%、agreement 73.3%。有效性判据：Sens+Spec > 1（bootstrap 2000 次，95% CI 1.2416–1.5227）。与 predicate 一致性：overall 94.5%、PPA 85.9% (55/64)、NPA 98.0% (150/153)。方法比对：n=141，斜率 1.01、截距 9.01、r 0.989（0.8–974 ng/mL）。
- K223921：方法比对 n=153（0.46–1071 ng/mL），Passing-Bablok 斜率 0.98 (0.97–0.99)、截距 0.06、R 1.00。临床 "Refer to K981985"。
- K231517：修改前后比对 n=110（0.56–396），加权 Deming 斜率 1.01 (0.997–1.012)、截距 0.106、R² 0.999。临床 "Refer to K041322"。
- K200215：Centaur vs ACS:180 重新分析 n=201（2.0–78.9），Deming 斜率 0.97、截距 0.11、r 1.00。

### 9. 厂家间差异与要点
- AMI 差异显著：Access 0.2–1,000；Vista 0.2–1,000；VITROS 0.31–400（自动稀释至 40,000）；Centaur 2.0–100 ng/mL（K200215 将 LoQ 上调至 2.0）。
- 溯源不统一（WHO 72/225 vs 内部标准 vs predicate 赋值），因此各家均强调结果不可互换。
- 只有 Vista（2008）以 RCV 形式给出显著变化定义；Ortho/Beckman/Siemens 近年修改型申报均引用原清关，未在本次文件中载明 % change。
- Ortho 与 Beckman（K223921、K240479、K240927）2022–2024 系列申报的共同主题：biotin 干扰改造（预结合 biotinylated 抗体）和 DxI 9000/Lumi-Phos PRO 平台迁移。

### 10. 来源
- K223921: https://www.accessdata.fda.gov/cdrh_docs/reviews/K223921.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K223921.pdf
- K231517: https://www.accessdata.fda.gov/cdrh_docs/reviews/K231517.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K231517.pdf
- K200215: https://www.accessdata.fda.gov/cdrh_docs/reviews/K200215.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K200215.pdf
- K071603: https://www.accessdata.fda.gov/cdrh_docs/reviews/K071603.pdf

---

## 糖类抗原 125（CA 125）

### 1. 法规定位
- Product code **LTK**（Test, epithelial ovarian tumor-associated antigen (CA 125)）；21 CFR 866.6010；Class II；Panel Immunology (82)。
- openFDA 累计清关数：**total = 26**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K240479 | 2024-05-10 | Beckman Coulter / Access OV Monitor（DxI 9000） | 顺磁微粒 CLIA 两位点夹心 | 血清、heparin 血浆 | aid in the management of ovarian cancer patients; serial testing…in conjunction with other clinical methods | 决策摘要 + summary |
| K221355 | 2022-12-12 | Ortho / VITROS CA 125 II Reagent Pack（5600） | 免疫化学发光（HRP），biotin 抗体预结合 | 血清、EDTA/heparin 血浆 | aid in monitoring response to therapy for patients with epithelial ovarian cancer | 决策摘要 + summary |
| K200199 | 2020-04-06 | Siemens / ADVIA Centaur CA 125II | 直接化学发光（M11 吖啶酯 + OC125 荧光素 / 抗荧光素捕获） | 血清、EDTA/Li-heparin 血浆（新增） | aid in monitoring patients previously treated for ovarian cancer…early detection of cancer recurrence…monitoring progression or regression…not intended for screening or diagnosis | 决策摘要 + summary |
| K143534 | 2015-08-06 | Roche / Elecsys CA 125 II（cobas e 411）+ CalCheck | ECLIA，biotin-M11 + Ru-OC125 | 血清、Li-hep、K2/K3-EDTA、Li-hep 分离胶管 | aid in the detection of residual or recurrent ovarian carcinoma; monitoring…disease progress or response to therapy | 决策摘要 + summary |
| K142895 | 2015-05-21 | Fujirebio / Lumipulse G CA125II + G1200 System | CLEIA（ALP–AMPPD），OC125/M11 | 血清、Na/Li-heparin、K2-EDTA 血浆 | aid in monitoring recurrence or progressive disease in patients with ovarian cancer | 决策摘要 + summary |

### 3. 预期用途与声明类型
- 均为 **监测（recurrence / progression / response to therapy）**，Rx only，非 POC。Siemens 版本最明确："not intended for screening or diagnosis of ovarian cancer"，并建议由"physician trained and experienced in the management of gynecological cancers"开具。
- Roche 版本描述测量物为 "OC 125 reactive determinants…associated with a high molecular weight glycoprotein…of women with primary epithelial invasive ovarian cancer (excluding those with cancer of low malignant potential)"。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 健康参考上限（ULN） | 建立方式 / 人群 | serial change 定义 |
|---|---|---|---|
| K142895 Lumipulse G CA125II | **ULN = 30.8 U/mL**（97.5th 百分位） | 表观健康女性 n=240（120 绝经前 + 120 绝经后，19–69 岁），97.5% 低于 ULN；良性妇科 n=260（86.5% <ULN）、良性非妇科 40、妊娠 40（72.5%）、CHF 40、高血压 40；未治疗卵巢癌 n=105（35.2% <ULN）；决策摘要注明干扰实验特意包含"接近常规正常上限 35 U/mL"的样本 | **≥20% 升高**：由 1.645×√2×总 CV（8.6%）≈ 20% 推导；"cut-off of >20% was chosen to balance sensitivity and specificity"；ROC AUC 0.727 (SE 0.047) |
| K143534 Elecsys CA 125 II | **ULN = 38.1 U/mL**（95th 百分位） | 健康女性 n=240（120/120，18–87 岁，白人 96%），均值 17.2、中位 14.1、5th 6.4；95% <ULN；非卵巢癌 199 例、良性 411 例、其他 80 例（共 690 份）：子宫内膜癌 70% <ULN、良性妇科 76%、妊娠 92% | 文件未载明（临床研究 "Not Applicable"，本次为标签现代化/LoB-LoD-LoQ 补充） |
| K221355 VITROS CA 125 II | **35 U/mL** | predicate 200 名正常女性 98.5% ≤35；本次 EP28-A3c 验证 120 名健康非吸烟者（30 女 <50 岁、30 女 >50 岁、60 男 <50 岁）3 lot，仅 1 例女性 >35 (98.3% <35)；标签分布表：卵巢癌 n=75 中 21 ≤35、40 >100；良性妇科 50 中 48 ≤35；妊娠 30 中 25 ≤35 | **"clinically meaningful change…at least 25% higher…when the VITROS CA 125 result was outside the normal range (>35 U/mL)"**（Refer to K983875） |
| K240479 Access OV Monitor | Refer to K023597 | — | 文件未载明（Refer to K023597） |
| K200199 ADVIA Centaur CA 125II | established in K020828 | — | 文件未载明（Refer to K020828） |

- 关键点：**35 U/mL 并非 ROC 最优，而是各家健康女性 95th–97.5th 百分位（Lumipulse 30.8、Roche 38.1、Ortho 35）**；serial change 阈值由分析不精密度推导（Lumipulse 20%）或 predicate 既定（Ortho 25%）。

### 5. 生物学 / 生理学依据
- 申报文件内：K142895 test principle 段落——CA 125 由 OC125 单抗识别（免疫原为卵巢浆液性囊腺瘤细胞），"可在 >80% 手术证实的上皮性卵巢癌中检出，其水平与临床病程显著相关"；第二代试剂使用 OC125（固相）+ M11（标记）。良性妇科病、妊娠、CHF、肝硬化/浆膜炎均可升高（K143534、K142895、K221355 分布表）。
- 背景（非申报文件）：CA 125 = MUC16 黏蛋白；GCIG 复发定义（CA125 ≥2×ULN 或 ≥2×最低点）与这里 20–25% 的 RCV 型定义不同。

### 6. 样本类型与样本要求
- K143534：血清、Li-heparin、K2/K3-EDTA、Li-heparin 分离胶管；矩阵比对 n=51–52（2.2–2980 U/mL），Passing-Bablok 斜率 0.98–1.00；predicate 允许 Na-heparin/枸橼酸（枸橼酸需 +10% 校正），本次删除；样本稳定 2–8°C 5 天、15–25°C 8 h、−20°C 24 周、2 次冻融。
- K142895：5 种管型（红头、SST、K2-EDTA、Li-hep、Na-hep），45 名供者 ×8 个加标浓度，加权 Deming 斜率 0.99–1.01；2–10°C 10 天、−20°C 10 天、≤6 次冻融；自动稀释仅允许 1:10（1:100 自动稀释未达标，须手工）。
- K221355：Li-hep/K2-EDTA vs 血清 n=49，斜率 0.984/0.990，r 1.00。
- K200199：K2-EDTA n=162（斜率 0.95）、Li-hep n=119（1.02）；抗凝剂 3×/5× 滴定回收 99–104%。
- K240479：血清/heparin 血浆，矩阵 Refer to K023597；样本量 30 µL（原 25）。
- 溯源：Roche 标准化至 Enzymun-Test CA125 II → Fujirebio CA125 II RIA；Fujirebio 溯源至内部参考校准品（对应其 CA125 II RIA）；Ortho 溯源至内部参考物（对应另一商业方法）；**"There are no reference standards for CA 125"**（K143534）。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果摘录 |
|---|---|---|
| K240479 | EP05-A3；EP06 2nd；EP09c；EP17-A2；EP34 | 精密度 8 水平 n=120，within-lab CV 2.6–6.1%；线性 0.5–6,162 U/mL；EMI 1:20 自动稀释至 100,000 U/mL（回收 94–99%）；LoB 0.5/LoD 0.7/LoQ 2.0；方法比对 n=152 斜率 0.98 |
| K221355 | EP05-A3；EP07-A3；EP09c；EP37；EP17-A2；EP06-Ed2；EP28-A3c；EP25-A | 精密度 3 lot n=240 总 CV 1.5–2.4%；线性 3.6–1288；干扰（Hb 750 mg/dL，1000 时 ≥10% 偏差；RF 975 U/mL；总蛋白 11.5 g/dL）；LoB 0.2/LoD 5.5/LoQ 5.5（claimed）；货架期仅 16 周；方法比对 n=146 斜率 1.02 |
| K200199 | EP05-A3；EP06-A；EP07 3rd；EP17-A2 | LoB 2.0/LoD 3.0/LoQ 3.0（AMI 3.0–600）；Centaur vs Immuno-1 n=224 斜率 1.03 |
| K143534 | EP05-A3；EP06-A；EP17-A2 | 精密度 n=80 总 CV 1.3–4.2%；线性 2.0–3000（含血浆样本斜率 0.88–0.89）；hook 50,000 U/mL；LoB 0.6/LoD 1.2/LoQ 2.0；干扰 biotin 70 ng/mL、HAMA 805 µg/mL、RF 1500 IU/mL、Hb 3200 mg/dL；16 常规药 + 23 抗癌药 |
| K142895 | EP05-A2；EP07-A2；C28-A3c；EP06-A；EP09-A3；FDA Tumor Associated Antigen 510(k) 指南；软件指南 | 精密度 20 天总 CV ≤2.6%，三站点 ≤8.4%；线性 2.5–1000；hook >200,000（>20,000 平台，标签警示）；LoB 0.1/LoD 0.5/LoQ 0.5；HAMA 1000 ng/mL、RF 1000、biotin 19.7 mg/dL |

### 8. 临床验证设计与结果
- **K142895（Lumipulse G CA125II）**：59 名卵巢癌女性（16–84 岁，67.8% 绝经后，72.9% 白人），348 份可评估样本、**289 对连续观察**（每人 3–20 次，中位 5.9 次；随访 29–2446 天）。医生依据临床/影像判定进展。以 ≥20% 升高：Progression 52 对中 35 阳性；No-progression 237 对中 57 阳性 → **Sens 67.31%、Spec 75.95%、Total concordance 74.39%、PPV 38.04%、NPV 91.37%**；预设成功标准 Sens+Spec >125%（实测 143.26）。连续比值中位数：NED 0.99、Stable 1.01、Responding 0.91、Progression 1.64。方法比对 vs ADVIA Centaur（predicate）n=102 加权 Deming 斜率 1.13 (1.06–1.20)——超出 0.9–1.1 验收标准，FDA 接受但要求标签警示不可互换。
- **K221355**："Refer to K983875"，标签给出 25% 规则；本次仅方法比对/参考区间验证。
- **K143534**：临床 "Not Applicable"；method comparison n=80（4.7–2680），Passing-Bablok y=0.98x+1.2，R² 0.99；新旧校准品 n=111 斜率 0.99。
- **K240479**：方法比对 n=152（2.4–5001 U/mL），斜率 0.98 (0.97–0.99)、截距 −0.14、R 1.00；临床 Refer to K023597。
- **K200199**：临床 Refer to K020828。

### 9. 厂家间差异与要点
- 抗体对：Fujirebio/Abbott/Roche/Ortho 均用 OC125 + M11（第二代 "II"）；Siemens 用 M11-吖啶酯 + OC125-荧光素 + 抗荧光素固相。
- AMI：Access 2.0–5,000；Lumipulse 2.5–1,000；VITROS 5.5–1,000；Elecsys 2.0–3,000（predicate 0.6–5,000）；Centaur 3.0–600。
- ULN：30.8 / 35 / 38.1 U/mL，取决于百分位选择（97.5th vs 95th）与人群；serial change 20%（Lumipulse，按 1.645×√2×CV 推导）vs 25%（Ortho）。
- 方法间斜率可达 1.13（Lumipulse vs Centaur），是"不可互换"警示的实证。
- Lumipulse 特别注明 1:100 自动稀释失败、>20,000 U/mL 信号平台。

### 10. 来源
- K240479: https://www.accessdata.fda.gov/cdrh_docs/reviews/K240479.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K240479.pdf
- K221355: https://www.accessdata.fda.gov/cdrh_docs/reviews/K221355.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K221355.pdf
- K200199: https://www.accessdata.fda.gov/cdrh_docs/reviews/K200199.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K200199.pdf
- K143534: https://www.accessdata.fda.gov/cdrh_docs/reviews/K143534.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K143534.pdf
- K142895: https://www.accessdata.fda.gov/cdrh_docs/reviews/K142895.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K142895.pdf

---

## 糖类抗原 19-9（CA 19-9）

### 1. 法规定位
- Product code **NIG**（System, Test, Carbohydrate Antigen (CA 19-9) for Monitoring and Management of Pancreatic Cancer）；21 CFR 866.6010；Class II；Panel Immunology (82)。
- openFDA 累计清关数：**total = 12**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K253528 | 2026-08-10 | Siemens / Atellica IM CA 19-9 II | 吖啶酯 CLIA 夹心（1116-NS-19-9 抗体，Fujirebio 授权） | 血清、EDTA/Li-heparin 血浆 | serial measurement…aid in the management of patients diagnosed with cancers of the exocrine pancreas…who have levels…exceeding the upper limit of normal | 决策摘要 + summary |
| K231525 | 2023-08-09 | Ortho / VITROS CA 19-9 Reagent Pack | 免疫化学发光（HRP），biotin 预结合 | 血清、EDTA/heparin 血浆 | aid in the management…monitor disease status in patients with confirmed pancreatic cancer who show measurable CA 19-9 values | 决策摘要 + summary |
| K191973（+K200997 special 510(k)） | 2019-10-22（2020-05-14） | Fujirebio / Lumipulse G CA19-9-N | CLEIA（ALP–AMPPD） | 血清、Na/Li-heparin、K2-EDTA 血浆 | aid in the management of patients diagnosed with cancer of the exocrine pancreas who have detectable levels of CA 19-9 at some point | 决策摘要 + summary |
| K100375 | 2011-04-06 | Siemens / Dimension Vista LOCI CA19-9 Flex + LOCI 7 Calibrator | 均相 LOCI 化学发光（1116-NS-19-9） | 血清、Li-heparin/EDTA 血浆 | serial measurement…aid in managing patients…confirmed pancreatic cancer who have levels…exceeding the median concentration determined for the apparently healthy cohort | 决策摘要 + summary |

### 3. 预期用途与声明类型
- 全部为**监测/管理已确诊外分泌胰腺癌患者**，Rx only，非 POC。Lumipulse 标签："should not be used for cancer screening or diagnosis"。
- 独特限定条件：仅适用于"病程中某时点 CA 19-9 曾高于健康人群 ULN（或中位数）/曾有可测值"的患者；并附 **Lewis 血型抗原阴性者无法产生 CA 19-9** 的警示（K253528、K191973）。
- K253528 明确说明抗体为 1116-NS-19-9（Fujirebio 授权），"assays using antibodies other than 1116-NS-19-9 may give different results"。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 健康参考上限 | 建立方式 / 人群 | serial change 定义 |
|---|---|---|---|
| K253528 Atellica CA 19-9 II | **35 U/mL**（ADVIA Centaur 既定 ULN） | 本次按 EP28-A3c 用 40 名健康者（20 男 20 女）验证，95% (38/40) 在范围内 | **>15% 升高** = positive change（沿用 predicate ADVIA Centaur 标签值） |
| K231525 VITROS CA 19-9 | **≤37 U/mL** | predicate K052889 建立：健康 200（100 男 100 女）中 194 ≤37、6 在 37.1–70、0 >70；本次 60 名健康者验证（57 ≤37、2、1 >70，>70 者经 predicate 复核一致）；标签分布：胰腺癌 50 中 33 >70；结直肠癌 100 中 43 >70；肝硬化 67 中 6 >70 | 文件未载明（Refer to K052889） |
| K191973 Lumipulse G CA19-9-N | 参考区间 **<0.7–50.0 U/mL**（2.5th–97.5th，n=240，22–93 岁；男 <0.7–40.1，女 <0.7–50.6）；88.8% ≤25、97.5% ≤50、100% ≤75 | 健康 240 + 良性胰腺炎 75 + 未治疗胰腺癌 120（21.7% ≤25）+ 胆道/乳腺/CRC/胆囊/肝/肺/卵巢/胃癌各 37–40 + 良性肺/肾/肝硬化/糖尿病/胆囊/肝炎/直肠息肉各 38–42 | **≥15% 升高**；ROC 表给出 10/20/30/40/50% 变化的 Sens/Spec（10%：67.1/55.6；20%：62.9/62.8；50%：44.3/76.6） |
| K100375 Vista LOCI CA19-9 | **≤37 U/mL**（300 名健康成人 98.7% ≤37） | 非参数分布；良性 200 例、恶性 398 例分布表（胰腺癌 105 中 30.5% ≤37、32.4% >1000） | **RCV = 84.7%**（CV_A 13.8% @13.9 U/mL，CV_I 27.2% 取文献）；predicate ADVIA Centaur 用 15%（取自其说明书） |

- 观察：ULN 35/37 U/mL 为各家健康人群百分位/分布（非 ROC）；15% 变化为 Siemens Centaur 传统值，被 Atellica 与 Lumipulse 沿用；Vista 用 RCV 公式得到 84.7%，代价是极低灵敏度（见第 8 节）。

### 5. 生物学 / 生理学依据
- 申报文件内：Lewis 抗原基因型阴性者不表达 CA 19-9，表型检测可能不足以识别（K253528、K191973）；良性胰腺炎、肝硬化、胆道疾病可升高（K191973 分布：胰腺炎 n=75 中位 18.5、最高 1966；胆道癌中位 73.3）；与 AFP/CA125/CA15-3/CA27.29/CEA/PSA 无交叉（K253528、K100375）。K253528 溯源注明"no internationally recognized reference method or material for CA 19-9"，按 ISO 17511:2020 赋值。
- 背景（非申报文件）：CA 19-9 为 sialyl-Lewis^a 糖抗原，胆道梗阻时非特异升高；约 5–10% 人群 Lewis 阴性。

### 6. 样本类型与样本要求
- K253528：血清、K2/K3-EDTA、Li-heparin 血浆 n=62 配对（4.0–674.7），加权 Deming 斜率 0.99–1.04，r ≥0.991。
- K231525：按 **CLSI EP35** 做矩阵等效，n=41（4 例超范围剔除），斜率 0.97/0.98，R² 0.99/1.00。
- K191973：50 名供者 5 种管型（红头、SST、K2-EDTA、Li-hep、Na-hep），斜率 0.989–1.009；样本稳定 2–10°C 4 天、−20°C 14 天。
- K100375：Li-hep n=60（y=1.013x−1.649）、EDTA n=63（y=0.9754x+0.4812）；加标回收在 Li-hep 中近 cutoff 处偏高（113.8–115.5%）。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果 |
|---|---|---|
| K253528 | EP05-A3；EP06 2nd；EP07 3rd；EP09c；EP17-A2；EP25 2nd；EP28-A3c；EP34；EP37；**ISO 17511:2020** | 精密度 5 水平 n=80，within-lab CV 3.5–8.2%；线性 1.2–732；hook 3,845,000 U/mL 无；自动稀释 1:10/100/200（EMI 至 140,000）；LoB 1.2/LoD 1.5/LoQ 2.0；HAMA 1000 µg/L、RF 1000 IU/mL、biotin 3500 ng/mL；交叉反应 6 种标志物；货架期 10 个月、在机 23 天 |
| K231525 | EP05-A3；EP06；EP07；EP17-A2；EP25-A；EP28-A3c；EP35；EP37 | 精密度 3 lot n=240 总 CV 3.5–6.8%；线性 0.9–1152（低端 12–14% 偏差）；**Hb 1000 mg/dL 偏差 190%（7.1 U/mL）/30.5%（45.2）、RF 1035 U/mL 偏差 27.4%** → 标签限制；货架期仅 20 周；LoD 1.4/LoQ 1.4 |
| K191973 | EP05-A3；EP6-A；EP07；EP14-A2；EP17-A2；EP25-A；EP28-A3c；EP34 | 精密度 n=80 总 CV 1.6–5.3%；线性 1.11–486；LoB 0.1/LoD 0.19/LoQ 0.57（15% CV）；稀释至 80,000；hook 200,000；干扰 22 种物质 |
| K100375 | EP05-A2；EP06-A；EP07-A2；EP09-A2；EP17-A | 精密度含站点/校准品 lot 变量；LoB 1.0/LoD 2.0；hook 1,230,509 U/mL；HAMA 至 327.1 mg/mL；62 种外源物 |

### 8. 临床验证设计与结果
- **K253528（Atellica）**：前瞻/回顾（"remnant sample documentation"）79 例外分泌胰腺癌（IA–IV 期），排除后 71 例（47–87 岁，均 65.3），413 份样本/342 次随访对；状态分类 NED/SD/RD/PD/RC。>15% 升高：Progression 63 中 43 阳、No-progression 279 中 111 阳 → **Sens 68.3% (59.5–78.4)、Spec 60.2% (54.7–65.2)、PPV 27.9%、NPV 89.4%**。与 predicate（ADVIA Centaur XPT）同一队列 n=321：Sens 66.7 vs 68.4%，Spec 60.9 vs 61.2%。比值中位数：NED 0.92、SD 1.06、RD 0.71、PD/RC 1.48。
- **K191973（Lumipulse）**：83 例（47–87 岁，中位 65；89.2% 白人），374 对（均 5.6 次/人）。≥15%：Progression 70 中 45 阳、No-progression 304 中 120 阳 → **Sens 64.29% (55.49–73.08)、Spec 60.53% (55.60–65.45)、Total concordance 61.23%、PPV 27.27%、NPV 88.04%**。同队列 predicate ARCHITECT CA 19-9XR（≥14%）n=301：Sens 64.71%、Spec 62.40%。方法比对 vs ARCHITECT n=84：Passing-Bablok 斜率 0.903 但 **R 仅 0.6628、平均偏差 164%**——文件如实记录，未做进一步解释。
- **K100375（Vista）**：38 回顾 + 34 前瞻 = 72 例（45.4–69.5 岁），189 对。>84.7%（RCV）升高：Progression 73 中 14 阳、No-progression 116 中 13 阳 → **Positive concordance 19.2% (10.9–30.1)、Negative concordance 88.8% (81.6–93.9)、Total 61.9%**；predicate（>15%）Positive 39.7%、Negative 69.8%、Total 58.2%。方法比对 vs ADVIA Centaur n=293，Passing-Bablok y=1.12x−5.71，R² 0.787。
- **K231525**：临床 "Refer to K052889"；方法比对 n=117/118（两 predicate lot），加权 Deming 斜率 0.96/0.97，R² 0.99。

### 9. 厂家间差异与要点
- 抗体均为 1116-NS-19-9（Fujirebio 授权），但不同平台间 R 可低至 0.66（Lumipulse vs ARCHITECT）、R² 0.787（Vista vs Centaur）——CA 19-9 是本组中方法间一致性最差的靶点。
- 变化阈值：15%（Siemens/Fujirebio）vs 84.7% RCV（Vista）→ Sens 64–68% vs 19%。
- AMI：Atellica 2.0–700；VITROS 1.4–1000；Lumipulse 0.7–500；Vista 2–1000 U/mL。
- Ortho VITROS 存在明显 Hb/RF 干扰（其它三家未报告类似问题）。
- 预期用途措辞差异："exceeding the upper limit of normal"（Atellica 2026）vs "exceeding the median concentration"（Vista 2011 / Centaur 2003）vs "detectable levels"（Lumipulse）。

### 10. 来源
- K253528: https://www.accessdata.fda.gov/cdrh_docs/reviews/K253528.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253528.pdf
- K231525: https://www.accessdata.fda.gov/cdrh_docs/reviews/K231525.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K231525.pdf
- K191973: https://www.accessdata.fda.gov/cdrh_docs/reviews/K191973.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K191973.pdf
- K200997: https://www.accessdata.fda.gov/cdrh_docs/reviews/K200997.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K200997.pdf
- K100375: https://www.accessdata.fda.gov/cdrh_docs/reviews/K100375.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf10/K100375.pdf

---

## 前列腺特异抗原用于前列腺癌管理/监测（Total PSA for Management of Prostate Cancer）

### 1. 法规定位
- Product code **LTJ**（Prostate Specific Antigen (PSA) for Management of Prostate Cancers）；21 CFR 866.6010；Class II；Panel Immunology (82)。
- openFDA 累计清关数：**total = 34**（多数为 1997–2003 年清关；2013 年后仅 NanoEnTek FREND 与 2026 年 Siemens Atellica tPSAII）。
- **与筛查/检测用途的界线（文件实证）**：K251630 的 predicate 为 **P950021/S015（PMA）**，其"aid in the detection of prostate cancer in conjunction with DRE in men aged 50 years and older"用途及全部分析性能均"Refer to P240021"（PMA）；510(k) K251630 仅承载**新增的 "management (monitoring)" 用途**及其临床监测研究。即：PSA 用于癌症检测/筛查 = Class III PMA（product code MTF），PSA 用于已确诊患者监测 = Class II 510(k)（LTJ）。同一试剂可同时具备两种用途，但分别走两条法规路径。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K251630 | 2026-01-05 | Siemens / Atellica IM total PSA II (tPSAII) | 吖啶酯 CLIA，三单抗（含未标记抗 fPSA 抗体，等摩尔设计） | 血清、EDTA/Li-heparin 血浆 | aid in the detection…（PMA 部分）+ **aid in the management (monitoring) of patients with prostate cancer**（本 510(k) 新增） | 决策摘要 + summary |
| K162378 | 2017-05-17 | NanoEnTek / FREND PSA Plus（修改：样本量 35 µL、LoQ 0.08） | 微流控卡匣荧光免疫（POC 型小型仪器） | 血清、Li-heparin、K3-EDTA 血浆 | serial measurement of total PSA…aid in the management of patients with prostate cancer | 决策摘要 + summary |
| K124056 | 2013-05-29 | NanoEnTek / FREND PSA Plus（原始清关） | 同上 | 血清、heparin/EDTA 血浆 | 同上；"for use in clinical laboratories upon prescription" | 决策摘要 + summary |

### 3. 预期用途与声明类型
- 声明类型：**监测（management/monitoring）已确诊前列腺癌患者**；Rx only。FREND 为台式小型荧光读卡仪（属 POC 形态），但标签限定 "for use in clinical laboratories upon prescription by the physician"（K124056），非 CLIA-waived。
- K251630 标签警示：不同方法 PSA 值不可互换，更换方法须补做序列检测确认基线。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 参考区间 | 建立方式 | serial change 定义 |
|---|---|---|---|
| K251630 Atellica tPSAII | Refer to P240021 | — | **RCV = 50.0%**（"increase in 50% RCV of the visit compared to previous visit"）；来源/推导文件未载明（仅说明用于候选与 predicate 两者） |
| K162378 FREND PSA Plus | 121 名 50–88 岁健康男性：均值 1.27、中位 1.01、**95th 2.92、97th 3.19、99th 4.59 ng/mL** | 表观健康男性百分位 | **PSA >1 ng/mL 时升高 >20%；PSA ≤1 ng/mL 时升高 >0.2 ng/mL**；"An absolute cutoff is not applicable for monitoring" |
| K124056 FREND PSA Plus | 196 名 ≥50 岁健康男性：**95th 1.68 ng/mL**（90% CI 1.28–2.03）；predicate（Tosoh ST AIA-PACK PA）1.52；">99% ≤4.0 ng/mL" | 百分位分布（0.5th–95th 列表） | **20%** = 2.5 × 总 CV（FREND 总 CV 8% → 20%；predicate Tosoh CV 3.4% → 8.5%），定义 v_ij=1 若 (x_j−x_i) ≥ d·x_i |

- 用户提到的 4 ng/mL：仅在 K124056 中以"健康男性 >99% ≤4.0 ng/mL"形式出现，且文件明确该值不是监测 cutoff。4.0 ng/mL 作为检测辅助阈值属 PMA 范畴（背景，非本组文件）。

### 5. 生物学 / 生理学依据
- 申报文件内：K124056 说明治疗（前列腺切除、放疗、激素）可使 PSA 下降而肿瘤仍在进展，故须结合 DRE/影像；K162378 干扰实验涵盖前列腺癌治疗药（flutamide、goserelin、leuprolide、finasteride、tamsulosin、docetaxel）及 PAP、kallikrein（K124056）。
- 背景（非申报文件）：PSA（hK3）为前列腺上皮丝氨酸蛋白酶；术后生化复发定义（如 ≥0.2 ng/mL 两次）与此处"≤1 ng/mL 时 +0.2 ng/mL"规则思路相近但非文件引用。

### 6. 样本类型与样本要求
- K251630：血清、EDTA/Li-heparin 血浆；矩阵比对 "Refer to P240021"；临床研究要求血清且在 SOC PSA 采血 ±28 天内。
- K162378：Li-hep/K3-EDTA vs 血清 n=40（0.08–23.59），Passing-Bablok 斜率 0.96 (0.89–1.03)/1.03 (0.99–1.09)，r 0.99。
- K124056：样本 2–8°C 一周，>1 周须 ≤−20°C（临床纳入标准）。
- 溯源：FREND 溯源 WHO 96/670（90:10 PSA-ACT:fPSA），校准品经 ARCHITECT total PSA (P910007)/Tosoh AIA PA (P910065) 确认；Atellica 溯源 "Refer to P240021"。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果 |
|---|---|---|
| K251630 | 仅列 **CLSI EP34 1st Ed**（其余 "Refer to P240021"） | 自动稀释 1:5/10/50/100/500（Multi-Diluent 2）；AMI 0.009–50.00 ng/mL |
| K162378 | EP05-A3；EP6-A；EP07-A2；EP09-A2；EP14-A3；EP17-A2；EP25-A；C28-A3 | 精密度 3 lot n=240：0.08 ng/mL 总 CV 14.6%、4.0 时 7.4%、21.4 时 6.1%；站点间 n=75 总 CV 7.9–11.0%；线性 0.10–25.53 与 0.05–1.00；LoB 0.02/LoD 0.03/LoQ 0.08；hook 1200 ng/mL 无；干扰 Hb 500 mg/dL、RF 1075 IU/mL、HAMA 70 ng/mL；vs ARCHITECT tPSA n=207 斜率 0.98 |
| K124056 | EP05-A2；EP06-A；EP07-A2；I/LA19-A；EP17-A | 精密度总 CV ~8%（@4 ng/mL）；LoB 0.04；方法比对 vs Tosoh ST AIA-PACK PA（n=236 监测样本 + 160 单点）；HAMA/PAP/kallikrein |

### 8. 临床验证设计与结果
- **K251630（Atellica tPSAII）**：前瞻性、6 个采集点，88 名 >50 岁确诊前列腺癌男性（54–86 岁，均 67.7；83% 白人、17% 非裔），随访 121–1014 天（中位 768），访视间隔中位 182 天；323 份（88 基线 + 235 随访）。状态 NED/Stable/Responding/Progression 由医生依据 DRE/活检/超声/MRI/CT/PET/骨扫描判定。RCV 50%：Progression 33 中 18 阳、No-progression 290 中 36 阳 → **Sens 54.5% (36.4–71.9)、Spec 87.6% (83.3–91.2)、PPV 33.3%、NPV 94.4%、LR+ 4.39、LR− 0.52**，患病率 10.22%。predicate（Atellica IM PSA，PMA）同队列：Sens 54.5%、Spec 89.3%。候选 vs predicate：**PPA 91.8% (45/49)、NPA 96.7% (265/274)**。
- **K162378（FREND 修改版）**：63 例（49 白人；I 期 3、II 期 25、III 期 22、IV 期 12；40 例接受手术），257 份（63 基线 + 194 随访，均 4.08 次/人），回顾性，样本须来自确诊后。规则（>20% 或 >0.2 ng/mL）：Progression 73 中 46 阳、No-progression 121 中 23 阳 → **Sens 63.0% (51.5–73.2)、Spec 81.0% (73.1–87.0)、PPV 66.7%、NPV 78.4%、LR+ 3.32**；分期/手术分层亦列出（手术组 Sens 60.5%/Spec 82.2%；非手术 66.7%/79.2%）。predicate（原 FREND）同队列 Sens 65.8%、Spec 81.0%、Total concordance 75.3%。方法比对 vs 原 FREND n=64 斜率 0.98。
- **K124056（FREND 原始）**：75 套回顾性 serial sets（311 份，3–9 次/人，均 4.2），236 对；Gleason 5–9，多种治疗。≥20%：Progression 108 中 84 阳、No-progression 128 中 46 阳 → **Sens 77.8% (69.1–84.6)、Spec 64.1% (55.5–71.9)、PPV 64.6%、NPV 77.4%、ROC AUC 0.759 (0.697–0.822)**。predicate Tosoh 用 8.5%。另有 85 例单点前列腺癌样本分布。

### 9. 厂家间差异与要点
- 显著变化定义三种：RCV 50%（Siemens，推导未载明）、20%/0.2 ng/mL 双层规则（FREND 2017）、2.5×CV=20%（FREND 2013）。对应 Sens 54.5% / 63.0% / 77.8%，Spec 87.6% / 81.0% / 64.1%。
- FREND 是 POC 型平台，低端精密度差（0.08 ng/mL 时 CV 14.6%），AMI 上限仅 25 ng/mL。
- Siemens 采用 PMA + 510(k) 双路径，510(k) 文件几乎不含分析性能。

### 10. 来源
- K251630: https://www.accessdata.fda.gov/cdrh_docs/reviews/K251630.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K251630.pdf
- K162378: https://www.accessdata.fda.gov/cdrh_docs/reviews/K162378.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162378.pdf
- K124056: https://www.accessdata.fda.gov/cdrh_docs/reviews/K124056.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K124056.pdf

---

## 甲胎蛋白用于睾丸癌/生殖细胞肿瘤（AFP for Testicular / Germ Cell Cancer）

### 1. 法规定位
- Product code **LOJ**（Kit, Test, Alpha-fetoprotein for testicular cancer）；21 CFR 866.6010；Class II；Panel Immunology (82)。
- openFDA 累计清关数：**total = 17**。
- **与神经管缺陷 AFP 的 Class III 区别（文件实证）**：K071597 的 predicate 为 **Abbott AxSYM AFP, P820060/S019（PMA）**，其差异表列出 predicate 另有"quantitative determination of AFP in amniotic fluid at 15–21 weeks gestation to aid in the detection of fetal open neural tube defects (NTD)"用途，而候选 Vista AFP 该项为 "Not applicable"。即 AFP 用于羊水/母血 NTD 检测属 Class III PMA（文件仅记录该用途归 PMA P820060，未载明其 product code），AFP 用于非精原细胞睾丸癌管理属 Class II 510(k)（LOJ）。K071597 差异表还显示 predicate 的样本类型含 "Serum, plasma and amniotic fluid"，候选仅血清。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K220176 | 2022-09-15 | Roche / Elecsys AFP（cobas e 601，biotin 改良） | ECLIA 一步夹心 | 血清、Li-hep、K2/K3-EDTA 血浆 | aid in the management of patients with non-seminomatous germ cell tumors | 决策摘要 + summary |
| K213626 | 2022-06-15 | Ortho / VITROS AFP Reagent Pack（biotin 改良 + 单位换算因子 1.04→1.21） | 免疫化学发光（羊抗 AFP 捕获 + HRP 鼠单抗） | 血清 | aid in the management of patients with non-seminomatous testicular cancer | 决策摘要 + summary |
| K071597 | 2008-06-04 | Dade Behring (Siemens) / Dimension Vista AFP Flex + LOCI 5 Calibrator | 均相 LOCI | 血清 | aid in managing non-seminomatous testicular cancer when used in conjunction with physical examination, histology/pathology and other clinical evaluation procedures | 决策摘要 |
| K090236 | 2009-04-27 | Siemens / Dimension Vista AFP（新增 Li-heparin 血浆） | 同上 | 血清、Li-heparin 血浆 | 同上 | 决策摘要 |

### 3. 预期用途与声明类型
- 均为**管理/监测非精原细胞性生殖细胞肿瘤**，Rx only，非 POC。Vista 版本强调须与体检、组织病理及其它临床评估联用。
- 三份近年/补充文件（K220176、K213626、K090236）临床均 "Refer to" 原清关（K981282、K983031、K071597）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 参考上限 | 建立方式 | serial change |
|---|---|---|---|
| K071597 Vista AFP | 健康男性（18–61 岁）n=231：**97.4% <8.0 ng/mL [6.6 IU/mL]** | 分布表（共 803 份：睾丸精原细胞瘤 10、肝细胞癌 123（44.7% ≤8）、胰腺 35、其他 GI 164、肝硬化 122（50.8% ≤8）、肝炎 118） | **RCV = 33.7%**（predicate AxSYM 37.7%）；CV_I 取文献 12%（Trapé 2003）；Clinical cut-off "Not applicable" |
| K213626 VITROS AFP | **<7.22 IU/mL（<8.74 ng/mL，新因子 1.21）**，408 名正常献血者上 97.5% | 本次 EP28-A3c 用 60 名健康非吸烟者（30 男 30 女）3 lot 验证：1 男 8.19–8.77、1 男 6.98–7.26，其余 <7.22；标签分布：原发性肝癌 36（24 >500）、睾丸非精原 117（41 >500、25 ≤7.22）、甲肝 50、健康男 210/女 198 | 文件未载明；标签要求 "Serial testing…reported and interpreted with the same conversion factor and cutoff values"（因换算因子变更使 ng/mL 区间上升 16%） |
| K220176 Elecsys AFP | 140 名正常者 **97% ≤6.90 IU/mL (≤8.7 ng/mL)**（K981282 建立） | Refer to K981282 | 文件未载明 |
| K090236 | Refer to K071597 | — | — |

### 5. 生物学 / 生理学依据
- 申报文件内：肝细胞癌与肝硬化 AFP 显著升高（Vista 分布表）；精原细胞瘤 100% ≤8 ng/mL（n=10）而非精原细胞瘤 35% >500 IU/mL（Ortho 表）；交叉反应物 hCG（1,000,000 mIU/mL）、α1-酸性糖蛋白、α1-抗胰蛋白酶、铜蓝蛋白、HPL、转铁蛋白、催乳素无交叉（K213626）；溯源 WHO 1st IRP 72/225（Roche、Ortho、Vista 均是），IU/mL↔ng/mL 换算因子各家不一（Ortho 由 1.04 改为 1.21 "to align…to other vendors"）。
- 背景（非申报文件）：AFP 由卵黄囊成分分泌，精原细胞瘤不产生 AFP；IGCCCG 分期用 AFP 绝对值（<1000 / 1000–10000 / >10000 ng/mL）。

### 6. 样本类型与样本要求
- K220176：血清、Li-hep、K2-EDTA、K3-EDTA；矩阵 n=52/48/53（1.63–927 IU/mL），Passing-Bablok 斜率 0.976–0.985，r ≥0.998。
- K213626：仅血清；矩阵 "Not applicable"。
- K090236：新增 Li-heparin 血浆，n=70 配对（0.9–999.5 ng/mL），Passing-Bablok 斜率 0.99 (0.96–1.00)、截距 −0.05、r 0.997；样本 −20°C 冻存 ≤60 天后检测。
- K071597：样本 7 天应激/−20°C 30 天/冻融稳定。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果 |
|---|---|---|
| K220176 | EP05-A3；EP06-A；EP17-A2；EP28-A3c | 精密度 n=84 within-lab 3.4–4.2%；3 lot 总 CV 3.5–9.1%；线性 R² 0.958（15 水平），AMI 1.5–1000 IU/mL；hook 1,000,000 IU/mL 无；biotin 无干扰至 1200 ng/mL（原 60）；Hb 2200 mg/dL、RF 1500；LoB 0.75/LoD 1.5/LoQ 1.5；货架期 21 个月；方法比对 3 lot n≈181，斜率 0.965–0.971 |
| K213626 | EP05-A3；EP09c；EP07-A3；EP37；EP17-A2；EP06-A；EP28-A3c；EP25-A | 精密度 3 lot 总 CV 2.4–5.6%；线性 0.33–723 IU/mL；自动稀释 1:400、手工 1:4000；干扰 HAMA 800 µg/L、Hb 500、RF 900；LoB 0.229/LoD 0.476/LoQ 0.800；货架期 32 周；方法比对 n=150 加权 Deming 斜率 0.99–1.01 |
| K071597 | FDA bundling/leftover 指南；EP9-A2；EP17-A（另引 EP5-A2、EP7-A2） | 精密度 within-lab ≤2.3%；线性 1.7–1014 ng/mL；LoB/LoD 支持 0.5 ng/mL；hook 2877 IU/mL；13 种化疗药；HAMA 3 份 <−4.8%；方法比对 vs AxSYM n=317 斜率 0.93、vs ADVIA Centaur n=84 斜率 0.97 |
| K090236 | EP5-A2 | 血浆精密度（249 ng/mL：within-lab 2.0%；7.5 ng/mL：1.7%） |

### 8. 临床验证设计与结果
- **K071597（Vista AFP，唯一含监测临床研究）**：74 套回顾性睾丸癌 serial sets，排除 4 例非非精原细胞瘤，70 例（1.1–53.7 岁，均 30.8；无非裔样本；I–IV 期）；医生依据体检/影像（CT/MRI/X 线/超声）/手术（穿刺、睾丸切除）判定状态。逐访视 n=244。RCV >33.7% 升高：Progression 59 中 16 阳、No-progression 185 中 24 阳 → **Sensitivity 27.1% (16.4–40.3)、Specificity 87.0% (81.3–91.5)、Accuracy 72.5%**；predicate（RCV 37.7%）Sens 27.1%、Spec 90.3%。Sens+Spec >1 bootstrap 95% CI 1.0062–1.2861（下限勉强 >1）。与 predicate 一致性 overall 97.5%、PPA 100% (34/34)、NPA 97.1%。注：三分类表显示 Progression 59 例中 21 例 AFP 反而下降 >33.7%（响应治疗中仍被判进展），解释了低灵敏度。
- K220176 / K213626 / K090236：临床 Refer to 原清关；仅方法比对与参考区间验证（见第 7 节）。

### 9. 厂家间差异与要点
- 单位：Roche 与 Ortho 报 IU/mL（AMI 1.5–1000 与 0.8–500），Vista 报 ng/mL（0.5–1000）；Ortho 换算因子改变导致 ng/mL 期望值整体 +16%，是标签层面需要关注的"cutoff 漂移"案例。
- 参考上限：8.0 ng/mL（97.4%，Vista）、7.22 IU/mL（97.5th，Ortho）、6.90 IU/mL（97%，Roche）。
- 仅 Vista 2008 给出 RCV（33.7%）；在睾丸癌监测中 RCV 型规则灵敏度极低（27%）。
- biotin 改良是 2022 年 Roche/Ortho 申报共同主题。

### 10. 来源
- K220176: https://www.accessdata.fda.gov/cdrh_docs/reviews/K220176.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220176.pdf
- K213626: https://www.accessdata.fda.gov/cdrh_docs/reviews/K213626.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K213626.pdf
- K071597: https://www.accessdata.fda.gov/cdrh_docs/reviews/K071597.pdf
- K090236: https://www.accessdata.fda.gov/cdrh_docs/reviews/K090236.pdf

---

## 人附睾蛋白 4（HE4）及 ROMA 算法

### 1. 法规定位
- HE4：Product code **OIU**（Test, Epithelial Ovarian Tumor Associated Antigen (HE4)）；21 CFR 866.6010；Class II；Panel Immunology (82)。openFDA 累计清关数：**total = 5**（含 1 件 Tumor Marker Control K103676）。
- ROMA：Product code **ONX**（Ovarian adnexal mass assessment score test system）；21 CFR **866.6050**；Class II（special controls，"Class II Special Controls Guidance Document: Ovarian Adnexal Mass Assessment Score Test System"）；Panel Immunology (82)。openFDA 累计清关数：**total = 6**（含 OVA1 DEN090004 de novo 与 OVA1 Next Generation K150588；ROMA 4 件）。原型为 **K103358 ROMA (HE4 EIA + ARCHITECT CA 125 II)**，其 predicate 为 OVA1（k081754）。

### 2. 代表性 510(k) 一览

**HE4（OIU）**

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K072939 | 2008-06-09 | Fujirebio / HE4 EIA（首个 HE4，predicate 为 ARCHITECT CA 125 II k042731） | 手工 EIA（2H5 捕获/3D8-HRP，TMB） | 血清 | aid in monitoring recurrence or progressive disease in patients with epithelial ovarian cancer; serial testing | 决策摘要 |
| K093957 | 2010-03-18 | Fujirebio / ARCHITECT HE4 + Calibrators + Controls | CMIA（Abbott ARCHITECT i） | 血清 | 同上 | 决策摘要 |
| K112624 | 2012-09-10 | Roche / Elecsys HE4 + CalSet + PreciControl + CalCheck 5 | ECLIA（12A2 biotin / 2H5-Ru） | 血清、K2/K3-EDTA、Li-hep 血浆 | 同上 | 决策摘要 + summary |
| K151378 | 2015-11-24 | Fujirebio / Lumipulse G HE4 + Calibrators | CLEIA（12A2/3D8） | 血清、Li-hep、K2-EDTA 血浆 | 同上 | 510(k) summary（**无决策摘要**） |

**ROMA（ONX）**

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K103358 | 2011-09-01 | Fujirebio / ROMA (HE4 EIA + ARCHITECT CA 125 II) | 软件算法 + 2 项免疫分析；绝经前/后两套 logistic 方程，输出 0.0–10.0 分 | 血清 | aid in assessing whether a premenopausal or postmenopausal woman who presents with an ovarian adnexal mass is at high or low likelihood of finding malignancy on surgery; >18 岁、计划手术、尚未转诊肿瘤科；not a screening or stand-alone diagnostic assay | 决策摘要 + summary |
| K151502 | 2016-04-28 | Fujirebio / ARCHITECT ROMA | ARCHITECT HE4 + ARCHITECT CA 125 II | 血清 | 同上 | 决策摘要 + summary |
| K160090 | 2016-05-16 | Fujirebio / Lumipulse G ROMA | Lumipulse G HE4 + CA125II | 血清、K2-EDTA、Li-hep 血浆 | 同上 | 决策摘要 + summary |
| K153607 | 2016-06-15 | Roche / ROMA Calculation Tool Using Elecsys Assays (RCTUEA) | Elecsys HE4 + Elecsys CA 125 II | 血清、K2/K3-EDTA、Li-hep 血浆 | 同上 | 决策摘要 + summary |

### 3. 预期用途与声明类型
- **HE4**：监测已确诊上皮性卵巢癌复发/进展（serial testing），Rx only，非 POC；决策摘要多次强调 "There is no assay cut-off for monitoring the progression"。
- **ROMA**：**术前风险分层（辅助评估附件包块恶性可能性）**，是本组中唯一的"单次检测、定性 high/low likelihood"声明；带 PRECAUTION 黑框："should not be used without an independent clinical/radiological evaluation…not intended to be a screening test or to determine whether a patient should proceed to surgery. Incorrect use…carries the risk of unnecessary testing, surgery, and/or delayed diagnosis"。Rx only。绝经状态须基于卵巢功能的临床判断；两种绝经状态的分数同时报告给医生。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**HE4 参考上限（健康女性百分位）**

| 产品 | ULN | 建立方式 |
|---|---|---|
| K072939 HE4 EIA | **150 pM**（推荐单一上限） | 204 名健康女性（76 绝经前、103 绝经后、3 未分类、22 妊娠），非参数单侧 95th 百分位：全体 153.4、绝经后 154.1、绝经前 119.1、绝经前+妊娠 136.7 pM；因绝经状态难确定、95% CI 102–205，取 150 |
| K093957 ARCHITECT HE4 | **绝经前 70 pmol/L、绝经后 140 pmol/L** | 400 名健康女性（210/190，14–93 岁，~95% 白人）：95th/97.5th 百分位 全体 98.2/145.2；绝经前 65.3/90.6；绝经后 125.9/163.0；厂家取整为 70/140；分布表：绝经前 95.7% ≤70，绝经后 81.6% ≤70、95.3% ≤140 |
| K112624 Elecsys HE4 | 未设固定 cutoff；按年龄分层 95th 百分位（<40 岁 61.24；40–49 63.55；50–59 105.20；60–69 96.46；≥70 96.89 pmol/L）；绝经前 95th 67.36、绝经后 96.89 | 345 名健康女性（198/147，20–79 岁，289 白人、52 非裔）；99.7% ≤140；文件以 140 pmol/L 作为 NED 分析阈值 |
| K151378 Lumipulse G HE4 | 参考区间 2.5th–97.5th：全体 32.5–108.4；绝经前 31.9–87.1；绝经后 32.5–112.2 pM；96% <90、99% <135 | 240 名健康女性（120/120） |

**HE4 serial change**

| 产品 | 定义 | 推导 |
|---|---|---|
| K072939 | ≥25% | 2.5 × 最大总 CV 9.3% ≈ 25%；同时提供 ROC 多 cutoff 表（−29.0%…+67.8%）供临床自选；无公认 cutoff |
| K093957 | ≥14% | 假设总 CV ≤10% → 两次差值变异 14%；亦给出 0/5/10/14/20/25/50/75/100% 的 Sens/Spec 表 |
| K112624 | ≥20% | 2 × 总 CV 上限 7%；同样给出多 cutoff 表（0–100%） |
| K151378 | ≥18% | "takes into account the variability of the assay"（具体推导未载明） |

**ROMA cutoff（分数 0–10）**

| 产品 | 绝经前 | 绝经后 | 建立方式 |
|---|---|---|---|
| K103358（原型）、K151502、K160090 | **≥1.31** 高可能 | **≥2.77** 高可能 | 训练集（两项 Pilot Study 合并）logistic 模型；cut-point 按 **specificity 75%** 设定，并要求合并人群 sensitivity ≥80% 作为接受阈值 |
| K153607 Roche RCTUEA | **≥1.14** | **≥2.99** | 沿用同一 ROMA 算法框架，但因 Elecsys HE4/CA125 II 标准化差异而重设；推导细节文件未载明 |

- 注：用户提及的 "7.4%/25.3%" 或 "11.4%/29.9%" 形式是以"恶性概率 %"表示的 ROMA 阈值（背景：见于美国以外标签/文献）；**FDA 清关文件一律采用 0.0–10.0 分制，未出现上述百分比数值**，两者为同一 logistic 概率的不同表达尺度（文件未载明换算关系）。
- 灰区：ROMA 无灰区；决策摘要另以 cutoff 上下的中位数再分组显示恶性率随分数递增（如 K103358 绝经后 0–1.50 分 4.5%、1.50–2.77 10.6%、2.77–6.16 32.6%、6.16–10 91.1%）。
- 干扰导致的判定限制：K103358 中 **RF ≥500 IU/mL 使 ROMA 分数偏差 12.6–28.2%**，标签规定 RF >250 IU/mL 样本不适用 ROMA。

### 5. 生物学 / 生理学依据
- 申报文件内：HE4（WFDC2 基因产物）在肺腺癌上调，肺癌灵敏度 42%（K072939）；CHF 与非妇科良性病中 HE4 意外升高（K072939 CHF 特异度 78.1%，K112624 CHF 46.7% >140 pmol/L；文件称"not associated with the known tissue expression of the HE4 gene"）；HE4 随年龄上升（K112624 分层表）；良性妇科病（子宫内膜异位、肌瘤、囊腺瘤）中 HE4 多低于 ULN（K093957 良性妇科 96% <140），这是 ROMA 相对 CA125 的设计依据（文件隐含）。Roche 交叉反应物 SLPI、Elafin/SKALP 无交叉。
- 背景（非申报文件）：HE4 为 WAP 四二硫核心域蛋白 2；肾功能不全时升高（文件仅提及 ESRD 个例）。

### 6. 样本类型与样本要求
- K072939 / K093957 / K103358 / K151502：仅血清（K093957 验证红头、SST、Hemogard 管型等效；4 天 2–8°C、24 h 室温、6 次冻融）。
- K112624：血清、Li-hep、K2/K3-EDTA，n=40 配对（17.1–1458 pmol/L），Passing-Bablok 斜率 0.980–1.010；稳定 15–25°C 5 h、2–8°C 2 天、−20°C 12 周、2 次冻融。
- K151378 / K160090：血清、Li-hep、K2-EDTA（SST/K2EDTA/Li-hep 斜率 95% CI 在 0.9–1.1 内）；ROMA 矩阵比对 n=86 配对 Deming 斜率 1.00。
- K153607：血清、K2/K3-EDTA、Li-hep；ROMA 血清 vs K2-EDTA n=89，Deming 斜率 1.00；另做"最坏情况"模拟矩阵效应。
- 溯源：HE4 无国际标准；Fujirebio 用重组 Ig-HE4 融合蛋白（CHO 表达，MALDI-TOF 定值）；Roche 用 OvCar-3 细胞来源 HE4、标准化至 Fujirebio HE4 EIA。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果 |
|---|---|---|
| K072939 | EP5-A2；EP6-A（另用 EP17-A、EP7-A） | 三站点精密度总 CV 3.3–6.5%（验收 ≤20%，声明 <15%）；线性 15–900 pM，回收 84.8–102.1%；hook 300,000 pM 无；LoD 1.1–2.2 pM、LoQ 3.85 pM；脂类 3 g/dL 回收 90.3%、蛋白 12 g/dL 113.1%、HAMA 400 ng/mL、RF 568 IU/mL、10 种化疗药 |
| K093957 | C28-A2；EP7-A2；EP5-A2；EP9-A2；EN 13640（稳定性）；EP6-A；EP17-A；FDA Tumor Associated Antigen 510(k) 指南；leftover specimen 指南 | 精密度 2 lot 三站点总 CV ≤5.4%（上 95% CL ≤6.1%）；线性 20–1500（FDA 重新分析偏差 −2.9%–7.3%）；1:10 稀释回收 91–106%；LoB 0.05/LoD 0.18/LoQ 1.0 pmol/L（FDA 按 30% TE 重算）；HAMA 45–155 ng/mL、RF 21–445 IU/mL；交叉 CA125/CA15-3/CA19-9/CEA/AFP 100×；方法比对 vs EIA n=390 Passing-Bablok 斜率 1.06 (1.03–1.08)、Spearman 0.97；单次复孔可行性 |
| K112624 | EP5-A2；EP6-A；EP17-A | 精密度 n=84 总 CV 2.7–4.3%；lot-to-lot n=126 斜率 0.958；线性 12.6–1510；hook 40,000 pmol/L；LoB 5/LoD 15/LoQ 20 pmol/L（标签声明，实测 0.23–4.4）；biotin 50 ng/mL、HAMA 805 ng/mL、RF 1500、IgG 70 g/L；18 常规药 + 14 抗癌药 |
| K151378 | ISO 17511:2003；EP5-A3；EP7-A2；C28-A3c；EP17-A2；EP6-A；EP9-A3；FDA TAA 指南 | 精密度总 CV ≤3.5%（三站点 ≤6.1%）；线性 20–1500；回收 91–107%；hook 300,000（>30,000 警示）；LoB 0.1/LoD 3.5/LoQ 3.5 pM；biotin 19.8 mg/dL、HAMA 1000 ng/mL、RF 1000；23 种药物；vs HE4 EIA n=143 加权 Deming 斜率 1.035 |
| K103358 | EP7-A；C28-A3（另用 EP5-A2） | ROMA 分数精密度：lot 间总 CV 绝经前 ≤7.72%、绝经后 ≤4.10%；三站点总 CV 绝经前 2.07–25.9%、绝经后 0.98–11.16%；模拟精密度覆盖 0–10；干扰（RF 见第 4 节） |
| K153607 | EP05-A3；EP09-A3；C28-A3c；ONX Special Controls 指南 | ROMA 精密度 n=84 总 CV ≤6.5%（绝经前）/≤4.7%（绝经后）；lot 间最高 12.2%（0.71 分处）；vs K103358 ROMA n=187 Deming 斜率 0.99–1.00 |
| K151502 / K160090 | EP5-A2/A3；EP7-A2；EP09-A3；C28-A3(c)/EP28-A3c；ONX Special Controls 指南 | vs K103358 ROMA：ARCHITECT 绝经前 Deming 0.98、绝经后 1.00；Lumipulse 绝经前 1.00、绝经后 1.00 |

### 8. 临床验证设计与结果

**HE4 监测研究（Progression vs No-progression，逐访视）**

| K 号 | 人群 | n（对） | 阈值 | Sens | Spec | 其他 |
|---|---|---|---|---|---|---|
| K072939 HE4 EIA | 80 名 EOC 女性（85% III/IV 期；66 绝经后），回顾性，来自单一癌症中心；非劣效于 ARCHITECT CA125 II | 354 | ≥25% | 60.3% (76/126) | 75.0% (171/228) | TPR−FPR 0.353 (0.251–0.455)；**ROC AUC HE4 0.725 (0.675–0.770) vs CA125 0.709 (0.659–0.756)**；CA125 ≥25%：Sens 69.1%、Spec 68.9%；按月变化率调整 AUC 0.73；HE4 ≤150 pM 判 NED 的 Sens 89%、Spec 64% |
| K093957 ARCHITECT HE4 | 76 名 EOC（均 54 岁，84% 白人，仅 2 绝经前；III 期 57%），506 份（3–31 次/人，随访中位 368 天，间隔中位 83 天） | 430 | ≥14% | **53.5% (43.2–63.6)** | **78.5% (73.7–82.9)** | PPV 42.7%、NPV 85.0%、总一致 72.8% (68.3–76.9)；成功标准：GEE 模型 log OR 检验总一致 >62.8%（EIA 既定值），通过；AUC 0.685 (0.618–0.750)；HE4 ≤140 pmol/L 对 NED Sens 97.8%、Spec 46.2%；NED 组 %change 中位 −1.4%，Progression +16.8%，Responding −15.6% |
| K112624 Elecsys HE4 | 80 名 EOC（20–85 岁，71 绝经后；51 有分期，43 例 III/IV），493 份 | 413 | ≥20% | **46.9% (37.3–56.3)** | **84.0% (81.7–86.3)** | AUC 0.699 (0.630–0.767)，与 ARCHITECT 统计等价；绝经后 Sens 44.4%/Spec 85%，绝经前 66.7%/76.3%；III 期 47.9%/83.0%；HE4 ≤140 对 NED Sens 99.4% (175/176)、Spec 34.2%；ROC 各 cutoff 表与 ARCHITECT 逐点对照（如 20%：46.9 vs 47.5；25%：39.5 vs 39.4） |
| K151378 Lumipulse G HE4 | 72 例，330 对（均 5.6 次/人） | 330 | ≥18% | 49.2% (27.8–70.6) | 79.6% (65.9–93.2) | 总一致 73.9%、PPV 35.3%、NPV 87.4% |

**ROMA 关键（pivotal）研究——四份文件共用同一前瞻性、多中心（13 站点）、盲法队列（512 例入组）**

- 设计：>18 岁、因卵巢囊肿/附件包块（单纯/复杂/实性）在综合或专科医院就诊于 generalist、已决定手术者；术前由非妇科肿瘤医师完成 Initial Cancer Risk Assessment（**ICRA**，76.1% 由普通妇产科医生完成）；所有患者手术，本院病理 + 独立病理学家（University of Maryland）复核；绝经状态：病历 → 年龄（≤49 绝经前、≥55 绝经后）→ 停经 1 年 → FSH。金标准：组织病理。评价 ICRA 单独、ROMA 单独、ICRA+ROMA 联合（任一阳性即阳性）。主要终点：联合使用较 ICRA 单独提高 NPV（差值 95% CI 不含 0）。
- 可评估例数：K103358 461（240 绝经前/221 绝经后；良性 81.3%、LMP 3.9%、EOC 10.4%）；K151502 459；K160090 450；K153607 455。

| 分析集：EOC+LMP，合并绝经状态 | K103358（HE4 EIA + ARCHITECT CA125） | K151502（ARCHITECT ROMA） | K160090（Lumipulse G ROMA） | K153607（Roche RCTUEA） |
|---|---|---|---|---|
| 患病率 | 15.0% (66/441) | 15.0% (66/440) | 15.1% (65/431) | 14.9% (65/436) |
| ICRA 单独 Sens/Spec | 77.3% / 84.3% | 77.3% / 84.2% | 76.9% / 84.2% | 76.9% / 84.4% |
| ROMA 单独 Sens/Spec | 87.9% / 75.5% | 87.9% / 86.1% | 87.7% / 76.0% | 86.2% / 79.5% |
| ICRA+ROMA Sens/Spec | 90.9% / 67.2% | 90.9% / 75.7% | 90.8% / 67.5% | 90.8% / 70.4% |
| ICRA+ROMA NPV（vs ICRA） | 97.7% vs 95.5%（Δ2.2%，CI 0.43–3.98） | 97.9% vs 95.5% | 97.6% vs 95.4% | 97.8% vs 95.4%（Δ2.3%，显著） |
| LR+（ICRA+ROMA 双阳） | 9.94 vs ICRA 4.91 | 13.88 vs 4.90 | 10.01 vs 4.85 | 11.18 vs 4.92 |

- 绝经分层（K103358，EOC+LMP）：绝经前患病率 6.8%，ICRA Sens 43.8% → ICRA+ROMA 81.3%，Spec 90.0% → 68.6%，NPV 95.7% → 98.1%（Δ2.4%，CI 0.07–4.73）；绝经后患病率 24.4%，Sens 88.0% → 94.0%，Spec 76.1% → 65.2%，NPV 95.2% → 97.1%（CI −0.75–4.66，不显著）。K153607 对应值：绝经前 Sens 43.8%→81.3%、NPV 95.8%→98.2%；绝经后 87.8%→93.9%、NPV 94.7%→97.1%（均报告显著）。
- "All cancers + LMP"分析集（K103358，n=461，患病率 18.7%）：ICRA+ROMA Sens 88.4%、Spec 67.2%、NPV 96.2%；联合使用较 ICRA 额外检出 13 例癌（含 LMP）。
- 分数与恶性率关系（K153607，绝经后）：0–1.50 分 4.6%、1.50–2.99 12.3%、2.99–6.57 44.7%、6.57–10 92.1%。
- 健康女性 ROMA 分布（各 n=240 左右）：K103358 绝经前 95th 2.36、绝经后 2.75（15.8%/5.0% 被判高可能）；K153607 1.70/2.58；K151502 1.65/2.02；K160090 1.73/2.03。良性妇科病 23–26% 判高可能（K103358）。

### 9. 厂家间差异与要点
- HE4 ULN：150 pM（EIA 单一）→ 70/140 pmol/L 绝经分层（ARCHITECT、Roche 分析亦用 140）→ Lumipulse 以 90/135 pM 为分布标记；serial change 14%/18%/20%/25%，均由分析 CV 倍数推导，无临床共识值，各文件均附多 cutoff ROC 表让临床自选。
- HE4 监测灵敏度普遍 47–60%、特异度 75–84%，与 CA125 相当（非劣效设计，K072939）；HE4 绝对值 ≤140 对 NED 灵敏度 98–99%，是文件强调的第二种判读方式。
- ROMA cutoff：Fujirebio/Abbott 系 1.31/2.77；Roche 1.14/2.99（因 HE4/CA125 定值差异）；四份 ROMA 文件复用同一 512 例队列（各自重测样本），联合 NPV 提升 ~2.2–2.4% 是共同主要终点。
- ROMA 特有限制：RF >250 IU/mL 不适用（K103358）；绝经状态由医生判定；两种分数同时报告。
- 文件类型缺口：K151378 仅有 510(k) summary，无决策摘要。

### 10. 来源
- K072939: https://www.accessdata.fda.gov/cdrh_docs/reviews/K072939.pdf
- K093957: https://www.accessdata.fda.gov/cdrh_docs/reviews/K093957.pdf
- K112624: https://www.accessdata.fda.gov/cdrh_docs/reviews/K112624.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K112624.pdf
- K151378: https://www.accessdata.fda.gov/cdrh_docs/pdf15/K151378.pdf（决策摘要 reviews/K151378.pdf 不可获得）
- K103358: https://www.accessdata.fda.gov/cdrh_docs/reviews/K103358.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf10/K103358.pdf
- K151502: https://www.accessdata.fda.gov/cdrh_docs/reviews/K151502.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K151502.pdf
- K160090: https://www.accessdata.fda.gov/cdrh_docs/reviews/K160090.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K160090.pdf
- K153607: https://www.accessdata.fda.gov/cdrh_docs/reviews/K153607.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153607.pdf

---

## 甲状腺球蛋白用于分化型甲状腺癌监测（Thyroglobulin, Tg）

### 1. 法规定位
- Product code **MSW**；21 CFR 866.6010；Class II；Panel Immunology (82 / IM)。
- openFDA 累计清关数：**total = 13**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K242981 | 2025-06-20 | Siemens / Atellica IM Thyroglobulin (Tg) | 吖啶酯 CLIA 夹心（双鼠单抗，biotin-streptavidin 预结合） | 血清、EDTA/Li-heparin 血浆 | aid in monitoring differentiated thyroid cancer patients who have undergone thyroidectomy with or without radioiodine ablation | 决策摘要 + summary |
| K221890 | 2023-09-30 | Roche / Elecsys Tg II（cobas e 411） | ECLIA 两步夹心 | 血清、Li-hep、K2/K3-EDTA 血浆 | aid in monitoring for the presence of persistent or recurrent/metastatic disease in patients who have DTC and have had thyroid surgery (with or without ablative therapy) | 决策摘要 + summary |
| K220972 | 2023-09-15 | Beckman Coulter / Access Thyroglobulin（biotin 改良） | 顺磁微粒 CLIA 一步夹心（4 种 biotin 化抗体预偶联） | 血清 | 同上 + "who lack serum thyroglobulin antibodies" | 决策摘要 + summary |
| K241423 | 2024-06-07 | Beckman / Access Thyroglobulin（新增 heparin 血浆） | 同上 | 血清、Li/Na-heparin 血浆 | 同上 | 决策摘要 + summary |
| K240927 | 2024-06-28 | Beckman / Access Thyroglobulin（DxI 9000, Lumi-Phos PRO） | 同上 | 血清、heparin 血浆 | 同上 | 决策摘要 + summary |

### 3. 预期用途与声明类型
- 全部为**术后监测**（persistent / recurrent / metastatic disease），Rx only，非 POC；Beckman 明确限定"lack serum thyroglobulin antibodies"，Siemens/Roche 以警示/黑框形式要求 anti-Tg 阳性样本不得检测 Tg。
- Roche 标签含 blackbox 声明：Tg 自身抗体可致假高/假低；不同方法结果不可直接比较，更换方法须平行检测。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
| 产品 | 临床 cutoff | 建立方式 | 健康参考区间 | 无病 DTC 患者分布 |
|---|---|---|---|---|
| K242981 Atellica Tg | **0.2 ng/mL**（阳性 = 术后复发证据）；分层：<0.2 Excellent、0.2–<1.0 Indeterminate、≥1.0 Biochemical Incomplete | **直接采用 2015 ATA 指南**"Excellent response"定义（非 TSH 刺激、无 TgAb 时 suppressed Tg <0.2） | 321 名健康者（157 男/164 女，22–80 岁）：中位 15.9，2.5th 2.4、**97.5th 74.9 ng/mL**（男 70.1、女 78.3） | 136 例术后 ≥4 年无病 DTC（排除 aTg >LoQ、妊娠、<22 岁）：95% ≤1.272 ng/mL |
| K221890 Elecsys Tg II | **0.2 ng/mL**；同样三分层 | 2015 ATA 指南 | 463 名健康者（244 男/219 女，22–79 岁）：中位 16.6，2.5th 3.6、**97.5th 77 ng/mL**（男 63.2、女 104） | 127 例 ≥4 年无病 DTC：0.1–11.6 ng/mL，80.3% <0.1，95% ≤0.786 |
| K220972 / K241423 / K240927 Access Tg | "Not applicable" / Refer to K002905 | — | predicate 152 名健康成人：1.15–130.77，中位 9.08，2.5th 1.59、**97.5th 50.03 ng/mL**；本次 28 名健康者 3 lot 验证 97.5th 44.76–46.44 | 文件未载明 |

- 要点：Tg 是本组中唯一**直接引用临床指南既定值（ATA 2015）**作为 cutoff 的靶点，且 cutoff 是绝对浓度而非 serial change；健康人群 97.5th（50–78 ng/mL）与临床 cutoff（0.2 ng/mL）相差两个数量级，因目标人群为甲状腺全切后。

### 5. 生物学 / 生理学依据
- 申报文件内：anti-Tg 自身抗体干扰（各文件警示；Access 用 80 份嗜异性血清测试，19 份有 HAMA 干扰；Atellica 6 份 HAMA 样本中 1 份 775 ng/mL 干扰）；RAI 治疗与非 RAI 患者分层的结构性疾病概率不同（Atellica：Tg ≥1.0 时 RAI 组 63.6%、非 RAI 31.6%；Roche：45%/22%）；交叉反应物 TBG、TSH、T3、T4、FSH、AFP、galectin-3、VEGF 无交叉。溯源统一至 **BCR CRM 457**。
- 背景（非申报文件）：Tg 仅由甲状腺滤泡细胞合成，全切+消融后应不可测；ATA 2015 将 Tg <0.2（非刺激）或 <1（刺激）定义为 excellent response。

### 6. 样本类型与样本要求
- K242981：血清（含 SST 胶管）、Li-hep（含 PST）、K2-EDTA，n=84–99 配对（0.05–145 ng/mL），Passing-Bablok 斜率 0.99–1.01；稳定性：血清 20–30°C 4 天、2–8°C 7 天；血浆 3 天/4 天；−20°C 12 个月、−70°C 24 个月；4 次冻融。
- K221890：Li-hep、K2/K3-EDTA n=65（分 0.197–19.9 与 20–490 两段），斜率 0.95–1.01；0.2 ng/mL 处 K3-EDTA 偏差 9.9%；稳定 2–8°C 与 15–25°C 14 天、−20°C 24 个月。
- K220972：仅血清；K241423 新增 Li-hep/Na-hep n=45（0.227–494），斜率 1.000/1.021，r 0.999（依 CLSI EP35）。
- 临床研究纳入：影像须在采血 ±30 天内；rhTSH 刺激 72 h 内排除；TgAb 阳性排除（Atellica 用 Beckman Access TgAb II >0.9 IU/mL；Roche 用 Elecsys anti-Tg ≥22 IU/mL）。

### 7. 分析性能验证所依据的标准
| K 号 | 标准 | 项目/结果 |
|---|---|---|
| K242981 | EP05-A3；EP06 2nd；EP07 3rd；EP09c；**EP12-A2**（临床定性性能）；EP17-A2；EP25-A；EP28-A3c；EP34；EP37；**I/LA30-A**（内源抗体干扰） | 精密度 7 水平 n=80，0.08 ng/mL 时 within-lab 9.0%，≥0.15 时 ≤3.3%；三站点 n=90 总 CV 1.9–5.8%；线性 0.04–163.8；稀释 1:10/20/50（EMI 至 7500）；hook 82,972 ng/mL 无；biotin 3510 ng/mL；LoB 0.039/LoD 0.044/LoQ 0.05；货架期 12 个月、在机 28 天；方法比对 vs Access Tg n=483，Passing-Bablok 斜率 **1.09** (1.06–1.12)、R 0.985 |
| K221890 | EP05-A3；EP06 2nd；EP07 3rd；EP09c；**EP12 3rd**；EP17-A2；EP37 | 精密度 n=84，0.13 ng/mL 时 within-lab 9.34%；三站点 3 lot 总 CV 6.3–13.6%；线性 0.005–578；hook 120,000；biotin 1200 ng/mL、HAMA 805 µg/L；LoB 0.02/LoD 0.04/LoQ 0.1；货架期 15 个月；方法比对 vs Access Tg n=126，**斜率 1.394、截距 −0.102**（系统性差异，补做 CRM 457 回收 91.7–110.5%） |
| K220972 | EP05-A3；EP06-Ed2；EP07-A3；EP09c；EP17-A2（稳定性用 EP25-A） | 精密度 0.11 ng/mL 时 within-lab 9.6%、0.17 时 15.9%；线性 0.05–563.7；hook 40,000；biotin 3510 ng/mL（原 10）；LoB 0.02/LoD 0.05/LoQ 0.05；方法比对 3 lot n=102 斜率 0.96–0.99 |
| K241423 | **EP35** | 矩阵等效（见第 6 节） |
| K240927 | EP05-A3；EP06-2nd；EP09c；EP17-A2 | 3 lot × 3 仪器精密度（0.3 ng/mL 时 within-lab 8.4%）；线性 0.03–642；LoB 0.03/LoD 0.05/LoQ 0.1；方法比对 vs Access 2 n=187 斜率 0.996、r 0.999 |

### 8. 临床验证设计与结果
- **K242981（Atellica Tg）**：前瞻性、美国 3 站点，407 例入组，排除 70 例，189 例（≥22 岁；全切/近全切 ± RAI 后 ≥6–12 周；ATA 初始风险低 37.6%/中 37.6%/高 24.9%；女性 70.9%），**291 次可评估访视**（影像在采血 30 天内）。金标准：结构性疾病（SD+）= 超声/CT/MRI 或 RAI 扫描/FDG-PET 阳性。Tg ≥0.2 判阳：SD+ 55 中 54 阳，SD− 236 中 110 阳 → **Sens 98.2% (94.6–100)、Spec 53.4% (47.8–58.0)**；按患病率 4.99% 调整 **NPV 99.8%、PPV 10.0%**（bootstrap 处理访视相关性，依 CLSI EP12-A2）。LR：RAI 组 Tg ≥1.0 LR 5.52，<0.2 LR 0.15；非 RAI 组 ≥1.0 LR 4.67。
- **K221890（Elecsys Tg II）**：前瞻性、美国 9 站点；纵向队列 219 例（术后 4–12 周入组，6/12/18/24 个月随访，5 次计划访视）+ 横断面队列 72 例（术后 >12 周且有结构性疾病，用于增加 SD+ 例数）；排除 242 份后 **530 份**（461 纵向 + 69 横断）；人群 68.5% 女性，中位 53 岁，RAI 20.4%，ATA 低/中/高 44.7/32.5/22.8%。Tg ≥0.2：SD+ 92 中 91 阳，SD− 438 中 204 阳 → **Sens 98.91% (94.10–99.81)、Spec 53.42% (48.74–58.05)**；以纵向队列真实患病率 4.99% (23/461) 调整 **NPV 99.89%、PPV 10.03%**。LR：RAI 组 ≥1.0 为 7.22、<0.2 为 0.03；非 RAI 组 ≥1.0 为 6.25。
- Access 系列（K220972/K241423/K240927）：临床 Refer to K002905（2000 年清关，文件未获取）。

### 9. 厂家间差异与要点
- 三家 AMI：Atellica 0.05–150（EMI 7500）；Elecsys 0.1–500（EMI 5000）；Access 0.1–500。LoQ：0.05 / 0.1 / 0.05–0.1 ng/mL——均满足 ATA 0.2 ng/mL 阈值所需的"高敏"要求。
- 尽管均溯源 CRM 457，Elecsys Tg II vs Access 斜率 1.394、Atellica vs Access 1.09，表明**方法间偏差仍显著**，FDA 接受了 CRM 回收实验作为补充证据。
- Siemens 与 Roche 2023–2025 两份新器械文件采用几乎相同的临床设计（前瞻、SD+ 影像金标准、ATA 三分层、按真实患病率调整 PPV/NPV），Roche 额外用横断面富集队列。
- Beckman 三份文件为改良型（biotin、血浆、DxI 9000），不含临床数据。

### 10. 来源
- K242981: https://www.accessdata.fda.gov/cdrh_docs/reviews/K242981.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K242981.pdf
- K221890: https://www.accessdata.fda.gov/cdrh_docs/reviews/K221890.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K221890.pdf
- K220972: https://www.accessdata.fda.gov/cdrh_docs/reviews/K220972.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220972.pdf
- K241423: https://www.accessdata.fda.gov/cdrh_docs/reviews/K241423.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K241423.pdf
- K240927: https://www.accessdata.fda.gov/cdrh_docs/reviews/K240927.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K240927.pdf

---

## 附：本组横向对比速览（均取自上述文件）

| 靶点 | 参考上限（典型） | 建立方式 | 显著变化定义 | 监测研究 Sens / Spec（代表） |
|---|---|---|---|---|
| CEA | 非吸烟 3.0 / 吸烟 5.0 ng/mL | 健康人群非参数百分位（n=347；n=250） | RCV 36.2%（Vista） | 55.2% / 83.6%（K071603） |
| CA 125 | 30.8–38.1 U/mL（35 传统） | 健康女性 95th–97.5th（n=240；n=200） | ≥20%（1.645√2·CV）或 ≥25%（>35 时） | 67.3% / 76.0%（K142895） |
| CA 19-9 | 35–37 U/mL；Lumipulse 97.5th 50 | 健康人群分布（n=40 验证；n=300；n=240） | >15%（Siemens/Fujirebio）；RCV 84.7%（Vista） | 68.3% / 60.2%（K253528）；19.2% / 88.8%（K100375） |
| PSA（监测） | 95th 1.68–2.92 ng/mL（≥50 岁健康男） | 百分位 | RCV 50%（Siemens）；>20% 或 >0.2 ng/mL（FREND） | 54.5% / 87.6%（K251630）；63.0% / 81.0%（K162378） |
| AFP（睾丸癌） | 6.9–8.0 IU/mL 或 ng/mL | 97–97.5th（n=140–408） | RCV 33.7%（Vista） | 27.1% / 87.0%（K071597） |
| HE4 | 70/140 pmol/L（绝经前/后）或 150 pM | 95th/97.5th（n=400；n=204；n=345） | ≥14/18/20/25%（2–2.5×CV） | 53.5% / 78.5%（K093957） |
| ROMA | 分数 1.31/2.77（Fujirebio/Abbott）或 1.14/2.99（Roche） | 训练集 logistic，specificity 75% 设点 | 不适用（单次定性） | 联合 ICRA：Sens 90.9%、NPV 97.7%（K103358） |
| Tg | 0.2 ng/mL（ATA 2015）；健康 97.5th 50–78 | 指南既定值 | 不适用（绝对阈值 + 三分层） | 98.2% / 53.4%（K242981）；98.9% / 53.4%（K221890） |
