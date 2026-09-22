# 靶点组 F：感染性疾病血清学抗体检测（FDA 510(k) 申报信息整理）

> 数据来源说明：以下所有申报数据均摘自 FDA 公开的 510(k) Decision Summary（`cdrh_docs/reviews/<K>.pdf`）或 510(k) Summary（`cdrh_docs/pdfYY/<K>.pdf`）。文件中未载明的信息一律标注「文件未载明」。凡来自笔者背景知识（指南、生理学解释）的内容均单独标注为「背景（非申报文件）」。openFDA 累计清关数为 `fda_fetch.py list <CODE>` 命令首行 `total=` 的数值（查询日期 2026-09-22）。

---

## 单纯疱疹病毒 2 型特异性 IgG（HSV-2 type-specific IgG；兼看 HSV-1 IgG）

### 1. 法规定位
- Product code **MYF**（Enzyme linked immunosorbent assay, Herpes Simplex Virus, HSV-2）；HSV-1 型特异性 IgG 为 **MXJ**；双抗原 immunoblot 型（HerpeSelect 1&2 Immunoblot，K000238）为 **LGC**。
- 21 CFR **866.3305** Herpes simplex virus serological assays；**Class II（Special Controls）**；Panel：Microbiology (83 / MI)。
- 适用 Special Controls 指导原则：*Class II Special Controls Guidance Document: Herpes Simplex Virus Types 1 and 2 Serological Assays*（K081687 引用 2007-04-03 版；K103603 引用 2010-09-28 版；K220924 引用 2011-08-09 版；K243575 引用同名文件未注明日期）。
- openFDA 累计清关数：**MYF total = 10**；MXJ total = 20。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K243575 | 2025-02-12 | Biokit S.A.（Abbott 分销）/ ARCHITECT HSV-2 IgG | CMIA（两步法，重组 gG2 抗原包被微粒，acridinium 标记抗人 IgG） | 血清、SST、K2-EDTA、Li-heparin、Li-heparin PST 血浆 | sexually active adults / expectant mothers；aid in presumptive diagnosis；预测值取决于流行率；未用于儿科、新生儿、免疫抑制；未清关用于献血者筛查 | 决策摘要 + 510(k) summary |
| K220924（原始清关 K121895，2012） | 2022-10-12 | Roche / Elecsys HSV-2 IgG | ECLIA 双抗原夹心（生物素化与钌标记重组 HSV-2 特异抗原） | 血清、Li-heparin、K2/K3-EDTA 血浆 | sexually active individuals / expectant mothers；不用于献血者筛查；儿科、新生儿、免疫抑制、POC 未建立 | 决策摘要 + 510(k) summary（K220924 为生物素耐受改良，性能引用 K121895） |
| K181334 | 2018-08-23 | Biokit S.A.（Siemens 分销）/ ADVIA Centaur Herpes-2 IgG | 间接化学发光两步夹心（重组 gG2 固相，acridinium ester 标记单抗抗人 IgG） | 血清、EDTA 及 Li-heparin 血浆 | sexually active adults / expectant mothers；不用于献血者筛查 | 仅 510(k) summary（决策摘要未获取到） |
| K081687 | 2008-11-10 | DiaSorin / LIAISON HSV-2 Type Specific IgG | 间接 CLIA（磁微粒包被重组 gG2，isoluminol 标记鼠单抗抗人 IgG） | 血清 | sexually active adults / expectant mothers；未用于儿科、新生儿、免疫抑制；不用于献血者 | 决策摘要 |
| K103603 | 2011-05-20 | ZEUS Scientific / ZEUS ELISA HSV gG-2 IgG Test System | 手工 ELISA（亲和纯化 **天然** gG-2 抗原，HRP 标记羊抗人 IgG Fc） | 血清 | sexually active individuals / pregnant women；不用于 donor screening 或自测 | 决策摘要 + 510(k) summary |
| K120959（原始 K090409） | 2012-07-25 | Bio-Rad / BioPlex 2200 HSV-1 & HSV-2 IgG | 多重微珠流式免疫分析（Luminex），重组 gG1（55 kD）与 gG2（31 kD）分别包被两种微珠 | 血清、EDTA 或肝素血浆 | 同时检测并区分 HSV-1/HSV-2 IgG；sexually active individuals / expectant mothers | 决策摘要 + 510(k) summary（本次为改良申报，临床性能引用 K090409） |

### 3. 预期用途与声明类型
- 全部为 **处方用（Rx only）定性辅助诊断**：措辞统一为 "aid in the presumptive diagnosis of HSV-2 infection"，目标人群限定为 **sexually active adults/individuals 与 expectant mothers/pregnant women**。
- 共同限定语：(1) 预测值取决于人群流行率与验前概率（"The predictive value of a positive or negative result depends on the population's prevalence and the pretest likelihood"）；(2) 未在儿科、新生儿筛查、免疫抑制/免疫功能低下患者中建立性能；(3) **not FDA cleared for screening blood or plasma donors**；(4) K243575、K220924 另加 "test results may not determine the state of active lesions or associated disease manifestations, particularly for primary infection"。
- K103603 额外声明 "not intended for donor screening or for self testing"。
- 无任何产品声明用于免疫状态判断或人群筛查；BioPlex（K120959）声明为 "qualitative detection and differentiation of IgG antibodies to HSV-1 and HSV-2"。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 结果单位 | 判定规则 | Equivocal 区 | cutoff 建立方式（文件所载） |
|---|---|---|---|---|
| ARCHITECT HSV-2 IgG（K243575） | S/CO | <1.00 Nonreactive；≥1.00 Reactive | **无** equivocal 区 | 505 份血清（271 反应性、228 非反应性、6 份 equivocal）做 **ROC 分析**，选定 1.00 S/CO，再在多中心临床一致性研究中验证 |
| Elecsys HSV-2 IgG（K220924/K121895） | COI（signal/cutoff） | <1.0 Non-reactive；≥1.0 Reactive | 无 | 最初用 **267 份血清**（sexually active adults 与 pregnant women 两队列）测定，将阳/阴性分布与 predicate（FDA 清关 immunoblot）比较后设定 cutoff；K121895 510(k) summary 载 LoB 0.11 COI、LoD 0.1244 COI、hook effect ≤10 COI 无 |
| ADVIA Centaur Herpes-2 IgG（K181334） | Index | 以 Index 报告；文件未载明具体判定阈值数值（510(k) summary 仅给出各 Index 区间的精密度设计要求 0.51–0.79、0.80–1.20 等） | 文件未载明 | 文件未载明 |
| LIAISON HSV-2 Type Specific IgG（K081687） | Index | cutoff Index 1.0 | **>0.90 且 <1.1 为 equivocal** | 决策摘要 "Assay cut-off: NA"，建立方法文件未载明 |
| ZEUS ELISA HSV gG-2 IgG（K103603） | Index Value / OD ratio | ≤0.90 Negative；≥1.10 Positive | **0.91–1.09 Equivocal**（须复测双份；仍 equivocal 则用 Western blot 等替代方法或 1–3 周后复采） | 用 **25 份确认阴性 + 9 份临床特征化阳性血清**，计算阴性人群均值与 SD，"using a mathematical calculation based on these data" 得理论 cutoff 并用特征化标本验证（具体倍数 SD 文件未载明） |
| BioPlex 2200 HSV-1 & HSV-2 IgG（K120959） | Antibody Index (AI) | 文件（改良申报）称 cutoff 保持不变，结果分 Positive/Equivocal/Negative；精密度表中 "Near Cutoff" 样本 AI≈1.0–1.2 | 有 equivocal 类别，具体区间文件未载明（引用 K090409） | 引用 K090409 |

- FDA 对 **低流行人群 PPV** 的关注在申报文件中的体现：
  - K081687（LIAISON）决策摘要列出 "HSV-2 Prevalence vs. Hypothetical Predictive Values" 表：以灵敏度/特异度 98.1%/98.0%（sexually active）计算，流行率 10% 时 PPV 仅 84.5%（expectant mothers 79.6%）。
  - K103603（ZEUS）同样列表：流行率 5% 时 sexually active PPV 92.9%、pregnant women PPV 82.5%；package insert 限制条款要求 "False positive test results may occur. Repeat testing or testing with a different device may be indicated in some settings e.g., patients with low likelihood of HSV infection"。
  - K243575、K103603、K081687、K121895 均专门设置 **Low Prevalence 人群**（16–19 岁 / 17–19 岁非 STD 场所个体，或非 STD 门诊患者）做特异度/NPA 研究（见第 8 节）。

### 5. 生物学 / 免疫学依据
- **抗原选择（申报文件）**：除 ZEUS（亲和纯化天然 gG-2）外，各厂家均使用 **重组 glycoprotein G-2（gG2）** 作为型特异性抗原；BioPlex 同时用重组 gG1（55 kD）与 gG2（31 kD）双微珠实现分型；Elecsys 用 "HSV-2-specific recombinant antigens"（大肠杆菌表达，具体抗原文件未细述）。predicate HerpeSelect immunoblot 则为 HSV-1/HSV-2 天然抗原 + 重组 gG1 + 重组 gG2 混合条带。
- **类特异性（申报文件）**：K243575 做了 class-specificity 研究证明仅检测 IgG；K103603 用去除总 IgG 的 HSV-2 IgM 阳性标本（0/10 阳性）证明不检测 IgM。
- **gG 缺陷株导致假阴性（申报文件）**：K121895 限制条款载 "False negative results may occur when the HSV virus is glycoprotein G (gG) deficient (0.2 % HSV isolates were gG deficient)"。
- **交叉反应（申报文件）**：主要评估 HSV-1 IgG 阳性标本（K081687：383 份 HSV-1 IgG 阳性，LIAISON 4 阳性/2 equivocal；K181334：117 份 HSV-1，Centaur 反应 29 vs 参考 31，即多为真阳性）；其他包括 EBV、CMV、VZV、HHV-6/8、Parvovirus B19、RF、ANA、Treponema pallidum、Chlamydia、Candida、Gardnerella、Toxoplasma 等。K243575 表中 HHV-8 IgG 2/10、Parvovirus B19 2/13、anti-dsDNA 1/8 反应（未说明是否真阳性）。K103603 注明 Candida albicans 未评估并写入限制条款。
- 背景（非申报文件）：gG1 与 gG2 的氨基酸序列同源性低，是唯一可靠区分 HSV-1/HSV-2 抗体的糖蛋白，故 CDC/ASM 推荐仅接受 gG-based 型特异性血清学；IgG 血清转换通常在原发感染后 2–12 周，故急性期单次阴性不能排除感染（申报文件限制条款亦有类似表述）。

### 6. 样本类型与样本要求
- **ARCHITECT（K243575）**：血清、SST、K2-EDTA、Li-heparin、Li-heparin PST；61 组配对基质回归 y=1.00–1.03x，r 0.996–0.998。样本稳定性：≤-20 °C 1 个月，2–8 °C 离凝块 14 天/在凝块 7 天，室温 48 h，3 次冻融。
- **Elecsys（K220924/K121895）**：血清、Li-heparin、K2/K3-EDTA；枸橼酸钠血浆未评估；样本稳定性研究仅用血清（K121895 限制条款）。
- **ADVIA Centaur（K181334）**：血清、SST、EDTA、Li-heparin；68 组配对 Deming 回归斜率 0.98–0.99，r 0.997–0.999。
- **LIAISON（K081687）、ZEUS（K103603）**：仅血清；ZEUS 干扰研究显示阴性样本在高胆红素（信号变化 +272–318%）、高血红蛋白等条件下信号变化 >20% 但仍低于 cutoff，未改变定性结果。
- **BioPlex（K120959）**：血清、EDTA、肝素钠血浆；44–47 组配对回归斜率 0.95–1.00，r ≥0.99。

### 7. 分析性能验证所依据的标准

| K 号 | 文件 "Standards/Guidance Documents Referenced" 所列 | 对应项目 |
|---|---|---|
| K243575 | CLSI EP05-A3 (2014)；CLSI EP07 3rd ed. (2018)；CLSI EP09c 3rd ed. (2018)；CLSI EP12-A2 (2008)；CLSI EP25-A (2009)；IEC/EN 61010-1:2010；CLSI I/LA21-A2 *Clinical Evaluation of Immunoassays*；HSV Special Controls Guidance | 精密度（20 天、12 天、3 中心重复性）；内源/药物干扰；基质比较；定性一致性；试剂/样本稳定性；仪器电气安全；临床评价 |
| K220924 | CLSI EP05-A3；CLSI EP07 3rd ed.；CLSI EP09c 3rd ed.；HSV Special Controls Guidance (2011-08-09) | 精密度（21 天，n=84）；生物素干扰（≤1200 ng/mL 偏差 ≤10%）；与原版本 206 份样本方法比较 |
| K181334（510k summary） | CLSI EP05-A3（精密度、多中心重复性）；CLSI EP7-A2（干扰） | 20 天精密度 80 重复；3 中心 5 天；干扰 ≤10% |
| K081687 | HSV Special Controls Guidance (2007-04-03)；*Guidance for the Content of Premarket Submissions for Software Contained in Medical Devices* (2005-05-11) | 3 中心重复性（8 份工程化血清，5 天 × 2 run × 4 rep）；交叉反应 840 份；干扰 |
| K103603 | HSV Special Controls Guidance (2010-09-28)；CLSI EP5 2nd ed.（EP5-A2）；CLSI EP7-A2 | 20 天精密度 80 重复；3 中心 180 重复；14 类交叉反应各 10 份 |
| K120959 | CLSI EP05-A2；CLSI EP07-A2；CLSI EP09-A2（用于基质比较）；CLSI EP12-A2；"EN 1360:200"（原文如此，疑为 EN 13640 稳定性标准，笔者推测） | 5 天精密度；14 种干扰物各 10 重复；基质比较；定性一致性；稳定性 |

### 8. 临床验证设计与结果

**ARCHITECT HSV-2 IgG（K243575）**——多中心（3 家美国外部实验室）前瞻性采集 915 份（sexually active + pregnant；女 670、男 245，年龄 14–99）；**参考方法为 3 种市售 anti-HSV-2 IgG 检测组成的复合比较法（2/3 多数规则）**（未用 Western blot）。
- Sexually active（n=618）：PPA 96.54%（223/231；95% CI 93.32–98.24），NPA 96.90%（375/387；94.66–98.22）。
- Pregnant（n=297）：PPA 95.12%（78/82；88.12–98.09），NPA 98.60%（212/215；95.98–99.52）。
- CDC panel（50 份血清 ×2 = 100 份盲样）：PPA 100%（30/30），NPA 97.14%（68/70）。
- Low prevalence（139 份 16–19 岁非 STD 场所）：NPA 98.48%（130/132；Wilson 94.64–99.58）。
- 观察到 HSV-2 阳性率 35%（女 36%、男 31%）。

**Elecsys HSV-2 IgG（K121895，经 K220924 引用）**——参考方法为 predicate immunoblot，并对 Western blot 计算相对灵敏度/特异度。
- Expectant mother（n=125）：PPA 97.83%（88.47–99.94），NPA 98.73%（93.15–99.97）。
- Sexually active（n=469）：PPA 93.63%（88.60–96.90），NPA 98.72%（96.75–99.65）。
- Low prevalence（n=200）：PPA 75.00%（19.41–99.37，阳性例数极少），NPA 98.47%（95.59–99.68）。
- 相对 Western blot：灵敏度 pregnant 100%、sexually active 99.3%、low prevalence 82.4%；特异度 93.1%、95.1%、100%。
- CDC panel 一致性 100%。K220924 本身为生物素改良，仅做 206 份样本与原版 100% 定性一致。

**ADVIA Centaur Herpes-2 IgG（K181334）**——864 份（≥18 岁，含 274 孕妇），3 家外部实验室；比较法为市售 anti-HSV-2 IgG immunoblot；比较法 equivocal（22 份）送 **University of Washington（Seattle）Western blot** 裁定（20 份阴性、2 份仍 equivocal）。
- 总体：灵敏度 95.3%（245/257；92.0–97.3），特异度 98.5%（598/607；97.2–99.2），总一致 97.6%。
- 孕妇（n=274）：灵敏度 100%（34/34；89.9–100），特异度 98.3%（236/240；95.8–99.4）。
- CDC panel 100 份 100% 一致；ZeptoMetrix ToRCH 混合盘 24 份 100% 一致。
- 交叉反应/特殊人群 522 份总一致 96.9%（506/522）。

**LIAISON HSV-2（K081687）**——951 份美国东北部样本，2 家外部实验室；比较法为 FDA 清关 immunoblot，predicate 重复 equivocal 送 "Reference Laboratory in the Pacific Northwest" 做 Western blot。
- Sexually active/STD 门诊（n=401）：灵敏度 98.1%（104/106；95.6–99.9），特异度 98.0%（289/295；96.0–99.1）。
- Expectant mothers（n=430）：灵敏度 94.8%（91/96；89.4–97.9），特异度 97.3%（325/334；95.3–98.6）。
- Low prevalence（n=120，非 STD 门诊）：灵敏度 100%（20/20），特异度 100%（100/100）。
- CDC panel（52% 阳性/48% 阴性）100% 一致（52/52，48/48）。

**ZEUS ELISA HSV gG-2 IgG（K103603）**——3 中心，788 份；比较法为市售 HSV-1/2 immunoblot；equivocal 不一致按不利于受试产品计。
- Sexually active（n=336）：PPA 93.3%（70/75；86.9–98.5），NPA 99.6%（260/261；97.9–100）。
- Pregnant（n=252；孕早/中/晚 121/64/67）：PPA 98.6%（70/71；92.5–100），NPA 98.9%（179/181；96.1–99.9）。
- Low prevalence（n=100，17–19 岁非 STD）：NPA 98.0%（98/100；2 份 equivocal），PPA 不适用。
- CDC panel（46 阳/54 阴，其中 24 份 HSV-1/2 双阳）：PPA 100%（46/46），NPA 96.9%（文件表格记 54/54 与 96.9% 并存，原文如此）。

**BioPlex 2200 HSV-1 & HSV-2 IgG（K120959）**——改良申报，仅做与已清关 BioPlex 的一致性（399 份 sexually active）：HSV-2 PPA 99.4%（166/167），NPA 100%（232/232）；CDC panel 80 份：HSV-2 PPA/NPA 均 100%（40/40，40/40），HSV-1 PPA 100%（42/42）、NPA 97.4%（37/38）。原始灵敏度/特异度需查 K090409（本次未获取）。

### 9. 厂家间差异与要点
1. **Equivocal 区策略分化**：老一代 ELISA/CLIA（ZEUS 0.91–1.09；LIAISON >0.90–<1.1）设 ±10% 灰区并要求复测/Western blot；新一代全自动 CMIA/ECLIA（ARCHITECT、Elecsys）以 ROC 或 predicate 对齐后**不设灰区**（ARCHITECT 精密度 Serum Panel 1 均值 0.95 S/CO，总 CV 以 SD 报告）。
2. **参考方法演变**：2008–2018 年清关产品以 FDA 清关 immunoblot 为比较法、以 UW Western blot 裁定 equivocal（K081687、K181334、K121895）；2025 年 ARCHITECT 改用 **3 种市售 IgG 检测 2/3 复合比较法**，并不再要求 Western blot。
3. **低流行人群**：所有决策摘要均包含 16–19 岁或非 STD 门诊人群 NPA（98.0–100%），并（K081687、K103603）附 PPV-流行率假设表，直接回应 Special Controls 指南对低流行人群假阳性的关注。
4. **抗原**：ZEUS 为唯一用天然亲和纯化 gG-2 者；其余均为重组 gG2；BioPlex 借多重微珠同步分型。
5. **样本基质**：仅血清（DiaSorin 2008、ZEUS）→ 血清+EDTA/肝素血浆（Bio-Rad、Siemens、Roche）→ 血清/SST/K2-EDTA/Li-heparin/PST 全覆盖并附 EP09c 回归（Abbott 2025）。
6. **性能差异**：sexually active 人群 PPA 范围 93.3%（ZEUS）–98.1%（LIAISON）；孕妇 PPA 94.8%–100%；NPA 普遍 ≥96.9%。

### 10. 来源
- K243575：https://www.accessdata.fda.gov/cdrh_docs/reviews/K243575.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K243575.pdf
- K220924：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220924.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220924.pdf
- K121895：https://www.accessdata.fda.gov/cdrh_docs/reviews/K121895.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K121895.pdf
- K181334：https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181334.pdf（决策摘要 reviews/K181334.pdf 未能下载）
- K081687：https://www.accessdata.fda.gov/cdrh_docs/reviews/K081687.pdf
- K103603：https://www.accessdata.fda.gov/cdrh_docs/reviews/K103603.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf10/K103603.pdf
- K120959：https://www.accessdata.fda.gov/cdrh_docs/reviews/K120959.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K120959.pdf

---

## 莱姆病伯氏疏螺旋体抗体（*Borrelia burgdorferi* antibodies, Lyme disease serology）

### 1. 法规定位
- Product code **LSR**（Reagent, Borrelia Serological Reagent）；控制品另编 QCH。
- 21 CFR **866.3830**（法规名称为 *Treponema pallidum treponemal test reagents*，Lyme 血清学试剂历史上归入此条）；**Class II**；Panel：Microbiology (83 / MI)。
- 相关 FDA 指导原则（K233367 引用）：*Establishing the Performance Characteristics of in Vitro Diagnostic Devices for the Detection of Antibodies to Borrelia burgdorferi – Guidance for Industry and FDA Staff*。
- openFDA 累计清关数：**LSR total = 112**。
- 注：用户提示的 K182192/K182193 经核查非 Lyme 产品（K182193 为正畸托槽，K182192 无文件）；ZEUS 的 MTTT 修改用途清关为 **K190907 / K191240 / K191398**（均 2019-07-29）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K183446 | 2019-03-12 | Bio-Rad / BioPlex 2200 Lyme Total | 多重微珠流式免疫（重组 p58、OspC type B、合成肽 FVlsE［FlaB+VlsE 序列］；抗人 IgG+IgM PE 标记） | 血清、EDTA/肝素血浆 | 总抗体一线筛查；reactive/equivocal 须做 second tier（Western blot） | 决策摘要 |
| K191240 | 2019-07-29 | ZEUS Scientific / ZEUS ELISA Borrelia VlsE1/pepC10 IgG/IgM + ZEUS ELISA B. burgdorferi IgM | 手工 ELISA；VlsE1+pepC10 肽 / 全菌抗原（IgM 版含抗人 IgG 沉淀去除 IgG 与 RF） | 血清 | **MTTT** 修改用途：两个 EIA 串联替代 Western blot | 决策摘要 |
| K190907 | 2019-07-29 | ZEUS / VlsE1/pepC10 IgG/IgM + B. burgdorferi IgG/IgM（全菌） | 同上 | 血清 | MTTT（IgG/IgM） | 决策摘要 |
| K202574（配套 K202573 IgM） | 2021-02-18 | DiaSorin / LIAISON Lyme IgG（+ LIAISON Lyme Total Antibody Plus 修改用途） | 间接 CLIA，磁微粒包被重组 VlsE（B. burgdorferi B31、B. garinii Pbi）；IgM 版为 OspC（B. afzelii pKo）+VlsE（B31） | 血清、SST、K2-EDTA、Li-heparin | 一线（STTT）或二线（MTTT，与 Lyme Total Antibody Plus 串联） | 决策摘要 |
| K203289 | 2021-03-22 | Gold Standard Diagnostics / Borrelia burgdorferi VlsE-OspC IgG/IgM ELISA（捆绑 IgG/IgM、IgG、IgM 全菌 ELISA 修改用途） | 手工 ELISA，重组 VlsE + OspC（B31，IMAC 纯化） | 血清 | 一线或二线（MTTT/STTT） | 决策摘要 |
| K173496 | 2018-08-30 | Quidel / Sofia 2 Lyme FIA | 免疫荧光侧向层析双向试条（IgM/IgG 分侧；重组蛋白+合成肽） | **指尖全血**（CLIA waived） | ≥2 周症状疑似者辅助诊断；阳性须 Western blot 确认 | 决策摘要 |
| K233367 | 2024-08-12 | ID-FISH Technology / iDart Lyme IgG ImmunoBlot | 线性重组抗原免疫印迹（P93、P41、P39、P23、P31、P66、P58、P45、P34、P30、P28、P18 + LSA［嵌合 VlsE 肽］） | 血清 | **MTTT 二线**免疫印迹；不用于无症状者筛查 | 决策摘要 |

### 3. 预期用途与声明类型
- 均为 **Rx 定性辅助诊断**，限用于 "patients with signs and symptoms consistent with / suspected of Lyme disease"；K233367 明确 "not intended as a screen for asymptomatic patients"。
- **两步法声明（申报文件）**：
  - 2019 年前（K183446、K173496）：所有 reactive/equivocal 须以 "second tier test such as Lyme IgG and IgM Western blot" 补充；"Non-reactive first tier or negative second tier results should not be used to exclude borreliosis"。
  - 2019 年起（K190907/K191240、K202574、K203289、K233367）：intended use 明文列出两条路径——**(1) STTT** "using IgG or IgM Western blot testing following current interpretation guidelines"；**(2) MTTT** 用指定的第二个 EIA/CLIA/免疫印迹。"Positive test results by either the STTT or MTTT methodology are supportive evidence for the presence of antibodies and exposure to B. burgdorferi"。
- 背景（非申报文件）：CDC 于 2019 年 MMWR（Mead P. et al., *MMWR* 68(32):703）更新推荐，承认以 FDA 清关的两个 EIA 串联（MTTT）作为 STTT 的可接受替代；上述 ZEUS 三项清关即为首批 MTTT 用途清关。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 单位 | 判定 | Equivocal | cutoff 建立（文件所载） |
|---|---|---|---|---|
| BioPlex 2200 Lyme Total（K183446） | AI | ≤0.8 Non-reactive；≥1.1 Reactive | **0.9–1.0 AI** | **1,372 份正常血清**按千分位百分位数扫描，取约 **第 98 百分位** RFI 定为 1.0 AI，再用临床诊断样本验证；灰区用于包夹 cutoff |
| ZEUS VlsE1/pepC10 IgG/IgM（K190907/K191240） | Index Value (IV) | 文件（修改用途申报）"Assay Cut-Off: N/A"；重复性样本 High− 0.74、Low+ 1.36 | 有 equivocal（临床数据分 positive/equivocal） | 引用原清关 K113397 / K885317 / K900196 |
| LIAISON Lyme IgG（K202574）/ IgM（K202573） | Index（量程 0.01–45.0） | cutoff 1.0；负值为 Negative，正值 Positive | **0.90–1.10 equivocal**（"to account for normal measurement imprecision"） | 用美国来源特征化 Lyme 患者样本 + 常规送检样本，结合其他已清关血清学结果分类为预期阳/阴，选 1.0 为 "best balance of sensitivity and specificity" |
| GSD VlsE-OspC IgG/IgM（K203289） | Units = (OD/cutoff)×10；cutoff = Calibrator OD × Correction Factor | <9.0 Negative；>11.0 Positive | **9.0–11.0 Units** | 238 份正常血清（101 流行区 + 95 非流行区，文件如此）取 **均值 + 3 SD**；再用 197 份（125 份各期 Lyme：79 份来自 Dr. Steere、46 份 CDC；72 份 CDC 交叉反应病例）做 **ROC** 确认 |
| Sofia 2 Lyme FIA（K173496） | 仪器判读 Positive/Negative（IgM、IgG 各一） | 无数值 | 无 | 用配对指尖全血/血清/血浆测定，使全血 cutoff 的 PPA/NPA 与已清关血清/血浆版本统计学相当，再在临床试验中验证 |
| iDart Lyme IgG ImmunoBlot（K233367） | 条带强度以 C2（Protein L）对照带为基准；≥C2 记 "+" | **Positive = LSA 带阳性 且 P93、P41、P39、P23、P31、P34 中至少两组各有 ≥1 带**；否则 Negative | 无 | 文件未载明条带规则建立过程 |

- **STTT Western blot 判读规则（申报文件）**：K190907 明确 "A positive B. burgdorferi IgM Western blot result was defined by at least **2 out of 3** positive antibody bands and IgG Western blot result was defined by at least **5 out of 10** positive antibody bands"。背景（非申报文件）：此即 CDC/ASTPHLD 1995 采纳的 Dressler/Engstrom 标准（IgM：23/24、39、41 kDa 中 2 条；IgG：18、21/23、28、30、39、41、45、58、66、93 kDa 中 5 条）。
- **MTTT 判定**：第一 EIA 阳性或 equivocal → 第二 EIA；第二 EIA 阳性**或 equivocal** 即判 MTTT 阳性（K190907、K191240、K202574、K203289 均如此规定）。

### 5. 生物学 / 免疫学依据
- **抗原设计（申报文件）**：
  - 一线高灵敏抗原：VlsE 或其 C6/pepC10/FVlsE 肽（Bio-Rad、ZEUS、DiaSorin、GSD、iDart LSA），常与 **OspC**（早期 IgM 主要靶抗原：Bio-Rad OspC type B、DiaSorin IgM 用 B. afzelii pKo OspC、GSD B31 OspC）或 p58 联用。
  - 二线：全菌超声裂解物 ELISA（ZEUS B31 whole cell、GSD B31+2591 株）或多抗原免疫印迹（iDart 12 种重组蛋白）。
  - 免疫印迹 K233367 中 LSA 带为嵌合 VlsE 肽，作为必要条件。
- **IgM/IgG 动力学（申报文件数据体现）**：各期分层灵敏度显示急性/EM 期一线阳性率仅 60–80%，晚期近 100%（见第 8 节）；ZEUS IgM ELISA 加入抗人 IgG 沉淀以去除 IgG 及 RF（防止 RF 假阳性与 IgG 竞争）；DiaSorin 用 **DTT 灭活 IgM** 证明 IgG 检测类特异性（IgM 信号降至 <0.13，IgG 变化 ≤9.7%）。
- **交叉反应机制（申报文件）**：多个文件显示与 **梅毒（螺旋体共同抗原）、其他蜱传病（Babesiosis、Ehrlichiosis、tick-borne relapsing fever）、Leptospirosis、EBV/传染性单核细胞增多症、H. pylori** 的交叉最常见（K203289：TBRF 3/4、Ehrlichiosis IgM 3/6、Babesiosis 6/16、Leptospirosis 2/10，且 predicate 亦阳性；K202574：Babesiosis 5/10；K183446：EBV NPA 91.8%、Ehrlichiosis 90.9%；K173496 Sofia IgM 在 Lupus 8/20、Syphilis 8/28、ANA 4/10 阳性）。背景（非申报文件）：p41 鞭毛蛋白与其他螺旋体/细菌鞭毛高度同源，是 IgM 假阳性主要来源，故 STTT 要求 IgM 至少 2 条带。

### 6. 样本类型与样本要求
- 血清为所有产品共有；BioPlex（EDTA、肝素钠/锂，配对 74–78 组，斜率 0.97–1.02，r 0.951–0.981）与 LIAISON（SST、K2-EDTA、Li-heparin，32 组 Passing-Bablok，比例偏差 0.96–1.04）另支持血浆。
- LIAISON 样本稳定性：室温 3 天、2–8 °C 14 天、-20 °C（3 个月数据，研究进行中）、5 次冻融。GSD：2–8 °C 7 天，3 次冻融。iDart：冻存（-20 °C 9–44 天）与新鲜结果 100% 一致（72 份）。
- Sofia 2 为**指尖全血**（CLIA-waived，CW170015），与血清/血浆 321 组配对阳性率相当（IgM 42.4/42.1/41.4%）。
- ZEUS、GSD、iDart 仅血清，不支持血浆（文件未做基质研究）。

### 7. 分析性能验证所依据的标准

| K 号 | 文件所列标准/指南 | 说明 |
|---|---|---|
| K183446 | "Standard/Guidance Document Referenced: N/A"；正文引用 CLSI EP5-A3（精密度 20 天）、EP15-A3（3 中心重复性）、EP07-A2（干扰） | 精密度 80 重复；120 重复/中心；10 种干扰物 |
| K191240 / K190907 | N/A（修改用途申报） | 引用原始 3 中心 5 天重复性、20 天重复性；13 类交叉反应各 10 份；6 种干扰物 |
| K202574 / K202573 | CLSI EP05-A3；CLSI EP15-A3；CLSI EP07-A3（文件写作 "Second Edition"，原文如此） | 12 天两批精密度；3 中心；222 份/22 种交叉反应；7 种干扰物；DTT 类特异性；冻融/室温/4 °C/-20 °C 稳定性 |
| K203289 | "None"；正文引用 CLSI EP7-A3 干扰浓度 | 12 天精密度；3 批批间；3 中心；238 份健康人；226 份交叉反应 |
| K173496 | N/A | 12 天精密度；3 中心；200 份健康人；17 类交叉反应；17 种干扰物（含抗生素） |
| K233367 | FDA *Establishing the Performance Characteristics of IVD Devices for the Detection of Antibodies to B. burgdorferi* 指南 | 3 中心 6 操作者 90 重复；376 份交叉反应；5 种内源干扰 |

### 8. 临床验证设计与结果

**统一设计要素**：(1) 流行区前瞻性送检人群（"test-ordered"）与 predicate 一线一致性 + Western blot 二线一致性；(2) **CDC Lyme Serum Repository 280 份盲样**（Stage I/早期 EM 急性 & 恢复期、Stage II 播散［心脏/神经］、Stage III 晚期关节炎/神经；90 份 look-alike 疾病［梅毒、RA、MS、纤维肌痛、传染性单核细胞增多症、重度牙周炎各 15］；健康对照 100 份［流行区/非流行区各 50］）；(3) 各分期临床特征化样本灵敏度；(4) 流行区/非流行区无症状人群特异度。

**BioPlex 2200 Lyme Total（K183446）**
- 前瞻 792 份（3 中心）vs 市售 IgM/IgG 免疫分析：一线 PPA 70.1%（68/97），NPA 95.1%（661/695）；二线 WB（IgG+IgM 合并）PPA 97.9%（46/47）。
- 105 份分期样本：Acute（<3 月）69.4%（50/72）vs 比较法 58.3%；Convalescent（<12 月）61.5%（16/26）；Late（>12 月）85.7%（6/7）；总 68.6%。
- CDC panel 280：Acute 84.6%（33/39）、Convalescent 93.5%、Late 100%（20/20）；look-alike 96.7%（1 反应+2 equivocal/90）；健康对照 97.0%。细分：早期 EM 急性 80.0%（vs 比较法 63.3%）、恢复期 93.3%、心脏/神经播散 100%、晚期关节炎/神经 100%、梅毒 15/15 阴性、传单 2 份 equivocal。
- 健康人 836 份（流行/非流行各 ~420）：NPA 96.4%/96.4%。

**ZEUS MTTT（K191240 IgM 路径；K190907 IgG/IgM 路径）**
- 回顾队列 356 份（CDC 280 + 46 Stage 2 + 30 Stage 3；共 166 例 LD：Stage 1/2/3 = 60/56/50；90 份其他疾病；100 份健康）。VlsE1/pepC10 一线 160 阳性 + 6 equivocal。
  - IgM 路径（K191240）灵敏度 STTT vs MTTT：Stage I 46.7% → **76.7%**；Stage II 50.0% → 75.0%；Stage III 16.0% → 72.0%；健康对照特异度 100%/100%；疾病对照 100%/97.8%。
  - IgG/IgM 路径（K190907）：Stage I 63.3% → **78.3%**；Stage II 60.7% → 66.1%；Stage III 100%/100%；健康 100%/100%；疾病对照 100%/97.8%。
- 前瞻队列 2,932 份（MA 900、WI 990、MN 1042）：一线 363 阳 + 58 equivocal。
  - IgM：MTTT vs WB-STTT PPA 96.2%（101/105；90.6–98.5），NPA 95.5%（2701/2827；94.7–96.2）；126 份 MTTT+/STTT− 中 5 份临床符合 Stage 1，91 份无临床资料。
  - IgG/IgM：PPA 93.3%（167/179；88.6–96.1），NPA 97.7%（2690/2753；97.1–98.2）；63 份 MTTT+/STTT− 中 4 份为确诊 Lyme（3 Stage 1、1 晚期）。

**LIAISON Lyme IgG（K202574）/ IgM（K202573）**
- 前瞻 2,621 份（14 州 5 区域，3 实验室）：
  - IgG 一线 vs ZEUS IgG ELISA：PPA 55.1%（166/301）、NPA 96.5%（2210/2320）；两者 WB IgG 阳性数相同（各 109）；predicate+/LIAISON− 的 135 份中 123（91.1%）WB 阴性。STTT 算法级 PPA 89.0%（97/109；81.7–93.6）、NPA 99.5%（2500/2512）。
  - IgG 作 MTTT 二线（Lyme Total Antibody Plus 一线阳/equivocal 225 份）vs WB-STTT：PPA 96.9%（93/96；91.2–98.9），NPA **30.2%**（39/129；23.0–38.6）。
  - IgM（K202573）：一线 PPA 56.5%（190/336）、NPA 96.5%（2206/2285）；STTT 级 PPA 91.6%（109/119）、NPA 99.4%（2486/2502）；MTTT vs STTT PPA 93.0%（93/100）、NPA 57.6%（72/125）。
- CDC panel（IgG 一线）：Acute 74.4%（vs ZEUS 61.5%）、Convalescent 93.5%、Late 100%、健康 97%（3/100 假阳）、疾病对照 90%（9 假阳：梅毒 4、传单 2、牙周炎 1、RA 1、纤维肌痛 1）。
- CDC panel MTTT vs WB-STTT（IgG）：Stage I **30% → 80%**；Stage II 50% → 90%；Stage III 100%/100%；健康 100%/100%；疾病对照 100%/97.8%。IgM（K202573）：Stage I 50% → 73.3%；Stage II 90%/90%；Stage III 35% → 45%；对照均 100%。
- 健康人 300 份：流行区阳性 6.7%（10/150），非流行区 1.3%（2/150）。

**GSD VlsE-OspC IgG/IgM ELISA（K203289）**
- 前瞻 481 份（3 中心）vs predicate 全菌 IgG/IgM ELISA：一线 PPA 90.7%（49/54）、NPA 96.7%（413/427）；STTT 算法级 PPA 100%（36/36）、NPA 99.6%（443/445）。
- MTTT（VlsE-OspC → IgG/IgM ELISA）vs WB-STTT：全样本 PPA 100%（36/36）、NPA 97.1%（432/445）；仅一线阳性 68 份内 NPA 40.0%（12/30）。IgG 路径 PPA 100%（23/23）、NPA 96.7%；IgM 路径 PPA 94.7%（18/19）、NPA 95.2%。
- 125 份分期样本（一线单独）：Early 62.9%（39/62）、Disseminated 100%（22/22）、Late 97.6%（40/41）；MTTT IgG 路径 vs WB-STTT IgG：Early 48.4% vs **8.1%**、Disseminated 81.8% vs 27.3%、Late 97.6% vs 95.1%；MTTT IgM 路径 vs WB-STTT IgM：Early 62.9% vs 58.1%、Late 82.9% vs 24.4%。
- CDC panel 一线：Early 81.7%（49/60；69.6–90.5）、Cardiac 66.7%（2/3）、Neurological 100%（7/7）、Late 100%（20/20）、健康 97%、疾病对照 95.6%（vs predicate 82.2%）。MTTT IgG/IgM：Early 76.7% vs STTT 61.7%；健康 100%；疾病对照 98.9%。
- 无症状人群 238 份：流行区特异度 97.0%（4/132 阳/equivocal），非流行区 99.1%（1/106）。

**Sofia 2 Lyme FIA（K173496）**
- 前瞻 327 例症状患者（11 流行区 POC 站点，CLIA-waived 操作者），324 可评估，配对指尖全血 vs VIDAS Lyme IgM/IgG（predicate）：IgM 一线 PPA 82.4%（75/91）、NPA 79.8%（186/233）；IgG PPA 88.9%（48/54）、NPA 85.9%（232/270）；二线 WB PPA：IgM 94.1%（48/51）、IgG 95.7%（22/23）。
- 95 份临床/培养确诊样本：IgM 总灵敏度 64.2%（Acute <1 月 EM 60.9%；Late 66.7%），IgG 总 80.0%（Acute 78.1%；Late 100%），均高于 predicate（58.9%/49.5%）。
- CDC panel 280：IgM 在阴性对照仅 82.6% 一致（33/190 假阳），早期 EM 81.7%，晚期 73.3%；IgG 阴性对照 86.8%（25/190 假阳），早期 81.7%，晚期 100%。
- 健康人 200 份：IgM 阴性率 89.5%（流行区 86.0%），IgG 96.5%。

**iDart Lyme IgG ImmunoBlot（K233367）**
- 768 份（Bay Area Lyme Foundation 前瞻性库存 290；IGeneX 送检 248 + 230）vs STTT（FDA 清关 EIA + immunoblot）：Cohort 1 PPA 95.00%（19/20）、NPA 86.67%（234/270）；Cohort 2 PPA 95.00%（114/120）、NPA 90.63%；Cohort 3 PPA 90.91%（10/11）、NPA 96.80%（212/219）。
- CDC panel：Stage I 灵敏度 **58.33%（35/60）vs STTT 30.00%**；Stage II 90%（9/10）/90%；Stage III 100%（20/20）；健康 100 份及疾病对照 90 份均 0 阳性。
- 健康人：流行区 313 份阳性 0.64%（2/313）；非流行区 112 份 0%。

### 9. 厂家间差异与要点
1. **STTT → MTTT 的转变**：2019 年 ZEUS 三项清关首次将 "两个 EIA 串联" 写入 intended use，并以回顾（CDC 盘）+ 前瞻（~2,900 份）两队列证明 MTTT 对 Stage I 灵敏度较 WB-STTT 提升 15–30 个百分点而对照特异度维持 ≥97.8%。DiaSorin（2021）、GSD（2021）、ID-FISH（2024）沿用同一证据框架。
2. **MTTT 二线的 NPA 天然偏低**（DiaSorin IgG 30.2%、GSD 40.0%、GSD IgM 53.2%）：因二线 EIA 比 WB 更敏感，FDA 决策摘要接受此现象并明确 "The above performance table artificially inflates the negative percent agreement…"，评价重点转向分期灵敏度与健康/疾病对照特异度。
3. **cutoff 建立方法多样**：Bio-Rad 用 1,372 份正常血清 98 百分位；GSD 用 238 份正常血清 均值+3SD 再 ROC；DiaSorin 以特征化样本平衡灵敏/特异；Quidel 以全血-血清 PPA/NPA 等效为目标。灰区宽度普遍 ±10%（0.9–1.1 Index 或 9–11 Units），Bio-Rad 为 0.9–1.0 AI。
4. **抗原组合决定分期表现**：含 OspC 的一线（Bio-Rad、GSD、DiaSorin IgM）在 CDC 早期 EM 组灵敏度 80–85%；纯 VlsE IgG（DiaSorin IgG）急性期 74.4%；全菌 predicate 61.5–63.3%。
5. **交叉反应差异**：微珠/CLIA 平台对梅毒、传单假阳性 1–4/15；POC 侧向层析（Sofia）IgM 在 SLE、梅毒、ANA、EBV 假阳性率明显更高（20–40%），CDC 盘阴性对照假阳性 17.4%。
6. **免疫印迹判读现代化**：iDart 用 Protein L 对照带（C2）作强度校准并要求 LSA 带 + 两组特异带，与传统 5/10 条带规则不同，Stage I 灵敏度 58% vs STTT 30%。

### 10. 来源
- K183446：https://www.accessdata.fda.gov/cdrh_docs/reviews/K183446.pdf
- K191240：https://www.accessdata.fda.gov/cdrh_docs/reviews/K191240.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K191240.pdf
- K190907：https://www.accessdata.fda.gov/cdrh_docs/reviews/K190907.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K190907.pdf
- K202574：https://www.accessdata.fda.gov/cdrh_docs/reviews/K202574.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K202574.pdf
- K202573：https://www.accessdata.fda.gov/cdrh_docs/reviews/K202573.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K202573.pdf
- K203289：https://www.accessdata.fda.gov/cdrh_docs/reviews/K203289.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K203289.pdf
- K173496：https://www.accessdata.fda.gov/cdrh_docs/reviews/K173496.pdf
- K233367：https://www.accessdata.fda.gov/cdrh_docs/reviews/K233367.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K233367.pdf

---

## 梅毒螺旋体抗体（*Treponema pallidum* treponemal antibodies）

### 1. 法规定位
- Product code **LIP**（Enzyme linked immunoabsorption assay, Treponema pallidum）；非密螺旋体（RPR）抗原为 GMQ（866.3820）；校准品 JIT、控制品 JJX。
- 21 CFR **866.3830** Treponema pallidum treponemal test reagents；**Class II**；Panel：Microbiology (83 / MI)。
- openFDA 累计清关数：**LIP total = 35**。
- 注：用户提及的 "BioPlex 2200 Syphilis Total & IgM" 未见于 LIP 清单；Bio-Rad 在 LIP 下的清关为 K063866（Syphilis IgG，2007）、K120439（EBV IgG + Syphilis IgG，2012）与 **K170413（Syphilis Total & RPR，2017）**，本节采用 K170413。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K153730 | 2016-06-15 | Abbott / ARCHITECT Syphilis TP | 两步 CMIA；重组 TpN15、TpN17、TpN47 包被微粒；acridinium 标记抗人 IgG+IgM | 血清、K2/K3-EDTA、Li/Na-heparin、Li-heparin PST、SST | 总抗体（IgG+IgM）；**initial diagnostic test 或与非密螺旋体试验联用**；不用于血液/血浆/组织供者筛查 | 决策摘要 + 510(k) summary |
| K170413 | 2017-05-11 | Bio-Rad / BioPlex 2200 Syphilis Total & RPR | 多重微珠：rTP47/rTP17 融合蛋白微珠 + cardiolipin 微珠；抗人 IgG+IgM PE 标记 | 血清、K2/K3-EDTA、Li/Na-heparin 血浆 | 总抗体 + RPR 定性/滴度；"may be used to supplement a previously determined reactive treponemal or non-treponemal test" | 决策摘要 + 510(k) summary |
| K112343 | 2012-01-20 | Siemens / ADVIA Centaur Syphilis | 直接双抗原夹心化学发光；重组 TpN17+TpN15 生物素化及 acridinium ester 标记 | 血清、EDTA/肝素/枸橼酸血浆 | aid in the diagnosis of syphilis；不用于 blood/tissue donor screening | 决策摘要 + 510(k) summary |
| K061247 | 2006-07-31 | DiaSorin / LIAISON Treponema | 一步双抗原夹心 CLIA；重组 Tp17（DNA-Tp17）磁微粒与 isoluminol-Tp17 示踪 | 血清 | 总抗体（IgG/IgM）；"in conjunction with non treponemal laboratory tests and clinical findings" | 决策摘要 |
| K241427 | 2024-09-06 | Beckman Coulter / Access Syphilis | 两步酶化学发光双抗原夹心；重组 Tp17、Tp47 顺磁微粒 + 生物素化 Tp17/Tp47 + ALP 结合物 | 血清、SST、Li-heparin、Li-heparin gel、K2/K3-EDTA、枸橼酸钠、CPD、ACD、CPDA | aid in the diagnosis of syphilis 或与非密螺旋体试验联用；不用于 donor screening | 决策摘要 + 510(k) summary |
| K160910（K211302 为生物素改良） | 2016-07-28 / 2021-07-20 | Roche / Elecsys Syphilis | 双抗原夹心 ECLIA；重组 TpN15、TpN17、TpN47（单体与多聚体） | 血清、Li/Na-heparin、K2/K3-EDTA、CPDA、枸橼酸钠血浆 | 总抗体；aid in the diagnosis in conjunction with clinical signs；不用于 blood/tissue donor | 决策摘要 + 510(k) summary |

### 3. 预期用途与声明类型
- 全部为 **Rx 定性、总抗体（IgG+IgM）或 "antibodies to T. pallidum"**（ADVIA Centaur 措辞因 "administrative error" 删去 IgG 字样）。
- **算法定位（申报文件）**：
  - 2006 年 LIAISON（K061247）：仅 "in conjunction with nontreponemal laboratory tests"，并附警示 "A positive result is not useful for establishing a diagnosis of syphilis. In most situations, such a result may reflect prior treated infection"；equivocal/阳性须补做定量 RPR/VDRL。FDA 在该文件中同时评价了 "clinical laboratory screen (treponemal test followed by a non treponemal test)" 与 "diagnostic confirmatory test（RPR/VDRL 阳性后）" 两种用法。
  - 2016 年起 ARCHITECT（K153730）、Access（K241427）明文 "**intended to be used as an initial diagnostic test** or in conjunction with a nontreponemal laboratory test"，即支持反向算法（treponemal 先行）。Elecsys 要求初始反应性样本双份复测，"Repeatedly reactive samples must be confirmed according to recommended confirmatory algorithms"。
- 共同限定："not intended for use in screening blood, plasma, or tissue donors"；结果须结合其他 treponemal/nontreponemal 结果、临床与病史 "to produce a diagnosis of syphilis by disease stage"。
- 背景（非申报文件）：CDC 传统算法为 RPR/VDRL 筛查 → treponemal 确认；反向算法（reverse sequence）为 treponemal EIA/CIA 筛查 → RPR 定量 → 不一致时 TP-PA 裁定（CDC MMWR 2011）。申报文件中 "2 out of 3（TP-CLIA、RPR、TP-PA）" 复合比较法即模拟该裁定逻辑。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 单位 | 判定 | Equivocal | cutoff 建立（文件所载） |
|---|---|---|---|---|
| ARCHITECT Syphilis TP（K153730） | S/CO；Cutoff RLU = 校准品均值 RLU × **0.2** | <1.00 Nonreactive；≥1.00 Reactive | **无**（文件明确 "without the need to assign an equivocal zone"） | **6,845 份**样本 ROC（409 抗体阳性；6,436 "normal" 供者+诊断人群），以与 **TP-PA** 最优一致为目标选定 cutoff multiplier 0.2，再在临床研究中验证 |
| BioPlex 2200 Syphilis Total（K170413） | AI | ≤0.8 Nonreactive；≥1.1 Reactive | **0.9–1.0 AI**（±10%，"based on the precision of the assay"）；RPR 无灰区（<1.0/≥1.0） | 可行性阶段用健康人、送检者、确诊者样本，以 predicate 结果为标准做 ROC，调整校准值使 cutoff = 1.0 AI |
| ADVIA Centaur Syphilis（K112343） | Index（两点校准） | <0.90 Nonreactive；≥1.10 Reactive | **≥0.90–<1.10 Equivocal**（±10%，"to encourage additional testing"） | 阴/阳样本评估显示 cutoff 1 清楚分离；临床试验中用 ROC（806 健康 + 285 确诊）验证，AUC 与 IMMULITE 2000 无显著差异（p=0.5282），cutoff 1.0 时 PPA 97.9%、NPA 99.4% |
| LIAISON Treponema（K061247） | Index（量程 0.1–70） | <0.9 Negative；≥1.1 Positive | **0.9–1.1 Equivocal**（复测；再 equivocal 则 ≥1 周后复采） | 欧洲 **3,625 份**（1,000 常规、131 特征化梅毒、2,494 供者）累积频率分布/ROC；美国验证 97% 阴性 <0.5、97.9% 阳性 >1.4 |
| Access Syphilis（K241427） | S/CO | <1.00 Nonreactive；≥1.00 Reactive | 无 | **301 阳性 + 3,047 阴性** ROC（DxI 9000 与 Access 2），以样本信号/阳性校准品信号归一化比值取最佳折中点 |
| Elecsys Syphilis（K160910/K211302） | COI | <1.00 Nonreactive；≥1.00 Reactive；初始反应须双份复测 | 无 | 原型批测定天然血清设初步 cutoff → 与市售检测比较 → 欧美多中心临床最终验证 |

### 5. 生物学 / 免疫学依据
- **抗原（申报文件）**：均为大肠杆菌重组 T. pallidum 脂蛋白：Tp17（TpN17）为所有产品共有；Tp47、Tp15 为 Abbott/Roche（三抗原）、Beckman（Tp17+Tp47）、Siemens（Tp17+Tp15）所加；Bio-Rad 用 rTP47/rTP17 融合蛋白；DiaSorin 仅 Tp17。predicate CAPTIA Syphilis-G 为 Nichols 株全菌超声抗原。
- **格式**：双抗原夹心（Siemens、DiaSorin、Roche、Beckman）可同时捕获 IgG/IgM；Abbott 用抗人 IgG+IgM 混合结合物；Bio-Rad 用抗 IgG+抗 IgM PE 结合物。
- **IgM 检测证明（申报文件）**：Abbott 用 Protein G 去 IgG 后以 IgM 特异结合物证明 IgM 信号增加 >2 倍；Bio-Rad 用 DTT 去 IgM + 9 成员 **血清转换盘**（Syphilis Total 在第 31 天 AI 1.3 阳性，早于 IgG 检测与 TP-PA［第 45 天］）。
- **Hook effect**：Abbott（≥20 S/CO 稀释）、Beckman、DiaSorin（Index 73–132 无 hook）均评估；Bio-Rad 称两步法带洗涤不适用。
- **交叉反应（申报文件）**：多数反应性标本经 TP-PA/RPR/FTA-ABS 证实为真阳性（Siemens：HAV 10/20、HSV 5/10、HAMA 2/10 等均 predicate 亦阳性；Roche 9/266 均经 Bio-Rad EIA 证实）。未被证实的假阳性零星见于 CMV IgG、HIV、hyper-IgG、E. coli（Abbott 各 1 份）、HTLV-II、myeloma、HBc IgM（Beckman 各 1 份）。DiaSorin 2006 年文件注明 HAV、yaws、pinta、leptospirosis 未评估并写入限制。Beckman 报告 γ-球蛋白 ≥60 g/L 有干扰（≤47.5 g/L 无）。
- 背景（非申报文件）：Tp15/Tp17/Tp47 为免疫优势外膜脂蛋白，treponemal 抗体一旦产生通常终身存在（治疗后不转阴），故治疗史患者阳性不代表活动感染——与 LIAISON 警示语一致；yaws/pinta 等地方性密螺旋体病抗体无法区分。

### 6. 样本类型与样本要求
- ARCHITECT：7 种管型（血清、SST、K2/K3-EDTA、Li/Na-heparin、Li-heparin PST），28 组配对，阴性差异 ≤0.20 S/CO、阳性 %差 ≥-20%；稳定性 室温 ≤72 h、2–8 °C ≤7 天、≤-10 °C ≤30 天、6 次冻融；仪器上 3 h。
- Access：10 种管型（含枸橼酸钠、CPD、ACD、CPDA），65 组配对；2–8 °C 7 天、20–25 °C 72 h、-20 °C 30 天、5 次冻融。
- BioPlex：血清、K2/K3-EDTA、Li/Na-heparin（82–113 组）；室温 3 天、-20 °C 6 个月、5 次冻融。
- ADVIA Centaur：血清（玻璃/SST）、K2-EDTA、Li/Na-heparin、枸橼酸钠；**ACD 管不推荐**（阴性样本 Index 轻微升高）。
- LIAISON（2006）：仅血清；2–8 °C ≤7 天，避免反复冻融。
- Elecsys：血清、Li/Na-heparin、K2/K3-EDTA、CPDA、枸橼酸钠。

### 7. 分析性能验证所依据的标准

| K 号 | 文件所列标准/指南 | 对应项目 |
|---|---|---|
| K153730 | "Not applicable" | 22 天 3 批 2 仪器精密度；3 中心 3 批重复性；7 种内源干扰；416 份 34 类交叉反应；hook；携带污染；IgG/IgM 反应性；样本稳定性 |
| K170413 | CLSI EP05-A3；EP07-A2；EP09-A3（仅基质比较）；EP12-A2；EP15-A3；EP25-A | 20 天精密度；3 中心 120 重复；RPR 滴度重复性；9 种干扰；332 份 28 类交叉反应；血清转换盘；稳定性 |
| K112343 | CLSI EP05-A2；CLSI GP10-A（ROC 评价临床准确度）；CLSI EP12-A2；FDA *Assayed and Unassayed Quality Control Material* 指南（2007-06-07） | 20 天精密度；3 中心 2 批 10 天；ROC；211 份特殊人群 + 265 份 20 类交叉反应；8 种干扰 |
| K061247 | "Not Applicable" | 3 中心 3 批 5 天重复性；hook；244 份 21 类交叉反应；3 种干扰；冻融/2–8 °C 稳定性；携带污染 |
| K241427 | CLSI EP05-A3；EP07 3rd ed.；EP12-A2；**EP24-A2**（ROC）；EP25-A；EP37；GP44-A4；GP41 7th ed.；ISO 7000、ISO 7010、ISO 15223-1、ISO 13485:2016 | 20 天 3 批双平台精密度；3 中心；6 种内源 + 7 种药物干扰（含生物素 0.351 mg/dL）；424 份 39 类交叉反应；hook；试剂/样本稳定性；标签符号 |
| K160910 / K211302 | K160910 "Not applicable"（正文引 CLSI EP5-A2 重复性）；K211302 引用 *Deciding When to Submit a 510(k) for a Change* 与 *Refuse to Accept* 指南 | 21 天精密度；生物素 ≤1200 ng/mL 偏差 ≤10%；266 份交叉反应 |

### 8. 临床验证设计与结果

**统一设计（2012 年后）**：前瞻性 intended-use 人群（常规送检、孕妇、HIV 阳性）+ 回顾性预选阳性 + 医学确诊分期（primary/secondary/latent，treated/untreated）+ 表观健康人（含儿科）；**参考方法为 TP-CLIA（predicate）+ RPR + TP-PA 三项 2/3 复合算法**（Abbott、Bio-Rad、Beckman、Roche），或 predicate 单一比较法（Siemens）、predicate + TP-PA/RPR 裁定（DiaSorin）。

**ARCHITECT Syphilis TP（K153730）**——2,220 份，2015 年 7–12 月，8 个美国站点。
- 前瞻 1,145 份（常规 442、孕妇 304、HIV+ 399；儿科 136）：PPA 96.2%（153/159；92.0–98.3），NPA 99.0%（976/986；98.1–99.4）。分层：常规 PPA 97.3%/NPA 99.5%；孕妇 NPA 99.7%（303/304，各孕期 99.2–100%）；HIV+ PPA 95.9%（117/122）/NPA 97.5%（270/277）。6 份假阴性均为 TP-CLIA 阴性但 TP-PA+RPR 阳性。
- 回顾预选 406 份：PPA 98.9%（376/380），NPA 92.3%（24/26）；20 份阳性孕妇 100%。
- 医学确诊 179 份：Primary treated 33/44 反应（11 份 nonreactive，其中 9 份 TP-CLIA 亦阴性）、Primary untreated 25/25、Secondary 56/56、Latent 54/54。
- 表观健康 480 份：反应率 3.1%（成人 15/367=4.1%，14 份 TP-CLIA 亦阳性；儿科 0/113）。前瞻人群阳性率 14.2%。

**BioPlex 2200 Syphilis Total & RPR（K170413）**——2,008 份，2016 年 11–12 月，3 站点（1 内部）。
- 前瞻 1,001 份（送检 401、孕妇 295、HIV+ 305）：Syphilis Total PPA 92.45%（147/159；87.27–95.63），NPA 97.86%（824/842）；RPR vs predicate RPR PPA 81.52%（75/92）、NPA 96.48%。
- 回顾 544 份：Total PPA 99.59%（486/488）、NPA 100%（56/56）；RPR PPA 98.14%、NPA 80.70%。
- 孕妇 372 份：Total PPA 100%（32/32）、NPA 98.83%（338/340）；HIV+ 362 份：PPA 93.33%（140/150）、NPA 93.87%（199/212）。
- 确诊 156 份：Primary untreated 96.2%（25/26）、treated 86.2%（25/29）；Secondary、Latent 均 100%；总 96.8% vs 复合比较法 97.4%。RPR 在治疗后分期敏感度明显下降（Primary treated 65.5%、Latent treated 66.7%）。
- 表观健康 301 份：NPA 99.0%。

**ADVIA Centaur Syphilis（K112343）**——2,108 份，3 站点，比较法为市售梅毒检测（predicate IMMULITE 2000）。
- 总体 PPA 97.9%（700/715；96.6–98.8），NPA 99.4%（1382/1391；98.8–99.7）。
- 表观健康 806 份（孕妇 332、儿科 75、其他 399）：反应 0.6%，NPA 98.8%（791/801）。
- 预期阳性 561 份（TPPA 反应 276、确诊 285）：PPA 99.4%（535/538）；确诊组 21/285 nonreactive（比较法亦阴性）。
- Intended-use 741 份（常规 + HIV+）：PPA 98.2%（160/163），NPA 98.4%（568/577）。

**LIAISON Treponema（K061247）**——3 站点冻存样本，predicate Trinity CAPTIA Syphilis-G。
- 美国确诊 51 份：PPA 97.9%（47/48）；欧洲 127 份 99.2%（118/119）。
- 送检 999 份（NY + 东南部）：PPA 55%（22/40）、NPA 98.9%（909/919）；12 份经 RPR/TP-PA 仍不能裁定。
- HIV+ 200 份：PPA 75.8%（69/91）、NPA 96.2%（100/104）；孕妇 200 份：PPA 100%（4/4）、NPA 100%（192/192）；表观健康 992 份：PPA 62.7%（54/86）、NPA 99.3%（881/887）。
- 作确认试验（RPR/VDRL 阳性 204 份）：PPA 99.5%（200/201）。
- 注：该文件的低 PPA 反映 predicate（全菌 IgG EIA）与 Tp17 重组抗原检测间不一致，非灵敏度不足；FDA 因此在 2016 年后转向三项复合比较法。

**Access Syphilis（K241427）**——1,910 份，2023 年 6–9 月，3 站点，双平台（Access 2 / DxI 9000）。
- 前瞻 1,104 份（常规 399、孕妇 405、HIV+ 300）：PPA **100%**（184/184；98.0–100），NPA 96.7%（890/920；95.4–97.7）。分层：常规 NPA 99.7%；孕妇 PPA 100%（6/6）/NPA 99.7%（398/399）；HIV+ NPA **84.6%**（154/182），28 份比较法阴性者 25（89.3%）经另一 TP-CLIA 亦反应。
- 回顾 402 + 高危 STI 50 份：PPA 100%（378/378；20/20）；高危组 NPA 80.0%（24/30）。
- 确诊 150 份：Primary treated 96.3%（26/27）、untreated 100%（49/49）；Secondary 100%；Latent treated 100%（25/25）、**untreated 53.8%（7/13，6 份另一 TP-CLIA 亦阴性）**；总 95.3%。
- 表观健康 204 份：成人反应 9.3%（18/194），<21 岁 0/10。前瞻阳性率 19.4%（男 38.8%、女 8.4%）。

**Elecsys Syphilis（K160910；K211302 引用）**——前瞻常规、301 孕妇、457 HIV+：PPA 100%（228/228；98.4–100），NPA 99.2%（2038/2054；98.7–99.6）；HIV+ NPA 95.6%（282/295）；孕妇 NPA 100%（301/301）。回顾 169 份：PPA 98.7%（155/157）、NPA 100%（12/12）。确诊：Primary treated 16/29 反应（13 份 nonreactive，其中 12 份复合算法亦阴性），Primary untreated 25/25，Secondary treated 24/25、untreated 25/25，Latent 50/50。表观健康 209 份：反应 9.6%（20/209；女 11.3%、男 8.5%）。K211302 仅做 232 份样本与原版 100% 定性一致。

### 9. 厂家间差异与要点
1. **灰区**：DiaSorin（0.9–1.1）、Siemens（0.9–<1.1）、Bio-Rad（0.9–1.0）保留 ±10% 灰区；Abbott、Beckman、Roche 以大样本 ROC（3,348–6,845 份）取消灰区，Abbott 明确以 TP-PA 一致性为优化目标。
2. **参考方法**：从 predicate 单一比较（2006、2012）演进为 **TP-CLIA + RPR + TP-PA 2/3 复合算法**（2016 起 Abbott、Bio-Rad、Roche、Beckman 一致），FTA-ABS 仅在 Abbott 交叉反应研究中作最后裁定。
3. **抗原数量与早期灵敏度**：三抗原（Tp15/17/47）产品在 untreated primary 100%；治疗后 primary 反应率 75–96%（Abbott 33/44、Roche 16/29、Bio-Rad 25/29、Beckman 26/27），反映治疗后抗体水平下降（背景推断）。Beckman 在 untreated latent 仅 53.8%，为文件中最低分期数据。
4. **HIV+ 人群 NPA 偏低**（Beckman 84.6%、Bio-Rad 93.9%、Roche 95.6%、Abbott 97.5%）：文件解释多为比较法阴性但另一 TP-CLIA 阳性，提示复合算法在高流行人群中低估真阳性。
5. **Bio-Rad 独有 RPR 微珠**：以自动化非密螺旋体检测（AI）+ 在线稀释滴度（1:4–1:64）与 treponemal 同管完成，可同时输出反向算法两步结果；RPR 对 SLE 交叉（11/12 NPA）写入限制。
6. **管型**：从仅血清（DiaSorin 2006）扩展到 10 种含枸橼酸/CPD/ACD 管型（Beckman 2024）；Siemens 排除 ACD。

### 10. 来源
- K153730：https://www.accessdata.fda.gov/cdrh_docs/reviews/K153730.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153730.pdf
- K170413：https://www.accessdata.fda.gov/cdrh_docs/reviews/K170413.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K170413.pdf
- K112343：https://www.accessdata.fda.gov/cdrh_docs/reviews/K112343.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K112343.pdf
- K061247：https://www.accessdata.fda.gov/cdrh_docs/reviews/K061247.pdf
- K241427：https://www.accessdata.fda.gov/cdrh_docs/reviews/K241427.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K241427.pdf
- K160910：https://www.accessdata.fda.gov/cdrh_docs/reviews/K160910.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K160910.pdf
- K211302：https://www.accessdata.fda.gov/cdrh_docs/reviews/K211302.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K211302.pdf

---

## 巨细胞病毒 IgG / IgM（Cytomegalovirus, CMV IgG / IgM）

### 1. 法规定位
- Product code **LFZ**（Cytomegalovirus serological reagents）；部分文件附 JJE（分析仪）、JJX（控制品）。
- 21 CFR **866.3175** Cytomegalovirus serological reagents；**Class II**；Panel：Microbiology (83 / MI)。
- openFDA 累计清关数：**LFZ total = 58**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K220949 | 2022-10-27 | Abbott / ARCHITECT CMV IgG | 两步 CMIA；**CMV 病毒裂解物（AD169 株）**包被微粒；acridinium 标记鼠抗人 IgG；6 点校准 0–250 AU/mL | 血清、SST、Li-heparin、Li-heparin PST、K3-EDTA | aid in the diagnosis + **determination of serological status**，含育龄妇女；未清关用于血液/血浆/组织供者筛查 | 决策摘要 + 510(k) summary |
| K181213 | 2018-07-30 | Siemens / ADVIA Centaur CMV IgG | 两步间接化学发光；生物素化 CMV 病毒裂解物抗原 + acridinium ester 抗人 IgG | 儿科及成人血清、K2-EDTA、Li-heparin 血浆 | serological status + aid in diagnosis，含孕妇；不用于 blood/tissue donor | 决策摘要 + 510(k) summary |
| K131605（K220911 为生物素改良） | 2014-02-28 / 2022-10-12 | Roche / Elecsys CMV IgG | 双抗原夹心 ECLIA；**重组** CMV 特异抗原（E. coli）生物素化/钌标记 | 血清、Li-heparin、K2/K3-EDTA | determination of serological status，含孕妇；免疫抑制者未评估；不用于新生儿筛查、POC、献血者 | 决策摘要 + 510(k) summary |
| K142133（K163569 为 cobas e 801 平台扩展） | 2014-10-24 / 2017-03-17 | Roche / Elecsys CMV IgM | **µ-capture** ECLIA：生物素化单抗抗人 IgM 捕获 + 钌标记重组 CMV 抗原 | 血清、Li-heparin、K2/K3-EDTA | aid in the diagnosis of **recent or current** infection，含孕妇；免疫抑制未评估；不用于新生儿筛查、POC、献血者 | 决策摘要 + 510(k) summary |
| K162969（原始 K040290，2005） | 2017-01-06 | DiaSorin / LIAISON CMV IgG + Serum Control Set | 间接 CLIA（Special 510(k)：仅控制品基质变更） | 血清 | aid in the determination of serological status to CMV | Special 510(k) 决策摘要（无性能数据）；K040290 文件未能获取 |

### 3. 预期用途与声明类型
- **IgG**：双重声明——(a) "aid in the determination of serological status"（免疫状态/既往暴露）；(b) "aid in the diagnosis of CMV infection"（Abbott、Siemens）。Roche IgG 仅声明 serological status。目标人群 "individuals for whom a CMV IgG test was ordered, including pregnant women / women of child-bearing age"（Abbott）。
- **IgM**：仅 "aid in the diagnosis of recent or current CMV infection"（Roche），不作免疫状态。
- 共同限定：**"not intended/cleared for use in screening blood, plasma, or tissue donors"**（全部）；Roche 另加免疫抑制者未评估、不用于新生儿筛查、不用于 POC；Siemens 儿科（2–21 岁）人群纳入 intended use 并有专门亚组数据；移植受者仅 Siemens 做了亚组（术前/术后）。
- 均为 Rx 用；定性（Abbott 以 AU/mL 报告但声明为定性，reportable range 0–250 AU/mL）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 单位 | 判定 | Equivocal | cutoff 建立（文件所载） |
|---|---|---|---|---|
| ARCHITECT CMV IgG（K220949） | AU/mL | <6.0 Nonreactive；≥15.0 Reactive | **6.0–<15.0 Grayzone/Equivocal**（建议 CMV IgM 复核或 ~2 周后复采） | 用 OUS 验证研究 530 份 × 3 批（1,590 结果）以 **AxSYM CMV IgG 为真值**，确认沿用 AxSYM cutoff 并设下限 6/上限 15 AU/mL 以保持灵敏度；美国临床研究验证 |
| ADVIA Centaur CMV IgG（K181213） | Index | <1.00 Nonreactive；≥1.00 Reactive | **无** | 389 份剩余临床样本（196 阴、186 阳、余 equivocal）与比较法比较，cutoff 设在 **PPA ≥95% 且 NPA ≥95%** 处（1.00 Index） |
| Elecsys CMV IgG（K131605/K220911） | COI | <0.5 Non-reactive；≥1.0 Reactive | **0.5–<1.0 Borderline**（复测、2 周内复采或替代方法） | 931 份样本 **ROC** 优化灵敏/特异，临床研究验证；hook ≤2,500 COI 无 |
| Elecsys CMV IgM（K142133/K163569） | COI | <0.7 Non-reactive；≥1.0 Reactive | **0.7–<1.0 Indeterminate (Border)**（复测；仍不定则 2 周内复采） | 用多种市售 IgG/IgM 检测特征化样本；**931 份低流行队列**评估特异度，**152 份原发感染不同阶段**样本优化灵敏度；外部临床验证。LoB 0.243 COI、LoD 0.267 COI（EP17-A2）；hook ≤18.7–32.9 COI 无 |
| LIAISON CMV IgG（K162969） | — | 文件（Special 510(k)）未载明 | 文件未载明 | 引用 K040290（未获取） |

### 5. 生物学 / 免疫学依据
- **抗原（申报文件）**：Abbott、Siemens 用 **AD169 株病毒裂解物**（多抗原天然混合）；Roche IgG/IgM 用 **大肠杆菌重组 CMV 特异抗原**（文件未指明具体蛋白）。predicate VIDAS CMV IgG（K920661）亦为 AD169 抗原（K131605 表格）。
- **IgM 捕获法 vs 间接法（申报文件）**：Elecsys CMV IgM 为 µ-capture，文件说明 "designed to specifically bind the anti-CMV IgM antibodies due to the µ-capture design and antigens minimizing interference from CMV IgG antibodies"；其 predicate Diamedix Is-CMV IgM Capture 亦为捕获 EIA。IgM 特异性验证：11 份 IgM 阳性样本 DTT 处理后 10 份转阴、1 份 COI 下降 82%。RF 干扰：K142133 声明 RF <2000 IU/mL 无干扰；K163569 用高 RF 血清 1:1 稀释实验声明 RF <899 IU/mL。本组文件中未包含间接法 IgM 产品（LIAISON CMV IgM 属 K040290，未获取）。背景（非申报文件）：间接法 IgM 易受 RF（IgM 类抗 IgG）与高浓度特异 IgG 竞争影响，故需 IgG 吸附；µ-capture 先捕获总 IgM 再加抗原，可避免上述干扰。
- **IgG 动力学（申报文件）**：Abbott 结果解释 "Reactivity for anti-CMV IgG indicates past or acute infection"；Nonreactive "presumed to be not infected with CMV and susceptible to primary infection"。Elecsys IgM 前瞻研究显示送检人群 IgM 阳性率仅 1.67%（10/599）。
- **交叉反应（申报文件）**：由于 CMV IgG 人群流行率高（Siemens 送检人群 56.7%、孕妇 82.8%、HIV+ 94–98%），交叉反应研究多以 "反应性 = 比较法亦阳性" 呈现（Siemens 297 份 88 反应，比较法 88；Roche IgG 249 份仅 HTLV 1 份不一致）。Abbott IgG 187 份：Toxoplasma IgG 1 份反应（比较法亦阳）、Parvovirus B19 1 份 equivocal（比较法阳）。Roche IgM 205 份：autoimmune 2/7 与 influenza vaccine 1/10 出现 border/反应（比较法阴性），文件注明 "Potential cross-reactivity with autoimmune markers and antibodies against influenza vaccination could not be ruled out"；Roche IgG 注明 E. coli 载体与自身免疫标志物交叉不能排除，VZV/Measles/Mumps/Parvo-B19 IgG 未评估。EBV：Abbott 0/10、Siemens EBV IgG 0/13、Roche IgM EBV 0/14。
- 背景（非申报文件）：EBV/HHV-6 等疱疹病毒急性感染可产生多克隆 IgM 反应导致 CMV IgM 假阳性；孕早期原发感染风险评估通常需 IgG avidity（本组 510(k) 文件均未涉及 avidity）。

### 6. 样本类型与样本要求
- ARCHITECT：血清、SST、Li-heparin、Li-heparin PST、K3-EDTA（≥40 供者配对，接受标准 <6 AU/mL 差 ≤0.6、6–<15 差 ≤1.5、≥15 %差 ≥-10%）；室温 3 天、2–8 °C 7 天、-20 °C 28 天。
- ADVIA Centaur：血清、K2-EDTA、Li-heparin（38 组，斜率 1.01，r 0.96–0.99）。
- Elecsys IgG：血清、Li-heparin、K2/K3-EDTA（38 组，阴性差 <0.2 COI、阳性 ±20%）；Elecsys IgM：加 SST；样本 2–8 °C 4 周、20–25 °C 7 天、-20 °C 3 个月、5 次冻融。
- LIAISON：仅血清（文件）。

### 7. 分析性能验证所依据的标准

| K 号 | 文件所列标准/指南 | 对应项目 |
|---|---|---|
| K220949 | CLSI EP05-A3 (2014)；EP07 3rd ed. (2018)；EP12-A2 (2008)；EP25-A (2009)；EP37 1st ed. (2018)；以及 FDA 510(k) 格式、eCopy、Replacement Reagent and Instrument Family Policy、软件、互操作性、网络安全、OTS 软件、Benefit-Risk 等指南 | 20 天 3 批精密度（校准曲线存储 31–35 天）；3 中心重复性；4 种内源 + 11 种药物干扰（含 biotin 3510 ng/mL、ganciclovir、valganciclovir、foscarnet、cidofovir）；187 份交叉反应；定性方法比较；试剂稳定性 |
| K181213 | "N/A"；正文引 CLSI EP05-A3（精密度、重复性）、EP07-A2（干扰） | 20 天 80 重复；3 中心 2 批；8 种干扰（biotin 4500 ng/mL）；297 份交叉反应 |
| K131605 / K220911 | K131605：CLSI EP5-A2、EP15-A2；K220911：CLSI EP05-A3、EP07 3rd ed.、EP37、EP12-A2 | 21 天精密度；3 中心；6 种干扰（biotin 100→1200 ng/mL）；249 份交叉反应；hook；280 份新旧版本一致性 |
| K142133 / K163569 | K142133：CLSI EP5-A2、EP7-A2、EP15-A2；K163569：CLSI EP5-A3、EP15-A2、**EP17-A2**（LoB/LoD） | 精密度；3 中心；4 种内源 + RF + 19 种药物；205 份交叉反应；DTT IgM 特异性；hook；LoB/LoD |
| K162969 | Special 510(k) 设计控制声明；稳定性、精密度（20 天）、基质效应（控制品） | 无临床性能 |

### 8. 临床验证设计与结果

**ARCHITECT CMV IgG（K220949）**——989 份（美国 791：常规 591、孕妇 200；OUS 198），3 站点，比较法为一种 FDA 清关市售 CMV IgG（equivocal 者再用 2 种检测 2/3 共识）。
- 常规送检：PPA 97.7%（514/526；96.1–98.7），NPA 99.2%（261/263；97.3–99.8）；7 份 ARCHITECT equivocal/比较法阳性计入分母。
- 孕妇：PPA 99.0%（98/99；94.5–99.8），NPA 100%（102/102）。
- CDC CMV IgG 盘 80 份：PPA 100%（91.59–100），NPA 92.11%（78.62–98.34），总 96.25%。
- 美国人群阳性率 64.2%，grayzone 1.3%。

**ADVIA Centaur CMV IgG（K181213）**——1,842 份（前瞻 1,699：送检 684、孕妇 348、儿科 229、HIV+ 44、移植 394；回顾 HIV+ 143），6 采集点，比较法为市售 CMV IgG（equivocal 再用 2 种检测共识）。
- 合并前瞻：PPA 99.6%（1037/1041；99.0–99.9），NPA 96.8%（637/658；95.2–98.0）。
- 送检 PPA 99.5%/NPA 95.8%；孕妇 PPA 100%（287/287）/NPA 98.4%（60/61）；儿科 PPA 98.8%/NPA 98.6%；HIV+ 100%/100%（1/1）；移植术前 100%/100%（151/151；87/87）；移植术后 PPA 99.0%/NPA 90.7%（49/54）；回顾 HIV+ 100%/100%。
- CDC 盘 80 份：100% 一致（39/39，41/41）。

**Elecsys CMV IgG（K131605）**——美国 605 份（送检 400、孕妇 205），2 外部 + 1 内部站点，比较法为市售 CMV IgG。
- 送检人群：PPA 98.94%（280/283；96.93–99.78），NPA 92.86%（299/322；89.47–95.42）（12 份 Elecsys+/比较法− 与 8 份 equivocal 计入）。
- 孕妇：PPA 96.55%（84/87；90.25–99.28），NPA 100%（118/118；96.92–100）。
- K220911（生物素改良）：280 份 3 批 vs 原版，PPA 96.5–100%、NPA 98.8–100%（保守法）。

**Elecsys CMV IgM（K142133）**——比较法为 2–3 种 FDA 清关 CMV IgM 的 2/3 共识；Elecsys equivocal 一律按不利计。
- 前瞻疑似感染 418 份：NPA 97.2%（384/395；95.1–98.6），PPA 87.0%（20/23；66.4–97.2）（3 份共识 equivocal/Elecsys 阴性按假阴性计）。
- 前瞻孕妇 199 份：NPA 98.5%（190/193；95.5–99.7），PPA 0/6（0.00–45.9；共识仅 1 阳 + 5 equivocal，全部按不利计）。
- 预选 IgM 阳性 134 份：PPA 100%（134/134；97.3–100）。
- 送检人群 IgM 反应率 1.67%、border 1.17%。K163569（cobas e 801）223 份 vs e 601：PPA 100%、NPA 100%，border 一致 75%（6/8）。

**LIAISON CMV IgG（K162969）**——Special 510(k)，仅控制品基质由 5% 血清改为 100% 血清/去纤维蛋白血浆、开瓶稳定 4→8 周；性能 "No Change"，引用 K040290（文件未能获取）。

### 9. 厂家间差异与要点
1. **IgG 灰区设计**：Abbott 设宽灰区 6–<15 AU/mL（约 2.5 倍跨度）并建议以 IgM 或复采解决；Roche 0.5–<1.0 COI；Siemens **无灰区**（以 PPA/NPA ≥95% 为 cutoff 设定准则）。对 "immune status" 用途而言，灰区处理直接影响孕妇/移植前的易感判定（背景推断）。
2. **IgM 判定**：Roche 三分类（<0.7 / 0.7–<1.0 / ≥1.0），cutoff 以低流行 931 份保特异 + 152 份原发感染各期保灵敏；FDA 要求 equivocal 全部按不利计，导致孕妇 PPA 0/6，决策摘要仍接受（预选阳性 134/134）。
3. **抗原**：病毒裂解物（Abbott、Siemens）vs 重组抗原（Roche）；Roche 因此对 E. coli 抗体交叉需加限制说明。
4. **特殊人群**：Siemens 唯一提供儿科（2–21 岁）与移植前/后亚组；Abbott 强调育龄妇女；Roche 明确免疫抑制者未评估。所有产品一致声明不用于供者筛查。
5. **CDC 盘**：Abbott IgG NPA 92.11%（3/38 假阳）、Siemens 100%；IgM 产品未见 CDC 盘数据。
6. **生物素耐受**：Roche 两次改良（IgG：100→1200 ng/mL；IgM 100 ng/mL）；Abbott 3510 ng/mL；Siemens 4500 ng/mL。

### 10. 来源
- K220949：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220949.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220949.pdf
- K181213：https://www.accessdata.fda.gov/cdrh_docs/reviews/K181213.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181213.pdf
- K131605：https://www.accessdata.fda.gov/cdrh_docs/reviews/K131605.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K131605.pdf
- K220911：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220911.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220911.pdf
- K142133：https://www.accessdata.fda.gov/cdrh_docs/reviews/K142133.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K142133.pdf
- K163569：https://www.accessdata.fda.gov/cdrh_docs/reviews/K163569.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K163569.pdf
- K162969：https://www.accessdata.fda.gov/cdrh_docs/reviews/K162969.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162969.pdf（原始 K040290 决策摘要/510(k) summary 均未能下载）
