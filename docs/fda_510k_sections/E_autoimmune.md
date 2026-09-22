# 靶点组 E（自身免疫血清学）：FDA 510(k) 申报信息整理

> 数据来源说明：本节全部申报数据（cutoff、人群、n、灵敏度/特异度、PPA/NPA、标准清单）均摘自 FDA 公开的 510(k) 决策摘要（Decision Summary，`cdrh_docs/reviews/<K>.pdf`）或 510(k) Summary（`cdrh_docs/pdfYY/<K>.pdf`），文本抽取件存于 `scratchpad/txt/<K>.txt`。文件未写明的内容一律标「文件未载明」；凡来自背景知识（指南、免疫学机理）的说明均单独标注「背景（非申报文件）」。openFDA 累计清关数以 `fda_fetch.py list <CODE>` 首行 `total=` 为准（查询日期 2026-09-22）。

---

## 抗环瓜氨酸肽抗体（anti-CCP / ACPA）

### 1. 法规定位
- Product code **NHX**（Antibodies, Anti-Cyclic Citrullinated Peptide (CCP)）；21 CFR **866.5775**（Rheumatoid factor immunological test system）；**Class II**；Panel: Immunology (82)。校准品/质控品另配 JIT/JIX（Class II 校准品）、JJX/JJY（Class I 质控品，Clinical Chemistry 75）。
- openFDA `product_code:NHX` 累计清关数 **total=21**（2002-04-29 K020414 QUANTA Lite CCP ELISA 至 2016-08-04 K153551）。注：2016 年以后未再见以 NHX 为主代码的新清关记录。
- 相关但独立的代码：OQZ（anti-mutated citrullinated vimentin, 866.5775, total=1）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K143754 | 2015-09-21 | Inova Diagnostics / QUANTA Flash CCP3（BIO-FLASH） | 顺磁微粒化学发光（CIA），半定量，CU 单位 | 血清 | IgG anti-CCP3；aid in the diagnosis of RA，与临床及其他实验室检查联合 | 决策摘要 + 510(k) summary |
| K153551 | 2016-08-04 | Axis-Shield（申请人）/ Siemens ADVIA Centaur Anti-CCP IgG (aCCP) | 吖啶酯化学发光（CMIA），半定量，U/mL | 血清、K2-EDTA 及 Li-heparin 血浆 | 同上；“autoantibody levels represent one parameter in a multi-criteria diagnostic process” | 决策摘要 + 510(k) summary |
| K121576 | 2013-04-04 | Siemens / IMMULITE 2000 Anti-CCP IgG | 固相化学发光免疫测定（碱性磷酸酶），半定量，U/mL | 血清、EDTA/Li-heparin 血浆 | 同上 | 决策摘要 + 510(k) summary |
| K093954 | 2010-08-24 | Bio-Rad / BioPlex 2200 Anti-CCP Kit | 多重流式微珠免疫分析（PE 荧光），半定量，U/mL | 血清、EDTA/钠肝素血浆 | 同上 | 决策摘要（无 510(k) summary 抽取件） |
| K083868 | 2009-09-25 | Axis-Shield（申请人）/ Abbott ARCHITECT Anti-CCP | CMIA，半定量，U/mL | 血清、SST、K2-EDTA、Li-heparin 血浆 | 同上 | 决策摘要 |
| K081338 | 2008-09-11 | Roche / Elecsys Anti-CCP（ECLIA） | 电化学发光，半定量，U/mL | 血清、K3-EDTA/Li-heparin 血浆 | 同上 | 决策摘要 |
| K061165 | 2006-08-03 | Phadia（Sweden Diagnostics）/ EliA CCP（ImmunoCAP 100/250） | 荧光酶免疫测定（FEIA），半定量，EliA U/mL | 血清、肝素/EDTA/枸橼酸血浆 | 同上 | 决策摘要 |
| K070789 | 2007-12-31 | Euroimmun / Anti-CCP ELISA (IgG) | 手工 ELISA，定性（OD ratio）/半定量（RU/mL） | 血清、EDTA/肝素/枸橼酸血浆 | 同上 | 决策摘要 |
| K110296 | 2011-08-18 | Axis-Shield / Axis-Shield Anti-CCP (FCCP600) | 手工 ELISA，定性/半定量，U/mL | 血清、SST、EDTA、Li-heparin、枸橼酸血浆 | 同上 | 决策摘要 + 510(k) summary |

### 3. 预期用途与声明类型
- 所有产品均为 **处方用（Rx only）**、**辅助诊断（aid in the diagnosis of rheumatoid arthritis）**，且明确要求 “in conjunction with clinical findings and other laboratory tests / other clinical information”。Siemens/Abbott 系列进一步写明 “autoantibody levels represent one parameter in a multi-criteria (multicriterion) diagnostic process, encompassing both clinical and laboratory-based assessments”。
- 分析物一律限定为 **IgG 类** anti-CCP（QUANTA Flash CCP3 的前代 QUANTA Lite CCP3.1 曾含 IgG/IgA，K072944，本次未取文件）。
- 测量类型：绝大多数为 **semi-quantitative**（无国际标准物质，单位为厂家任意单位）；Euroimmun 与 Axis-Shield ELISA 另提供 qualitative（OD ratio）判读方案。
- 无任何产品声明用于 RA 预后分层、治疗监测或“早期 RA”单独适应症（Axis-Shield K110296 仅在临床数据中分列 early RA / established RA 亚组）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 (K 号) | Cutoff 及判读 | 建立方式 | 参考人群 |
|---|---|---|---|
| QUANTA Flash CCP3 (K143754) | ≥20 CU 阳性，<20 CU 阴性；无灰区。510(k) summary 载明 20 CU 经 IUIS/CDC ACPA 参考试剂（IS2723，赋值 100 U/mL，本法测得 379.5 CU）换算 ≈ **5.3 U/mL**；软件可自动输出 U/mL | **参考人群第 99 百分位** | 210 例（170 表观健康献血者 + 40 例其他疾病对照） |
| ADVIA Centaur aCCP (K153551) | ≥5 U/mL 阳性，<5 U/mL 阴性；无灰区 | **与 predicate ARCHITECT 对齐**（“Assay cut-off 5.0 U/mL — Same”）。文件引用 EP24-A2 但决策摘要未给 ROC 数据 | 期望值：127 例健康者 99th percentile 1.27 U/mL（远低于 cutoff） |
| IMMULITE 2000 (K121576) | ≥4.00 U/mL Reactive，<4 U/mL Non-reactive | **ROC 分析**，“balanced consideration of sensitivity and specificity” | 388 例：212 表观健康 + 106 RA 阳性 + 70 例 predicate 法 anti-CCP 阳性；另 200 例健康者非参数 99th percentile = 4.06 U/mL，1.5% 阳性 |
| BioPlex 2200 Anti-CCP (K093954) | 3.0 U/mL（<3.0 阴性，>3.0 阳性）；无灰区 | **以 predicate（Axis-Shield DIASTAT）结果为金标准做一致性 + ROC**，选 3.0 U/mL 使 PPA 92.9%、NPA 98.2% | 1394 例：177 正常、504 例 RF 已测或阳性、82 例 >70 岁、287 RA、344 非 RA；该批样本未再用于其他研究 |
| ARCHITECT Anti-CCP (K083868) | ≥5.0 U/mL 阳性 | **为与 predicate AxSYM 及其他 anti-CCP 方法保持一致而选定**；vs AxSYM 的 ROC AUC 0.873 (0.849–0.897) | 199 例健康者 <0.5–2.5 U/mL |
| Elecsys Anti-CCP (K081338) | ≥17 U/mL 阳性（predicate IMMUNOSCAN 为 25 U/mL） | **ROC 最优 cutoff**，AUC 0.85；该点灵敏度 67.7%、特异度 97.0% | 792 例确诊 RA + 420 例无症状健康者 + 907 例其他风湿/非风湿病 |
| EliA CCP (K061165) | <7 U/mL 阴性；**7–10 U/mL equivocal**；>10 U/mL 阳性 | 预先定义的 cutoff，用健康人群分布**确认**（99th percentile 6.2 U/mL、mean+3SD 10.1 U/mL 均低于/落在灰区） | 400 例表观健康高加索成人献血者，性别年龄均衡 |
| Euroimmun Anti-CCP ELISA (K070789) | 半定量：>5 RU/mL 阳性；定性：OD ratio >1.0 阳性 | **ROC 分析**，5 RU/mL 处灵敏度 78.5%、特异度 98.2% | 419 RA + 1144 对照（400 献血者、28 PsA、35 其他关节炎、108 SLE、106 SjS、98 SSc、159 自身免疫甲状腺炎、25 Wegener、126 抗 parvovirus B19 阳性、54 病毒性肝炎、5 抗 HIV 阳性） |
| Axis-Shield FCCP600 (K110296) | 半定量：>5 U/mL 阳性；定性：OD ratio <0.95 阴性、0.95–1.0 borderline（建议复测）、>1.0 阳性 | **沿用 predicate DIASTAT (K023285) cutoff**；510(k) summary 载 ROC AUC 0.910 (0.881–0.940) vs DIASTAT 0.903 | 514 例（229 RA、150 健康、135 非 RA） |

要点：
- 各厂家单位（CU、U/mL、RU/mL、EliA U/mL）互不可换算；文件均声明 “no recognized standard or reference material for anti-CCP”。仅 Inova（K143754）报告以 IUIS/CDC ACPA 参考试剂建立 CU↔U/mL 换算。
- cutoff 建立路径有三类：健康/疾病对照百分位（Inova 99th）、ROC（Siemens IMMULITE、Roche、Euroimmun）、与 predicate 对齐（Abbott、Siemens ADVIA、Bio-Rad 以 predicate 为标准的 ROC、Axis-Shield）。
- 灰区：仅 Phadia EliA（7–10 U/mL）与手工 ELISA 定性方案（Axis-Shield ratio 0.95–1.0）设 equivocal/borderline；自动化 CIA/CMIA/ECLIA 产品均为单一 cutoff。

### 5. 生物学 / 免疫学依据
- 申报文件所载：所有产品抗原均为 **synthetic cyclic citrullinated peptide（第二代，“CCP, second generation”）**，Inova 标为 CCP3；检测 IgG 类自身抗体；申报文件将 anti-CCP 归入 21 CFR 866.5775（RF 试验系统）并强调 RA 诊断为多参数过程。Bio-Rad K093954 引用文献说明 anti-CCP 亦可见于原发性干燥综合征（非侵蚀性滑膜炎）并提示 Centromere B 阳性（23%）、SS-A 阳性（12%）、骨髓瘤 IgG（30%）样本可能出现交叉阳性；Siemens IMMULITE K121576 引用 Ann Intern Med 2007/2010 两篇 Meta 分析佐证其灵敏度/特异度水平。
- **与 RF 的关系（申报文件）**：(1) 法规上 anti-CCP 与 RF 同属 866.5775。(2) 各产品均做 RF 干扰试验：ADVIA/IMMULITE/BioPlex/ARCHITECT 至 200 IU/mL 无干扰、Elecsys <150 IU/mL、EliA IgM-RF 550 IU/mL、Euroimmun 用 10 例 IgM-RF 阳性 RA 血清做交叉反应“无阳性”；IMMULITE 说明书加注 “RF may have unpredictable effects on patient samples with low concentrations of anti-CCP”。(3) BioPlex cutoff 研究纳入 504 例 “RF tested or positive” 患者以挑战特异性；IMMULITE 临床研究说明入组仅依诊断、predicate anti-CCP 与 RF 状态对申办方设盲。(4) **RF 与 anti-CCP 联合判读的诊断性能**：所有文件均未载明。
- 背景（非申报文件）：ACPA 靶向 PAD 酶介导的瓜氨酸化蛋白（filaggrin、fibrinogen、vimentin、α-enolase 等），CCP 为合成环肽模拟表位；ACPA 与 RF 一同列入 2010 ACR/EULAR RA 分类标准血清学项（高滴度 >3×ULN 计 3 分），可在临床发病前数年出现，与侵蚀性病程相关。

### 6. 样本类型与样本要求
- QUANTA Flash CCP3：**仅血清**（predicate QUANTA Lite 允许枸橼酸/EDTA 血浆）；稳定性 RT 48 h、2–8 °C 21 d、≤−20 °C 3 次冻融。
- ADVIA Centaur：血清、SST、K2-EDTA、Li-heparin（各约 50 对，Passing-Bablok 斜率 0.98–1.01）；样本稳定性沿用 K083868（RT 22 h、2–8 °C 14 d、2 次冻融）。
- IMMULITE 2000：血清、SST、EDTA、Li-heparin（39 组匹配，Deming 斜率 1.01–1.02）；20–25 °C 2 d、2–8 °C 7 d、血清 −20 °C 6 个月、血浆 4 个月；短抽血 5 倍抗凝剂过量 8 例中 7 例回收 ±10%。
- BioPlex 2200：血清、EDTA、肝素血浆（41 组匹配，EP9-A2；斜率 0.96）。
- ARCHITECT：血清、SST、K2-EDTA、Li-heparin（18 组加标）；2–8 °C 7 d 或 30 °C 22 h，长期 ≤−20 °C。
- Elecsys：血清、K3-EDTA（y=1.10x−0.98，n=31）、Li-heparin（y=1.03x−0.14，n=33）。
- EliA CCP：血清、肝素/EDTA/枸橼酸血浆（50 组，r ≥0.996）。
- Euroimmun：血清、EDTA/肝素/枸橼酸血浆（4 对系列稀释，Passing-Bablok）。
- Axis-Shield FCCP600：血清、SST、K-EDTA、Li-heparin、Na-citrate（19 组匹配 + 16 组加标）。

### 7. 分析性能验证所依据的标准（文件 “Standards/Guidance Documents Referenced” 原文）
| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K143754 QUANTA Flash CCP3 | EP05-A2（精密度）、EP06-A（线性）、EP07-A2（干扰）、EP09-A3（方法比对）、EP17-A2（LoB/LoD）、C28-A3（参考区间/cutoff） | 20 天×2 次×2 复精密度；3 中心再现性；线性 4.6–2776.8 CU；LoB 410 RLU；胆红素/血红蛋白/甘油三酯干扰；CDC ANA 12 份参考血清交叉反应 |
| K153551 ADVIA Centaur | EP5-A2、EP6-A、EP7-A2、EP9-A3、EP17-A2、**EP24-A2（ROC 诊断准确性）**、EP25-A（稳定性） | 精密度、线性 0.40–200 U/mL、LoB 0.24/LoD 0.40 U/mL、干扰（含 biotin 500 ng/dL、RF 200 IU/mL）、11 种自身抗体交叉反应 |
| K121576 IMMULITE 2000 | EP5-A2、EP06-A、EP17-A | 2 仪器×2 批×20 天精密度（n=320/样本）；线性 1.5–200 U/mL；LoB 0.26、LoD 1.46、功能灵敏度 2.34 U/mL |
| K093954 BioPlex 2200 | EP5-A、EP7-A、EP6-A、EP12-A2（定性一致性）、EP14-A2（基质效应）、EP15-A2（用户验证）、EP17-A | 三基质精密度、线性 0–300 U/mL、LoD 0.2 U/mL、13 种干扰物、163 例 ANA 阳性样本交叉反应 |
| K083868 ARCHITECT | EP5-A2、EP6-A、EP7-A2、EP9-A2、EP17-A | 三型号仪器精密度、线性 0.1–246.5 U/mL、LoB 0.02/LoD 0.11 U/mL、ORDAC 自动稀释、20 份自身抗体加标交叉反应 |
| K081338 Elecsys | EP17-A、EP5-A2 | LoB ≤7、LoD ≤8、LoQ 8 U/mL；21 天精密度；无钩效应至 7000 U/mL |
| K061165 EliA CCP | 仅引用 FDA 软件指南（Guidance for the Content of Premarket Submissions for Software Contained in Medical Devices） | 6 批次×3 孔批精密度；1:2–1:16 稀释线性（2 个比值超标，说明书注明并非所有血清可线性稀释）；校准品溯源 WHO IRP 67/86 |
| K070789 Euroimmun | None referenced | 20 复孔批内、4 次×3 天批间、3 批次；LoB 0.4/LoD 2.5 RU/mL；干扰回收 85–115% |
| K110296 Axis-Shield FCCP600 | EP05-A2、EP06-A、EP07-A2、EP09-A2、EP17-A | LoD 1.04 U/mL；多聚回归线性；Passing-Bablok 基质比对 |

### 8. 临床验证设计与结果
| 产品 | RA 组（分类标准） | 对照组构成 | 临床灵敏度 (95% CI) | 临床特异度 (95% CI) | 与 predicate/参考法一致性 |
|---|---|---|---|---|---|
| QUANTA Flash CCP3 (K143754) | 352 例 “diagnosed with RA”（分类标准：文件未载明） | 376 例：AS 13、OA 49、PMR 20、PsA 14、SLE 53、SjS 20、AIH 19、UC 11、乳糜泻 20、Lyme 25、parvovirus 12、沙门菌 10、HBV 22、HCV 13、其他 75 | **70.7%** (65.8–75.2) | **96.5%** (94.2–98.0)；对照阳性率 3.5%（PMR 10%、PsA 7.1%、SLE 3.8%） | vs QUANTA Lite CCP3 ELISA（n=420 在测量范围内）：PPA 94.2% (90.3–96.6)，NPA 90.3% (85.4–93.7) |
| ADVIA Centaur (K153551) | 307 例，**2010 ACR 分类标准** | 460 例：AS 26、OA 74、PMR 13、PM 3、PsA 27、反应性关节炎 24、SLE 92、SSc 8、SjS 22、桥本 10、血管炎 6、纤维肌痛 4、痛风 2、骨髓瘤/肿瘤 23、Crohn 13、UC 12、Lyme 7、West Nile 10、Chikungunya 15、登革 13、EBV 15、HCV 41 | **68.08%** (62.54–73.26) | **97.17%** (95.22–98.49)；对照阳性 2.8%（PsA 7.4%、OA 6.8%、SLE 4.3%、SjS 4.5%） | vs ARCHITECT（n=253）：PPA 96.80% (92.06–98.75)，NPA 96.88% (92.24–98.78) |
| IMMULITE 2000 (K121576) | 1048 例；两临床中心由风湿科医师按 **1987 ACR 标准**前瞻入组（791 例，其中 204 例病程 <2 年，多数在用抗风湿药）；供应商样本“as diagnosed” | 464 例：AS 31、桥本 16、Crohn 9、DM 6、EBV 5、痛风 10、OA 93、PMR 26、PsA 49、SSc 15、SjS 18、SLE 79、UC 8、Wegener 6、其他 78 等 | **63.6%** (60.0–66.6) | **97.0%** (95.0–98.3) | 外部中心 1255 例同测 predicate：predicate 63.8%/96.1% vs IMMULITE 58.9%/97.0%；范围内 255 例 PPA 86.9% (81.0–91.5)、NPA 46.2% (28.8–64.5)（阴性样本数少所致） |
| BioPlex 2200 (K093954) | 496 例 “previously diagnosed RA”（标准未载明） | 300 健康献血者 + 201 例其他风湿病（Crohn 18、纤维肌痛 15、痛风 19、炎性关节炎 16、OA 8、SSc 18、SjS 21、SLE 31、UC 17、CREST 2、Wegener 2 等） | **83.1%** (79.5–86.1) | **97.8%** (96.1–98.8) | vs DIASTAT（n=822）：PPA 97.5% (95.4–98.7)，NPA 91.4% (88.5–93.7) |
| ARCHITECT (K083868) | 496 例 confirmed RA（标准未载明） | 200 健康 + 299 例其他疾病（AS、桥本、Crohn、DM、EBV、Lyme、OA、PMR、PM、PsA、反应性关节炎、SSc、SjS、SLE、UC） | **70.6%**（文件写 66.3–77.4）；predicate AxSYM 70.8% (66.5–74.7) | 健康 99.5%（1/200 阳性）、疾病对照 97.3%、合计 **98.2%** | vs AxSYM（n=995）：PPA 98.9%，NPA 99.5%；AUC 0.873 vs 0.872 |
| Elecsys (K081338) | ROC/临床：792 例 confirmed RA；方法比对：614 例按 **ARA 1987 修订标准** | 420 无症状健康 + 907 例其他风湿/非风湿病；方法比对 319 健康 + 673 非 RA | **67.7%** | **97.0%**；AUC 0.85 | vs IMMUNOSCAN（n=1606）：PPA 94.3% (91.7–96.2)，NPA 98.4% (97.5–99.1) |
| EliA CCP (K061165) | 82 例 clinically confirmed RA | 452 例：CTD 120（UCTD 11、SLE 30、SSc 39、SjS 20、MCTD 10、肌炎 10）、感染 185（CMV 25、EBV 30、HBV 10、HCV 10、HIV 15、其他 95）、其他 147（OA 35、血管炎 21、AITD 20、Crohn 20、UC 20、杂项 31） | **87.8%** (72/82) | **96.7%** (437/452)；SSc/SjS 各 3/39、3/20 阳性 | vs DIASTAT（n=80）：PPA 98.1%，NPA 89.3% |
| Euroimmun ELISA (K070789) | 419 例 RA（方法比对 259 例 “according to ACR criteria”） | 1144 例（见第 4 节） | **78.5%** (74.3–82.4) | **98.2%**；献血者 99.5%、PsA 100%、SLE 97.2%、SjS 98.1%、SSc 96.9% | vs IMMUNOSCAN（259 RA）：PPA 98.6%，NPA 90.7% |
| Axis-Shield FCCP600 (K110296) | 229 例，**ACR 1987** | 150 健康 + 135 非 RA（炎性多关节炎 41、EBV 18、桥本 17、SjS 16、SLE 16 等） | Early RA 43 例 **63%**；established RA 186 例 **82%**；全部 **78%** | 健康 99.3%、非 RA 97.8%、合计 98.6% | vs DIASTAT（n=514）：PPA 99.4%，NPA 98.8%；AUC 0.910 |

### 9. 厂家间差异与要点
1. **灵敏度差异主要源于 RA 人群构成**：前瞻性、治疗中、含早期 RA 的队列（IMMULITE 63.6%、Axis-Shield 早期 RA 63%）显著低于 “previously diagnosed” 库存样本（BioPlex 83.1%、EliA 87.8%）；特异度普遍 96–98%。
2. **cutoff 数值不可跨平台比较**（3、4、5、17、20 单位不等），但各家都用 predicate 一致性 + 健康人群百分位/ROC 两条证据链；FDA 接受“沿用 predicate cutoff”（Abbott、ADVIA、Axis-Shield）而不强制重新做 ROC。
3. **灰区**只见于 Phadia EliA 和手工 ELISA 定性方案，自动化平台一律单点 cutoff。
4. **样本基质**：Inova 自动化 CIA 只申报血清；其余大多数含 EDTA/肝素血浆，Euroimmun/EliA/Axis-Shield 还含枸橼酸血浆。
5. **交叉反应**：Bio-Rad 报告 SS-A、centromere B、骨髓瘤 IgG 样本可能出现假阳性；Abbott/Siemens 加标 11 种自身抗体无干扰；QUANTA Flash 对照组 PMR 10%、PsA 7.1% 阳性。
6. **溯源**：仅 Inova 报告与 IUIS/CDC ACPA 参考试剂关联（20 CU≈5.3 U/mL）；Bio-Rad 校准品赋值来自 predicate DIASTAT；其余为厂内任意单位。
7. **RF 状态**从未作为入组/分层变量报告；文件中 RF 仅出现在干扰试验与 BioPlex cutoff 人群。

### 10. 来源
- K143754：https://www.accessdata.fda.gov/cdrh_docs/reviews/K143754.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K143754.pdf
- K153551：https://www.accessdata.fda.gov/cdrh_docs/reviews/K153551.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153551.pdf
- K121576：https://www.accessdata.fda.gov/cdrh_docs/reviews/K121576.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K121576.pdf
- K093954：https://www.accessdata.fda.gov/cdrh_docs/reviews/K093954.pdf
- K083868：https://www.accessdata.fda.gov/cdrh_docs/reviews/K083868.pdf
- K081338：https://www.accessdata.fda.gov/cdrh_docs/reviews/K081338.pdf
- K061165：https://www.accessdata.fda.gov/cdrh_docs/reviews/K061165.pdf
- K070789：https://www.accessdata.fda.gov/cdrh_docs/reviews/K070789.pdf
- K110296：https://www.accessdata.fda.gov/cdrh_docs/reviews/K110296.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K110296.pdf

---

## 抗核抗体筛查（ANA screen：IFA HEp-2 与固相法）

### 1. 法规定位
- 21 CFR **866.5100**（Antinuclear antibody immunological test system），**Class II**，Panel: Immunology (82)。
- Product codes 与 openFDA 累计清关数：
  - **DHN** Antinuclear Antibody, Indirect Immunofluorescent, Antigen, Control（IFA HEp-2）：**total=110**（含大量 Bio-Rad Liquichek 质控品记录）。
  - **LKJ** Antinuclear Antibody, Antigen, Control（多重微珠/CLIA 等）：**total=50**。
  - **LJM** Antinuclear Antibody (Enzyme-Labeled), Antigen, Controls（ELISA 筛查）：**total=112**。
  - 相关：PIV（Automated indirect immunofluorescence microscope and software-assisted system，21 CFR 866.4750，NOVA View 经 DEN140039 建立）、NVI（BioPlex Medical Decision Support Software，21 CFR 862.3100，total=1）。
- FDA 专门指南：*Guidance for Industry and FDA Staff: Recommendations for Anti-Nuclear Antibody (ANA) Test System Premarket (510(k)) Submissions*（2009-01-22），被 K150155、K131791、K131185、K131330、K112996 明确引用。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K150155 | 2015-04-08 | Inova / NOVA Lite DAPI ANA Kit（手工或 NOVA View 自动荧光显微镜） | IFA HEp-2，定性/半定量滴度；DAPI 复染 | 血清 | aid in the diagnosis of SLE and other systemic rheumatic diseases；软件结果须经培训操作者确认 | 决策摘要 + 510(k) summary |
| K131791 | 2014-02-26 | Euroimmun / IFA 40: HEp-20-10（BIOCHIP/TITERPLANE） | IFA，HEp-20-10 细胞 | 血清 | aid in the diagnosis of systemic rheumatic diseases | 决策摘要 + 510(k) summary |
| K172745 | 2018-06-05 | Immco / ImmuGlo HEp-2 Elite IFA（标准 HEp-2 + PSIP1/DFS70 敲除 HEp-2 按 9:1 混合） | IFA | 血清 | 同上；ICAP AC-1~AC-28 核型 | 决策摘要 + 510(k) summary |
| K201956 | 2022-04-29 | Zeus Scientific / Zeus IFA ANA HEp-2 Test System + Zeus dIFine（自动读片） | IFA + 软件判读 | 血清 | 同上 | 决策摘要 + 510(k) summary（Statement） |
| K153117 | 2016-07-28 | Aesku / AESKUSLIDES ANA HEp-2-Gamma + HELIOS 全自动 IFA 系统 | IFA（Fc-γ 特异性结合物）+ 软件判读 | 血清 | 同上 | 决策摘要 + 510(k) summary |
| K041658 | 2004-12-20 | Bio-Rad / BioPlex 2200 ANA Screen（13 微珠：dsDNA、chromatin、Rib-P、SS-A 60/52、SS-B、Sm、Sm/RNP、RNP-A/68、Scl-70、Jo-1、CENP-B）；K113610（2012，Special 510(k)，加 MDSS 及 QC 频次变更） | 多重流式微珠免疫分析 | 血清、EDTA/肝素血浆 | qualitative ANA screen + 定量 dsDNA + 10 项半定量；aid in the diagnosis of systemic rheumatic diseases | 决策摘要（K041658）；K113610 仅 ODE 备忘 + 510(k) summary |
| K131185 | 2013-07-15 | Euroimmun / ANA Screen ELISA (IgG)（dsDNA、组蛋白、Rib-P、nRNP/Sm、Sm、SS-A、SS-B、Scl-70、Jo-1、centromere 混合抗原） | ELISA，定性 ratio | 血清、EDTA/Li-heparin/枸橼酸血浆 | aid in the diagnosis of MCTD、SLE、SjS、PSS、PM/DM | 决策摘要 + 510(k) summary |
| K131330 | 2014-01-28 | Gold Standard Diagnostics / ANA Screen ELISA（HEp-2 裂解物 + 纯化/重组 ENA） | ELISA，定性 ratio，含 equivocal | 血清 | 同上 | 决策摘要 + 510(k) summary |
| K083188 | 2009-03-13 | Phadia / Varelisa ReCombi ANA Screen（8 抗原；Sm 抗原改为 16 肽） | EIA，定性（cut-off calibrator） | 血清、枸橼酸/EDTA 血浆（肝素禁用） | aid in the diagnosis of SLE、SSc、MCTD、SS、PM/DM | 决策摘要（device modification） |

### 3. 预期用途与声明类型
- 均为 **Rx only**、**qualitative and/or semi-quantitative（IFA 滴度）** 的 **辅助诊断** 声明，措辞为 “aid in the diagnosis of systemic lupus erythematosus and other systemic rheumatic diseases / systemic autoimmune diseases, in conjunction with other serological tests and clinical findings”。固相法多列出目标疾病清单（SLE、MCTD、SjS、SSc/PSS、PM/DM；BioPlex 另含 UCTD、RA、CREST、Raynaud）。
- 自动读片系统（NOVA View、dIFine、HELIOS）附加特殊条件：“All software-aided results must be confirmed by the trained operator”。
- BioPlex MDSS（K043341/K113610）为可选的 kNN 模式识别工具，输出 Negative / No Association / Association with Disease（≤2 个疾病分类），不改变各分析物判读。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**IFA HEp-2**

| 产品 | 筛查稀释度 | “阳性”定义 | cutoff 依据 |
|---|---|---|---|
| NOVA Lite DAPI ANA (K150155) | **1:80**（predicate NOVA Lite HEp-2 为 1:40） | 1:80 有特异荧光并报核型（homogeneous、speckled、centromere、nucleolar、nuclear dots；NOVA View 另有 unrecognized）；强度 1+~4+；阳性者以 1:80→1:2560 倍比稀释定终点滴度 | “1:80 selected to provide optimal clinical sensitivity and specificity”；NOVA View **LIU cutoff=48**：120 例健康献血者 LIU 非参数 **90th percentile 48.8**（90% CI 31–78），以手工读片特异度为参照；cutoff 固化于软件不可修改 |
| Euroimmun IFA 40 (K131791) | **1:40** | ≥1+ 荧光且有可辨核型；0 = 阴性；滴度分级低 1:40–1:80、中 1:160–1:320、高 ≥1:640 | 1:40 “determined from the literature”（Tan et al., Arthritis Rheum 1997 健康人 ANA 范围；Arch Pathol Lab Med 2000 临床应用指南） |
| ImmuGlo HEp-2 Elite (K172745) | **1:40** | 1+~4+ 且可识别核型；ICAP 命名（AC-1 均质、AC-2/4/5 颗粒、AC-3 着丝点、AC-6/7 核点、AC-8/9/10 核仁；胞浆 AC-15~23）；~90% 细胞为 PSIP1 敲除，DFS70 阳性仅在 ~10% 野生型细胞显色，用于区分 AC-2 与均质/颗粒 | “established from literature”；建议实验室自建正常值 |
| Zeus IFA ANA HEp-2 + dIFine (K201956) | **1:40** | 任何核型的苹果绿荧光，1+~4+；dIFine 识别 8 种核型（均质、颗粒、着丝点、核仁、核点、核膜、胞浆核糖体、胞浆线粒体）；软件可输出 Uncertain | “Please refer to K781516”；1:40 以下视为阴性 |
| AESKUSLIDES HEp-2-Gamma + HELIOS (K153117) | **1:40 或 1:80**（由实验室选择） | ≥1+ 且可辨核型；低/中/高滴度分级同上 | 厂家建议实验室按自身人群与仪器自定筛查稀释度 |

**固相法**

| 产品 | cutoff / 判读 | 建立方式 | 参考人群 |
|---|---|---|---|
| BioPlex 2200 ANA Screen (K041658) | 12 项半定量：Antibody Index (AI) ≥1.0 阳性（测量范围 0.2–8.0 AI）；dsDNA 定量 IU/mL：≤4 阴性、**5–9 indeterminate**、≥10 阳性；ANA Screen 定性 = 任一项阳性；SS-A 60/52 或 RNP-A/68 任一阳性即报 SS-A/RNP 阳性 | **AI 1.0 ≈ 非患病人群第 99 百分位**；cutoff 研究 719 例（285 例 predicate 筛查阴性正常人 + 434 例 predicate EIA 单项阳性），做 ROC，各项灵敏度 67–100%、特异度 94–100%（dsDNA 67%/99%、Chromatin 94/100、RNP 96/94、SS-B 100/98、SS-A 99/98、Scl-70 98/98、Sm 96/100、CENP 98/98、Sm/RNP 96/99、Rib-P 82/98、Jo-1 95/100）；dsDNA 溯源 WHO Wo/80 | 222 例美国正常献血者：ANA Screen 6.8% 阳性、dsDNA 1.4% 阳性 + 3.6% indeterminate |
| Euroimmun ANA Screen ELISA (K131185) | **ratio ≥1.0 阳性**（样本 OD / 校准品 OD）；无灰区 | 沿用 predicate Aeskulisa（ratio 1.0）；校准品逐批用 ≥8 阳性 + 2 阴性参考血清调校 | 200 例健康献血者：3.0% 阳性，ratio 均值 0.2 |
| Gold Standard ANA Screen ELISA (K131330) | ratio units：**<0.83 阴性、0.83–1.2 equivocal、>1.2 阳性** | **99 例正常献血者 + 3 例已知阴性的 mean + 3SD** 定 cutoff；equivocal = cutoff ±20%，目的“encourage repeat testing” | 99 例：2 阳性、7 equivocal、prevalence 7.0%；95th percentile 男 1.09 / 女 0.82 |
| Varelisa ReCombi ANA Screen (K083188) | 与 cut-off calibrator 比较，报 positive / equivocal / negative（具体比值：文件未载明） | 在 K050625 建立；本次用 460 例健康高加索献血者**确认** | 460 例（年龄性别均衡） |

**IFA 1:40 vs 1:80 的实测差异（K150155，同一 410 例队列）**：健康对照 150 例阳性率 1:40 为 27.3%、1:80 为 11.3%；SLE 85% vs 80%；SjS 76.7% vs 70%；SSc 66.7% vs 50%；RA 66.7% vs 56.7%。1:80 vs 1:40 PPA 79.9% (74.1–85.0)、NPA 97.3% (91.5–100)。Aesku K153117 同样报告 1:80 较 1:40 SLE 灵敏度下降（A 法 86.4%→78.1%）、非健康对照特异度上升（58.0%→68.3%）。

### 5. 生物学 / 免疫学依据
- 申报文件所载：IFA 以 HEp-2（人喉癌上皮）细胞为基质，检测 IgG 类抗核抗体；Euroimmun 用 HEp-20-10、Immco 用 PSIP1 敲除 HEp-2 消除 LEDGF/DFS70 抗原，文件引用 ICAP 2014–2015 共识并称 “due to the high prevalence of AC-2 patterns in an ANA screening population and relatively low association with autoimmune rheumatic diseases, ICAP recommends all clinical labs report the DFS70 pattern”；Euroimmun IFA 决策摘要引用 “As per the American College of Rheumatology a prevalence of ANAs in healthy individuals is about 3.0–15.0%”，NOVA Lite 摘要称健康人 “10–20% of positivity may be seen in reference subjects according to published literature”。固相筛查以有限抗原组合（dsDNA、组蛋白、ENA 等）替代全细胞，各文件均用 CDC/IUIS ANA 参考血清验证：Euroimmun ELISA 对 CDC #6（U3-RNP/fibrillarin）和 #11（PM-Scl）阴性，Euroimmun ENA pool 对 CDC #1（dsDNA）、#6、#8（centromere）、#11 阴性——即固相法漏检其抗原谱之外的特异性。
- 背景（非申报文件）：ACR 2009 年立场声明（*Methodology of testing for antinuclear antibodies*）认定 **IFA HEp-2 为 ANA 筛查的金标准/参考方法**，要求实验室报告方法学；2019 EULAR/ACR SLE 分类标准把 HEp-2 IFA ANA ≥1:80 作为入门标准。固相法（ELISA/CLIA/多重微珠）因抗原有限、对核仁型/DFS70 等敏感度不同，与 IFA 一致性有限，这正是 K041658 中 BioPlex vs EIA 筛查 PPA 仅 70% 的机理基础。

### 6. 样本类型与样本要求
- 所有 IFA 产品：**仅血清**（NOVA Lite、Euroimmun、Immco、Zeus、Aesku 均 “Serum”，Matrix comparison: Not applicable）。样本稳定性各 IFA 文件未载明（NOVA Lite 仅给试剂/结合物稳定性）。
- BioPlex 2200：血清、EDTA、钠肝素（214 例匹配，EDTA 回收 95.3–100.2%、钠肝素 90.5–103.6%）；锂肝素在再现性研究中验证。
- Euroimmun ANA Screen ELISA：血清、EDTA/Li-heparin/枸橼酸血浆（各 12 对，Passing-Bablok 斜率 0.99–1.00）。
- Gold Standard：仅血清。
- Varelisa ReCombi：血清、枸橼酸、EDTA 血浆；**肝素血浆干扰，说明书警告不得使用**；样本稳定性依 CLSI H18-A3。

### 7. 分析性能验证所依据的标准
| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K150155 NOVA Lite DAPI | FDA ANA 指南（2009）；C28-A3（LIU cutoff 参考区间）；EP07-A2（干扰：胆红素 10 mg/dL、Hb 200、TG 1000、胆固醇 224、RF IgM 56 AU）；EP09-A2IR（基质比对，实际未做） | 3 项重复性（LIU 与读片一致性）、3 中心×2 操作者再现性、3 批次 |
| K131791 Euroimmun IFA | FDA ANA 指南（2009）；DIN EN 13640:2002（稳定性） | 29 样本 10 次批内/20 次批间/3 批/双观察者（±1 强度级）；半定量 1:40→1:10240 |
| K172745 Immco | EP05-A3（精密度，120 复/样本）、EP06-A（稀释线性）、EP07-A2（干扰含 naproxen、ibuprofen、prednisone、MMF、HCQ、CTX、rituximab、belimumab）、**EP12-A2（定性比对）** | 3 中心再现性（核型/滴度 ±1）、3 批次、3 读片者研究 |
| K201956 Zeus | EP05-A3、EP06（2nd）、EP07-A2、EP17-A2、**EP28-A3c**、IEC 61010-1 | 10 天精密度、5 天 3 中心再现性、3 批次；23 份 ICAP 参考样本；干扰 |
| K153117 Aesku/HELIOS | EP07-A2、EP17-A2、EP25-A、EP28-A3c、EP06-A、ISO 14971 | 1:40 与 1:80 双稀释度精密度/再现性/批间/仪器间；携带污染 |
| K041658 BioPlex | NCCLS EP5-A（3 中心精密度）、EP7-A（干扰） | 交叉反应（心磷脂、AMA-M2、ssDNA、组蛋白、MPO、PR3、gliadin、SMA、TPO、Tg、RF 各 10 例）；dsDNA 线性（WHO Wo/80） |
| K131185 Euroimmun ELISA | FDA ANA 指南（2009） | 20 复孔批内、10 次批间、3 批；干扰 Hb 1000、TG 2000、胆红素 40 mg/dL、RF 500 IU/mL |
| K131330 Gold Standard | EP07-A2；FDA ANA 指南 | 重复性、3 中心再现性、3 批；CDC + AMLI 参考血清；HAMA 干扰 |
| K083188 Varelisa | H18-A3（样本处理） | 5 次×8 复精密度；CDC/AMLI 参考血清 |

### 8. 临床验证设计与结果

**IFA**

| 产品 | 目标疾病组（分类标准） | 对照组 | 灵敏度 (95% CI) | 特异度 (95% CI) | 与 predicate/参考法 |
|---|---|---|---|---|---|
| NOVA Lite DAPI 1:80 (K150155) | 方法比对队列：SLE 100、SjS 30、SSc 30、AIM 10、MCTD 20（“clinically characterized”，标准未载明）；临床队列 463 例：SLE 75、SSc 20、SS 20、AIL 20、MCTD 21、AIM 26、DIL 25 等 | 健康 150（方法比对）；临床队列健康 75 + HBV 20、HCV 5、HIV 5、梅毒 5、RA 20、纤维肌痛 25、ANCA 26、IBD 20、AITD 24、乳糜泻 24、其他 7 | 方法比对：SLE **80.8%** (71.7–88.0)，SARD **71.4%** (64.4–77.8)（1:40 predicate 85.9%/78.8%）；3 中心手工读片 SLE 72.0/70.7/82.7%，CTD+AIL 62.9/65.6/71.0% | 方法比对（RA+感染 n=60）65.0% (51.6–76.9)（1:40 为 56.7%）；3 中心手工（非健康对照 n=174）74.1/67.2/67.2% | 1:80 vs 1:40：PPA 79.9%、NPA 97.3%；含/不含 DAPI 结合物（同 1:80，n=407）PPA 98.6%、NPA 94.3%；NOVA View vs 手工总一致 87.0–89.8%，核型一致 72.7–86.3% |
| Euroimmun IFA 40 (K131791) | 无临床灵敏度/特异度研究（“Not applicable”） | — | — | — | vs HEp-2000 ANA-Ro（n=200 常规送检）PPA 92.4% (83.2–97.5)、NPA 95.5% (90.5–98.3)；第二队列 n=156 PPA/NPA 100%；核型一致 91.4%；200 健康者 3.5% >1:40；138 例常规筛查 11.6% 阳性 |
| ImmuGlo HEp-2 Elite (K172745) | SARD 256：SLE 78、SjS 65、SSc 27、肌炎 36、MCTD 10、DIL 12、APS 15、AIH 13（标准未载明） | 497：非 SARD 自身免疫 254（RA 87、PBC 25、Crohn 25、UC 23、ANCA 20、AITD 20、恶性贫血 20、银屑病 20、乳糜泻 10、Raynaud 4）、炎症 13、肿瘤 57、感染 173 | **78.1%** (72.7–82.8)；SLE 94.9%、SjS 95.3%、SSc 63.0% | **75.1%** (71.1–78.7)；RA 59.8%、PBC 68% 阳性；感染 16.8% | vs predicate ImmuGlo ANA（n=591）PPA 99.7% (98.2–99.9)、NPA 98.3% (96.1–99.3)；DFS70 富集队列 121 例 PPA 100%；128 例正常 5.5% 阳性 |
| Zeus IFA + dIFine (K201956) | 190 例：SLE 40、SS 30、SSc 20、AM 20、MCTD 20、CREST 20、AIH 20、DIL 20（商业来源，标准未载明） | 190（中心 2 为 160）：乳糜泻 22、ANCA 血管炎 28、Crohn 10、RA 30、IBD 10、AITD 30、UC 10、肿瘤 20、纤维肌痛 10、感染 10 | SLE **45.0–57.5%**；CTD+ANA 相关病 **51.6–60.0%**（3 中心×2 技师×3 方法） | **67.4–83.1%** | 手工 A vs 数字图像 B（2268 次读片）PPA 96.4% (95.0–97.5)、NPA 97.0% (96.0–97.8)；180 例健康者 A/B 10.6% 阳性，dIFine 7.8% 阳性 + 7.8% uncertain |
| AESKUSLIDES HEp-2-Gamma/HELIOS (K153117) | 中心 A/B 460 例：SLE 90、AIH 18、CREST 4、MCTD 8、肌炎 35、SSc 39、SjS 64、UCTD 6（标准未载明） | APS 5、乳糜泻 20、HBV 30、HCV 28、RA 60、SpA 35、血管炎 18 | 1:40 合并：SLE **86.4%** (76.1–92.7)，CTD+AIL 84.8% (78.8–89.3)；1:80：78.1%/78.3% | 1:40 **58.0%** (50.6–65.1)；1:80 68.3% (61.1–74.8) | 手工 D vs predicate E 总一致 ≥93%（1:40）；软件 A vs E 核型一致 92.6–97.5% |

**固相法**

| 产品 | 目标疾病组 | 对照组 | 灵敏度 (95% CI) | 特异度 (95% CI) | 与 predicate / IFA |
|---|---|---|---|---|---|
| BioPlex 2200 ANA Screen (K041658) | 908 例风湿科前瞻样本；**Criteria-diagnosed targeted CTD n=408**（ACR 或文献分类标准：SLE、SSc、Raynaud、CREST、SjS、MCTD、UCTD、DM/PM）；physician-diagnosed n=413 | 500 例无目标 CTD 诊断 | **67.9%** (63.2–72.6) | **77.4%** (73.6–81.2) | vs Helix EIA ANA Screen（n=908）PPA 70.2% (66.1–74.4)、NPA 90.4% (87.5–93.4)；criteria-CTD 内 PPA 80.0%、NPA 87.7%；单项 vs EIA PPA 58.6（chromatin）–98.5%（CENP）、NPA 91.6–100%；**与 IFA 的一致性：文件未载明** |
| Euroimmun ANA Screen ELISA (K131185) | 429 例：MCTD 21、SLE 213、PM/DM 26、SSc 81、SjS 88（标准未载明） | 309：乳糜泻 10、Wegener 17、RA 203、其他自身免疫 63、感染 16 | **72.5%** (68.0–76.7)：MCTD 95.2%、SLE 73.2%、SjS 81.8%、SSc 72.8%、PM/DM 15.4% | **95.8%** (92.9–97.7)；RA 94.1% | vs Aeskulisa ANA HEp-2（n=290）PPA 96.5% (92.0–98.8)、NPA 98.0% (94.2–99.6)；**与 IFA：文件未载明** |
| Gold Standard ANA Screen ELISA (K131330) | 510 例 CTD：SLE 322、SSc 40、PM 12、DM 15、PM/DM 重叠 5、肌炎 10、MCTD 28、UCTD 3、SjS 75（商业血清库，标准未载明） | 201 非 CTD（RA 100、OA 20、PBC 8、AIH 2、桥本 17、Graves 17、UC 5、乳糜泻 5、PAPS 22、GPA 5）+ 42 感染/RF | equivocal 计阳性 **83.5%** (79.1–86.0)；计阴性 79.2% (75.4–82.7)；SLE 82.9%、MCTD 100% | **96.7%** (93.6–98.6)；RA 2% 阳性 | vs Aeskulisa（n=848 常规送检）equivocal 计阳性 PPA 94.9% (91.5–97.2)、NPA 92.2% (89.7–94.2)；55 例单项阳性样本同测 HEp-2 IFA：dsDNA/SS-A/SS-B/Sm/Sm-RNP/Scl-70/centromere 组 HEp-2 100% 阳性，Jo-1 组 20%、Ro52 与 Rib-P 组 80% |
| Varelisa ReCombi (K083188) | 无临床灵敏度/特异度（“Not applicable”） | — | — | — | 新旧版本比对（253 例：感染 50、非 CTD 50、CTD 103、RA 50 + 92 献血者）PPA 99.0% (94.7–100)、NPA 96.7% (92.4–98.9) |

### 9. 厂家间差异与要点
1. **筛查稀释度**：多数 IFA 沿用文献/predicate 的 1:40；Inova 改为 1:80 并用同一队列证明健康人阳性率由 27.3% 降至 11.3%、SLE 灵敏度由 85.9% 降至 80.8%；Aesku 同时申报 1:40/1:80 由实验室选择。
2. **自动读片系统的“cutoff”是荧光强度阈值**（NOVA View LIU 48 = 健康人 90th percentile；dIFine 输出 uncertain 类），临床判读仍以稀释度 + 核型为准并须人工确认。
3. **固相法 cutoff** 均为比值/指数：BioPlex AI 1.0（≈99th percentile，配合 predicate-EIA-阳性样本 ROC）；ELISA 用 cutoff 校准品 ratio 1.0，Gold Standard 独有 mean+3SD 及 ±20% 灰区。
4. **临床性能**：IFA 灵敏度（SLE 72–95%）通常高于固相筛查（SLE 73–83%），但 IFA 在非 SARD 自身免疫/感染对照中阳性率高（特异度 58–83%），固相法特异度 95–97%（对照多为 RA/器官特异性自身免疫/感染）。
5. **与 IFA 的 PPA/NPA**：本次所取固相法文件均以另一 ELISA/EIA 为 predicate，未系统报告对 IFA 的一致性（Gold Standard 仅 55 例小样本）。
6. **分类标准**：仅 BioPlex 明确按 ACR/文献分类标准划分 “criteria-diagnosed targeted CTD”；其余多为 “clinically characterized”，具体标准文件未载明。
7. **样本类型**：IFA 一律血清；固相法部分允许 EDTA/肝素/枸橼酸血浆，但 Varelisa 明确肝素血浆干扰。

### 10. 来源
- K150155：https://www.accessdata.fda.gov/cdrh_docs/reviews/K150155.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K150155.pdf
- K131791：https://www.accessdata.fda.gov/cdrh_docs/reviews/K131791.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K131791.pdf
- K172745：https://www.accessdata.fda.gov/cdrh_docs/reviews/K172745.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K172745.pdf
- K201956：https://www.accessdata.fda.gov/cdrh_docs/reviews/K201956.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K201956.pdf
- K153117：https://www.accessdata.fda.gov/cdrh_docs/reviews/K153117.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153117.pdf
- K041658：https://www.accessdata.fda.gov/cdrh_docs/reviews/K041658.pdf ；K113610：https://www.accessdata.fda.gov/cdrh_docs/reviews/K113610.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K113610.pdf
- K131185：https://www.accessdata.fda.gov/cdrh_docs/reviews/K131185.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K131185.pdf
- K131330：https://www.accessdata.fda.gov/cdrh_docs/reviews/K131330.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K131330.pdf
- K083188：https://www.accessdata.fda.gov/cdrh_docs/reviews/K083188.pdf

---

## 可提取核抗原抗体与抗 dsDNA（ENA panel：dsDNA、Sm、RNP、SS-A/Ro、SS-B/La、Scl-70、Jo-1、centromere）

### 1. 法规定位
- 21 CFR **866.5100**，**Class II**，Panel: Immunology (82)。
- Product codes 与 openFDA 累计清关数：
  - **LLL** Extractable Antinuclear Antibody, Antigen and Control：**total=182**（ENA 单项/组合、CTD screen 多在此码下）。
  - **LSW** Anti-DNA Antibody, Antigen and Control：**total=8**（QUANTA Flash dsDNA、EliA dsDNA、Aeskulisa/AESKUSLIDES nDNA、FIDIS dsDNA）。
  - **LRM** Anti-DNA Antibody (Enzyme-Labeled), Antigen, Control：**total=29**（dsDNA ELISA：Euroimmun NcX、FARRZYME、Immco 等）。
  - **KTL**（Crithidia luciliae IFA / anti-nDNA）：**total=19**（NOVA Lite DAPI dsDNA CLIFT、Euroimmun CLIFT EUROPattern、Zeus/Immuno Concepts nDNA）。
  - 其他相关：LKO（Anti-RNP）、LKJ（EliA Ro52/Ro60 用此码）、MQA（Anti-Ribosomal P）、OBE/PET（IDS 产品同时列出）。
- 提示：openFDA 仅按主代码索引，多重/多项产品（如 BioPlex ANA Screen 列 LKJ/LRM/MQA/LKO/LJM/LLL）只计入一个代码。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K152013 | 2016-04-11 | Inova / QUANTA Flash dsDNA（BIO-FLASH） | CIA，**定量** IU/mL（溯源 WHO Wo/80） | 血清 | IgG anti-dsDNA；aid in the diagnosis of SLE | 决策摘要 + summary |
| K072393 | 2007-12-07 | Phadia / EliA dsDNA（ImmunoCAP 100/250）+ EliA ANA Control | FEIA，定量 IU/mL（Wo/80） | 血清、肝素/EDTA/枸橼酸血浆 | aid in the diagnosis of SLE | 决策摘要 |
| K083381 | 2009-04-15 | Euroimmun / Anti-dsDNA-NcX ELISA (IgG)（dsDNA-核小体复合物抗原） | ELISA，定量/半定量 IU/mL（Wo/80） | 血清、EDTA/枸橼酸血浆 | 同上 | 决策摘要 |
| K062183 | 2006-11-21 | The Binding Site / FARRZYME Human High Avidity anti-dsDNA EIA | ELISA（高严格洗涤，仅测高亲合力抗体），定量 IU/mL | 血清 | 同上（“high avidity IgG”） | 决策摘要 |
| K192916 | 2020-12-11 | Inova / NOVA Lite DAPI dsDNA Crithidia luciliae Kit（手工或 NOVA View） | **CLIFT** IFA，定性/半定量 | 血清 | 同上 | 决策摘要 + summary |
| K213403 | 2023-09-29 | Inova / Aptiva CTD Essential Reagent（dsDNA、RNP、Sm、Ro52、Ro60、SS-B、Scl-70、Jo-1、Centromere、Ribo-P 十项） | 粒子多重荧光免疫分析（PMAT），FLU 半定量 + dsDNA 定量 | 血清 | 各项分别 aid in diagnosis of SLE / MCTD / SjS / SSc / IIM | 决策摘要（46 页）+ summary |
| K253367 | 2026-06-25 | Phadia / EliA CTD 13 Screen（Phadia 250；13 抗原含 RNA Pol III、Rib-P、dsDNA） | FEIA，定性 ratio | 血清 | aid in the diagnosis of SLE、SS、SSc、PM、DM、MCTD | 决策摘要 + summary |
| K190710 | 2019-11-29 | Phadia / EliA SymphonyS（Phadia 250 与 2500/5000） | FEIA，定性 ratio | 血清、Li-heparin/EDTA 血浆 | 同上（无 dsDNA） | 决策摘要 + summary |
| K182353 | 2018-11-27 | Phadia / EliA CENP、EliA U1RNP、EliA RNP70（迁移至 Phadia 2500/5000） | FEIA，半定量 EliA U/mL | 血清、Li-heparin/EDTA 血浆 | CENP→scleroderma (CREST)；U1RNP/RNP70→MCTD 与 SLE | 决策摘要 |
| K210902 | 2022-07-27 | Phadia / EliA Ro52、EliA Ro60 | FEIA，半定量 | 血清 | Ro52→SLE、SS、IIM、SSc；Ro60→SLE、SS | 决策摘要 + summary |
| K152635 | 2016-06-01 | Inova / QUANTA Flash Scl-70 | CIA，半定量 CU | 血清 | aid in the diagnosis of systemic sclerosis | 决策摘要 + summary |
| K151429 | 2016-02-12 | Inova / QUANTA Flash Jo-1 | CIA，半定量 CU | 血清 | aid in the diagnosis of IIM | 决策摘要 + summary |
| K141328 | 2015-02-12 | Inova / QUANTA Flash Ro60 | CIA，半定量 CU | 血清 | aid in the diagnosis of SLE 与 SjS | 决策摘要 + summary |
| K141210 | 2015-01-29 | Inova / QUANTA Flash SS-B | CIA，半定量 CU | 血清 | 同上 | 决策摘要 + summary |
| K112996 | 2013-04-09 | Euroimmun / Anti-ENA Pool ELISA (IgG)（nRNP/Sm、Sm、SS-A 60/Ro-52、SS-B、Scl-70、Rib-P） | ELISA，定性 ratio | 血清 | aid in the diagnosis of MCTD、SLE、SjS、PSS | 决策摘要 + summary |
| K250666 | 2025-10-22 | Zeus / Alegria Flash CTD Screen（SSA-60、SSA-52、SS-B、Jo-1、Scl-70、SmRNP、Sm、dsDNA、Rib-P、核小体、CENP-B） | 磁微粒 CLIA，定性 CLIA Units | 血清 | ANA screen（aid in the diagnosis of ANA-associated CTD） | 决策摘要 + summary（Statement） |
| K253817 | 2026-08-21 | Immunodiagnostic Systems / IDS SS-A/Ro、IDS SS-A/Ro 60 kDa、IDS SS-B/La | CLIA，半定量 AU/mL | 血清 | SS-A/Ro→SLE、SjS、SSc；Ro60/SS-B→SLE、SjS | 决策摘要 + summary |
| K041658 | 2004-12-20 | Bio-Rad / BioPlex 2200 ANA Screen（各单项，见 ANA 节） | 多重微珠 | 血清、EDTA/肝素血浆 | — | 决策摘要 |

### 3. 预期用途与声明类型
- 全部为 **Rx only**、**辅助诊断**。单项抗体声明与疾病绑定：dsDNA→SLE；Sm→SLE；RNP/U1RNP/RNP70→MCTD 与 SLE；Ro60/SS-B→SLE 与 SjS；Ro52→SLE、SjS、IIM、SSc；Scl-70→SSc；centromere/CENP-B→SSc（CREST）；Jo-1→IIM/肌炎；Rib-P→SLE；RNA Pol III→弥漫型 SSc（EliA RNA Pol III K202541，为 CTD 13 的 predicate）。
- 测量类型：dsDNA 多为 **quantitative（IU/mL，溯源 WHO Wo/80）**；其余 ENA 单项 **semi-quantitative**（CU、FLU、EliA U/mL、AU/mL）；筛查组合为 **qualitative**（ratio/CLIA Units）。
- 组合筛查（EliA CTD 13、SymphonyS、Alegria、Euroimmun ENA pool、BioPlex）不报告单项特异性（BioPlex 除外，其同时输出 11 项结果）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**anti-dsDNA**

| 产品 | cutoff | 建立方式 | 参考人群 |
|---|---|---|---|
| QUANTA Flash dsDNA (K152013) | **<27 IU/mL 阴性；27–35 indeterminate；>35 阳性** | 参考人群 **95th percentile = 27**（cutoff）、**99th percentile = 35**（灰区上限） | 170 例（121 健康献血者 + 49 疾病对照）；期望值另测 300 健康者：2.0% 灰区、1.3% 阳性 |
| Aptiva dsDNA (K213403) | 同上（27–35 IU/mL，“Same” as QUANTA Flash） | 120 例疾病对照初定 + 预期用途人群验证 | 乳糜泻 16、桥本 18、感染 25、PBC 7、PBC/AIH 2、PSC 8、PSC/AIH 1、RA 30、Lyme 13 |
| EliA dsDNA (K072393) | <10 阴性；**10–15 IU/mL equivocal**；>15 阳性 | 预定义，用健康人群分布确认（95th percentile 6.3、99th 17.3 IU/mL） | 400 例健康高加索献血者 |
| Euroimmun NcX ELISA (K083381) | **100 IU/mL**（predicate FARRZYME 为 30） | 预定义，健康人群确认（均值 6.8、范围 0.6–72 IU/mL） | 400 例健康献血者 |
| FARRZYME (K062183) | **>30 IU/mL 阳性** | 150 例健康献血者全部 <17 IU/mL（97% <12.3） | 150 例 |
| BioPlex dsDNA (K041658) | ≤4 阴性、5–9 indeterminate、≥10 IU/mL 阳性 | 与其他项一起在 719 例 predicate 定性样本上做 ROC（灵敏度 67%、特异度 99%） | 222 例正常：1.4% 阳性、3.6% 灰区 |
| NOVA Lite DAPI dsDNA CLIFT (K192916) | 筛查稀释 **1:10**；阳性 = 动基体（kinetoplast）±核着色强于阴性对照；NOVA View 报 LIU 及 negative/positive/indeterminate，SWT 软件估算终点滴度 | 与 predicate 同（1:10）；建议实验室自建滴定方案 | 120 例健康者：手工 3.3%、NOVA View 9.2%、数字读片 0.8% 阳性 |

**ENA 单项 / 组合**

| 产品 | cutoff | 建立方式 | 参考人群 |
|---|---|---|---|
| QUANTA Flash Scl-70 (K152635) | **≥20 CU** | 254 例（疾病对照 + 健康）按 C28-A 取 **99th percentile** → 7,387 RLU（=20 CU）；再用 19 例 predicate 阳性 SSc 样本调整，**cutoff 上调至 15,000 RLU 并重新赋值为 20 CU** | 期望值 100 例健康者 95th percentile 1.4 CU |
| QUANTA Flash Jo-1 (K151429) | **≥20 CU（10,000 RLU）** | 207 例疾病对照（SSc 31、SLE 30、Crohn 21、MS 19、HCV 19、UC 18、PsA 13、梅毒 10、健康 10、PMR 9、RA 8、SpA 5、SjS 3、其他 11）**99th percentile** 9,151 RLU → 取整 10,000 RLU（1 例 SSc 阳性） | 400 例健康者 0% 阳性 |
| QUANTA Flash Ro60 (K141328) | **≥20 CU** | 156 例（健康 115、病毒性肝炎 8、梅毒 5、RA 23、HIV 5）非参数 **99th percentile** 8,193.6 RLU；结合 12 份 CAP/UKNEQAS 能力验证样本上调至 **12,000 RLU = 20 CU** | 98 例健康者 2% 阳性 |
| QUANTA Flash SS-B (K141210) | **≥20 CU** | 99th percentile 4,097 RLU（1 例参考样本阳性）；用 32 例 IIF/FIA/ELISA 确认的 SS-B 阳性样本（37,743–942,882 RLU）上调至 **12,000 RLU = 20 CU**（510(k) summary） | 138 例健康者 0.7% 阳性 |
| Aptiva CTD Essential 其余 9 项 (K213403) | **≥5.00 FLU 阳性** | 同上 120 例疾病对照初定 + 验证 | 115 例健康者：RNP 1.7%、Ro60 1.7%、Ro52/Scl-70/centromere 各 0.9% 阳性 |
| EliA CENP / RNP70 (K182353) | <7 阴性、**7–10 equivocal**、>10 EliA U/mL 阳性；EliA U1RNP：<5、**5–10**、>10 | 沿用 predicate K082759/K083117 | 400 例健康者：99th percentile CENP 1.3、U1RNP 3.0、RNP70 2.5 U/mL，无一进入灰区 |
| EliA Ro52 / Ro60 (K210902) | <7 阴性、**7–10 equivocal**、>10 阳性；灰区建议 8–12 周后复测 | Ro52：69 健康 + 19 SLE + 9 SS，再用 10 IIM + 14 SSc 验证；Ro60：70 健康 + 22 SLE + 6 SS | 同左 |
| EliA CTD 13 Screen (K253367) | ratio **<0.7 阴性、0.7–1.0 equivocal、>1.0 阳性** | 70 健康献血者 + 4 MCTD + 8 IIM + 7 SSc + 6 SLE + 5 SjS | 330 例健康者：95th percentile 0.5、97th 0.8 |
| EliA SymphonyS (K190710) | 同上 ratio 三分段 | 70 健康 + 30 例已知单项阳性 SARD；健康 95th percentile 低于灰区 | 558 例健康者：99th percentile 1.1；1.3% 阳性、0.5% 灰区 |
| Euroimmun Anti-ENA Pool (K112996) | ratio **≥1.0 阳性** | 沿用 predicate Aeskulisa | 200 例健康者 1.5% 阳性 |
| Alegria Flash CTD Screen (K250666) | **100 CLIA Units** | 用 28 例 ANA 阳性 + 25 例健康者**验证**（建立方法未载明） | 200 例健康者 5.0% 阳性 |
| IDS SS-A/Ro、Ro60、SS-B/La (K253817) | **≥10.0 AU/mL 阳性**（决策摘要一处写 “100 AU/mL” 系笔误，判读表为 10.0） | 建立方法：决策摘要未载明；summary 建议实验室确认 | 112 例健康者 0.9% 阳性（SS-A/Ro） |
| BioPlex 各项 (K041658) | AI ≥1.0（≈99th percentile） | 见 ANA 节 | — |

### 5. 生物学 / 免疫学依据
- 申报文件所载：
  - **dsDNA 亲合力/抗原形式**决定不同方法一致性：FARRZYME 以更严格洗涤只检测高亲合力抗体，故与常规 ELISA（PPA 58.8%）和 CLIFT（PPA 67.3%）一致性低，但与 Farr RIA 高（PPA 89.8%）；Bio-Rad 解释 dsDNA 方法差异源于“捕获抗原类型、是否同时检出 ssDNA 抗体、对高/低亲合力抗体的捕获能力”；Phadia 指出 AMLI-J（低亲合力 dsDNA）在 EliA 阴性可能因严格洗涤；Euroimmun 用 dsDNA–核小体复合物（NcX）抗原。
  - 参考血清：各产品用 CDC/IUIS ANA 参考血清与 AMLI 共识盘验证抗原特异性（如 QUANTA Flash dsDNA 对 CDC #1 50.7 IU/mL 阳性、#5 (Sm) 22.7、#7 (SS-A) 14.8 IU/mL 阴性）。
  - Euroimmun ENA pool 发现 CDC #10（Jo-1）阳性系样本共存 Ro-52 抗体所致；Aptiva 注明 MCTD 样本可因疾病关联而 Ro52 阳性。
- 背景（非申报文件）：抗 dsDNA 与抗 Sm 为 SLE 分类标准（ACR 1997、SLICC 2012、EULAR/ACR 2019）免疫学项，抗 dsDNA 滴度随疾病活动/狼疮肾炎波动；Farr 法（放射免疫沉淀）偏向高亲合力抗体、CLIFT 特异性高、ELISA 灵敏度高但可检出低亲合力抗体——这是 FDA 文件中“以 Farr/CLIFT 作比较法”并出现低 PPA 的免疫学背景。Ro52（TRIM21）与 Ro60（hY-RNA 结合蛋白）为不同抗原，Ro52 与肌炎/SSc 相关性更强，因此新一代平台分别检测。

### 6. 样本类型与样本要求
- QUANTA Flash dsDNA/Scl-70/Jo-1/Ro60/SS-B、Aptiva、EliA Ro52/Ro60、EliA CTD 13、IDS、Euroimmun ENA pool：**仅血清**。QUANTA Flash dsDNA 样本稳定性 RT 48 h、2–8 °C 10 d、3 次冻融；Scl-70/Jo-1 RT 48 h、2–8 °C 21 d；EliA CTD 13 支持冷冻 24.5 个月；IDS 2 次冻融。
- EliA dsDNA：血清、肝素/EDTA/枸橼酸血浆（46 组，斜率 1.015–1.067）。
- Euroimmun NcX：血清、EDTA/枸橼酸血浆（35 组）；**肝素血浆不适用**。
- FARRZYME、CLIFT：血清（CLIFT 样本 RT 48 h、2–8 °C 7 d、3 次冻融）。
- EliA SymphonyS、EliA CENP/U1RNP/RNP70（Phadia 2500/5000）：血清、Li-heparin、EDTA（**枸橼酸血浆在新平台被取消**）；SymphonyS 62 组匹配无阴阳转换。
- BioPlex：血清、EDTA、钠肝素。
- Alegria Flash CTD Screen：仅血清（“Sample matrix: Serum — Same” as predicate）；样本稳定性文件未载明。

### 7. 分析性能验证所依据的标准
| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K152013 QF dsDNA | EP05-A2、EP6-A、EP07-A2、EP09-A3、EP17-A2、C28-A3 | 精密度（9 样本 20–21 天）、3 中心再现性、线性 9.8–666.9 IU/mL、LoB/LoD/LoQ（RLU）、干扰含 RF IgM 947 IU/mL、CDC 参考血清、cutoff 百分位 |
| K152635 / K151429 / K141328 / K141210 QF ENA | EP05-A2、EP6-A、EP07-A2、EP09-A3、EP17-A2（Ro60/SS-B 为 EP17-A）、C28-A3（Ro60 写 C28-A3c；Jo-1 summary 写 EP28-A3c）、EP09-A2IR（基质）、FDA 校准品 abbreviated 510(k) 指南 | 同上；Scl-70 干扰含 IgG 70 mg/mL、prednisone、naproxen |
| K213403 Aptiva | EP05-A3、EP06-Ed2、EP07-A2、EP17-A2、EP28-A3c | 十项精密度/线性/LoB-LoD-LoQ；干扰（rituximab 未评估）；CDC 参考盘 |
| K253367 EliA CTD 13 | EP07（3rd）、**EP12-Ed3（定性二元性能）**、EP17-A2、EP28-A3c、EP34（扩展测量区间）、EP37（干扰补充表）、GP44-A4（样本处理） | 按 EP12 设计的定性精密度（20 天、3 中心、3 批）；23 种内外源干扰物（belimumab、rituximab、infliximab、MTX、MMF、HCQ 等） |
| K190710 EliA SymphonyS | EP05-A3、EP7-Ed3、EP17-A、H18-A4 | 精密度、干扰、检测限；基质比对 |
| K182353 EliA CENP/U1RNP/RNP70 | EP05-A3、EP6-A、EP17-A | 3 仪器×7 天精密度；LoB/LoD/LoQ；Passing-Bablok 仪器比对 |
| K210902 EliA Ro52/Ro60 | EP05-A3、EP06-Ed2、EP07（3rd）、EP09c（3rd）、EP17-A2、EP28-A3c、EP25-A、EP37 | 双仪器（Phadia 250 / 2500E）LoB/LoD、线性、干扰、稳定性 |
| K072393 EliA dsDNA | None referenced | 3 仪器精密度；灰区样本（~15 IU/mL）精密度；干扰（RF IgM 500 IU/mL）；CDC/AMLI |
| K083381 Euroimmun NcX | None referenced（LoD/LoQ 按 ICH Q2B；样本稳定性按 CLSI H18-A3） | 13 样本精密度；线性 10–800 IU/mL；CDC/AMLI；Scl-70 与 CCP 阳性血清交叉反应阴性 |
| K062183 FARRZYME | None referenced | 精密度、线性、Kokusai 干扰盘、6 例 IgG 骨髓瘤阴性 |
| K192916 CLIFT | **EP12-A2**、EP07-A2 | 手工/数字/NOVA View 三种读法精密度；1:10→1:5120 滴度线性；内源 + 9 种药物（rituximab、belimumab、MTX、AZA、MMF、HCQ、CTX 等）干扰；CDC 参考盘 |
| K112996 Euroimmun ENA pool | FDA ANA 指南（2009）；DIN EN 13640:2002 | 批内/批间/3 批（100% 判读一致）；干扰 RF 500 IU/mL；6 种单项 ELISA 比对 |
| K250666 Alegria | EP07-Ed3、**EP12-Ed3**、EP17-A2、EP25-Ed2、EP28-A3c、GP44-A4 | 定性精密度、稳定性 |
| K253817 IDS | EP05（3rd）、EP06（2nd）、EP07（3rd）、EP17（2nd）、EP25-A2、EP28-A3c、EP34、EP37 | 精密度、线性、LoB 0.0/LoD 1.2/LoQ 5.3 AU/mL、干扰 ≤15% |

### 8. 临床验证设计与结果

**anti-dsDNA（含 Farr / CLIFT 比较）**

| 产品 | SLE 组（标准） | 对照组 | 灵敏度 (95% CI) | 特异度 (95% CI) | 方法比对 |
|---|---|---|---|---|---|
| QUANTA Flash dsDNA (K152013) | 644 SLE（标准未载明）；另 22 例药物性狼疮排除 | 465：pAPS 20、SjS 50、乳糜泻 20、SSc 40、IIM 20、MCTD 20、Crohn 20、Graves 21、桥本 61、RA 101、血管炎 34、感染 58 | 灰区计阴性 **42.7%** (38.9–46.6)；计阳性 49.1% (45.2–52.9) | 计阴性 **94.4%** (91.9–96.3)；计阳性 87.5% (84.2–90.2)；对照阳性率 5.7%（MCTD 15%、pAPS 15%、RA 7.9%） | vs QUANTA Lite dsDNA SC ELISA（n=481）：灰区计阴性 PPA 80.3% (74.6–84.9)、NPA 79.1% (73.6–83.6)；计阳性 68.5%/70.0% |
| Aptiva dsDNA (K213403) | 230 SLE | 1039（SjS 141、SSc 217、MCTD 91、IIM 200 及 26 种其他） | 计阴性 **46.1%** (39.8–52.5)；计阳性 50.4% (44.0–56.8) | 计阴性 **92.9%** (91.2–94.3)；计阳性 90.1% | vs QUANTA Flash dsDNA（n=428）：计阴性 PPA 78.9% (70.8–85.1)、NPA 91.8% (88.2–94.4) |
| EliA dsDNA (K072393) | 无临床灵敏度/特异度（“None provided”） | — | — | — | **vs DPC Anti-DNA RIA**（n=97，60 SLE + 46 非 SLE，9 例 equivocal 剔除）：PPA 90.6% (79.3–96.9)、NPA 100% (92–100)、总一致 94.8%；ImmunoCAP 100 vs 250 r=0.986 |
| Euroimmun NcX (K083381) | 213 SLE（573 例冻存血清） | 360：RA 165、SjS 88、PSS 81、肌炎 26 | **59.6%** (52.7–66.3) | **97.2%** (95.0–98.7)；RA 4.2% 阳性 | vs FARRZYME（n=214：93 SLE、54 RA、26 ANA+ CTD、41 健康）：PPA 89.3% (71.8–97.7)、**NPA 81.2%** (74.8–86.5)（NcX 阳性 35 例 predicate 阴性） |
| FARRZYME (K062183) | 224 SLE（未选活动度） | 150 健康献血者；疾病对照 UC 13、Crohn 16、白喉毒素送检 9 及 Scl-70/Jo-1/PR3/MPO/GBM/gliadin/TPO/ASCA 阳性血清均阴性 | SLE **36%** 阳性（文献 28–63%） | 健康 100%（全部 <17 IU/mL） | **vs Farr RIA**（n=166）：PPA 89.8% (44/49)、NPA 94.0% (110/117)、总 92.8%；**vs Crithidia luciliae IFA**（n=252）：PPA 67.3% (37/55)、NPA 92.9%；vs BINDAZYME ELISA：PPA 58.8%、NPA 97.7%（文件解释：ELISA/CLIFT 检出低中高亲合力，FARRZYME 仅高亲合力） |
| NOVA Lite DAPI dsDNA CLIFT (K192916) | 391 SLE（766 例队列）；另 3 中心 100 SLE | 375：DIL 20、感染 60、血管炎 30、pAPS 20、SjS 30、乳糜泻 20、SSc 30、IIM 20、MCTD 20、Crohn 20、Graves 20、桥本 30、RA 35、AIH 20 | 手工/数字 **48.1%** (43.2–53.0)；NOVA View 57.0%；3 中心手工 32.7% (27.6–38.2)、NOVA View 40.0% | 手工 **91.2%** (87.9–93.7)、数字 92.3%、NOVA View 88.8%；3 中心手工 95.5%、NOVA View 85.4%；SSc 组 NOVA View 特异度仅 55–70% | vs predicate NOVA Lite dsDNA 手工（n=744）：PPA 87.8% (82.7–91.5)、NPA 96.0% (94.0–97.4)；grade ±2 级一致 99.6% |
| BioPlex dsDNA (K041658) | 前瞻 907 例 | — | ROC vs predicate 定性：67%/99% | — | vs Varelisa dsDNA EIA 定量相关性差：<25 IU/mL r=0.536；全范围 r=0.544；BioPlex 校准品在 Varelisa 回收 36–92%；vs EIA 定性 PPA 78.5%、NPA 96.7%（n=955，92 例灰区剔除） |

**ENA 单项（逐项 PPA/NPA 与临床性能）**

| 产品 / 分析物 | 目标疾病（n，标准） | 对照（n） | 灵敏度 (95% CI) | 特异度 (95% CI) | vs predicate PPA / NPA (n) |
|---|---|---|---|---|---|
| Aptiva RNP | SLE 230；MCTD 91 | 948 | SLE 37.4% (31.4–43.8)；MCTD 68.1% (58.0–76.8) | 94.8% (93.2–96.1) | vs QF RNP：97.9% / 90.6% (480) |
| Aptiva Sm | SLE 230 | 1039 | 10.4% (7.1–15.1) | 99.6% (99.0–99.9) | vs Orgentec anti-Sm ELISA：85.7% / 98.2% (418) |
| Aptiva Ro52 | SLE 230、SjS 141、SSc 217、IIM 200 | 481 | SLE 24.3%、SjS 60.3%、SSc 15.2%、IIM 19.0% | 93.6% (91.0–95.4) | vs QF Ro52：97.6% / 97.8% (1028) |
| Aptiva Ro60 | SLE 230、SjS 141 | 898 | SLE 52.2% (45.7–58.5)、SjS 66.7% (58.5–73.9) | 88.8% (86.5–90.7)；SSc 16.6%、MCTD 20.9% 阳性 | vs QF Ro60：98.5% / 91.7% (551) |
| Aptiva SS-B | SLE 230、SjS 141 | 898 | SLE 15.7%、SjS 46.8% (38.8–55.0) | 96.5% (95.1–97.6) | vs QF SS-B：96.7% / 96.3% (550) |
| Aptiva Scl-70 | SSc 217 | 1052 | 30.4% (24.7–36.8) | 94.3% (92.7–95.5) | vs QF Scl-70：95.5% / 97.6% (435) |
| Aptiva Jo-1 | IIM 200 | 1069 | 11.5% (7.8–16.7) | 99.3% (98.7–99.7) | vs QF Jo-1：96.0% / 99.7% (416) |
| Aptiva Centromere | SSc 217 | 1052 | 47.0% (40.5–53.6) | 97.0% (95.7–97.8) | vs QF Centromere：96.2% / 98.3% (449) |
| Aptiva Ribo-P | SLE 230 | 1039 | 11.7% (8.2–16.5) | 99.7% (99.2–99.9) | vs QUANTA Lite Ribosomal P：100% / 99.7% (387) |
| QF Scl-70 (K152635) | SSc 123 | 375（SLE 32、RA 31、IIM 25、MCTD 25、乳糜泻 25、AITD 25、SjS 20、Crohn 54、血管炎 15、OA 28 等） | **42.3%** (33.9–51.1) | **98.7%** (96.9–99.4) | vs QUANTA Lite Scl-70：95.5% (84.9–98.7) / 93.3% (88.1–96.3) (193，含 41 例配制样本) |
| QF Jo-1 (K151429) | IIM 206（DM 95、PM 71、JDM 7、其他 5、未分类 28） | 281（SjS 15、SLE 41、SSc 44、RA 59、MCTD 103、败血症 19） | **11.7%** (8.0–16.7) | **99.3%** (97.4–99.8) | vs FIDIS Connective 10（n=105，borderline 计阴性）：90.9% (72.2–97.5) / 86.7% (77.8–92.4) |
| QF Ro60 (K141328) | SLE 150；SS 39 | 271（sAPS 15 例排除） | SLE **24.0%** (17.4–31.6)；SS **66.7%** (49.8–80.9)；合并 32.8% | **97.4%** (94.8–99.0) | vs Hycor AUTOSTAT SS-A Ro ELISA（n=63 在测量范围内；equivocal 计阴性）：93.0% (80.9–98.5) / 80.0% (56.3–94.4) |
| QF SS-B (K141210) | SLE 290；SS 40 | 416 | SLE **13.1%** (9.4–17.5)；SS **35.0%** (20.6–51.7) | 98.0% / 97.8% | vs QUANTA Lite SS-B ELISA（n=142）：88.9% (77.4–95.8) / 92.0% (84.3–96.7) |
| EliA Ro52 (K210902) | SLE 120、SS 60、IIM 94、SSc 91（“clinically and ethnically defined”） | 390 | equivocal 计阴/阳：SLE 47.5% (38.3–56.8) / 50.8%；SS 50.0 / 55.0%；IIM 36.2 / 37.2%；SSc 20.9 / 26.4% | 96.9% (94.7–98.4) / 95.6% | vs QUANTA Flash Ro52（n=181，计阴性）：80.8% (67.5–90.4) / 98.4% (94.5–99.8) |
| EliA Ro60 (K210902) | SLE、SS（n 同上体系） | — | SLE 48.3% (39.1–57.6) / 50.8%；SS 68.3% (55.0–79.7) / 71.7% | 98.6% (97.1–99.5) / 98.4% | vs QUANTA Flash Ro60（n=104）：93.9% (85.2–98.3) / 92.1% (78.6–98.3) |
| EliA CENP / U1RNP / RNP70 (K182353) | 临床性能引用 K082759 / K083117（本次未取） | — | — | — | 仪器比对 Phadia 250 vs 2500/5000（n=100–113）：PPA 96.8–100%、NPA 80.8–100%，无阴阳互转 |
| IDS SS-A/Ro (K253817) | SLE 228、SjS 95、SSc 86 | 415 | SLE 46.9% (40.6–53.4)、SjS 80.0% (70.9–86.8)、SSc 17.4% (10.9–26.8) | 96.4% (94.1–97.8) | vs Orgentec SS-A ELISA（n=132，equivocal 计阴性）：100% / 91.1% |
| IDS SS-A/Ro 60 kDa | SLE 228、SjS 95 | 501 | 46.1% / 74.7% | 98.2% (96.6–99.1) | 89.9% (80.5–95.0) / 97.5% (148) |
| IDS SS-B/La | SLE 228、SjS 95 | 501 | 19.3% / 42.1% | 98.6% (97.1–99.3) | 100% / 98.4% (136) |
| BioPlex 单项 vs EIA (K041658，n≈1047) | — | — | — | — | PPA / NPA：dsDNA 78.5/96.7；Chromatin 58.6/91.6；Rib-P 84.2/97.8；SSA 92.3/97.8；SSB 88.9/97.9；Sm 80/97.4；Sm/RNP 83.7/98.5；RNP 85.3/96.9；Scl-70 64.7/98.2；Jo-1 87.8/100；CENP 98.5/99 |

**组合筛查**

| 产品 | 目标疾病组 | 对照 | 灵敏度 (95% CI) | 特异度 (95% CI) | vs predicate |
|---|---|---|---|---|---|
| EliA CTD 13 Screen (K253367) | 575：SLE 200（+狼疮肾炎 11）、PM 45、DM 51、MCTD 53、SS 97、SSc 57、lSSc 22、Raynaud 39 | 452：APS 47、AIH 20、细菌感染 35、乳糜泻 28、Crohn 21、Graves 27、HIV 24、桥本 10、HCV 19、白血病/淋巴瘤 40、PMR 25、PBC 20、PSC 24、RA 40、UC 22、血管炎 32、其他病毒 18 | equivocal 计阴性 **86.8%** (83.7–89.4)；计阳性 90.4%；SLE 100%、SS 100%、SSc 87.7%、MCTD 54.7%、DM 47.1% | 计阴性 **81.2%** (77.3–84.7)；计阳性 72.6%；**RA 62.5%、APS 57.4%、AIH 90% 阳性** | vs SymphonyS+dsDNA+Rib-P+RNA Pol III 组合（n=1082，计阴性）：PPA 98.3% (96.9–99.2)、NPA 88.2% (85.0–90.9) |
| EliA SymphonyS (K190710) | 194：SLE 83、狼疮肾炎 19、SS 28、SSc 32、MCTD 22、PM/DM 10 | 250：RA 40、pAPS 28、Graves 28、桥本 10、乳糜泻 28、Crohn 22、UC 22、AIH 6、PBC 5、GPA 1、细菌 30、病毒 30 | 计阴性 **52.6%** (45.3–59.8)；计阳性 55.2% | **94.8%** (91.3–97.2) / 94.4% | vs EliA Symphony（n=633：SLE 97、SS 96、SSc 87、PM/DM 78、MCTD 46 + RA 85、HBV 36、HCV 36、HIV 27、肿瘤 25、细菌 20）：计阴性 PPA 97.6% (95.0–99.0)、NPA 98.3% (96.3–99.4) |
| Euroimmun Anti-ENA Pool (K112996) | 261：MCTD 44、SLE 85、SSc 66、SjS 66 | 171：PM/DM 26、乳糜泻 21、Wegener 17、RA 39、其他自身免疫 52、感染 16 | **60.9%** (54.7–66.9)：MCTD 100%、SjS 72.7%、SLE 55.3%、SSc 45.5% | **96.5%** (92.5–98.7) | vs Aeskulisa ANA HEp-2（n=278）：97.8% (93.8–99.5) / 99.3% (96.1–100)；vs 6 种单项 ELISA 合并 PPA 100%、NPA 91.7% |
| Alegria Flash CTD Screen (K250666) | 589：SjS 135、SSc 55、CREST 55、PM 51、DM 53、MCTD 60、SLE 180 | 354：AIH 30、APS 20、肿瘤 20、乳糜泻 24、DIL 5、纤维肌痛 20、Crohn 31、UC 29、HBV 18、HCV 10、HIV 12、HSV 10、PBC 30、RA 20、血管炎 16、萎缩性胃炎 20、Graves 20、桥本 19 | **70.8%** (67.0–74.3)：SjS 83.7%、SLE 77.2%、MCTD 65.0%、DM 41.5% | **90.1%** (86.6–92.8) | vs Zeus ANA Screen ELISA (K940362)（n=943，indeterminate 计阴性）：PPA 84.2% (80.6–87.2)、NPA 87.8% (84.5–90.4) |

### 9. 厂家间差异与要点
1. **dsDNA cutoff 与单位**：虽都溯源 WHO Wo/80（IU/mL），厂家 cutoff 差异极大（BioPlex 10、EliA 15、FARRZYME 30、QUANTA Flash/Aptiva 35、Euroimmun NcX 100 IU/mL），K041658 明确指出校准品在他家平台回收仅 27–121%，故 IU/mL 不可跨平台比较。
2. **dsDNA 灰区**：QUANTA Flash/Aptiva（27–35，95th→99th percentile）、EliA（10–15）、BioPlex（5–9）设 indeterminate；FARRZYME 与 Euroimmun 无灰区。灰区计阴/阳会使灵敏度/特异度各移动 5–7 个百分点（K152013）。
3. **Farr/CLIFT 作为比较法**只见于 FARRZYME（K062183，PPA 89.8% vs Farr、67.3% vs CLIFT）和 EliA dsDNA（K072393，vs DPC RIA PPA 90.6%）；新一代 CIA/PMAT 产品均以 ELISA 或同厂 CIA 为 predicate。CLIFT 自身作为申报产品（K192916）灵敏度 32.7–48.1%、特异度 91–96%。
4. **CIA 单项 cutoff 的“二次调整”**：Inova QUANTA Flash 系列先按 CLSI C28/EP28 取 99th percentile，再用已知阳性/能力验证样本上调 RLU 阈值并统一赋值 20 CU（Scl-70 7,387→15,000 RLU；Ro60 8,194→12,000；SS-B 4,097→12,000；Jo-1 9,151→10,000），FDA 接受此做法。
5. **多重平台逐项一致性**：Aptiva vs QUANTA Flash PPA 95.5–100%、NPA 90.6–99.7%（RNP、Ro60 NPA 偏低）；EliA Ro52 vs QF Ro52 PPA 仅 80.8%；BioPlex vs 单项 EIA PPA 58.6–98.5%（chromatin、Scl-70 最低）。EliA CTD 13 相对 predicate 组合 NPA 88.2%，且在 RA（62.5%）、APS、AIH 对照中阳性率高，故整体特异度 81.2%。
6. **单项临床灵敏度**普遍偏低且依赖疾病定义：Sm 10.4%、Jo-1 11.5–11.7%、Rib-P 11.7%、SS-B（SLE）13–20%、Scl-70 30–42%、centromere 47%、Ro60（SLE）24–52%、Ro52（SjS）60%；特异度 88.8%（Ro60）–99.7%。
7. **分类标准**：绝大多数文件仅写 “clinically characterized/defined”，未指明 ACR/EULAR 标准（BioPlex 例外）。
8. **样本基质**趋势：Inova/IDS/新版 Phadia 单项限血清；Phadia 平台迁移时取消枸橼酸血浆；Euroimmun NcX 明确肝素血浆不适用。
9. **干扰研究**新近扩展到生物制剂/免疫抑制剂（belimumab、rituximab、infliximab、MTX、MMF、HCQ、AZA、CTX），见 K253367、K192916、K213403（rituximab 未评估）。

### 10. 来源
- K152013：https://www.accessdata.fda.gov/cdrh_docs/reviews/K152013.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K152013.pdf
- K072393：https://www.accessdata.fda.gov/cdrh_docs/reviews/K072393.pdf
- K083381：https://www.accessdata.fda.gov/cdrh_docs/reviews/K083381.pdf
- K062183：https://www.accessdata.fda.gov/cdrh_docs/reviews/K062183.pdf
- K192916：https://www.accessdata.fda.gov/cdrh_docs/reviews/K192916.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K192916.pdf
- K213403：https://www.accessdata.fda.gov/cdrh_docs/reviews/K213403.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K213403.pdf
- K253367：https://www.accessdata.fda.gov/cdrh_docs/reviews/K253367.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253367.pdf
- K190710：https://www.accessdata.fda.gov/cdrh_docs/reviews/K190710.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K190710.pdf
- K182353：https://www.accessdata.fda.gov/cdrh_docs/reviews/K182353.pdf
- K210902：https://www.accessdata.fda.gov/cdrh_docs/reviews/K210902.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K210902.pdf
- K152635：https://www.accessdata.fda.gov/cdrh_docs/reviews/K152635.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K152635.pdf
- K151429：https://www.accessdata.fda.gov/cdrh_docs/reviews/K151429.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K151429.pdf
- K141328：https://www.accessdata.fda.gov/cdrh_docs/reviews/K141328.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K141328.pdf
- K141210：https://www.accessdata.fda.gov/cdrh_docs/reviews/K141210.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K141210.pdf
- K112996：https://www.accessdata.fda.gov/cdrh_docs/reviews/K112996.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K112996.pdf
- K250666：https://www.accessdata.fda.gov/cdrh_docs/reviews/K250666.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K250666.pdf
- K253817：https://www.accessdata.fda.gov/cdrh_docs/reviews/K253817.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253817.pdf
- K041658：https://www.accessdata.fda.gov/cdrh_docs/reviews/K041658.pdf
