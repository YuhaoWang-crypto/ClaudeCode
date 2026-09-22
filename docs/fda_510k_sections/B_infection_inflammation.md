# 靶点组 B：感染 / 炎症（PCT、粪便钙卫蛋白、D-dimer）

> 数据来源说明：本节所有申报数据均摘自 FDA 公开的 510(k) 决策摘要（decision summary，`cdrh_docs/reviews/<K>.pdf`）或 510(k) summary（`cdrh_docs/pdfYY/<K>.pdf`），文本已抽取至 `scratchpad/txt/<K>.txt`。文件中未出现的信息一律标注「文件未载明」；来自申报文件之外的临床指南 / 生理学解释单独标注为「背景（非申报文件）」。openFDA 累计清关数来自 `python3 fda_fetch.py list <CODE>` 首行 `total=`（检索日期 2026-09-22）。

---

## 降钙素原（Procalcitonin, PCT）

### 1. 法规定位

| 项目 | 内容 |
|---|---|
| 21 CFR | **866.3215** — *Device to detect and measure non-microbial analyte(s) in human clinical specimens to aid in assessment of patients with suspected sepsis*（DEN150009 De Novo 于 2016 年建立该法规；K070310 决策摘要中写的是「21 CFR 866.3610, Endotoxin activity」+ product code NTM，属 De Novo 之前的旧归类） |
| Class | Class II (Special Controls)；DEN150009 列出 7 项特殊控制（含 IFU 声明、分析性能、前瞻性临床研究或等效样本集、无症状人群 PCT 水平评估、用户培训、结果解释与样本局限性标注） |
| Product code | **PMT**（BRAHMS PCT sensitive KRYPTOR 及 De Novo 后同类）、**PRI**（扩展用途/抗生素决策类）、**NTM**（2016 年前旧代码）、**PTF**（试剂类，如 Diazyme、DiaSys、Beckman Access PCT） |
| Panel | 83 – Microbiology |
| openFDA 累计清关数 | PMT total=5；PRI total=8；NTM total=3；PTF total=4（另按 device_name 检索「procalcitonin」共 5 条，含 PTF 代码的 DiaSys K242294、Diazyme K162297） |

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K070310 | 2008-03-31 | B·R·A·H·M·S AG / **B·R·A·H·M·S PCT sensitive KRYPTOR**（原型） | TRACE 时间分辨免疫荧光（cryptate/XL665） | 血清、血浆（肝素、EDTA；柠檬酸不推荐） | ICU 首日风险评估：进展为 severe sepsis / septic shock | 决策摘要（无 510k summary） |
| DEN150009 | 2016-02-20 | B·R·A·H·M·S GmbH (Thermo Fisher) / PCT sensitive KRYPTOR | 同上 | 血清、EDTA 或肝素血浆 | De Novo：新增 ΔPCT 28 天全因死亡风险；建立 866.3215 | De Novo 决策摘要 |
| K160911 | 2016-06-28 | bioMérieux / **VIDAS B·R·A·H·M·S PCT** | ELFA（酶联荧光，碱性磷酸酶/4-MUP） | 血清、Li-heparin 血浆（EDTA 不可用） | ICU 首日风险 + ΔPCT 28 天死亡风险（MOSES 样本再检） | 决策摘要 + 510k summary |
| K160729 | 2016-06-13 | Roche / **Elecsys BRAHMS PCT**（cobas e 411） | ECLIA 电化学发光 | 血清、K2/K3-EDTA、Li-heparin 血浆 | ΔPCT 28 天死亡风险（本次仅此项） | 决策摘要 + 510k summary |
| K162827 | 2017-02-23 | bioMérieux / VIDAS B·R·A·H·M·S PCT | ELFA | 血清、Li-heparin 血浆 | 新增：LRTI 抗生素启停决策、脓毒症抗生素停用（PMT+PRI） | 决策摘要 + 510k summary |
| K170652 | 2017-06-01 | Fisher Diagnostics (Thermo Fisher) 申报 / **ARCHITECT B·R·A·H·M·S PCT**（Abbott ARCHITECT i2000SR 平台） | CMIA 化学发光微粒子 | 血清、Li-heparin、K2-EDTA 血浆 | 四项声明（ICU 风险、28 天死亡、LRTI、脓毒症停药）；PRI/PMT/PTF | 决策摘要 + 510k summary |
| K181002 | 2018-07-16 | Siemens / **Atellica IM BRAHMS PCT** | 直接化学发光（吖啶酯，3 株鼠单抗） | 血清、EDTA、Li-/Na-heparin 血浆 | 四项声明 | 决策摘要 + 510k summary |
| K162297 | 2017-04-18 | Diazyme / **Diazyme PCT Assay**（Olympus/Beckman AU400） | 乳胶增强免疫比浊 PETIA | 血清、EDTA、Li-heparin 血浆 | 仅 ICU 首日风险评估；明确「不用于抗生素决策」；PTF | 决策摘要 + 510k summary |
| K192271 | 2019-11-26 | Beckman Coulter / **Access PCT** | 顺磁微粒化学发光 | 血清、Li-heparin、EDTA 血浆 | 仅 ICU 首日风险评估；PTF | 决策摘要 + 510k summary |
| K220262 | 2022-08-26 | Siemens / **Dimension EXL LOCI BRAHMS PCT** | LOCI 均相化学发光 | 血清（SST/RST）、Li-/Na-heparin、K2/K3-EDTA | 四项声明；PRI | 决策摘要 + 510k summary |
| K242294 | 2025-05-09 | DiaSys / **DiaSys Procalcitonin FS**（Abbott ARCHITECT c8000） | PETIA（羊多抗包被聚苯乙烯颗粒，660 nm） | 血清、Li-heparin 血浆 | 仅 ICU 首日风险评估；PTF | **仅 510k summary**（决策摘要缺失） |

### 3. 预期用途与声明类型

- **风险分层（prognostic）而非诊断**：所有产品均声明「in conjunction with other laboratory findings and clinical assessments」，警示「not indicated to be used as a stand-alone diagnostic assay」。
- 四类声明（以 K162827 / K170652 / K181002 / K220262 为完整版）：
  1. **ICU 首日风险评估**：进展为 severe sepsis / septic shock（K070310 原型声明）。
  2. **28 天全因死亡累积风险**：severe sepsis / septic shock 患者，用 ΔPCT（Day 0 或 Day 1 → Day 4）（DEN150009 首次建立）。
  3. **LRTI（CAP、急性支气管炎、AECOPD）抗生素治疗决策**：住院或急诊场景（K162827 首次建立；明确排除门诊，因门诊数据「not generalizable」）。
  4. **疑似/确诊脓毒症抗生素停用决策**（K162827）。
- 「借用型」产品（Diazyme K162297、Beckman K192271、DiaSys K242294，均 PTF）只申报第 1 项，且 Diazyme/Beckman 标签明确写入「not indicated to be used as an aid in decision making on antibiotic therapy」。
- 全部为 **prescription use only**；均为中心实验室自动化平台，**无 POC 声明**（文件未载明任何 POC/CLIA-waived 用途）。
- 警示语（K162827 起标准化）：肾功能不全可能影响 PCT；多发伤、烧伤、大手术、心源性休克可升高 PCT；非典型病原体（*Chlamydophila pneumoniae*、*Mycoplasma pneumoniae*）PCT 可不升高；<17 岁、孕妇、免疫抑制人群未在支持性 RCT 中正式分析。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**(a) ICU 首日进展风险（K070310 原型 → 各产品一致沿用）**

| 结果 | 解释 |
|---|---|
| PCT > 2.0 ng/mL (µg/L) | 进展为 severe sepsis / septic shock 高风险 |
| PCT < 0.5 ng/mL | 低风险；但「不排除感染」——局部感染或系统感染 <6 h 时可仍低 |
| 0.5–2.0 ng/mL | 灰区：结合临床背景解释，建议 6–24 h 内复测（K160911、K162297 标签） |

- K070310 决策摘要仅写「Data support the following interpretative risk assessment criteria」，**临床敏感性/特异性栏为 N/A**——即原型的 0.5/2.0 cutoff 没有在决策摘要中给出建立方法（文件未载明）。「Assay cut-off」栏另写 0.02 ng/mL（即分析灵敏度）。
- 参考人群：K070310 美国正常人群 151 例，146 例 <0.1 ng/mL（无年龄性别分布）。

**(b) LRTI 抗生素启停 & 脓毒症停药（K162827 建立，K170652/K181002/K220262 逐字照抄）**

| 场景 | PCT | 建议 |
|---|---|---|
| LRTI 启动 | < 0.10 ng/mL | 强烈不建议抗生素 |
| | 0.10–0.25 ng/mL | 不建议 |
| | 0.26–0.50 ng/mL | 建议 |
| | > 0.50 ng/mL | 强烈建议 |
| LRTI 停药 | ≤ 0.25 ng/mL **或** ΔPCT > 80%（较峰值下降） | 可停用抗生素 |
| 脓毒症停药 | ΔPCT > 80% **或** PCT ≤ 0.50 ng/mL | 可停用抗生素 |

- 决策摘要注明这些 cutoff 「unchanged from prior 510(k) submissions」并直接来自纳入 meta 分析的 RCT 所用算法（文件未逐条列出各 RCT 名称；方法比对样本来自「ProRESP trial bank」）。K162827 原始摘要曾漏掉 LRTI 停药 cutoff，2017-05-24 更正。
- 背景（非申报文件）：0.1/0.25/0.5 与 80% 下降规则即 Schuetz/Christ-Crain 系列 RCT（ProRESP 2004、ProCAP 2006、ProHOSP 2009 等）使用的算法；ProACT（2018）发表于 K162827 清关之后，未在本文件中出现。

**(c) 28 天全因死亡风险（DEN150009 → K160911 / K160729 → 后续产品）**

- ΔPCT = (PCT_Day0 − PCT_Day4) / PCT_Day0；**ΔPCT ≤ 80% = 阳性（高风险）**，**ΔPCT > 80% = 阴性（低风险）**；等价于 PCT ratio (Day4/Day0) > 0.2 vs ≤ 0.2。Day 0 缺失时可用 Day 1。
- 分层附加信息：初始 PCT ≤ 2.0 vs > 2.0 ng/mL 与 Day 4 是否仍在 ICU 联合解释（K160911 标签）。
- Cutoff 建立方式：来自前瞻性 MOSES 试验（NCT01523717）预设分析，未在决策摘要中给出 ROC 优化过程（文件未载明 80% 的推导细节，仅给出 ≤80%/＞80% 二分与 Cox 风险比）。
- 厂家提供在线「Change in Procalcitonin Calculator」（www.BRAHMS-PCT-Calculator.com），FDA 审查了其软件危害分析。

**参考人群（各产品「Expected values」）**

| 产品 | n | 定义 | 结果 |
|---|---|---|---|
| K070310 KRYPTOR | 151 | 美国正常人群 | 146/151 < 0.1 ng/mL |
| K160911/K162827 VIDAS | 200（98 男/102 女） | 表观健康；<60 岁 189，>60 岁 11 | 95th < 0.05；99th 0.09 ng/mL |
| K160729 Elecsys | 282 | 自报健康；<60 岁 276 | 95th 上限 0.08 ng/mL |
| K162297 Diazyme | 216 | 21–65 岁健康成人，CLSI C28-A3 | 参考区间 0.02–0.30 ng/mL |
| K181002 Atellica | 144 | 自报健康，EP28-A3c | 99th < 0.05 ng/mL |
| K220262 Dimension | 33（18 男/15 女） | 表观健康，6 种管型 | 全部 <0.05–0.07 ng/mL，0% 超出 KRYPTOR 参考范围 <0.1 |
| K242294 DiaSys | — | — | 文件未载明 |

**cutoff 一致性验证（「借用」逻辑）**：后续厂家不重做临床结局研究，而是（i）与 predicate（KRYPTOR 或 VIDAS）做方法比对并在 0.1/0.25/0.5/2.0 ng/mL 各决策点计算 PPA/NPA；（ii）四项声明产品用 MOSES 试验的 banked 样本再检测，计算与 KRYPTOR 的临床一致性（见第 8 节）；（iii）在 cutoff 浓度处专门验证 LoQ / 总误差与干扰（K162827 在 0.10、0.20、2.0、150 ng/mL 四个水平做干扰；Diazyme 在 0.2/0.5/2.0 报告 TE）。

### 5. 生物学 / 生理学依据

- 申报文件（K160911、K160729 标签）：「Procalcitonin (PCT) is a biomarker associated with the inflammatory response to bacterial infection」；ΔPCT 下降幅度反映感染控制与预后。K162827 引用文献：病毒/非典型病原体感染 106 例 PCT 中位数 0.09 ng/mL (IQR 0.06–0.17)；重度肾衰可能升高 PCT（「remains an area for further investigation」）。
- K070310 分析特异性：降钙素、katacalcin、α-CGRP、β-CGRP、鲑鱼/鳗鱼降钙素无交叉反应（提示抗体针对 PCT 而非其裂解产物）。
- 背景（非申报文件）：PCT 为降钙素前体（116 aa），细菌感染时由甲状腺外多种组织在 IL-6/TNF-α 及内毒素刺激下大量表达，半衰期约 24 h，感染后 6–12 h 上升；病毒感染时 IFN-γ 抑制其表达。上述解释未出现在申报文件中。

### 6. 样本类型与样本要求

| 产品 | 管型 | 基质等效性研究 | 稳定性 |
|---|---|---|---|
| K070310 KRYPTOR | 血清、肝素血浆、EDTA 血浆；**柠檬酸血浆不推荐**（PCT 下降） | 10 份加标样本、4 种管型、三重复 | 文件未载明样本稳定性 |
| K160911/K162827 VIDAS | 血清、**Li-heparin 血浆**；**EDTA 不可用**（测值下降）；同一患者需用同一管型 | 见 K071146 | 3 次冻融（K160911）；其余见 K071146 |
| K160729 Elecsys | 血清、SST、PST、Li-heparin、K2-/K3-EDTA | 53–55 对样本各管型 vs 血清，Passing-Bablok 无显著差异 | 2–8°C 48 h；15–25°C 24 h；−20°C 12 个月；1 次冻融 |
| K170652 ARCHITECT | 血清、SST、Li-heparin、K2-EDTA | 4 管型配对（瑞士两家医院 ED/ICU/住院），以 K2-EDTA 为基准，偏差 ≤10% | 室温 24 h；机上 3 h（53 例 0.02–62.28 ng/mL） |
| K181002 Atellica | 血清、EDTA、Li-/Na-heparin | 51 组匹配（商业来源，加标），斜率 vs 血清 1.03–1.05，r 1.00 | 文件未载明（决策摘要） |
| K162297 Diazyme | 血清、K3-EDTA、Li-heparin | 40 组配对加标，斜率 1.0±0.1，R² ≥0.95 | 决策摘要仅载校准品/质控稳定性；样本稳定性文件未载明 |
| K220262 Dimension | SST、RST、Li-/Na-heparin、K2/K3-EDTA | 76 组匹配加标，Passing-Bablok 斜率 0.98–1.00 | 文件未载明（决策摘要） |
| K242294 DiaSys | 血清、Li-heparin | ≥20 例配对 Passing-Bablok，斜率 0.85–1.15 | 20–25°C 24 h；2–8°C 5 d；−20°C 14 d；2 次冻融 |

### 7. 分析性能验证所依据的标准

| K 号 | 文件 "Standards/Guidance Referenced" | 对应项目 |
|---|---|---|
| K070310 | CLSI EP5-A2（精密度）；EP6-A（线性）；EP9-A（方法比对）；EP17-A（LoB/LoD/LoQ） | 内/外部精密度 3 站点、线性 0.02–5000、LoD 0.0227、LoQ 0.075 ng/mL（TE ≤30%） |
| K160911 | EP17-A2；EP5-A3；EP6-A；EP9-A3 | LoD 0.03、LoQ 0.05；VIDAS vs VIDAS 3 Deming 比对 |
| K160729 | EP05-A3；EP06-A；EP17-A2 | LoD 0.0181（声明 0.02）；LoQ 0.045（声明 0.06） |
| K162827 | EP05-A3；EP07-A2；EP17-A2 | 低端精密度 0.05–0.25 ng/mL 补充研究；干扰（4 种内源 + 36 种药物）；LoQ 在 0.10/0.25 处 TE |
| K170652 | EP05-A3；EP06-A；EP07-A2；EP09-A3；EP15-A3；EP17-A2；EP25-A；EP28-A3c | 全套；LoB 0.0004 / LoD 0.0018 / LoQ 0.0077 ng/mL |
| K181002 | EP05-A3；EP06-A；EP07-A2；EP09-A3；EP17-A2；EP28-A3c；ANSI/AAMI/ISO 14971 | LoD 0.03、LoQ 0.04 ng/mL；参考区间 144 例 |
| K162297 | EP5-A2；EP6-A；EP7-A2；EP9-A2；EP17-A2；C28-A3（另在携带研究引用 EP10-A2） | LoB 0.060 / LoD 0.16 / LoQ 0.20；AMR 0.20–52 ng/mL |
| K220262 | EP05-A3；EP07 3rd；EP09c 3rd；EP17-A2；EP25-A；EP28-A3c；EP34（扩展测量区间）；EP06 2nd | LoD 0.04、LoQ 0.05；AMI 0.05–50，EMI 至 1000 |
| K242294 | EP05-A3、EP15-A3；EP06-A；EP17-A2；EP07-A3；EP25-A；EP09-A3；ISO 14971:2019 | LoB 0.081；LoD = LoQ 0.23 ng/mL；稳定性 24 个月 |
| DEN150009 | EP05-A2；EP06-A；EP07-A2；EP09-A2；EP25-A；EP17-A | De Novo 分析包 |

（未见任何 PCT 文件引用 EP12 或 EP24。）

### 8. 临床验证设计与结果

**K070310（原型）**：无临床敏感性/特异性研究（N/A）；仅方法比对 vs BRAHMS PCT LIA（K040887），3 站点 184 份样本，Passing-Bablok y = 0.95x + 0.03，r² = 0.98。

**K160911 / K160729 / DEN150009 — MOSES 研究（28 天死亡）**：前瞻性、美国 13 中心，858 例成人 severe sepsis / septic shock 患者，Day 0、1、4 测 PCT；per-protocol 598 例（44% 女，平均 64 岁；severe sepsis 51% / shock 49%；社区获得性 91%；Day 4 仍在 ICU 44%）；PP 死亡率 16.8%（全人群 22%）。主要分析：Fisher 精确检验 ΔPCT（≤80% vs >80%）vs Day 28 存活；Cox 回归。
- VIDAS 3（K160911）：ΔPCT Day0→4 ≤80% vs >80% 的 HR 2.27 (1.41–3.63, p=0.0007)；VIDAS：HR 2.05 (1.30–3.23)。Day 4 仍在 ICU 者：ΔPCT>80% 死亡率 18.4% vs ≤80% 31.3%；预后敏感性 78.4%、特异性 35.7%。非 ICU：5.4% vs 11.4%。校正 APACHE/SOFA/年龄/Day 4 位置后 ΔPCT 仍独立（HR 1.6–2.1）。
- Elecsys（K160729）：ICU 组 Day0→4：ΔPCT>80% 死亡率 22.1% (13.3–31.0) vs ≤80% 29.6% (22.9–36.4)；敏感性 73.4%、特异性 35.0%；非 ICU 组 5.6% vs 11.0%，敏感性 72.3%、特异性 44.4%。基线 PCT 分层死亡率：<0.5 15.2%、0.5–2.0 12.5%、>2.0 19.5%。
- 对比：单独「PCT Day 0 >2 vs ≤2」HR 仅 1.38 (0.89–2.14, p=0.149)——即绝对值分层在此人群对死亡的预测弱于 ΔPCT。

**K162827 — 抗生素决策（meta 分析替代诊断准确性）**：系统文献回顾 + 研究级与患者级 meta 分析（PubMed/Cochrane，Cochrane 偏倚工具，漏斗图，固定/随机效应）。2016-11-10 Microbiology Devices Panel 会议审议。
- LRTI：研究级 11 RCT / 4090 例；患者级 13 RCT / 3142 例。抗生素启动 OR 0.26 (0.13–0.52)（研究级）；患者级 PCT 组 71.4% vs 标准组 88.4%（adj OR 0.27）；疗程中位 7 vs 10 d（−2.87 d）；总暴露 5 vs 9 d（−3.60 d）；30 天死亡 6.7% vs 7.4%（OR 0.95, 0.77–1.16）；并发症 18.0% vs 21.1%（OR 0.82, 0.68–0.99）；住院日 −0.18 d。按初始 PCT 分层：<0.10 组启动率 34.5% vs 71.4%；0.10–0.25 组 57.2% vs 87.4%；>0.50 组两组均 ≈99%。算法依从率 59–91%。DOOR/RADAR 分析支持 PCT 指导组更优。
- 脓毒症停药：研究级 10 RCT / 3489 例；患者级 5 RCT / 598 例；疗程 −1.49 d (−2.27, −0.71)；死亡 RR 0.90 (0.79–1.03)；患者级总暴露 8 vs 12 d（−3.20 d）；30 天死亡 19.9% vs 23.8%（OR 0.87）。依从 47–91%。
- 方法比对 vs KRYPTOR（ProRESP bank 203 份）：0.10 ng/mL PPA 83.7% / NPA 89.2% / OPA 86.7% / κ 0.731（未达预设 κ 标准但被接受，因不一致不跨越两级）；0.25：94.6/98.6/97.5/κ 0.938；0.50：100/98.8/99.0/κ 0.971；2.00：100/97.3/97.5/κ 0.870。

**K170652 ARCHITECT（一致性验证型）**：
- 方法比对 EP09-A3：ICU 血清 142 份（130 天然 + 12 人工）vs KRYPTOR：Weighted Deming ARCHITECT = −0.02 + 1.00×KRYPTOR，r 0.99；0.5 ng/mL：NPA 97.96% (48/49)、PPA 95.70% (89/93)；2.0：NPA 100% (77/77)、PPA 100% (65/65)。
- 临床一致性：MOSES banked 样本 n=2331（<0.1 µg/L 84 例；<0.25 351；<0.5 594；<2.0 1091）：PPA/NPA/总一致/κ — 0.10：96.1%/95.2%/96.1%/0.619；0.25：96.8/95.2/96.5/0.871；0.50：96.9/96.8/96.9/0.920；2.00：96.9/98.4/97.6/0.951；Passing-Bablok 斜率 0.95 (0.94–0.96)，截距 −0.04；范围 0.02–862.43 µg/L。

**K181002 Atellica**：方法比对 623 份血清（KRYPTOR 赋值 0–1000 ng/mL，4 台仪器 4 批试剂），n=522（0.06–49.20）回归：Weighted Deming 斜率 1.02 (0.99–1.05)、截距 −0.02，r 0.98；Passing-Bablok 斜率 1.06、截距 −0.04。决策点一致性：0.1：PPA 99.3% (539/543)、NPA 95.0% (76/80)；0.25：99.0/94.6；0.5：96.7/97.4；2.0：97.2/97.6。临床一致性 MOSES banked n=2285：0.10：PPA 98.6%、NPA 77.6%、总 97.8%、κ 0.710；0.25：98.7/83.9/96.4/0.854；0.50：97.5/92.8/96.3/0.903；2.00：95.3/97.6/96.4/0.927；Passing-Bablok 斜率 0.941 (0.927–0.954)，截距 0.017。

**K162297 Diazyme**：方法比对 vs VIDAS（K071146）3 站点 219 份 ICU 血清：0.5：NPA 85.3% (29/34)、PPA 96.2% (178/185)；2.0：NPA 97.1% (136/140)、PPA 91.1% (72/79)；Passing-Bablok 总斜率 0.944 (0.882–0.984)，截距 0.001；Weighted Deming 斜率 0.866。**独立临床研究**：116 例连续入 medical ICU 患者（21–93 岁）首日血清，ACCP/SCCM 分类（无感染 4、SIRS 26、sepsis 18、severe sepsis 36、septic shock 32）：以 severe sepsis/shock 为阳性，0.5 ng/mL：敏感性 100% (68/68)、特异性 43.8% (21/48)；2.0 ng/mL：敏感性 97.1% (66/68)、特异性 91.7% (44/48)。

**K192271 Beckman Access PCT**：vs VIDAS 207 份：0.5：NPA 94.0% (78/83)、PPA 100% (124/124)；2.0：NPA 98.4%、PPA 97.6%；Weighted Deming 斜率 0.96、截距 0.02、r 0.99。

**K220262 Dimension EXL LOCI**：595 份天然血清（0.05–1000 ng/mL，含 ≥60 岁 402 例）vs KRYPTOR；预设标准：Passing-Bablok 斜率 1.00±0.10、截距 ≤LoQ (0.05)、r ≥0.950；PPA/NPA ≥75%（0.10、0.25）、≥85%（0.50、2.00）。结果（lot FB1218）：斜率 1.07 (1.05–1.09)、截距 −0.01、r 0.96；0.10：PPA 97.96%、NPA 89.09%；0.25：97.78/92.36；0.50：96.53/95.58；2.00：97.85/97.51。**临床研究栏「Not applicable」**——即 2022 年后四项声明可仅凭方法比对借用，无需 MOSES 样本再检。

**K242294 DiaSys**（仅 510k summary）：两次方法比对 vs VIDAS：2021 年 n=120，Passing-Bablok 斜率 1.08、截距 0.105、r 0.991；0.5：NPA 93.9% (79.8–99.3)、PPA 100%；2.0：NPA 93.3%、PPA 97.8%。2025 年 n=210，斜率 0.940、截距 0.017、r 0.965；0.5：NPA 98%、PPA 100%；2.0：NPA 100%、PPA 94%。无独立临床结局研究。

### 9. 厂家间差异与要点

1. **声明范围两极分化**：BRAHMS 授权体系（KRYPTOR、VIDAS、Elecsys、ARCHITECT、Atellica、Dimension、LIAISON、Lumipulse 等以 "BRAHMS" 冠名）可拿到四项声明；非授权 PETIA 产品（Diazyme、DiaSys）及 Beckman Access 仅 ICU 风险评估，并被要求写明「不用于抗生素决策」。
2. **借用 cutoff 的证据链随时间放宽**：2017–2018（ARCHITECT、Atellica）需用 MOSES banked 样本 n≈2300 做四决策点一致性；2022（Dimension）仅需 595 份方法比对 + 预设 PPA/NPA 阈值，临床研究 N/A。
3. **低端一致性最弱**：0.10 ng/mL 处 NPA 最低（VIDAS vs KRYPTOR NPA 89.2%、κ 0.73；Atellica MOSES 样本 NPA 77.6%），因此 FDA 要求在 0.10/0.25 处单独证明 LoQ/TE（VIDAS：0.10 处 TE ≤30%、0.25 处 ≤20%）。PETIA 平台 LoQ 0.20–0.23 ng/mL，**无法覆盖 0.10/0.25 决策点**，是其只能申报 ICU 声明的分析学原因之一（文件未明言，属推断）。
4. **抗凝剂差异**：VIDAS 禁用 EDTA；KRYPTOR 禁用柠檬酸；Elecsys/ARCHITECT/Atellica/Dimension 支持 EDTA。
5. **参考区间**：各家健康人群 95th/99th 百分位均 <0.1 ng/mL；Diazyme 报告 0.02–0.30 ng/mL（较宽，与其 LoQ 0.20 相关）。
6. **软件/报告要求**：四项声明产品的标签须建议实验室报告绝对值并附 ΔPCT (≤80%/>80%) 及计算器链接（K160729、K170652）。

### 10. 来源

- K070310：https://www.accessdata.fda.gov/cdrh_docs/reviews/K070310.pdf
- DEN150009：https://www.accessdata.fda.gov/cdrh_docs/reviews/DEN150009.pdf
- K160911：https://www.accessdata.fda.gov/cdrh_docs/reviews/K160911.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K160911.pdf
- K160729：https://www.accessdata.fda.gov/cdrh_docs/reviews/K160729.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K160729.pdf
- K162827：https://www.accessdata.fda.gov/cdrh_docs/reviews/K162827.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162827.pdf
- K170652：https://www.accessdata.fda.gov/cdrh_docs/reviews/K170652.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K170652.pdf
- K181002：https://www.accessdata.fda.gov/cdrh_docs/reviews/K181002.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181002.pdf
- K162297：https://www.accessdata.fda.gov/cdrh_docs/reviews/K162297.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162297.pdf
- K192271：https://www.accessdata.fda.gov/cdrh_docs/reviews/K192271.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K192271.pdf
- K220262：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220262.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220262.pdf
- K242294：https://www.accessdata.fda.gov/cdrh_docs/pdf24/K242294.pdf（决策摘要 reviews/K242294.pdf 未获取到）

---

## 粪便钙卫蛋白（Fecal Calprotectin, fCAL）

### 1. 法规定位

| 项目 | 内容 |
|---|---|
| 21 CFR | **866.5180** — *Fecal calprotectin immunological test system*（K050007 De Novo 建立；特殊控制指南 *Class II Special Controls Guidance Document: Fecal Calprotectin Immunological Test Systems*, July 2006） |
| Class | Class II (Special Controls) |
| Product code | **NXO** – Calprotectin, Fecal（提取装置也归 NXO；校准品 JIT） |
| Panel | 82 – Immunology |
| openFDA 累计清关数 | NXO total=15（含 DEN060001 PhiCal；device_name 含「calprotectin」检索到 6 条） |

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K050007 (DEN060001) | 2006-04-26 | Genova Diagnostics / **PhiCal Test**（原型，De Novo） | 手工 ELISA（兔多抗 / 碱性磷酸酶，405 nm） | 粪便，0.1 g + 5 mL 提取液，终稀释 1:2500 | 辅助诊断 IBD（CD、UC）并与 IBS 鉴别 | 决策摘要 |
| K130945 | 2014-01-16 | Eurospital / **Calprest**（PhiCal 同一生产商） | ELISA，6 点校准 | 粪便，1:50 提取 → 1:2500 | 同上 | 决策摘要 + 510k summary |
| K170993 | 2017-12-22 | Inova / **QUANTA Flash Calprotectin**（BIO-FLASH） | 顺磁微粒化学发光（兔多抗捕获/单抗示踪） | 粪便提取液（手工称重） | 同上 | 决策摘要 + 510k summary |
| K180971 | 2018-10-16 | Inova / **Fecal Extraction Device (FED)** | 附件（体积法取样 ~56 mg） | 粪便 | 同上（附件） | 决策摘要 + 510k summary |
| K181012 | 2018-06-04 | BÜHLMANN / **fCAL ELISA** | ELISA（单抗/HRP/TMB，450 nm） | 粪便 <1 g，1:7500 | 同上 | 决策摘要 + 510k summary |
| K182698 | 2018-12-26 | DiaSorin / **LIAISON Calprotectin + Q.S.E.T. Device** | 顺磁微粒 CLIA（异鲁米诺） | 粪便；手工称重或 Q.S.E.T. 装置（~10.5 mg） | 同上 | 决策摘要 + 510k summary |
| K190784 / K191718 / K232057 | 2019-06-25 / 2019-09-24 / 2024-02-06 | BÜHLMANN / **fCAL turbo + CALEX Cap**（Roche cobas c501/c502） | PETIA（禽多抗聚苯乙烯纳米颗粒，546 nm） | 粪便；手工 1:50 或 CALEX Cap 1:500 | 同上；K232057 为 Special 510(k) 延长提取液稳定性 | 决策摘要 + 510k summary |
| K213858 | 2022-07-26 | DiaSorin / **LIAISON Q.S.E.T. Device Plus** | 附件（预装 6 mL 缓冲液；LIAISON XL/XS） | 粪便 BSFS 2–7 | 同上（附件） | 决策摘要 + 510k summary |
| K220763 | 2023-04-13 | ALPCO / **Calprotectin Immunoturbidimetric Assay**（Beckman AU680） | PETIA（鼠单抗乳胶） | 粪便；Easy Stool Extraction Device 或手工 50–100 mg + 99× 缓冲液 | 同上 | 决策摘要 + 510k summary（openFDA 标注 Statement，但决策摘要可获取） |

### 3. 预期用途与声明类型

- 统一措辞：「aid in the diagnosis of inflammatory bowel disease (IBD), specifically Crohn's disease (CD) and ulcerative colitis (UC), and aid in the differentiation of IBD from irritable bowel syndrome (IBS), in conjunction with other laboratory and clinical findings」——属**辅助诊断 / 鉴别诊断**，非筛查、非监测（文件中无疾病活动度监测或复发预测声明）。
- 定量报告（µg/g 或 mg/kg），但临床解释为三分类（正常 / 灰区 / 升高）。
- 全部 prescription use；均为实验室检测（手工 ELISA、免疫分析仪或生化分析仪），**无 POC 或家用声明**。
- 年龄：K182698 纳入 ≥4 岁；K181012/K190784 纳入 4–21 岁儿科亚组；K220763 仅 ≥22 岁；K050007 文件未载明年龄范围。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 正常 | 灰区 / 临界 | 升高 | 灰区处理 |
|---|---|---|---|---|
| PhiCal K050007 | <15.6–50 µg/g | 50–120 µg/g | >120 µg/g | 4–6 周复测 |
| Calprest K130945 | <50 | 50–120 | >120 | 4–6 周复测（决策摘要「Assay cut-off」表；同文件差异表中另写 50–100/>100，前后不一致） |
| QUANTA Flash K170993 | <50 mg/kg | ≥50–<120 | ≥120 | 文件未载明复测建议 |
| LIAISON K182698 | <50 | 50–120 | >120 | — |
| fCAL ELISA K181012 / turbo K190784 | **<80** | **80–160** | **>160** | 4–6 周内复测 |
| ALPCO K220763 | <50 | **50–100** | **>100** | 4–6 周复测 |

**cutoff 建立方式**

- **K050007（原型）**：以 124 例正常样本（文件写「normal sera」）均值 40 µg/g、SD 40 → 上限 120 µg/g；另在 908 例总人群做 ROC，120 µg/g 对应特异性 95%、敏感性 70%；50 µg/g 为健康人群「expected value <50」。即**「均值+2SD」+ ROC 双重建立，50/120 双 cutoff 形成灰区**。
- **K170993 Inova**：61 例参考人群（53 健康 + 8 例息肉/鳞癌等）非参数 95th 百分位 137.8 mg/kg (90% CI 125.0–148.1)，再用 31 例已知 IBD（无一 <120）把 cutoff 定在 120 mg/kg，「以尽量减少假阴性」；依据 CLSI EP28-A3c。
- **BÜHLMANN 80/160**：K181012/K190784 决策摘要未给出 80/160 的推导过程，「Clinical cut-off: Not applicable」（文件未载明），仅在临床研究中同时报告 ≥80 与 ≥160 两个 cutoff 的敏感性/特异性。
- **ALPCO 50/100**：沿用其 predicate K191807（文件未载明推导）。
- **Predicate 对齐**：Calprest（PhiCal 同厂）、LIAISON（predicate PhiCal）直接沿用 50/120。
- 灰区一律双算：「borderline 视为阳性」与「borderline 视为阴性」分别报告敏感性/特异性。

**参考人群**

| 产品 | n / 定义 | 结果 |
|---|---|---|
| K050007 | 161 例正常（在 908 例队列内） | 预期值 <50 µg/g |
| K130945 | 90 例健康 | 85/90 (92.2%) 阴性；5 例 51.8–60.4 mg/kg |
| K170993 | 164 例表观健康（95 女/70 男，17–89 岁，均值 44.9 岁） | 全部 <120；均值 19.0，范围 <16.1–49.2 mg/kg |
| K182698 | 112 成人 + 15 儿童（3–21 岁），排除慢性腹部主诉、IBD/IBS 史、NSAIDs | 均值 23.5、中位 15.0、范围 5.0–103 µg/g；90% 中心区间 5.0–79.8；0% 升高、11.8% 临界 |
| K190784 | 141 例健康成人 ≥21 岁 | <80：75.2%；80–160：12.8%；>160：12.0% |
| K220763 | 120 例无症状（65 女/55 男，22–82 岁，均值 43.8），EP28-A3c | 104/120 <11 µg/g；中位 3.9；上限 67.5 (90% CI 20.0–130.3)；2 例升高（111.2、130.3）、1 例灰区 |
| K181012 | — | 「Not applicable」（文件未载明） |

### 5. 生物学 / 生理学依据

- 申报文件：「fecal calprotectin, a neutrophilic protein that is a marker of (intestinal) mucosal inflammation」（K050007、K182698、K190784 IFU）；K190784 注明分析物为 MRP8/14（S100A8/A9）二聚体；K050007 注明消化道出血 100 mL/天 仅使 fCAL 升高 ≤15 µg/g；K050007 称 IBD 患者水平可达健康上限的「五倍至数千倍」，IBS 一般不升高但可重叠。
- 干扰研究显示常见肠道菌（*E. coli*、*Salmonella*、*Shigella*、*Yersinia*、*Klebsiella*、*Citrobacter*）、5-ASA、硫唑嘌呤、泼尼松、PPI、抗生素、铁剂、血红蛋白等不干扰测定（K050007、K130945、K181012、K182698、K220763）。
- 背景（非申报文件）：钙卫蛋白占中性粒细胞胞浆蛋白约 60%，肠黏膜炎症时中性粒细胞迁移入肠腔而释放；在粪便中室温稳定数日；NSAIDs、感染性肠炎、结直肠肿瘤亦可升高（申报文件仅以排除 NSAIDs 使用者的入组标准间接体现）。

### 6. 样本类型与样本要求

| 产品 | 提取方式 | 样本量 / 稀释 | 稳定性 |
|---|---|---|---|
| K050007 PhiCal | 手工称重 ~0.1 g + 5 mL 提取液，振荡 30 min，离心，上清 1:50 | 终稀释 1:2500 | 原便 2–8°C 或变温条件下 11 天内稳定；提取液 −20°C 11 天；7 天内检测或 −20°C 冻存 |
| K130945 Calprest | 1:50 提取 → 1:2500 | 1–5 g 送检 | 原便 2–8°C 4 天（±5% 标准）；提取液 −20°C 3 个月；提取重复性 CV 7.0–13.6% |
| K181012 fCAL ELISA | 手工提取 | <1 g；1:7500 | 原便 2–8°C 3 天；提取液 2–8°C 7 天、−20°C 3 年、3 次冻融；提取重复性总 CV 9.5–20.5% |
| K190784/K191718 fCAL turbo + CALEX Cap | 手工 1:50 或 CALEX Cap 1:500 | — | 手工提取液 2–8°C 11 天、−20°C 2 个月；CALEX 提取液 2–8°C 3.5 天（K191718）→ 2–8°C 15 天、18–28°C 48 h（K232057）；CALEX Cap 装置 24 个月，运输 37°C 4 周 |
| K182698 LIAISON + Q.S.E.T. | 手工称重或 Q.S.E.T.（平均 10.5 mg，CV 4.9%，BSFS 2–6） | AMR 5–800 µg/g，自动 1:10 至 8000 | 原便 2–8°C 72 h、−20°C 16 周、3 次冻融；手工提取液 2–8°C 7 天；Q.S.E.T. 未离心 2–8°C 6 h、离心后 8 天 |
| K213858 Q.S.E.T. Plus | 预装 6 mL；BSFS 2–7（液体便可直接移液 12 µL） | 平均 11.3 mg，CV 9.4% | 未离心 −20°C 16 周、2–8°C 72 h、20–24°C 6 h；提取后 −20°C 7 天、1 次冻融 |
| K170993/K180971 QUANTA Flash + FED | 手工称重或 FED（平均 56 mg，CV 6.9%，BSFS 2–5） | AMR 16.1–3500 mg/kg | FED 提取液室温 72 h、2–8°C 14 天、4 次冻融；标签警示 FED 精密度较手工差 |
| K220763 ALPCO | Easy Stool Extraction Device（不用于水样便）或手工 50–100 mg + 99× 缓冲液，涡旋 30 min，3000 g 5 min | 1:100 | 提取液 2–8°C 3 天、−80°C 14 天；原便稳定性见 K191807 |

基质比对：各文件均注明「Stool is the only matrix」（不适用）；粪便稠度（Bristol Stool Form Scale）被作为提取装置精密度/方法比对的自变量（K180971、K182698、K191718、K213858）。

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 对应项目 |
|---|---|---|
| K050007 | **None referenced** | 精密度（20 重复；6 池 20 次提取）、线性（水相/基质）、加标回收 99–119%、功能灵敏度（20% CV → 6.25 ng/mL = 15.6 mg/kg） |
| K130945 | EP05-A2；EP06-A；EP7-A2；EP17-A；Special Controls Guidance (2006) | LoB 3.04、LoD 3.98 mg/kg；线性 15.6–500；干扰 |
| K170993 | Special Controls Guidance；EP05-A3；EP17-A2；EP06-A；EP28-A3c | 精密度、LoB/LoD、线性、参考区间（非参数 95th） |
| K180971 | EP05-A3 | FED 取样重量精密度、提取重复性 |
| K181012 | EP5-A3；EP6-A；EP7-A2；EP17-A2；EP25-A；EP28-A3c；X5-R（计量溯源） | LoB 8.3 / LoD 12.6 / LoQ 30 µg/g；线性 30–1800；稳定性 24 个月 |
| K182698 | EP05-A3；EP06-A；EP07-A2；**EP12-A2**（定性一致性）；EP15-A3；EP17-A2；EP28-A3c | LoB 0.107 / LoD 0.395 / LoQ 0.400 µg/g；hook 至 100,000 µg/g；PPA/NPA |
| K190784 / K191718 | Special Controls Guidance；EP05-A3；EP6-A；EP07-A2；EP10-A2；EP17-A2；EP25-A | LoB 16.7 / LoD 23.7 / LoQ 30 µg/g；携带污染；CALEX 稳定性 |
| K213858 | EP05-A3（另引用 EP25-A 做装置稳定性） | Q.S.E.T. Plus 取样精密度、稳定性 |
| K220763 | EP05-A3；EP06-Ed2；EP07 3rd；EP17-A2；EP28-A3c；Special Controls Guidance | LoB 2.2 / LoD 3.9 / LoQ 11.0 µg/g；线性 11–1000；hook 至 18,000 |
| K232057 | EP25（Special 510(k)） | CALEX 提取液稳定性延长 |

（未见引用 EP09（方法比对均用 Passing-Bablok/Deming 但未列标准）、EP24；溯源：各家均声明「无国际参考物质」，溯源至内部标准——PhiCal 溯源至 Fagerhol 实验室「Gold Standard」制剂，Calprest 溯源至 K050007 内部参考，fCAL 溯源至重组或人血清钙卫蛋白，ALPCO 溯源至其 CL-ELISA predicate。）

### 8. 临床验证设计与结果

**K050007 PhiCal（原型）**：908 例（IBD 255、IBS 410、其他肠病 82、正常 161）；金标准文件未载明（仅按诊断分组）；三分类：>120 µg/g 254 IBD / 29 非 IBD；50–120 43/89；<50 42/451。排除 132 例 borderline 后：敏感性 86% (254/296)、特异性 94% (451/480)。40 份样本与挪威 Fagerhol 实验室分样比对 y = 0.9603x + 9.7691，r² 0.9618。

**K130945 Calprest**：138 例（IBD 98：CD、UC、IC；非 IBD 40：IBS、慢性腹泻、RAP、乳糜泻），IBD「由临床发现和/或结肠镜确认」。borderline 视为异常：敏感性 96.9% (91.3–99.4)、特异性 85.0% (70.2–94.3)、PPV 94.1%、NPV 91.9%；borderline 视为正常：敏感性 79.6%、特异性 92.5%、NPV 64.9%。方法比对 vs PhiCal n=131，Deming y = 0.98x − 1.85（斜率 0.96–1.01）；定性一致 PPA 94.9%、NPA 98.1%。

**K170993 QUANTA Flash**：165 例经**回肠结肠镜**（所有患者）+ 组织学 + Rome III 诊断：IBD 58（CD 31、UC 26、淋巴细胞性结肠炎 1），非 IBD 107（IBS 75、乳糜泻 6、慢性腹泻 10、胃炎 5、RAP 10、小肠梗阻 1）；排除 <14 岁、既往 IBD 诊断、过多黏液。intermediate 视为阴性：敏感性 89.5% (78.9–95.1)、特异性 86.1% (78.3–91.4)、PPV 77.3%、NPV 93.9%；视为阳性：96.5% / 74.1% / 66.3% / 97.6%。vs QUANTA Lite ER：n=77（AMR 内）斜率 1.10、r 0.956；定性 PPA 98.5%、NPA 94.4%。

**K181012 fCAL ELISA**：美国 13 中心前瞻性，478 例（415 成人 + 63 例 4–21 岁），有 IBD/IBS 症状并转诊内镜者（含已知 IBD 疑似复发），粪便在内镜前 ≥1 天或内镜后 ≤3 天采集；最终诊断由研究者依据**内镜 + 组织学 + Rome III**；可评估 337 例：IBD 135（成人 102/儿科 33）、IBS 130、其他 GI 72。分布：IBD <80 6.7%、80–160 8.9%、>160 84.4%；IBS 72.3%/13.1%/14.6%。IBD vs 非 IBD：≥80 敏感性 93.3% (87.7–96.9)、特异性 70.3% (63.5–76.5)；≥160：84.4% / 83.7%。IBD vs IBS：≥80：93.3% / 72.3%；≥160：84.4% / 85.4%。vs Calprest NG n=371：80 处 PPA 96.2%、NPA 82.5%；160 处 PPA 98.8%、NPA 85.1%。

**K190784 fCAL turbo**：同一 337 例队列：IBD vs 非 IBD ≥80：敏感性 91.1% (85.0–95.3)、特异性 74.3% (67.7–80.1)；≥160：80.0% / 85.1%；IBD vs IBS ≥80：91.1% / 76.2%；≥160：80.0% / 87.7%。vs fCAL ELISA n=220：Passing-Bablok 斜率 1.025 (0.990–1.058)、截距 −4.5、r 0.972；80 处偏差 −3.1%、160 处 −0.3%；定性 n=248：PPA80 93.6%、NPA80 91.3%、PPA160 93.9%、NPA160 95.3%。**K191718 CALEX vs 手工**（n=202）：斜率 1.149 (1.100–1.201)、截距 −8.3、r 0.921；80 处偏差 +4.6%、160 处 +9.7%；定性 n=241：PPA80 98.1%、NPA80 89.9%、PPA160 97.6%、NPA160 96.6%。

**K182698 LIAISON**：美国 14 中心前瞻性，411 例入组 → 240 例可评估（160 女/80 男；19 例 <22 岁）；纳入 ≥4 岁、有 IBD/IBS 症状、经结肠镜确诊；排除肠切除/转流、结肠镜前 7 天内 NSAIDs、6 个月内免疫调节剂/生物制剂、孕哺；IBD 102（CD 用 SES-CD、UC/IC 用 Mayo 内镜评分）、IBS 67（1 年内阴性结肠镜 + Rome III）、其他 GI 71（憩室病 35、慢性腹泻 16、RAP 11、艰难梭菌 3 等）。IBD vs 非 IBD：borderline 视为升高——敏感性 98.0% (100/102; 93.1–99.8)、特异性 66.8% (95/138; 60.4–76.7，文件原值)；视为正常——88.2% (90/102)、90.6% (125/138)。IBD vs IBS：98.0% / 65.7%；88.2% / 88.1%。vs PhiCal n=164：Passing-Bablok 斜率 0.97 (0.89–1.00)、截距 1.50、R² 0.96；PPA 96.9%/NPA 88.9%（borderline=阳性）。Q.S.E.T. vs 手工 n=128：斜率 0.96、50 处偏差 −2.91、120 处 −5.43 µg/g。

**K213858 Q.S.E.T. Plus vs Q.S.E.T.**：n=159，斜率 0.88 (0.84–0.92)、120 处偏差 −14.1%、R² 0.967；PPA 97.8%/NPA 100%（equivocal=阴性）。

**K220763 ALPCO PETIA**：15 中心前瞻性 349 例 ≥22 岁拟行结肠镜者；IBD 依内镜和/或活检组织学；IBS 依 Rome IV + 结肠及末端回肠内镜阴性；排除已知肠癌、化疗/免疫抑制、已诊断且用药的 IBD、肠道感染、上消化道病、2 周内 PPI/H2RA、2 周内 ≥7 天 NSAIDs。IBD 63（UC 29、CD 24、未定 10）、IBS 105、其他 181。IBD vs 非 IBD：equivocal 视为升高（≥50）——敏感性 90.5% (80.7–95.6)、特异性 93.4% (89.9–95.7)、PPV 75.0%、NPV 97.8%；视为正常（>100）——76.2% / 97.6% / 87.3% / 94.9%。IBD vs IBS（n=168）：≥50：90.5% / 90.5%；>100：76.2% / 96.2%。vs ALPCO CL-ELISA n=168：斜率 0.98 (0.95–1.02)、R² 0.99；equivocal 视为升高 NPA 80.0%、PPA 93.2%。

（本组文件均未报告 AUC；评价指标为固定 cutoff 下的敏感性/特异性/PPV/NPV 与 predicate 的 PPA/NPA/Passing-Bablok。）

### 9. 厂家间差异与要点

1. **cutoff 体系分裂**：PhiCal 系（PhiCal、Calprest、LIAISON、Inova）用 50/120；BÜHLMANN 用 80/160；ALPCO 用 50/100。因无国际参考物质，各家溯源至自家内部标准，跨平台数值不可直接互换（fCAL ELISA vs Calprest NG 在 80 处 NPA 仅 82.5%）。
2. **敏感性-特异性权衡**：灰区归阳性时敏感性 90–98%、特异性 66–93%；归阴性时敏感性 76–88%、特异性 88–98%。特异性最高的是 ALPCO（其他 GI 疾病占比大 181/349，且排除标准最严）。
3. **金标准演进**：原型 K050007 未描述金标准；2017 年后（Inova、BÜHLMANN、DiaSorin、ALPCO）均要求全员内镜 ± 组织学，IBS 需内镜阴性 + Rome III/IV；DiaSorin 进一步用 SES-CD / Mayo 评分定义 IBD 状态。
4. **提取装置成为独立 510(k)**：FED（K180971）、Q.S.E.T./Q.S.E.T. Plus（K182698/K213858）、CALEX Cap（K191718）、Easy Stool Extraction Device（K191807）。体积法取样重量 10–56 mg 不等；与手工称重相比斜率 0.88–1.15、120/160 处偏差 −14% 至 +10%，FDA 要求标签注明精密度差异。
5. **平台趋势**：手工 ELISA（2006–2018）→ 化学发光免疫分析仪（2017–2018）→ 生化分析仪 PETIA（2019 BÜHLMANN、2023 ALPCO），LoQ 11–30 µg/g 均低于最低 cutoff。
6. **儿科**：BÜHLMANN、DiaSorin 纳入 ≥4 岁；ALPCO 仅成人；PhiCal/Calprest 文件未载明。

### 10. 来源

- K050007 (DEN060001)：https://www.accessdata.fda.gov/cdrh_docs/reviews/K050007.pdf
- K130945：https://www.accessdata.fda.gov/cdrh_docs/reviews/K130945.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K130945.pdf
- K170993：https://www.accessdata.fda.gov/cdrh_docs/reviews/K170993.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K170993.pdf
- K180971：https://www.accessdata.fda.gov/cdrh_docs/reviews/K180971.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K180971.pdf
- K181012：https://www.accessdata.fda.gov/cdrh_docs/reviews/K181012.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181012.pdf
- K182698：https://www.accessdata.fda.gov/cdrh_docs/reviews/K182698.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K182698.pdf
- K190784：https://www.accessdata.fda.gov/cdrh_docs/reviews/K190784.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K190784.pdf
- K191718：https://www.accessdata.fda.gov/cdrh_docs/reviews/K191718.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K191718.pdf
- K213858：https://www.accessdata.fda.gov/cdrh_docs/reviews/K213858.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K213858.pdf
- K220763：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220763.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220763.pdf
- K232057：https://www.accessdata.fda.gov/cdrh_docs/reviews/K232057.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K232057.pdf

---

## D-二聚体（D-dimer）

### 1. 法规定位

| 项目 | 内容 |
|---|---|
| 21 CFR | **864.7320** — *Fibrinogen/fibrin degradation products assay*（K090264 另列 864.5425 用于质控品） |
| Class | Class II |
| Product code | **DAP** – Fibrinogen/Fibrin Split Products, Antigen, Antiserum, Control（含 VTE 排除声明的主流产品）；**GHH** – Fibrin Split Products（多为早期 FDP/D-dimer 及无排除声明产品，如 Diazyme K112120，其同时列 DAP、JIT）；相关：GGN（质控/线性物质）、JPA、NBC、LYR（多分析物或仪器组合） |
| Panel | 81 – Hematology |
| openFDA 累计清关数 | DAP total=58；GHH total=28；按 device_name 含「d-dimer」检索共 65 条（含 GGN 质控品 8 条、JPA 1、NBC 1、LYR 1） |

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K040882 | 2004-06-10 | bioMérieux / **VIDAS D-dimer Exclusion**（原型，在 K030328 DVT 排除基础上新增 PE） | ELFA | 柠檬酸血浆 | 与 PTP 模型联合**排除** DVT 和 PE（门诊） | 决策摘要 |
| K081732 → K091916 → K093626 | 2008-10-24 → 2009-10-29 → 2010-11-29 | Siemens / **INNOVANCE D-Dimer**（BCS/BCS XP） | 乳胶免疫比浊（单抗 8D3） | 血浆（柠檬酸） | K081732：「aid in diagnosis of VTE」；K091916 加 PE 排除；K093626 加 DVT 排除（non-high PTP） | 决策摘要（K081732、K093626） |
| K112818 | 2012-07-31 | bioMérieux / **VIDAS D-Dimer Exclusion II (DEX2)** | ELFA（20 min） | 柠檬酸 3.2%/3.8%、CTAD 血浆 | 排除 DVT/PE（门诊） | 决策摘要 + 510k summary |
| K090264 / K070927 / K151534 / K160885 / K172903 | 2010-02-05 / 2007 / 2015-07-06 / 2016-10-27 / 2017-11-22 | Instrumentation Laboratory / **HemosIL D-Dimer HS 500** 与 **HemosIL D-Dimer HS**（ACL TOP） | 乳胶免疫比浊（F(ab')2 单抗） | 柠檬酸血浆 | 排除 VTE（门诊）；HS cutoff 230 ng/mL、HS 500 cutoff 500 ng/mL；Special 510(k) 加入年龄校正文献说明 | K090264 决策摘要；其余 Special 510(k) 备忘 + summary |
| K110303 | 2011-05-16 | Siemens / **Stratus CS Acute Care D-dimer (DDMR)** | 固相 RPIA 荧光 | 柠檬酸或肝素血浆 | non-high PTP 排除 PE + VTE 辅助诊断；**临床实验室和 POC 场景** | 决策摘要 + 510k summary |
| K112120 | 2013-01-24 | Diazyme / **Diazyme D-Dimer**（Roche Modular P） | 乳胶免疫比浊 | 柠檬酸血浆 | **仅**「aid in detecting the presence of intravascular coagulation and fibrinolysis」，明确「Not for exclusion of DVT and PE」；GHH | 决策摘要 + 510k summary |
| K162227 | 2016-12-10 | Diagnostica Stago / **STA-Liatest D-Di**（STA-R、STA Compact、STA Satellite） | 微乳胶免疫比浊（8D2 + 2.1.16 双单抗，540 nm） | 3.2% 柠檬酸静脉血浆 | 新增 DVT 排除（原有 PE 排除） | 决策摘要 + 510k summary |
| K253985 | 2026-08-11 | Siemens / **INNOVANCE D-Dimer 2.0**（CS-5100） | 乳胶免疫比浊 | 3.2% 柠檬酸血浆 | 排除 DVT 和 PE（门诊，与 PTP 联合） | 决策摘要 + 510k summary |

### 3. 预期用途与声明类型

- **排除（exclusion）声明**：「for use in conjunction with a clinical pretest probability (PTP) assessment model to exclude DVT and PE in outpatients suspected of DVT or PE」——限定：(i) 门诊/急诊；(ii) 必须与 PTP 联合；(iii) Siemens/Stratus 写明「**non-high** clinical PTP」。这是 D-dimer 的核心声明，属「rule-out」型。
- **辅助诊断（aid in diagnosis）**：K081732 初次清关只写「aid in the diagnosis of VTE」；K110303 同时保留「aid in the diagnosis of VTE」+ PE 排除。
- **DIC/纤溶监测型**：K112120 Diazyme 仅「aid in detecting the presence of intravascular coagulation and fibrinolysis」，其 predicate Roche Tina-quant（K062203）声明含 DIC 监测 + non-high PTP 下 <0.5 µg FEU/mL 排除 DVT/PE，但 Diazyme 主动放弃排除声明，因此无临床研究、无临床 cutoff。
- **POC**：仅 Stratus CS（K110303）声明「clinical laboratory and point of care (POC) settings」，由受训医护人员使用；其余均为实验室凝血分析仪/免疫分析仪。
- 全部 prescription use。K253985 新增儿科限制：<10 岁性能未建立；<18 岁排除声明未验证。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

| 产品 | 临床 cutoff | 单位 | 建立方式（文件所述） |
|---|---|---|---|
| VIDAS Exclusion K040882 | **500 ng/mL FEU** | FEU | 决策摘要仅写「Assay cut-off: 500 ng/mL (FEU)」；用于前瞻性管理研究（<500 阴性、≥500 阳性）；建立过程文件未载明 |
| VIDAS DEX2 K112818 | 500 ng/mL FEU | FEU | 「transferred from the predicate device (K040882)」，在既往管理研究冻存样本中验证 |
| INNOVANCE K081732/K093626 | 0.50 mg/L FEU | FEU | 单中心前瞻性 359 例门诊疑似 VTE「derivation study」（敏感性 98%、特异性 38%、NPV 98%）确定 0.5；再在 902 例验证 |
| INNOVANCE 2.0 K253985 | 0.500 mg/L FEU | FEU | 与 predicate 相同；方法比对在 MDL 0.500 处偏差 0.0010 mg/L |
| HemosIL D-Dimer HS 500 K090264 | 500 ng/mL（FEU） | FEU（「Optimization of the House Standard to a cut-off value at 500 ng/mL FEU」） | ROC：295 例急诊冻存样本（75 VTE：47 PE、28 DVT）+ 100 例盲法研究；再由 4 医院 747 例管理研究验证 |
| HemosIL D-Dimer HS K070927/K151534/K160885 | **230 ng/mL** | 文件未标注单位（仅写「Cut-off 230 ng/mL」）；LoD 21→137 ng/mL；线性 150–69000 | 文件未载明（K070927 原文件未获取）。背景（非申报文件）：HemosIL D-Dimer HS 以 D-dimer units (DDU) 报告，230 ng/mL DDU ≈ 460 ng/mL FEU；DDU→FEU 换算系数约 2 |
| Stratus CS K110303 | **450 ng/mL [µg/L] FEU** | FEU | 「cutoff previously established and evaluated in K051597」 |
| STA-Liatest D-Di K162227 | 0.50 µg/mL FEU | FEU | 「The clinical cut-off was established at 0.50 µg/mL (FEU)」；沿用 K964728 |
| Diazyme K112120 | **无临床 cutoff**（Not applicable） | µg/mL FEU | 参考范围 <0.5 µg/mL FEU（90% 健康人） |

- **判读规则**：结果 < cutoff = 阴性（排除）；≥ cutoff = 阳性（需影像）。VIDAS 原型和 HemosIL 管理研究中「阴性 + 低/中 PTP → 不做进一步检查，随访 3 个月」；「阳性或高 PTP → 影像」。
- **年龄校正 cutoff**：K151534（HS）和 K172903（HS 500）为 Special 510(k)，仅在说明书「Summary and Principle」加入文献陈述：「D-Dimer levels also rise with age… age-adjusted cut-off values for DVT and PE suspicion have been shown to increase specificity… in patient populations greater than 50 years」，并在 Limitations 加入「**The performance of this assay has not been validated for use with age-adjusted cut-off values**」——即 FDA 允许提及文献但不允许声明。K253985（2026）文件未载明任何年龄校正内容。背景（非申报文件）：年龄校正公式为 年龄×10 µg/L FEU（>50 岁），来自 ADJUST-PE 等研究，未被任何本组 510(k) 采纳为声明。
- **参考人群（expected values）**

| 产品 | n / 定义 | 结果 |
|---|---|---|
| K040882 | 200 例献血者 | 96% <500 ng/mL FEU |
| K112818 | 215 例献血者，4 台 VIDAS、3 批试剂 | 90% (193/215) <500 ng/mL FEU |
| K093626 | 健康人 | <0.59 mg/L（n 文件未载明） |
| K253985 | 153 例 ≥18 岁健康成人（70 男/83 女，3 中心）+ 24 例儿童 | 非参数 EP28-A3c：<0.190–0.942 mg/L FEU；90th 百分位 0.468；儿科对照文献验证 |
| K090264 | 140 例健康献血者，C28-A2 | 95% 参考区间 29.1–500.1 ng/mL |
| K110303 | 131 例健康（81 男/50 女，18–59 岁） | 柠檬酸 <552 ng/mL FEU（38–804，均值 258）；Li-heparin <682（37–971，均值 304） |
| K162227 | K964728 建立 | 成人 <0.50 µg/mL FEU |
| K112120 | 120 例健康成人，C28-A3 | 90% <0.5 µg/mL FEU |

### 5. 生物学 / 生理学依据

- 申报文件：D-dimer 为「cross-linked fibrin degradation products」；INNOVANCE 说明 D-dimer 交联区具「stereosymmetrical structure」，单抗 8D3 表位出现两次，故单一抗体即可触发乳胶凝集（K081732、K253985）；Diazyme 510k summary：「Thrombus formation is normally followed by an immediate fibrinolytic response. The resultant generation of plasmin causes the release of fibrin degradation products (predominantly containing D-Dimer) into the circulation」；HemosIL/Roche 说明书：DVT、PE、DIC 升高，妊娠期升高，随年龄升高（K151534/K172903）。K040882 注明纤维蛋白原及 FDP 不干扰（比例 1:100 内），但「specificity… not tested against FDPE, cross-reactivity cannot be ruled out」。
- 背景（非申报文件）：D-dimer 由 FXIIIa 交联后的纤维蛋白经纤溶酶降解产生，是体内凝血激活并伴纤溶的标志；其排除价值基于高敏感性/高 NPV，而阳性预测值低（炎症、肿瘤、手术、妊娠、高龄均升高）；FEU 与 DDU 的差别源于校准物（纤维蛋白原当量 vs 纯 D-dimer 质量），FEU ≈ 2×DDU。

### 6. 样本类型与样本要求

| 产品 | 抗凝剂 | 基质研究 | 稳定性 |
|---|---|---|---|
| K040882 VIDAS | 柠檬酸血浆 | 文件未载明 | 文件未载明 |
| K112818 DEX2 | 3.2% 或 3.8% 柠檬酸；CTAD | 柠檬酸 vs CTAD：124 例献血者 + 18 人工高值，500 处偏差 −0.8% (−4.2%; 2.6%)；3.2% vs 3.8%：91 例 + 18 加标，500 处偏差 −1.7% | 试剂 12 个月；样本稳定性以冻存 vs 新鲜比对支持临床研究 |
| K081732/K093626 INNOVANCE | 柠檬酸血浆 | 文件未载明 | 文件未载明 |
| K253985 INNOVANCE 2.0 | **3.2% 柠檬酸静脉血浆** | Matrix comparison「Not applicable」 | 38 份样本 2 中心：15–25°C 4 h（原管或分杯）；2–8°C 24 h；≤−18°C 3 个月；≤−74°C 6 个月；冻融 1 次后 15–25°C 4 h；试剂 18 个月、开瓶 31 天、机上 120 h |
| K090264 HemosIL HS 500 | 柠檬酸血浆 | 新鲜 vs 冻存 57 例 r 0.999；3.2% vs 3.8% 柠檬酸 20 对，斜率 0.998、偏差 −2.3% | 试剂 21 个月（进行中） |
| K110303 Stratus CS | 柠檬酸或肝素血浆 | 见 K022976；两种基质各做独立临床研究 | 文件未载明 |
| K162227 Stago | 3.2% 柠檬酸 | 见 K964728 | 见 K964728 |
| K112120 Diazyme | 柠檬酸血浆 | 不适用 | 2–8°C 4 天；−20°C 3 个月；1 次冻融（10 例 0.32–5.01 µg/mL） |

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 对应项目 |
|---|---|---|
| K040882 | （空白，未列） | 精密度（CV 3.9–7.1%）、线性 45–10,000 ng/mL FEU、LoD ≤45、纤维蛋白原/FDP 特异性 |
| K081732 | EP7-A2；EP5-A2；EP9-A2；EP6-P2；NCCLS/CLSI EP17-A；FDA *Assayed and Unassayed Quality Control Material* 指南 | LoB 0.02、LoD 0.05 mg/L；线性 0.12–4.46；干扰（含 RF 3640 IU/mL，回收 92–101%） |
| K093626 | EP5A2；EP6A2；EP7A2；EP9A2；EP17A | 分析性能引用 K081732 |
| K112818 | EP5-A2；EP6-A；EP7-A；EP9-A2；**EP12-A2**；EP17-A | 3 中心精密度（CV ≤7%）；LoB 7.4 / LoD 16.5 / LoQ 29 ng/mL FEU；47 种药物干扰 |
| K090264 | EP05-A2；EP07-A2；C28-A2；EP09-A2；EP17-A | LoB 146 / LoD 203 / LoQ 208 ng/mL；线性 215–128,000（auto-rerun）；干扰 Hb 500 mg/dL、胆红素 18、TG 1327、RF 1400 IU/mL、FDP 10 µg/mL、HAMA |
| K110303 | Not applicable | 分析性能引用 K022976/K051597/K063356 |
| K112120 | EP5-A2；EP6-A；EP7-A；EP9-A2；EP17-A；C28-A3 | LoB 0.06 / LoD 0.09 / LoQ 0.15 µg/mL FEU；线性至 8.6；hook 至 30 µg/mL；3 中心再现性 |
| K162227 | **CLSI H59-A** *Quantitative D-dimer for the Exclusion of Venous Thromboembolic Disease* | 仅临床研究；分析性能引用 K964728/K983460/K082248 |
| K253985 | EP05-A3；EP06-A2 (2020)；EP07 (2018)；EP09c (2018)；**EP12-Ed3 (2023)**；EP14-A3（互通性）；EP17-A2；EP25-A2 (2023)；EP28-A3c；EP34（扩展测量区间）；EP39（替代样本）；**H21-A5**（凝血样本采集处理）；**H59-A** | LoB 0.0277 / LoD 0.0416 / LoQ 0.135 mg/L FEU；AMR 0.190–80.000 |
| K151534 / K160885 / K172903 | ISO 14971；EP17-A2（K160885 LoD 更新） | Special 510(k)，无新性能研究 |

### 8. 临床验证设计与结果

**K040882 VIDAS Exclusion（原型）**
- DVT：3 中心前瞻性队列冻存样本，556 例连续门诊首次疑似 DVT，Wells 模型分低/中/高 PTP；阴性 + 低/中 PTP 不再检查、随访 3 个月；阳性或高 PTP 做系列加压超声。n=555：敏感性 100% (93.6–100)、特异性 32.9% (28.8–37.2)、NPV 100% (97.8–100)；低 PTP 295：特异性 39.7%、NPV 100% (96.7–100)；中 PTP 189：特异性 26.7%、NPV 100% (92.3–100)；高 PTP 71：特异性 16.0%、NPV 100% (63.1–100)。
- PE：3 中心前瞻性，965 例急诊疑似 PE，新鲜样本，「WICKIE model」（文件原文；背景：即 Wicki/Geneva 评分）分层；阴性者不治疗不检查；阳性者超声/螺旋 CT/血管造影，随访 3 个月。全人群：敏感性 100% (98.4–100)、特异性 37.7%、NPV 100% (98.7–100)、PPV 32.4%；低+中 PTP n=891：敏感性 100% (97.7–100)、NPV 100% (98.7–100)、PPV 25.8%；高 PTP 74：PPV 91.3%。

**K081732 INNOVANCE（辅助诊断 → 排除的基础数据）**：902 例急诊门诊疑似 VTE（533 女/369 男，18–95 岁）冻存样本，2 采集中心 + 2 检测中心，Wells 评分，标准客观检查确诊，阴性者随访 4 个月；VTE 频率 22.0%；cutoff 0.5 mg/L：敏感性 97%、特异性 42%、NPV 98%。方法比对 vs Stratus CS：2 中心 318 份，Passing-Bablok 斜率 0.951、截距 0.059、r 0.97。

**K093626 INNOVANCE DVT 排除**：多中心前瞻性，455 例连续急诊疑似 DVT（29 例排除 → 426），Wells likely/unlikely；阳性者加压超声/静脉造影，阴性者及影像阴性者随访 3 个月；DVT 患病率 21.8% (93/426)。全体：敏感性 100% (96.1–100)、特异性 34.5% (29.4–39.9)、NPV 100% (96.8–100)；unlikely PTP n=267：敏感性 100% (83.9–100)、特异性 37.0%、NPV 100% (96.0–100)。方法比对 vs VIDAS n=265（0.17–4.17 mg/L FEU）：斜率 1.11、截距 −0.075、r 0.96。

**K253985 INNOVANCE 2.0（2026）**：冻存样本，急诊 ≥18 岁，1246 例疑似 PE + 1263 例疑似 DVT，Wells likely/unlikely；阳性者影像，阴性者随访 3 个月。DVT 患病率 6.2% (78/1263)：敏感性 97.4% (LCL 92.1)、特异性 48.4% (LCL 46.0)、NPV 99.7% (LCL 99.0)、PPV 11.1%（2 例假阴性）。PE 患病率 7.8% (97/1246)：敏感性 99.0% (LCL 95.2)、特异性 57.4% (LCL 55.0)、NPV 99.8% (LCL 99.3)、PPV 16.4%（1 例假阴性）。方法比对 vs INNOVANCE：377 例（13 操作者，3 美国 + 1 境外中心；59.4% 新鲜/40.6% 冻存；8% 为 10–18 岁儿科；含疑似 PE/DVT/DIC），Passing-Bablok 斜率 1.036 (1.024–1.052)、截距 −0.017、r 0.995、MDL 偏差 0.0010 mg/L FEU。

**K112818 VIDAS DEX2**：既往 VTE 管理研究的冻存样本 n=315（VTE 患病率 23.5%），客观检查 + 3 个月随访确诊；低+中 PTP n=303：DEX2 敏感性 100% (94.2–100)、特异性 35.7% (29.6–42.1)、NPV 100% (95.8–100)、PPV 28.6%；旧 VIDAS 同批：特异性 37.8%。方法比对 vs 旧 VIDAS n=326：250 处偏差 +5.0%、500 处 +12.2% (8.5–16.0)、3000 处 +18.2%。

**K090264 HemosIL D-Dimer HS 500**：
- cutoff 建立：295 例急诊冻存样本 ROC（VTE 75 例）；盲法研究 100 例（HS 500 敏感性 100% (28/28)、特异性 31.9%；VIDAS 特异性 25.0%）。
- 多中心管理研究（4 医院，747 例连续急诊：DVT 疑似 401 例，平均 65.9 岁；PE 疑似 346 例，平均 51 岁）：Wells；阴性 + 低 PTP 不再检查、3 个月随访；阴性 + 中 PTP 由医生决定随访或影像；阳性或高 PTP 影像。DVT 患病率 22.4%：全体敏感性 100% (90/90; 96–100)、特异性 42.1% (36.6–47.8)、NPV 100% (97.2–100)；低+中 PTP 322：敏感性 100% (45/45)、特异性 43.3%、NPV 100% (97–100)。PE 患病率 15%：全体敏感性 100% (52/52; 93.2–100)、特异性 48.3%、NPV 100% (97.4–100)；低+中 PTP 322：敏感性 100% (43/43)、特异性 49.1%、NPV 100% (97.3–100)。3 个月随访中无 D-dimer 阴性者发生 VTE。
- 结局研究 295 例：HS 500 敏感性 100% (95.2–100)、特异性 42.3%、NPV 100%；VIDAS 特异性 35%。

**K162227 STA-Liatest D-Di（DVT 排除）**：前瞻性，16 中心（美国 9，法/意/西/加 7），1219 例急诊或门诊疑似 VTE；Wells 评分；**高 PTP 158 例排除**、其他排除 96 → 980 例分析；影像 + 3 个月随访。DVT 患病率 8.4% (85/980)（美国 6.0%，欧加 10.3%）；敏感性 100% (85/85; LCL 95.8%)、特异性 55.2% (LCL 51.9%)、NPV 100% (494/494; LCL 99.3%)、PPV 17.5%。美国亚组 369：敏感性 100% (LCL 84.6)、特异性 58.2%。510k summary 引述既往 PE 排除研究：低+中 PTP n=1130，敏感性 97.0% (91.6–99.4)、NPV 99.7% (99.2–100)。

**K110303 Stratus CS（POC，PE 排除）**：
- 柠檬酸血浆：730 例连续急诊疑似 PE（排除 75 → 655），Wells 高/中/低；PE 患病率 14.0%；cutoff 450 ng/mL FEU：敏感性 98.9% (94.1–100)、特异性 42.5% (38.3–46.7)、NPV 99.6% (97.7–100)；低+中 PTP 625：98.7% / 43.0% / 99.6%。
- 肝素血浆：468 → 427 例；患病率 14.1%；敏感性 98.3% (91.1–100)、特异性 29.7%、NPV 99.1% (95.04–100)；低+中 401：97.9% / 29.9% / 99.1%。
- 方法比对 vs INNOVANCE n=396（54–4506 ng/mL）：斜率 0.950、截距 −12.32、r 0.938。

**K112120 Diazyme（DIC 监测型）**：临床敏感性/特异性、临床 cutoff 均「Not applicable」；仅方法比对 vs Roche Tina-quant：88 例独特样本/128 次测定（ICU、创伤、术后患者），bootstrap 分析，斜率 1.0±0.1、r² >0.90、截距 ±0.15 达标。

### 9. 厂家间差异与要点

1. **排除声明 vs DIC 声明的分水岭**：有排除声明者必须提交前瞻性（或前瞻性队列冻存样本）管理研究——连续急诊/门诊疑似 VTE 患者、Wells（或 Geneva/Wicki）PTP 分层、影像金标准 + 3 个月（K081732 为 4 个月）随访、报告敏感性/NPV 及其 95% CI 下限；无排除声明者（Diazyme K112120）只需方法比对 + 参考区间，且被要求在标签中写明「Not for exclusion of DVT and PE」。
2. **NPV 要求**：文件未写明固定的 NPV 接受阈值（文件未载明），但所有清关的排除声明研究 NPV ≥99.1%（2011 Stratus 肝素）～100%，敏感性 CI 下限 ≥83.9%（小亚组）；2026 年 K253985 报告 LCL（单侧下限）而非双侧 CI。
3. **单位与 cutoff**：主流为 500 ng/mL (0.5 mg/L, 0.5 µg/mL) FEU；Stratus CS 450 ng/mL FEU；HemosIL D-Dimer HS 230 ng/mL（文件未标单位；背景：DDU）。IL 公司同时维持 HS（230）和 HS 500（500 FEU）两条产品线，K090264 明确将 HS 500 的内部标准「optimized… to 500 ng/mL FEU」以与 VIDAS 对齐，说明 cutoff 一致性通过校准物赋值而非改变临床阈值实现。
4. **患病率随时间下降**：2004–2010 研究 VTE 患病率 14–25%；2016 Stago 8.4%；2026 Siemens 6.2%（DVT）/7.8%（PE）→ 特异性从 33–42% 升至 48–57%，NPV 维持 ≥99.7%。
5. **PTP 限定措辞**：bioMérieux/IL/Stago 写「in conjunction with a PTP assessment model」；Siemens 写「**non-high** PTP」；Stago 研究直接排除高 PTP 患者；VIDAS 原型和 HemosIL 研究纳入高 PTP 但分层报告。
6. **年龄校正**：仅 IL 通过 Special 510(k) 在说明书中提及文献，并声明「not validated」；无任何产品获得年龄校正 cutoff 声明（截至本组文件）。
7. **抗凝剂**：几乎全部限定柠檬酸（3.2%，部分兼容 3.8%、CTAD）；Stratus CS 与 Roche Tina-quant 兼容肝素血浆（Stratus 肝素基质特异性明显低于柠檬酸：29.7% vs 42.5%）。
8. **标准演进**：2004 年无标准引用 → 2008–2013 EP5/6/7/9/17 → 2016 起引用 CLSI **H59-A**（D-dimer 排除 VTE 专用指南）→ 2026 年扩展至 EP12-Ed3、EP14、EP25-A2、EP34、EP39、H21-A5。

### 10. 来源

- K040882：https://www.accessdata.fda.gov/cdrh_docs/reviews/K040882.pdf
- K081732：https://www.accessdata.fda.gov/cdrh_docs/reviews/K081732.pdf
- K093626：https://www.accessdata.fda.gov/cdrh_docs/reviews/K093626.pdf
- K112818：https://www.accessdata.fda.gov/cdrh_docs/reviews/K112818.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K112818.pdf
- K090264：https://www.accessdata.fda.gov/cdrh_docs/reviews/K090264.pdf
- K151534：https://www.accessdata.fda.gov/cdrh_docs/reviews/K151534.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K151534.pdf
- K160885：https://www.accessdata.fda.gov/cdrh_docs/reviews/K160885.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K160885.pdf
- K172903：https://www.accessdata.fda.gov/cdrh_docs/reviews/K172903.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K172903.pdf
- K110303：https://www.accessdata.fda.gov/cdrh_docs/reviews/K110303.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K110303.pdf
- K112120：https://www.accessdata.fda.gov/cdrh_docs/reviews/K112120.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K112120.pdf
- K162227：https://www.accessdata.fda.gov/cdrh_docs/reviews/K162227.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162227.pdf
- K253985：https://www.accessdata.fda.gov/cdrh_docs/reviews/K253985.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253985.pdf
