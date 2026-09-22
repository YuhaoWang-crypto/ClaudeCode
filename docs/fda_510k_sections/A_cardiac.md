# 靶点组 A：心脏标志物（cTnI/cTnT、NT-proBNP/BNP、hsCRP）FDA 510(k) 申报信息整理

> 数据来源说明：本文所有申报数据均取自 FDA 决策摘要（decision summary，`cdrh_docs/reviews/<K>.pdf`）及 510(k) summary（`cdrh_docs/pdfYY/<K>.pdf`）原文；openFDA `device/510k` API 用于清单与累计数。文件未写明的信息一律标注「文件未载明」。凡来自本人背景知识（指南、生理学）的内容均以「背景（非申报文件）」标出。
> 检索日期：2026-09-22。

---

## 心肌肌钙蛋白 I / T（cardiac troponin, cTnI / cTnT；hs-cTn）

### 1. 法规定位
- Product code：**MMI**（"Immunoassay method, troponin subunit"）
- 21 CFR：**862.1215**（Creatine phosphokinase/creatine kinase or isoenzymes test system）
- Class：**II**；Panel：Clinical Chemistry (75 / CH)
- openFDA 该代码累计清关数：**total = 95**（`list MMI` 首行）
- 配套 calibrator 走 JIT (862.1150)，control 走 JJX/JJY (862.1660)（K162895 载明）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K253051 | 2026-06-11 | Abbott Point of Care / i-STAT hs-TnI cartridge + i-STAT Alinity | 卡式 ELISA + 电化学检测（≈15 min，22 µL） | Li-heparin 全血 / 血浆（及无抗凝全血，3 min 内） | aid in the diagnosis of MI；POC 或临床实验室 | 决策摘要 + 510(k) summary |
| K252393 | 2025-10-29 | Ortho / VITROS Immunodiagnostic Products hs Troponin I（VITROS 5600） | 化学发光（HRP–luminol）免疫测定 | 肝素（Li-heparin）血浆 | aid in the diagnosis of MI | 决策摘要 + 510(k) summary |
| K231974 | 2024-03-20 | PHC / PATHFAST hs-cTnI-II | 化学发光酶免疫（Magtration，POC 台式） | 肝素或 EDTA 全血 / 血浆 | aid in the diagnosis of AMI；临床实验室或 POC | 决策摘要（极简，多处"见 K100130"）+ 510(k) summary |
| K191595 | 2019-09-13 | Abbott / ARCHITECT STAT High Sensitivity Troponin-I（i2000SR） | CMIA（两步法） | K2-EDTA 血浆 | aid in the diagnosis of MI | 决策摘要 + 510(k) summary |
| K172783 | 2018-06-12 | Beckman Coulter / Access hsTnI（UniCel DxI 800） | 顺磁微粒化学发光（ALP） | 血清、Li-heparin 血浆 | aid in the diagnosis of MI | 决策摘要 + 510(k) summary |
| K171274 | 2018-07-12 | Siemens / ADVIA Centaur High-Sensitivity Troponin I (TNIH) | 三位点夹心直接化学发光（acridinium ester） | 血清、Li-heparin 血浆 | aid in the diagnosis of AMI | 决策摘要 + 510(k) summary |
| K162895 | 2017-01-18 | Roche / Elecsys Troponin T Gen 5 STAT（cobas e 411/601） | ECLIA（9 min STAT） | Li-heparin 血浆（仅此一种） | aid in the diagnosis of MI；明确**不再**声明 ACS 风险分层/慢性肾衰风险 | 决策摘要 + 510(k) summary |
| K201441 | 2021-09-21 | Roche / Elecsys Troponin T Gen 5（cobas e 801，加入 biotin scavenger 抗体） | ECLIA | Li-heparin 血浆 | aid in the diagnosis of MI（cutoff 与 K162895 相同） | 决策摘要 + 510(k) summary |
| K063243 | 2007-12-14 | bioMérieux / VIDAS Troponin I Ultra（"原型"，非 hs） | ELFA（酶联荧光） | 血清、Li-heparin 血浆 | aid in diagnosis of MI（去掉了 predicate 的 ACS 风险分层声明） | 仅决策摘要（pdf07 无 summary） |

（另下载了 K240984 i-STAT 1 hs-TnI（K253051 的 predicate）与 K121790 Access AccuTnI，未单独展开。）

### 3. 预期用途与声明类型
- 全部为**辅助诊断 MI/AMI**（"aid in the diagnosis of myocardial infarction"），**处方用**（Rx only），体外诊断用。
- 明确的**声明收窄**：Roche K162895 与 predicate K051752（4 代 TnT STAT）比较表写明，新一代 hs 产品"Not indicated for these uses"——即不再声明 ACS 患者风险分层、慢性肾衰心脏风险、强化治疗选择。VIDAS TnI Ultra（K063243）也去掉了 predicate（Dade Dimension cTnI K010313）的"risk stratification ... relative risk of mortality"声明。PATHFAST 510(k) summary 直接注明"Not for risk stratification"。
- **POC 声明**：i-STAT hs-TnI（K253051/K240984）与 PATHFAST hs-cTnI-II（K231974）标注"for use in clinical laboratory or point of care (POC) settings"，样本含全血。其他为中央实验室自动化分析仪。
- 定量报告：新产品统一以 ng/L 报告（PATHFAST K231974 的申报目的之一就是把单位从 ng/mL 改为 ng/L 并把可报告范围下限降到 4.1 ng/L）。
- 说明书通用限制（多份文件重复出现）：troponin 对心肌坏死特异但对 MI 不特异，需结合 ECG、症状及**系列采样的升/降趋势**（Universal Definition of MI）解读；troponin 自身抗体可见于 10–20% 的 ED 患者并可致假性偏低（K191595、K172783、K171274、K253051）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**hs 产品的 cutoff 全部 = 表面健康人群的 99th percentile URL**（总体 + 性别特异三套值并列写入说明书）。逐产品：

| 产品 | 女性 | 男性 | 总体 | 参考人群（n、定义）| 统计方法 |
|---|---|---|---|---|---|
| Roche Elecsys TnT Gen 5 STAT（K162895/K201441） | 14 ng/L | 22 ng/L | 19 ng/L | n=1301（645 男/656 女），4 个外部采集点；排除癌症、任何慢性病或服用慢性病处方药、高血压、ACS 史、3 个月内住院、妊娠或产后 6 周内 | 决策摘要仅写"99th percentile upper reference limit"，具体统计方法（参数/非参数）文件未载明；CI 文件未载明 |
| Abbott ARCHITECT hs-TnI（K191595） | 17 [90% CI 14–20] | 35 [27–44] | 28 [22–33] ng/L | n=1531（765 女/766 男），美国人群，21–75 岁；入组生物标志物筛查：BNP ≤25（男）/≤40（女）pg/mL，HbA1c ≤6%，GFR ≥60 mL/min/1.73 m²；冻存后单次测定 | CLSI EP28-A3c 的 **robust statistical method**（注意：非非参数法） |
| Beckman Access hsTnI（K172783） | 血浆 14.9 (95% CI 10.1–27.1)；血清 13.6 (10.0–25.6) pg/mL | 血浆 19.8 (15.9–38.4)；血清 19.8 (15.4–44.8) | 血浆 17.9 (14.7–27.1)；血清 18.1 (14.3–25.6) | 1088 份 Li-heparin 血浆（593 女/495 男）+1085 份血清，5 个地理分散的采集点，"no known diseases of the cardiovascular system or other serious acute or chronic diseases or infections" | non-parametric empirical univariate distribution function |
| Siemens ADVIA Centaur TNIH（K171274） | Li-hep 36.99；血清 39.59 pg/mL | Li-hep 57.27；血清 58.05 | Li-hep 47.34；血清 46.47；**说明书总体值取 47.34** | n=1990（1006 女/984 男），美国，22–91 岁，同上健康定义；两名女性（≈400、≈5000 pg/mL）作为 outlier 剔除，申办方证明剔除不改变 cutoff | non-parametric empirical univariate distribution function；CI 文件未载明 |
| Ortho VITROS hs Troponin I（K252393） | 9 (90% CI 3.9–17.5) | 12 (8.8–20.9) | 11 (8.2–14.3) ng/L | n=952（486 女/466 男），22–91 岁，59% ≥50 岁；排除肾病、糖尿病、心脏病、癌症、肺病、甲状腺病、卒中史，高血压/高脂，3 个月内肌肉骨骼损伤或手术，吸烟，妊娠；加测 HbA1c ≥6.5%、NT-proBNP >125（<75 岁）/ >450（≥75 岁）pg/mL、eGFR <60 者排除 | non-parametric（EP28-A3c） |
| Abbott i-STAT hs-TnI（K253051，沿用 K240984 值） | 13 (90% CI 10–17) | 28 (19–58)* | 21 (14–30) ng/L | n=896 入组/895 分析（490 女/404 男），美国；筛查：NT-proBNP <125（<75 岁）/ <450（≥75 岁）pg/mL，HbA1c ≤6.5%，GFR ≥60；排除 BMI >35 或 <16、1/2 型糖尿病、3 个月内住院、心血管病史（含需服药的高血压、MI、心绞痛）、PCI/CABG/血管手术、6 个月内他汀、妊娠/产后 6 周 | nonparametric method（EP28-A3c）；*男性值剔除 1 例 outlier；在 Alinity 上实测 F 12 / M 29 / 总 21，判定与 i-STAT 1 相似故沿用 |
| PHC PATHFAST hs-cTnI-II（K231974） | 文件未载明 | 文件未载明 | 29 ng/L（K100130 建立，非参数） | 参考人群 n、纳排：本文件未载明（"See K100130"）；文件称 >50% 健康人可测出 >LoD | non-parametric |
| bioMérieux VIDAS TnI Ultra（K063243，非 hs） | — | — | **cutoff 0.11 ng/mL = 10% CV 处 LoQ**，而非 99th | 399 份美国健康人肝素血浆 99th = 0.01 ng/mL（即 LoD 水平，低于 cutoff） | 文件未载明方法 |

**cutoff 建立逻辑的演变（申报文件可见）**：
- 旧一代（VIDAS 2007）：健康人 99th 落在 LoD 附近无法用，故取 **10% CV 浓度（功能灵敏度）**作为 cutoff；Roche 4 代 TnT（K051752，作为 predicate 被引述）cutoff 0.1 ng/mL 基于 **ROC**。
- hs 一代（2017 起）：cutoff = 99th URL；FDA 要求同时评估**性别特异**与**总体** cutoff 在临床队列中的敏感性/特异性/PPV/NPV，并在说明书写入"女性用较低的性别特异 cutoff 时 PPV 更低、阳性中非 MI 比例可高达 xx%"（K162895：69/82/78%；K191595：87.6%；K172783：75%；K171274：71%）以及"男性用较高的性别特异 cutoff 会增加假阴性"（K171274：最多 +2.9%；K253051：+2.00%）。
- 灰区：troponin 申报文件未设置灰区；解读依赖系列采样的 rise/fall（说明书语言）。
- **hs 定义**（K252393、K171274 510(k) summary 引 IFCC TF-CB）：① 99th 处 %CV ≤10%；② ≥50% 健康男性及 ≥50% 健康女性可测出 >LoD 的浓度。K252393 用 LoQ 研究的不精密度曲线回归估算 99th 处 CV：女 4.7%、男 4.5%、总体 4.5%。

### 5. 生物学 / 生理学依据
- **申报文件所载**：troponin 在心肌细胞坏死时释放，"cardiac specific but not specific for MI"，可见于心律失常、主动脉综合征、急性心衰、高血压危象、心肌炎、心包炎、肺栓塞、Takotsubo（K162895 说明书节选）；Ortho K252393 临床试验中 15.79% (311/1970) 非 MI 者至少一次结果 >11 ng/L，其中 91.96% 有 eGFR <60 或心绞痛/房颤/心肌病/CAD/心衰/心动过速。判定标准为 Universal Definition of MI（2007 版：K162895、K191595、K172783、K171274；第三版：K171274 510(k) summary；第四版：K252393、K253051，要求 rise and/or fall + 缺血证据 + 系列采样）。
- **背景（非申报文件）**：cTnI/cTnT 为心肌肌钙蛋白复合体亚单位，胞质游离池早期释放、结构池持续释放；hs 检测把可检出下限推到健康人分布内，使 0/1 h、0/2 h 快速排除策略成为可能（ESC 指南），但 FDA 510(k) 说明书仅按时间窗报告诊断性能，不清关具体的 0/1 h 算法阈值。

### 6. 样本类型与样本要求
| 产品 | 样本 | 稳定性 / 基质研究（文件所载） |
|---|---|---|
| Roche TnT Gen 5（K162895） | **仅 Li-heparin 血浆**（predicate 允许血清+血浆） | 2–8 °C 24 h；−20 °C 12 个月；仅冻融 1 次；Hb >0.1 g/dL 假性降低；高剂量 biotin（>5 mg/d）停药 8 h；无基质比对（"intended for lithium heparin only"） |
| Abbott ARCHITECT hs-TnI（K191595） | 仅 K2-EDTA 血浆 | 室温 8 h（因此 20 天精密度用质控物完成）；fibrinogen 1000 mg/dL 干扰 +11.6%；总蛋白 12.4 g/dL −12~−18% |
| Beckman Access hsTnI（K172783） | 血清 + Li-heparin 血浆 | 室温 4 h、2–8 °C 48 h、−20 °C 3 天；两种基质各自完成精密度/线性/LoQ/临床/99th（未做单独基质等效研究）；温度补偿算法，18–30 °C 残余偏倚 ≤8%；ALP >400 U/L 假阳性、asfotase alfa 患者禁用 |
| Siemens Centaur TNIH（K171274） | 血清 + Li-heparin 血浆 | 所有性能在两种基质分别完成，无独立基质比对；样本稳定性文件未载明 |
| Ortho VITROS hs TnI（K252393） | 仅 Li-heparin 血浆（predicate 含血清/EDTA） | 临床样本 −20 °C 冻存并提供稳定性证据；biotin 0.351 mg/dL 无干扰；样本稳定性细节文件未载明 |
| Abbott i-STAT hs-TnI（K253051） | Li-heparin 全血、血浆；无抗凝全血 | Li-hep 全血/血浆室温 4 h；无抗凝全血 **3 min**；全血 vs 血浆基质等效按 **CLSI EP35**（临床点配对样本，双份测定）；无抗凝 vs Li-hep 全血 n=86；HCT 15–55% 准确，≥55% 不精密度 >10%、偏倚 ±10%，HCT >55% 只测血浆；海拔 ≈10,000 ft 研究 |
| PATHFAST hs-cTnI-II（K231974） | EDTA / Li-heparin 全血与血浆 | 四种基质分别测 LoB/LoD/LoQ（LoQ 均 4.1 ng/L）；稳定性"unchanged from K100130" |
| VIDAS TnI Ultra（K063243） | 血清、Li-heparin 血浆 | 37 对样本 Passing-Bablok：serum = 0.95 [0.91;0.97]·plasma − 0.01；说明书建议系列采样用同一管型 |

### 7. 分析性能验证所依据的标准（文件 "Standards/Guidance Documents Referenced" 逐项）

| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K162895 (Roche 2017) | CLSI EP5-A2；EP6-A；EP17-A2 | 精密度（21 天、3 批、2 台）；线性（21 点，≤6.3%/12.4% 偏离）；LoB/LoD/LoQ（LoQ = ≤20% CV 精密度曲线；另提供 11 ng/L 处 CV=10% 数据） |
| K201441 (Roche 2021) | 文件"Standards"栏未逐项读取；目的为降低 biotin 干扰 | 精密度、biotin 干扰（血清 biotin 可达 355 ng/mL）、临床 |
| K191595 (Abbott 2019) | EP05-A2；EP06-A；EP07-A2；EP17-A2；EP28-A3c | 精密度（EP05 20 天用质控）；线性 10 样本 ≤14.1%；内源/药物/交叉反应干扰；LoB/LoD/LoQ（parametric LoB/LoD，LoQ 精密度曲线 ≤20%）；参考区间（robust 法） |
| K172783 (Beckman 2018) | EP05-A3；EP06-A；EP07-A2；EP17-A2；**EP25-A**；EP28-A3c | 精密度（含 3 温度 × 3 温度校准矩阵、多点重现性、模拟现场精密度）；线性（加权回归）；干扰（判定 ±10% 或 ±2.3 pg/mL @≤11.5 pg/mL）；检出限（LoB 非参数）；试剂稳定性；参考区间 |
| K171274 (Siemens 2018) | EP05-A3；EP06-A；EP07-A2；EP17-A2；EP28-A3c | 精密度（20 天，2 台，80 重复；3 批；3 点重现性）；线性（0–150、2–25,000 两段）；干扰含 HAMA/RF；LoB 非参数 240 测定/基质/批，LoD 640 测定；参考区间 |
| K252393 (Ortho 2025) | EP05-A3；EP17-A2；EP06 2nd；EP07 3rd；EP37；EP28-A3c；**ISO 17511:2021** | 精密度；检出限（LoB 非参数、LoD 参数、LoQ 幂函数拟合 20% CV）；线性；干扰；参考区间；溯源（内部标准） |
| K253051 (Abbott POC 2026) | EP05-A3；EP06 2nd；EP07 3rd；**EP09c**；**EP12**；EP17 2nd；EP28-A3c；**EP35**；EP37；IEC 61326-1 | 精密度（含 POC 操作者多点研究）；线性（全血/血浆各 11 样本）；干扰（全血/血浆）；方法学比对（vs i-STAT 1）；定性一致性；LoB/LoD（非参数）/LoQ（≤20% CV 精密度曲线）；参考区间；样本类型等效（全血 vs 血浆）；EMC |
| K231974 (PATHFAST 2024) | EP6 2nd；EP17-A2 | 线性（12 水平，≤9.7% 偏离）；LoB 非参数、LoD 参数、LoQ <20% CV |
| K063243 (VIDAS 2007) | EP5-A2；EP6-A（列两次）；LoB/LoD"protocol similar to EP17-A" | 精密度（10 天、2 批、3 点，n=240）、LoQ@10% CV；回收率法线性；LoB = 零标准均值 + 2SD |

### 8. 临床验证设计与结果（逐产品给数字）

- **Roche K162895**：① APACE（国际多中心前瞻性，ED 胸痛/心绞痛，症状高峰 ≤12 h，唯一排除 = 透析；独立裁定委员会含心脏科医生，2007 UDMI，60 天随访；n=718，3 h 515 例，6 h 310 例）。19 ng/L：基线 Sens 93.5% (115/123, 87.6–97.2)，Spec 86.4% (514/595)，PPV 58.7%，NPV 98.5%；3 h Sens 98.3%，NPV 99.7%；6 h Sens 100%，NPV 100%。女 14 ng/L 基线 Sens 97.1% (34/35)，Spec 77.4%，PPV 41.5%；男 22 ng/L 基线 Sens 90.9% (80/88)，Spec 89.3%，PPV 66.1%。② 美国多中心 n=1679（173 例裁定 MI；排除 3 个月内 MI/手术/PCI 等，FDA 注明特异性和 PPV 可能被高估）：cobas e 411，19 ng/L 基线 Sens 86.7% (143/165)，Spec 87.8%，NPV 98.3%；3 h Sens 94.4%，NPV 99.3%。
- **Roche K201441**（cobas e 801 重做美国 ED 队列，症状 ≤12 h，冻存样本）：19 ng/L 男性 <1.5 h Sens 81.5%、Spec 85.7%、NPV 94.7%；>2.5–3.5 h Sens 98.5%、NPV 99.7%；女性 <1.5 h Sens 74.3%、Spec 86.4%、NPV 96.3%。
- **Abbott K191595**：11 家 ED，n=1065，MI 10.8% (116/1065)，3 名认证心脏科医生按 2007 ACC/AHA/ESC 裁定，盲法，冻存样本。总体 28 ng/L：女基线 Sens 91.7% (22/24)、Spec 92.0%、NPV 99.4%；男基线 Sens 81.8% (54/66)、Spec 81.5%、NPV 96.9%；男 2–4 h Sens 91.7%。性别特异：女 17 ng/L 基线 Sens 95.8%、Spec 87.6%、PPV 32.4%；男 35 ng/L 基线 Sens 78.8%、Spec 84.5%，4–9 h Sens 93.7%。说明书：PPV 下限女 12.4%、男 33.7%。
- **Beckman K172783**：2010–2011 年多中心，ED 胸痛/缺血等症状 n=1929 入组、1854 可测，2007 UDMI 独立专家裁定，盲法，0/1–3/3–6/6–9 h。血浆总体 17.9 pg/mL：基线 Sens 88.1% (89/101)、Spec 88.5%、PPV 57.8%、NPV 97.7%；1–3 h 94.1%/89.8%/99.2%；3–6 h 94.1%/89.8%/99.1%；6–9 h 98.6%/85.1%/99.7%。血清总体 18.1：基线 87.3%/89.3%，6–9 h 97.1%/85.0%。女性血清 13.6：6–9 h Sens 100% (20/20)，PPV 37.7%。
- **Siemens K171274**：n=2495（29 个美国采集点、3 个检测点），MI 13% (329/2495)，独立专家裁定，0–<1.5 h 至 ≥24 h 共 8 个窗。血浆总体 47.34 pg/mL：0–<1.5 h Sens 78.0% (70.5–84.1)、Spec 92.8%、NPV 96.6%；≥1.5–<2.5 h 89.5%/90.7%/98.3%；≥2.5–<3.5 h 92.9%/90.4%/98.9%；≥4.5–<6 h 95.2%/89.5%/99.3%；≥6–<9 h 92.8%/88.1%/98.2%。女 36.99：≥2.5–<3.5 h Sens 95.8%、Spec 91.9%；男 57.27：0–<1.5 h Sens 74.0%。
- **Ortho K252393**：24 家美国 ED，n=2145（999 女/1146 男），≥22 岁，冻存 Li-hep 血浆，心脏科医生小组按第四版 UDMI 裁定，MI 患病率 8.16%（女 6.21%、男 9.86%），窗 0–2/≥2–4/≥4–6/≥6–11 h。逐窗表格文件有但本文未逐格转录；按受试者的假阴性率：女 9 ng/L 3.2% (2/62) vs 总体 11 ng/L 6.5% (4/62)；男 12 ng/L 8.0% (9/113) vs 11 ng/L 7.1% (8/113)。PPV 95% CI 下限：女 9 ng/L 28.54%（≥4–6 h）、男 12 ng/L 29.42%（≥2–4 h）。
- **Abbott i-STAT K253051**：28 个美国站点，n=3586（2296 女/1290 男），Li-hep 静脉全血与离心血浆平行测定，心脏科+急诊科医生小组按第四版 UDMI 盲法裁定，MI 患病率女 6.8%、男 11.6%。全血、性别特异 cutoff（女 13/男 28）：女 0–<1 h Sens 92.19%（单侧 97.5% CI 下限 86.22）、Spec 83.12%；>1–3 h 96.64%/82.82%；>3–6 h 97.14%/77.90%；男 0–<1 h 79.23%/84.30%；>1–3 h 90.68%/84.09%；>3–6 h 94.20%/82.16%。血浆结果几乎一致（女 0–<1 h 92.19%/82.43%）。女性用 13 ng/L 比 21 ng/L 假阴性率低 ≤3.82%；男性用 28 比 21 高 ≤2.00%。
- **PATHFAST K231974**：临床性能"See K100130"（本文件未载明）；申报论证：降低 LoQ 不改变 99th 与灵敏度/特异性。
- **VIDAS K063243**（原型）：美国社区医院 ED n=302（36 MI）：入院 0–6 h Sens 69.44% (52.80–82.20)、Spec 96.24%；首采后 4–12 h Sens 97.22%、Spec 96.21%。法国心脏 ICU n=243（116 MI）：0–6 h 78.45%/90.55%；4–12 h（n=153）98.70%/92.11%。与商品化 TnI 方法比对 n=534，Passing-Bablok slope 0.42 (0.38–0.44)，r 0.97；一致率 87.2–95.3%。
- **方法学比对 vs predicate**：hs 产品普遍写"Not applicable"，以临床研究替代（Beckman、Siemens、Abbott、Ortho）；i-STAT Alinity 对 i-STAT 1 做 EP09c 比对（数值本文未转录）。

### 9. 厂家间差异与要点
1. **99th 值差异极大**（Ortho 11 / i-STAT 21 / Roche 19 / Beckman 18 / Abbott 28 / PATHFAST 29 / Siemens 47 ng/L），与抗体表位、标准化（Abbott/i-STAT/VIDAS 溯源 NIST SRM 2921；Beckman/Siemens 溯源"商品化 troponin 方法"；Ortho 溯源内部标准；Roche 溯源前一代 TnT）以及参考人群筛查严格度相关——Abbott、Ortho、i-STAT 用 BNP/NT-proBNP + HbA1c + eGFR 生物标志物筛查，Roche/Beckman/Siemens 用问卷式健康定义；Siemens 允许剔除 2 名 outlier 女性。
2. **统计方法**：多数用非参数法；Abbott ARCHITECT 用 EP28-A3c robust 法；Roche 2017 文件未写明。
3. **样本类型策略**：Roche/Ortho 只留 Li-heparin 血浆，Abbott ARCHITECT 只留 K2-EDTA，Beckman/Siemens 血清+Li-hep 并分别给出 99th（血清/血浆值略不同），POC 产品覆盖全血。
4. **LoQ 定义**：新产品统一 20% CV（EP17-A2 精密度曲线）；同时 FDA 关注 99th 处 CV ≤10%（Roche 给 11 ng/L=10% CV；Beckman 给 10% CV 处 4.6 pg/mL；Ortho 给 99th 处 4.5–4.7%）。Abbott ARCHITECT、Siemens、i-STAT 文件未直接给出 10% CV 浓度（Siemens 仅声明满足 IFCC 标准）。
5. **临床试验共性**：ED 疑似 ACS 前瞻多中心、独立盲法裁定（UDMI）、按采样时间窗报告 Sens/Spec/PPV/NPV、总体与性别 cutoff 双套分析、说明书写入 PPV 下限与假阴性率差异。冻存样本需额外提供稳定性证据（Abbott、Ortho、Roche 2021、i-STAT 未提冻存）。
6. **特殊干扰**：Beckman（ALP 标记）受内源 ALP >400 U/L 与 asfotase alfa 影响；Roche 受 biotin（2021 版加 scavenger）；Abbott 受 fibrinogen；Siemens 写入自身抗体致负偏倚并建议首采阴性者再采两次；i-STAT 受 HCT >55% 影响并声明 <1% 卡式假高/假低。

### 10. 来源
- K253051：https://www.accessdata.fda.gov/cdrh_docs/reviews/K253051.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253051.pdf
- K252393：https://www.accessdata.fda.gov/cdrh_docs/reviews/K252393.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K252393.pdf
- K231974：https://www.accessdata.fda.gov/cdrh_docs/reviews/K231974.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K231974.pdf
- K191595：https://www.accessdata.fda.gov/cdrh_docs/reviews/K191595.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K191595.pdf
- K172783：https://www.accessdata.fda.gov/cdrh_docs/reviews/K172783.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K172783.pdf
- K171274：https://www.accessdata.fda.gov/cdrh_docs/reviews/K171274.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K171274.pdf
- K162895：https://www.accessdata.fda.gov/cdrh_docs/reviews/K162895.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162895.pdf
- K201441：https://www.accessdata.fda.gov/cdrh_docs/reviews/K201441.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K201441.pdf
- K063243：https://www.accessdata.fda.gov/cdrh_docs/reviews/K063243.pdf

---

## 利钠肽（NT-proBNP / BNP；B-type natriuretic peptide test system）

### 1. 法规定位
- Product code：**NBC**（B-type natriuretic peptide test system）
- 21 CFR：**862.1117**；Class **II（special controls）**；Panel：Clinical Chemistry (75)
- Special control：*Class II Special Controls Guidance Document for B-Type Natriuretic Peptide Premarket Notifications; Final Guidance for Industry and FDA Reviewers*（2000-11-30）——K032646、K072437、K073091、K232164、K223637 均列入引用。
- openFDA 该代码累计清关数：**total = 45**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K252169 | 2026-03-24 | Beckman Coulter / Access BNP II（DxI 9000） | 顺磁微粒化学发光 | K2-EDTA 血浆（唯一） | ED 疑似新发/失代偿/加重 HF 的辅助诊断；ACS 风险分层；HF 风险分层 | 决策摘要 + summary |
| K253539 | 2026-02-18 | Axis-Shield（Abbott）/ Alere NT-proBNP for Alinity i | CMIA | 血清、血浆 | ED 确诊新发/加重 HF 患者的 **12 个月预后风险分层**（全因死亡、心脏死亡、心脏相关住院）；ACS 风险分层 | 决策摘要 + summary |
| K232164 | 2024-04-12 | Beckman Coulter / Access NT-proBNP（DxI 9000） | 顺磁微粒化学发光（兔抗人 NT-proBNP） | 血清、K2-EDTA、Li-heparin 血浆 | ED 急性 HF 辅助诊断；HF 严重度评估；HF 风险分层；ACS 风险分层 | 决策摘要 + summary |
| K223637 | 2023-07-21 | Roche / Elecsys proBNP II & proBNP II STAT（cobas e 601） | ECLIA（18 / 9 min） | 血清、血浆 | HF 辅助诊断；**新增 ED ADHF 年龄分层 cutoff**；ACS/HF 风险分层；稳定 CAD 心血管事件风险 | 决策摘要 + summary |
| K220265 | 2023-09-24 | Siemens / ADVIA Centaur NT-proBNPII (PBNPII) | 化学发光 | 血清、EDTA / Li-heparin 血浆 | ED 与门诊（OP）新发或加重 HF 辅助诊断 | 决策摘要 + summary |
| K201312 | 2021-10-04 | Ortho / VITROS NT-proBNP II（VITROS 3600） | 化学发光免疫（HRP） | 血清、K2-EDTA、Li-heparin 血浆 | HF 辅助诊断；HF 严重度评估 | 决策摘要 + summary |
| K192380 | 2020-08-24 | Fujirebio (Tosoh) / ST AIA-PACK BNP | 荧光酶免疫（AIA） | K2-EDTA 血浆 | ED 新发/失代偿/加重 HF 辅助诊断 | 决策摘要 + summary |
| K072437 | 2008-02-05 | Roche / Elecsys proBNP II（多抗→单抗） | ECLIA | 血清、Li/NH4-heparin、K2/K3-EDTA 血浆 | CHF 辅助诊断；ACS/CHF 风险分层；稳定 CAD 事件/死亡风险 | 仅决策摘要 |
| K073091 | 2008-02-29 | bioMérieux / VIDAS NT-proBNP | ELFA | 血清、Li-heparin 血浆（EDTA 不推荐） | 疑似 CHF 辅助诊断 | 仅决策摘要 |
| K053597 | 2006-07-21 | i-STAT Corp（Abbott）/ i-STAT BNP（POC） | 卡式 ELISA 电化学 | EDTA 全血或血浆 | CHF 辅助诊断及严重度评估 | 仅决策摘要 |
| K032646 | 2003-11-12 | Roche / Elecsys proBNP（"原型"，多克隆羊抗体） | ECLIA | 血清、血浆 | CHF 辅助诊断 + 新增 ACS/CHF 风险分层 | 仅决策摘要 |

### 3. 预期用途与声明类型
- 核心声明：**aid in the diagnosis of (congestive) heart failure**；2020 年后新清关产品把人群细化为"ED 就诊、临床怀疑新发、急性失代偿或加重 HF"（K192380、K252169、K232164、K223637、K201312），Siemens K220265 与 Ortho K201312 另加**门诊（OP）**人群。
- 附加声明：① **HF 严重度评估**（i-STAT K053597、Ortho K201312、Beckman NT-proBNP K232164）；② **ACS 风险分层 / HF 风险分层**（Roche 2003 起，Beckman 2024/2026，Alere 2026）——文件显示均以**同行评议文献**支持（Roche/Beckman NT-proBNP：GUSTO IV 亚组 James 2003、Jernberg JACC 2002、Fisher Heart 2003；Beckman BNP：deLemos NEJM 2001、Doust BMJ 2005、Vrtovec 2003、Harrison 2002、Logeart 2004、2014 AHA/ACC NSTE-ACS 指南）；③ **12 个月预后**（Alere K253539，前瞻随访自有队列）；④ 稳定 CAD 患者心血管事件与死亡风险（Roche K072437 起）。
- 全部处方用；i-STAT 为 POC（全血）。
- 排除 / rule-out 逻辑：ED 用 NT-proBNP 产品设置年龄无关阴性 cutoff（300 pg/mL）与年龄分层阳性 cutoff，中间为灰区（indeterminate）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**NT-proBNP**
| 产品 | cutoff（文件原值） | 建立方式 |
|---|---|---|
| Roche Elecsys proBNP（K032646, 2003）/ proBNP II（K072437, 2008） | 125 pg/mL（<75 岁）；450 pg/mL（≥75 岁） | K022516 建立（本文件"previously established"）；未说明为 ROC 还是参考区间；K072437 用 NYHA I–IV 及参考组做 cutoff 一致性 |
| VIDAS NT-proBNP（K073091） | 125 / 450（同上） | "based on the Roche Elecsys proBNP assay which the VIDAS assay is traceable to"；年龄匹配 ROC 总 AUC 0.965 |
| Roche proBNP II（K223637, 2023） | ED-ADHF：<300 阴性（全年龄）；阳性 >450（<50 岁）/ >900（50–75）/ >1800（>75）；其间灰区。保留 125/450（predicate 值） | 申办方以美国 17 站点 ED 队列（ICON-Reloaded 队列，文件第 985 行提及）验证既定 cutoff；文件未称由 ROC 重新推导 |
| Beckman Access NT-proBNP（K232164） | <300 ng/L 不太可能；≥450（22–<50 岁）/ ≥900（50–≤75）/ ≥1800（>75）可能；其间灰区 | 与已清关 predicate 的 cutoff 对齐，用自有 2384 例 ED 队列验证；健康人 675 例参考区间另给（EP28-A3） |
| Siemens Centaur PBNPII（K220265） | ED：<300 阴性；灰区 ≥300–≤450/≤900/≤1800；阳性 >450/>900/>1800。OP：<125 阴性、≥125 阳性 | 两个前瞻队列验证；OP 队列 ROC AUC 0.839 (0.804–0.868) |
| Ortho VITROS NT-proBNP II（K201312） | ED：rule-out 300；rule-in 450/900/1800。OP：<125 阴性；≥125"考虑 HF 及其他升高原因" | 同上；predicate 值 125/450 |
| Alere NT-proBNP Alinity（K253539） | 预后用**人群四分位**：≤1214.7 / >1214.7–3261.3 / >3261.3–7031.3 / >7031.3 pg/mL | 861 例随访队列的四分位（非固定 cutoff） |

**BNP**
| 产品 | cutoff | 建立方式 |
|---|---|---|
| i-STAT BNP（K053597, 2006） | >100 pg/mL 异常 | 与 Abbott ARCHITECT/AxSYM 对齐（一次性赋值校准，slope 0.97）；AxSYM 队列 ROC：100 pg/mL 处 Sens 74.2%、Spec 91.5%，AUC 0.90 (0.86–0.92)；年龄匹配后 AUC 0.87 (0.85–0.90) |
| Tosoh ST AIA-PACK BNP（K192380, 2020） | 100 pg/mL | "As recommended in the AHA guidelines for acute heart failure（Weintraub 2010 Circulation 122:1975）"；430 例健康人参考分布（95th 143.2 pg/mL；91.4% <100） |
| Beckman Access BNP II（K252169, 2026） | ≤100 正常；>100 异常提示 HF | "Same as described in K033383"（Triage BNP for Beckman）；用 1323 例 ED 队列验证 |

**参考人群（健康人）定义与数值（文件所载）**
- Beckman NT-proBNP K232164：675 例（369 女/306 男），排除心血管病、未控高血压（≥140/85）、服心血管药（降压药除外）、BMI ≥30、糖尿病、CKD、其他严重慢性病、急性感染；Li-hep 女性 95th：<50 岁 152、50–75 岁 282、>75 岁 443 ng/L；男性 84/188/243；<125 ng/L 比例女 67.2%、男 86.6%。
- Siemens K220265：723 例自报无 HF（362 女/361 男），EP28-A3c 非参数；95th：男 <50 124、50–75 322、>75 154；女 133/192/178；总体 163 pg/mL。
- Ortho K201312：血清 385 女/374 男，排除吸烟、心脏病、高血压、肾病、糖尿病、5 年内癌症/卒中/肺病、高脂、甲状腺病、妊娠，加 troponin >99th（VITROS TnI ES 0.034 ng/mL）、HbA1c >6.5%、eGFR ≤60；95% 参考区间上限：女 22–<50 95.3、50–<75 221、≥75 296；男 125/299/326；总体 217 pg/mL。
- i-STAT BNP K053597（AxSYM 数据转移）：890 例非 HF（含非透析肾病、糖尿病、高血压、COPD 门诊患者，BNP 与健康人无差异）：全体 95th 135 pg/mL，91.5% <100；≥75 岁 95th 254。

### 5. 生物学 / 生理学依据
- **申报文件所载**：利钠肽因心肌壁应力升高而释放；慢性 HF 患者即使临床稳定也可持续升高（K252169 引 AHA/ACC/HFSA Stage C 定义）；BNP 与 BMI 负相关（K252169：31 例假阴性中 74.2% BMI ≥30；K220265：44 例假阴性中 93% BMI ≥30；K192380：BMI ≥37 Sens 仅 76%）；随肾功能恶化升高（K220265：eGFR <60 者假阳性 45.9% vs 17.6%；K192380 eGFR<60 Spec 59%）；随年龄升高（K192380 ≥75 岁 Spec 53%）；有 HF 史者假阳性高（K220265 46.1% vs 16.1%）；透析、nesiritide 输注患者排除（Natrecor 不与 NT-proBNP 检测交叉，K032646）。Roche K223637 说明书：eGFR <60 者假阳性率更高，<50 及 50–75 岁假阴性率更高，"use caution"。
- **背景（非申报文件）**：proBNP(1-108) 裂解为有活性的 BNP(1-32) 与无活性的 NT-proBNP(1-76)；NT-proBNP 半衰期更长、更依赖肾清除，故年龄/肾功能依赖性更强；ED 年龄分层 450/900/1800 与 rule-out 300 源自 ICON 合并分析（Januzzi 2006），BNP 100 pg/mL 源自 Breathing Not Properly（Maisel 2002）；上述文献名称在申报文件中未直接出现（Roche 文件仅提 ICON-Reloaded 队列）。

### 6. 样本类型与样本要求
| 产品 | 样本 | 基质研究 / 稳定性（文件所载） |
|---|---|---|
| Roche proBNP（K032646） | 血清、肝素血浆 | EDTA 血浆结果**约低 10%**；接受标准回收 90–110% 或 slope 0.9–1.1 + r>0.95 |
| Roche proBNP II（K072437） | 血清（含分离胶）、Li-/NH4-heparin、K2-/K3-EDTA | Passing-Bablok vs 血清：Li-hep slope 1.000 (n=30)、NH4-hep 0.9838、K2-EDTA 0.9986 (n=58)、K3-EDTA 0.994；回收 97.1–100.9% |
| VIDAS NT-proBNP（K073091） | 血清（含胶）、Li-heparin | 63 例：Li-hep slope 0.98、胶管 1.00；**EDTA 不符合（nonconformity），说明书不推荐** |
| i-STAT BNP（K053597） | EDTA 全血/血浆 | 仅 EDTA；毛细血管/指尖血未验证 |
| Tosoh BNP（K192380）/ Beckman BNP II（K252169） | 仅 K2-EDTA 血浆 | Beckman：分离后 EDTA 血浆室温 2 h（文件第 509 行）；无基质研究（单一基质） |
| Beckman NT-proBNP（K232164） | 血清、K2-EDTA、Li-hep | 68 例供者配对，Passing-Bablok，三重复；精密度/线性/临床均三基质完成 |
| Siemens PBNPII（K220265） | 血清（无胶/SST/RST）、K2-EDTA、Li-hep | 50 例配对 vs 血清：K2-EDTA slope 1.00 截距 1（50–30,134 pg/mL，r 1.000）；Li-hep 1.00/4；SST 1.00/−2；RST 1.00/4 |
| Ortho NT-proBNP II（K201312） | 血清、K2-EDTA、Li-hep | 160 例配对 Passing-Bablok（数值本文未转录）；稳定性条款在文件第 525–527 行 |
| Roche proBNP II（K223637） | 同 K210546 | 临床样本冻存并提供稳定性证据 |

### 7. 分析性能验证所依据的标准

| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K032646 (2003) | NCCLS EP-5A；BNP special controls guidance (2000) | 精密度（6 次/天 ×10 天）；分析灵敏度 5 pg/mL，功能灵敏度（20% CV）<50 pg/mL |
| K072437 (2008) | BNP guidance (2000)；EP6-A；EP05-A2；EP17-A | 线性（%回收标度，≤10–15% 偏离）；精密度（含 450 cutoff 两侧样本 CV 2.3%/3.4%）；LoB/LoD（claim 5）/LoQ（20% CV，claim 50 pg/mL）；EP9 型方法比对 n=1551 |
| K073091 (2008) | EP5-A2；EP17-A；EP9-A2；EP6-A；BNP guidance；leftover specimen 知情同意指南 (2006) | 3 点 ×2 批精密度（n=240/样本）；LoB 3.4 / LoD 6.7（claim <20）/ LoQ 21.9 pg/mL；方法比对 n=713（slope 0.905, r 0.989；PPA 99.48%、NPA 97.12%）；线性/回收 |
| K053597 (2006) | EP7-A；EP9-A2；**C28-A2**（参考区间转移） | 干扰；vs ARCHITECT 比对（Deming，slope 0.97）；参考范围可转移性论证 |
| K192380 (2020) | EP28-A3c；EP06-A；EP07 3rd；EP37；EP17-A；EP05-A3 | 参考区间 430 例；线性 14 水平；干扰/交叉；LoB 0.9/LoD 1.9/LoQ 3.5 pg/mL；精密度 |
| K201312 (2021) | EP25-A；EP05-A3；EP17-A2；EP06-A；EP07 3rd；EP37；EP28-A3c | 试剂稳定性；精密度；LoB 0.034/LoD 0.49/LoQ claim 20.0 pg/mL（观测 0.46）；线性 20–30,000；干扰；参考区间 |
| K220265 (2023) | EP05-A3；EP06 2nd；EP07 3rd；EP17-A2；EP28-A3c；EP37；（summary 另引 EP09c-ed3 基质等效） | LoB 13 / LoD 20 / LoQ 35 pg/mL（测量范围 35–35,000）；参考区间；基质比对 |
| K223637 (2023) | EP05-A3；EP06 2nd；EP17-A2；BNP guidance | 精密度 8 血清 +2 质控 21 天 3 批（例：68.3 pg/mL 中间精密度 4.8%）；LoB 1.48 / LoD 2.57 pg/mL（e 601），LoQ claim 36 pg/mL；e 411 vs e 601 比对 n=157 PB slope 0.988、截距 1.03、r 0.999 |
| K232164 (2024) | EP-09c；EP05-A3；EP17-A2；EP07 3rd；EP25-A；EP06 2nd；BNP guidance | LoB 1.1；LoD 血清 4.8 / K2-EDTA 3.2 / Li-hep 4.3 ng/L；LoQ 4.8 ng/L（≤20% CV）；无 hook 至 400,000；范围 10–35,000 ng/L |
| K252169 (2026) | EP09c；EP05-A3；EP17-A2；EP07 3rd；EP06 2nd；**EP24-A2**（ROC）；EP28-A3c | LoB 0.3 / LoD 1 / LoQ 2 pg/mL；vs Access BNP 144 例加权 Deming；ROC 评价 |
| K253539 (2026) | EP32-R（溯源）；EP09c 3rd | 器件不变（K241176），本次仅新增预后声明 |

### 8. 临床验证设计与结果

**NT-proBNP，ED 急性 HF/ADHF 队列（设计高度一致：≥22 岁、ED 呼吸困难/怀疑 HF、排除透析或 eGFR<15 / CKD 4–5、排除非 HF 明确原因的呼吸困难与胸外伤；独立中心裁定小组（心脏科/HF 专科或急诊科医生）；以年龄分层阳性 + 300 阴性 cutoff 报告患病率、验后概率、LR）**
- **Roche K223637**：17 家美国站点，n=1485（741 女/744 男，22–97 岁），STAT 版全体：<50 岁患病率 7.5% (36/478)，阳性验后概率 53.4% (31/58)，阴性验后"无 ADHF"98.7% (392/397)，LR+ 14.10、LR− 0.16；50–75 岁患病率 21.6%，阳性 56.2% (146/260)，灰区 17.6%，阴性 97.5% (429/440)，LR+ 4.64、LR− 0.09；>75 岁患病率 34.5%，阳性 58.1% (43/74)，灰区 21.2%，阴性 100% (25/25)，LR+ 2.63。
- **Siemens K220265**：30 站点，n=3128（1148 急性 HF/1980 非），全体：<50 岁验前 25.3%，>450 验后 63.0% (209/332)，LR+ 5.01，<300 者 HF 仅 2.4% (13/549)，LR− 0.07；50–75 岁验前 37.0%，>900 验后 68.5%，<300 者 4.8%，LR− 0.09；>75 岁验前 47.6%，>1800 验后 68.3%，<300 者 5.4%，LR− 0.06。**门诊队列**：28 站点 n=1033（185 新发 HF），cutoff 125：男 ≤75 Sens 85.5% (53/62)、Spec 68.9%、NPV 95.2%；女 ≤75 79.2%/70.8%/95.7%；男 >75 90.5%/37.2%；女 >75 100%/20.8%；AUC 0.839。
- **Ortho K201312**：ED 20 站点 n=2200（1095 HF/1105 非，22–106 岁，冻存）；交叉表：阳性 949 HF/196 非，灰区 133/310，阴性 13/599；22–<50 岁阳性验后 84.7%，阴性验后非 HF 96.5%，LR+ 6.84、LR− 0.05；50–<75 岁 80.4%/98.1%，LR+ 4.81、LR− 0.02。**门诊** n=777（10 站点），cutoff 125：Sens 91.7% (44/48)、Spec 67.2% (490/729)、NPV 99.2%、PPV 15.6%。
- **Beckman K232164**：17 采集点 +3 检测点，n=2384（1059 HF/1325 非；<50 岁 879 例 HF 35.3%，50–75 岁 915 例 48.7%，>75 岁 590 例 51.4%），Li-hep：阳性 899 HF/326 非、灰区 118/253、阴性 42/746（全体）。
- **VIDAS K073091**（2008，病例–对照）：407 例确诊 CHF（3 站点）vs 411 例无 CHF/心血管病史参考者，125/450：欧洲站 1 Sens 94.63%/Spec 97.39%；欧洲站 2 94.96%/96.69%；美国站 84.87%/81.31%；男性合计 92.92%/95.67%；女性 90.42%/89.44%；AUC 0.965。
- **Roche K032646/K072437**（K022516 数据）：按年龄/性别 Sens/Spec，男 <75 岁 Sens 89.0% (85.95–91.58)、Spec 90.0%；女 <75 岁 Sens 90.6%、Spec 76.7%；≥75 岁男 86.5%/88.9%、女 81.8%/87.9%；NPV 96.8–100%。K072437 vs proBNP 比对 n=1551 PB y=0.979x−0.314、τ 0.932、r 0.996；cutoff 一致率参考组 97.6–99.0%、NYHA I–IV 99–100%。
- **Alere K253539**（预后）：K241176 队列 880 例裁定中 861 例随访 6±1、12±2 个月（495 男/366 女，平均 62.0 岁）；12 个月复合终点绝对风险 Q1 44.03% (36.86–50.39) → Q4 63.08% (55.99–69.02)；全因死亡 10.21% → 25.50%；心脏死亡 2.03% → 15.78%；心脏相关住院 39.69% → 53.41%；KM 与 Cox（单/多变量）结果文件有，HR 数值本文未转录。

**BNP，ED 队列（单一 cutoff 100 pg/mL）**
- **Beckman K252169**：18 采集点，n=1323（572 女/751 男；449 HF/874 非；冻存），Sens 93.1% (418/449, 90.4–95.1)，Spec 65.6% (573/874, 62.3–68.6)，PPV 58.1%，NPV 94.9%；<50 岁 Sens 93.7%/Spec 82.0%；50–75 岁 91.8%/66.4%；>75 岁 96.1%/38.1%（NPV 93.0%）。假阴性 6.9%，假阳性 34.4%（57.1% 有 HF 史，36.9% eGFR<60，28.6% >75 岁）。NYHA I/II/III/IV 中位数 900/661/712/808 pg/mL。
- **Tosoh K192380**：8 家 ED，825 份/724 例分析（329 HF）：Sens 88.4% (84.5–91.5)、Spec 70.6% (66.0–74.9)、PPV 71.5%、NPV 88.0%、LR+ 3.012、LR− 0.164；一致率 78.7%。
- **i-STAT K053597**（AxSYM 数据转移，非自有队列）：100 pg/mL 男性总 Sens 71.0% (328/462)、Spec 94.8%；女性 80.5%/88.4%；≥75 岁男 86.1%/89.5%；HF 组 693 例（NYHA I–IV 中位数 133/266/335/1531 pg/mL）。参考范围转移依据 C28-A2 + slope 0.97。

### 9. 厂家间差异与要点
1. **NT-proBNP 三个时代的 cutoff**：2003–2008（125/450 两档，Roche 起源，VIDAS 通过溯源直接借用）→ 2021 起 ED 用 300 + 450/900/1800 + 灰区（Ortho 首先在本清单中写入，Roche 2023 才把自家产品说明书更新为同一套），门诊仍用 125。各厂无一在文件中重新用 ROC 推导 cutoff，均为"验证既定 cutoff"。
2. **BNP 阵营**统一 100 pg/mL，但依据不同：i-STAT 用校准对齐 + AxSYM ROC，Tosoh 引 AHA 2010 科学声明，Beckman 引 predicate K033383。
3. **风险分层声明**几乎全靠文献（FDA 接受），Alere 2026 是清单中唯一用自有前瞻随访数据建立 12 个月预后声明的产品，且用四分位而非固定 cutoff。
4. **性能报告口径**：新一代 ED 研究不再报告总体 Sens/Spec/AUC，而是按年龄层报告验前/验后概率与 LR（Roche、Siemens、Ortho、Beckman NT-proBNP）；BNP 产品仍报告 Sens/Spec/PPV/NPV。说明书必须写入 BMI、eGFR、HF 史、年龄对假阴/假阳的影响。
5. **参考区间 vs cutoff 分离**：健康人 95th/97.5th（Beckman 女 >75 岁 95th 443 ng/L、Siemens 总体 163 pg/mL、Ortho 总体 217 pg/mL）远低于 ED 阳性 cutoff，说明书同时列出。
6. **样本**：BNP 类（i-STAT、Tosoh、Beckman）只允许 EDTA；NT-proBNP 允许血清/肝素/EDTA（早期 Roche 提示 EDTA 低 10%，VIDAS 拒绝 EDTA）。
7. **LoQ 口径**：从"功能灵敏度 <50 pg/mL"（Roche 2003/2008 claim 50）到 EP17-A2 精密度曲线（Ortho 20、Roche 2023 36、Siemens 35、Beckman 4.8 ng/L）。

### 10. 来源
- K252169：https://www.accessdata.fda.gov/cdrh_docs/reviews/K252169.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K252169.pdf
- K253539：https://www.accessdata.fda.gov/cdrh_docs/reviews/K253539.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K253539.pdf
- K232164：https://www.accessdata.fda.gov/cdrh_docs/reviews/K232164.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K232164.pdf
- K223637：https://www.accessdata.fda.gov/cdrh_docs/reviews/K223637.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K223637.pdf
- K220265：https://www.accessdata.fda.gov/cdrh_docs/reviews/K220265.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K220265.pdf
- K201312：https://www.accessdata.fda.gov/cdrh_docs/reviews/K201312.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K201312.pdf
- K192380：https://www.accessdata.fda.gov/cdrh_docs/reviews/K192380.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K192380.pdf
- K072437：https://www.accessdata.fda.gov/cdrh_docs/reviews/K072437.pdf
- K073091：https://www.accessdata.fda.gov/cdrh_docs/reviews/K073091.pdf
- K053597：https://www.accessdata.fda.gov/cdrh_docs/reviews/K053597.pdf
- K032646：https://www.accessdata.fda.gov/cdrh_docs/reviews/K032646.pdf

---

## 高敏 / 心脏 C 反应蛋白（hsCRP / cardiac CRP, cCRP）

### 1. 法规定位
- Product code：**NQD**（"Cardiac C-Reactive Protein, Antigen, Antiserum, and Control"）
- 21 CFR：**866.5270**（C-reactive protein immunological test system）；Class **II**；Panel：Immunology (82 / IM)（Roche K042485 标注 Chemistry 75）
- openFDA 该代码累计清关数：**total = 19**。
- 关键指导文件（文件引用）：*Guidance for Industry and FDA Staff: Review Criteria for Assessment of C-Reactive Protein (CRP), High Sensitivity C-Reactive Protein (hsCRP) and Cardiac C-Reactive Protein (cCRP) Assays*，2005-09-22（K173833、K233242、K260026 引用）；*AHA/CDC Scientific Statement, Circulation 2003;107:499-511*（K033908 直接列为"standard"）。
- 与常规 CRP 的区别（申报文件可见）：NQD 承载"aid in identification/stratification of individuals at risk for future cardiovascular disease"及"independent marker of prognosis for recurrent events in stable coronary disease or ACS"两句声明；常规 CRP 声明仅为"evaluation of infection, tissue injury, inflammatory disorders"。Roche K260026 的最新 cardiac hsCRP III 更明确写"**Not intended for the detection and evaluation of inflammatory disorders**"（predicate K053603 两者兼有）。（背景（非申报文件）：常规 CRP 与 hsCRP 炎症用途产品的 product code 为 DCK/DCN，按用户提示注明，未在本批文件中核实。）

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K260026 | 2026-04-02 | Roche / Tina-quant Cardiac high sensitivity CRP III（cobas c 503） | 乳胶增强免疫比浊（PETIA） | 血清、Li-heparin、K2-EDTA 血浆 | 未来 CVD 风险评估；稳定 CAD/ACS 复发事件预后；**不用于炎症评估** | 决策摘要 + summary |
| K233242 | 2024-01-18 | Siemens / Atellica CH High Sensitivity CRP 2 (hCRP2) | 颗粒增强免疫散射比浊（PENIA） | 血清、Li-/Na-heparin、K2-EDTA 血浆 | 未来 CVD 风险个体识别；ACS 复发预后 | 决策摘要 + summary |
| K212559 | 2022-12-16 | Siemens / CardioPhase hsCRP（即原 N High Sensitivity CRP 更名） | PENIA（BN 系统） | 同 K033908 | 不变（Special 510(k)：溯源 ERM-DA470→ERM-DA474；CRP1 范围 3.1–200→3.1–100 mg/L） | 仅决策摘要（2 页） |
| K173833 | 2018-09-27 | Sentinel CH / CRP Vario（在 Abbott ARCHITECT c8000 上） | 乳胶免疫比浊（兔多抗） | 血清、肝素/EDTA 血浆 | cCRP 风险识别/分层；ACS 复发预后 | 决策摘要 + summary |
| K070626 | 2007-05-04 | Beckman Coulter / Synchron Systems High Sensitivity Cardiac CRP (CRPH) | 近红外颗粒免疫比浊速率法 | 血清、肝素、EDTA 血浆 | 炎症评估 + 心血管风险识别/分层 + 预后（在 K010597 上加心脏声明） | 仅决策摘要 |
| K042485 | 2004-10-29 | Roche / Tina-quant CRP (Latex) HS | PETIA | 血清、血浆 | 炎症 + 冠心病风险 + 复发预后（仅改声明，性能"subject of k003400"） | 仅决策摘要 |
| K041799 | 2004-08-25 | Ortho / VITROS Chemistry Products hsCRP（VITROS 5,1 FS） | 乳胶免疫比浊 | 血清、肝素血浆 | "CRP is used to evaluate the risk of developing CHD. The risk of CHD increases with values of CRP that exceed 3 mg/L." | 仅决策摘要 |
| K033908 | 2004-01-22 | Dade Behring / N High Sensitivity CRP（BN Systems）——"原型" | PENIA | 血清、肝素、EDTA 血浆 | 在 K991385 基础上新增心脏风险声明（文献 + AHA/CDC 声明支持） | 仅决策摘要 |

### 3. 预期用途与声明类型
- 声明类型：**风险分层 / 风险识别**（primary prevention）+ **预后**（secondary：stable CAD 或 ACS 的复发事件独立标志物），均为"aid"性质；无诊断、无监测声明。Roche K042485 文件转录了 AHA/CDC 使用建议：有感染/全身炎症/外伤时不做风险评估；持续 >10 mg/L 应查非心血管原因；最好取间隔 2 周两次结果的平均值；不推荐全民筛查；不能替代传统危险因素；不应据此单独决定 ACS 处理或二级预防；**不用系列 hsCRP 监测治疗**。
- 处方用；均为中央实验室化学/免疫分析仪，本组无 POC 产品（NQD 清单中 Dade Stratus CS K062924/K060369 为急诊台式，未展开）。
- 2004–2007 年的清关多为"给已清关炎症 CRP 加心脏声明"（K033908、K042485、K070626），2018 年后为新试剂或改版，仍**无临床研究**，依靠分析等效 + 文献。

### 4. 阳性 / 阴性判定与 cutoff 逻辑
- 全部产品采用 **AHA/CDC 2003 三分位**（Pearson et al., Circulation 2003;107:499–511；部分另引 Ridker Circulation 2003;107:363–369）：
  - Low <1.0 mg/L；Average 1.0–3.0 mg/L；High >3.0 mg/L（K033908、K070626、K173833、K233242、K260026）。
  - Ortho K041799 加一档：High 3.0–10.0；**Indeterminate >10.0 mg/L**（提示其他炎症/感染来源）。
  - Sentinel K173833 说明书表中 Low 误印为"<0.1"（文件原样）。
- cutoff 的建立：**不是由申办方数据建立**。K042485 决策摘要明文："The cutoff was established previously in the literature not with this device"；"The expected range was established in the literature not with this device"。K041799 论证逻辑：CDC/AHA 分类基于 Dade Behring BN ProSpec 比浊法建立 → 195 例（101 女/94 男）血清与 BN ProSpec 比对等效 → 已发表分类适用于本产品。K070626 同样以 Deming 回归对 Dade CardioPhase 的等效为桥。
- 产品自有参考区间（补充信息，非 cutoff）：Beckman K070626：551 例南加州血库健康非吸烟成人 95% <7.48 mg/L；Dade CardioPhase（在 K070626 比较表引述）2147 例健康人 90th 1.69、95th 2.87 mg/L，并注明"each laboratory should determine its own reference interval"。Roche K260026、Siemens K233242 未提供自有健康人分布，仅列三分位。
- 灰区：无；仅 >10 mg/L 作为"另有炎症来源"提示（Ortho 表格、Roche 文件的 AHA/CDC 文字）。

### 5. 生物学 / 生理学依据
- **申报文件所载**（Roche K042485 器件描述）：CRP 为经典急性期蛋白，肝脏合成，五个相同多肽链组成五聚环，分子量 120,000；与配体复合后经 C1q 激活补体，启动调理与吞噬；主要功能是结合并解毒组织损伤产生的内源性毒性物质；高敏测定用于"apparently healthy persons"冠心病风险预测与复发事件预后。Beckman K070626、Siemens K233242 等直接引 AHA/CDC 2003 声明。
- **背景（非申报文件）**：hsCRP 是 IL-6 驱动的下游炎症标志物，反映动脉粥样硬化的低度慢性炎症；"hs"指测量下限 ≤0.3 mg/L 以区分 1/3 mg/L 三分位，而非分析原理不同。

### 6. 样本类型与样本要求
| 产品 | 样本 | 基质研究 / 稳定性（文件所载） |
|---|---|---|
| Roche K260026 | 血清、Li-hep、K2-EDTA | 70 例配对 Passing-Bablok：K2-EDTA slope 1.023 截距 −0.029 r 0.999；Li-hep 1.028/−0.027/1.000（0.185–9.23 mg/L）。稳定性：血清/Li-hep 15–25 °C 14 天、2–8 °C 28 天、−20 °C 12 个月；K2-EDTA 室温仅 2 天；可冻融 4 次 |
| Siemens K233242 | 血清、Na-hep、Li-hep、K2-EDTA | 60 例（55 例在范围内）Deming：Na-hep slope 1.070、K2-EDTA 0.970、Li-hep 1.070（0.22–7.40 mg/L）。稳定性：室温 26 h、2–8 °C 188 h、−20 °C 2 年（predicate 为 2–8 °C 8 天、冻 8 个月） |
| Sentinel K173833 | 血清（胶/无胶）、EDTA、Li-hep 胶、Na-hep | 每管型 40 例 vs 无胶血清 Passing-Bablok：血清胶 1.01、EDTA 0.99、Li-hep 胶 0.96、Na-hep 0.97 |
| Beckman K070626 | 血清、Na-/Li-heparin、EDTA | 47 例健康人 Deming：Na-hep 0.984x+0.049；Li-hep 1.017x−0.012；EDTA 0.982x+0.021；脂血/混浊样本不可用 |
| Ortho K041799 | 血清、Li-heparin | 血清 vs SST vs Li-hep 无显著差异（数值文件未载明） |
| Dade K033908 / Siemens K212559 | 血清、肝素、EDTA | "By reference to K991385" |
| Roche K042485 | 血清、血浆 | "subject of k003400"（文件未载明数值） |

### 7. 分析性能验证所依据的标准

| K 号 | 引用标准 | 对应项目 |
|---|---|---|
| K033908 (2004) | AHA/CDC Scientific Statement (Circulation 2003) | 仅新增声明；性能"By reference to K991385" |
| K042485 (2004) | "Not Applicable subject of k003400" | 无新数据 |
| K041799 (2004) | NCCLS EP9-A2；EP5-A；EP6-A；EP7-A | 方法比对 n=110（y=0.954x+0.063，r 0.993）；20 天精密度（0.60 mg/L 5.0% CV；1.84 2.1%；10.64 1.8%）；线性 0.10–15.00 mg/L；干扰（Hb 1000 mg/dL、bilirubin 50、Intralipid 450、RF 1200 IU/mL、18 种药物）；LoD 0.1 mg/L（=20% CV 功能灵敏度） |
| K070626 (2007) | EP5-A2；EP6-A；EP9-A2 | 精密度（3 水平 20 天：0.630 mg/L 总 5.3%）；线性（Y=1.0051X+0.1415；低端 0.9861X+0.0941）；分析灵敏度 0.11（LoB，20 重复）；功能灵敏度 0.18 mg/L（20% CV）；Deming vs Dade n=269（1.048X+0.024, r 0.9899），心脏范围 0.2–10 n=149（1.030X−0.008） |
| K173833 (2018) | FDA 1994 Points to Consider；2005 510(k) 格式指南；**2005 CRP/hsCRP/cCRP review criteria**；EP17-A2；EP07-A2；EP6-A；EP05-A3；**EP21**（总误差） | 精密度（8 水平 20 天，围绕 1 mg/L 决策点：0.98 mg/L 总 2.6%，0.34 mg/L 9.9%）；线性 0.30–10.00（slope 1.006, r 0.9998，≤4%）；LoB 0.02 / LoD 0.04 / LoQ 观测 0.07、claim 0.30 mg/L；干扰（Hb 1984 mg/dL 等，≤±10%）；hook 至 1500 mg/L；vs CardioPhase n=115 PB slope 1.019、截距 0.061、r 0.996 |
| K233242 (2024) | EP05-A3；EP17-A2；EP06-ED2:2020；EP07 3rd；EP37；2005 review criteria guidance | 精密度（1.06 mg/L 内部 1.6%）+ 3 台 ×3 批重现性（≤2.5%）；线性 0.16–9.50；干扰（Hb 1000、bilirubin 40、Intralipid 3000 mg/dL、RF 500 IU/mL）；prozone：1300 mg/L 仍报 >9.50；LoB 0.06 / LoD 0.11 / LoQ 0.16 mg/L；EP09c 加权 Deming vs CardioPhase n=100（slope 0.96、截距 0.03、r 0.999，0.26–9.41） |
| K260026 (2026) | EP05-A3；EP17-A2；EP06-ED2:2020；2005 review criteria guidance；（干扰按 EP07 Ed 03） | 精密度（0.429 mg/L 中间 1.9%）；线性 0.150–10.0（>1 mg/L 偏离 <4%，≤1 mg/L 绝对差 <0.03）；干扰 11 浓度梯度（RF 520 IU/mL、Intralipid 2000 mg/dL、IgG 79.2 g/L）+ 40 种药物（<±8.1%）；hook 至 1000 mg/L；LoB 0.100 / LoD 0.150 / LoQ 0.150 mg/L；vs K053603 n=104 PB slope 1.068、截距 0.0302、r 0.999，医学决策水平预测偏倚 1 mg/L +3.8% (2.1–5.4)、3 mg/L +5.8% (4.9–6.5)；并证明对 K033908 的偏倚相似 |
| K212559 (2022) | Special 510(k)（New 510(k) Paradigm 设计控制） | 溯源材料变更 ERM-DA470→ERM-DA474 |

溯源材料演变：CRM 470/RPPHS（K033908、K041799、K070626）→ ERM-DA472（K173833）→ ERM-DA474/IFCC（K212559、K233242、K260026）。

### 8. 临床验证设计与结果
- **本组 8 份文件全部写明"Clinical sensitivity / specificity: Not applicable"或 N/A**，没有申办方自行开展的临床结局研究。
- 支持心脏声明的证据链：K033908（Dade）提供"substantial peer-reviewed literature using the Dade test"+ AHA/CDC 专家小组声明 → FDA 认为心脏声明恰当；后续产品（Roche K042485、Ortho K041799、Beckman K070626、Sentinel K173833、Siemens K233242、Roche K260026）通过与 Dade/Siemens 比浊法在 0.2–10 mg/L 心脏区间的**方法比对等效**（Deming / Passing-Bablok，n=100–269，slope 0.954–1.068，r 0.990–0.999）"分析桥接"到该文献证据。
- 评价指标因此为：slope/截距/r，以及 Roche K260026 新增的**医学决策水平（1、3 mg/L）预测偏倚及 95% CI**；Sentinel 在 1 mg/L 决策点设计精密度；Siemens 报告 3 台 3 批总重现性。无 Sens/Spec/AUC/PPA/NPA。

### 9. 厂家间差异与要点
1. **声明措辞差异**：Beckman/Roche 2004/Dade 版本同时保留炎症用途 + 心脏用途；Ortho 2004 只写冠心病风险（>3 mg/L）；Sentinel/Siemens 2024 只写心脏用途；Roche 2026 明确排除炎症用途——反映 FDA 2005 review criteria 把 CRP / hsCRP / cCRP 三类分开处理的思路（背景（非申报文件）：该指南按预期用途区分产品代码与所需证据）。
2. **测量范围收窄**：早期 0.2–80 mg/L（Beckman，含 ORDAC 扩展 60–380）→ 0.10–15（Ortho）→ 0.3–10（Sentinel）→ 0.16–9.50（Siemens）→ 0.150–10.0（Roche 2026）：cardiac 专用产品把上限压到 10 mg/L 左右，超出即提示炎症。
3. **检出限**：新产品 LoQ 0.15–0.30 mg/L（EP17-A2，20% CV），旧产品用"功能灵敏度"0.1–0.18 mg/L；所有产品在 1 mg/L 附近 CV <3%（Sentinel 2.6%、Siemens 1.6%、Roche 1.9%）。
4. **标准化**：从 CRM 470 迁移到 ERM-DA474/IFCC 是 2022–2026 年多次申报的主要变更（Siemens Special 510(k)、Roche III 代），FDA 关注对 1/3 mg/L 决策点偏倚的影响（Roche 给出 +3.8%/+5.8%）。
5. **方法学**：Siemens/Dade 沿用散射比浊（PENIA），Roche/Ortho/Sentinel/Beckman 为透射比浊（PETIA/近红外速率），互为 predicate 被 FDA 接受。
6. 无 POC、无全血产品；EDTA 血浆稳定性明显短于血清/肝素（Roche：2 天 vs 14 天）。

### 10. 来源
- K260026：https://www.accessdata.fda.gov/cdrh_docs/reviews/K260026.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf26/K260026.pdf
- K233242：https://www.accessdata.fda.gov/cdrh_docs/reviews/K233242.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K233242.pdf
- K212559：https://www.accessdata.fda.gov/cdrh_docs/reviews/K212559.pdf
- K173833：https://www.accessdata.fda.gov/cdrh_docs/reviews/K173833.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K173833.pdf
- K070626：https://www.accessdata.fda.gov/cdrh_docs/reviews/K070626.pdf
- K042485：https://www.accessdata.fda.gov/cdrh_docs/reviews/K042485.pdf
- K041799：https://www.accessdata.fda.gov/cdrh_docs/reviews/K041799.pdf
- K033908：https://www.accessdata.fda.gov/cdrh_docs/reviews/K033908.pdf

---

### 附：本节未能从文件中获得的信息（汇总）
- Troponin：Roche K162895 参考区间的统计方法与 CI；Abbott ARCHITECT / Siemens / i-STAT 的"10% CV 处浓度"具体数值；PATHFAST K231974 的参考人群与临床性能（决策摘要全部转引 K100130，未下载）；Ortho K252393 逐时间窗 Sens/Spec 表格（文件有，未逐格转录）；Siemens K171274、Ortho K252393 的样本稳定性时限。
- NT-proBNP/BNP：Roche 125/450 cutoff 的原始建立方法（在 K022516，未下载）；Beckman NT-proBNP、Ortho 的逐年龄层 LR/验后概率完整表格（文件有，本文只转录部分）；Alere K253539 的 Cox HR 数值；i-STAT BNP 方法比对回归参数（决策摘要中表格为图片，文本抽取缺失）；Beckman BNP II vs predicate 回归参数（未转录）。
- hsCRP：Roche K042485、Dade K033908、Siemens K212559 的全部分析性能（转引 predicate，未下载 K003400/K991385）；Ortho K041799 的基质比对数值；Roche/Siemens 新产品的自有健康人分布（文件未提供）。
