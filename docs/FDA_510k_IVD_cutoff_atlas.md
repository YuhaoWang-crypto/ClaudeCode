# FDA 510(k) IVD 试剂盒按靶点整理：验证标准、阳性/阴性判定与 cutoff 逻辑图谱

> 整理日期：2026-09-22。数据来源：openFDA `device/510k` 与 `device/classification` 接口（数据更新至 2026-09-14）、FDA 510(k) 决策摘要（Decision Summary，`accessdata.fda.gov/cdrh_docs/reviews/K*.pdf`）与 510(k) Summary（`cdrh_docs/pdfYY/K*.pdf`）、FDA 指南文件。
>
> **严格性标注约定**（全篇通用）：
>
> - 「申报文件」：数字或表述直接取自该 K 号的决策摘要 / 510(k) summary，可按第 10 小节的 URL 复核。
> - 「背景（非申报文件）」：来自临床指南、生理学或文献，用于解释 cutoff 的来源，不是 FDA 文件原文。
> - 「文件未载明」：FDA 公开文件里没有这项信息（常见于早期只有 510(k) summary 的申报，或 “performance previously established in K…” 的修改型申报）。

---

## 0. 这份整理回答什么问题

对每一个靶点（analyte），把市面上不同厂家 510(k) 申报时的以下要素并排放在一起：

| 问题 | 对应小节 |
|---|---|
| 这个靶点在 FDA 归哪一类、哪个产品代码、有多少个已清关产品 | §1 法规定位 |
| 有哪些代表性产品、什么方法学、什么样本 | §2 代表性 510(k) 一览 |
| 申报的预期用途是“辅助诊断”“风险分层”“监测”还是“排除” | §3 预期用途 |
| 阳性/阴性怎么判、cutoff 是多少、cutoff 是怎么定出来的 | §4 判定与 cutoff 逻辑 |
| 为什么这个 cutoff 在生物学/生理学上说得通 | §5 依据 |
| 用什么样本、抗凝剂、稳定性、基质等效 | §6 样本 |
| 分析性能验证依据哪些 CLSI / ISO / 溯源标准 | §7 标准 |
| 临床验证怎么设计、金标准是什么、结果数字 | §8 临床验证 |
| 厂家之间的差异 | §9 |
| 原始文件链接 | §10 |

覆盖的靶点分六组（共 25 个靶点；全文提及 274 个不同的 K 号，含作为 predicate 被引用者，见目录）。选择原则：都是走 510(k)（Class II）的定量或定性 IVD，且 cutoff 逻辑彼此不同，能覆盖下面 §2 里的全部八种原型。每组由一名研究助手独立通读 FDA 文件后按统一模板撰写，关键数字由我逐组抽样回查原文核对（抽查项见 §2.1 表中加粗者）。

---

## 1. 先划边界：哪些靶点不走 510(k)

按靶点整理时，第一步是确认该靶点的预期用途在 FDA 分类里落在 Class II（510(k)）还是 Class III（PMA），同一个分子在不同用途下可以分属不同类别：

| 靶点 | Class II / 510(k)（本文覆盖） | Class III / PMA 或 De Novo（本文不覆盖，只做提示） |
|---|---|---|
| PSA | 产品代码 LTJ：已确诊前列腺癌患者的管理/监测 | 产品代码 MTF/MTG：50 岁以上男性癌症检出（筛查）、free/total PSA 鉴别，均为 PMA |
| AFP | 产品代码 LOJ：睾丸癌/生殖细胞肿瘤监测；NSF：AFP-L3% 肝癌风险 | 产品代码 LOK/LTQ：孕中期母血神经管缺陷筛查，Class III |
| CA 125 | 产品代码 LTK：卵巢癌治疗反应/复发监测 | 卵巢癌筛查声明无 510(k) 途径 |
| ROMA（HE4+CA125） | 产品代码 ONX（21 CFR 866.6050），K103358 及后续同类走 510(k) | 866.6050 条款本身由 OVA1（K081754）经 De Novo 建立 |
| PCT | 产品代码 NTM/PMT/PRI：2008 年 K070310 以 510(k) 清关；2016 年 DEN150009 经 De Novo 建立 21 CFR 866.3215 与 7 项特殊控制，之后新声明走 510(k) | — |
| HBsAg、抗-HCV、抗-HIV | — | 血液筛查用途归 CBER，PMA |
| KRAS / EGFR / BRAF 等肿瘤基因伴随诊断 | — | PMA（伴随诊断）。这也是本仓库 `grn_pipeline` 里讨论的 KRAS G12C 入组 biomarker 所在的法规类别 |
| HbA1c | 产品代码 LCP：监测；2010 年后诊断/筛查声明也走 510(k) | — |
| 肌钙蛋白 | 产品代码 MMI：AMI 辅助诊断 | — |

结论：本文覆盖的全部靶点都在 Class II，依据是 21 CFR 862 / 864 / 866 各条款的“特殊控制（special controls）”。特殊控制通常就是（a）一份 FDA 指南（如 CRP 的 2005 年 Review Criteria、肿瘤标志物的 Tumor Associated Antigen 510(k) 指南）加（b）CLSI 共识标准。这决定了 §7 里看到的标准清单为什么高度趋同。

---

## 2. cutoff 逻辑的八种原型（跨靶点抽象）

把 25 个靶点的 §4 放在一起看，申报文件里 cutoff 的“来源”只有八种写法。这个分类是本文的主线，后面每个靶点都会标注属于哪一种。

| 原型 | cutoff 怎么来 | 申报时要证明什么 | 典型靶点 |
|---|---|---|---|
| **A. 健康参考人群统计百分位** | 在“表观健康”人群中测定分布，取 97.5th（双侧参考区间上限）、95th/99th（单侧决策限）百分位；CLSI EP28-A3c 非参数法，通常 n ≥ 120 每组 | 参考人群的纳入/排除定义（如 hs-troponin 要求排除 NT-proBNP 升高、eGFR 低、HbA1c 高者）、样本量、是否分性别/年龄；决策限处的 CV | hs-cTnI/T（99th URL，性别特异）、TSH、ferritin、CEA、CA 19-9、CA 125、HE4（95th，绝经分层）、anti-CCP 部分厂家（99th）、ANA 固相法 |
| **B. 临床结局驱动的固定决策点** | cutoff 先由指南或 RCT 确定（NT-proBNP 300/450/900/1800 pg/mL，PCT 0.10/0.25/0.5 ng/mL 与 Δ80%，hsCRP 1/3 mg/L，HbA1c 6.5%，Tg 0.2 ng/mL），厂家不再重新推导 | “cutoff 迁移”研究：在预期用途人群中验证该固定 cutoff 下的敏感性/特异性/NPV，或与已建立临床效用的 predicate 做方法比对 + 在 cutoff 处的一致性 | NT-proBNP/BNP、PCT（三类声明）、hsCRP（心血管风险）、HbA1c（诊断）、D-dimer（VTE 排除 500 ng/mL FEU）、Tg |
| **C. ROC 最优点** | 用确诊病例 vs 疾病对照 + 健康人做 ROC，取 Youden 或预设特异度（如 ≥95–98%）对应的浓度 | 训练集与验证集分开；疾病对照谱系（spectrum）要覆盖交叉反应疾病；EP24 | anti-CCP 部分厂家（Siemens IMMULITE、Roche、Euroimmun）、粪便钙卫蛋白原型（均值+2SD 再 ROC）、HemosIL D-dimer HS |
| **D. 阴性人群分布 + 灰区** | 血清学 index/ratio：以阴性人群均值 + k·SD 或 predicate 对齐定 1.0，设 equivocal 区（如 0.9–1.1） | 与参考方法（Western blot、TPPA、FTA-ABS、CDC 血清盘）的 PPA/NPA；equivocal 的比例和处理；交叉反应盘 | HSV-2 IgG、Lyme EIA、梅毒 treponemal、CMV IgG/IgM、ENA 单项 |
| **E. 算法 / 多参数** | 多个测量值合成一个分数或串联判读 | 每个成分的性能 + 算法整体的临床性能；分层（绝经前/后、疾病分期） | ROMA（HE4+CA125，0–10 分制 1.31/2.77 或 1.14/2.99）、Lyme 两步法（STTT / MTTT）、梅毒反向算法 |
| **F. 无临床 cutoff、量值溯源型** | 不声明诊断 cutoff，只给参考区间；申报重心是校准溯源（NIST SRM、WHO IS、ERM-DA471、NGSP/IFCC）和与参考方法的偏倚 | ID-LC-MS/MS 或参考实验室比对的 Passing-Bablok 斜率/截距、总误差 | 25-OH 维生素 D、cystatin C、ferritin、TSH（参考区间型） |
| **G. 监测型 serial change** | 参考区间上限 + “有临床意义的变化”定义（由分析 CV 推导，或 RCV = 1.96·√2·√(CVa²+CVi²)），预期用途限定为已确诊患者的监测 | 在监测人群中，标志物变化方向与疾病状态变化（进展/缓解/稳定，由影像或临床判定）的一致性 %；不做筛查性能 | CEA、CA 125、CA 19-9、PSA（LTJ）、AFP（LOJ）、HE4、Tg |
| **H. 与 predicate 对齐（派生型）** | 新平台直接沿用 predicate 的 cutoff，只做方法比对 + cutoff 处分类一致性 | 方法比对回归、PPA/NPA；FDA 2005 CRP 指南明文允许“文献 + 桥接研究” | 几乎每个靶点的后续厂家：hsCRP 全部产品、Abbott/Siemens anti-CCP、PCT 借用四项声明的各平台 |

一个靶点常同时用两种原型：hs-troponin 的 cutoff 是 A 型（99th URL），但多数申报同时报告 B 型的 0/1h 或 0/2h 算法性能；CA 125 是 A 型上限 + G 型 serial change；anti-CCP 三种路径并存（A、C、H）。

### 2.1 靶点 × 原型 × 代表性 cutoff 速览（数字均取自分节 §4，加粗者为本人回查原文核对过）

| 靶点 | 原型 | 代表性 cutoff（申报文件） | 备注 |
|---|---|---|---|
| hs-cTnI / hs-cTnT | A (+B) | 总体 99th URL：Ortho 11、Roche TnT 19、Beckman 18、i-STAT 21、Abbott 28、PATHFAST 29、**Siemens 47.34** ng/L；均另给性别特异值 | 参考人群 n≈900–2000，部分加 NT-proBNP/HbA1c/eGFR 筛查；LoQ 统一按 20% CV |
| NT-proBNP | B | ED：**<300 rule-out；450/900/1800** 按 <50/50–75/>75 岁 rule-in，中间为灰区；门诊 125 pg/mL | 2021 Ortho 起写入灰区，Roche 2023 跟进 |
| BNP | B | 100 pg/mL（各厂家一致） | POC 与实验室平台相同 |
| hsCRP / cCRP | B + H | <1 / 1–3 / >3 mg/L（AHA/CDC 2003）；FDA 指南要求 1.0 mg/L 处 CV ≤ 10% | 8 份文件无一有自有临床研究 |
| PCT | B + H | ICU 风险 0.5 / 2.0 ng/mL；LRTI/脓毒症抗生素决策 **0.10 / 0.25 / 0.50 ng/mL**；28 天死亡 **ΔPCT ≤ 80%** | 0.5/2.0 原型 K070310 未载建立方法；抗生素决策来自 RCT meta 分析 |
| 粪便钙卫蛋白 | C (+A) | 原型 **<50 阴性 / 50–120 灰区 / >120 µg/g 阳性**（n=908，120 处特异度 95%、灵敏度 70%）；BÜHLMANN 80/160；ALPCO 50/100 | 后续厂家多写 “clinical cut-off not applicable”，靠对 predicate 一致性 |
| D-dimer | B / C | VTE 排除 **500 ng/mL FEU**（VIDAS 原型）；HemosIL HS 230 ng/mL（ROC 建立）；Stratus CS 450 ng/mL FEU | 无任何产品获年龄校正声明 |
| CEA | A + G | 参考上限 3 / 5 ng/mL（非吸烟/吸烟）；Vista RCV **36.2%** | 后续修改型申报转引 1998–2004 原清关 |
| CA 125 | A + G | 30.8 / 35 / 38.1 U/mL（按厂家）；serial change 20%（= 1.645·√2·CV）或 25% | — |
| CA 19-9 | A + G | 35–37 U/mL；Vista RCV 84.7% | — |
| PSA（管理） | A + G | 4 ng/mL 参考上限；Siemens RCV 50%；FREND 20% | 筛查用途走 PMA |
| AFP（睾丸癌） | A + G | 参考上限按厂家；Vista RCV 33.7% | NTD 用途 Class III |
| HE4 | A + G | 绝经前 70 / 绝经后 140 pmol/L（95th）；serial change 14–25% | — |
| ROMA | E | **1.31 / 2.77**（Fujirebio/Abbott/Lumipulse，0–10 分制）；1.14 / 2.99（Roche） | 常见的 7.4%/25.3% 是非 FDA 标签表达 |
| Tg | B + G | **0.2 ng/mL** 三分层（ATA 2015）；Siemens 以 structural disease 为金标准 | 唯一直接采用指南值的肿瘤标志物 |
| HbA1c（诊断） | B + F | 6.5%；21 CFR 862.1373 特殊控制：TE ≤ 6%、5.0/6.5/8.0/12% 四水平、≥120 样本、HbF 黑框 | NGSP + IFCC 双溯源 |
| TSH | F (+A) | 参考区间约 0.3–4.x mIU/L，孕期分层；WHO **81/565** 溯源 | 排除 TPOAb/TgAb 阳性者 |
| 25-OH VitD | F | 无 cutoff；VDSP / Ghent RMP / NIST **SRM 972a** 溯源 | 20/30 ng/mL 属背景指南 |
| Ferritin | F | 性别分层参考区间；WHO 94/572 或 80/602 | 钩状效应必报 |
| Cystatin C | F | 无 cutoff；ERM-DA471 溯源；文件均无 eGFR 方程声明 | — |
| anti-CCP | A / C / H | Inova **20 CU**（99th）；IMMULITE 4、**Roche 17**、Euroimmun 5 RU/mL（ROC）；Abbott/ADVIA 5、BioPlex 3 U/mL（predicate） | 灵敏度 63.6–87.8%、特异度 96.5–98.6% |
| ANA | A / D | IFA 1:40 vs 1:80 同队列实测；NOVA View **48 LIU**（健康人 90th）；固相 ratio/AI ≈1.0 | 固相法以 ELISA 为 predicate，对 IFA 的 PPA/NPA 多未报 |
| ENA / dsDNA | D | 单项 index/IU/mL，dsDNA 10–100 IU/mL 各平台不可换算（校准品回收 27–121%） | dsDNA 以 Farr / CLIFT 为比对 |
| HSV-2 IgG | D | index 1.0，equivocal 0.9–1.1（各厂家略异） | 参考为 UW Western blot |
| Lyme | D + E | EIA index 1.0/1.1；WB **IgM 2/3、IgG 5/10** 条带；MTTT 双 EIA 串联 | 用 CDC 血清盘分期评价 |
| 梅毒 treponemal | D + E | **index 1.0**（Abbott S/CO、Bio-Rad AI）；TPPA/FTA-ABS 为参考 | 反向算法为背景 |
| CMV IgG / IgM | D | IgG index/AU 加 equivocal；IgM µ-capture | 限定“不用于血液/器官供者筛查” |

---

## 3. 分析性能验证标准对照表

FDA 决策摘要 “Standards/Guidance Documents Referenced” 一节里反复出现的标准，和它们各自对应的性能项目。CLSI 版本号以各 K 号文件为准（早期为 NCCLS / -A、-A2；近年为 -A3、Ed3）。

| 性能项目 | 标准 | 申报文件里的典型写法 |
|---|---|---|
| 精密度（重复性、实验室内） | CLSI EP05-A2 / A3 | 20 天 × 2 轮 × 2 重复，≥3–5 个浓度，含 cutoff 附近的浓度 |
| 线性 / 可报告范围 | CLSI EP06-A / EP06 Ed2 | 高低样本按比例混合 ≥ 9–11 级 |
| 干扰物 | CLSI EP07-A2 / Ed3 | 血红蛋白、胆红素、脂血、生物素、RF、HAMA、常用药物；判定标准通常为偏倚 ≤ 10% |
| 方法比对 / 与 predicate 对比 | CLSI EP09-A2 / A3 / c | Passing-Bablok 或 Deming 回归；n 通常 ≥ 100–200，覆盖测量范围 |
| 定性方法的一致性 | CLSI EP12-A / A2 | PPA / NPA / 总一致性 + 95% CI；cutoff ±20% 处的 C5/C95 |
| LoB / LoD / LoQ | CLSI EP17-A / A2 | hs-troponin 强调 LoD 和 “10% CV / 20% CV 处浓度”；EP17 定义 LoQ 的 CV 目标 |
| ROC / 定 cutoff | CLSI EP24-A2（原 GP10） | 自身免疫、钙卫蛋白 |
| 稳定性 | CLSI EP25-A；EN 13640 | 实时 + 加速；开瓶、校准周期 |
| 参考区间 | CLSI EP28-A3c（原 C28-A2/A3） | 非参数 2.5–97.5 百分位，n ≥ 120；hs-troponin 用 99th |
| 定性检测的性能 | CLSI EP12；FDA 2007《Statistical Guidance on Reporting Results from Studies Evaluating Diagnostic Tests》 | 比较方法不是金标准时用 PPA/NPA 而非敏感性/特异性；不允许 discrepant resolution；报告 95% CI；cutoff 须在临床研究前预设 |
| 校准溯源 | ISO 17511；具体标准品：WHO IS（TSH 81/565；ferritin 94/572）、IFCC/BCR/CAP CRM 470 / ERM-DA470k（CRP）、ERM-DA471（cystatin C）、NIST SRM 972a（25-OH VitD）、NGSP + IFCC（HbA1c）、NIST SRM 2921（cTnI） | 决策摘要 “Traceability” 小节 |
| 基质等效（血清 vs 各类血浆 vs 全血） | 无单独标准，按 EP09 做配对比对 | 常作为 §6 样本类型声明的依据 |
| 钩状效应（高剂量） | 无单独 CLSI 标准 | 肿瘤标志物、hCG、ferritin 必报 |
| CLIA waiver | FDA《Recommendations for Clinical Laboratory Improvement Amendments of 1988 (CLIA) Waiver Applications》（2008/2020） | POC 产品：非专业操作者研究、flex 研究 |

---

## 4. 临床验证设计的常见范式

决策摘要 “Clinical Studies” 小节的写法按 §2 原型分成四种模板：

1. **诊断准确度型**（B、C、D 原型）：前瞻性或回顾性收集预期用途人群（如急诊呼吸困难患者、疑似 RA 患者、STD 门诊），金标准是独立的临床判定或参考方法；报告敏感性/特异性/PPV/NPV 及 95% CI，按亚组（年龄、性别、分期、绝经状态）分层。若比较方法不是金标准，改报 PPA/NPA（FDA 2007 统计指南）。
2. **参考区间型**（A、F 原型）：不做患者研究，报告表观健康人群的百分位及其 90% CI，说明纳入/排除条件与样本量。
3. **cutoff 迁移 / 桥接型**（B 原型的后续厂家）：和已建立临床效用的 predicate 做方法比对，再在预设 cutoff 处报告分类一致性；FDA 2005 CRP 指南明确允许“文献 + 器械特异桥接研究”。
4. **监测一致性型**（G 原型）：在已确诊患者的连续样本中，把标志物变化（上升/下降/不变，按预设 % 变化）与临床状态变化（进展/缓解/无变化）做 3×3 或 2×2 表，报告一致性 %。

---

## 5. 样本类型（概览）

- **心脏标志物**：血清、Li-肝素血浆、EDTA 血浆三者常见；POC 产品（i-STAT、PATHFAST、Stratus CS）另含肝素/EDTA 全血。hs-cTn 的性别特异 99th 常按基质分别给出（如 Siemens 血清与肝素血浆各一套）。
- **PCT**：血清、EDTA/肝素血浆；**D-dimer**：仅 3.2% 枸橼酸血浆（Stratus CS 另验证肝素）；**钙卫蛋白**：粪便，各厂家提取装置不同（CALEX Cap、Q.S.E.T.、FED），提取步骤本身是申报验证项目。
- **肿瘤标志物**：血清为主，部分扩展到 EDTA/肝素血浆并做基质等效。
- **内分泌/代谢**：HbA1c 用 EDTA 全血与指尖血（POC/waived 产品）；TSH、VitD、ferritin、cystatin C 用血清与血浆。
- **自身免疫**：多数仅血清；少数验证 EDTA/肝素血浆。
- **感染血清学**：血清为主，多数验证 EDTA/肝素/枸橼酸血浆；Quidel Sofia Lyme 为指尖全血（CLIA waived）。

各靶点具体到抗凝剂、稳定性和基质等效研究见分节 §6。

---
## 目录（分节）

- **A. 心脏标志物**
  - 心肌肌钙蛋白 I / T（cardiac troponin, cTnI / cTnT；hs-cTn）
  - 利钠肽（NT-proBNP / BNP；B-type natriuretic peptide test system）
  - 高敏 / 心脏 C 反应蛋白（hsCRP / cardiac CRP, cCRP）
- **B. 感染 / 炎症**
  - 降钙素原（Procalcitonin, PCT）
  - 粪便钙卫蛋白（Fecal Calprotectin, fCAL）
  - D-二聚体（D-dimer）
- **C. 肿瘤标志物**
  - 癌胚抗原（CEA, Carcinoembryonic Antigen）
  - 糖类抗原 125（CA 125）
  - 糖类抗原 19-9（CA 19-9）
  - 前列腺特异抗原用于前列腺癌管理/监测（Total PSA for Management of Prostate Cancer）
  - 甲胎蛋白用于睾丸癌/生殖细胞肿瘤（AFP for Testicular / Germ Cell Cancer）
  - 人附睾蛋白 4（HE4）及 ROMA 算法
  - 甲状腺球蛋白用于分化型甲状腺癌监测（Thyroglobulin, Tg）
  - 附：本组横向对比速览（均取自上述文件）
- **D. 内分泌 / 代谢**
  - 糖化血红蛋白（HbA1c, Hemoglobin A1c / Glycated hemoglobin）
  - 促甲状腺激素（TSH, Thyroid-Stimulating Hormone / Thyrotropin）
  - 25-羟维生素 D（25-OH Vitamin D, 25(OH)D）
  - 铁蛋白（Ferritin）
  - 胱抑素 C（Cystatin C）
- **E. 自身免疫**
  - 抗环瓜氨酸肽抗体（anti-CCP / ACPA）
  - 抗核抗体筛查（ANA screen：IFA HEp-2 与固相法）
  - 可提取核抗原抗体与抗 dsDNA（ENA panel：dsDNA、Sm、RNP、SS-A/Ro、SS-B/La、Scl-70、Jo-1、centromere）
- **F. 感染性疾病血清学**
  - 单纯疱疹病毒 2 型特异性 IgG（HSV-2 type-specific IgG；兼看 HSV-1 IgG）
  - 莱姆病伯氏疏螺旋体抗体（*Borrelia burgdorferi* antibodies, Lyme disease serology）
  - 梅毒螺旋体抗体（*Treponema pallidum* treponemal antibodies）
  - 巨细胞病毒 IgG / IgM（Cytomegalovirus, CMV IgG / IgM）


# 第 A 组 · 心脏标志物

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

# 第 B 组 · 感染 / 炎症

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

# 第 C 组 · 肿瘤标志物

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

# 第 D 组 · 内分泌 / 代谢

> 数据来源说明：本节全部申报数据均摘自 FDA 公开的 510(k) / De Novo 决策摘要（decision summary，`cdrh_docs/reviews/`）或申请人 510(k) summary（`cdrh_docs/pdfYY/`）。文件中没写的信息一律标「文件未载明」。凡来自背景知识（临床指南、生理学）的解释，单独标注「背景（非申报文件）」。openFDA 清关计数为 `python3 fda_fetch.py list <CODE>` 首行 `total=`（查询日期 2026-09-22）。

---

## 糖化血红蛋白（HbA1c, Hemoglobin A1c / Glycated hemoglobin）

### 1. 法规定位

HbA1c 在 FDA 体系中有**两条平行的法规路径**，这是本靶点最核心的要点：

| 预期用途 | Product code | 21 CFR | Class | Panel | openFDA 累计清关 |
|---|---|---|---|---|---|
| 糖尿病**监测**（long-term glycemic control） | **LCP** – Glycosylated Hemoglobin Assay | **864.7470** | II | HE – Hematology (81)（早期个别文件写 Chemistry 75，如 K110313） | **249** |
| 糖尿病**诊断 / 风险识别**（aid in diagnosis; aid in identifying patients at risk） | **PDJ** – Hemoglobin A1c Test System | **862.1373**（Class II, special controls；2014-08-25 生效，源自 Roche De Novo DEN130002） | II | CH – Chemistry (75) | **26** |

- 21 CFR 862.1373 的特殊控制（special controls）原文要点（来自 eCFR 全文）：
  1. 须获得 FDA 认可的糖化血红蛋白标准化组织的**初始及每年一次**的标准化验证（实际即 NGSP 认证）；
  2. 510(k) 须包含精密度、准确度、线性、干扰的性能测试：
     - (i) 精密度至少用 **≈5.0%、6.5%、8.0%、12% HbA1c** 四个水平的血样，≥20 天、≥3 个试剂批、≥3 台仪器；
     - (ii) 准确度至少 **120 份**跨测量区间的血样，与标准化方法（standardized test method）比对，须显示“little or no bias”；
     - (iii) **总误差（total error）以单次测量对标准化方法评价，须 ≤ 6%**；
     - (iv) 须证明对常见血红蛋白变异体 **HbC、HbD、HbE、HbA2、HbS** 几乎无干扰；
  3. 若观察到 **HbF** 或低频变异体的干扰，须在所有标签材料中以**黑框警示（black box）**说明干扰及受影响人群。
- 与之对应，LCP 监测类产品的决策摘要一律写「Assay cut-off: Not applicable」「Clinical cut-off: Not applicable」，并在特殊使用条件中明示「not intended for use in the diagnosis of or screening for diabetes」（K221326、K163633、K121842、K140829）。
- POC 产品常同时携带仪器代码 JQT（21 CFR 862.2400 反射式密度计，Class I，如 Afinion、Allegro）或 JJE（862.2160 离散光度分析仪，如 cobas b 101、DxI 9000）。
- 注：DEN130002（Roche Tina-quant HbA1c Gen.2，2013-05-23，PDJ 首个产品）在 openFDA 有记录，但 `cdrh_docs/reviews/DEN130002.pdf` 返回 404，本节未能读取其正文；诊断用途的规则内容改用 K121610（Roche cobas c 501 Tina-quant HbA1cDx Gen.3，2013-08-08，PDJ 第 2 个清关）与 eCFR 条文。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K121610 | 2013-08-08 | Roche / cobas c 501 Tina-quant HbA1cDx Gen.3 | 浊度抑制免疫比浊 (TINIA)，TTAB 溶血 | 全血：Li-heparin、K2/K3-EDTA、KF/Na2-EDTA、Na-heparin、NaF/K-oxalate、NaF/Na2-EDTA | **aid in diagnosis of diabetes；identify at-risk**（PDJ, 862.1373） | 决策摘要 + 510(k) summary |
| K180296 | 2018-05-07 | Alere (Abbott) / Afinion HbA1c Dx (AS100) | 硼酸亲和法 + 反射光度 | K2-EDTA 静脉全血、指尖毛细血管全血 | **aid in diagnosis；at-risk；POC（moderate complexity）**（PDJ） | 决策摘要 + 510(k) summary |
| K214117 (+CW210007) | 2023-09-27 | Abbott Diagnostics Technologies AS / Afinion HbA1c, Afinion 2 & AS100 | 硼酸亲和法 | 静脉及毛细血管全血 | **monitoring**（marker of long-term metabolic control）；**Dual 510(k)+CLIA waiver**；Rx；POC | 决策摘要 + 510(k) summary |
| K221326 (+CW220003) | 2024-11-22 | Nova Biomedical / Nova Allegro HbA1c Assay + Allegro Analyzer | 乳胶增强免疫比浊 | **仅指尖毛细血管全血**（1.5 µL） | **monitoring**；**Dual 510(k)+CLIA waiver**；明示不用于诊断/筛查；不用于新生儿 | 决策摘要 + 510(k) summary |
| K163633 | 2017-07-28 | Roche / cobas HbA1c Test, cobas b 101 | 乳胶凝集抑制免疫法（单次使用 disc） | 指尖或静脉全血：EDTA (K2/K3)、Li-heparin | monitoring；POC + 实验室；不用于筛查/诊断/新生儿 | 决策摘要 + 510(k) summary |
| K121842 | 2012-12-12 | Axis-Shield (Abbott) / ARCHITECT HbA1c (CMIA) | 化学发光微粒子免疫法（磁性硅微粒） | 全血：K2-EDTA、NaF/K-EDTA、NaF/Na-EDTA、fluoride oxalate（Li-heparin 不可） | monitoring；「Should not be used for the diagnosis of diabetes mellitus」 | 决策摘要 + 510(k) summary |
| K140829 | 2014-07-14 | Beckman Coulter / UniCel DxC Synchron HbA1c3 Reagent | 免疫比浊抑制（A1c3）+ 比色总 Hb（Hb3） | 全血：K2/K3-EDTA、Li-/Na-heparin | monitoring；「Do not use for screening or diagnosis」 | 决策摘要 + 510(k) summary |
| K110313 | 2011-12-23 | Roche / Tina-quant HbA1c Gen.2 (Integra 800) | TINIA；全血或溶血液两种应用 | 全血/溶血液，7 种抗凝剂 | monitoring；本次为增加 3 种抗凝剂 + LoB/LoD 声明 | 决策摘要 + 510(k) summary |
| K071132 | 2008-09-16 | Tosoh / G8 Automated Glycohemoglobin Analyzer HLC-723G8 | 阳离子交换 HPLC | EDTA 全血 | monitoring（clinical management of diabetes） | 决策摘要（无 summary） |
| K130990 / K151809 | 2013-05-09 / 2015-09-25 | Bio-Rad VARIANT II TURBO 2.0 / Alere Afinion（取样器玻璃→塑料） | HPLC / 硼酸亲和 | — | Special 510(k) 器械修改，无性能数据 | 决策摘要（仅审查备忘录） |

（未读取但已在列表中确认的其他 PDJ 诊断类清关：Bio-Rad D-100 K151321、Tosoh G8 K131580/K200904、Sebia CAPILLARYS K171861、Siemens ADVIA A1c_E K171771、Beckman HbA1c Advanced K182651、Abbott K130255/K140654、Ortho VITROS K142595 等。）

### 3. 预期用途与声明类型

- **监测（LCP）**：典型措辞「The measurement of % HbA1c is recommended as a marker of long-term metabolic control in persons with diabetes mellitus」（K214117）、「used for the monitoring of long-term blood glucose/metabolic control in individuals with diabetes mellitus」（K221326）、「monitoring long term glycemic control in diabetic patients」（K121842）。附带的限制语句：不用于诊断/筛查；不用于评估每日血糖控制、不能替代家庭血糖/尿糖监测；不用于红细胞寿命缩短者（溶血病、妊娠、显著急/慢性失血；K121842 还列出总 Hb <7 或 >21 g/dL 不可用）。
- **诊断（PDJ）**：「to be used as an aid in the diagnosis of diabetes and as an aid in identifying patients who may be at risk for developing diabetes」（K121610、K180296）。诊断类文件附带更长的禁忌清单（K180296）：不用于妊娠期糖尿病、HbF >10%（HPFH）、血红蛋白病但红细胞周转正常者（如镰状细胞性状）、红细胞周转异常者（溶血/缺铁性贫血）、遗传性球形红细胞增多症、恶性肿瘤、严重慢性肝肾病、3 周内输血或化疗；快速进展的 1 型糖尿病中 HbA1c 升高滞后于血糖，须以血糖诊断；不能替代 1 型、儿童、孕妇的血糖检测。
- **POC / CLIA waived**：K214117（Abbott Afinion）与 K221326（Nova Allegro）均为 **Dual 510(k) + CLIA Waiver by Application**（CW210007 / CW220003）。K180296 Afinion HbA1c Dx 为 POC 但在 **moderate complexity** 实验室场景验证（与预定 Roche Dx Gen.3「clinical laboratories only」相比，差异项即「intended use sites: point of care settings」）。K163633 cobas b 101 为「professional use in a clinical laboratory setting or PoC」。
- 全部为 **Rx（处方用）**；K221326 为 Rx 且「Intended Users: Prescription use. Waived users」。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

**A. 监测类（LCP）——无 cutoff，仅给期望值/治疗目标**

| 产品 | 「Assay cut-off」 | 期望值/参考区间（申报文件原文） | 建立方式 |
|---|---|---|---|
| K214117 Afinion | Not applicable | 标签引用 ADA：非妊娠成人合理目标 HbA1c <7.0% (53 mmol/mol)，目标须个体化（ADA Diabetes Care 2019;42 Suppl 1:S61–S70） | 引用临床指南，非自建 |
| K221326 Nova Allegro | Not applicable | 标签引用 ADA 2023：<7% (53 mmol/mol)；宽松目标 <8% (64 mmol/mol)（Diabetes Care 2023;46 Suppl 1:S97–S110）；另计算 eAG (mg/dL)=28.7×A1C−46.7（NGSP 公式） | 引用指南 |
| K163633 cobas b 101 | Not applicable | 引用 ADA 2016：<7%；严格 6.5% (48)；宽松 8% (64) | 引用指南 |
| K121842 ARCHITECT | Not applicable | ADA 2012：<8% 宽松、<7% 一般、<6.5% 严格；**5.7–6.4% (39–46 mmol/mol) 为糖尿病风险增加**；predicate AxSYM 的参考区间研究中央 95% 为 4.6–6.0% | 引用指南 + predicate 参考区间 |
| K140829 Beckman DxC | Not applicable | 成人正常 HbA1c 4.0–6.0% (NGSP)；20–42 mmol/mol (IFCC)，引用 Panteghini 2007、Tietz、Henry's；建议各实验室自建 | 文献 |
| K110313 Roche Gen.2 | Not applicable | IFCC 29–42 mmol/mol；DCCT/NGSP 4.8–5.9%（Junge W 等，EUROMEDLAB 2003 海报） | 文献 |
| K071132 Tosoh G8 | Not applicable | 146 名非糖尿病健康成人（代表美国人群）EDTA 全血：4.4–6.1% | 自建（n=146；纳入排除细节文件未载明） |

**B. 诊断类（PDJ）——以临床指南既定值为 cutoff**

- K180296 Afinion HbA1c Dx：「The diagnostic cut-off is **6.5% HbA1c**. Patients with HbA1c values in the range **5.7–6.4%** are identified as having an increased risk for developing diabetes」，引用 International Expert Committee 2009 与 ADA 2018。决策摘要「Clinical cut-off: See expected values」。
- K121610 Roche Dx Gen.3：参考区间 IFCC 20–42 mmol/mol、DCCT/NGSP 4.0–6.0%；「According to the recommendations of the ADA values above 48 mmol/mol (6.5%) are suitable for the diagnosis of diabetes mellitus. Patients with 39–46 mmol/mol (5.7–6.4%) may be at a risk of developing diabetes」；治疗行动建议 >64 mmol/mol (8%)；<53 mmol/mol (7%) 达 ADA 目标。
- 逻辑：cutoff **不是**由申报方自建 ROC/参考人群得出，而是直接采用 ADA/IEC 既定值；申报方需要证明的是**在 cutoff 附近（5.0/6.5/8.0/12%）的精密度、偏倚与总误差 ≤6%**（见第 8 节），并通过 NGSP 认证/IFCC 溯源保证与 DCCT 建立的 cutoff 可比。
- 单位与换算：K180296 结果以 IFCC mmol/mol 计算再用 IFCC Master Equation 换算：NGSP-HbA1c (%) = 0.09148 × IFCC (mmol/mol) + 2.152；K121610 写为 IFCC = (NGSP − 2.15)/0.092。

**C. 溯源**

- 所有产品声明 **NGSP 认证（一年有效，须每年更新）+ IFCC 参考方法溯源**（K214117、K221326、K163633、K121842、K140829、K180296、K121610、K071132）。K180296 另说明「Traceability is verified annually through the IFCC Whole Blood Monitoring Program」。K071132 Tosoh 通过日本 ReCCS 参考物质赋 JDS 值，再经 IFCC-JDS 关系与 IFCC-NGSP Master Equation 对齐。K140829 校准品定值使用 IFCC calibrators/controls 及 NGSP 溯源样本。

### 5. 生物学 / 生理学依据

- 申报文件中的表述：HbA1c 为 β 链 N 端缬氨酸糖化的血红蛋白（K163633：糖化 HbF 不含表征 HbA1c 的糖化 β 链故不被检出，但 HbF 计入总 Hb，导致 HbF >10% 时 %HbA1c 偏低）；HbA1c 反映既往 2–3 个月血糖（K121610：高于参考区间上限提示既往 2–3 个月高血糖；低于参考区间可能提示近期低血糖、Hb 变异体或红细胞寿命缩短）。红细胞寿命缩短（溶血病、失血、缺铁、妊娠）会使 HbA1c 不能反映血糖（多个文件的限制语句）。K140829 提示红细胞沉降（ESR 升高）与未混匀全血可致错误结果。
- 背景（非申报文件）：ADA 采用 6.5% 诊断阈值系基于视网膜病变患病率拐点的流行病学数据（IEC 2009）；NGSP 通过 DCCT/UKPDS 结果与 %HbA1c 的关系保证临床解释一致；IFCC 参考方法（HPLC–ESI-MS / HPLC–CE）测量 β-N-1-deoxyfructosyl-Hb 六肽，报告 mmol/mol。

### 6. 样本类型与样本要求

| 产品 | 样本类型 | 抗凝剂 / 基质等效性研究 | 稳定性 / 其他 |
|---|---|---|---|
| K180296 Afinion Dx | K2-EDTA 静脉全血；指尖毛细血管全血 | 「Matrix comparison: Not applicable. Venous K2-EDTA and capillary are the only matrices claimed」；方法比对同时用指尖血与配对静脉血 | 凝固或溶血样本不可用；>14% 溶血 (2000 mg/dL) 报信息码；Hb <6 或 >20 g/dL 不出结果 |
| K221326 Nova Allegro | **仅指尖毛细血管全血**（1.5 µL，采样后 1 min 内上机） | 精密度以 K2-EDTA 静脉血作为验证过的替代样本，并提供数据证明静脉血与毛细血管血行为一致；总 Hb 6–20 g/dL 无干扰，超出报「HbA1c Bad Sample」 | 海拔 12,000 ft 无影响 |
| K163633 cobas b 101 | 指尖血、静脉全血 K2/K3-EDTA、Li-heparin | 91 份 K2 vs K3-EDTA 配对：y=1.03x−0.00, r=0.99；方法比对分别用指尖血、K2-EDTA、Li-heparin | 「Use fresh whole blood only. Do not use plasma and serum」；disc 开袋后 20 min 内使用 |
| K121842 ARCHITECT | K2-EDTA（对照）、NaF/K-EDTA、NaF/Na-EDTA、fluoride oxalate | 20 非糖尿病 + 20 糖尿病配对管：斜率 1.03–1.07，r=0.99；**Li-heparin 不可用**（标签明示） | — |
| K140829 Beckman | K2-EDTA（参考）、K3-EDTA、Li-/Na-heparin | Deming：K3-EDTA y=1.012x−0.043；Li-hep y=1.007x−0.034；Na-hep y=1.007x−0.034（n=61–62, r=0.999） | 校准品复溶后 15–25 °C 8 h / 2–8 °C 48 h；分装冻存 60 天 |
| K110313 Roche Gen.2 | 7 种抗凝剂（新增 Na-heparin、NaF/K-oxalate、NaF/Na2-EDTA），全管与半满管均验证 | Passing-Bablok 相对 K2-EDTA 斜率 0.994–1.010（n=132, 4.97–15.97%） | 全血 8 h 内检测，否则制备溶血液：15–25 °C 4 h、2–8 °C 24 h、−15~−25 °C 6 个月；新抗凝剂冻存稳定性未确定 |
| K121610 Roche Dx Gen.3 | 7 种抗凝剂（沿用 k102914） | — | TTAB 溶血剂不裂解白细胞，无需去除不稳定 HbA1c 预处理 |
| K071132 Tosoh G8 | EDTA 全血（4 µL） | Not applicable | — |

### 7. 分析性能验证所依据的标准

| K 号 | 文件列出的标准（Standards/Guidance Documents Referenced） | 对应性能项目 |
|---|---|---|
| K180296 | CLSI EP5-A3（精密度）、EP7-A2（干扰）、ISO 14971:2007 / EN ISO 14971:2012（风险管理） | 精密度按 EP05-A3（内部 + 3 个 POC 站点）；干扰按 EP7-A2；总误差按 862.1373 特殊控制公式 |
| K121610 | CLSI EP5-A2（精密度）、EP6-A（线性）、EP17-A（LoB/LoD） | 精密度 3 仪器 × 3 批 × 21 天；线性 20 级稀释系列；LoB/LoD |
| K163633 | CLSI EP05-A3、EP06-A、EP7-A2、EP9-A3、EP25-A（稳定性）、IEC 62304、IEC 61010-1、IEC 60601-1-2 | 精密度、线性、干扰、方法比对、试剂稳定性、软件/电气安全 |
| K121842 | CLSI EP9-A2、EP6-A、EP5-A2、EP7-A、「EN 13460 06: Stability Testing of In Vitro Diagnostic Reagents」（文件原文；对应标准应为 EN 13640）；正文另引用 EP17A（LoB/LoD/LoQ） | 方法比对、线性、精密度、干扰、稳定性、检出限 |
| K140829 | CLSI EP05-A2、EP06-A、EP07-A2、EP09-A2、EP17-A2、**EP28-A3**（参考区间） | 精密度、线性、干扰、方法比对、检出限、参考区间 |
| K110313 | CLSI EP17-A | LoB/LoD |
| K214117 | ISO 14971（仅风险管理；性能沿用 K050574） | — |
| K221326 | 决策摘要「None referenced」；510(k) summary 提及 CLSI EP5-A3 | 精密度 |
| K071132 | UL 61010-1、IEC 60601-1-2、EN 980、EN 375（标签/电气） | 无 CLSI 标准列出 |
| 溯源 | 全部：NGSP 认证 + IFCC 参考方法（K180296 另引 IFCC Whole Blood Monitoring Program；K071132 引 JDS/ReCCS） | 溯源 |

### 8. 临床/准确度验证设计与结果

**K180296 Afinion HbA1c Dx（PDJ，POC 诊断）**
- 精密度（EP05-A3，特殊控制要求的 4 个水平）：内部 4 份 K2-EDTA 静脉血（≈5, 6.5, 8, 12%），11 名操作者，3 批试剂，9 台 AS100，20 天，每份 720 次：合并总 CV **1.64% (5.14%)、1.51% (6.55%)、1.31% (8.06%)、1.42% (11.26%)**。外部 3 个 POC 中等复杂度站点、6 台仪器、3 批、10 天，每站点每份 240 次：总 CV 1.78 / 1.48 / 1.36 / 1.22%。指尖血精密度由三项研究合成：总 CV 2.03% (低)、1.58% (阈值 6.5%)、1.49% (中)、1.30% (高)。
- 准确度：120 份配对指尖血 + K2-EDTA 静脉血（4.6–11.4%），2 批试剂，3 个 POC 站点；样本分布刻意围绕决策点（6.1–6.5% 占 25%，6.6–7% 占 25%）；比对方法 **NGSP comparator Tosoh G8 HPLC**（静脉血双份均值）。指尖血：加权 Deming 斜率 0.997 (95% CI 0.966–1.027)、截距 0.000；Passing-Bablok 斜率 1.000、截距 −0.040。静脉血：Deming 斜率 0.991、截距 0.053；PB 斜率 1.000、截距 −0.030。
- 决策水平偏倚：指尖血 Deming 在 5.0/6.5/8.0/12.0% 处 %bias 均约 −0.33%；静脉血 +0.20/−0.05/−0.21/−0.43%。
- **总误差**：%TE = |%Bias| + 1.96 × %CV × (1 + %Bias/100)。指尖血 Deming：**4.30% (5.0%)、3.42% (6.5%)、3.24% (8.0%)、2.87% (12%)**；PB：4.75/3.69/3.40/2.87%。静脉血（外部精密度）：3.69/2.95/2.87/2.81%。全部 <6%。
- Hb 变异体：221 份 K2-EDTA 样本（HbA2 26、HbS 21、HbC 25、HbE 20、HbD 21、HbF 121），比对方法 Premier Hb9210（HbA2/S/E）与 Tosoh G8（HbF）；显著干扰阈值 ≤7% 平均相对偏差。平均 %bias（~6.5% / ~8.5%）：HbA2 −3.4/−2.5；HbS −4.1/−1.0；HbC −5.5/−1.8；HbE +3.5/+3.7；HbD −2.2/−2.9。**HbF：10.4% 为无显著干扰的最高浓度，超过则负干扰** → 黑框警示。

**K121610 Roche cobas c 501 Tina-quant HbA1cDx Gen.3（PDJ，实验室诊断）**
- 精密度：4 份静脉全血（5, 6.5, 8, 12%）+ 2 控制品，3 台 c 501 × 3 批 × 21 天；合并总 CV 1.9% (5.05%)、1.7% (6.40%)、2.1% (7.99%)、2.0% (11.34%)。
- 准确度：**141 份无变异体样本**（4.7–12.2%），3 天单次测定，比对 **NGSP 二级参考实验室 HPLC**；样本分布 6–6.5% 占 19.9%、6.5–7% 占 23.4%。Deming 斜率 1.002 (0.974–1.0285)、截距 −0.114；PB 斜率 1.006、截距 −0.135。决策水平偏倚：5.2% −1.98%；6.5% −1.45%；8.0% −1.06%。
- **总误差**：5.2% → **6.0%**（%CV 2.07）；6.5% → 4.7%；8.0% → 5.1%（此文件公式写作 %TE = |%Bias| + 1.96×%CV×(1+%Bias)）。
- Hb 变异体：116 份（HbS 20、HbC 19、HbE 20、HbD 20、HbF 20、HbA2 17）；相对偏倚（~6% / ~9%）：C −3.07/−0.35；S +2.17/+3.42；E −1.58/+3.46；D −2.30/+3.35；A2 −5.73/−4.12；**HbF >7% 时偏倚超过 −7%，负偏倚与 %HbF 成正比** → 黑框警示 + 「HbSS/HbCC/HbSC 疑影响时不得用于诊断」。
- 线性 4.2–20.2%；LoB 2.3%、LoD 2.5%；干扰接受标准 ±7%（脂血 600 mg/dL、胆红素 60 mg/dL、RF 750 IU/mL、葡萄糖 1000 mg/dL、总蛋白 21 g/dL 无干扰）。

**K221326 Nova Allegro（LCP，CLIA waived，指尖血）**
- 方法比对：4 个医师办公室、15 名代表 CLIA waived 用户的操作者、533 名患者指尖血 vs 配对 K2-EDTA 静脉血在 **NGSP 二级参考实验室**的 FDA 已清关比对法；8 台仪器、3 批；PB（n=526，4.4–13.8%）合并斜率 0.972、截距 0.217、r 0.993；各站点斜率 0.968–0.983。
- 精密度（waived 用户）：20 天控制品 4 站点合并再现性 CV 2.00% (5.67%)、2.55% (9.47%)；5 天静脉血单站点 within-lab 1.55–2.71%；指尖血重复性（524 例第二次采血）1.50% (4.0–6.0%)、1.51%、1.92%、2.38% (10.1–14.0%)；合成再现性 2.31/1.87/2.45/2.91%。
- 线性：11 份 4.0–15.4% 与 Tosoh G8 比：斜率 1.016、截距 −0.14、r 0.999；声明 4–14%。
- 变异体：119 份（HbC/E/D/F 各 20、HbS 25、HbA2 14）；平均偏倚 C +1.9%、E +4.4%、D +3.0%、A2 +0.6%、S −4.5%；**HbF 仅 5.4% 以内无显著干扰（<10%）** → 标签警示「HbA1c results are invalid for patients with abnormal amounts of HbF」。干扰接受标准为 ±10%（与诊断类 ±7% 不同）。

**K163633 Roche cobas b 101（LCP，POC）**
- 精密度（EP05-A3）：内部 1 台 × 2 批 × 21 天 n=84，总 CV 1.2–2.9%；外部 3 个 POC 站点、6 台、3 批，合并总 CV 1.3–2.7%（全血）、3.1–3.7%（控制品）。
- 方法比对（EP09-A3）：3 个 POC 站点、n=379 前瞻性样本，指尖血/K2-EDTA/Li-heparin 分别 vs 配对 K2-EDTA 在 **NGSP 二级参考实验室 Tosoh G8**：各站点 PB y=1.00x−0.10~−0.20（K2-EDTA 站点 3 为 0.97x−0.04），r 0.99。
- 变异体：130 份（S/C/D/E/F 各 20、A2 10），比对 Tina-quant Gen.3 (k121610) 与 Tosoh G8；接受标准 ±10%；HbS ≤41%、HbC 36%、HbD 42%、HbE 27%、HbA2 6.2% 无干扰；HbF >10% 致偏低。
- 线性 3.6–12.9%：斜率 0.996、截距 −0.014、r 0.9961；声明 4–12%。

**K121842 Abbott ARCHITECT HbA1c**
- 精密度（EP5-A2）2 台 × 2 批 × 20 天 n=80：总 CV 2.8–3.8%（5.1–9.9%）。
- 方法比对 vs predicate AxSYM（EP9-A2）：n=127 EDTA 全血 4.07–13.61%，r 0.95 (0.93–0.96)，斜率 1.04 (0.97–1.12)，截距 −0.07。
- LoB 2.7%、LoD 2.8%、LoQ 4.0%（EP17-A）；线性 y=0.96x−0.09, r²=0.98（4.0–15.8%）。
- 变异体：HbA2/C/D/E/F/J/S 样本 4.7–11.3% 全部 <10% 偏倚，但标签仍写 **HbD、HbE、HbF (>9%)、HbS 会干扰**；总 Hb 7–20 g/dL 无干扰。

**K140829 Beckman DxC HbA1c3**
- 精密度（EP05-A2）DxC 600/800 各 n=80：总 CV 1.34–1.93%。
- 方法比对 vs predicate（Deming）：DxC 600 n=119（4.4–16.6%）r 0.998、斜率 1.031、截距 −0.267；DxC 800 n=118 r 0.999、斜率 1.031、截距 −0.294。
- 变异体 vs HPLC 参考：C ≤36.4%、D ≤35.1%、E ≤21.7%、S ≤31.9% 无干扰（接受 ±7%）；HbF >10% 致偏低（接受 ±10%）；内源干扰接受标准 ±6%。

**K071132 Tosoh G8（HPLC）**
- 精密度 within-run/between-run CV 0.39–1.30%；线性 4.0–16.9%（回收 97.9–103.2%）；方法比对 vs G7 n=114：斜率 1.020、截距 −0.16、r 0.998。
- 变异体：HbAD/AS/AC 在 A0 峰后洗脱并被扣除，≤30% 无干扰；HbF ≤10% 无干扰；**HbAE 无法与其他峰区分，故干扰 sA1c**。

**K214117 Afinion（2023 标签更新）**：精密度/线性/方法比对均「Previously established in K050574」；本次新增干扰数据（≤7% 标准）：HbA2 ≤5.7%、HbS ≤42%、HbC ≤36%、HbE ≤26%、HbD ≤42%、HbF ≤10.4% 无干扰，>10.4% HbF 负干扰；乙酰化 Hb 4.6 mg/mL、氨甲酰化 Hb 13.8 mg/mL、不稳定 HbA1c 11.4 mg/mL、糖化白蛋白 7.7 mg/mL 无交叉反应。

### 9. 厂家间差异与要点

1. **监测 vs 诊断的分水岭在于样本设计与总误差**：诊断类（PDJ）必须用 5.0/6.5/8.0/12% 四个精密度水平、≥120 份围绕 6.5% 富集的样本、与 NGSP 二级参考实验室 HPLC 比对、并计算 %TE ≤ 6%；监测类（LCP）只需与 predicate 或 NGSP 参考法回归比对，决策摘要不计算 TE。
2. **干扰接受标准**：诊断类文件多采用 ±7%（K180296、K121610、K214117 因沿用 Dx 数据），监测类多为 ±10%（K221326、K163633、K121842），Beckman 内源干扰用 ±6%、变异体 ±7%。
3. **HbF 阈值差异显著**：Nova Allegro 5.4%、Roche Dx Gen.3 7%、ARCHITECT 9%、Afinion 10.4%、cobas b 101/Beckman/Tosoh G8 10%。免疫法（抗体针对糖化 β 链 N 端）与硼酸亲和法对 HbF 均为负干扰（HbF 计入总 Hb 分母）。
4. **HbE**：Tosoh G8 HPLC 明确写 HbAE 干扰；ARCHITECT 标签写 HbD/E/S 干扰；硼酸亲和法与 TINIA 对 S/C/D/E/A2 均声明无干扰。
5. **样本**：Nova Allegro 仅声明指尖血（静脉血作为精密度替代）；Afinion 同时声明静脉 K2-EDTA + 指尖血；Roche/Beckman 实验室型声明多种抗凝剂并做管型等效性；ARCHITECT 明确排除 Li-heparin。
6. **CLIA waiver 研究设计**（K221326）：4 个 POC 站点、每站点 ≥3 名未经训练/waived 操作者、≥90 例/站点、≥20 天、8 台仪器、3 批，指尖血与中心实验室 NGSP 方法比对 + 第二次指尖采血评估重复性；20 天控制品精密度由 waived 操作者完成。
7. **参考区间**：监测类产品普遍不自建，而引用 ADA 目标或文献；自建者仅 Tosoh G8（n=146）与 AxSYM predicate（中央 95%）。

### 10. 来源

- K121610: https://www.accessdata.fda.gov/cdrh_docs/reviews/K121610.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K121610.pdf
- K180296: https://www.accessdata.fda.gov/cdrh_docs/reviews/K180296.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K180296.pdf
- K214117: https://www.accessdata.fda.gov/cdrh_docs/reviews/K214117.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K214117.pdf
- K221326: https://www.accessdata.fda.gov/cdrh_docs/reviews/K221326.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K221326.pdf
- K163633: https://www.accessdata.fda.gov/cdrh_docs/reviews/K163633.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K163633.pdf
- K121842: https://www.accessdata.fda.gov/cdrh_docs/reviews/K121842.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K121842.pdf
- K140829: https://www.accessdata.fda.gov/cdrh_docs/reviews/K140829.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K140829.pdf
- K110313: https://www.accessdata.fda.gov/cdrh_docs/reviews/K110313.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf11/K110313.pdf
- K071132: https://www.accessdata.fda.gov/cdrh_docs/reviews/K071132.pdf
- K130990: https://www.accessdata.fda.gov/cdrh_docs/reviews/K130990.pdf ；K151809: https://www.accessdata.fda.gov/cdrh_docs/reviews/K151809.pdf
- 21 CFR 862.1373 全文：https://www.ecfr.gov/current/title-21/section-862.1373 （DEN130002 决策摘要 URL 404，未能读取）

---

## 促甲状腺激素（TSH, Thyroid-Stimulating Hormone / Thyrotropin）

### 1. 法规定位

- Product code **JLW**（Thyroid stimulating hormone test system），**21 CFR 862.1690**，Class II，Panel CH – Clinical Chemistry (75)。
- openFDA 累计清关 **total = 242**。
- 常见伴随代码：JIS/JIT（校准品，862.1150）、JJX（质控品，862.1660）、JJE（分析仪，862.2160）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K251543 | 2026-02-06 | Siemens / Atellica IM TSH3-Ultra II (TSH3ULII) | 夹心化学发光（吖啶酯；FITC/抗 FITC 磁微粒） | 血清、血浆（EDTA、Li-heparin） | diagnosis of thyroid or pituitary disorders；「third-generation assay」；Rx | 决策摘要 + 510(k) summary |
| K233050 | 2024-04-04 | Siemens / ADVIA Centaur TSH3-Ultra II | 同上 | 血清、血浆（EDTA、Li-heparin） | 同上；含儿科参考区间 | 决策摘要 + 510(k) summary |
| K221225 | 2022-11-10 | Beckman Coulter / Access TSH (3rd IS) on DxI 9000 | 双位点免疫酶法（碱性磷酸酶 + Lumi-Phos PRO） | 血清、血浆（Li-heparin） | 「capable of providing 3rd generation TSH results」；Rx | 决策摘要 + 510(k) summary |
| K153651 | 2016-08-18 | Beckman Coulter / Access TSH (3rd IS) + Calibrators (DxI 800) | 同上（Lumi-Phos 530） | 血清、血浆（Li-heparin；血清含/不含凝胶） | 同上；本次重点为改用 WHO 3rd IS 81/565 标准化 | 决策摘要 + 510(k) summary |
| K162606 | 2017-01-23 | Roche / Elecsys TSH on cobas e 801 | 电化学发光 ECLIA（生物素-链霉亲和素，钌标记） | 血清（含分离胶）、Li-heparin、K2/K3-EDTA | diagnosis of thyroid or pituitary disorders | 决策摘要 + 510(k) summary |
| K130469 | 2013-04-05 | DiaSorin / LIAISON TSH + Control Thyroid 1/2/3 (LIAISON XL) | 夹心化学发光（异鲁米诺） | **仅血清** | 同上 | Triage 决策摘要（无内容）+ 510(k) summary |

（K190773 Roche Elecsys TSH 2019：决策摘要为 Triage-Quick Review，510(k) summary 为 PDF portfolio 无法抽取文本，未采用。K083173 Ortho VITROS TSH 2008 为 Special 510(k) 加装 VITROS 3600，无性能数据。）

### 3. 预期用途与声明类型

- 统一措辞：「Measurements of thyroid stimulating hormone produced by the anterior pituitary are used in the diagnosis of thyroid or pituitary disorders」（K251543、K233050、K221225、K162606、K130469）→ **辅助诊断**（甲状腺/垂体疾病），非筛查、非监测专用声明。
- **“第三代”声明**：Beckman 在 Indications 中直接写「This assay is capable of providing 3rd generation TSH results」（K153651、K221225），其 predicate Access HYPERsensitive hTSH 曾声明「3rd generation (HYPERsensitive) and/or 2nd generation (Fast hTSH)」；Siemens 在原理段落写「This assay is a third-generation assay」（K233050、K251543）；Roche、DiaSorin 文件未使用“代”的措辞。
- 全部 Rx，实验室用；无 POC/waived 产品在本次样本中（列表中有 Boditech AFIAS、Nanoentek FREND、Qualigen FastPack 等 POC 型，未读取）。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

- 所有文件「Assay cut-off: Not applicable」「Clinical cut-off: Not applicable」；判定依赖**参考区间**。

| 产品 | 参考区间（申报文件） | 人群与建立方法 |
|---|---|---|
| K153651 Beckman Access TSH (3rd IS)（K221225 沿用） | **一般成人 0.45–5.33 µIU/mL**（中位 1.48，观察范围 0.32–7.08）；孕早期 0.05–3.70；孕中期 0.31–4.35；孕晚期 0.41–5.18 | 每人群约 367 例（一般人群男女约各半、21–88 岁；孕妇按三孕期均匀分布，n=318/362/335）；**先以 Beckman TPOAb 与 TgAb II 筛查，阳性者（约 10%）排除**；2 个外部站点、1 批试剂、单次测定；非参数法，按 **CLSI C28-A3** 取中央 97.5%（文件原文如此），并给双侧非参数 95% CI |
| K233050 Siemens Centaur TSH3-Ultra II（K251543 沿用「unchanged」） | **甲状腺功能正常成人 0.55–4.78 µIU/mL (n=229)**；婴儿 1–23 月 0.87–6.15 (n=94)；儿童 2–12 岁 0.67–4.16 (n=198)；青少年 0.48–4.17 (n=150) | 「sponsor provided information to support」；纳入排除细节文件未载明；标准列表含 CLSI EP28-A3c |
| K162606 Roche Elecsys TSH | **0.270–4.20 µIU/mL** | 516 名健康受试者的 2.5–97.5 百分位；纳入排除、年龄性别文件未载明；建议各实验室验证可转移性 |
| K130469 DiaSorin LIAISON TSH | **0.357–4.789 mIU/L**（中位 1.438） | 130 名「apparently healthy」受试者血清，观察 95% 正常范围；纳入排除未载明 |

- 逻辑：参考区间由厂家自建（健康人群百分位），单位 µIU/mL = mIU/L；不同厂家上限 4.20–5.33 的差异与人群筛选（是否排除甲状腺自身抗体阳性）、标准品代次（WHO 2nd IRP 80/558 vs 3rd IS 81/565）及抗体特异性有关（申报文件未直接讨论原因，此为背景推断）。
- **溯源**：Beckman K153651 改为 **WHO 3rd International Reference Preparation, NIBSC 81/565**（predicate 为 2nd IRP 80/558），溯源过程「based on EN ISO 17511」，一级工作校准品由 WHO 81/565 制备；Siemens TSH3-Ultra II 溯源 **WHO 3rd IRP 81/565**；Roche Elecsys TSH 溯源 **2nd IRP WHO 80/558**；DiaSorin 溯源 **2nd IRP WHO 80/558**。

### 5. 生物学 / 生理学依据

- 申报文件：TSH 由垂体前叶产生（Indications 原文）；Beckman 参考区间研究把妊娠各孕期单列并给出孕早期下限低至 0.05 µIU/mL（申报数据）；交叉反应研究针对结构相近的糖蛋白激素 hCG、FSH、LH、hGH（K162606、K153651、K233050、K130469）；Siemens 列出荧光素干扰（视网膜荧光血管造影后 48–72 h 内样本 TSH 假性降低）。
- 背景（非申报文件）：TSH 与游离 T4 呈对数-线性负反馈，是原发性甲状腺功能异常最敏感的一线指标；“第三代”指功能灵敏度 ≤0.01–0.02 mIU/L，可区分抑制的 TSH（甲亢）与正常低值；孕早期 hCG 的促甲状腺作用使 TSH 生理性下降，这是 Beckman 分孕期给区间的生理基础。

### 6. 样本类型与样本要求

| 产品 | 样本 | 基质等效性研究 |
|---|---|---|
| K251543 / K233050 Siemens | 血清、K2-EDTA 血浆、Li-heparin 血浆 | K233050：EDTA vs 血清 y=0.99x−0.019 (n=52, r 0.999)；Li-hep y=1.01x−0.034 (n=57, r 0.990)；K251543 沿用 |
| K153651 / K221225 Beckman | 血清（含/不含凝胶）、Li-heparin 血浆 | 79 组配对 + 9 份脱 TSH 后加标（0.01–0.2 µIU/mL）：血清(无胶) vs 血清(胶) y=1.00x−0.0002, r 1.00；血清(胶) vs Li-hep y=1.00x+0.0002, r 0.99 |
| K162606 Roche | 血清、血清分离胶、Li-heparin、K2-EDTA、K3-EDTA（predicate 另含枸橼酸钠、NaF/K-oxalate，本次未声明） | 56 组配对（0.0147–92.5 µIU/mL）PB：Li-hep 斜率 1.001；K2-EDTA 0.994；K3-EDTA 0.983，r 0.999 |
| K130469 DiaSorin | 仅血清（200 µL） | 未做血浆 |
| 稳定性 | 样本稳定性各文件未载明；试剂/校准稳定性：Roche 16 周开瓶；Siemens 63 天在机；Beckman 校准品开瓶 90 天 | — |

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 项目 |
|---|---|---|
| K251543 | CLSI EP05-A3、EP09c (3rd ed.)、EP17-A2、EP06 (2nd ed.) | 精密度、方法比对、LoB/LoD/LoQ、线性 |
| K233050 | CLSI EP05-A3、EP17-A2、EP28-A3c、EP07 (3rd ed.)、EP06 (2nd ed.)、EP25-A | 精密度、检出限、参考区间、干扰/交叉反应、线性、稳定性 |
| K221225 | CLSI EP05-A3、EP06-Ed2、EP25-A、IEC 62304、IEC 62366-1、ISO 15223-1 | 精密度、线性、稳定性、软件、可用性、符号 |
| K153651 | CLSI EP05-A3、EP6-A、EP7-A2、EP09-A3、EP17-A2、EP25-A、EP28-A3c；溯源 EN ISO 17511 | 精密度、线性、干扰、方法比对、检出限、稳定性、参考区间、溯源 |
| K162606 | CLSI EP5-A2 (2005)、EP6-A、EP9-A3、EP17-A2、IEC 61010-2-101；交叉反应引用 EP7-A2 | 精密度、线性、方法比对、检出限、电气安全 |
| K130469 | CLSI EP9-A2、EP5-A2、EP6-A、EP17-A2、EP7-A2 | 方法比对、精密度、线性、检出限、干扰 |

### 8. 临床/准确度验证设计与结果

**功能灵敏度 / 检出限（“第三代”如何验证）**

| 产品 | LoB | LoD | LoQ（定义） | 测量区间 |
|---|---|---|---|---|
| K251543 Atellica TSH3-Ultra II | 0.001 | 0.003 | **0.004 µIU/mL**（within-lab CV ≤20%） | 0.008–150 |
| K233050 Centaur TSH3-Ultra II | 0.005 | 0.008 | **0.010 µIU/mL**（within-lab CV ≤20%） | 0.010–150 |
| K221225 Beckman DxI 9000 | 0.002（3 批 95% 非参数） | 0.003 | **0.003 µIU/mL**（13 份低值血清、3 台 × 3 批 × 5 天 × 9 复；log-log 二次精密度曲线取 **between-run CV 10%**） | 0.01–50.0 |
| K153651 Beckman DxI 800 | 0.0004 | 0.001 | **0.001 µIU/mL**（7 份低值血清各 135 复，between-run CV 10%） | 0.01–50.0 |
| K162606 Roche e 801 | 0.0025（60 复，95 百分位） | 0.005（LoB+1.653 SD） | **0.005 µIU/mL**（5 份低值 0.0078–0.064，**functional sensitivity = within-lab CV 20%**） | 0.005–100 |
| K130469 DiaSorin | 0.014 mIU/L | 0.02 | **0.02 mIU/L**（inter-assay CV <20%，6 份低值 × 72 次） | 0.02–90 |

要点：各厂家均按 CLSI EP17-A2 做 LoB/LoD/LoQ，“功能灵敏度”实质上即 LoQ，但 CV 判据不同（Roche/Siemens/DiaSorin 20%，Beckman 10% between-run）；声明的测量区间下限（0.005–0.02）普遍高于 LoQ。

**精密度（低端水平尤其关键）**
- K251543：血清 0.091 µIU/mL 重复性 2.2%、within-lab 4.2%；0.205 → 4.1%；再现性（3 台 × 3 批 × 5 天）0.087 → 4.43%。
- K233050：血清 0.088 → 重复性 2.5%、within-lab 3.6%。
- K221225：0.022 µIU/mL 总 CV 6.4%（80 复）；再现性 0.024 → 4.0%。
- K153651：0.02 → 总 CV 4.4%；0.37 → 3.5%；38.76 → 5.9%。
- K162606（EP5-A3，n=84）：Serum 1 **0.00851 µIU/mL** 总 CV 11.3%；0.209 → 2.5%；PreciControl TS 0.184 → 2.3%。
- K130469：0.2660 mIU/L 总 CV 5.5%（2 批 × 20 天 n=160）。

**方法比对**
- K251543 vs ADVIA Centaur TSH3-Ultra II（EP09c）：n=323 血清 0.011–147.2，PB 斜率 0.97、截距 −0.006、r 0.998。
- K233050 vs predicate TSH3-Ultra (K083844)：n=404，0.028–147.337，y=0.95x−0.016，r 0.999。
- K221225 vs Access 2：n=111，0.01–47，PB 斜率 1.06 (1.04–1.07)、截距 −0.019、r 1.00。
- K153651 vs HYPERsensitive hTSH（2nd IRP 标化）：n=155，PB 斜率 **0.940** (0.92–0.97)、截距 −0.02、r 0.98（反映 2nd→3rd IS 换标后的系统差）。
- K162606 vs Elecsys 2010：n=130，0.005–85.9，PB y=0.936x−0.003，r 0.999。
- K130469 vs 市售免疫法（EP9-A2）：n=181，0.0257–59.56 mIU/L，PB y=1.005x−0.0030，斜率 CI 0.988–1.026。

**线性 / 高剂量钩状效应 / 交叉反应**
- 线性：K251543 15 级 0.004–158，最大偏离 ±6.5%；K221225 两项研究（全量程 9 级 0.002–58.5 与低端 0.002–4.898），≤0.02 µIU/mL 时偏离 ≤0.000001，>0.02 时 ≤5–6%；K162606 12 级 0–102，斜率 0.952、r 0.9986；K130469 7 级 0.015–91.3，y=0.9807x+0.0013。
- 钩状效应：Roche ≥1,000 µIU/mL 无；Beckman 至 1000 µIU/mL 无（WHO 81/565 加标）；Siemens 至 3000 µIU/mL 无（>150 报「>150」）。
- 交叉反应：Roche LH 10,000 mU/mL <0.038%、FSH 10,000 <0.080%、hGH 1,000 0%、hCG 50,000 0%；Beckman hCG 1,000,000 mIU/mL <0.010%、hFSH 1,000 <0.10%、hLH 3,000 <0.10%；Siemens hCG 200,000、FSH 1,500、LH 600 mIU/mL ≤5%；DiaSorin LH 1000、FSH 5000、hGH 100 ng/mL、hCG 200,000、HAMA 1268 ng/mL 无干扰。
- 生物素：Roche 56.0 ng/mL 无干扰，标签要求高剂量生物素（>5 mg/day）停用 8 h 后采样；Siemens 0.35 mg/dL。
- 异嗜性抗体：Beckman K153651 明示「the study showed interference in some cases」并加限制语句；Siemens 加荧光素警示（>0.24 µg/mL 使结果降低）。

### 9. 厂家间差异与要点

1. **标准品代次**：Beckman 2016 年以后与 Siemens TSH3-Ultra II 溯源 WHO 3rd IS 81/565；Roche Elecsys TSH 与 DiaSorin 仍溯源 2nd IRP 80/558。Beckman 换标后与旧方法斜率 0.94，提示不同代标准品间约 6% 的系统差（申报数据）。
2. **参考区间人群定义严格程度不同**：Beckman 是本组唯一在决策摘要中写明**排除 TPOAb/TgAb 阳性**并分孕期建立区间的产品；Siemens 提供分年龄段儿科区间；Roche 仅 n=516 健康人 2.5–97.5 百分位；DiaSorin n=130。
3. **功能灵敏度判据**：Beckman 用 10% between-run CV（更严格），其余用 20%。测量区间下限多设在 LoQ 的 2–10 倍。
4. **样本**：DiaSorin 仅血清；Roche 声明最广（含 K2/K3-EDTA）；Beckman 未声明 EDTA。
5. 声明“第三代”并不是 FDA 定义的法规术语，决策摘要中只以 LoQ/功能灵敏度数据支持。

### 10. 来源

- K251543: https://www.accessdata.fda.gov/cdrh_docs/reviews/K251543.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf25/K251543.pdf
- K233050: https://www.accessdata.fda.gov/cdrh_docs/reviews/K233050.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf23/K233050.pdf
- K221225: https://www.accessdata.fda.gov/cdrh_docs/reviews/K221225.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K221225.pdf
- K153651: https://www.accessdata.fda.gov/cdrh_docs/reviews/K153651.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153651.pdf
- K162606: https://www.accessdata.fda.gov/cdrh_docs/reviews/K162606.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K162606.pdf
- K130469: https://www.accessdata.fda.gov/cdrh_docs/reviews/K130469.pdf （Triage）；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K130469.pdf
- K190773（未采用）: https://www.accessdata.fda.gov/cdrh_docs/reviews/K190773.pdf ；K083173: https://www.accessdata.fda.gov/cdrh_docs/reviews/K083173.pdf

---

## 25-羟维生素 D（25-OH Vitamin D, 25(OH)D）

### 1. 法规定位

| 路径 | Product code | 21 CFR | Class | Panel | openFDA 累计 |
|---|---|---|---|---|---|
| 免疫/结合法 | **MRG** – Vitamin D test system | **862.1825** | II | CH (75) | **46** |
| LC-MS/MS | **PSL** – Vitamin D mass spectrometry test system | **862.1840**（Class II, special controls，源自 De Novo DEN170019） | II | CH (75) | **1**（DEN170019 AB Sciex Vitamin D 200M Assay for the Topaz System，2017-05-18） |

- 862.1840 特殊控制（DEN170019 决策摘要末页可读部分；条文前半段在文本抽取中缺失）：(2)(i) 精密度须用预期样本类型在**医学决策点**浓度，至少一份未修饰患者样本，按 FDA 认可标准评价重复性与再现性；(ii) 准确度 ≥115 份血清/血浆，与参考方法或合法上市的标准化质谱法比对；(iii) 须在标签描述维生素 D2/D3、1-OH-D2/D3、3-epi-25OH-D2/D3、1,25(OH)2-D2/D3、3-epi-1,25(OH)2-D2/D3、25,26-(OH)2-D3、24(R),25-(OH)2-D3、23(R),25-(OH)2-D3 的干扰；(3) 参考区间研究须用美国≥3 个气候区、不同季节、≥21 岁、人种性别代表美国人口的健康成人；(4) 结果**只能报告总 25-OH 维生素 D**。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K210901 | 2021-09-01 | Roche / Elecsys Vitamin D total III (cobas e 601) | 竞争性 ECLIA，钌标记 VDBP，DTT/NaOH 预处理 | 血清、Li-heparin、K2/K3-EDTA 血浆、SST | aid in the assessment of vitamin D sufficiency **in adults**；Rx | 决策摘要 + summary |
| K200509 | 2020-05-29 | Siemens / ADVIA Centaur Vitamin D Total (VitD) | 竞争性化学发光（吖啶酯抗体，FITC-维生素 D 类似物） | 血清、SST、Na-/Li-heparin、K2/K3-EDTA | aid in the determination of vitamin D sufficiency；本次为恢复血浆声明 | 决策摘要 + summary |
| K142373 | 2014-12-22 | Beckman Coulter / Access 25(OH) Vitamin D Total (Access 2) | 两步竞争性免疫酶法（羊单抗，ALP 结合物） | 血清（含/不含胶）、Li-heparin 血浆 | aid in the assessment of vitamin D sufficiency | 决策摘要 + summary |
| K223503 | 2023-01-19 | Beckman Coulter / Access 25(OH) Vitamin D Total (DxI 9000) | 同上 | 同上 | 同上；声明「Traceable to NIST-Ghent ID-LC-MS/MS」 | Triage 决策摘要 + 510(k) summary |
| K153375 | 2016-08-12 | Abbott / ARCHITECT 25-OH Vitamin D 5P02 (+Calibrators, Controls) | 延迟一步竞争 CMIA（兔单抗） | 血清、SST、K2/K3-EDTA、Na-/Li-heparin、PST | aid in the assessment of vitamin D sufficiency；标化 NIST SRM 2972 | **仅 510(k) summary**（无决策摘要） |
| K132492 | 2013-09-05 | DiaSorin / LIAISON 25 TOTAL-D (+Control Set, Calibration Verifiers) | 直接竞争 CLIA（羊多抗，异鲁米诺） | 血清、EDTA 及 Li-heparin 血浆 | 「25-hydroxyvitamin D and other hydroxylated vitamin D metabolites」；aid in the assessment … in adults | Triage 决策摘要 + 510(k) summary |
| K121608 | 2013-02-15 | Ortho / VITROS 25-OH Vitamin D Total Reagent Pack + Calibrators | 直接竞争化学发光（HRP 结合物，羊单抗） | **仅血清** | assessment of Vitamin D sufficiency；成人；儿科性能未建立 | 决策摘要 + summary |
| DEN170019 | 2017-05-18 | AB Sciex / Vitamin D 200M Assay for the Topaz System | **LC-MS/MS**（蛋白沉淀，同位素内标，分别定量 25-OH-D3 与 D2 求和） | 仅血清（红头管） | assessment of vitamin D sufficiency in adults；trained lab professional | De Novo 决策摘要 |

### 3. 预期用途与声明类型

- 统一为「aid in the assessment (determination) of vitamin D sufficiency」——既非“诊断维生素 D 缺乏症”，也非监测；Roche、DiaSorin、Ortho、AB Sciex 限定**成人**（Ortho：「performance characteristics … have not been established in a pediatric population」；Siemens 反而在标签中给出儿科 12 月–21 岁参考区间）。
- DiaSorin 措辞「25-hydroxyvitamin D **and other hydroxylated vitamin D metabolites**」反映其抗体对 1,25(OH)2D 等有交叉反应（见第 8 节）。
- 全部 Rx、实验室用；本组无 POC/waived（列表中 Nanoentek FREND K162754、Immunostics ALFIS K221817 为小型仪器，未读取）。
- LC-MS/MS 产品要求「by a trained laboratory professional in a clinical laboratory」。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

- 所有文件「Assay cut-off: Not applicable」「Clinical cut-off: Not applicable」（K142373 补充「this is a quantitative assay」）。
- 申报方给出的是**健康人群参考区间**（多为 2.5–97.5 百分位），而 20/30 ng/mL 的“缺乏/不足/充足”分层**仅以引用文献形式**出现在标签中（不是申报的性能声明）：

| 产品 | 参考区间（申报数据） | 人群设计 | 标签中引用的分层 |
|---|---|---|---|
| K210901 Roche | 全部 n=463：2.5–97.5 百分位 **10.2–49.4 ng/mL**（中位 25.7）；夏 n=245：12.5–52.4；冬 n=218：9.4–44.1 | 美国南/中/北 3 地，夏冬两季，男女约各半，约 30% 深肤色，22–79 岁；EP28-A3c | 文件未载明分层引用 |
| K200509 Siemens | 成人中位 22.5，2.5–97.5 百分位 **7.4–44.0 ng/mL**；儿科（12 月–21 岁）中位 23.8，11.4–45.8 | 本次仅 20 份健康人样本按 EP28-A3c **验证**原区间可转移 | 标签引用 Holick 2007 等：缺乏 <20；不足 20–<30；充足 30–100 ng/mL（儿科 <15 / 15–<20 / 20–100） |
| K142373 Beckman | n=367，中位 24.9，**11.9–43.6 ng/mL** | 21–89 岁健康成人，美国多地域，性别人种代表美国人口，冷暖季节，20% 服补充剂；纳入要求 Ca、Mg、P、PTH、TSH 正常 | 文件未载明 |
| K153375 Abbott | 合并 n=283：**6.6–49.9 ng/mL**（均值 ±1.96 SD）；冬 n=129 6.2–45.5；夏 n=154 7.0–53.2 | ≥21 岁，男 142 女 141，北/南/中 3 地，≥30% 深肤色与 ≥30% 浅肤色，非裔/西裔/白人；排除 Ca（2.15–2.50/2.55 mmol/L）、TSH（0.35–4.94）、iPTH（15.0–68.3 pg/mL）异常者；≤50% 服补充剂；C28-A3c | 控制品靶值 20 与 40 ng/mL 「bracket the lower and upper medical decision points of 20.0 and 30.0 ng/mL」 |
| K132492 DiaSorin | n=395，中位 22.9，中央 95% **6.8–54.2 ng/mL** | 4 个美国站点（北/南/中），21–90 岁，夏冬，深浅肤色；健康定义：总钙、iPTH、TSH 正常；无肾/胃肠/肝/甲状腺/甲状旁腺病史、无癫痫、无减重手术、非孕/哺乳；≥50% 不服补充剂，服者 <2000 IU/day；C28-A2 | 文件未载明 |
| K121608 Ortho | n=399，中位 33.4，2.5–97.5 百分位 **14.7–68.3 ng/mL**（三系统均值） | 21–79 岁（招募 21–90），男女各半，北/南/中，夏冬各半，≥30% 深肤 ≥30% 浅肤；排除甲状旁腺/钙调节病、肾/胃肠/肝病、减重手术；不服 >1000 IU/day | 标签引用 Endocrine Society 2011（Holick）：缺乏 <20；不足 20–29；充足 30–100；潜在毒性 >100 ng/mL |
| DEN170019 AB Sciex | n=404，内 95% **8.6–49 ng/mL** | ≥21 岁健康男女，3 个地理区；EP28-A3 | 文件未载明 |

- 逻辑：因维生素 D 状态受季节/纬度/肤色/补充剂影响，FDA（自 2013 前后）要求参考区间研究覆盖美国≥3 个气候区、两季、不同肤色，并在 862.1840 中固化为特殊控制；免疫法产品虽无特殊控制，但 Roche/Abbott/DiaSorin/Ortho/Beckman 的设计高度一致。
- **溯源**：Roche → ID-LC-MS/MS 参考测量程序 → NIST SRM 2972；Siemens → NIST SRM 972a；Beckman → JCTLM 认可的 Ghent 大学 ID-LC-MS/MS RMP → NIST SRM 2972，并按 **VDSP** 用 40 份 CDC/VDSCP 样本验证（Deming y=1.01x−2.87, r 0.99）；Abbott → NIST SRM 2972；Ortho → 内部参考校准品由溯源 LC-MS 值的患者样本定值；DiaSorin（2013）→ 内部标准品 UV 分光定值（无 NIST/VDSP 溯源）；AB Sciex → **CDC VDSCP 认证**。

### 5. 生物学 / 生理学依据

- 申报文件：25-OH-D 与维生素 D 结合蛋白（DBP）结合，免疫法需预处理释放（Roche DTT+NaOH；Siemens 释放剂；Beckman DBP releasing agent；Abbott ANSA/甲醇）；Abbott 文件引用 Heijboer 2012 说明**血液透析患者**因 DBP 浓度差异使自动免疫法相对 LC-MS/MS 出现负偏倚（Abbott 实测 −15.3%），孕妇三孕期无显著偏倚（4.5/−2.2/0.1%）；Beckman 引用 Carter 2007 说明外源加标 25-OH-D 在竞争性免疫法中回收不足，故交叉反应需按 25-OH-D3 归一化。
- 背景（非申报文件）：25-OH-D 是循环中半衰期最长（2–3 周）的维生素 D 代谢物，是营养状态指标；D2 来自植物/处方补充剂，D3 来自皮肤合成和动物食物，故“总 25-OH-D”要求 D2/D3 等摩尔识别；C3-epimer 在婴儿中比例高，在成人通常 <5 ng/mL（K142373 引 Keevil 2012 给出成人最高 4.9 ng/mL 的申报引用）。

### 6. 样本类型与样本要求

| 产品 | 样本 | 基质等效性 | 稳定性 |
|---|---|---|---|
| K210901 Roche | 血清、SST、Li-hep、K2-EDTA、K3-EDTA | 46 组配对 PB：Li-hep 斜率 1.03；K2-EDTA 0.973；K3-EDTA 0.92；SST 1.02 | 文件未载明 |
| K200509 Siemens | 血清、SST、Na-hep、Li-hep、K2-、K3-EDTA | 66 天然 + 8 人工配对（13.0–142.9），加权 Deming vs 普通血清：SST 0.97；Na-hep 1.02；Li-hep 0.99；K2 0.97；K3 0.96 | 文件未载明 |
| K142373 Beckman | 血清（胶/无胶）、Li-hep | 45 组：血清(胶) 斜率 1.01；Li-hep 1.05 (1.02–1.10) | 校准品开瓶 56 天 |
| K153375 Abbott | 血清、SST、K2/K3-EDTA、Na-hep、Li-hep 粉、PST | 51 组 PB：SST 0.99x+0.04；K2-EDTA 0.93x+0.76；K3 0.92x+0.83；Na-hep 0.94x+0.84；Li-hep 0.94x+0.90；PST 0.93x+1.09 | 2–8 °C ≥12 天（含/不含血凝块）；22–30 °C ≥72 h；≥4 次冻融 |
| K132492 DiaSorin | 血清、SST、EDTA、Li-hep | 64 组（4.0–139.3）：SST 斜率 0.9942；EDTA 1.0092；Li-hep 1.0007 | 文件未载明 |
| K121608 Ortho | 仅血清 | Not applicable | 文件未载明 |
| DEN170019 AB Sciex | 仅血清（红头促凝管） | Not applicable | 文件未载明 |

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 项目 |
|---|---|---|
| K210901 | CLSI EP05-A3、EP06-A、EP17-A2、EP28-A3c | 精密度、线性、检出限、参考区间 |
| K200509 | CLSI EP17-A2、EP06-A、EP05-A3、**EP34**（扩展测量区间/稀释）、EP28-A3c、EP25-A | 检出限、线性、精密度、稀释回收、参考区间验证、稳定性 |
| K142373 | CLSI EP5-A2、EP7-A2、EP17-A2、EP9-A2、EP25-A、EP6-A；溯源 JCTLM/VDSP | 精密度、干扰、检出限、方法比对、稳定性、线性 |
| K223503 | CLSI EP05-A3、EP06-2nd、EP17-A2、EP09c-A3 | 精密度、线性、检出限、方法比对 |
| K153375 | EP5-A2、EP6-A、EP17-A2、EP7-A2、EP09-A3、C28-A3c（summary 正文引用） | 精密度、线性、检出限、干扰/交叉、方法比对与病理人群、参考区间 |
| K132492 | CLSI EP9-A2、C28-A2、EP5-A2、EP6-A、EP17-A2、EP7-A2 | 方法比对、参考区间、精密度、线性、检出限、干扰 |
| K121608 | EN 13640、C28-A2、EP7-A、EP9-A2、EP10-A2、EP5-A2；正文引用 EP6-A、EP17 | 稳定性、参考区间、干扰、方法比对、初步评价、精密度、线性、检出限 |
| DEN170019 | CLSI EP25-A、EP28-A3c、EP17-A2、EP9-A3、EP07-A2、EP06-A、EP05-A3、**C62-A**（LC-MS 方法）、IEC 61010-1/-2-101/-2-061 | 稳定性、参考区间、检出限、方法比对、干扰、线性、精密度、LC-MS 方法学、电气安全 |

### 8. 临床/准确度验证设计与结果

**准确度（与参考方法比对）**

| 产品 | 比对方法 | n / 范围 | 回归 | 决策点偏倚 |
|---|---|---|---|---|
| K210901 Roche | **CDC 维生素 D 参考实验室 ID-LC-MS/MS 定值的 VDSCP 验证样本** | 157 份天然单供者血清，5.6–118.4 | Deming y=0.981x+0.795, r 0.982 | **30 ng/mL 处预测偏倚 0.8%** |
| K142373 Beckman | Ghent 大学 ID-LC-MS/MS RMP（EP9-A3） | 109 份，8.0–109.4 | PB 斜率 1.02 (0.95–1.12)、截距 −3.57 (−6.15~−1.30)、r 0.95；另 VDSP 40 份 Deming y=1.01x−2.87 | 文件未载明 |
| K153375 Abbott | 外部参考实验室 ID-LC-MS/MS（EP09-A3） | ≥100 份，4.0–153.2（≤10% 加标） | PB 斜率 1.02、截距 −0.99、r 0.99 | 文件未载明（病理人群：孕妇 +4.5/−2.2/+0.1%，透析 −15.3%） |
| DEN170019 AB Sciex | CDC VDSCP 定值样本 | 118 份，5.6–133（12 份人工） | PB 斜率 1.008、截距 −0.3949、r 0.991 | 文件未载明 |
| K200509 Siemens | predicate（自身旧配方） | 126 份（118 天然 + 8 加标），5.93–130.85 | Deming 两批：斜率 1.03/1.04，截距 0.85/1.74，r 0.99 | — |
| K223503 Beckman DxI 9000 | predicate Access 2 | 150 份，7.0–120 | 斜率 1.05 (0.99–1.10)、截距 0.94、r 0.97（接受 R²≥0.90，斜率 1.00±0.10） | — |
| K132492 DiaSorin | predicate LIAISON 25 OH Vitamin D TOTAL | 391 份（403 中 12 份 <4.0 排除），4.18–135 | PB 斜率 0.99 (0.97–1.01)、截距 −0.22、r 0.990 | — |
| K121608 Ortho | predicate IDS-iSYS | n=102–103/系统，12.8–126 | 5600: 0.99x−5.12, r 0.92；3600: 1.08x−7.87, r 0.93；ECi: 0.96x−9.07, r 0.94 | — |

**精密度（决策点附近）**
- Roche：21 天 n=84，12.3 ng/mL within-lab 9.8%；28.7 → 6.0%；33.0 → 5.6%；52.3 → 4.0%。
- Siemens：20 天 n=80，21.29 → within-lab 9.6%；32.16 → 7.4%；65.47 → 5.5%。
- Beckman K142373：3 台 × 3 批 × 20 天，24.6 → 总 7.5%；49.8 → 7.3%；K223503 DxI 9000：7.0 → 18.6%（SD 1.3）；28 → 6.8%；71 → 4.3%。
- Abbott：2 台 × 3 批 × 20 天 n≈358，21.1 → 3.2%；30.5 → 3.1%；5.3 → 6.9%（SD 0.37）；接受标准 <8 ng/mL 时 SD ≤0.8，≥8 时 CV ≤10%。
- DiaSorin：1 批 20 天 n=80，17.7 → 总 10.3%；28.2 → 8.4%；34.3 → 8.6%。
- Ortho：多系统多批，22.5 → within-lab 15.3%；31.1 → 13.3%；≈120 → 5.5%。
- AB Sciex LC-MS/MS：3 站点 × 5 天 × 5 复，31.0 → 再现性 6.5%；未修饰患者样本 28.4 → 7.4%。

**LoB / LoD / LoQ 与测量区间**
- Roche 2.0 / 3.0 / 6.0（20% CV），6–120 ng/mL（可 1:2 稀释至 240）；Siemens 1.7 / 3.2 / 4.2，4.2–150；Beckman Access 2 0.55 / 1.0 / 3.0，7.0–120（DxI 9000：2.5/4.5/7.0）；Abbott 1.6 / 2.2 / 2.4，3.4–155.9；DiaSorin <2.03 / 2.91 / 4.0，4–150；Ortho 4.34 / 8.64 / 12.8，12.8–126；AB Sciex LLMI 2.9（<20% bias 且 <20% CV），4–140。

**D2/D3 交叉反应与 C3-epimer（关键差异）**

| 产品 | 25-OH-D2（相对 D3） | 3-epi-25-OH-D3 | 1,25(OH)2-D | 24,25(OH)2-D3 | 其他 |
|---|---|---|---|---|---|
| Roche K210901 | 103.3%（归一化） | **121.6%**（归一化；3-epi-D2 102.7%） | D3 未检出；D2 0.9% | 8.1% | 生物素 >600 ng/mL 才干扰（≤10%） |
| Beckman K142373 | 96–116%（归一化） | **54–70%**（归一化） | D2 1043–2262%、D3 212–324%（归一化；测试浓度为内源 125–375 倍） | 2–9% | 帕立骨化醇 264–288%，标签要求用药 24 h 内不测；Hb >50 mg/dL 假性升高 |
| Abbott K153375 | 80.5%、82.4%（内源 D2 样本 vs LC-MS/MS） | **1.3%**（100 ng/mL） | 0.1% / −0.4% | **101.9–189.2%**（D3）、71.4–114.2%（D2） | 帕立骨化醇 0.6%；TG >500 mg/dL 负干扰（800 mg/dL 时 −10~−17.5%） |
| DiaSorin K132492 | 93% | 1.9% | D3 17.1%、D2 27.1% | 文件未载明 | 维生素 D3 3.6%、D2 1.9% |
| Ortho K121608 | 104.9% | **37.4%** | D2 >100%（0.2 ng/mL 时偏倚 1.7 ng/mL）、D3 −5.0% | 34.3–34.8% | 帕立骨化醇有正干扰（标签限制） |
| AB Sciex DEN170019 | 分别定量 D2/D3 求和 | 3-epi-D2/D3 100 ng/mL 无干扰（<10%） | 10,000 pg/mL 无干扰 | 150 ng/mL 无干扰 | 79 种物质均 <10% |

### 9. 厂家间差异与要点

1. **溯源与准确度基础**：2014 年后的主流免疫法（Roche、Beckman、Abbott、Siemens）均声明 NIST SRM 972/972a/2972 或 Ghent ID-LC-MS/MS RMP 溯源，并以 ID-LC-MS/MS（CDC/VDSCP 或 Ghent 定值样本）作为方法比对的比较方法；2013 年的 DiaSorin/Ortho 仍以 predicate 免疫法比对、内部 UV 定值。LC-MS/MS 产品（PSL）则要求 CDC VDSCP 认证与 ≥115 份样本比对。
2. **C3-epimer**：Roche 完全交叉（121.6%）、Beckman 54–70%、Ortho 37%、Abbott 1.3%、DiaSorin 1.9%；LC-MS/MS 通过色谱分离而无干扰。这直接影响婴幼儿样本的适用性（各文件成人限定的部分原因，申报文件未明说）。
3. **1,25(OH)2-D 与帕立骨化醇**：Beckman/DiaSorin/Ortho 交叉显著（Beckman、Ortho 标签限制帕立骨化醇），Roche/Abbott 近零。
4. **24,25(OH)2-D3**：Abbott 100–189% 完全交叉，Roche 8%，Beckman 2–9%。
5. **参考区间**均按“美国 3 气候区 + 两季 + 肤色配额 + 钙/PTH/TSH 正常 + 补充剂限制”设计；结果差异（下限 6.6–14.7，上限 43.6–68.3 ng/mL）反映方法偏倚与季节构成。
6. **“充足/不足”分层**只以文献引用出现（Siemens、Ortho），FDA 决策摘要均写 clinical cut-off 不适用。
7. **特殊人群**：Abbott 是本组唯一按 EP09-A3 做孕妇与透析患者偏倚研究的产品，并明确透析人群 −15.3% 偏倚。

### 10. 来源

- K210901: https://www.accessdata.fda.gov/cdrh_docs/reviews/K210901.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf21/K210901.pdf
- K200509: https://www.accessdata.fda.gov/cdrh_docs/reviews/K200509.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf20/K200509.pdf
- K142373: https://www.accessdata.fda.gov/cdrh_docs/reviews/K142373.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K142373.pdf
- K223503: https://www.accessdata.fda.gov/cdrh_docs/reviews/K223503.pdf （Triage）；https://www.accessdata.fda.gov/cdrh_docs/pdf22/K223503.pdf
- K153375: https://www.accessdata.fda.gov/cdrh_docs/pdf15/K153375.pdf （决策摘要 URL 无文件）
- K132492: https://www.accessdata.fda.gov/cdrh_docs/reviews/K132492.pdf （Triage）；https://www.accessdata.fda.gov/cdrh_docs/pdf13/K132492.pdf
- K121608: https://www.accessdata.fda.gov/cdrh_docs/reviews/K121608.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf12/K121608.pdf
- DEN170019: https://www.accessdata.fda.gov/cdrh_docs/reviews/DEN170019.pdf

---

## 铁蛋白（Ferritin）

### 1. 法规定位

- Product code **DBF**（Ferritin, antigen, antiserum, control），**21 CFR 866.5340** Ferritin immunological test system，Class II，Panel **IM – Immunology (82)**。
- openFDA 累计清关 **total = 87**。
- 伴随代码：JIX（多分析物校准品，862.1150）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K193650 | 2021-09-14 | DiaSorin / LIAISON Ferritin (LIAISON XL) | 夹心 CLIA（鼠单抗磁微粒 + 异鲁米诺单抗） | 血清、SST、Li-heparin 血浆 | aid in the diagnosis of iron deficiency anemia and iron overload；Rx | 决策摘要 + summary |
| K171642 | 2017-08-31 | Siemens / Atellica IM Ferritin (Fer) | 两位点夹心化学发光（羊多抗吖啶酯 + 鼠单抗磁微粒） | 血清、EDTA 及 Li-heparin 血浆 | aid in the diagnosis of iron deficiency anemia and iron overload | 决策摘要 + summary |
| K100538 | 2010-06-22 | Roche / Tina-quant Ferritin Gen.4 (Roche/Hitachi 902/912/917/Modular P) | 乳胶增强免疫比浊（兔抗） | 血清、Li-heparin、EDTA 血浆 | aid of diagnosis of diseases affecting iron metabolism | 决策摘要 + summary |
| K110736 | 2011-08-17 | Siemens / ADVIA Chemistry Ferritin (FRT) + Liquid Specific Protein Calibrator | 乳胶免疫比浊（658 nm） | 血清、肝素血浆、EDTA 血浆 | aid in the diagnosis of diseases affecting iron metabolism, such as hemochromatosis (iron overload) and iron deficiency anemia | **仅 510(k) summary** |
| K191562 | 2020-03-06 | Horiba ABX / Yumizen C1200 Ferritin（同批含 Transferrin、RF） | 乳胶增强免疫比浊（兔抗） | **仅血清** | 同上 | 决策摘要 + summary |

### 3. 预期用途与声明类型

- 两类措辞：(a) 「aid in the diagnosis of iron deficiency anemia and iron overload」（DiaSorin、Siemens Atellica，源自 ADVIA Centaur predicate）；(b) 「aid in the diagnosis of diseases affecting iron metabolism, such as hemochromatosis (iron overload) and iron deficiency anemia」（Roche、Siemens ADVIA Chemistry、Horiba，Beckman predicate 另写「Serum ferritin is an indicator of body iron stores: it has been shown to correlate with stainable bone marrow iron」）。
- 均为**辅助诊断**、Rx、实验室用；无 POC/waived。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

- 所有文件「Assay cut-off: Not applicable」；以**性别分层参考区间**为判定依据。

| 产品 | 参考区间（申报文件） | 建立/验证方式 |
|---|---|---|
| K171642 Siemens Atellica | **女 (n=275, 16–94 岁) 中位 46.4，95 百分位区间 7.3–270.7 ng/mL；男 (n=179, 15–95 岁) 中位 55.9，10.5–307.3 ng/mL** | 自建，CLSI EP28-A3c；健康人，肝功能酶、胆红素、血清铁正常 |
| K193650 DiaSorin | 女 (n=39, 18–92 岁) 中位 52.2，6.3–317.1；男 (n=39, 21–82 岁) 中位 124.3，3.6–362.4 ng/mL | 「Previously established ferritin reference ranges were verified」——78 份血清（美国 5 地，按 census.gov 人种分布选样）**验证**；原区间来源文件未载明 |
| K100538 Roche | 男 (20–60 岁) 30–400 ng/mL；女 (17–60 岁) 15–150 ng/mL | 沿用 predicate；来源文件未载明 |
| K110736 Siemens ADVIA Chem | 男 20–250 ng/mL；女 10–120 ng/mL（predicate N Latex：男 n=216 20–290；绝经前女 n=193 4.5–170；绝经后女 n=47 24–260 µg/L） | 文件未载明来源与 n |
| K191562 Horiba | 女 10–120、男 20–250 ng/mL | 引用 Tietz 第 4 版（Roberts 等 2006），用血库样本女 50、男 95 份**验证** |

- **溯源（WHO 标准品）**：DiaSorin → **WHO NIBSC 94/572**（3rd IS，重组）；Horiba → 3rd IS 94/572；Siemens ADVIA Chemistry 校准品 → **WHO 3rd IS 94/572**；Siemens Atellica IM → **WHO 2nd IS 80/578**（与 ADVIA Centaur 相同）；Roche Tina-quant Gen.4 → 对 Elecsys Ferritin 标化，Elecsys 溯源 **WHO IS 94/572、80/578、80/602**；Elecsys predicate（K971833）在 DiaSorin 对比表中标为 80/602。
- 逻辑：无临床 cutoff；缺铁/铁过载判断由临床按参考区间与文献阈值（背景：WHO 缺铁 <15 µg/L 成人；炎症状态另议——此为背景非申报）进行。

### 5. 生物学 / 生理学依据

- 申报文件：Roche 抗体「specific for ferritin from human liver and recognize ferritin from human spleen … no cross reactivity to the human ferritin H subunit」（识别 L 亚基为主）；Beckman predicate 表述血清铁蛋白与骨髓可染铁相关；Siemens Atellica 警示「grossly hemolyzed samples … release of intracellular ferritin can cause elevated results」。
- 背景（非申报文件）：血清铁蛋白与体内储存铁成比例，是缺铁最早的实验室指标，同时为急性期反应蛋白，炎症/肝病/恶性肿瘤时升高，故预期用途限于“辅助”。

### 6. 样本类型与样本要求

| 产品 | 样本 | 基质等效性 | 稳定性 |
|---|---|---|---|
| K193650 DiaSorin | 血清、SST、Li-hep | 43 组（37 天然 + 6 人工，3.8–1946）PB：SST 斜率 1.00；Li-hep 0.98 (0.94–0.99)，截距 −1.73 | 20–25 °C 24 h；2–8 °C 7 天；−20 °C 1 个月 |
| K171642 Siemens Atellica | 血清、EDTA、Li-hep | 56 组（7 人工，2.5–1440.7）加权 Deming：Li-hep y=0.95x+1.6；K-EDTA y=0.96x+0.1，r 0.997–0.998 | 文件未载明 |
| K100538 Roche | 血清、Li-hep、K2-EDTA、K3-EDTA（predicate 含枸橼酸，本次删除） | n=94–95 PB：Li-hep 1.020x−1.259；K2-EDTA 0.999x−0.877；K3-EDTA 0.985x−1.800 | 文件未载明 |
| K110736 Siemens ADVIA Chem | 血清、肝素、EDTA 血浆 | 文件未载明 | 文件未载明 |
| K191562 Horiba | 仅血清 | 未做 | 20–25 °C 7 天；2–8 °C 7 天；−20 °C 1 年（引用文献） |

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 项目 |
|---|---|---|
| K193650 | CLSI EP5-A3、EP15-A3（用户验证精密度/偏倚）、EP07-A、EP06-A、EP17-A2、EP28-A3c；正文另引用 EP07-A2、**EP37-A**（干扰补充表） | 精密度、线性、干扰、检出限、参考区间验证 |
| K171642 | CLSI EP05-A3、EP06-A、EP17-A2、EP28-A3c、EP7-A2、EP09-A3 | 精密度、线性、检出限、参考区间、干扰、方法比对 |
| K100538 | CLSI EP17-A、EP5-A2 | 检出限、精密度 |
| K191562 | CLSI EP05-A3、EP17-A2、EP06-A、C28-A3、EP25-A | 精密度、检出限、线性、参考区间、稳定性 |
| K110736 | 文件未载明 | — |
| 溯源 | WHO NIBSC 94/572、80/578、80/602（各文件） | 溯源 |

### 8. 临床/准确度验证设计与结果

- **K193650 DiaSorin**：精密度 2 台 × 2 批 × 20 天 n=320：5.8 ng/mL 总 CV 5.6%；18.3 → 4.6%；178 → 4.1%；1883 → 6.4%。再现性（5 天 × 6 复）5.6 → 9.4%。线性 4 份样本各 7 级，斜率 0.99–1.02，回收 90–110%；扩展区间 1:50 自动稀释至 100,000 ng/mL，回收 99–108%。LoB 0.004、LoD 0.073、LoQ 0.461（20% CV）；测量区间 0.46–2200。方法比对 vs Elecsys Ferritin：n=173（2.0–1932）PB 斜率 0.96 (0.95–0.98)、截距 −1.12、R² 0.995。干扰（±10%，18 复）：白蛋白 60 mg/mL、TG 30 mg/mL、Hb 10 mg/mL、胆红素 0.2 mg/mL、HAMA 54.5 µg/mL、RF 33.1 IU/mL 及 16 种药物无干扰。
- **K171642 Siemens Atellica IM**：精密度 1 台 × 20 天 n=80：8.2 → within-lab 5.4%；41.9 → 4.2%；118.3 → 4.0%；1453.6 → 6.3%。线性 0.9–1650，加权 Deming y=0.914x−0.888, R² 0.999；钩状效应至 80,000 ng/mL 无。LoB 0.3、LoD 0.7、LoQ 0.9（≤20% CV）。方法比对 vs ADVIA Centaur Ferritin：126 份中 107 份在范围内（3.4–1641.4），r 0.99，加权 Deming y=1.03x−0.5（斜率 CI 1.01–1.04）。干扰（≤10%）：Hb 900 mg/dL、胆红素 60、Intralipid 2000 mg/dL、抗坏血酸 176 mg/dL 等。
- **K100538 Roche Gen.4**：Hitachi 917 精密度 n=84：8.48 ng/mL 重复性 7.2%、中间精密度 9.9%；235 → 0.9/1.8%。线性 5–1000（917），PB y=1.0098x−0.0646；902 声明 5–800。LoB 3、LoD 5、LoQ 7（between-run ≤20%）。方法比对 vs predicate：917 n=94（15.0–775.3）PB y=0.987x+0.040；902 n=84 y=0.979x−1.188。干扰接受：≤40 ng/mL 时 ±4 ng/mL，>40 时 ±10%；胆红素 60 mg/dL、Hb 500 mg/dL、Intralipid 1000 mg/dL、RF 1200 IU/mL 无干扰；钩状效应至 80,000 无。
- **K110736 Siemens ADVIA Chemistry**：方法比对 vs N Latex Ferritin n=47，斜率 1.00 (0.97–1.03)、截距 0.00 (−3.4~3.4)；测量区间 6–(450–500) ng/mL；钩状至 40,000。精密度、线性、干扰数据 summary 未载明。
- **K191562 Horiba**：精密度 3 台 × 20 天 n=240：29.6 → 总 7.9%；172.6 → 1.9%；3 批再现性 19.1 → 11.8%。线性 13.3–426.6 y=0.99x+3.16；LoB 4.21、LoD 6.30、LoQ 9.39（TE ≤20%）；测量区间 10–450（1:5 重测至 2250）；前带效应 ≥5043 报警。方法比对 vs Beckman OSR61203：n=103（12.2–440.7）PB 斜率 **0.91** (0.91–0.92)、截距 2.10、r 0.999。

### 9. 厂家间差异与要点

1. **WHO 标准品代次不统一**：94/572（3rd IS，重组；DiaSorin、Horiba、Siemens ADVIA Chem）vs 80/578（2nd IS；Siemens Atellica/Centaur）vs Roche 多重溯源（94/572、80/578、80/602）。不同标准品对齐可能是各厂家参考区间上限差异（150/250/270/307/362 ng/mL）的原因之一（申报文件未讨论）。
2. **参考区间建立方式**：Siemens Atellica 自建（n=454，按 EP28-A3c，附肝功能/胆红素/血清铁正常筛选）；DiaSorin、Horiba 仅“验证”既有区间（n=78 / 145）；Roche、ADVIA Chem 沿用文献/predicate。仅 Siemens N Latex predicate 有绝经前/后分层。
3. **测量区间**：化学发光法宽（0.46–2200、0.9–1650）；比浊法窄（5–1000、10–450、6–500），需稀释重测；前带/钩状效应阈值差异大（Horiba 5043 vs Roche/Siemens 80,000）。
4. **样本**：Horiba 仅血清；其余含 Li-heparin，Siemens/Roche 另含 EDTA；Roche Gen.4 删除枸橼酸血浆。
5. 干扰接受标准：DiaSorin/Siemens/Horiba ±10%；Roche 低浓度用绝对值 ±4 ng/mL。

### 10. 来源

- K193650: https://www.accessdata.fda.gov/cdrh_docs/reviews/K193650.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K193650.pdf
- K171642: https://www.accessdata.fda.gov/cdrh_docs/reviews/K171642.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K171642.pdf
- K100538: https://www.accessdata.fda.gov/cdrh_docs/reviews/K100538.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf10/K100538.pdf
- K110736: https://www.accessdata.fda.gov/cdrh_docs/pdf11/K110736.pdf （决策摘要 URL 无文件）
- K191562: https://www.accessdata.fda.gov/cdrh_docs/reviews/K191562.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf19/K191562.pdf

---

## 胱抑素 C（Cystatin C）

### 1. 法规定位

- Product code **NDY**（Test, Cystatin C），归类于 **21 CFR 862.1225 Creatinine test system**（无独立法规条款），Class II，Panel CH – Clinical Chemistry (75)。
- openFDA 累计清关 **total = 22**（本代码下清关数最少的靶点之一）。
- 伴随代码：JIT（二级校准品）、JIX（多分析物校准品）、JJX（质控品）。

### 2. 代表性 510(k) 一览

| K 号 | 决定日期 | 厂家 / 产品 | 方法学 | 样本类型 | 预期用途关键词 | 文件类型 |
|---|---|---|---|---|---|---|
| K242585 | 2025-05-16 | SENTINEL CH. SpA / Cystatin C（Alinity c 系统） | 颗粒增强免疫比浊 PETIA（兔 IgG 乳胶） | 血清、血浆（SST、Li-hep±胶、Na-hep、K2/K3-EDTA） | aids in the diagnosis and treatment of renal diseases；lab professional；Rx | 决策摘要 + summary |
| K181082 | 2018-05-25 | Siemens / ADVIA Chemistry Cystatin C_2 (CYSC_2) | PETIA（571/805 nm） | 血清、Li-hep、K-EDTA 血浆 | aids in the diagnosis and treatment of renal disease | Triage 决策摘要 + 510(k) summary |
| K171072 | 2017-05-12 | Siemens Healthcare Diagnostics Products GmbH / N Latex Cystatin C + N Protein Standard UY (BN II / BN ProSpec) | 颗粒增强免疫**散射**比浊 PENIA（兔抗） | 血清、Li-heparin 血浆 | used in the diagnosis and treatment of renal diseases；本次目的：**重新标化至 ERM-DA471/IFCC** | Triage 决策摘要 + 510(k) summary |
| K141143 | 2014-07-17 | Roche / Tina-quant Cystatin C Gen.2 + C.f.a.s. Cystatin C + Control Set Gen.2 (cobas c 501) | PETIA（546 nm，兔抗） | 血清、Li-hep、K2/K3-EDTA 血浆 | aid in the diagnosis and treatment of renal diseases | 决策摘要 + summary |
| K161817 | 2016-07-27 | Roche / Tina-quant Cystatin C Gen.2（Special 510(k)） | 同上 | 同上 | 增加样本稳定性声明；增加 HARA 干扰警示；钩状效应限值 20→12 mg/L | Special 510(k) 审查备忘录 + summary |

### 3. 预期用途与声明类型

- 统一措辞「Cystatin C measurements are used as an aid in the diagnosis and treatment of renal diseases」——辅助诊断 + 治疗（监测），无筛查声明；均 Rx、实验室用。**申报文件中均未出现 eGFR 方程（CKD-EPI 2012 cystatin 等）的声明**，只有 K242585 在参考区间人群定义中使用 eGFR >80 mL/min/1.73 m² 作为纳入标准。
- 列表中另有 Diazyme Cystatin C POC Test Kit（K111664, 2012），未读取。

### 4. 阳性 / 阴性判定与 cutoff 逻辑

- 所有文件「Assay cut-off: Not applicable」「Clinical cut-off: Not applicable」；无固定 cutoff，仅参考区间。N Latex（K171072）提到精密度样本覆盖「medical decision point (MDP) of 1.1 mg/L」——这是文件中唯一出现的决策点数值（未作为声明）。

| 产品 | 参考区间（申报） | 人群与方法 |
|---|---|---|
| K242585 Sentinel/Alinity c | **0.59–1.28 mg/L** | 美国健康人 105 女 + 145 男（n=250），**eGFR >80 mL/min/1.73 m²**，18–69 岁；2 台 × 2 批；非参数 2.5–97.5 百分位；CLSI EP28-A3c |
| K181082 Siemens ADVIA CYSC_2 | **0.64–1.23 mg/L** | 208 份「apparently healthy」样本，2.5–97.5 百分位；EP28-A3c；年龄性别未载明 |
| K171072 Siemens N Latex（ERM 重标化后） | **0.49–1.19 mg/L**（21–72 岁） | 203 份美国健康成人或非肾病患者血清，性别均衡，2.5–97.5 百分位；EP28-A3c |
| K141143 Roche Gen.2 | **0.61–0.95 mg/L** | 273 名健康受试者，性别均衡，21–77 岁，2.5–97.5 百分位 |

- **溯源**：全部声明 **ERM-DA471/IFCC**（Sentinel：校准品重量法制备并对 ERM-DA471/IFCC 定值，标准列表含 ISO 17511 2nd ed.；Siemens N Latex：本次 510(k) 的唯一目的即把校准品定值改为溯源 ERM-DA471/IFCC，「materially the same as the assay cleared in K041878」；ADVIA CYSC_2、Roche Gen.2 均溯源 ERM-DA471/IFCC）。
- 逻辑：参考区间由健康人群百分位建立；各厂家上限 0.95–1.28 mg/L 差异较大。Roche 0.61–0.95 与其他 ≈1.2 的差异申报文件未解释（背景：可能与人群、抗体及标化年代有关）。

### 5. 生物学 / 生理学依据

- 申报文件：Cystatin C 用于肾脏疾病的诊断与治疗（Indications）；干扰警示：使用兔抗的方法在接受兔抗体治疗或产生人抗兔抗体（HARA）者可得异常值（K242585 引 Kricka 1999；K161817 因客户投诉新增此警示）；RF 可与试剂免疫球蛋白反应（K242585 引 Boscato 1988）；总蛋白 15 g/dL 引起假性升高（K242585：0.79 mg/L 时 +0.27 mg/L，3.95 mg/L 时 +18.5%）。
- 背景（非申报文件）：Cystatin C 为 13 kDa 半胱氨酸蛋白酶抑制物，所有有核细胞恒定产生，经肾小球自由滤过、近曲小管重吸收并降解，不受肌肉量影响，故用于 eGFR（CKD-EPI 2012 cystatin C、2021 CKD-EPI creatinine-cystatin）；ERM-DA471/IFCC（2010）为国际认证参考物质，标化后各方法结果可用同一 eGFR 方程。

### 6. 样本类型与样本要求

| 产品 | 样本 | 基质等效性 | 稳定性 |
|---|---|---|---|
| K242585 Sentinel | 血清（普通/SST）、Li-hep（±胶）、Na-hep、K2-、K3-EDTA | CLSI **EP35**，50 组（0.6–9.315）PB vs 血清：SST 1.031；Li-hep 1.067；Li-hep 胶 1.067；Na-hep 1.033；K2-EDTA 1.016；K3-EDTA 0.969，r 1.000 | 血清/Li-hep：15–25 °C 2 天；2–8 °C 7 天；−20 °C 30 天；≤2 次冻融 |
| K181082 Siemens ADVIA | 血清、Li-hep、K-EDTA | summary 未载明数据 | 未载明 |
| K171072 Siemens N Latex | 血清、Li-hep | 「unchanged since K041878」 | 未载明 |
| K141143 Roche | 血清、Li-hep、K2-、K3-EDTA、分离胶管 | 57 管/抗凝剂 + 半满管；PB：Li-hep 1.010x+0.020；K2-EDTA 1.020x−0.010；K3-EDTA 1.030x+0.000；胶管 1.000x−0.010，r 1.000 | K161817 新增：血清/血浆 15–25 °C 7 天、2–8 °C 7 天 |

### 7. 分析性能验证所依据的标准

| K 号 | 标准 | 项目 |
|---|---|---|
| K242585 | CLSI EP05-A3、EP06 2nd、EP07 3rd、EP09c 3rd、EP17-A2、EP25 2nd、EP28-A3c、**EP34**（扩展测量区间）、**EP35**（样本类型等效）、**EP37**（干扰补充表）、**ISO 17511 2nd**（溯源） | 精密度、线性、干扰、方法比对、检出限、稳定性、参考区间、EMI、基质、溯源 |
| K181082 | CLSI EP09-A3、EP6-A、EP07-A2、EP28-A3c（summary 正文） | 方法比对、线性、干扰、参考区间 |
| K171072 | CLSI EP05-A2、EP06-A、EP17-A2、EP28-A3c、EP09-A3 | 精密度、线性、LoQ、参考区间、方法比对 |
| K141143 | CLSI EP5-A2、EP6-A、EP17-A2 | 精密度、线性、检出限 |
| K161817 | 无（Special 510(k)） | — |

### 8. 临床/准确度验证设计与结果

- **K242585 Sentinel（Alinity c）**：20 天精密度 n=80：0.49 mg/L within-lab 1.8%；0.92 → 0.9%；5.89 → 0.6%；8.95 → 1.0%。3 站点 × 2 批再现性 n=240：0.49 → 5.4%；0.93 → 1.6%。线性 3 组 × 13 级 0.13–10.76，偏离 ≤9.7%；AMI 0.30–10.00 mg/L，EMI 10–40（1:4 稀释）；前带至 40 mg/L 无。LoB 0.03、LoD 0.05、LoQ 0.30（TE ≤25%，实测 23.2%）。方法比对 vs Roche Tina-quant Gen.2（EP09c）：n=161 血清 0.60–8.20，PB 斜率 1.03 (1.02–1.05)、截距 −0.07、r 1.00。干扰（±10%）：胆红素 60 mg/dL、RF 550 IU/mL、Hb 1000 mg/dL、总蛋白 10.2 g/dL、TG 1500 mg/dL、生物素 3510 ng/mL 等无干扰；RF 1200 IU/mL 在 0.79 mg/L 时 +0.19 mg/L；总蛋白 15 g/dL 显著正干扰（标签限制）。
- **K181082 Siemens ADVIA CYSC_2**：方法比对 vs N Latex（BN ProSpec，EP09-A3）：n=158 冷冻血清 0.44–8.07，PB 斜率 1.008 (0.998–1.022)、截距 0.101 (0.086–0.113)、r 1.000；加权最小二乘斜率 1.021、截距 0.079。精密度 20 天 n=80：0.52 → within-lab 3.7%；1.09 → 1.2%；8.92 → 2.4%。线性 0.25–8.93（斜率 1.087，>2.5 mg/L 偏离 <10%，≤2.5 时 ≤0.25 mg/L）。LoB 0.10、LoD 0.25、LoQ 0.25（<10% CV）。
- **K171072 Siemens N Latex（重标化）**：精密度 20 天 n=80，BN ProSpec：0.465 → 3.85%；1.142 → 1.41%；8.425 → 1.62%。线性 0.27–10.3；LoQ 0.08（<8% CV）；测量区间 0.27–9.4；钩状至 42.91 mg/L 无。方法比对 vs Roche Gen.2（EP09-A3，1 个美国外部站点）：n=186 血清，BN ProSpec r 0.988、斜率 0.982、截距 −0.078；BN II r 0.998、斜率 0.981、截距 −0.079。
- **K141143 Roche Gen.2**：精密度 21 天 n=80：0.560 → within-run 1.8%、between-run 2.0%；2.80 → 0.6/1.3%；6.39 → 0.6/1.1%。线性血清 21 级 0.1–7.6 y=1.001x−0.0057；血浆 13 级 y=1.006x−0.0178；测量区间 0.4–6.8（1:1.5 重测）；钩状至 20 mg/L 无（K161817 改为 12）。LoB 0.3、LoD 0.4、LoQ 0.4（13.3% CV）。方法比对 vs Diazyme（k093680）：n=103（0.42–6.21）PB y=0.997x−0.064, r 0.937；线性回归 y=1.031x−0.153, r 0.988。干扰（≤10%）：胆红素 60 I-index、Hb 1000 H-index、脂血 1000 L-index、RF 1200 IU/mL 及 16 种药物无干扰。

### 9. 厂家间差异与要点

1. **标化是 2014–2017 年这一代 510(k) 的主线**：Siemens N Latex K171072 整份申报只为改用 ERM-DA471/IFCC 定值；Roche Gen.2、ADVIA CYSC_2、Sentinel 均声明 ERM-DA471/IFCC。方法比对均以 Roche Tina-quant Gen.2 或 N Latex 为比较方法，斜率 0.98–1.03。
2. **参考区间**上限差异（0.95 / 1.19 / 1.23 / 1.28 mg/L）；Sentinel 是唯一以 eGFR >80 定义健康人的产品，Siemens N Latex 允许纳入非肾病患者。
3. **eGFR 方程**：任何文件均未把 CKD-EPI cystatin 方程写入预期用途或标签（决策摘要中无 eGFR 声明）；这与肌酐类产品一致，方程由实验室/临床应用。
4. **测量区间**：Sentinel 0.30–10（EMI 至 40）最宽；Roche 0.4–6.8 最窄；PENIA（N Latex）0.27–9.4。
5. **干扰**：兔抗体系（Roche、Siemens、Sentinel）均警示 HARA；Sentinel 额外发现高总蛋白正干扰；Roche 2016 年将钩状效应限值下调至 12 mg/L。
6. **样本**：N Latex 仅血清 + Li-heparin；其余含 EDTA；Sentinel 按 EP35 做 6 种管型等效。

### 10. 来源

- K242585: https://www.accessdata.fda.gov/cdrh_docs/reviews/K242585.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf24/K242585.pdf
- K181082: https://www.accessdata.fda.gov/cdrh_docs/reviews/K181082.pdf （Triage）；https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181082.pdf
- K171072: https://www.accessdata.fda.gov/cdrh_docs/reviews/K171072.pdf （Triage）；https://www.accessdata.fda.gov/cdrh_docs/pdf17/K171072.pdf
- K141143: https://www.accessdata.fda.gov/cdrh_docs/reviews/K141143.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf14/K141143.pdf
- K161817: https://www.accessdata.fda.gov/cdrh_docs/reviews/K161817.pdf ；https://www.accessdata.fda.gov/cdrh_docs/pdf16/K161817.pdf

# 第 E 组 · 自身免疫

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

# 第 F 组 · 感染性疾病血清学

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

---

## 附录 A. 怎样复现或扩展到新的靶点

仓库 `docs/tools/fda_510k_fetch.py` 是本文所有 K 号文本的抓取脚本（只依赖 `pypdf` 和 `curl`）：

```bash
pip install pypdf fonttools
# 1) 查产品代码：用 openFDA classification 接口按 device_name 关键词搜
curl 'https://api.fda.gov/device/classification.json?search=device_name:"troponin"&limit=10'
# 2) 列出该产品代码最近的 510(k)
python3 docs/tools/fda_510k_fetch.py list MMI 40
# 3) 下载决策摘要 + 510(k) summary 并抽成文本
python3 docs/tools/fda_510k_fetch.py get K231974 K252393
# 4) 阅读 txt/K231974.txt，按本文 §1–§10 模板填写
```

注意事项：
- `accessdata.fda.gov` 对无 User-Agent 的请求返回反爬页面，脚本已带浏览器 UA。
- 决策摘要（`cdrh_docs/reviews/`）从约 2003 年起才有，且并非每个 K 号都有；更早或缺失的只能用 510(k) summary（`cdrh_docs/pdfYY/`），信息量明显少。
- 修改型申报（Special 510(k) / 换平台）常写 “previously established in Kxxxxxx”，需要沿 predicate 链回溯到原始申报才能拿到 cutoff 的建立过程。

## 附录 B. 引用的 FDA 指南

- Statistical Guidance on Reporting Results from Studies Evaluating Diagnostic Tests（2007-03-13）。
- Review Criteria for Assessment of C-Reactive Protein (CRP), High Sensitivity C-Reactive Protein (hsCRP) and Cardiac C-Reactive Protein (cCRP) Assays（2005-09-22）：明确 cCRP 决策点 ≤ 1.0 mg/L、AHA/CDC 低风险 cutoff 1.0 mg/L 处 CV ≤ 10%、标准化到 IFCC/BCR/CAP CRM 470、允许“文献 + 桥接研究”支持风险分层声明。
- Guidance Document for the Submission of Tumor Associated Antigen Premarket Notifications (510(k)s)（21 CFR 866.6010 特殊控制）。
- Recommendations for CLIA Waiver Applications（POC 产品）。
- 21 CFR 862.1117（natriuretic peptide）、862.1215（cTn，挂靠 CK 条款）、862.1225（cystatin C）、862.1690（TSH）、862.1825/1840（vitamin D）、864.7320（fibrin degradation products）、864.7470（HbA1c）、866.3175/3305/3830（CMV、HSV、Borrelia/Treponema 血清学）、866.3215（PCT）、866.5100（ANA）、866.5180（钙卫蛋白）、866.5270（CRP）、866.5340（ferritin）、866.5775（anti-CCP）、866.6010（肿瘤标志物）、866.6050（卵巢附件肿块评估）。
