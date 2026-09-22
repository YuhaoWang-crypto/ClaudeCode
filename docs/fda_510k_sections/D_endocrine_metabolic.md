# 靶点组 D：内分泌 / 代谢（HbA1c、TSH、25-OH 维生素 D、铁蛋白、胱抑素 C）

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
