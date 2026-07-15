# G034 项目 — 抗原鉴定与靶点评估报告

> **一句话结论:** G034 提供的抗原序列 = **人 α‑taxilin(TXLNA,UniProt P40222)第 210–546 位残基**(BLAST 100% 一致)。该蛋白的"第二身份"正是 **白介素‑14α(IL‑14α)/ 高分子量 B 细胞生长因子(HMW‑BCGF)**——一个在干燥综合征(Sjögren's disease, SjD)动物模型和患者血清中都有验证的 **B 细胞高活化驱动因子**。靶点**新颖度极高、专利空白极大**(与背景"希望实现专利突破"高度契合);但作为"阻断型全长抗体药"的靶点,存在一个**必须先回答的根本性成药学问题**:α‑taxilin 经典身份是**胞内**囊泡运输蛋白、无经典信号肽,而抗体只能作用于**胞外/可溶**的部分。整个项目成立与否,取决于"可溶性 IL‑14α 是否真实、可被中和"这一命题能否被验证。

本报告严格区分 **✅ 已确证** / **⚠️ 假说或有争议** / **❗ 关键风险**,逐条回应甲方 4 项需求(见文末《需求对照表》)。

---

## 0. 抗原鉴定(需求 2:抗原信息)

### 0.1 序列与翻译

- 提供的 DNA:**1014 nt**,ATG 起始、TAG 终止,可读框完整(1014 / 3 = 338 密码子)。
- 翻译产物:**337 aa**(去除终止子)。

```
MKLLQKKQSQLVQEKDHLRGEHSKAVLARSKLESLCRELQRHNRSLKEEGVQRAREEEEKRKEVTSHFQVTLND
IQLQMEQHNERNSKLRQENMELAERLKKLIEQYELREEHIDKVFKHKDLQQQLVDAKLQQAQEMLKEAEERHQR
EKDFLLKEAVESQRMCELMKQQETHLKQQLALYTEKFEEFQNTLSKSSEVFTTFKQEMEKMTKKIKKLEKETTM
YRSRWESSNKALLEMAEEKTVRDKELEGLQVKIQRLEKLCRALQTERNDLNKRVQDLSAGGQGSLTDSGPERRP
EGPGAQAPSSPRVTEAPCYPGAPSTEASGQTGPQEPTSARA
```

### 0.2 BLAST 鉴定结果(✅ 已确证)

| 项目 | 结果 |
|---|---|
| 最佳匹配 | **Alpha‑taxilin,Homo sapiens**(NP_001363786 / **UniProt P40222**;基因 **TXLNA**) |
| 一致度 | **337 / 337 = 100%**,E‑value = 0.0,无 gap |
| 比对位置 | 覆盖全长 546 aa 蛋白的 **第 210–546 位**(即 C 端约 62%) |
| 跨物种 | 与猕猴、猩猩、黑猩猩等灵长类 100% 一致 → **高度保守**(❗ 见 5.4 免疫耐受/安全性) |

> 换言之,G034 抗原是 **α‑taxilin / IL‑14α 的 C 端片段**,涵盖了功能性的 **coiled‑coil(卷曲螺旋)syntaxin 结合域** + C 端富脯氨酸低复杂度尾巴。结构上:约前 275 aa 为规则七肽重复卷曲螺旋,C 端约 60 aa 为无序富脯氨酸区(PRD)。

### 0.3 关键的"双重身份"(✅ 数据库一致标注)

TXLNA / α‑taxilin 在 HGNC / UniProt / ChEMBL / Reactome 中同时携带以下别名与注释:

- **别名:** Interleukin‑14(IL‑14)、IL‑14α、HMW‑BCGF(高分子量 B 细胞生长因子)、TXLN‑α。
- **GO 注释(ChEMBL/UniProt,并存):** 一方面 `cytoplasm / cytosol / syntaxin binding / exocytosis`(胞内运输身份);另一方面 `extracellular region / B cell activation / "Other interleukin signaling"(Reactome R‑HSA‑449836)`(胞外细胞因子身份)。
- Reactome 明确记载:*"Interleukin‑14, renamed alpha‑taxilin (TXLNA)…"*。

**这一"胞内运输蛋白 ↔ 分泌型 B 细胞生长因子"的双重身份,是本项目全部机会与风险的根源。** 下文分别展开。

---

## 1. 靶点生物学(需求 1 + 项目背景:干燥症抗体药)

### 1.1 身份 A — α‑taxilin:胞内 syntaxin 结合 / 囊泡运输蛋白(✅ 机制清楚)

- α‑taxilin 通过 C 端卷曲螺旋结合 syntaxin‑1/3/4,参与 SNARE 复合体调控下的 **胞吐(exocytosis)与内体循环**(如转铁蛋白受体经 SNX4 的循环通路);在神经内分泌细胞中参与 **Ca²⁺ 依赖性胞吐**。
- **泛表达**、定位于胞质;taxilin 家族含 α/β/γ 三成员。
- 该身份下 α‑taxilin **无经典信号肽、非分泌**——这是后文成药性讨论的核心约束。

> 备注(与疾病的潜在机制联系,⚠️ 推测):泪腺/唾液腺的液体分泌本质上依赖 SNARE 介导的胞吐(水通道蛋白转运、分泌颗粒融合)。因此"胞内 α‑taxilin 功能异常 → 腺体分泌障碍"在机制上并非不可能,但这属于**胞内机制**,**无法用抗体药干预**。

### 1.2 身份 B — IL‑14α:B 细胞生长因子 / 干燥症驱动因子(✅ 动物+人证据,⚠️ 分子机制有争议)

历史脉络:HMW‑BCGF/IL‑14 由 Ambrus 等在 1980s–90s 描述为一种促记忆 B 细胞增殖、上调自身受体、抑制抗体分泌的高分子量 B 细胞生长因子(据 PubMed,DOI [10.1016/S0021-9258(19)67851-8](https://doi.org/10.1016/S0021-9258(19)67851-8);DOI [10.1097/00045391-199512000-00006](https://doi.org/10.1097/00045391-199512000-00006))。后被归入 il14 基因,IL‑14α 由 plus 链外显子 3–10 编码,并与 taxilin 命名合并。**IL‑14 从未进入正式白介素命名体系的"主流",其作为独立分泌型细胞因子的分子身份在领域内一直存在争议**(⚠️)。

尽管分子机制有争议,**功能性证据链却相当完整且可重复**:

| 证据 | 内容 | 强度 |
|---|---|---|
| **IL‑14α 转基因鼠(IL‑14αTG)** | 组成型过表达人 IL‑14α 的小鼠,按与患者相同的时间进程依次出现:高丙种球蛋白血症 → 自身抗体 → 唾液腺功能下降 + 淋巴细胞浸润 → 累及颌下腺/腮腺/**泪腺** + 间质性肺病 → 晚期 **大 B 细胞淋巴瘤**。**完整复现原发性干燥综合征(pSS)全过程。** 据 PubMed,DOI [10.4049/jimmunol.177.8.5676](https://doi.org/10.4049/jimmunol.177.8.5676);DOI [10.1016/j.clim.2008.10.006](https://doi.org/10.1016/j.clim.2008.10.006) | ✅ 强(公认的 SjD 模型) |
| **遗传上位实验** | 在 IL‑14αTG 鼠中特异性敲除边缘区 B 细胞(MZB,删 RBP‑J)→ **全部 SjD 表型消失**(唾液分泌正常、自身抗体阴性、腺体组织学正常);而删 B1 细胞(btk)无效 → 证明 IL‑14α 通过 MZB 驱动疾病。据 PubMed,DOI [10.1016/j.clim.2016.04.008](https://doi.org/10.1016/j.clim.2016.04.008) | ✅ 强(因果) |
| **人血清生物标志物** | pSS 患者血清 IL‑14α(Western blot)显著高于非 SS 干眼、类风关、健康对照;与抗 SSA/Ro、抗 SSB/La 相关;与 BAFF 并列作为干眼分层标志物。据 PubMed,DOI [10.3389/fimmu.2021.673658](https://doi.org/10.3389/fimmu.2021.673658) | ✅ 中(人相关性,非因果) |
| **2025 最新** | 唾液腺转录组 + 免疫表型研究仍以 IL‑14αTG 为标准 SjD 模型;当前治疗讨论集中于抗 CD20 等 B 细胞疗法,**尚无直接靶向 IL‑14α 的疗法** → 说明该治疗方向仍是**空白**。 | ✅ |

**小结(靶点‑疾病相关性):** 在"IL‑14α 是 pSS/干燥症关键上游驱动因子"这一层面,**有转基因模型的因果证据 + 人血清相关性证据**,验证程度**高于绝大多数"全新靶点"**。这对甲方"用于干燥症的抗体药"背景是强正向信号。

---

## 2. 文献调研 A:靶点新颖性(需求 1)

| 维度 | 评估 | 依据 |
|---|---|---|
| **作为药物靶点的新颖性** | **极高(近乎全新)** | ClinicalTrials.gov 检索"taxilin/IL‑14"干预性试验 = **0**;无任何抗 IL‑14 抗体进入临床;2025 年综述仍不把 IL‑14α 列为在研治疗靶点。 |
| **疾病生物学验证度** | **中‑高** | 转基因鼠因果证据 + MZB 上位实验 + 人血清标志物(见 1.2)。 |
| **分子机制清晰度** | **低‑中(有争议)** | "分泌型细胞因子 IL‑14α"与"胞内运输蛋白 α‑taxilin"身份未统一;受体、信号通路未被现代分子手段确证(⚠️)。 |
| **竞争格局** | **几乎无竞争** | 干燥症在研生物药集中在 BAFF/BAFF‑R(ianalumab)、BAFF+APRIL(telitacicept)、CD40/CD40L、抗 CD20、LTβR、IL‑17 等成熟通路;**IL‑14α 方向无企业管线**。 |

> **新颖性的两面性:** 新颖 = 专利/竞争空白大(见第 3 节),但也 = **靶点未被行业充分验证**,机制风险需自行承担。这是典型的"first‑in‑class 高风险高回报"靶点。

---

## 3. 文献调研 B:专利态势与"专利突破"空间(需求 1)

甲方明确"重点……以及希望实现专利突破"。结论:**这是本靶点最有吸引力的一面。**

### 3.1 现有相关专利 / 现有技术(prior art)

| 类型 | 关键条目 | 对 FTO / 突破的影响 |
|---|---|---|
| **基础专利** | **US 7,622,574** — *"IL‑14α RNA inhibitors and antibodies to IL‑14α for treatment of autoimmune diseases and lymphomas"*(Ambrus / University at Buffalo 体系)。权利要求覆盖:抗 IL‑14α 多肽的单抗 + 反义寡核苷酸,用于 SLE、Sjögren、淋巴瘤。 | ❗ 需重点规避的**方法/用途**核心现有技术。但该专利优先权约在 2000s、授权 2009,**20 年期限已到期或临近到期**——意味着**用途保护正在或即将失效**,后来者反而**更易获得自由操作空间(FTO)**。需由专利律师核实法律状态与同族。 |
| **公开的抗体现有技术** | Peng 等 2009:用合成 **IL‑14α‑C 多肽**免疫 BALB/c 得到鼠抗 IL‑14α‑C 单抗(IgG2a/κ,Kaff≈1×10⁸ M⁻¹≈10 nM)。据 PubMed,DOI [10.1089/hyb.2009.0007](https://doi.org/10.1089/hyb.2009.0007) | ⚠️ 注意:该抗体针对的 **"IL‑14α‑C" 正是 C 端片段,与 G034 抗原(210–546)高度重叠**。这是"能否成药"的正面证据(C 端可成功免疫出特异抗体),但也构成**表位层面的现有技术**——新抗体应通过独特 CDR / 表位实现区分。 |
| **诊断专利** | University at Buffalo 另有 SjD 诊断(IL‑14α 检测)相关专利。 | 与治疗性抗体权利要求不同域,冲突有限。 |
| **商业研究抗体** | Proteintech、R&D、OriGene、Abbexa 等大量 **科研级** 抗 TXLNA 多/单抗(WB/ELISA/IHC 用)。 | 均非治疗级、无序列专利,不构成 composition‑of‑matter 障碍;可用作检测/对照工具。 |

### 3.2 专利突破可行性(强)

- **Composition‑of‑matter(序列/CDR)空白:** 目前**没有**针对 IL‑14α 的**人源化治疗性抗体的序列(CDR)专利**在有效期内占据。自研人源化全长抗体 → **可获得全新的物质专利**,这是最干净的突破点。
- **用途专利突破:** 基础用途专利(US7622574)若已到期/临近到期,"抗 IL‑14α 抗体治疗干燥症/干眼"可通过**新适应症细分、新给药途径(如局部滴眼)、新剂型、联用**等重构用途权利要求。
- **规避策略:** 选择与 IL‑14α‑C 现有抗体**不同的表位**(例如靶向 coiled‑coil 中段特定构象表位,而非 C 端富脯氨酸尾),既提升新颖性又利于差异化功能。

> **专利小结:** ❗ 需专业检索确认 US7622574 及其同族的**法律状态与到期日**;若确认到期,则本靶点 = **专利空白 + 可建物质专利**,与"专利突破"诉求**高度契合**。

---

## 4. 抗体形式与用途可行性(需求 3:阻断抗体 / 成药 / 人源化 / 全长)

甲方要求:**阻断抗体、需成药、人源化、全长(带 Fc)**。逐项评估:

### 4.1 ❗ 核心成药学问题:靶点可及性(make‑or‑break)

- 全长 IgG **不能进入细胞质**。若 IL‑14α 的致病作用发生在**胞内(α‑taxilin 运输功能)**,则**任何阻断抗体都无法奏效**。
- 项目成立的**必要前提**是:存在**足量、可及的胞外/可溶性 IL‑14α**,且其 B 细胞刺激活性可被中和。
- **支持前提成立的证据:** ① 患者**血清**可检出 IL‑14α(Liang 2021,WB);② IL‑14 历史上被描述为由滤泡树突细胞/活化 T 细胞**分泌**;③ 转基因**过表达分泌型构建**即致病 → 提示胞外 IL‑14α 有生物学活性。
- **削弱前提的因素(⚠️):** ① 经典 α‑taxilin **无信号肽**,分泌机制若存在则为**非经典分泌**,尚未被现代方法定量确证;② 血清中被 WB 识别的"IL‑14α"分子种类/来源仍需严格定性;③ 领域对"IL‑14 是独立分泌细胞因子"仍有保留。

> **建议(见第 6 节):在投入完整抗体发现前,先用 1 个"靶点可及性验证包"回答此问题。** 这是决定 go/no‑go 的第一道闸门。

### 4.2 抗体形式建议

| 需求 | 评估与建议 |
|---|---|
| **阻断/中和** | 需要**功能性中和 assay**(如 IL‑14α 诱导的 B 细胞增殖 / MZB 扩增被抗体抑制)来定义"阻断"。抗原 C 端富脯氨酸尾为无序区,免疫原性好但未必是功能表位;**优先筛选靶向 coiled‑coil 功能域(syntaxin/受体结合面)的中和抗体**。 |
| **全长(带 Fc)** | 与"阻断可溶性配体"目标一致(长半衰期、便于全身给药)。若担心 Fc 效应功能带来的耗竭/ADCC(此处目标是中和配体、非杀细胞),可选 **IgG1 LALA / IgG4** 等**效应减弱 Fc**。若未来考虑局部滴眼,再评估片段化(但甲方明确要"全长",故以全长为主线)。 |
| **人源化** | 现有现有技术抗体均为**鼠源**(Peng 2009 等)→ 自研**人源化/全人**抗体既满足成药要求,又是**新颖性与专利突破**的核心。路线:杂交瘤/免疫文库 → 人源化 CDR 移植,或直接用全人转基因鼠/噬菌体人源文库。 |
| **靶点保守性带来的免疫难点(❗)** | α‑taxilin 跨物种 100% 保守(见 0.2)。用野生型小鼠免疫**人 IL‑14α‑C** 时,若鼠源同源蛋白序列高度一致 → 可能**免疫耐受、滴度低**。缓解:用 **KLH 偶联多肽 + 强佐剂**、选择人鼠差异表位、或用 **TXLNA‑KO 小鼠**免疫;噬菌体全人文库可绕开耐受。 |

### 4.3 用途:治疗 vs 检测

- 甲方本项目 = **阻断(治疗)抗体**。
- 顺带机会:同一抗原可衍生**伴随诊断**(血清 IL‑14α 分层干眼/pSS,Liang 2021 已示范)——利于整体 IP 组合与临床分层,但非本次交付主线。

---

## 5. 交付标准 / 开发性评估(需求 4:免疫原性 / 稳定性 / 可合成性,"都看")

> 说明:开发性(developability)评估的对象是**未来产出的抗体分子**,当前尚无候选序列,故以下为**评估框架 + 针对本靶点/本抗原的具体注意点**;抗体序列产出后再跑定量计算(如 TAP、AggScore、Therapeutic Antibody Profiler 等)。

### 5.1 免疫原性(anti‑drug antibody, ADA 风险)

- 抗体侧:走**人源化/全人 + 低 T 细胞表位**路线(EpiVax/IEDB/NetMHCII in‑silico 去免疫原性);全长人 IgG1 骨架 ADA 风险可控。
- 抗原侧(与本任务直接相关的一点):❗ 靶点是**人自身蛋白**且**高度保守**。这带来两点——(a)免疫产抗难(见 4.2);(b)理论上,治疗性抗体若与患者自身广泛表达的 α‑taxilin 交叉,需评估**打破自身耐受/免疫复合物**风险。**注:** 我们**未**对本抗原运行 EDEN 免疫原性预测器,因为该模型要求**完整天然 CDS**、明确不适用于**片段(partial)序列**(本抗原为 210–546 片段,属分布外),强行预测会给出误导性数值——此处如实弃用,建议后续用**全长天然 CDS**再评估。

### 5.2 稳定性(stability)

- 分子层面:抗体产出后评估热稳定(Tm/DSF)、胶体稳定(自相互作用、pI)、化学降解热点(CDR 的 Asn 脱酰胺 NG/NS、Asp 异构 DG、Met/Trp 氧化、游离 Cys)。
- 与本抗原相关点:抗原 C 端含 **1 个游离 Cys(…VTEAP**C**YPGA…,对应 P40222 的 C531)** 及富脯氨酸区 —— 免疫/重组抗原制备时需注意该 Cys 的**错配二硫/聚集**倾向,建议对抗原做 Cys→Ser 突变或加帽以获得均一免疫原(不影响筛抗体)。

### 5.3 可合成性 / 可生产性(manufacturability / synthesizability)

- 标准 CHO 表达全长 IgG 的可生产性主要由序列决定(表达量、聚集、疏水性、CDR liabilities),产出候选后用 **Therapeutic Antibody Profiler(TAP)/ Developability Index** 打分即可。
- 抗原生产(用于免疫与筛选):C 端片段(210–546)为**卷曲螺旋 + 无序尾**,原核表达易形成包涵体/聚集;建议 **哺乳或杆状病毒表达 + 融合标签(Fc/His/SUMO)**,或分段/加可溶化标签,并处理 5.2 的游离 Cys。

### 5.4 ❗ 安全性(甲方未单列,但"成药/都看"应涵盖)

- α‑taxilin **泛表达且为基础囊泡运输机器**;虽抗体不入胞,但需评估:(a)是否存在**细胞表面/膜结合池**导致 on‑target 结合正常组织;(b)中和 IL‑14α 对**正常记忆 B 细胞维持**的影响(免疫抑制/感染风险)。这些应纳入早期毒理与靶点安全性评估。

---

## 6. 关键风险与建议(整体 go / no‑go)

### 6.1 风险矩阵

| 风险 | 等级 | 说明 |
|---|---|---|
| **靶点可及性(胞外可溶 IL‑14α 是否足量、可中和)** | ❗❗❗ 最高 | 决定"阻断全长抗体"路线是否根本成立(4.1)。 |
| **分子机制/受体未确证** | ❗❗ 高 | 难以建立干净的机制性中和 assay、难做 PK/PD。 |
| **自身抗原免疫耐受(产抗难)** | ❗ 中 | 影响发现效率;可用文库/KO 鼠缓解(4.2)。 |
| **on‑target 安全性(泛表达)** | ❗ 中 | 需早期评估(5.4)。 |
| **专利法律状态未核实** | ❗ 中 | US7622574 到期日需律师确认(3.1)。 |
| 靶点新颖性 / 专利空白 | ✅ 正向 | 竞争空白、可建物质专利(2、3)。 |
| 疾病相关性验证 | ✅ 正向 | 转基因鼠因果 + 人血清证据(1.2)。 |

### 6.2 建议的分阶段路线(先验证、后投入)

**Gate 1 — 靶点可及性验证包(最小投入、最高价值,建议先做):**
1. 定性患者血清/唾液中的 IL‑14α 分子种类(免疫沉淀 + 质谱),确认可溶species 与来源;
2. 重组表达可溶 IL‑14α(全长 + 210–546),建立 **B 细胞/MZB 增殖中和 assay**;
3. 用现有鼠抗 IL‑14α‑C(或快速多抗)做**概念验证中和**——能否抑制 IL‑14α 诱导的 B 细胞活化。
→ **若 Gate 1 阳性**:靶点可及性成立,进入正式抗体发现;**若阴性**:及时止损或转向胞内策略(非抗体)。

**Gate 2 — 抗体发现与人源化(满足甲方需求 3):**
- 全人噬菌体/酵母文库 或 免疫(KO 鼠/差异表位)→ 杂交瘤 → 人源化;
- 优先 coiled‑coil 功能域中和表位,避开 IL‑14α‑C 现有表位以利专利;
- 输出:人源化 **全长 IgG(效应减弱 Fc 可选)**,附亲和力 + 功能中和数据。

**Gate 3 — 开发性 + 专利(满足需求 4 + "专利突破"):**
- TAP/Developability 打分 + 去免疫原性 + 稳定性液相条件;
- 同步提交**物质专利(CDR)**;核实并绕开 US7622574。

### 6.3 一句话建议

> **靶点科学吸引力与专利空白俱佳,但成药性风险集中于一个可先验证的命题。** 建议**先花小钱做 Gate 1 靶点可及性验证**,以其结果决定是否全力投入抗体发现;不建议在未验证可溶 IL‑14α 之前直接启动完整人源化抗体 campaign。

---

## 7. 补充问答:成药前提的实验证据 + 若为胞内的替代模态(答甲方追问)

### 7.1 「可溶性 IL‑14α 存在且可被中和」是否有实验证据?

**有,且比预期更强——但存在一条未闭合的分子鸿沟。**

**✅ 支持证据:**

| 证据 | 内容 | 意义 |
|---|---|---|
| Ambrus 1985(PMID 3876402,*J Immunol*) | 从 **T 细胞系上清中把 HMW‑BCGF 纯化到均一**,证明**特异结合活化 B 细胞**,并制得**能从上清吸附/中和其活性的单抗** | **最硬证据**:存在一个分泌的、有受体、可被抗体中和的 B 细胞生长因子 |
| Blood 1995(PMID 7795235) | 侵袭性 B 细胞淋巴瘤患者**积液中检出 IL‑14/HMW‑BCGF** | 体液天然存在 |
| Liang 2021(DOI [10.3389/fimmu.2021.673658](https://doi.org/10.3389/fimmu.2021.673658));anti‑PD‑1 2023 | pSS **血清** IL‑14α 升高;肿瘤血清 IL‑14α 预测 PD‑1 疗效 | 循环可检测 |

**⚠️ 未闭合的鸿沟(须如实告知):**

- 1993 年最初的 IL‑14 cDNA(PMID 8327514,PNAS)预测的是**带信号肽、53 kDa、3 个 N‑糖位点的分泌蛋白**;而如今等同为 IL‑14α 的 **α‑taxilin(P40222,~62 kDa,无信号肽、无 N‑糖)不是显然同一分子**。"纯化到的分泌活性 = taxilin 基因产物"这一分子桥**从未用现代方法干净证明**。
- α‑taxilin **无经典信号肽** → 若分泌则为**非经典/无前导分泌**(如 IL‑1/FGF2/HMGB1),或血清中只是**细胞裂解被动释放**;二者未被实验区分。
- **IL‑14 受体从未被分子克隆鉴定。**

> 结论:做阻断抗体的**必要前提有正向实验证据**(存在可中和的分泌活性),但"该活性 = α‑taxilin 本身"及其分泌途径**尚未证明**——正是 §6.2 Gate‑1(IP‑质谱定性 + 中和 assay)要闭合的缺口。

### 7.2 阻断 α‑taxilin 与「什么」结合?(靶点相互作用)

本抗原 210–546 = C 端 coiled‑coil,**正是相互作用结构域**。

- **胞内(运输轴,已验证):** 结合**游离 Syntaxin‑1A/3A/4A**(不结合已成 SNARE 者)→ **阻止 t‑SNARE 形成、抑制胞吐**(Reactome "TXLNA(IL14)binds syntaxin3");另结合 **SNX4**(转铁蛋白受体循环)、**NAC 复合体**、**聚合微管蛋白**。❗方向性:α‑taxilin 是胞吐**抑制因子**,清除它反而**增强**胞吐。
- **胞外(IL‑14 轴,治疗真正想打):** 结合**活化 B 细胞上推定的"IL‑14R"**(功能上有、分子身份未鉴定)。抗体的目标 = 阻断 **IL‑14α ↔ IL‑14R** → 抑制 MZB/B 细胞过度活化。

### 7.3 若确为胞内——更适合的"阻断/清除"模态

| 形式 | 适配 | 说明 |
|---|---|---|
| **ASO / siRNA(敲低)** | ⭐⭐⭐ 首选 | 对无口袋胞内靶点最现实;❗ **US7622574 已 claim "IL‑14α RNA inhibitors"**,需新序列/化学/递送差异化;难点=淋巴/腺体递送 |
| **降解剂 PROTAC / 分子胶** | ⭐⭐ | 见 7.4;机制契合但卡在无配体 |
| **胞内抗体 intrabody / VHH / scFv** | ⭐ | 需基因/mRNA 递送,近研究工具 |
| **多肽 / 订书肽 PPI 抑制剂** | ⭐ | 模拟 syntaxin 螺旋竞争阻断;入胞是瓶颈 |

> 对胞内靶点,**"降解/敲低"优于"阻断界面"**:一次清除全部功能。现实首选 **ASO/siRNA**。

### 7.4 PROTAC 思路可行性

**方向正确(胞质蛋白→泛素‑蛋白酶体就在旁边,正好解决抗体的区室难题;催化性、清除全部功能),但当前被"无配体"卡死:**

- ❗ ChEMBL 全库 α‑taxilin **仅 1 条小分子结合记录**:激酶抑制剂样化合物,**Kd ≈ 3.5 µM**,且只是 **kinobead 化学蛋白组学脱靶命中**(ChemMedChem 2018)——弱、偶然、未优化,非可成药口袋。
- 结构 = 长 coiled‑coil + 无序 PRD,典型小分子"难成药"拓扑。
- → **现成 PROTAC 不可行,须先做苗头发现。**

**具体抓手:** 抗原含 **4 个 Cys(P40222:C245 / C373 / C471 / C523,后者在 C 端)** → **针对 Cys 的共价片段筛选**可获共价锚点 → 搭 **共价 PROTAC / 共价分子胶**;3.5 µM kinobead 命中亦为弱起点。属真实但**数年级**化学工程。
**其他降解:** 分子胶(更难理性设计);若存在分泌/表面池 → **LYTAC / 抗体基降解剂**可**复用本项目要做的抗体**降解胞外 α‑taxilin(收敛点)。

### 7.5 模态选择取决于 Gate‑1

| Gate‑1 结果 | 推荐模态 |
|---|---|
| 可溶 IL‑14α 真实且可中和 | ✅ **阻断抗体**(甲方首选)→ 打 IL‑14α–IL‑14R 轴 |
| 纯胞内 | **ASO/siRNA(主)** 或 **降解剂(需配体 campaign)**;抗体无效 |
| 不确定/对冲 | 抗体 + 敲低**双轨** |

---

## 需求对照表(逐条回应甲方 4 点 + 背景)

| 甲方需求 | 本报告回应 | 位置 |
|---|---|---|
| 项目背景:开发用于**干燥症**的抗体药 | 靶点 = IL‑14α,pSS/干燥症上游 B 细胞驱动因子,有转基因鼠 + 人血清验证 | §1.2 |
| 需求 1:文献调研(靶点新颖性 **+** 抗体专利 / 专利突破)——"都看" | 新颖性极高、竞争空白;基础专利 US7622574 可能到期→物质专利突破空间大 | §2、§3 |
| 需求 2:抗原信息("见文件") | = α‑taxilin/IL‑14α **210–546**,BLAST 100% 确证,附结构解析 | §0 |
| 需求 3:阻断抗体 / 成药 / 人源化 / 全长 | 全长带 Fc 与"中和可溶配体"一致;需人源化(现有技术均鼠源);**核心前提 = 胞外可及性** | §4 |
| 需求 4:开发性评估(免疫原性 / 稳定性 / 可合成性)——"都看" | 给出评估框架 + 本抗原具体注意点(游离 Cys、保守性耐受、片段表达);安全性补充 | §5 |

---

## 参考文献(据 PubMed / ClinicalTrials.gov / Google Patents / UniProt / ChEMBL)

> 本报告文献部分信息来自 **PubMed**,按其要求附 DOI 链接。

1. Shen L, *et al.* Development of autoimmunity in IL‑14alpha‑transgenic mice. *J Immunol.* 2006. DOI [10.4049/jimmunol.177.8.5676](https://doi.org/10.4049/jimmunol.177.8.5676) (PMID 17015757)
2. Shen L, *et al.* IL‑14 alpha, the nexus for primary Sjögren's disease in mice and humans. *Clin Immunol.* 2008. DOI [10.1016/j.clim.2008.10.006](https://doi.org/10.1016/j.clim.2008.10.006) (PMID 19038581)
3. Shen L, *et al.* Central role for marginal zone B cells in an animal model of Sjögren's syndrome. *Clin Immunol.* 2016. DOI [10.1016/j.clim.2016.04.008](https://doi.org/10.1016/j.clim.2016.04.008) (PMID 27140729)
4. Liang Y, *et al.* IL‑14α as a Putative Biomarker for Stratification of Dry Eye in Primary Sjögren's Syndrome. *Front Immunol.* 2021. DOI [10.3389/fimmu.2021.673658](https://doi.org/10.3389/fimmu.2021.673658) (PMID 34012457)
5. Peng X, *et al.* Characteristics of a novel monoclonal antibody against interleukin‑14alpha. *Hybridoma.* 2009. DOI [10.1089/hyb.2009.0007](https://doi.org/10.1089/hyb.2009.0007) (PMID 19663695)
6. Ambrus JL, *et al.* Intracellular signaling…HMW‑BCGF. *J Biol Chem.* 1991. DOI [10.1016/S0021-9258(19)67851-8](https://doi.org/10.1016/S0021-9258(19)67851-8) (PMID 1847385)
7. Ambrus JL, *et al.* A Potential Role for PGE and IL‑14 (HMW‑BCGF) in B‑Cell Hyperactivity of SLE. *Am J Ther.* 1995. DOI [10.1097/00045391-199512000-00006](https://doi.org/10.1097/00045391-199512000-00006) (PMID 11854811)
8. 专利 **US 7,622,574** — *IL‑14α RNA inhibitors and antibodies to IL‑14α for treatment of autoimmune diseases and lymphomas*(需核实法律状态/到期日)。
9. 专利 **WO2021091706A1** — *Treatment for Sjögren's syndrome*(Google Patents,内容待进一步核实)。
10. UniProt **P40222**(Alpha‑taxilin);ChEMBL Target **CHEMBL6066423**;基因 **TXLNA**(HGNC:30685)。
11. Frontiers in Dental Medicine 2025 — Salivary gland transcriptomic analysis and immunophenotyping in the IL‑14α transgenic mouse model of Sjögren's disease(2025 年最新,确认 IL‑14αTG 仍为标准模型,尚无直接抗 IL‑14α 疗法)。
12. Ambrus JL, *et al.* Purification to homogeneity of a high molecular weight human B cell growth factor; specific binding to activated B cells; and a monoclonal antibody to the factor. *J Immunol.* 1985. (PMID 3876402) — 分泌 BCGF 纯化 + 中和抗体的直接实验证据。
13. Ford R, Ambrus JL, *et al.* Identification of a cDNA for a human high‑molecular‑weight B‑cell growth factor (proposed IL‑14). *PNAS.* 1993;90:6330. (PMID 8327514) — 最初预测"带信号肽 53 kDa 分泌蛋白"的 cDNA(与现 α‑taxilin 分子身份存在鸿沟)。
14. Ambrus JL, *et al.* Identification of B‑cell growth factors (IL‑14; HMW‑BCGF) in effusion fluids from aggressive B‑cell lymphomas. *Blood.* 1995;86:283. (PMID 7795235) — 体液中检出。
15. α‑taxilin 的 syntaxin/SNX4 结合:Nogami *et al.* Taxilin, a novel syntaxin‑binding protein (PMID 12558796);α‑Taxilin–SNX4 与转铁蛋白受体循环,*PLoS One* 2014,DOI [10.1371/journal.pone.0093509](https://doi.org/10.1371/journal.pone.0093509)(PMID 24690921);Reactome R‑HSA‑9014052 "TXLNA(IL14)binds syntaxin3"。
16. 小分子可配体性:ChEMBL Target CHEMBL6066423 仅 1 条结合记录(molecule CHEMBL5653589,Kd≈3.5 µM,kinobead 脱靶,*ChemMedChem* 2018)——提示 α‑taxilin 目前无可成药口袋,PROTAC 需先做共价苗头发现(靶向 C245/C373/C471/C523)。

---

*方法说明:抗原经标准密码子表翻译为 337 aa,提交 NCBI blastp(nr 库)鉴定;靶点/疾病/专利信息经 PubMed、ClinicalTrials.gov、ChEMBL、UniProt、Google Patents / USPTO 检索交叉核对。全部结论按 ✅ 已确证 / ⚠️ 争议或假说 / ❗ 风险 标注,未做无依据的断言。*
