# 攻克 KRAS：从"不可成药"到隐匿口袋与分子胶三元复合物

> 一份把 KRAS 成药史与现代计算药物发现流程串起来的技术笔记，并附 4 个基于 Boltz‑2 的可复现 demo（见 [`demo/`](demo/)）。
>
> 📘 **想要一份带目录、串起全部内容（科学背景 + 计算 demo + 设计 playbook + 4 个设计 pilot）的总报告？见 [`REPORT.md`](REPORT.md)。** 本文件是其中"科学背景"的详版。

---

## 目录

1. [概述：为什么 KRAS 曾是"圣杯"级难题](#1-概述)
2. [三重壁垒：传统药物为何"抓不住"KRAS](#2-三重壁垒)
3. [破局的物理基础：构象呼吸](#3-破局的物理基础构象呼吸)
4. [突破口一：二硫键栓系钓出 Switch‑II 隐匿口袋（G12C 共价）](#4-突破口一)
5. [突破口二：从共价走向非共价（G12D/G12V）](#5-突破口二)
6. [OFF 态药物的软肋与 ON 态难题](#6-off-态药物的软肋)
7. [分子胶三元复合物：化学重塑伴侣蛋白](#7-分子胶三元复合物)
8. [现代药物发现流程：从动态采样到自由能金标准](#8-现代药物发现流程)
9. [可信 vs 不可信的计算判据](#9-可信-vs-不可信的计算判据)
10. [Demo：用 Boltz‑2 复现"两大突破口"](#10-demo)
11. [参考文献](#11-参考文献)

---

## 1. 概述

KRAS 曾被视为肿瘤药理学的"圣杯"，也是"**不可成药（Undruggable）**"靶点的代名词，直到近年通过**变构隐匿口袋（Cryptic Pocket）捕获**与**分子胶三元复合物（Tri‑complex Molecular Glue）**才实现破局。[^stt]

KRAS 是一种**小 G 蛋白分子开关**，在结合 GTP（激活态，ON）与结合 GDP（失活态，OFF）之间循环。癌症突变（主要集中在 **G12、G13、Q61** 位点）破坏了 GTP 酶激活蛋白（GAP）介导的 GTP 水解，使 KRAS 长期锁定在持续激活的 GTP 结合态，驱动 RAF–MEK–ERK 等下游通路持续开火。[^cbm]

| 组织来源 | 高频 KRAS 突变 | 侧链性质 | 代表药物 | 机制 |
|---|---|---|---|---|
| 肺腺癌（NSCLC） | **G12C** | 半胱氨酸（游离 ‑SH，亲核） | Sotorasib, Adagrasib | 共价锁定 OFF 态 |
| 胰腺癌（PDAC）/结直肠癌 | **G12D** | 天冬氨酸（负电羧基） | MRTX1133 | 非共价盐桥锁定 OFF 态 |
| 多种实体瘤 | **G12V** | 缬氨酸（疏水） | （在研） | 非共价 |
| 泛 RAS | 多位点（ON 态） | — | RMC‑6236（Daraxonrasib） | 分子胶三元复合物 |

---

## 2. 三重壁垒

传统药物之所以"抓不住"KRAS，源于三重壁垒：

### 2.1 竞争性小分子壁垒：皮摩尔亲和力 + 毫摩尔底物浓度

激酶抑制剂之所以成功，是因为激酶对 ATP 的亲和力仅在**微摩尔（μM）**级别；而 KRAS 对内源性 GTP/GDP 的亲和力高达 **皮摩尔（pM）** 级别，且细胞内 GTP 浓度高达 **毫摩尔（~0.5 mM）** 级别。[^jci] 任何试图直接去核苷酸口袋"硬抢位置"的小分子，在热力学和浓度梯度上都无法与内源 GTP 竞争。[^stt]

> **量级对比**：pM 亲和力（KRAS·GTP）比 μM（激酶·ATP）高约 10⁶ 倍；再叠加 mM 级胞内 GTP 浓度 → 正交竞争几乎不可能。

### 2.2 变构小分子壁垒："台球"般光滑的表面

除了唯一的核苷酸结合口袋外，静态晶体结构下的 KRAS 近似一个**紧实的球体**，缺乏深邃、疏水的凹槽来提供足够的范德华力与疏水相互作用能（ΔG），导致常规高通量筛选（HTS）命中的小分子"**结合即脱离**"。[^mb]

### 2.3 单克隆抗体壁垒：细胞内定位与空间屏蔽

KRAS 通过 C 端**法尼基化**锚定在细胞质膜内侧，而常规大分子单抗（~150 kDa）**无法穿透**疏水的细胞膜磷脂双分子层。[^cells] 即便针对突变表位设计，也需借助胞内抗体（Intrabodies）、穿膜纳米抗体，或靶向 MHC‑I 呈递新抗原肽段的 TCR‑mimic 抗体——常规胞外抗体对其束手无策。[^bba]

---

## 3. 破局的物理基础：构象呼吸

破局的关键在于打破"蛋白是静态刚体"的思维定势。KRAS 在溶液中存在剧烈的**构象呼吸（Conformational Breathing）**，尤其是控制信号输出的 **Switch‑I** 和 **Switch‑II** 环区具有高度柔性。[^jcim2017] 正是这种柔性，让平时"闭合"的隐匿口袋得以在瞬时打开——只要有分子恰好在它张开的那一刻楔入并稳定它。

---

## 4. 突破口一

### 利用二硫键栓系（Disulfide Tethering）钓出 Switch‑II 隐匿口袋（S‑IIP）

2013 年，UCSF 的 **Kevan Shokat** 团队（Ostrem 等）没有直接去筛非共价强结合剂，而是抓住了肺癌高频突变 **KRAS G12C** 中第 12 位半胱氨酸（Cys12）游离巯基（‑SH）的强亲核反应性。[^nature2013]

**实验逻辑（"钓鱼"策略）**：

1. 构建含二硫键的**分子片段库**（Disulfide Fragment Library）；
2. 在弱还原剂（2‑巯基乙醇）存在下与 GDP 态 KRAS G12C 孵育；[^shokatpdf]
3. 只有当某个片段恰好嵌入 Cys12 附近**瞬时张开**的凹槽、通过非共价作用把自己稳住时，它与 Cys12 形成的二硫键才**不会被还原剂轻易洗脱**（动力学"自筛选"）；
4. 用完整蛋白高分辨质谱（Intact Protein MS）直接检测分子量增加，再解析共晶结构。

**发现**：在 Switch‑II 环区下方，竟被撑开了一个平时闭合的 **Switch‑II 隐匿口袋（S‑IIP）**。[^rg]

**从探针到药物**：随后将二硫键替换为不可逆**共价弹头（丙烯酰胺）**；在优化中（安进研发 **Sotorasib / AMG 510**）发现 **His95 侧链翻转**会进一步露出一个次级隐匿沟槽（**His95 Cryptic Groove**），用芳香环填充后大幅提升亲和力，把 KRAS 死死锁在失活的 **GDP（OFF）态**。[^jho]

> ✅ Demo 1/2 复现了这一突破口：见 §10。

---

## 5. 突破口二

### 从共价走向非共价（攻克无半胱氨酸的 G12D/G12V）

胰腺癌和结直肠癌中最常见的是 **KRAS G12D**（天冬氨酸，带负电羧基）和 **G12V**（缬氨酸，疏水侧链），**没有半胱氨酸**可供共价反应。[^cr]

研发团队（Mirati 研发 **MRTX1133**）直接借用已探明的 Switch‑II 口袋骨架（**吡啶并嘧啶母核**），把共价弹头替换为带正电的**双环二胺**基团（[3.2.1] bicyclic diamino），与 G12D 的 **Asp12 羧基负离子**形成高强度**非共价盐桥（Salt Bridge）**与氢键网络；同时深入挖掘口袋周边的 **Glu62、His95** 等位点，使非共价亲和力直接推升至 **皮摩尔（pM）** 级。这证明了：**只要非共价相互作用足够精准，没有共价键也能撑开并稳住 Switch‑II 裂缝**。[^tlcr]

> ✅ Demo 3/4 复现了这一突破口，并验证了"**非共价精准 ⇒ 突变选择性**"：见 §10。

---

## 6. OFF 态药物的软肋

锁定失活态（OFF）的药物存在一个致命软肋：许多肿瘤细胞中，突变 KRAS 长期处于 **GTP 结合的激活态（ON）**。此时：

- Switch‑II 口袋构象**紧闭**；
- 小分子单独结合 KRAS(ON) 依然**抓不牢**；
- 无法完全阻断下游效应蛋白（如 RAF）的结合。[^science]

这就需要一种全新的范式——不是去核苷酸口袋硬抢，也不是等 Switch‑II 张开，而是**制造一个全新的结合界面**。

---

## 7. 分子胶三元复合物

### 7.1 灵感来源：自然界微生物的"借刀杀人"

新药思路源自天然免疫抑制剂——**环孢素 A（Cyclosporin A）、他克莫司（FK506）、雷帕霉素（Rapamycin）**。这类天然大环内酯本身尺寸有限，但它们先**高亲和力抓取胞内丰度极高的伴侣蛋白**（亲环蛋白 A / **CYPA**，或 **FKBP12**），形成"伴侣蛋白–药物"二元复合物；这个复合体表面形成一片**全新的化学修饰界面**，再去结合原本表面平坦的靶蛋白（钙调神经磷酸酶或 mTOR）。[^nature2024]

### 7.2 化学重塑伴侣蛋白（Chemical Remodeling of a Chaperone）

**Revolution Medicines** 团队（Schulze、Holderfield 等）受此启发，以天然产物 **Sanglifehrin A**（已知高亲和力结合 CYPA 的大环内酯）为锚点彻底改造，设计出：

- **RMC‑6291（Elironrasib）**：靶向 KRAS **G12C(ON)**；[^science]
- **RMC‑6236（Daraxonrasib）**：泛 **RAS(ON)** 非共价抑制剂。[^nature2024]

### 7.3 三元复合物凭什么能稳定结合

```
   小分子(glue)                 CYPA(高丰度伴侣)
        │  ①先高亲和力劫持            │
        └──────►  CYPA : glue  ◄──────┘   二元复合物
                      │
                      │ ② 靠近 KRAS(ON)，诱导新生 PPI 界面
                      ▼
          CYPA : glue : KRAS‑GTP   三元复合物
                      │
                      │ ③ 正协同 α≫1 + 立体位阻盖住效应器面
                      ▼
            k_off↓↓  →  驻留时间↑↑  →  阻断 RAF 结合
```

1. **预先劫持高丰度内源伴侣**：小分子入胞后首先以极高亲和力结合大量存在的 **CYPA**；此时小分子与 KRAS(ON) 单独的二元亲和力其实很弱、甚至测不出。[^science]
2. **诱导新生蛋白–蛋白界面（Neo‑PPI Interface）**：当"CYPA:小分子"靠近 GTP 态 KRAS 时，不仅小分子外露端（Warhead/Linker）接触 Switch‑I/II，更**迫使 CYPA 表面与 KRAS 表面大面积直接接触**，埋藏数百 Å² 的溶剂可及表面积（BSA），形成大量额外的蛋白–蛋白氢键、盐桥与疏水堆积。[^rsc]
3. **巨大热力学正协同（Positive Cooperativity, α≫1）+ 立体位阻封锁**：新生界面提供巨大的额外结合自由能补偿（ΔΔG_coop < 0），使三元复合物的解离速率（k_off）**数量级下降**；同时庞大的 CYPA 像盾牌一样直接盖在 KRAS 的效应器结合面上，通过物理位阻彻底阻断 RAF 等下游蛋白的靠近。[^pdb9bfx]

> 📌 结构实证：**PDB 9BFX**（CYPA : 分子胶 : KRAS 三元复合物）；最初解析的"初始三元"是 **PDB 8G9Q**（CypA : Compound‑1 : KRAS G12C, 1.40 Å, Schulze 2023）。
> ✅ 本节已做计算 demo（Boltz‑2 复现 9BFX + 有/无胶对照）：见 §10.2 与 [`demo/tricomplex/RESULTS.md`](demo/tricomplex/RESULTS.md)。

### 7.4 研发与设计思路（这套"招募伴侣造界面"怎么被设计出来）

- **思想源头**：Gregory Verdine "Drugging the undruggable"——选一个**伴侣蛋白（presenter）**，设计能与它**共同识别靶点**的分子；想法写进专利 **US20160199506A1**（Warp Drive Bio），在 Third Rock Ventures 孵化后成为 Revolution Medicines。
- **为什么是 CypA**：①胞内高丰度广泛表达；②有化学起点（天然产物 **Sanglifehrin A**，对 CypA 亲和力约环孢素 A 的 10–20×）；③与 KRAS 效应面有**局部电荷互补**（FKBP12 不具同样互补）。
- **工具分子 Compound‑1**：双功能（抓 CypA 端 + Cys 反应抓手），用**共价栓系**先锁住弱三元、拿到初始三元晶体（8G9Q）作 SAR 起点——与 2013 年钓 Switch‑II 隐匿口袋是**同一套 tethering 思路**。
- **工具分子→候选药的 4 个杠杆**：①骨架**大环化**（预组织、降熵罚）；②**共价/非共价并行**；③**更换弹头**（丙烯酰胺→炔酰胺 ynamide，三元 IC50 52→28 nM，结合速率 ×26）；④**优化连接段构象**（RMC‑4998，调整弹头姿势，速率再 ×14）。小分子结合后**CypA 表面也诱导重排**，共建识别面。
- **一个平台、三条 RAS 分支**：**elironrasib/RMC‑6291**（G12C，共价）、**daraxonrasib/RMC‑6236**（泛 RAS，非共价）、**zoldonrasib/RMC‑9805**（G12D，共价——用**温和弹头氮丙啶 aziridine** + 界面几何攻克弱亲核的 Asp12，"关键不只在弹头，也在界面"）。

> 📘 完整的**可迁移设计准则 + 对设计新分子的具体启发 + 计算验证方案**整理在 [`design-notes/molecular-glue-design-playbook.md`](design-notes/molecular-glue-design-playbook.md)。

---

## 8. 现代药物发现流程

要从零挖掘隐匿裂缝、并筛选出能稳定二元/三元复合物的药物，现代流程是一套**从动态采样到自由能计算**的严密管线：

### 步骤 1｜增强采样 + 混合溶剂 MD —— 探寻隐匿裂缝

对静态晶体结构运行**微秒级全原子 MD**、**加权系综（Weighted Ensemble）**，或加入苯/乙醇/氙气等疏水探针的**混合溶剂分子动力学（MixMD / Cosolvent MD）**。[^jcim2022] 疏水探针会挤入瞬时张开的缝隙并阻止其闭合；再用**马尔可夫状态模型（MSM）**重构构象自由能地景，识别具有热力学亚稳态特征的**隐匿口袋**。[^jctc_msm]

### 步骤 2｜共价片段栓系 / 生物物理初筛 —— 捕捉微弱苗头

实验端：对候选区域引入定点半胱氨酸突变做**二硫键栓系质谱筛选**；或用 **¹⁹F‑NMR / HSQC‑NMR**（监测残基化学位移扰动）与 **SPR** 寻找微弱结合的片段命中。[^nature2013] 对三元复合物，则在**固定饱和浓度伴侣蛋白（CYPA）**存在下进行配体筛选。[^science]

### 步骤 3｜动力学稳定性 + 水热力学 —— 过滤假阳性姿势

将命中结构置入**显式水盒子**做多重复无偏 MD 与**元动力学（Metadynamics）**：
- 关键氢键/盐桥**占据率 > 70–80%**；
- 口袋残基 **RMSF 显著降低**；
- 用 **WaterMap / GIST** 确认配体是否置换了口袋内的**高能非稳态水（"Unhappy Waters"）**。[^acsomega]

### 步骤 4｜炼金术自由能微扰（FEP+）—— 金标准定量

对通过动力学筛选的分子，构建**热力学循环**运行 **Alchemical FEP+**：[^jctc_fep]
- **二元复合物**：计算相对结合自由能 ΔΔG_bind；
- **分子胶三元复合物**：分别在"**有伴侣蛋白**"与"**无伴侣蛋白**"两个体系中执行配体/蛋白突变的 FEP，求算**协同自由能 ΔΔG_coop**（协同系数 α = exp(−ΔΔG_coop/RT)），筛选 **ΔΔG_coop 显著为负**的强正协同分子。[^jcim_coop]

### 步骤 5｜生物物理动力学（k_off）+ 高分辨共晶 —— 最终闭环

用 **SPR / BLI** 测定解离速率常数 k_off（算驻留时间 τ = 1/k_off）；用**氢氘交换质谱（HDX‑MS）**验证 Switch 区被锁定；用 **X 射线晶体学 / 冷冻电镜** 解析高分辨结构完成闭环。[^natcomm]

```
静态晶体 → ①MixMD/WE/MSM 找裂缝 → ②栓系/NMR/SPR 初筛
        → ③显式水 MD + WaterMap 过滤 → ④FEP+ 定量(ΔΔG / α)
        → ⑤SPR·k_off + HDX‑MS + Cryo‑EM 闭环
```

---

## 9. 可信 vs 不可信的计算判据

在评估像 KRAS 这种**高柔性隐匿口袋**或**分子胶三元复合物**是否真正"稳定结合"时，计算方法可靠性存在严格分水岭：

| 评估维度 | ❌ 不可信（高假阳性误区） | ✅ 可信（真实稳定性准则） |
|---|---|---|
| **初始结构与采样** | 单一静态 Apo 晶体**刚性对接**：Switch‑II 处于闭合态，直接对接会漏掉**所有**活性分子 | **系综对接 + 诱导契合/混合溶剂采样**：先用 MixMD / WE / MSM 采出口袋张开的亚稳态构象簇，再做柔性诱导契合 |
| **打分 / 亲和力** | 仅凭**对接打分**（Vina/Glide）或**单帧 MM‑GBSA**：忽略构象熵损失与显式水效应，与实验相关性极差 | **炼金术 FEP / TI**：基于完整统计力学配分函数与热力学循环，显式含溶剂重组与熵变，误差通常可控在 ~1 kcal/mol |
| **轨迹稳定性判据** | 仅看短时（10–50 ns）**蛋白整体 RMSD 平稳**：核心 β 折叠本就稳，掩盖不了配体在浅口袋里的滑移 | **配体相对口袋重原子 RMSD + 关键作用存活率 + 解离路径采样**：多条 μs 级独立轨迹，配体 RMSD 稳定、盐桥/氢键占据率高，用伞形采样/tMD 估解离能垒 |
| **三元复合物（分子胶）** | **刚性蛋白–蛋白对接**，或只算小分子对单蛋白亲和力：忽视界面构象重排与去水合代价 | **热力学循环协同自由能 ΔΔG_coop + 界面 BSA 分布**：双体系 FEP 算三元相对二元的自由能增益，结合新生 PPI 界面的埋藏面积与静电互补性判稳 |

**一句话总结**：*可信的结论来自"动态采样 + 显式溶剂 + 自由能"三件套；单一静态结构的刚性对接打分，几乎必然高假阳性。*

---

## 10. Demo

> 完整说明、脚本与结果见 [`demo/`](demo/)。§10.1 对应突破口一/二（二元复合物），§10.2 对应 §7 分子胶三元复合物。

### 10.1 二元复合物：突破口一/二

> 下表为运行配置；结果见 [`demo/RESULTS.md`](demo/RESULTS.md)。

我们用 **Boltz‑2**（结构 + 亲和力联合预测）复现"两大突破口"的核心主张，并与 **ChEMBL 实验数据**对照。构建体统一采用人源 KRAS 催化域（残基 1–169，UniProt P01116），在第 12 位点引入对应突变。

| # | 体系 | 突破口 / 目的 | 药物（来源） | 实验锚点（ChEMBL） |
|---|---|---|---|---|
| 1 | KRAS **G12C** + **Sotorasib**(AMG510) | 突破口一：S‑IIP 共价 | 已上市，CHEMBL4535757 | KRAS IC50 = 68 nM；G12C 细胞 5 nM |
| 2 | KRAS **G12C** + **Adagrasib**(MRTX849) | 突破口一：交叉验证 | 已上市，CHEMBL4594350 | 第二个获批 G12C 药 |
| 3 | KRAS **G12D** + **MRTX1133** | 突破口二：非共价盐桥 | 临床，CHEMBL4858364 | G12D(1–169) IC50 = **0.4 nM** |
| 4 | KRAS **WT** + **MRTX1133** | 对照：非共价精准 ⇒ 选择性 | 同上 | WT 细胞 IC50 **> 3000 nM** |

**验证逻辑**：
- Demo 1/2 应将配体准确放入 **Switch‑II 隐匿口袋**（而非核苷酸口袋），佐证"隐匿口袋可被小分子占据"；
- Demo 3 应在 **Asp12 ↔ 双环二胺**间形成盐桥，佐证"非共价也能撑住裂缝"；
- Demo 3 vs 4（G12D 强、WT 弱）应复现 **MRTX1133 的突变选择性**。

> ⚠️ **诚实边界**：Boltz‑2 是单结构 + 亲和力的深度学习预测，**不等于**第 8 节的 MD/FEP 金标准流程。它适合做快速假设生成与姿势合理性检查，但对"稳定结合"的最终判定仍需 §8 的动态采样 + 自由能闭环。本 demo 的定位是**用可获得的工具快速演示文中思路**，而非给出监管级结论。

### 10.2 分子胶三元复合物：CYPA : 胶 : KRAS，对照 PDB 9BFX

复现 §7：分子胶诱导的 **CYPA–KRAS 新生界面（neo‑PPI）**。参照 **PDB 9BFX**（Elironrasib/RMC‑6291 + KRAS G12C + CypA，ON 态）。
三个体系：**T1 非共价**（CYPA+KRAS+RMC‑6291+GppNHp+Mg）、**T1cov**（同 T1 但显式声明 Cys12 SG↔胶弹头共价键）、**T2 无胶对照**（CYPA+KRAS）。详见 [`demo/tricomplex/RESULTS.md`](demo/tricomplex/RESULTS.md)。

| 指标 | 9BFX 实验 | T1 非共价 | T1cov 共价 | T2 无胶 |
|---|---|---|---|---|
| binding_confidence | — | 0.903 | **0.931** | **0.317** |
| CYPA–KRAS 界面置信 protein_iptm | — | 0.819 | 0.754 | 0.21–0.68（波动） |
| 胶↔Cys12 (Å) | 1.68（共价） | 4.66 | **1.96** | — |
| 蛋白–蛋白 BSA (Å²，近似) | ~425 | ~337 | **~417** | ~580\* |
| **KRAS 界面残基 vs 9BFX** | — | 100% / 70% | 89% / 80% | 28%（非天然） |
| **CYPA 界面残基 vs 9BFX** | — | 67% / 55% | **90% / 82%** | 64% / 64% |

**结论**：含胶时 Boltz **高置信、可复现地复现了 9BFX 的天然 CYPA–KRAS 界面**；去掉胶后置信度崩塌、界面变成非天然硬对接——定性印证"**分子胶诱导新生界面 + 正协同**"。进一步**显式声明 Cys12 共价键（T1cov）**后，胶↔Cys12 收紧到 1.96 Å、BSA 回到 ~417 Å²（≈实验）、CYPA 界面回收率 55%→**82%**——约束越贴合真实化学，预测越贴近晶体。
⚠️ 这仍是协同性的**定性代理**，非 §8 的 ΔΔG_coop 金标准定量。

---

## 11. 参考文献

[^stt]: Signal Transduction and Targeted Therapy — <https://pmc.ncbi.nlm.nih.gov/articles/PMC8591115/>
[^cbm]: Cancer Biology & Medicine — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12302276/>
[^jci]: Journal of Clinical Investigation — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11789494/>
[^mb]: Molecular Biomedicine — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12462266/>
[^cells]: Cells — <https://pmc.ncbi.nlm.nih.gov/articles/PMC9954350/>
[^bba]: Biochimica et Biophysica Acta — <https://pmc.ncbi.nlm.nih.gov/articles/PMC9627686/>
[^jcim2017]: Journal of Chemical Information and Modeling — <https://pubs.acs.org/doi/10.1021/acs.jcim.7b00399>
[^nature2013]: Ostrem et al., Nature 2013 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC4274051/>
[^shokatpdf]: Shokat Lab — <https://shokatlab.ucsf.edu/pdfs/29320178.pdf>
[^rg]: Discovery of switch‑II pocket binder (figure) — <https://www.researchgate.net/figure/Discovery-of-switch-II-pocket-binder-a-Chemical-structures-of-the-screening-hits-2E07_fig2_341805423>
[^jho]: Journal of Hematology & Oncology — <https://pmc.ncbi.nlm.nih.gov/articles/PMC8767686/>
[^cr]: Cancer Research — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10922474/>
[^tlcr]: Translational Lung Cancer Research — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11257988/>
[^science]: RMC‑6291 / tri‑complex, Science — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10474815/>
[^nature2024]: RMC‑6236 pan‑RAS(ON), Nature — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11111408/>
[^rsc]: RSC Chemical Biology — <https://pmc.ncbi.nlm.nih.gov/articles/PMC9994104/>
[^pdb9bfx]: PDB 9BFX — <https://www.rcsb.org/structure/9BFX>
[^jcim2022]: Cosolvent / MixMD, JCIM — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11558672/>
[^jctc_msm]: Markov State Models, JCTC — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10373493/>
[^acsomega]: WaterMap / GIST, ACS Omega — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12959236/>
[^jctc_fep]: Alchemical FEP+, JCTC — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12159975/>
[^jcim_coop]: Cooperativity FEP, JCIM — <https://pmc.ncbi.nlm.nih.gov/articles/PMC9869335/>
[^natcomm]: Nature Communications — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10344917/>

---

*本文档为技术学习笔记整理，引用均来自公开文献；KRAS 构象呼吸、隐匿口袋与分子胶机制的定量结论请以原始论文与监管级实验/计算为准。*
