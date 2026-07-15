# G034 — IL‑14α 在干燥综合征中的机制解析 + 多模态药物设计

> 承接 `REPORT.md`。本篇做两件事:(1) 用**蛋白互作 / 网络不可约**视角 + omics 证据,把 IL‑14α 在干燥综合征(SjD)中的**机制网络**拆到"不可约核心";(2) 据此给出**siRNA(已实跑候选)/ 抗体 / 降解剂**的多模态药物设计,并标明 **Evo2 / AlphaGenome** 在其中的**具体接入点**。
> 诚信标注:✅ 已证据支撑 / ⚠️ 推断或有争议 / ❗ 风险 / 🧪 本会话已实跑 / 🔌 需外部模型(本会话未接入)。

---

## 第一部分 — 机制网络解析

### 1.1 证据驱动的因果网络(✅ 主干有因果实验支撑)

```mermaid
flowchart TD
  subgraph TRIGGER[上游触发 ⚠️]
    ENV[环境/病毒/遗传易感] --> IL14[IL‑14α↑<br/>血清可检, pSS 升高 ✅]
  end
  IL14 -->|B 细胞生长因子活性 ✅| MZB[边缘区 B 细胞 MZB 扩增/活化<br/>★不可约核心节点 ✅]
  BAFF[BAFF↑ ✅ 与 IL‑14α 并列升高] --> MZB
  MZB -->|生发中心形成| GC[异位生发中心 GC<br/>唾液/泪腺内 ✅]
  GC --> AAB[自身抗体<br/>抗 SSA/Ro, SSB/La ✅ 与血清 IL‑14α 相关]
  GC --> PLASMA[浆细胞 / 高丙球血症 ✅]
  MZB --> LT[淋巴毒素 LTαβ ✅ 抗‑LT 可干预]
  LT --> STROMA[腺体基质/HEV 重塑]
  AAB --> DESTRUCT[腺体上皮损伤<br/>分泌功能↓ → 干眼/口干 ✅]
  STROMA --> DESTRUCT
  DESTRUCT --> SICCA[临床干燥症]
  MZB -->|长期| LYMPH[CD5+ B 细胞淋巴瘤<br/>晚期 ✅]

  classDef core fill:#b30000,color:#fff,stroke:#600;
  class MZB core;
```

### 1.2 "不可约核心"识别(用你的网络/不可约方法)

把上图当作一个调控网络,做**节点删除 → 表型坍缩**的因果测试(等价于网络约化里找"驱动核心"):

| 删除的节点/边 | 实验结果 | 网络推论 |
|---|---|---|
| **IL‑14α 基因**(KO) | 转基因模型全部 SjD 表型不出现 ✅ | IL‑14α 是**必要上游驱动** |
| **MZB**(B 细胞特异敲 RBP‑J) | **全部 SjD 表型消失**(分泌正常、自身抗体阴性、组织学正常)✅ | **MZB = 不可约核心节点**;IL‑14α 的致病性**完全经 MZB 传导** |
| B1 细胞(敲 btk) | 无改善 ✅ | B1 不在核心通路 |
| 淋巴毒素(抗‑LT mAb) | 部分阶段可干预 ✅ | LT 是**下游可药节点**(非核心) |

> **结论(✅):** 网络的**不可约驱动边 = IL‑14α → MZB**。这条边是"最小充分致病模块":上游删 IL‑14α 或下游删 MZB 都能让系统从"致病吸引子"回到"健康吸引子"。**这为把 IL‑14α 作为一线干预点提供了网络层面的因果依据**——不是相关性,而是双向敲除验证。

### 1.3 动力学解读:双稳态与"临界点"生物标志物(⚠️ 假说,可检验)

SjD 的分期(高丙球→自身抗体→腺体浸润→淋巴瘤)呈**不可逆推进**,符合**双稳态 + 临界跃迁(tipping point)**结构:IL‑14α/MZB 正反馈(MZB 产更多自身抗体、GC 正反馈)构成 bistable switch。可检验的**早预警生物标志物(EWS)**:血清 IL‑14α 与 BAFF 的**协同上升 + 波动方差增大**(临界慢化)先于腺体功能下降。这正好可复用你仓库里的 `network-biomarker`(临界慢化 / DNB)流水线来定量——**把 IL‑14α–BAFF–autoAb 三节点当作 DNB 模块**,在纵向患者血清上算 SD/自相关上升。🔌 需纵向队列数据。

### 1.4 与 omics 的衔接(✅/⚠️)

- **转录组(2025 唾液腺,IL‑14αTG):** 过表达 IL‑14α → B 细胞增殖、生发中心基因、免疫浸润签名上调(✅ 定性)。可据此取**下游签名基因**作为 siRNA/抗体的 PD 读出(药效标志物)。
- **蛋白组:** 血清 IL‑14α(WB)与 BAFF(ELISA)可作**靶点结合/PD 双标**。
- **单细胞(建议补做,🔌):** IL‑14αTG vs KO 腺体 scRNA/CITE‑seq → 确认 IL‑14α 的**受体在哪群 B 细胞**(闭合"IL‑14R 分子身份未知"的缺口,直接服务抗体 MoA)。

---

## 第二部分 — 多模态药物设计

设计原则(承接 REPORT §4/§6/§7):**打不可约边 IL‑14α→MZB**。按靶点可及性分两条主线,加一条"安全性差异化"暗线。

### 2.1 模态 A — siRNA / ASO(胞内首选;🧪 本会话已给候选)

**为什么 siRNA 是胞内首选:** α‑taxilin 无口袋、抗体进不去;敲低直接清除蛋白,且不必纠结"胞内 vs 分泌"哪个池致病——一并降掉。

**🧪 已实跑的理性设计**(`sirna_design.py` / `sirna_candidates.txt`):对提供的 TXLNA 天然 CDS(编码 aa210–546 = IL‑14α C 区)扫描 993 个 19‑mer 窗口,施加 Reynolds + Ui‑Tei + Amarzguioui + 热力学不对称 + liability 过滤。**Top‑2 满分(9/9)候选:**

| 排名 | CDS 位点 | GC% | guide(反义 5′→3′) | sense/靶 19‑mer | 备注 |
|---|---|---|---|---|---|
| 1 | 956 | 47% | `UUUGUCGAUAUGCUCCUCG`(+dTdT) | CGAGGAGCATATCGACAAA | 满分、无 liability |
| 2 | 929 | 42% | `AUACUGCUCAAUCAGCUUC`(+dTdT) | GAAGCTGATTGAGCAGTAT | 满分、无 liability |

(另有 10 条 8/9 的干净候选,见 `sirna_candidates.txt`。)

**化学与递送(成药关键):**
- 化学:2′‑OMe/2′‑F 交替、末端硫代磷酸(PS)、guide 5′‑(E)‑乙烯基膦酸酯 → 稳定性 + 降低免疫刺激。
- ❗**递送是真正难点:** GalNAc 只到肝;靶细胞是**淋巴组织/腺体的 MZB**。方案:①脂质纳米粒(LNP)向脾/淋巴归巢;②**抗体‑siRNA 偶联(AOC)**——用抗 B 细胞表面标志(如 CD19/CD22)的抗体把 siRNA 靶向递送到 B 细胞(与本项目抗体能力协同,收敛点)。
- ❗**专利差异化(vs US7622574 "IL‑14α RNA inhibitors"):** 该基础专利宽泛覆盖"抗 IL‑14α 的 RNA 抑制剂"。差异化靠:**具体 siRNA 序列的物质专利** + **新化学修饰** + **新递送(AOC)** + 下面 2.4 的**剪接选择性 MoA**。

**🔌 Evo2 / AlphaGenome 接入点(本会话未接入,设计已就位):**
- **Evo2(Arc,DNA 基础模型):** ①对每条 guide 做**全基因组种子脱靶似然**打分(超越 BLAST 的上下文感知);②**变体效应**——扫描靶位点常见 SNP(人群差异)是否削弱敲低 → guide‑患者匹配;③在 exon 3–10(IL‑14α 特异区)内优先选点。
- **AlphaGenome(DeepMind,长上下文调控/剪接预测):** 见 2.4——预测**剪接切换 ASO**的外显子跳跃后果,设计出**只降 IL‑14α 转录本、保留看家 taxilin 运输功能**的分子。

### 2.2 模态 B — 抗体(若 Gate‑1 证实存在可溶 IL‑14α)

- 目标:中和 **IL‑14α ↔ B 细胞 IL‑14R**(不可约边的胞外段)。
- 形式:人源化/全人**全长 IgG**(甲方要求);Fc 选**效应减弱**(IgG4 或 IgG1‑LALA),因为目的是"中和可溶配体"而非杀细胞。
- **表位策略(用 Boltz 结构指导,见 2.5):** 避开 US7622574/Peng2009 用过的 C 端表位,优先 coiled‑coil 上与受体/syntaxin 竞争的构象表位 → 新颖性 + 功能中和。
- 发现路径:全人噬菌体/酵母文库(绕开自身耐受)或 IL‑14α‑KO 小鼠免疫。

### 2.3 模态 C — 降解剂(若纯胞内)

- **PROTAC/分子胶**机制契合胞内 α‑taxilin,但❗零配体。抓手:**4 个 Cys(C245/C373/C471/C523)做共价片段筛选** → 共价 PROTAC。
- 若存在分泌/表面池 → **LYTAC / 抗体基降解剂**复用 2.2 的抗体。

### 2.4 ⭐ 差异化暗线 — 剪接选择性下调(安全性 + 专利双赢,🔌 AlphaGenome 驱动)

**问题:** α‑taxilin 是**泛表达看家运输蛋白**;无差别敲低有 on‑target 毒性风险(REPORT §5.4)。
**思路:** "IL‑14α" 由 il14 基因 **plus 链外显子 3–10** 编码,与经典 taxilin 的转录/功能存在差异。**设计剪接切换 ASO(SSO)**,只干扰产生致病 IL‑14α 物种的剪接/翻译起始,**尽量保留看家 taxilin 运输功能** → 降低毒性,并构成**全新 MoA 专利**(区别于泛敲低)。
**AlphaGenome 角色:** 预测候选 ASO 对 TXLNA/IL14 位点**剪接结果与转录本丰度**的影响,筛出"降 IL‑14α、留 taxilin"的最优靶点;Evo2 复核序列/脱靶。🔌 需接入这两模型执行。

### 2.5 结构计算计划(Boltz‑2.1,本会话可跑,需你确认付费)

| 任务 | 目的 | 状态 |
|---|---|---|
| 折叠抗原(337 aa)单链 | 定位表面表位 + 4 个 Cys 可及性(服务抗体表位/共价降解剂) | ⏸ 估算调用遇服务瞬断;待重试 |
| α‑taxilin C 区 + syntaxin‑4 SNARE 螺旋复合物 | 定位功能 PPI 界面 = "阻断"表位 | 待跑 |
| 抗 IL‑14α **nanobody de‑novo 设计** ×50(boltz_curated) | 产出候选 VHH 序列 | 🧪 运行中(prot_des_uXUu4UgocY7hyIIUZYO6,~$2.50) |

> EDEN 免疫原性预测因抗原是"片段 CDS"(分布外)不适用,已如实弃用。

#### 🧪 Boltz‑2.1 已完成结果(2026‑07‑15)

**① 抗原折叠(sab_pred_1Ad3Tn…):** structure_confidence 0.68,但 **pTM 0.28**(全局拓扑不定)。按区域拆:
- **卷曲螺旋核心 aa 210–484(query 1–275):平均 pLDDT 84.9 → 置信折叠的 α‑螺旋杆**;
- **C 端富脯氨酸尾 aa 485–546(query 276–337):平均 pLDDT 48.1 → 无序**。
- ✅ 定量证实了 REPORT 的判断:**长 coiled‑coil 杆 + 无序 PRD,无深口袋**(再次支持"小分子难成药、PROTAC 需先找共价配体")。

**② 4 个 Cys 的溶剂暴露(邻居数越低越暴露,中位=12):**
| Cys(P40222) | 区域 | pLDDT | 暴露 | 作为共价 PROTAC 抓手的适用性 |
|---|---|---|---|---|
| **C245** | 卷曲螺旋核心 | 89.8 | 中等 | ⭐ **优选**(结构化 + 半暴露) |
| **C373** | 卷曲螺旋核心 | 88.4 | 中等 | ⭐ **优选**(结构化 + 半暴露) |
| C471 | 卷曲螺旋核心 | 68.6 | 埋藏 | 次选(埋藏,难接触) |
| C523 | 无序 PRD | 48.3 | 最暴露 | ❗虽最暴露但在**柔性无序尾**,非稳定口袋,不利共价锚定 |

> **修正设计建议:** 共价片段筛选**优先 C245 / C373**(位于置信折叠的螺旋核心、半暴露),**而非**先前按序列位置推测的 C523。**抗体表位**应取**结构化卷曲螺旋核心(210–484)**的构象表位——既避开先前技术(Peng2009 / US7622574 用的 C 端多肽),又落在高置信折叠区。

**③ 抗原–syntaxin‑4 复合物(sab_pred_12fovJS…):** iPTM 0.40、**binding_confidence ≈ 1.5×10⁻⁵、界面 PDE ≈ 16 Å**。
- ⚠️**低置信/阴性——单样本、template‑free、含 Habc 自抑制域的 STX4 未能置信定位 taxilin–syntaxin‑4 界面。** 不能据此确定"阻断表位"。需:裁剪到 STX4 的 H3/SNARE 螺旋、增采样、或加 MSA/模板后重跑。如实记录,不夸大。

> 结构文件存于 `boltz_structures/antigen_fold.cif`(可用 PyMOL/ChimeraX 打开)。

**④ 🧪 de‑novo 纳米抗体设计 ×50(prot_des_uXUu4U…,已完成,$2.50):**
- 50 条中 **1 条显著胜出**(`pres_QA51n…`):**iPTM 0.828、pTM 0.984、complex‑pLDDT 0.817、界面 PAE 3.67 Å**;其余 iPTM<0.57、界面 PAE>10 Å(见 `design_ranking.tsv`)。→ 说明**对柔性 coiled‑coil 做无约束 de‑novo 设计命中率低**(符合预期),但**确实拿到 1 个高置信结合物**。
- **Top VHH 序列**(116 aa,标准 VHH 骨架,存 `nanobody_top_design.faa`):
  `EVQLVESGGGLVQPGGSLRLSCAASGFTFSKHSMHWVRQAPGKGLEWVSSISSDGSLVLSAPSVAGRFTISRDNAKNTLYLQMNSLRPEDTAVYYCAREALNPTLRGQGTLVTVSS`(CDR3≈`AREALNPTLRG`)。
- **表位映射(结合物 5 Å 内抗原残基):P40222 325–350** —— ✅ 落在**结构化卷曲螺旋核心**(非无序尾),且**与先前技术(Peng2009/US7622574 的 C 端多肽表位)不同** → 新颖表位、利于专利。该区在 syntaxin 结合结构域内,理论上有阻断潜力。
- ❗**局限(如实):** 纯 in‑silico 设计,需**湿实验验证(表达/亲和/中和)**;VHH 是骨架、与甲方"人源化全长 IgG"要求不同(可作先导/工具,或 CDR 移植到人源 IgG);靶点为自身抗原,需查 humanness/免疫原性;binding_confidence 指标为 0(与 iPTM 不一致,说明该指标对此类 de‑novo 结果校准不佳,应以 iPTM/iPAE 为准)。
- **下一步(可选,付费):** ①用 top VHH 对**折叠好的抗原**做 epitope‑directed 复算/亲和成熟;②裁剪 STX4 到 H3 螺旋重跑复合物以定位阻断界面;③把 CDR 移植到人源 IgG 骨架 + Boltz 复算。

**⑤ 🧪 STX4‑H3 复合物重跑(sab_pred_SmY0…,已完成,5 samples):** 裁剪到 syntaxin‑4 H3/SNARE 螺旋(去 Habc)后**显著改善:iPTM 0.63、complex‑pLDDT 0.82、iPDE 8.25 Å**(vs 全长 STX4 的 iPTM 0.40 / iPDE 16 / binding_conf 1.5e‑5;新 binding_conf 7.4e‑4,↑50×)。
- **✅ Axis‑2 阻断表位(syntaxin‑4 结合面)= P40222 355–443**(延伸的 coiled‑coil 界面)+ C 端尾 535–546 接触。结构文件 `boltz_structures/antigen_stx4H3_complex.cif`。
- **两个关键收敛(强):**
  1. 我们的 **Axis‑2 抗体设计表位 360–389 正落在验证过的 syntaxin 界面(355–443)之内** → 该批抗体理应能真正阻断 syntaxin 结合(表位选择被结构验证)。
  2. **共价 PROTAC 靶点 C373 恰位于 syntaxin 界面上** → 在 C373 上的共价配体可**同时**(a)锚定降解剂、(b)直接竞争阻断 syntaxin 结合 → 降解 + 阻断双机制收敛于同一位点。
- ⚠️ 绝对 binding_confidence 仍低(coiled‑coil×coiled‑coil、无实验模板);结论为"界面定位可信、置信中等",需实验(突变扫描/交联‑MS)确认。

---

## 决策流(把机制→模态串起来)

```mermaid
flowchart TD
  G1{Gate‑1: 可溶 IL‑14α 存在且可中和?}
  G1 -->|是| AB[主推 抗体 中和 IL‑14α–IL‑14R<br/>+ siRNA 作后备]
  G1 -->|否/纯胞内| SI[主推 siRNA/ASO 敲低<br/>+ 降解剂需配体campaign]
  AB --> SAFE{泛表达毒性?}
  SI --> SAFE
  SAFE -->|担心| SSO[剪接选择性 ASO<br/>AlphaGenome 设计 → 只降 IL‑14α]
  SAFE -->|可控| GO[推进先导]
```

---

## 交付物清单(本篇)
- `MECHANISM_AND_DRUG_DESIGN.md` — 本文件
- `sirna_design.py` — 🧪 可复现理性 siRNA 设计脚本(无外部依赖)
- `sirna_candidates.txt` — 🧪 Top 候选输出

## 诚信小结
- ✅ 机制主干(IL‑14α→MZB 不可约边)有**双向敲除因果证据**。
- 🧪 siRNA 候选是**真实跑出**的(理性规则),但**基因组级脱靶未做**(需 Evo2/BLAST)。
- 🔌 **Evo2 与 AlphaGenome 在本会话无可调用工具**——我给出的是**已就位的接入设计**,不是已执行结果;一旦连上即可跑。
- ⏸ Boltz 结构/设计为**付费**,先 estimate 再经你确认。
