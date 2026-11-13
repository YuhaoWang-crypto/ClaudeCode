# 总报告：从 KRAS "不可成药"到隐匿口袋、分子胶与可计算的分子胶设计

> 一份把**科学背景 + 计算 demo + 可迁移设计方法**串起来的总报告。
> 全部计算用 **Boltz‑2.1**（结构+亲和力/共价/三元）+ 自写的纯 numpy 结构分析器完成；实验锚点取自 **ChEMBL v34** 与 **RCSB PDB**。
> 📏 **诚实标注贯穿全文**：✅=基于真实结构/实验的可信结论；⚠️=深度学习预测的假设生成，非金标准定量。
> 📂 每节都链接到可复现的详细文件与脚本。运行日期 2026‑11。

---

## 目录

- [第一部分 · 科学背景：KRAS 为何曾"不可成药"](#第一部分--科学背景kras-为何曾不可成药)
  - [1.1 三重壁垒](#11-三重壁垒)
  - [1.2 破局物理基础：构象呼吸](#12-破局物理基础构象呼吸)
  - [1.3 突破口一：G12C 共价隐匿口袋](#13-突破口一g12c-共价隐匿口袋)
  - [1.4 突破口二：G12D 非共价盐桥](#14-突破口二g12d-非共价盐桥)
  - [1.5 ON 态软肋与分子胶三元复合物](#15-on-态软肋与分子胶三元复合物)
  - [1.6 现代发现流程与可信/不可信判据](#16-现代发现流程与可信不可信判据)
- [第二部分 · 计算 Demo：复现两大突破口与分子胶](#第二部分--计算-demo复现两大突破口与分子胶)
  - [2.1 二元复合物：共价 G12C / 非共价 G12D / 选择性对照](#21-二元复合物共价-g12c--非共价-g12d--选择性对照)
  - [2.2 分子胶三元复合物：对照 PDB 9BFX + 显式共价键](#22-分子胶三元复合物对照-pdb-9bfx--显式共价键)
- [第三部分 · 设计 Playbook：从案例提炼的可迁移准则](#第三部分--设计-playbook从案例提炼的可迁移准则)
- [第四部分 · 设计 Pilots：把设计决策变成可计算判据](#第四部分--设计-pilots把设计决策变成可计算判据)
  - [Pilot 1 · presenter 选择（CypA vs FKBP12）](#pilot-1--presenter-选择cypa-vs-fkbp12)
  - [Pilot-flex · presenter 是刚是柔？](#pilot-flex--presenter-是刚是柔)
  - [Pilot 2 · 平台外推 STAT3：判定走哪条分支](#pilot-2--平台外推-stat3判定走哪条分支)
  - [Pilot 3 · STAT3 降解分支：CRBN 诱导邻近](#pilot-3--stat3-降解分支crbn-诱导邻近)
- [第五部分 · 总结果表 + 设计决策流程](#第五部分--总结果表--设计决策流程)
- [第六部分 · 方法学诚实边界总览](#第六部分--方法学诚实边界总览)
- [第七部分 · 文件索引与复现](#第七部分--文件索引与复现)
- [参考文献](#参考文献)

> 📘 更详细的分节原文：[`README.md`](README.md)（科学背景全文）、[`demo/RESULTS.md`](demo/RESULTS.md)、[`demo/tricomplex/RESULTS.md`](demo/tricomplex/RESULTS.md)、[`design-notes/molecular-glue-design-playbook.md`](design-notes/molecular-glue-design-playbook.md) 及各 pilot 的 RESULTS。

---

## 第一部分 · 科学背景：KRAS 为何曾"不可成药"

KRAS 是小 G 蛋白分子开关，在 GTP(ON)/GDP(OFF) 间循环；癌症突变（G12/G13/Q61）破坏 GAP 介导的 GTP 水解，使其锁死在激活态。它曾是"不可成药"的代名词，直到**变构隐匿口袋**与**分子胶三元复合物**破局。

### 1.1 三重壁垒
1. **竞争性壁垒**：KRAS 对内源 GTP/GDP 亲和力达 **皮摩尔(pM)**，胞内 GTP 高达 **~0.5 mM**——小分子无法去核苷酸口袋硬抢（激酶抑制剂成功是因为 ATP 亲和力才 μM）。
2. **变构壁垒**：静态晶体下 KRAS 表面像"台球"般光滑，缺乏深疏水凹槽，HTS 命中"结合即脱离"。
3. **抗体壁垒**：KRAS 经法尼基化锚在质膜内侧，常规 150 kDa 单抗无法穿膜。

### 1.2 破局物理基础：构象呼吸
KRAS 在溶液中存在剧烈**构象呼吸**，Switch‑I/II 环区高度柔性——平时闭合的隐匿口袋会瞬时张开。

### 1.3 突破口一：G12C 共价隐匿口袋
2013 Shokat 组用**二硫键栓系**抓住 G12C 的 Cys12 游离巯基，钓出 Switch‑II 下方闭合的 **S‑IIP 隐匿口袋**；换不可逆弹头（丙烯酰胺）+ 挖 His95 次级沟槽 → **Sotorasib/AMG510**，把 KRAS 锁死在 OFF 态。

### 1.4 突破口二：G12D 非共价盐桥
G12D/G12V 无半胱氨酸。Mirati 借 Switch‑II 骨架，用带正电**双环二胺**与 Asp12 羧基形成**非共价盐桥** → **MRTX1133**，亲和力推到 pM 级，证明"非共价够精准也能撑住裂缝"。

### 1.5 ON 态软肋与分子胶三元复合物
OFF 态药物的软肋：很多肿瘤里突变 KRAS 长期处于 ON 态、Switch‑II 紧闭。**Revolution Medicines** 受天然免疫抑制剂（环孢素/FK506/雷帕霉素借伴侣蛋白）启发，以 **Sanglifehrin A** 为锚改造出**招募亲环蛋白 CYPA** 的分子胶：先劫持高丰度 CYPA，再诱导一个原本不存在的 **CYPA–KRAS 新生界面(neo‑PPI)**，用正协同 + 立体位阻把 KRAS(ON) 钉住。平台三分支：**RMC‑6291(G12C,共价) / RMC‑6236(泛RAS,非共价) / RMC‑9805(G12D,共价)**。

### 1.6 现代发现流程与可信/不可信判据
5 步：①增强采样/混合溶剂 MD 找隐匿裂缝 → ②栓系/NMR/SPR 初筛 → ③显式水 MD + WaterMap 过滤 → ④**FEP 自由能金标准**（二元 ΔΔG / 三元协同 ΔΔG_coop）→ ⑤SPR k_off + HDX‑MS + 晶体/冷冻电镜闭环。
**核心判据**：可信结论来自"动态采样 + 显式溶剂 + 自由能"；单一静态刚性对接打分几乎必然高假阳性。（全文见 [`README.md`](README.md) §8–§9。）

---

## 第二部分 · 计算 Demo：复现两大突破口与分子胶

> 统一构建体：人源 KRAS 催化域 1–169（UniProt P01116），第 12 位引入对应突变。详见 [`demo/RESULTS.md`](demo/RESULTS.md)、[`demo/tricomplex/RESULTS.md`](demo/tricomplex/RESULTS.md)。

### 2.1 二元复合物：共价 G12C / 非共价 G12D / 选择性对照

4 个体系 Boltz 结构+亲和力预测 + 自写口袋接触分析（[`demo/analyze_pocket.py`](demo/analyze_pocket.py)）：

| # | 体系 | 药物(ChEMBL) | 结构验证关键发现 | binding_confidence |
|---|---|---|---|---|
| 1 | G12C+**Sotorasib** | 已上市 | 弹头落 **Cys12 SG 2.75 Å**，嵌 Switch‑II(8残基)+His95，核苷酸口袋零接触 | 0.925 |
| 2 | G12C+**Adagrasib** | 已上市 | Cys12 SG 2.73 Å，交叉印证 | 0.971 |
| 3 | G12D+**MRTX1133** | 临床 | **Asp12 羧基↔配体N = 3.0 Å 盐桥**，Switch‑II 9残基+His95 | 0.995 |
| 4 | WT+MRTX1133（对照） | — | 同口袋但第12位 **Gly12（4.14 Å，无盐桥）** | 0.988 |

**✅ 三条可直接对应机制的发现：**(a) 配体占据 **Switch‑II 隐匿口袋**而非核苷酸口袋（绕开竞争壁垒）；(b) 共价弹头精准对准 Cys12 并接触 His95；(c) **Asp12 盐桥是 G12D 选择性的结构根源**（WT Gly12 无此桥，实验 G12D 0.4 nM vs WT 细胞 >3000 nM）。
**⚠️ 诚实点**：Boltz 的标量 optimization_score 几乎没区分 G12D(0.76)/WT(0.70)，但**相互作用层面的结构分析**清楚揭示了盐桥有无——印证"不能只看一个分数"（§1.6）。

### 2.2 分子胶三元复合物：对照 PDB 9BFX + 显式共价键

对照 **PDB 9BFX**（Elironrasib/RMC‑6291 + KRAS G12C + CypA, ON 态）。自写界面/BSA/桥接分析（[`demo/tricomplex/tri_analyze.py`](demo/tricomplex/tri_analyze.py)，含纯 numpy Shrake–Rupley SASA）：

| 指标 | 9BFX 实验 | T1 非共价 | **T1cov 显式共价** | T2 无胶对照 |
|---|---|---|---|---|
| binding_confidence | — | 0.903 | **0.931** | 0.317 |
| CYPA–KRAS 界面置信 protein_iptm | — | 0.819 | 0.754 | 0.21–0.68(波动) |
| 胶↔Cys12 (Å) | 1.68(共价) | 4.66 | **1.96** | — |
| 蛋白–蛋白 BSA (Å²,近似) | ~425 | ~337 | **~417** | ~580\* |
| **KRAS 界面残基 vs 9BFX** | — | 100%/70% | 89%/80% | 28%/50% |
| **CYPA 界面残基 vs 9BFX** | — | 67%/55% | **90%/82%** | 64%/64% |

**✅ 结论**：含胶预测**复现 9BFX 天然 CYPA–KRAS 界面**（KRAS 侧界面残基 100% precision）；无胶对照界面置信度崩塌、变成非天然硬对接——印证"分子胶诱导新生界面"。进一步**显式声明 Cys12 共价键**后，胶↔Cys12 收紧到 1.96 Å、BSA≈实验、CYPA 回收率升到 82%——**约束越贴合真实化学，预测越贴近晶体**。
**⚠️ 诚实点**：这是协同性的**定性代理**，非 FEP 的 ΔΔG_coop 金标准；BSA 为近似 SASA。

---

## 第三部分 · 设计 Playbook：从案例提炼的可迁移准则

> 全文见 [`design-notes/molecular-glue-design-playbook.md`](design-notes/molecular-glue-design-playbook.md)。一句话：**靶点没口袋时，别找口袋——招募高丰度伴侣蛋白(presenter)，设计小分子让两者共建一个原本不存在的识别界面(neomorphic interface)。**

- **准则 A｜范式**：平坦 PPI/持续活化开关 → 造界面而非找口袋。
- **准则 B｜presenter 三选一**：①胞内高丰度广泛表达；②有化学起点（已知高亲和天然配体，如 sanglifehrin A↔CypA）；③与靶面有潜在**表面互补**。
- **准则 C｜工具分子先行**：双功能最小体（presenter 端 + 靶点栓系端），常用**共价栓系**先锁住弱三元、拿可解析晶体（8G9Q）作 SAR 起点。
- **准则 D｜四个优化杠杆**：①大环化/预组织（降熵罚）②弹头工程（共价/非共价、反应性调谐）③连接段构象（调弹头姿势，kobs 可数量级提升）④界面优化（诱导 presenter 表面重排）。
- **准则 E｜弹头‑残基匹配**：选择性来自**界面几何**不只弹头反应性（弱亲核 Asp12 用温和弹头 aziridine + 界面把弹头摆到位）。
- **准则 F｜平台化**：同一 presenter + 同一原理 → 多靶点/多突变分支。
- **反直觉推论**：presenter **不是靠"柔"通吃**——它是**刚性锚 + 本就带的潜在互补面 + 少量 rim 可塑**；真正"编程新界面"的是小分子；且真实campaign是把胶**刚化**来预付熵（"越柔越好"是错的）。

---

## 第四部分 · 设计 Pilots：把设计决策变成可计算判据

### Pilot 1 · presenter 选择（CypA vs FKBP12）
> 详见 [`design-notes/pilot1/RESULTS.md`](design-notes/pilot1/RESULTS.md)。把"为什么选 CypA"量化成**界面静电互补**（[`electro_complementarity.py`](design-notes/pilot1/electro_complementarity.py)）。

- **✅ 真实 9BFX**：CypA 面偏碱(+5) ↔ KRAS 效应面偏酸(−5)，**8 盐桥、0 冲突**，EC=+1.15——复现"局部电荷互补"。
- **⚠️ FKBP12 对照**（apo 对接）：找不到干净互补界面，EC+0.11、**首现 2 个同电荷冲突**、binding_confidence 最低(0.171)——方向上支持"选 CypA 不选 FKBP12"。
- **意义**：presenter 选择可被量化成**筛选器**。

### Pilot-flex · presenter 是刚是柔？
> 详见 [`design-notes/pilot_flex/README.md`](design-notes/pilot_flex/README.md)。apo CypA(2CPL) vs 结合态 CypA(9BFX) Kabsch 叠合。

- **全局 Cα RMSD = 0.41 Å**（晶体噪声级）→ CypA **整体刚性**，不是柔性适配器。
- 界面 rim 仅个别残基动：**Asn71 2.15 Å、Arg69 1.39 Å** = 局部诱导契合。
- **意义**：修正"presenter 靠柔性适配"的直觉——适配新界面靠**重设计胶**，不靠 presenter 变柔。

### Pilot 2 · 平台外推 STAT3：判定走哪条分支
> 详见 [`design-notes/pilot2/RESULTS.md`](design-notes/pilot2/RESULTS.md)。STAT3 = 转录因子、直接抑制失败的典型。

- **✅ 功能面刻画**（真实 1BG1 活化态二聚体）：STAT3 二聚界面近中性(+1)，核心是 **pTyr705(−2)↔SH2 Arg609(+) 的磷酸"锁孔"**——与 KRAS 的宽阔酸性面(−5)**性质完全不同**。
- **⚠️ 遮盖型 presenter 预筛**：CypA/FKBP12 对 STAT3‑SH2 都**抓不住**——binding_confidence **0.0002 / 0.019**，比 CypA–KRAS apo(0.317) 低约 1000×/17×。
- **方法学教训**：CypA–STAT3 的 **EC(+0.64) 反而高于** CypA–KRAS apo(+0.38)——任意对接能凑电荷，**EC 单看会误导**；判别要靠 **binding_confidence + 靶面生物学**（再证 §1.6）。
- **判定**：STAT3‑SH2 **不适合 CypA 式遮盖**，应走**降解分支(E3)或口袋直给**。

### Pilot 3 · STAT3 降解分支：CRBN 诱导邻近
> 详见 [`design-notes/pilot3/RESULTS.md`](design-notes/pilot3/RESULTS.md)。换分支的判据（[`proximity_analyze.py`](design-notes/pilot3/proximity_analyze.py)）：不看静电遮盖，看**三元 + 赖氨酸递呈**。

- 真实 STAT3 PROTAC **SD‑36**（已验证强效降解剂）招募 CRBN：三元 **binding_confidence 0.871**、protein_iptm 0.856、ligand_iptm 0.97——与遮盖分支的 ~0 天壤之别。
- **PROTAC 双端锚对**：谷氨酰亚胺端进 **CRBN tri‑Trp 口袋 2.84 Å**，pTyr 模拟膦酸端锚 STAT3‑SH2（2.47 Å）。
- **赖氨酸递呈**：STAT3 **K707/K709 被递到 CRBN 面前 6.1 Å**，25 Å 内还有 K658/K685 等——充足泛素化把手。
- **结论**：STAT3 能被这套范式攻克，但走**降解分支(PROTAC+CRBN)**——与 SD‑36 真实有效完全一致，闭环了"先判分支、再按分支判据建模"。
- **⚠️ 诚实点**：只建 CRBN‑TBD（无 DDB1/CUL4/E2~Ub 全酶），可及区是几何代理；PROTAC 三元预测是公认难题。

---

## 第五部分 · 总结果表 + 设计决策流程

### 全部 Boltz 体系一览

| 类别 | 体系 | binding_confidence | 关键读出 |
|---|---|---|---|
| 二元·共价 | KRAS G12C + Sotorasib | 0.925 | 弹头@Cys12 2.75Å, 占 S‑IIP |
| 二元·共价 | KRAS G12C + Adagrasib | 0.971 | 交叉印证 |
| 二元·非共价 | KRAS G12D + MRTX1133 | 0.995 | Asp12 盐桥 3.0Å |
| 二元·对照 | KRAS WT + MRTX1133 | 0.988 | Gly12 4.14Å 无盐桥 |
| 三元·胶(非共价) | CYPA:RMC‑6291:KRAS | 0.903 | 复现 9BFX 界面 |
| 三元·胶(共价) | 同上 + 显式 Cys12 键 | **0.931** | Cys12 1.96Å, CYPA回收82% |
| 三元·无胶对照 | CYPA + KRAS | 0.317 | 界面崩塌 |
| 预筛·遮盖 | CypA + KRAS(apo) | 0.317 | 基线 |
| 预筛·遮盖 | FKBP12 + KRAS(apo) | 0.171 | 互补差 |
| 预筛·遮盖 | CypA + STAT3‑SH2 | **0.0002** | 遮盖分支否定 |
| 预筛·遮盖 | FKBP12 + STAT3‑SH2 | **0.019** | 遮盖分支否定 |
| 诱导邻近·降解 | CRBN:SD‑36:STAT3‑SH2 | **0.871** | 降解分支证实, Lys递呈 |

### 可计算的设计决策流程（本工作的方法论产出）

```
新的不可成药靶点
    │
    ├─[1] 刻画靶"功能面"性质（真实结构/预测）
    │      宽阔带电平面? → 偏遮盖分支    磷酸/口袋热点? → 偏口袋/降解分支
    │
    ├─[2] 分支预筛（按分支选判据）
    │      遮盖分支: presenter×靶 的 EC+盐桥/冲突+binding_confidence   (Pilot 1/2)
    │      降解分支: E3+PROTAC 三元 + 赖氨酸递呈(可及区)                (Pilot 3)
    │      ※ 教训: 别只看单一分数(EC)，联合 binding_confidence 与生物学
    │
    ├─[3] 命中分支内优化（playbook 四杠杆）
    │      大环化/预组织 · 弹头工程 · 连接段姿势 · 界面/共价锚定
    │      计算读出: protein_iptm/BSA/界面回收率/共价距离(已在 demo 验证)
    │
    └─[4] 金标准闭环(未做): MD/FEP 协同自由能 + 实验 k_off/HDX/结构
```

---

## 第六部分 · 方法学诚实边界总览

- **Boltz‑2 ≠ 金标准**：单结构 + 标量/校准置信度的深度学习预测，适合**快速假设生成与位姿/几何合理性检查**，不替代 MD/FEP 对"稳定结合、协同自由能"的定量。
- **apo 蛋白‑蛋白预测本就低置信**：有意义的是**同口径的数量级对比**，不是绝对值。
- **静电互补 EC 是概念性相对度量**（非泊松–玻尔兹曼），且**单看会误导**——须与 binding_confidence + 生物学联合判读。
- **BSA 为近似 Shrake–Rupley SASA**，只做相对比较。
- **PROTAC/三元预测是公认难题**；共价与降解的"可及区"是几何代理，非催化胜任度证明。
- **构件为结构域片段**（KRAS 1–169 / STAT3‑SH2 / CRBN‑TBD），非全长语境。
- 所有机制性数字以**原始论文/实验**为准；时效性信息（如某药 FDA 批准日期）转述自来源、未独立核实。
- **专利/FTO**：Warp Drive/RevMed 的 US20160199506A1 等覆盖"大环 + presenter + 靶"三元范式，新设计需做规避分析。

---

## 第七部分 · 文件索引与复现

```
kras-druggability/
├── REPORT.md                     ← 本总报告
├── README.md                     ← 科学背景全文（第一部分详版）
├── demo/
│   ├── RESULTS.md                ← 二元复合物 demo 结果
│   ├── analyze_pocket.py         ← 口袋接触分析（纯 numpy）
│   ├── inputs.md                 ← 可复现输入（序列/SMILES/ID）
│   ├── structures/*.cif          ← 4 个二元预测结构
│   └── tricomplex/
│       ├── RESULTS.md            ← 分子胶三元 demo（对照 9BFX + 共价）
│       ├── tri_analyze.py        ← 界面/BSA/桥接分析
│       └── structures/*.cif      ← 9BFX + T1/T1cov/T2 结构
└── design-notes/
    ├── molecular-glue-design-playbook.md   ← 设计 playbook（第三部分详版）
    ├── pilot1/   (presenter 选择, electro_complementarity.py)
    ├── pilot_flex/ (CypA 刚性, rmsd.py)
    ├── pilot2/   (STAT3 分支判定)
    └── pilot3/   (STAT3 降解三元, proximity_analyze.py)
```

**复现要点**：所有分析脚本纯 numpy、无第三方依赖，`python3 <script>.py <cif> --args` 即可重跑；Boltz 预测的 prediction_id 与输入规格记录在各 pilot 的 RESULTS.md 与 `demo/inputs.md`。

---

## 参考文献

关键原始来源（完整列表见 [`README.md`](README.md) §11 与各 RESULTS）：
- Ostrem/Shokat et al., *Nature* 2013（G12C 二硫键栓系 S‑IIP）
- Schulze et al., *Science* 2023（CYPA 分子胶三元, 初始三元 PDB 8G9Q）；Weller et al., *Science* 2025（G12D 共价 aziridine）
- Revolution Medicines RMC‑6291/6236/9805；专利 US20160199506A1（Warp Drive Bio / Verdine）
- Becker et al., *Nature* 1998（STAT3 PDB 1BG1）；Zhou et al. 2019（STAT3 降解剂 SD‑36）
- 结构/数据库：RCSB PDB（9BFX / 8G9Q / 1BG1 / 2CPL / 5FQD）、ChEMBL v34、UniProt
- 计算：Boltz‑2.1（本工作全部预测）

---

*本报告为技术学习与方法探索整理；所有 ⚠️ 标注内容为假设生成层面，非监管级结论。定量结论请以原始论文与金标准实验/计算为准。*
