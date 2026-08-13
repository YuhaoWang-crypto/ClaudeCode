# ssDNA aptamer library candidates for EAB sensors against NPY and PP (PPY)

**In-silico prioritisation only. Nothing here has been experimentally tested.**
Every number is labelled `[MEASURED]` (from a cited publication),
`[DERIVED]` (computed deterministically from measured data),
`[PREDICTED]` (model output) or `[DESIGN]` (a proposal with no support yet).

Reproduce with `python3 -m aptamer_eab.run_all`.

---

## 摘要 (TL;DR)

1. **NPY 已经有一个能用的 ssDNA aptamer，而且已经做成了 EAB 传感器。** ssDNA
   aptamer **4.31**（Mendonsa & Bowser, CE-SELEX, Kd ≈ 0.3 µM）被 Seibold 等人
   直接用作 E-AB 识别元件：signal-ON，100 Hz SWV，信号变化 230 %，缓冲液中
   Kd = 385 nM、血清中 56 nM，LOD 162 pM。**先做这条，不要从零开始。**
2. **关键实验事实：把 4.31 截短成 40 nt 的核心区（去掉引物臂）后，E-AB 完全没有
   信号。** 引物臂不是"多余的扩增把手"：折叠分解显示它同时承担 **5′臂↔3′臂 9 bp
   的长程闭合茎**（开关的力学耦合，占全部配对 37.5 %）和 **5′臂↔核心 5 bp 的双链**
   （正交接触概率指认的界面，占 20.8 %）。所以本 library 默认保留 80 nt 架构。
3. **NPY / PP / PYY 的 C 端几乎完全一样**（36 位中 14 位三者全同，其中 8 位集中在
   24–36 段）。**针对 C 端片段选出来的 aptamer 必然全家族交叉反应。**
   要做"短片段"，应选 **N 端**：NPY 用 **5–16 (PDNPGEDAPAED)**，PP 用
   **1–12 (APLEPVYPGDNA)**。
4. **PP 最经济的起点是 4.31 本身。** 4.31 对 NPY 相对 hPP 只有 42 倍选择性 —— 也就是
   说它**确实结合 PP**（约 12 µM）。这是一个真实的结合先验；用高 doping (30 %) 的
   doped library + 对 NPY 反向筛选，就有机会把选择性反转过来。
5. 交付物：**58 条候选序列**（`aptamer_eab/output/ORDER_PANEL.csv`）+ **4 条可直接
   下单的简并寡核苷酸**（naive N40 库 / NPY 成熟库 / PP 重定向库 / 引物）。

---

## 1. 先看已有的东西（这一步省掉最多钱）

| 试剂 | 化学 | 长度 | Kd | 来源 | 标签 |
|---|---|---|---|---|---|
| **aptamer 4.31** | ssDNA | 80 nt | 0.3 ± 0.2 µM (NPY) | Mendonsa & Bowser, *JACS* 2005 | `[MEASURED]` |
| 4.31 做成 E-AB | ssDNA, 5′-SH / 3′-MB | 80 nt | 385 nM (buffer), 56 nM (serum) | Seibold et al., *ACS Sens* 2023 | `[MEASURED]` |
| **DP3** | 2′-NH₂-嘧啶 RNA | 95 nt | 370 nM | Proske et al., *JBC* 2002 | `[MEASURED]` |
| PP aptamers | ssDNA | — | — | Ali et al., *ACS Omega* 2019 | `[MEASURED]`，**序列未取到** |

**4.31 序列** `[MEASURED]`（5′→3′，下划线为固定引物臂）:

```
AGCAGCACAGAGGTCAGATG CAAACCACAGCCTGAGTGGTTAGCGTATGTCATTTACGGA CCTATGCGTGCTACCGTGAA
└──── 5′ arm 20 nt ──┘└──────── selected core 40 nt ─────────────┘└──── 3′ arm 20 nt ──┘
```

E-AB 构型：`5'-HS-(CH2)6-[80 nt]-MB-3'`，signal-ON，SWV @ 100 Hz，
动态范围 20–600 nM，最大信号变化 230 ± 15 %，PYY 干扰仅 4 ± 1 % SC。

### 三件必须先知道的事

**(a) 截短会杀死信号。** 同一篇论文报告 40 nt 的裸核心（STLA）"no appreciable
signal"。这直接决定了 library 的架构：核心区两侧保留引物臂。本流程算出的二级结构
给出了机制解释 —— 80 nt 构型的配对分解是 core↔core 10 bp / **5′臂↔3′臂 9 bp** /
**5′臂↔核心 5 bp**；裸核心只有 10 bp 且完全没有臂，是两个互不相干的局部发夹。
删掉臂等于同时失去长程闭合茎和臂–核心双链。详见 §4b。`[PREDICTED]`

**(b) 已发表传感器的动态范围可能不匹配你的样本。** 20–600 nM 的窗口远高于文献报道
的血浆 NPY 水平（通常为低 pM，即数十 pg/mL 量级）。若目标是**血浆**定量，亲和力
成熟不是可选项而是必需项；若目标是**局部释放/微透析/组织间液**，nM 窗口通常够用。
请先用你们自己的参考区间核对。`[MEASURED 文献值 / 推断]`

**(c) PP 的 aptamer 序列没拿到。** Ali et al. *ACS Omega* 2019 是开放获取的，但
`pubs.acs.org` 在本沙箱返回 403，Europe PMC 也没有镜像。**这是目前 PP 方向最高价值
的缺口** —— 请从 PDF/SI 里把序列取出来填进 `aptamer_eab/known.py` 的 `PP_KNOWN`，
之后 C 库可以直接以它们为种子重跑。

---

## 2. 靶点分析：为什么"短片段"必须选 N 端

`[DERIVED from UniProt P01303 / P01298 / P10082]`

```
pos      1234567890123456789012345678901234567890
NPY      YPSKPDNPGEDAPAEDMARYYSALRHYINLITRQRY
PP       APLEPVYPGDNATPEQMAQYAADLRRYINMLTRPRY
PYY      YPIKPEAPREDASPEELNRYYASLRHYLNLVTRQRY
all-same  *  *  *   *  *    *   ** * *  ** **
```

| 比较 | 一致性 |
|---|---|
| NPY vs PP | **50.0 %** (18/36) |
| NPY vs PYY | 63.9 % (23/36) |
| PP vs PYY | 44.4 % (16/36) |

三者全同的 14 个位置：2, 5, 8, 12, 15, 20, 24, 25, 27, 29, 32, 33, 35, 36 ——
**其中 8 个挤在 24–36 段**。这是 PP-fold 家族共同的 C 端 `…RY-NH₂` 受体结合面。

### 12-mer 滑窗判别力扫描 `[DERIVED]`

判别分 = (对**两个**旁系同源都不同的位点数 + 0.5 × 其中氨基酸类别也改变的位点数) / 窗口长度。
只有对 **PP 和 PYY 都不同**的位点才计分 —— 只能区分 NPY/PP 却区分不了 PYY 的表位，
对家族检测没有意义。

| 靶点 | 最佳窗口 | 序列 | 判别分 | 最差窗口（不要选） |
|---|---|---|---|---|
| **NPY** | 3–14 / **5–16** / 6–17 | SKPDNPGEDAPA / **PDNPGEDAPAED** | 0.625 | 24–35 `LRHYINLITRQR` → **0.08** |
| **PP** | **1–12** / 2–13 / 3–14 | **APLEPVYPGDNA** | 0.833 | 22–33 `ADLRRYINMLTR` → 0.38 |

**推荐**

- **NPY → 5–16 `PDNPGEDAPAED`**。判别分与 3–14 并列最高，但完全避开 Tyr1-Pro2。
  DPP-4 在体内把 NPY(1-36) 切成 NPY(3-36)，所以避开 1–2 位的表位对两种形式都响应
  ——若你要测 **total NPY**，这是正确选择。
  反过来，如果你**想区分** NPY(1-36) 与 NPY(3-36)（两者受体选择性不同：Y1 vs Y2），
  就应该刻意把 Tyr1-Pro2 放进表位，用 1–12 窗口。这是一个产品定义问题，需要你们定。
- **PP → 1–12 `APLEPVYPGDNA`**，判别分 0.833，是全家族最容易区分的一段。
  注意 PP 的 N 端也是 Xaa-Pro（Ala1-Pro2），即 DPP-4 的识别基序，同样的完整性问题
  可能存在，建议在你们的基质里先验证。`[推断，需实验确认]`

### 一个对 PP 不利的结构事实 `[MEASURED，文献]`

- **NPY** 在水相缓冲液中大体无序/部分螺旋，只有结合膜/胶束时才折叠 ——
  **线性短肽是自由 NPY 的合理替身**，用片段做 SELEX 风险较低。
- **PP** 折叠成紧密的 PP-fold（聚脯氨酸 II 螺旋贴在 α-螺旋上，疏水核心；
  典型晶体结构 PDB 1PPT）。**线性 N 端片段并不还原完整 PP 呈现的表面。**

因此本设计的立场是：**NPY 可以用片段选；PP 应当用完整肽选，片段只用于反向筛选的
概念界定。** 如果 PP 也必须用片段，就必须在结合验证阶段用完整 PP 回验，否则很可能
选出只认线性肽、不认天然激素的 aptamer。

---

## 3. 打分函数：先说它错在哪

E-AB 的信号完全来自结合诱导的构象变化 —— 所以要选的不是"折叠得最稳"的序列，而是
**可被靶点重塑**的序列。

**第一版打分函数是错的，必须说明。** 它用了原始 `P(MFE)` 和原始 MFE，两者都强烈依赖
长度；在 80 nt 上塌缩，结果**把已发表、能工作的 4.31（230 % 信号变化）排在它自己的
scramble 对照之下**。一个与唯一已知真值反相关的打分函数比没有打分函数更糟。

现版本只用**长度归一化**的描述量，并对本靶点仅有的两个已发表 E-AB 结果做锚定检查：

| 锚点 | 实验事实 | `eab_score` | terminal_bp | MFE/nt | ens_div/nt |
|---|---|---|---|---|---|
| 4.31 完整 80 nt | `[MEASURED]` 有效，230 % SC | **80.9** | 2 | −0.150 | 0.169 |
| 4.31 裸核心 40 nt | `[MEASURED]` 无信号 | 48.7 | **0** | −0.147 | 0.072 |
| scramble 对照 80 nt | 预期无效 | 52.8 | 2 | −0.076 | **0.344** |

**ANCHOR CHECK: PASS**（有效构型需高出两个阴性对照 >10 分）。`run_all` 在这一步
失败时会直接中止，不会拿一个坏打分函数去出订单。

> **必须直说的话：只有 2 个已知结果，而且窗口是在看过这两条序列之后设的。这是
> "锚定"，不是"验证"。两个点无法验证一个五项的启发式。** 它只用来给合成清单排序，
> 绝不是亲和力预测，也不是信号变化预测。
>
> 由此还有一个推论：panel 里 50 条"分数高于 4.31"**不构成任何结合更强的证据** ——
> 打分函数本来就是照着 4.31 的结构锚定的，高分只意味着"更像 4.31 的结构类别"。

打分项（权重）：终端茎 `terminal_bp` (0.30) · MFE/nt (0.25) · 归一化系综多样性 (0.20)
· 配对比例 (0.15) · 3′ 端自由碱基数 (0.10)；GC 与均聚物越界扣分。
**`terminal_bp` 的作用范围已按 §4b 更正**：它区分的是**架构**（framed vs 无臂），
在 A–D 库内部近乎常数（固定臂人人相同），库内区分力实际来自 `mfe_per_nt` 与
`norm_ens_div`。把这 0.30 的权重读作"强制架构"，不要读作"给核心排序"。
`self_dimer_dg` **只报告不计分** —— 有效的 4.31 自二聚（−25.8 kcal/mol）比它的
scramble（−12.2）更强，把它计入惩罚会直接把唯一已知有效的序列打下去。

---

## 4. 六个子库

| 库 | 靶点 | 策略 | 数量 | 备注 |
|---|---|---|---|---|
| **A** | NPY | 4.31 核心 doped 15 % + 保留引物臂 | 401 | 先验最强、风险最低 |
| **B** | NPY | de novo 结构池（茎环/双发夹/三向连接/随机） | 631 | 逃离 4.31 表位用 |
| **C** | PP | 4.31 核心 doped **30 %**，配对进 5′ 臂的 5 个位点只 doped 5 % | 400 | 利用 4.31 对 hPP 的残余亲和力；两级 doping 的理由见 §4b |
| **D** | PP | de novo 结构池（独立随机流） | 642 | |
| **S** | 通用 | 28–40 nt 无引物臂短开关 | 891 | 风险较高（见截短实验），但电子转移快 |
| **G** | 通用 | G-四链体，15–32 nt | 701 | **单独排序**，见下 |

### 关于 G4 轨道：两次踩坑，都记录在代码里

G4 值得单列，因为 E-AB 方法的原型传感器（凝血酶/TBA）就是一个 15 nt 的 G4，而且
G4 DNA 紧凑（电子转移快、相对信号变化大）、天然抗核酸酶 —— 对血清/体内运行很重要。

- **坑 1：ViennaRNA 的 DNA 参数里根本没有 G4 项。** 真实 G4 会被当成无规卷曲，
  所以任何 G4 候选在 `eab_score` 上必然垫底。这是建模盲区，不是生物学结论。
- **坑 2：直接最大化 G4Hunter 分数是错的目标。** 第一版这么做，池子塌缩成近乎
  纯 G 的均聚物（GC 0.86–0.96，G 连续段 10 个以上）—— 这类序列形成分子间 G-wire
  和聚集体，合成困难，成不了单层膜。**用 TBA 校准就能看出来：TBA 的 G4Hunter 只有
  1.13，低于通常 1.2 的"可能成 G4"阈值。** 现在改为按 aptamer 导向的
  `g4_quality` 排序（TBA 式窗口 + 恰好四段 G-tract + GC/均聚物硬约束），TBA 本身作为
  参照行留在库里（`g4_quality = 70.6`）。

---

## 4b. 正交验证的结果（HDOCK + 接触概率）与它改变了什么

一轮独立的正交验证（HDOCK 刚体对接 9 组交叉 + Boltz-2 distogram 接触概率 5 个复合物）
的报告并入本项目。结论要分成两半读。

### 不能用的那一半：三个打分判据全部没通过 decoy 对照

| 方法 | NPY−PP 分离度（信号） | 4.31−打乱序列（噪声） | 信号/噪声 |
|---|---|---|---|
| Boltz-2 ipTM | 0.054（符号还反了） | 0.185 | 0.29 |
| 接触概率 | 1.10 | 1.26 | 0.87 |
| HDOCK | 3.92 | 18.03 | 0.22 |

**三者全部 < 1**：把靶标从 NPY 换成 PP（真实可测的 42 倍差异），对分数的影响比把
aptamer 换成随机序列还小。更糟的是，HDOCK 和接触概率**都把实验已确认无 EAB 信号的
40-nt 截短版排在全长 4.31 前面**（HDOCK −278.3 vs −235.6；接触概率 13.78 vs 13.23）。
刚体对接在跟踪形状互补与静电，接触概率的判别项贡献近乎为零
（`Spearman(cp_on, 特异性得分) = 0.9933`）。

**这与本项目的立场一致，且明确排除了一条路：不要用共折叠分数、接触概率或对接分数给
SELEX 候选排序。** 本报告第 7 节最后一条已经说明为什么没有做共折叠打分，这轮验证把它
从"预判"变成了"实测结论"。

**关于本项目 `eab_score` 的诚实说明：它同样没有通过独立验证。** 它确实把 40-nt 截短版
排在全长之下（48.7 vs 80.9），而 HDOCK 和接触概率都排反了 —— 但这**不构成本方法更好的
证据**：`eab_score` 的窗口正是在看过这一对序列之后设的，排对是构造使然，不是发现。
三种方法现在都属于"未验证"，只是失效模式不同。

### 能用的那一半：固定臂在界面上 —— 两条独立证据吻合

接触概率给出的**区段分解**（比排序弱得多的一个主张）与本项目从二级结构独立算出的结果
高度一致：

| 证据来源 | 结论 |
|---|---|
| 接触概率（Boltz-2 distogram） | 5′ 固定臂承担 **32.1 %**（NPY×4.31）的接触；最强接触核苷酸 **T14、C15、A16** |
| 本项目二级结构（ViennaRNA，**完全没有蛋白**） | 5′ 臂位点 **14、15、18、19、20** 配对进入核心区（5 bp），核心区对应位点 44、45、46、59、60 |

**T14 和 C15 精确落在两法共同指认的位点上**（A16 位于 14-15 与 18-20 两段配对之间的
环内）。一个是蛋白–DNA 距离分布图，一个是无蛋白的核酸二级结构，两者独立指向同一段
5′ 臂 —— 这比任何一个排序分数都更值得相信，因为它是**位置**主张而不是**优劣**主张。

同时，本项目的 4.31 折叠给出完整的配对分解：

| 配对类型 | bp | 占比 |
|---|---|---|
| core ↔ core | 10 | 41.7 % |
| **5′arm ↔ 3′arm** | **9** | **37.5 %** |
| **5′arm ↔ core** | **5** | **20.8 %** |

而 40-nt 裸核心只有 5+4+1 = 10 bp，且完全没有臂。**这就是截短失效的结构解释**：
删掉固定臂，同时失去了 37.5 % 的长程闭合茎（开关的力学耦合）和 20.8 % 的臂–核心双链
（接触概率指认的界面）。

### 由此做出的两处实际修改

**1. 修正一个我说过头的说法。** 原先把 `terminal_bp` 描述成"区分能工作与不能工作的
单一特征"。这是错的：**固定臂在每一条 framed 候选里都完全相同**，所以
`terminal_bp`、`arm_arm_bp`、`arm_core_bp` 在 A–D 库内部近乎常数 ——
组成匹配的 scramble 对照同样有 `arm_arm_bp = 9`、`arm_core_bp = 4`。
`terminal_bp` 实际区分的是**架构**（framed vs 无臂），这正是 anchor check 需要的；
库内部的区分力来自 `mfe_per_nt` 和 `norm_ens_div`（4.31 = −0.150 / 0.169，
scramble = −0.076 / 0.344）。代码注释与权重解读已按此更正。

**2. C 库（PP 重定向）改用两级 doping。** 检查发现平铺 30 % doping 会把臂–核心双链
打掉：面板中 12 条只有 9 条保住（`arm_core_bp ≥ 4`）。把配对进入 5′ 臂的 5 个核心位点
（核心内 1-based **24、25、26、39、40**，即构建体 44、45、46、59、60）降到 5 % doping 后，
**保住比例升到 11/12，中位数从 4 bp 回到 6 bp**。A 库（15 % doping）本来就是 12/12，
不需要改。合成上这只是同一条简并寡核苷酸用两个 doping 水平，任何合成仪都支持。

`ORDER_PANEL.csv` 现在带 `arm_arm_bp` / `arm_core_bp` / `keeps_arm_core_duplex` 三列，
可据此对 framed 候选再筛一道。

### 关于 3′ 端（MB 位）

接触概率给出 3′ 臂只占 **7.2 %** 的接触。这与已发表的 3′-MB 端标记产生 230 % 信号变化
自洽：3′ 端**不直接参与结合**，所以它的位移是构象重排的**读出**，而不是结合本身的一部分
—— 这恰是一个好的 E-AB 报告基团位置应有的性质。本项目 `tail3_free` 项（奖励 3′ 端游离）
与之方向一致，4.31 的最后 5 nt 全部未配对。这是弱的支持性证据，不改变任何设计决定。

---

## 5. 交付物

### 5.1 定义好的候选面板 — `aptamer_eab/output/ORDER_PANEL.csv`（58 条）

22 条 NPY · 22 条 PP · 14 条与靶点无关的骨架（6 条短开关 + 8 条 G4）。
面板内任意两条核心区的编辑距离 ≥ 6，避免把同一条序列订 16 遍。

**NPY 前 5** `[PREDICTED / DESIGN]`

| id | 来源 | 架构 | nt | MFE | score | 5′臂↔核心 | 序列 (5′→3′) |
|---|---|---|---|---|---|---|---|
| A01 | A 种子 | doped 4.31 (6 mut) | 80 | -12.8 | 88.4 | **6 bp 保** | `AGCAGCACAGAGGTCAGATGCCAACCTCCGGCTGAGTGGTTGGCGTATGTCATTTCCGGACCTATGCGTGCTACCGTGAA` |
| B13 | B de novo | dual hairpin | 77 | -12.6 | 86.9 | 3 bp 破 | `AGCAGCACAGAGGTCAGATGTGTCCCGACATCGAAGTCGGGTTCTGTATCAACACAGCCTATGCGTGCTACCGTGAA` |
| B14 | B de novo | dual hairpin | 80 | -13.3 | 86.7 | 0 bp 破 | `AGCAGCACAGAGGTCAGATGTCACGATGGATTTTGCTCCATCACGACGGCGTAACCGTCGCCTATGCGTGCTACCGTGAA` |
| A02 | A 种子 | doped 4.31 (2 mut) | 80 | -12.7 | 86.1 | **7 bp 保** | `AGCAGCACAGAGGTCAGATGCAAACCACAGCCTAAGTGGTTAGCGTATCTCATTTACGGACCTATGCGTGCTACCGTGAA` |
| A03 | A 种子 | doped 4.31 (5 mut) | 80 | -13.6 | 85.8 | **6 bp 保** | `AGCAGCACAGAGGTCAGATGCAAACCACATCAGGAGTGGTTAGCATATGTCATTTGCGGACCTATGCGTGCTACCGTGAA` |

> A02 只有 2 个突变，是最保守的一步试探 —— 如果它比野生型 4.31 好，说明这条路走得通。

**PP 前 5** `[PREDICTED / DESIGN]`

| id | 来源 | 架构 | nt | MFE | score | 5′臂↔核心 | 序列 (5′→3′) |
|---|---|---|---|---|---|---|---|
| C23 | C 交叉种子 | doped 4.31 (13 mut) | 80 | -12.9 | 89.5 | **6 bp 保** | `AGCAGCACAGAGGTCAGATGCGAACGAAACCTTGAGGCGTTCGCGTATAACATTAGGGGACCTATGCGTGCTACCGTGAA` |
| D35 | D de novo | stem-loop | 76 | -12.2 | 88.7 | 3 bp 破 | `AGCAGCACAGAGGTCAGATGTTAGAACCTTCGGTACACCGGCGACCGAAGGCTGTTCCTATGCGTGCTACCGTGAA` |
| D36 | D de novo | dual hairpin | 77 | -12.2 | 87.4 | 0 bp 破 | `AGCAGCACAGAGGTCAGATGTTAGTGGTATATACCACCGTCCAGCCATCAGGGCTGGCCTATGCGTGCTACCGTGAA` |
| D37 | D de novo | stem-loop | 76 | -12.2 | 87.2 | **5 bp 保** | `AGCAGCACAGAGGTCAGATGGACCGCAACCACCTTAGCGGTGAGGTGGTTACAATACCTATGCGTGCTACCGTGAA` |
| C24 | C 交叉种子 | doped 4.31 (9 mut) | 80 | -12.7 | 86.5 | **8 bp 保** | `AGCAGCACAGAGGTCAGATGCAAACAACCGCCTGACTAGTTTGCATGTGTCCTTTACGGCCCTATGCGTGCTACCGTGAA` |

> `5′臂↔核心` 列 = 该候选保留了几个碱基对进入 5′ 固定臂（4.31 为 5 bp）。
> D 库是 de novo 核心，本来就没有保留这段双链的义务；C 库经两级 doping 后 11/12 保住。

**短开关轨道（S，无引物臂，30–38 nt）** — 电子转移更快、相对信号更大，但对本家族
风险更高（见 4.31 截短实验）：

| id | 架构 | nt | MFE | 序列 |
|---|---|---|---|---|
| S45 | 三向连接 | 35 | −5.4 | `TGGCCTGCTACTAGCATATGGACATAGGCCATCAA` |
| S46 | 茎环 | 30 | −4.9 | `ATAAGAGGAAGTAGGGTATCCTCTTCGCCT` |
| S47 | random | 38 | −6.5 | `ACAGGATACCGCAGGCTCTTACGCGGTTCCATCCACAG` |

**G4 轨道（G，20–25 nt）** — 按 `g4_quality` 排序，TBA 参照 = 70.6 / G4Hunter 1.13：

| id | nt | GC | G4Hunter | 序列 |
|---|---|---|---|---|
| G51 | 20 | 0.70 | 1.40 | `TGGGTTTCTGGGCGGGCTGG` |
| G52 | 20 | 0.70 | 1.40 | `CGGCTGGGAGGGAACGGGAA` |
| G53 | 20 | 0.70 | 1.40 | `CACGGGTCGGGTGGGAGGTT` |

### 5.2 可直接下单的简并寡核苷酸 — `aptamer_eab/output/order_specs.csv`

doped library 在湿实验里是**一条简并寡核苷酸**（合成仪每个位点按比例混单体），不是
一份序列清单。CSV 里的序列列表是同一空间的 in-silico 采样，供偏好定义阵列的人使用。

| 名称 | 规格 | 用途 |
|---|---|---|
| **E1** naive N40 | `5'-AGCAGCACAGAGGTCAGATG-(N)40-CCTATGCGTGCTACCGTGAA-3'` | 两个靶点的 naive 第一轮；沿用已发表选择的引物臂，工作流与传感器可直接对比 |
| **E2** NPY 成熟 | 同上，核心区按 **15 % doping**（85 % 野生型 / 各 5 % 其它三种） | 4.31 亲和力+开关幅度成熟；升压力 + PP/PYY 反向筛选 |
| **E3** PP 重定向 | 同上，核心区按 **30 % doping**，但核心位点 **24,25,26,39,40**（构建体 44,45,46,59,60）只按 **5 % doping** | 把 4.31 骨架重定向到 PP；**从第 1 轮起就对 NPY 和 PYY 反向筛选**。两级 doping 的理由见 §4b |
| **E4** 引物 | FWD `AGCAGCACAGAGGTCAGATG` / REV `CCTATGCGTGCTACCGTGAA` | 反向引物建议 5′-磷酸化或加 poly-A 尾以再生 ssDNA |

### 5.3 E-AB 构型（沿用已发表的成功配置）

```
5'-HS-(CH2)6-[candidate]-MB-3'
```
金电极自组装单层 + 6-巯基己醇封闭；方波伏安法，扫描频率 60–500 Hz 内实测挑增益最优点
（已发表 NPY 传感器用 100 Hz）。若打算长期/体内运行，3′ 端加 inverted-dT 抗外切酶。

### 5.4 必须一起订的对照 `[MEASURED / DESIGN]`

| 对照 | 序列 | 作用 |
|---|---|---|
| **阳性参照** | 4.31 完整 80 nt | 唯一有真实数据的基准；每块芯片都应带 |
| **带引文的阴性** | 4.31 裸核心 40 nt | 已发表在 E-AB 里无信号 —— 比随机对照更有价值 |
| scramble | 组成匹配打乱 | 常规阴性 |

---

## 6. 建议的推进路线

**先做这一步（最高性价比）：** 直接合成 4.31 完整 80 nt 的 E-AB 构型，在你们的缓冲液
和基质里复现 230 % 信号变化。这一步验证的是**你们的电极工艺**，不是 aptamer ——
在自制传感器上重现一个已发表的结果，是后面所有比较的前提。同时用它测 PP 的响应：
文献只给了 42 倍选择性这一个数（而且是溶液相 CE 测的），**E-AB 格式下 4.31 对 PP 的
实际响应曲线是本项目最有价值的一个未知数**。如果 4.31 对 PP 的响应可测，C 库的
重定向策略就有了实验依据；如果完全无响应，就应该改走 D 库的 de novo 路线。

之后：

1. **NPY** — E2 doped 库做 3–5 轮亲和力成熟 SELEX，逐轮加压，**每一轮**都用 PP 和
   PYY 反向筛选。若目标是血浆，把终点亲和力目标定在低 nM 到 pM。
2. **PP** — 先把 Ali et al. 的序列补进 `PP_KNOWN`；有了它们就以它们为种子重跑 C 库。
   若拿不到，走 E3（4.31 重定向，对 NPY/PYY 全程反向筛选）+ D 库 de novo 双路并行。
3. **特异性必须实测，不能假设。** 交叉反应矩阵至少覆盖 NPY(1-36)、NPY(3-36)、PP、
   PYY(1-36)、PYY(3-36)，在目标基质里做。家族同源度这么高，特异性只能靠反向筛选
   造出来再测出来。
4. 结合验证用与 E-AB 正交的方法（SPR/BLI/MST），再回到电极上确认开关幅度。
5. 最终优化：截短扫描（保留终端茎！）、MB 位置（末端 vs 内部，已发表工作两种都试过）、
   扫描频率、单层密度。

---

## 7. 必须随交付一起说明的局限

- 全部为 in-silico 设计，**panel 中没有任何一条序列被证明能结合 NPY 或 PP**。
- `eab_score` 是二级结构上的工程启发式，**不是亲和力预测，也不是信号变化预测**，
  且只有 2 个锚点、窗口是看过锚点后设的。
- ViennaRNA 的 DNA 参数不正确处理 Mg²⁺/Na⁺，**完全不建模 G-四链体和任何三级结构**。
  G4 候选必须用实验定构象（CD：平行 G4 约 264 nm 正峰 / 240 nm 负峰；K⁺ vs Li⁺
  依赖的熔解曲线）。
- 只对 MFE 单一构象打分，离子依赖的构象系综未建模。
- NPY/PP/PYY 相互一致性 44–64 %，C 端近乎相同：**家族特异性不能假设**。
- Library 是 SELEX / 筛选的起始池，不是 binder。
- **三个正交打分判据（ipTM / 接触概率 / HDOCK）均未通过 decoy 对照**，其中两个把
  实验无信号的 40-nt 截短版排在全长前面。详见 §4b。本项目的 `eab_score` 同样未经
  独立验证，只是失效模式不同。
- **本项目未做共折叠（Boltz-2）打分。** 该 skill 自己的交叉验证结论是：对柔性
  ssDNA，ipTM 偏乐观且不具区分度，唯一可信的信号是结构 footprint 重叠 —— 而 NPY 是
  一条 36 aa 的柔性肽，没有可用的 aptamer 复合物结构来定义 footprint。在这种情况下
  跑共折叠只会产出一个看起来精确、实际无判别力的分数，所以没有做。真正的判别力
  必须来自湿实验的反向筛选。

---

## 8. 需要你们决定/补充的两件事

1. **NPY 的产品定义：测 total NPY，还是要区分 NPY(1-36) 与 NPY(3-36)？**
   前者 → 表位选 5–16（避开 DPP-4 切点）；后者 → 表位必须包含 Tyr1-Pro2，选 1–12。
   这直接改变 SELEX 时用哪条合成肽做靶点。
2. **Ali et al. *ACS Omega* 2019 的 PP aptamer 序列**（`pubs.acs.org` 在本环境 403）。
   拿到后填进 `aptamer_eab/known.py` 的 `PP_KNOWN`，C 库可以立刻以真实 PP binder
   为种子重跑，先验强度会有量级差别。

---

## 参考文献

1. Mendonsa SD, Bowser MT. *In vitro selection of aptamers with affinity for
   neuropeptide Y using capillary electrophoresis.* J Am Chem Soc 2005;127(26):9382-3.
   [doi:10.1021/ja052406n](https://doi.org/10.1021/ja052406n) (PMID 15984861)
2. Seibold JM, Abeykoon SW, Ross AE, White RJ. *Development of an Electrochemical,
   Aptamer-Based Sensor for Dynamic Detection of Neuropeptide Y.* ACS Sens 2023;8(12).
   [doi:10.1021/acssensors.3c00855](https://doi.org/10.1021/acssensors.3c00855)
   (PMID 38033269, [PMC11214579](https://pmc.ncbi.nlm.nih.gov/articles/PMC11214579/))
3. Proske D, Höfliger M, Söll RM, Beck-Sickinger AG, Famulok M. *A Y2 receptor mimetic
   aptamer directed against neuropeptide Y.* J Biol Chem 2002;277(13):11416-22.
4. Ali ASM, El-Halawany MS, Ibrahim SA, Plückthun O, Khalil ASG, Mayer G.
   *Aptasensor for Quantifying Pancreatic Polypeptide.* ACS Omega 2019;4(2):2948-56.
   [doi:10.1021/acsomega.8b03131](https://doi.org/10.1021/acsomega.8b03131)
5. Bedrat A, Lacroix L, Mergny J-L. *Re-evaluation of G-quadruplex propensity with
   G4Hunter.* Nucleic Acids Res 2016;44(4):1746-59.
6. UniProt: [P01303](https://www.uniprot.org/uniprotkb/P01303) (NPY),
   [P01298](https://www.uniprot.org/uniprotkb/P01298) (PPY),
   [P10082](https://www.uniprot.org/uniprotkb/P10082) (PYY).
