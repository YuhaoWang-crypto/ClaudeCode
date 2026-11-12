# Demo 结果：分子胶三元复合物 CYPA : 胶 : KRAS —— 对照 PDB 9BFX

> 配套正文见 [`../../README.md`](../../README.md) §7（分子胶三元复合物）。
> 结构/亲和力来自 **Boltz‑2.1**；参照为 **PDB 9BFX**（*Tri‑complex of Elironrasib (RMC‑6291), KRAS G12C, and CypA*）。运行日期 2026‑11‑12。

---

## 1. 目标与设计

复现正文 §7 的核心主张——**分子胶先劫持高丰度伴侣蛋白 CYPA，再诱导一个原本不存在的 CYPA–KRAS 新生蛋白界面（neo‑PPI），把 KRAS(ON) 钉住**——并与真实晶体结构 **9BFX** 对照。

| 体系 | 组成 | binding 模式 | prediction_id |
|---|---|---|---|
| **T1 三元（非共价）** | CYPA + KRAS(G12C) + **RMC‑6291** + GppNHp + Mg²⁺ | ligand_protein_binding（胶为 binder） | `sab_pred_gpde8tALI5RaqTWFI6Nq` |
| **T1cov 三元（显式共价）** | 同 T1，另声明 **Cys12 SG ↔ 胶弹头β‑碳** 共价键 | ligand_protein_binding + `bonds` | `sab_pred_28wbJ0CBYwkBZccMRzm8` |
| **T2 对照（无胶）** | CYPA + KRAS(G12C)（apo） | protein_protein_binding（KRAS 为 binder） | `sab_pred_dIK0dZok68JJwAlw1jwU` |
| **参照** | PDB **9BFX**（实验，ON 态 GppNHp） | — | 链 A=KRAS, D=CypA, 胶=A1AOD |

- 胶分子 **RMC‑6291（Elironrasib）= PDB 配体 `A1AOD`**，C₅₅H₈₀FN₉O₈（73 重原子），由 9BFX 直接取得。
- KRAS 用规范 1–169 G12C、CypA 用规范 1–165（与前序 demo 编号一致）；ON 态用 GppNHp(`GNP`)+Mg²⁺ 还原。
- 成本：T1≈$0.20，T2≈$0.13。
- 界面分析脚本 [`tri_analyze.py`](tri_analyze.py)（纯 numpy，自实现 Shrake–Rupley SASA 算 BSA）；完整输出见 [`tri_analysis.txt`](tri_analysis.txt)。

---

## 2. 结果一：Boltz 置信度——胶在，界面才"可信且可复现" ✅

| 体系 | protein_iptm（PPI 界面置信） | ligand_iptm | **binding_confidence** | complex_ipde（越低越准） | 3 样本 protein_iptm |
|---|---|---|---|---|---|
| **T1 非共价** | 0.819 | 0.935 | 0.903 | 0.97 | 0.819 / 0.824 / 0.789 |
| **T1cov 共价** | 0.754 | 0.916 | **0.931** | **0.99** | 0.754 / 0.769 / 0.618 |
| **T2 无胶** | 0.644 | — | **0.317** | 5.69 | 0.644 / 0.677 / **0.207** |

- 含胶时 CYPA–KRAS 界面**置信度高（protein_iptm 0.82）、三样本稳定、界面误差低（ipde≈1）**。
- 去掉胶后 **binding_confidence 从 0.90 跌到 0.32**，界面误差暴涨（ipde 5.7–11.2），且三样本剧烈波动（0.21–0.68）——Boltz 对"两蛋白到底怎么贴"**没有把握、给不出一致答案**。这正对应实验事实：CYPA 与 KRAS 本身几乎无亲和力，**必须靠胶**。

> ⚠️ 注意：protein‑protein 模式下 Boltz **总会**摆出某个相对位置（它无法输出"不结合"），所以"无胶"不等于"零界面"。真正的判别信号是**置信度 + 下面的界面残基是否天然**，而不是"有没有接触"。这本身是正文 §9 的实证。

---

## 3. 结果二：含胶预测几乎复现 9BFX 的真实界面 ✅

界面接触判据：残基最近重原子 ≤ 4.5 Å；BSA 为自实现 SASA 估算（绝对值仅供比较，约 ±10–15%）。

| 指标 | **9BFX 实验** | **T1 非共价** | **T1cov 共价** | **T2 无胶** |
|---|---|---|---|---|
| CYPA–KRAS 接触原子对 | 79 | 36 | 54 | 169 |
| KRAS 界面残基数 | 10 | 7 | 9 | 18 |
| CYPA 界面残基数 | 11 | 9 | 10 | 11 |
| 蛋白–蛋白 BSA (Å²) | ~425 | ~337 | **~417** | ~580\* |
| 胶↔Cys12 最近 (Å) | 1.68（共价） | 4.66 | **1.96** | — |
| 胶同时桥接两蛋白 | ✓ | ✓ | ✓（含 Cys12） | — |
| **KRAS 界面残基 vs 9BFX** | — | precision 100% / recall 70% | **89% / 80%** | precision 28% / recall 50% |
| **CYPA 界面残基 vs 9BFX** | — | precision 67% / recall 55% | **90% / 82%** | precision 64% / recall 64% |

**关键判读（KRAS 侧界面是"胜负手"）：**

- **T1（含胶）预测的 KRAS 界面残基 {Tyr32, Asp33, Pro34, Ile36, Asp38, Glu63, Tyr64} —— 100% 命中 9BFX**，全部落在 Switch‑I（31–38）/ Switch‑II（63–64）区，正是胶模板化出来的天然界面。
- **T2（无胶）** 虽然 BSA 看起来更大（580 Å²），但那是 Boltz 把两蛋白硬贴出的**非天然**界面：KRAS 侧动用了 N 端 Lys5–Val7、Asp54–Leu56、Arg73–Gly75 等 **9BFX 中根本不参与**的残基，precision 仅 28%。
- CypA 侧两者都与 9BFX 重叠较多——因为 CypA 的结合凹槽（Arg55/Ile57/Trp121/His126/Arg148…）是其**固有 PPIase 口袋**，胶就锚在那里；区别主要体现在 KRAS 侧。

**胶–蛋白接触残基（T1 vs 9BFX，高度一致）：**

| 胶接触 | 9BFX 实验 | T1 预测 |
|---|---|---|
| ↔ KRAS | Cys12, Tyr32, Pro34, Thr35, Ile36, Glu37, Ala59, Gly60, Gln61, Tyr64, Met67, Tyr71 | Tyr32, Pro34, Thr35, Ile36, Glu37, Ala59, Gln61, Tyr64, Met67 |
| ↔ CYPA | Arg55, Ile57, Phe60, Met61, Gln63, …Trp121, Leu122, His126, Arg148 | Arg55, Ile57, Phe60, Met61, Gln63, …Trp121, Leu122, His126 |

---

## 3b. 结果三：显式声明共价键，三元排布更贴近晶体 ✅

RMC‑6291 的弹头是 α,β‑不饱和酰胺（Michael 受体），对 **KRAS G12C 的 Cys12 共价加成**（9BFX 中 SG↔胶 C15 = 1.68 Å）。在 Boltz 输入 `bonds` 里声明 **Cys12 SG ↔ 胶弹头 β‑碳** 后（**T1cov**）：

- 胶↔Cys12 从 **4.66 Å（非共价）收紧到 1.96 Å**（≈ 共价 C–S 键长），胶↔KRAS 接触残基重新纳入 **Cys12、Ala11**（非共价版丢失的共价锚点）；
- 蛋白–蛋白 **BSA 337 → 417 Å²**，几乎等于实验 **425 Å²**；
- **CYPA 侧界面回收率大幅提升**：precision 67%→**90%**、recall 55%→**82%**（KRAS 侧维持在 89%/80%）；
- binding_confidence 0.903 → **0.931**。

> 直觉：共价锚点把胶的一端钉死在 Cys12，等于给三元复合物加了一个正确的"铆钉"，使 CYPA 被胶拉到的相对取向更接近天然——**约束越符合真实化学，预测越贴近晶体**。

---

## 4. 一句话结论

> **胶把一个"测不出的二元亲和力"变成了一个"高置信、可复现、与晶体一致的三元界面"。** Boltz 在含胶时稳定复现 9BFX 的 CYPA–KRAS neo‑PPI（KRAS 界面残基 100% 命中），去掉胶后界面置信度崩塌且变成非天然接触——正是正文 §7 所说的"分子胶诱导新生蛋白界面 + 正协同"的计算侧印证。

---

## 5. 诚实边界 ⚠️

- **这不是协同自由能（ΔΔG_coop）的定量。** 置信度/界面回收率能定性区分"胶模板化的天然界面" vs "硬对接的非天然界面"，但要给出正文 §8 要求的**协同系数 α**，仍需双体系 FEP（有/无伴侣蛋白）。本 demo 的 T1/T2 对照是**协同性的定性代理**，非金标准定量。
- **共价键已显式建模（T1cov）。** RMC‑6291 对 KRAS **G12C 是共价**（9BFX 中胶↔Cys12 = 1.68 Å）。T1（非共价）胶↔Cys12 落在 4.66 Å，故 KRAS 侧界面/BSA 略小于实验；在输入 `bonds` 里声明 **Cys12 SG ↔ 胶弹头 β‑碳**（加成后饱和形式，β‑碳标 `[CH:1]`）后的 **T1cov**，胶↔Cys12 收紧到 **1.96 Å**（≈ C–S 键），BSA 回到 **417 Å²（≈ 实验 425）**，CYPA 界面回收率从 55%→**82%**。见下节「结果三」。
- **BSA 为近似 SASA**（自实现 Shrake–Rupley，192 点），绝对值有系统误差；文中只用其做**相对比较**与趋势，不作绝对结论。
- 构建体用规范 KRAS/CypA（非 9BFX 的 cysteine‑light 结晶突变体），编号与前序 demo 一致；不影响界面定性结论。

---

## 6. 复现

```bash
ls structures/   # 9BFX.cif  pred_T1_ternary.cif  pred_T1cov_ternary.cif  pred_T2_noglue.cif
# 界面 + BSA + 胶桥接分析（完整输出见 tri_analysis.txt）
python3 tri_analyze.py structures/9BFX.cif             --kras A --cypa D --glue-comp A1AOD
python3 tri_analyze.py structures/pred_T1_ternary.cif    --kras A --cypa C --glue L
python3 tri_analyze.py structures/pred_T1cov_ternary.cif --kras A --cypa C --glue L
python3 tri_analyze.py structures/pred_T2_noglue.cif     --kras A --cypa C
```

Boltz 输入：`structure_and_binding`。
- **T1** = 2 蛋白 + 胶(SMILES) + GNP(ccd) + MG(ccd)，binding=ligand_protein_binding(binder=胶)。
- **T1cov** = 同 T1，胶弹头 β‑碳在 SMILES 中标 `[CH:1]`（加成后饱和形式），并加 `bonds: [{atom1: polymer_atom(chain A, residue_index 11=Cys12, atom_name "SG"), atom2: smiles_atom(chain L, atom_map 1)}]`。
- **T2** = 2 蛋白，binding=protein_protein_binding(binder=KRAS)。

胶 SMILES 见 `chem_A1AOD.json`（来自 RCSB，PDB 配体 A1AOD）；共价版 SMILES 把 `...C(=O)/C=C/C(C)(C)N(C)C` 改为 `...C(=O)C[CH:1]C(C)(C)N(C)C`。
