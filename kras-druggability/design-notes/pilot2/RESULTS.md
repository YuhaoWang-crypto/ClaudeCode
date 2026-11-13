# Pilot 2 结果：把这套思路试在 STAT3——presenter 预筛 + 判定"适配哪条分支"

> 配套：playbook [`../molecular-glue-design-playbook.md`](../molecular-glue-design-playbook.md) §3B / §5。
> 目标：STAT3（转录因子、SH2 介导二聚化、直接抑制失败的典型）能否套用"招募 presenter 造界面"？
> 工具：真实结构刻画（PDB **1BG1**，活化态 pTyr705 STAT3 二聚体）+ Boltz‑2.1 对接 + [`../pilot1/electro_complementarity.py`](../pilot1/electro_complementarity.py)。运行日期 2026‑11‑13。

---

## 1. 严谨层：STAT3 的"功能面"是什么性质？

从真实 STAT3 二聚体（1BG1）解析二聚界面：

- 每侧界面 **26 残基、净电荷仅 +1（近中性）**——**不是** KRAS 那种宽阔酸性效应面（净 −5）。
- 界面核心是 **pTyr705（−2）↔ 对侧 SH2 Arg609（+）的磷酸"锁孔"**（两者都确认在界面上）。
- 即 STAT3 的功能面是 **"单磷酸热点 + SH2 口袋"** 型，而非"大片带电平面"。

> 设计含义：这种"热点+口袋"界面，**天然适合磷酸肽模拟物/口袋结合剂**，而**不适合刚性 presenter 的"宽面遮盖"**（CypA 对 KRAS 靠的是宽面电荷/形状互补）。

## 2. 预筛层：CypA / FKBP12 能抓住 STAT3‑SH2 吗？（Boltz 对接，同 pilot1 口径）

| 体系 | binding_confidence | 界面净电荷 presenter/靶 | 盐桥 | 冲突 | EC |
|---|---|---|---|---|---|
| **基线** CypA–KRAS 共价三元（可成胶） | **0.93** | +4 / −4 | 3 | 0 | +1.06 |
| 基线 CypA–KRAS apo（无胶） | 0.317 | +2 / −3 | 2 | 0 | +0.38 |
| **CypA–STAT3‑SH2 apo** | **0.0002** | +4 / −3 | 5 | 0 | +0.64 |
| FKBP12–KRAS apo（pilot1） | 0.171 | +2 / −1 | 3 | 2 | +0.11 |
| **FKBP12–STAT3‑SH2 apo** | **0.019** | +2 / −1 | 2 | 1 | +0.10 |

**读出：**
- **没有一个 presenter 对 STAT3‑SH2 形成可信复合物**：binding_confidence **0.0002（CypA）/ 0.019（FKBP12）**——比 CypA–KRAS **apo** 的 0.317 还低约 **1000× / 17×**，比可成胶的三元（0.93）低约 4600× / 48×。Boltz 实质判定"抓不住"。

## 3. 一个重要的方法学教训：别只看 EC（呼应正文 §9）

- CypA–STAT3 的 **EC=+0.64 竟高于 CypA–KRAS apo 的 +0.38**——因为**任意一个对接姿势都可能凑出几对电荷**，EC 单独看会**误导**。
- 真正的判别信号是 **binding_confidence（校准过的"是否真结合"）骤降到 ~0** + **靶面生物学（口袋热点 vs 宽平面）**。
- → **presenter 筛选要"binding_confidence + 界面几何/生物学"联合判读，不能只凭一个静电分**。这正是正文 §9"不可信的单一打分 vs 可信的联合判据"的又一实证。

## 4. 设计判定（Pilot 2 的真正价值：告诉你该走哪条分支）

**STAT3‑SH2 不适合"CypA 式刚性 presenter 遮盖"分支。** 原因：功能面是磷酸锁孔热点、近中性宽面，刚性 presenter 无宽阔互补面可抓（预筛 binding_confidence ~0 实证）。

STAT3 更适合这套范式的**另外两条路**：
1. **降解分支（E3 presenter）**：招募 **CRBN / VHL** 等 E3 → 泛素化降解 STAT3。诱导邻近做降解**不要求宽面互补**，只要能形成生产性三元并把一个 Lys 递给 E2——对"热点+口袋"型靶点更宽容（STAT3 降解剂本就是活跃方向）。⚠️ 注意：本 pilot 的"遮盖/静电互补"读出**不适用**于降解分支（几何判据不同），需换 E3 三元 + 赖氨酸可及性的建模。
2. **口袋直给**：SH2 磷酸口袋的**磷酸肽模拟物/口袋结合剂**（直接堵 pTyr705 识别），或**定制/工程化 presenter**（表面专为 STAT3 面设计）。

> 一句话：**预筛不是告诉你"能不能成药"，而是告诉你"这套范式里哪条分支适合这个靶"**——对 STAT3，答案是"走降解或口袋，不走 CypA 遮盖"。

## 5. 诚实边界 ⚠️

- apo 蛋白‑蛋白预测本就低置信（无胶时 CypA–KRAS 也才 0.317）；**有意义的是同口径的数量级对比**（CypA–KRAS apo 0.317 vs CypA–STAT3 0.0002 = 1000×），不是绝对值。
- binding_confidence 在这种极低区间可能有噪声，按**数量级方向**解读。
- 用的是 STAT3 **SH2 域片段（人源 580–715）**，非全长语境；1BG1 为活化态二聚体。
- 只测了**遮盖型** presenter（CypA/FKBP12）；**降解分支（E3）未建模**，因为静电/遮盖读出不是其正确判据。

## 6. 复现

```bash
# 真实 STAT3 二聚界面刻画：见本目录脚本注释；结构 structures/1BG1_assembly.pdb
# presenter 预筛（EC + 界面）：
python3 ../pilot1/electro_complementarity.py structures/pred_cypa_stat3.cif   --presenter P --kras S
python3 ../pilot1/electro_complementarity.py structures/pred_fkbp12_stat3.cif --presenter P --kras S
# 完整输出见 complementarity_analysis.txt
```

Boltz 输入：presenter(CypA P62937 1–165 / FKBP12 P62942 1–108) + STAT3‑SH2(人源 P40763, 580–715) + protein_protein_binding(binder=STAT3)，num_samples=3。
prediction_id：CypA `sab_pred_qTKAGBbu8dNTs53f5upg`；FKBP12 `sab_pred_RaO4d5papanJyJLW9rco`。
