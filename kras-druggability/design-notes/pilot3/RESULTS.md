# Pilot 3 结果：STAT3 的"降解分支"——CRBN 诱导邻近三元（看赖氨酸递呈）

> 承接 Pilot 2：遮盖分支（CypA/FKBP12）抓不住 STAT3（binding_conf ~0）。本 pilot 验证**降解分支**是否成立。
> 判据与遮盖分支**完全不同**：不看静电宽面互补，而看**三元能否形成 + STAT3 赖氨酸是否被递到泛素化可及区**。
> 工具：Boltz‑2.1 + [`proximity_analyze.py`](proximity_analyze.py)。运行日期 2026‑11‑13。

---

## 1. 体系

真实的 STAT3 PROTAC 降解剂 **SD‑36**（CHEMBL4452527；Bai 组 2019 已验证的强效 STAT3 降解剂）——双功能：
- 一端 = **pTyr 模拟的二氟膦酸**（结合 STAT3 SH2 的 Arg609 磷酸袋，即 SI‑109 弹头）
- 中间 = 烷基/炔连接段
- 另一端 = **来那度胺谷氨酰亚胺**（结合 CRBN 的 tri‑Trp 口袋）

Boltz 三元：**CRBN‑TBD(人源 Q96SW2, 319–442) + SD‑36 + STAT3‑SH2(P40763, 580–715)**，ligand_protein_binding(binder=SD‑36)，num_samples=3，prediction_id `sab_pred_8EIucgoc9aBesC5IRQQi`。

## 2. 结果：诱导邻近成立，且置信度高

**置信度（与 Pilot 2 遮盖分支对比，差异巨大）：**

| 体系 | 分支 | binding_confidence | protein_iptm |
|---|---|---|---|
| CypA–STAT3‑SH2（Pilot 2） | 遮盖 | **0.0002** | 0.78(波动) |
| FKBP12–STAT3‑SH2（Pilot 2） | 遮盖 | **0.019** | 0.79(波动) |
| **CRBN : SD‑36 : STAT3‑SH2（本 pilot）** | **降解** | **0.871** | **0.856** |

> 给 STAT3 一个**双功能 PROTAC**招募 CRBN 后，三元置信度从 ~0 跳到 **0.87**（ligand_iptm 0.97、ipde 0.73）——**降解分支对 STAT3 成立，遮盖分支不成立**，与 Pilot 2 的判定完全一致。

**[1] PROTAC 双功能桥（两端都锚对）：**
- 谷氨酰亚胺端 → **CRBN tri‑Trp 口袋 2.84 Å**（IMiD 识别袋，锚对 ✓）
- 膦酸(pTyr模拟)端 → STAT3‑SH2：25 残基接触、最近 2.47 Å
- CRBN↔PROTAC 16 残基；STAT3↔PROTAC 25 残基 → **双端同时锚定 ✓**

**[2] 诱导出新的 CRBN–STAT3 界面：** 33 / 26 残基接触，284 原子对，最近 2.43 Å。

**[3] STAT3 赖氨酸递呈（降解分支的核心判据）：**
STAT3‑SH2 构件含 11 个 Lys；朝 CRBN 一侧、落在泛素化可及区的最近两个是
**Lys707 / Lys709（距 CRBN 仅 6.1 Å）**，25 Å 内还有 K658(8.7)、K591(10.2)、K631(11.9)、K679(14.5)、K626(14.6)、K642(16.1)、**K685**(24.4) 等。
> 即三元几何**把 STAT3 的多个赖氨酸（含 pTyr705 尾部的 K707/K709）递到了 CRBN 面前**——有充足的泛素化"把手"。
> （构件按 580–715，Boltz 重编号为 1–136；真实 STAT3 号 = 构件号 + 579。）

## 3. 结论

**STAT3 确实能被这套"招募 presenter 诱导邻近"范式攻克——但走的是降解分支（PROTAC + E3/CRBN），不是 CypA 式遮盖分支。**
- Pilot 2 用遮盖判据（静电互补）正确地**否掉**了遮盖分支；
- Pilot 3 用降解判据（三元 + 赖氨酸递呈）**证实**了降解分支：高置信三元、PROTAC 双端锚对、STAT3 赖氨酸递呈到 CRBN（K707/K709 @6 Å）。
- 这正是 SD‑36 在实验上确为强效 STAT3 降解剂的结构解释，也闭环了"**预筛先判该靶适配哪条分支，再按该分支的判据建模**"的设计流程。

## 4. 诚实边界 ⚠️

- **只建了 CRBN‑TBD，没有 DDB1/CUL4/RBX1/E2~Ub 全酶**，所以"可及区<25Å"是**几何代理**，不是催化胜任度的证明；真正可降解还要求赖氨酸能够到达全 CRL4 组装中 E2~Ub 的活性位（需全复合物/MD/实际降解数据）。SD‑36 的实际强效降解活性提供了这一层的实验背书。
- **PROTAC 三元预测本身是公认难题**（连接段构象系综、协同性）；这里是单一高置信姿势，方向性强但非终局；连接段长度/刚性是真实设计变量。
- STAT3 用 **SH2 域片段**（非全长）；CRBN 用 TBD（非全长）。构件重编号见上。
- 赖氨酸递呈是**诱导邻近的必要几何条件之一**，非充分条件。

## 5. 复现

```bash
python3 proximity_analyze.py structures/pred_crbn_sd36_stat3.cif --crbn C --target S --protac L
```
Boltz 输入：CRBN‑TBD(Q96SW2 319–442) + SD‑36(CHEMBL4452527 SMILES) + STAT3‑SH2(P40763 580–715)，ligand_protein_binding(binder=L)。
