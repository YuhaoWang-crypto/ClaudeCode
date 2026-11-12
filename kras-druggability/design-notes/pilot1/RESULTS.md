# Pilot 1 结果：复现"为什么选 CypA"——presenter 界面静电互补的量化

> 配套：设计 playbook [`../molecular-glue-design-playbook.md`](../molecular-glue-design-playbook.md) §3.1 / §5。
> 目标：把视频里"选 CypA 的一个具体理由：局部表面电荷互补；FKBP12 不具同样互补"这个**设计决策**，用我们的结构+界面框架**量化、可复现**。
> 工具：[`electro_complementarity.py`](electro_complementarity.py)（纯 numpy）；Boltz‑2.1。运行日期 2026‑11‑12。

---

## 1. 做了什么

- **严谨层**：在**真实 9BFX** 与我们**预测的 T1cov** 的 CypA–KRAS 界面上，量化电荷互补（界面残基净电荷、跨界面盐桥/同电荷冲突、静电互补分 EC）。
- **对照层（有确认偏差，已标注）**：新跑 **FKBP12 + KRAS G12C（apo 二元）** Boltz 预测，和已有 **CypA + KRAS（apo，T2）** 同口径对比。
  - ⚠️ 为什么不能直接拿胶试 FKBP12：RMC‑6291 的 presenter 端源自 sanglifehrin A，**只结合 CypA 不结合 FKBP12**，带胶会变成"胶不结合"的混淆。故对照层用 **apo 二元**，只比**界面互补性**，不比结合与否。

## 2. 结果

| 体系 | binding_confidence | 界面净电荷 presenter / KRAS | 跨界面盐桥(<4Å) | 同电荷冲突(<4Å) | 静电互补 EC | 严谨度 |
|---|---|---|---|---|---|---|
| **9BFX CypA–KRAS（真实）** | — | **+5 / −5** | **8** | **0** | **+1.15** | ✅ 真实结构 |
| **T1cov CypA–KRAS（预测三元）** | 0.931 | +4 / −4 | 3 | 0 | +1.06 | ✅ 预测复现 |
| CypA–KRAS（apo, T2） | 0.317 | +2 / −3 | 2 | 0 | +0.38 | ⚠️ 低置信 apo |
| **FKBP12–KRAS（apo）** | **0.171** | +2 / **−1** | 3 | **2** | **+0.11** | ⚠️ 低置信 apo |

## 3. 判读

**✅ 严谨结论（基于真实结构）——视频的核心论断成立且可量化：**
- 天然 CypA–KRAS 界面是**强电荷互补**的：**KRAS 效应面偏酸（界面净电荷 −5）↔ CypA 面偏碱（+5）**，形成 **8 个跨界面盐桥、0 个同电荷冲突**，EC=+1.15。
- CypA 侧的碱性界面残基正是 **Arg55 / Arg69 / Lys125 / Arg148 / Lys151**（与我们前面界面分析一致），对上 KRAS Switch‑I/II 的酸性簇（Asp33/Glu37/Asp38/Glu63…）。
- 我们**预测的三元（T1cov）复现了这一互补**（+4/−4，EC+1.06，0 冲突）——说明该度量在预测结构上也稳。

**⚠️ 对照结论（方向性，非证明）——与"为什么不是 FKBP12"一致：**
- FKBP12 被迫与 KRAS apo 对接时：① **binding_confidence 最低（0.171 < CypA apo 0.317）**；② 界面被摊成一大片低质量接触（25/22 残基 vs 天然 ~10/~11）；③ **首次出现同电荷冲突（2 个）**；④ 静电互补几乎为零（EC +0.11）、KRAS 侧只凑到 −1。
- 即：**在相同口径下，FKBP12 找不到像 CypA 那样干净的电荷互补界面**——方向上支持"选 CypA 不选 FKBP12"。

> 📌 对设计的意义：**"presenter 选择"这一步可以被量化**——把候选 presenter 对目标靶面的"界面净电荷互补 + 盐桥/冲突计数 + EC"算出来排序，而不用凭感觉。这就把 playbook 准则 B③（界面可诱导互补）变成一个可计算的筛选器。

## 4. 诚实边界 ⚠️

- **最可信的是"真实 CypA 界面强互补"**（来自 9BFX 晶体）；FKBP12 一腿是**低置信 apo 对接上的度量**，位姿本身不确定，只能作**方向性**证据，不是严格证明。要严格比较需 FKBP12 的三元几何（现实中不存在）或实验。
- EC 是**概念性部分电荷 + 距离加权**的相对度量，非精确泊松–玻尔兹曼静电；只用于**同口径相对比较**。
- 视频原图亦标注"电荷分布为概念示意"、"电荷互补 ≠ 两者会自行稳定结合"——与我们结论一致：互补是**可成胶的必要有利条件**，真正把三元钉住仍要靠胶（见 `../../demo/tricomplex/`：无胶时界面置信度崩塌）。

## 5. 复现

```bash
ls structures/   # pred_fkbp12_kras_apo.cif  (CypA 结构复用 ../../demo/tricomplex/structures/)
S=../../demo/tricomplex/structures
python3 electro_complementarity.py $S/9BFX.cif             --presenter D --kras A
python3 electro_complementarity.py $S/pred_T1cov_ternary.cif --presenter C --kras A
python3 electro_complementarity.py $S/pred_T2_noglue.cif     --presenter C --kras A
python3 electro_complementarity.py structures/pred_fkbp12_kras_apo.cif --presenter F --kras A
```

Boltz 输入：FKBP12(人源 P62942, 1–108) + KRAS G12C(1–169)，protein_protein_binding(binder=KRAS)，num_samples=3，prediction_id `sab_pred_DWhyDbVoUtC481zPDkAJ`。
