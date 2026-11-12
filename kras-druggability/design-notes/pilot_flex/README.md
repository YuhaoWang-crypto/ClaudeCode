# Pilot-flex：CypA 在三元形成中到底动了多少？

测 apo CypA（PDB **2CPL**）vs 结合态 CypA（**9BFX** 链 D）的构象变化，回答"presenter 是靠柔性适配吗"。

**结果**（`python3 rmsd.py`）：
- 全局 Cα RMSD = **0.41 Å**（晶体噪声级）→ CypA 整体**刚性**，不是柔性适配器。
- 界面 rim 偏差均值 0.58 Å，个别残基在动：**Asn71 2.15 Å、Arg69 1.39 Å**，其余 <0.7 Å → **局部诱导契合**，非整体变形。

**结论**：适配新界面主要靠**重新设计小分子（胶）**，presenter 提供的是"刚性锚 + 潜在互补面 + 少量 rim 可塑"。详见 [`../molecular-glue-design-playbook.md`](../molecular-glue-design-playbook.md) §3B。

> ⚠️ 单结构叠合反映的是**两个晶体态之差**，非真实运动轨迹/时间尺度；要谈动力学需 HDX‑MS / MD（视频原图亦如此标注）。
