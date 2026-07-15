# G034 — 共价 PROTAC 设计(靶向 α‑taxilin C245/C373 + CRBN/VHL)

> 目标:α‑taxilin 无口袋、无配体(ChEMBL 仅 1 条 3.5 µM kinobead 脱靶),经典可逆 PROTAC 不可行。**策略 = 共价电泳体优先(electrophile‑first)**:用 Cys 反应弹头共价锚定 Boltz 判定的**高置信半胱氨酸 C245 / C373**(位于折叠螺旋核心、半暴露),再经 linker 连 E3 配体招募 CRBN 或 VHL → 泛素化降解。🧪 SMILES 已 RDKit 验证。
> 诚信:✅ 有依据 / ⚠️ 需实验验证 / 🧪 本会话已算。

## 1. 组件

| 模块 | 选择 | 说明 |
|---|---|---|
| **锚定位点** | **C245 或 C373**(P40222) | Boltz 折叠:两者 pLDDT 88–90、半暴露、位于结构化 coiled‑coil 核心(优于埋藏的 C471 与柔性尾的 C523) |
| **弹头 warhead** | **丙烯酰胺**(Michael 受体,温和、可逆性可调)/ **氯乙酰胺**(SN2,更强) | 靶向 Cys 硫醇;丙烯酰胺是临床首选(如 KRAS G12C、BTK) |
| **linker** | PEG2/PEG3 或 C4 烷基 | 长度/柔性影响三元复合物几何,需 SAR 优化 |
| **E3 配体** | **CRBN = 泊马度胺**(MW273)/ **VHL = VH032**(MW430) | 两条并行,泛素连接酶组织表达不同,择优 |

## 2. 候选分子(🧪 RDKit 验证通过)

| ID | E3 | 弹头 | linker | MW | cLogP | HBD/HBA | TPSA | 可旋转键 |
|---|---|---|---|---|---|---|---|---|
| **G034‑PT1** | CRBN | 丙烯酰胺 | PEG2 | 472.5 | −0.64 | 3/8 | 160 | 11 |
| **G034‑PT2** | CRBN | 氯乙酰胺 | C4 | 448.9 | 0.55 | 3/6 | 142 | 8 |
| **G034‑PT3** | CRBN | 丙烯酰胺 | PEG3 | 516.5 | −0.62 | 3/9 | 169 | 14 |
| **G034‑PT4** | VHL | 丙烯酰胺 | PEG2 | 629.8 | 1.56 | 4/9 | 159 | 15 |

**SMILES:**
```
G034-PT1  C=CC(=O)NCCOCCOCC(=O)Nc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O
G034-PT2  ClCC(=O)NCCCCC(=O)Nc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O
G034-PT3  C=CC(=O)NCCOCCOCCOCC(=O)Nc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O
G034-PT4  C=CC(=O)NCCOCCOCC(=O)N[C@@H](C(C)(C)C)C(=O)N1C[C@H](O)C[C@H]1C(=O)NCc1ccc(-c2scnc2C)cc1
```
> 说明:CRBN 版 MW 449–517、TPSA 高(泊马度胺酰亚胺贡献),cLogP 偏低 → 溶解性好但被动膜透性需关注;VHL 版 MW 630、cLogP 1.6,更典型 PROTAC 空间。PROTAC 属 "beyond‑Rule‑of‑5",上表仅供 linker/E3 取舍参考,非 Ro5 合格判定。

## 3. 计算/实验验证路径

1. 🔌/付费:**共价对接** —— 因当前 Boltz ligand 建模对**共价键 + SMILES 弹头**支持有限(atom‑level 共价仅支持 CCD 配体),建议用 covalent docking(如 CovDock/ICM)把弹头锚到 C245/C373,评估弹头朝向与可及性。
2. 付费可跑:**Boltz 小分子 ADME(adme‑v1)** 对 4 个 SMILES 出 Tier‑1 ADME;**Inductive Bio** 出 logD/pKa。
3. 湿实验:①重组 α‑taxilin 与弹头片段做**质谱共价加合**(确认 C245/C373 选择性、动力学 k_inact/K_I);②三元复合物(α‑taxilin–PROTAC–CRBN/VHL)pull‑down;③细胞内**降解 DC50/Dmax**、Hook effect;④选择性蛋白组(避免脱靶 Cys 蛋白)。
4. ❗**安全性**:α‑taxilin 泛表达且为看家运输蛋白 → 全身降解有 on‑target 毒性风险;考虑**组织/细胞靶向递送**或与 2.4 剪接选择性策略并行。

## 4. 与专利
共价 PROTAC + 特定 Cys(C245/C373)锚定 + 特定 E3 → **全新物质与 MoA 专利空间**,且绕开 US7622574(仅覆盖抗体 + RNA 抑制剂,不含降解剂/小分子)。
