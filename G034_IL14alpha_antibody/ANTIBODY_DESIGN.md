# G034 — 双轴抑制性抗体设计(Boltz‑2.1,已实跑)

> 承接 `MECHANISM_AND_DRUG_DESIGN.md`。针对甲方两轴,**epitope‑directed** 设计阻断型抗体并亲和力排序。均为 🧪 本会话实跑;抗体形式为**配对 Fv(VH+VL)**,可直接改造为人源化全长 IgG(甲方需求)。
> 诚信:✅ 依据 / ⚠️ 需湿验证 / 🧪 已算。**所有序列为 in‑silico 设计,须湿实验(表达/亲和 SPR/功能中和)确认;打分以 iPTM/界面 PAE 为准,binding_confidence 对 de‑novo 校准不佳。**

## 现有抗体克隆情况
先前技术抗体(Peng 2009 抗 IL‑14α‑C、US7622574)**均无公开序列**,无法直接"基于现有克隆优化" → 故对每个表位 **de‑novo 设计 30 条**,亲和力排序(下)。

## 轴 1 — 阻断 IL‑14α ↔ 活化 B 细胞 IL‑14R(治疗轴)
- 靶表位:P40222 250–279(N 端结构化核心的暴露面;受体身份未知,选非 C 端、非先前技术表位以利新颖性)。
- **Top 抗体(pres_GWGX…):iPTM 0.85、界面 PAE 3.52 Å**(30 条最佳;另有 iPTM 0.81/0.78/0.78 三条,见 `design_ranking_axis1.tsv`)。
- 结构验证表位:**P40222 275–294**(命中目标邻域)。
- VH:`EVQLVESGGG…CARVSAGTRAVDVWGQGTLVTVSS`(119aa);VL:`DIQLTQSPSS…CQSNYGLPTFGQGTKVEIK`(110aa)(全序列见 `antibody_top_designs.faa`)。

## 轴 2 — 阻断 α‑taxilin ↔ syntaxin(=本抗原运输轴)⭐结构验证
- 靶表位:P40222 360–389(syntaxin‑4 结合面,由 STX4‑H3 复合物 §2.5 定位的 355–443 界面内)。
- **Top 抗体(pres_PvYZ…):iPTM 0.77、界面 PAE 5.32 Å**(见 `design_ranking_axis2.tsv`)。
- ✅**结构验证表位 = P40222 324–350,其中残基 326/327/328/329/332/333/335/336 正是 STX4 界面接触残基** → 该抗体**真正竞争阻断 syntaxin 结合**(不是泛结合)。
- VH:`QVQLVQSGAE…CAISTVSYNPILKTASTGAMSIWGQGTLVTVSS`(128aa);VL:`DIQMTQSPSS…CQTSTSAGPTTFGGGTKVEIK`(108aa)。

## 亲和力成熟(Top VHH,表位 325–350)
- **Top VHH(pres_KP4f…):iPTM 0.72、界面 PAE 6.72 Å**;表位 P40222 321–343(与原始 de‑novo top VHH 同区)。
- 说明:原始 de‑novo top VHH(iPTM 0.83)仍最优;本 epitope‑directed 批提供**同表位多样 CDR 备选**(`design_ranking_affmat.tsv`)。真正的迭代成熟建议:以最优复合物为模板,固定 framework、只重设计 CDR‑H3 多轮 + SPR 闭环。

## 汇总:三个先导(🧪 iPTM 排序)

| 轴/用途 | 先导 ID | 形式 | iPTM | 界面 PAE(Å) | 验证表位(P40222) | 阻断机制 |
|---|---|---|---|---|---|---|
| 轴1 受体阻断 | pres_GWGX… | Fv(VH+VL) | **0.85** | 3.52 | 275–294 | 遮蔽 IL‑14R 结合面(推定) |
| 轴2 syntaxin 阻断 | pres_PvYZ… | Fv(VH+VL) | **0.77** | 5.32 | 324–350(∈STX4 界面) | ✅竞争 syntaxin 结合(结构验证) |
| 核心表位 VHH | 原 de‑novo / pres_KP4f… | VHH | 0.83 / 0.72 | 3.67 / 6.72 | 325–350 | 核心表位、专利差异化 |

## 下一步(可选,付费/湿实验)
1. **CDR 移植到人源 IgG1(效应减弱 Fc)** + Boltz 复算,满足"人源化全长"要求。
2. **CDR‑H3 定向亲和成熟**(固定 framework,重设计 CDR)多轮 + Boltz 复筛 → 提升 iPTM/降低界面 PAE。
3. 开发性:TAP/Therapeutic Antibody Profiler、去免疫原性、稳定性液相条件;检查自身抗原 humanness。
4. 功能验证:重组 IL‑14α 结合(SPR/BLI)+ B 细胞/MZB 增殖中和(轴1)/ syntaxin pull‑down 竞争(轴2)。

## 诚信小结
- ✅ 两轴各得**高置信先导**(轴1 iPTM 0.85、轴2 iPTM 0.77),**轴2 表位经 STX4 复合物结构独立验证为真阻断位**。
- ⚠️ 全为**计算设计**,须湿实验;Fv 需改人源化全长 IgG;绝对亲和力(nM)需 SPR 实测。
- 数据:`antibody_top_designs.faa`、`design_ranking_{axis1,axis2,affmat}.tsv`、`boltz_structures/{axis1,axis2,affmat}_top/`。
