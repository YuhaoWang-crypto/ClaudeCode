# VC · 虚拟细胞与体外替身 — 产品家族规格 v2

把 [`VIRTUAL_CELL.md`](../VIRTUAL_CELL.md) 的虚拟细胞分层与实测能力边界，落成
**QureGen BioDecision** 产品体系里的一个新家族：按**细胞系**和**适应症**组织，
可直接并入现有站点。

v2 依据《虚拟细胞模型与代表性细胞系实施方案》技术报告做了七处修订，见 §8。

| | |
|---|---|
| 家族代号 | `VC` |
| 家族名称 | 虚拟细胞与体外替身 / Virtual Cell & In-Vitro Twins |
| 产品线 | 7 条 |
| 已定义方案 | 38 项 |
| 覆盖细胞系 / 模型 | 13 |
| 组织轴 | 细胞系 × 适应症 |
| 语言 | 中文 / English（同一数据源） |

机器可读规格：[`vc_family.json`](vc_family.json)（`schema: quregen-biodecision/family@2`，
所有文案字段均为 `{zh, en}` 双语对象）。

---

## 1. 五模块架构

产品边界来自模块划分。最关键的一条：**生物状态 ≠ 仪器读数**。预测"细胞发生了什么"和
预测"CellTiter-Glo 最终读到多少"是两个模型，中间那层必须用配对实验标签单独校准
——这就是 AssayEmul-AI 独立成线的原因。

| 模块 | 职责 | 产品线 |
|---|---|---|
| **Cell Passport** | 细胞系身份、基因组、基础转录状态、表观调控、基因依赖、药物历史、形态表型、实验上下文 | `VC-TWN` |
| **Perturbation Encoder** | 药物结构与靶点、CRISPRi / KO / CRISPRa、细胞因子、剂量、时间 | `VC-PRT` · `VC-PHE` |
| **State Transition** | 扰动后单细胞或群体状态如何迁移；DEG、通路分数、亚群组成 | `VC-PRT` · `VC-PHE` · `VC-RES` |
| **Assay Emulator** | 把生物状态映射到活力、荧光、图像、蛋白、分泌 readout——必须用配对标签校准 | `VC-ASY` |
| **Uncertainty & Domain Check** | 训练覆盖不足、跨细胞系外推、异常组合的识别；置信区间与适用域告警 | 内建于每条线 |

---

## 2. 细胞系底座矩阵

选细胞系的第一个问题不是"它像不像病人"，是"**它有没有能训练和验收的底座数据**"。
没有底座就没有噪声天花板，没有噪声天花板就没有验收标准，所有模型分数都无法解释。

- **扰动底座** — 系内或多背景 Perturb-seq
- **化合物底座** — L1000 签名 + 药敏（GDSC / PRISM）
- **形态底座** — Cell Painting（JUMP / BBBC021）

| 细胞系 | 背景 / 适应症 | 扰动 | 化合物 | 形态 | 阶段 | 可支撑的产品线 |
|---|---|---|---|---|---|---|
| A549 | 肺腺癌 · KRAS 相关背景 | 全 | 全 | 全 | 1 | TWN · PRT · PHE · ASY · CMB |
| MCF7 | 乳腺癌 ER+ | 部分 | 全 | 部分 | 1 | TWN · PHE · ASY · RES · CMB |
| HT29 | 结直肠癌 · 上皮 / EMT | 全 | 全 | 部分 | 1 | TWN · PRT · PHE · RES |
| BxPC-3 | 胰腺癌 · 上皮与应激 | 全 | 部分 | 无 | 1 | TWN · PRT · PHE |
| K562 | 髓系白血病 · 悬浮 | 全 | 部分 | 无 | 1 | TWN · PRT · TOX · RES · CMB |
| HAP1 | 近单倍体遗传筛选模型 | 全 | 无 | 无 | 1 | TWN · PRT |
| HepG2 | 肝细胞 · 肝毒性代理 | 部分 | 全 | 无 | 2 | TWN · TOX · ASY |
| RPE1 | 相对正常 · p53 完整上皮 | 全 | 无 | 无 | 2 | TWN · PRT · TOX（选择性窗口） |
| Jurkat | T 细胞背景 | 部分 | 无 | 无 | 2 | PRT · TOX |
| U2OS | 成骨肉瘤 · Cell Painting 主力 | 无 | 无 | 全 | 2 | PHE · ASY |
| H1 hESC | 人胚胎干细胞 | 部分 | 无 | 无 | 2 | TWN · RES（分化与命运） |
| PC9 / HCC827 | 肺腺癌 · EGFR 敏感突变 | 无 | 部分 | 无 | 扩展 | RES |
| 客户自有系 | 按项目 | 待评估 | 待评估 | 待评估 | — | 先做 CellTwin 可行性判定 |

> **两条不能混的线。** 跨背景 Perturb-seq 覆盖 A549 / MCF7 / HT29 / BxPC-3 / K562 / HAP1
> 六条系（GEO GSE281048），这是训练"相同基因扰动在不同背景下为何结果不同"的语料；
> 而化合物与形态底座最厚的是 MCF7 / A549 / U2OS / HepG2。两条线的重叠只有 A549 和 MCF7
> ——所以基因扰动产品和化合物表型产品各自的主力细胞系不同，**A549 是唯一四种底座都齐的系**。

---

## 3. 七条产品线

沿用现有体系的五步骨架（标准输入 → 候选生成 → 证据校准 → 组合选择 → 实验输出）
与七阶段项目流程。

| 代号 | 产品线 | 解决的问题 | 组织轴 | 方案 | 项目范围 |
|---|---|---|---|---|---|
| `VC-TWN` | **CellTwin-AI**｜Cell Passport 与可建模余量判定 | 这条系的基线、可建模余量和噪声天花板是多少 | 按细胞系 | 8 | 6万–14万 |
| `VC-PRT` | **PerturbSim-AI**｜基因扰动虚拟预筛 | 敲低哪些基因值得进湿实验 | 按通路 / 适应症 | 6 | 8万–18万 |
| `VC-PHE` | **PhenoScreen-AI**｜化合物表型虚拟预筛 | 这批化合物属于哪类 MoA、先测哪 384 孔 | 按筛选场景 | 5 | 8万–20万 |
| `VC-ASY` | **AssayEmul-AI**｜检测 readout 虚拟化 | 生物状态怎么变成仪器真正读到的那个数 | 按检测类型 | 5 | 10万–20万 |
| `VC-TOX` | **ToxTwin-AI**｜器官毒性虚拟前哨 | 哪些毒性表型要提前设计对照 | 按器官 | 5 | 8万–16万 |
| `VC-RES` | **ResistSim-AI**｜耐药、隐性状态与状态转换 | 离开关翻转还有多远；处理前是否已有隐性耐药克隆 | 按适应症 | 5 | 10万–22万 |
| `VC-CMB` | **CombiRank-AI**｜组合用药虚拟排序 | n² 组合里先测哪 96 孔 | 按组合场景 | 4 | 8万–18万 |

### 方案清单

**VC-TWN · CellTwin-AI**（地基，其余六条线的前置）
`001` A549 · `002` MCF7 · `003` HT29 · `004` BxPC-3 · `005` K562 · `006` HAP1 ·
`007` 二期扩展包（HepG2 / RPE1 / Jurkat / U2OS / H1）· `008` 客户自有细胞系

**VC-PRT · PerturbSim-AI**
`001` KRAS-MAPK 合成致死预筛 · `002` 双基因遗传相互作用（K562）·
`003` HAP1 KO 功能与逃逸机制 · `004` DDR / PARP 合成致死预筛 ·
`005` T 细胞激活通路敲低（Jurkat）· `006` 因果调控网络虚拟敲低

**VC-PHE · PhenoScreen-AI**
`001` MoA 归类与去重 · `002` 疾病签名逆转检索 · `003` 跨细胞选择性排序 ·
`004` 剂量-时间响应插值 · `005` 天然产物 / 提取物 MoA 分诊

**VC-ASY · AssayEmul-AI**
`001` 活力 readout 校准（CellTiter-Glo / MTT）· `002` 转录 readout 与通路分数 ·
`003` Cell Painting 形态 readout · `004` 流式表面标志物（需 CITE-seq 配对）·
`005` ELISA / 细胞因子分泌（需分泌组配对）

**VC-TOX · ToxTwin-AI**
`001` 肝毒性（HepG2）· `002` 心脏毒性（hiPSC-CM）· `003` 骨髓抑制（K562 / HL-60）·
`004` 肾小管毒性（HK-2 / RPTEC）· `005` 癌 / 正常选择性窗口（RPE1 对照）

**VC-RES · ResistSim-AI**（接 `grn_pipeline` 动力学引擎）
`001` EGFR-TKI 耐药（PC9 / HCC827）· `002` KRAS G12C 适应性耐药 ·
`003` 隐性耐药状态与克隆选择（谱系条码）· `004` EMT 可逆性与 hysteresis 判定 ·
`005` ER+ 内分泌治疗耐药（MCF7）

**VC-CMB · CombiRank-AI**
`001` KRAS G12C 上下游组合 · `002` CDK4/6 + 内分泌组合 ·
`003` DDR 合成致死组合 · `004` 双基因组合靶点（K562）

---

## 4. 能力边界（合同级）

这是本家族与市面其他"虚拟细胞"产品的主要区别：**第三档写进合同附件作为排除条款**。
划分依据见 [`VIRTUAL_CELL.md` §3](../VIRTUAL_CELL.md)。

### ✅ 直接替代
- 靶点在该系的表达、拷贝数、突变与融合状态（CCLE 直接查）
- 该靶点在该系是否为必需基因（DepMap CRISPR 直接查，比任何模型都准）
- 化合物结构警报：PAINS / Brenk / NIH 干扰子与类药性
- sgRNA 全基因组脱靶枚举与 CFD 特异性打分（穷举，不是预测）
- 细胞系身份与交叉污染核对（STR / SNP 指纹）

### 🟡 排序后测
- 化合物的转录响应方向与 MoA 归类（须先打赢 L1000 检索基线）
- 单基因敲低的下游通路——仅限该系有扰动底座覆盖的靶点及邻域
- **同一模态内的 CRISPRa 预测**（该系自有 CRISPRa 数据，如 K562 双基因集）
- 化合物–靶点结合排序（`opt_score` ρ≈+0.6；`binding_confidence` ρ≈−0.2）
- n² 组合空间的优先级排序与非加和效应候选
- 剂量与时间的插值（在已测剂量-时间网格内部）
- 网络失稳早期预警：方差与 lag-1 自相关在分岔点附近上升

### ❌ 不承诺
- 外推到没有扰动底座的新细胞系——oracle 上下文特征比不给还差
- **用敲低数据训练的模型跨模态预测 CRISPRa / 过表达**（`r = 0.009`，等于掷骰子）
- 绝对 IC50 与完整剂量–响应曲线，只排序不给数
- **未经配对标签校准的蛋白、分泌与活力 readout（RNA ≠ 功能）**
- **时间外推：6 小时转录响应不代表 72 小时活力或长期耐药**
- 从细胞系外推到患者响应；细胞系不等于患者
- 监管级毒理放行结论

---

## 5. 验证协议：七种数据拆分

随机 cell split **不能作为主要证据**——同一扰动和批次会泄漏，分数虚高。

| 数据拆分 | 优先级 | 回答的问题 |
|---|---|---|
| 随机 cell split | 不推荐 | 同一扰动和批次泄漏，容易虚高 |
| Leave-perturbation-out | 必做 | 测试未见基因或药物 |
| Leave-drug-out | 必做 | 确保测试集不含同一化合物 |
| Leave-cell-line-out | 必做 | 测试跨细胞背景泛化——本家族已知最弱的一轴 |
| Leave-study-out | 必做 | 检验跨实验室、平台和批次迁移 |
| Unseen dose / time | 推荐 | 检验剂量和时间的插值与外推 |
| Unseen combination | 推荐 | 检验药物或基因组合的非加和效应 |
| External wet-lab set | 上线前必做 | 用与训练数据独立的真实实验评估 |

**指标。** 不能只看整体表达 MSE 或 Pearson。至少同时报告：DEG 方向一致率、DE overlap、
Pearson delta、扰动辨别分数、通路 AUROC、细胞亚群距离、IC50/AUC 误差、top-k 富集、
置信区间覆盖率，以及实验重复给出的性能上限。

**最低上线原则。** 若深度模型未稳定超过 no-change、平均扰动和线性 / 岭回归基线，
该任务只可作为探索性分析，不能作为减少湿实验的主要依据。

---

## 6. 结算口径

不报相关系数，报**富集倍数**和**每个真阳性的实验成本**。行业可达区间是把大规模组合
压缩到优先验证的 **5%–20%**。示意算例：

```
初筛库容量                     10,000
历史真实命中率                 0.5%   → 50 个真阳性
虚拟预筛取头部                 500    (5%)
头部命中率（8× 富集）          4%     → 20 个真阳性
──────────────────────────────────────────────
实验量 ↓ 20×   │   召回 40%（丢掉 30 个真阳性）
```

划不划算取决于**一个 hit 的边际价值对一个孔的边际成本**。靶点富集的小库不值得做预筛；
十万级表型库、或动物实验前的排序，杠杆才起来。

---

## 7. 三阶段与停止条件

| 阶段 | 周期 / 算力 | 交付 | 停止条件 |
|---|---|---|---|
| **Phase 0** 盘点与地基 | 2–4 周，零 GPU | 八层 Cell Passport、噪声天花板表、五级基线阶梯、污染审计判定、两条检索基线 | 噪声天花板低到基线已贴顶 → 没有可建模余量，换 readout 或换系，不硬训 |
| **Phase 1** 主力模块 | 1–2 月，少量 GPU | 系内扰动模型（分层报未见靶点与未见细胞系）、化合物表型检索、Assay Emulator 首个 readout head、实验减量对账表 | 打不赢 no-change / 平均扰动 / 线性 / 最近邻检索 → 改架构，不调参硬凑，不换指标重报 |
| **Phase 2** 整合与闭环 | 2–3 月 | 顺式/反式并列输出、动力学预警引擎、不确定性分层与适用域告警、主动学习闭环、外部湿实验验证集 | 外部湿实验集上未能复现富集倍数 → 这层不上线，回到 Phase 1 |

推荐入口：任选一条细胞系的 **CellTwin-AI**。最低承诺、最快出结论，且结论无论正负都可复用。

---

## 8. v2 相对 v1 的七处修订

依据《虚拟细胞模型与代表性细胞系实施方案》技术报告（2026-08-31）：

| # | 修订 | 原因 |
|---|---|---|
| 1 | **细胞系矩阵修正**：新增 HT29、BxPC-3、HAP1、H1 hESC，并把扰动底座从"仅 K562/RPE1"更正为**六系跨背景 Perturb-seq**（GEO GSE281048） | v1 低估了公开扰动数据的细胞系覆盖，把 PerturbSim 的可服务适应症限窄了 |
| 2 | **新增 `VC-ASY` AssayEmul-AI 产品线** | v1 把"生物状态"和"仪器读数"混在一层。要模拟 MTT / CellTiter-Glo / ELISA / 流式，必须单独学习试剂、批次、时间与仪器的映射，且需要配对标签 |
| 3 | **采用五模块架构词汇**（Cell Passport / Perturbation Encoder / State Transition / Assay Emulator / Uncertainty & Domain Check） | 比 L0–L7 分层更适合作为商业产品分解；CellTwin 的交付物直接改名为 Cell Passport（八层） |
| 4 | **CRISPRa 边界精确化** | v1 笼统写"不承诺 CRISPRa"。实际失败的是**跨模态迁移**（CRISPRi 训练 → CRISPRa 预测，r=0.009）；该系自有 CRISPRa 数据时（如 K562 双基因集）属于可预筛档 |
| 5 | **验证协议从 3 种拆分扩到 7 种** + 明确随机 cell split 为反模式 | 补上 leave-drug-out、leave-study-out、unseen dose/time、unseen combination 与外部湿实验集 |
| 6 | **ResistSim 新增两项方案**：隐性耐药状态与克隆选择（谱系条码）、EMT 可逆性与 hysteresis 判定 | 谱系条码回答的是"耐药是演化出来的还是处理前就存在被选择的"，与动力学推演互补；EMT hysteresis 直接对接 `grn_pipeline` 的双稳/滞回引擎 |
| 7 | **新增边界条款**：RNA ≠ 功能、时间不可外推、许可与合规审查 | 商业交付前需完成模型与数据集许可证审查；部分模型或 API 仅允许非商业研究使用 |

---

## 9. 与现有产品体系的关系

| 现有产品 | 关系 |
|---|---|
| `SC-MOD-*` ModelMatch-AI｜细胞系与类器官选择 | **上游**。ModelMatch 选出模型组合后，CellTwin 判定这些模型有没有可建模余量 |
| `SC-AUD-003` 扰动预测模型审计 | **姊妹件**。ModelAudit 审别人的模型，PerturbSim 交付自己的——共用同一套基线阶梯与七种拆分 |
| `SC-SAF-*` SafetyMap-AI｜靶点正常组织安全图谱 | **互补**。SafetyMap 看组织表达，ToxTwin 看体外毒性表型 |
| `MOL-FOC-*` FocusLibrary-AI｜靶点聚焦采购库 | **下游**。PhenoScreen 的富集结果可直接驱动聚焦库采购 |
| `IMG-*` 生物影像家族 | **共用 Assay Emulator**。Cell Painting readout head 两边都用 |
| `grn_pipeline`（本仓库） | **ResistSim-AI 的引擎**。CRNT 亏格判开关能力、Markevich 双稳窗口、临界慢化预警与 M18 滴定阳性对照 |

---

## 10. 文件

| 文件 | 内容 |
|---|---|
| `vc_family.json` | 机器可读双语规格（`{zh, en}`），可直接喂站点生成器 |
| `vc_family_zh.html` | 离线单文件产品图谱 · 中文版（50 路由，含语言切换） |
| `vc_family_en.html` | 离线单文件产品图谱 · English（同一数据源） |
| `README.md` | 本文档 |
| [`../VIRTUAL_CELL.md`](../VIRTUAL_CELL.md) | 底层技术依据：七层模型地图、公开数据盘点、12 条实测能力边界 |

两个 HTML 版本共用同一份数据结构，仅 `LANG_DEFAULT` 常量不同；页面右上角可切换语言。
`vc_family.json` 的每条产品线包含 `steps`（五步骨架）、`phase_text`（七阶段项目流程的
可填充文本）、`engagement`（合作交付表，含三档能力边界）与 `solutions`（已定义方案）。

## 参考资料

- Arc Institute — STATE / Stack / 2026 Virtual Cell Challenge
- Tahoe Therapeutics — Tahoe-100M
- Roohani Y, et al. GEARS. *Nat Biotechnol*
- Ahlmann-Eltze C, et al. Deep-learning perturbation models versus linear baselines. *Nat Methods*
- Viñas Torné R, et al. Systema evaluation framework. *Nat Biotechnol*
- DepMap · CMap/LINCS · JUMP Cell Painting · CZ CELLxGENE Census · AlphaGenome
- NCBI GEO **GSE281048** — Multi-context Perturb-seq

完整链接见 `vc_family.json` 的 `references` 字段与图谱页面的参考资料页。
