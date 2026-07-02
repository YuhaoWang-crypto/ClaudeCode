# 面向单一适应症的新 Biomarker 挖掘策略

> 目标：针对**一个适应症**，整合 [NVIDIA BioNeMo recipes](https://github.com/NVIDIA-BioNeMo/bionemo-recipes) 的生物基础模型、蛋白组/蛋白层建模能力（下称 PROTO）、本环境已连接的生物医学数据源，以及 Claude 的编排与推理能力，系统性地挖掘**新的蛋白 / 基因 biomarker**，覆盖：疾病风险预测、诊断、疗效评估、治疗程度监测、入组富集、伴随诊断（CDx）、疾病进程分析。

---

## 1. 先把"biomarker"讲清楚：按临床用途分类（BEST 框架）

你列的每一个用途，都对应 FDA-NIH **BEST**（Biomarkers, EndpointS, and other Tools）里的一个标准类别。先对齐语言，后面所有工作都挂在这张表上：

| 你的说法 | BEST 类别 | 回答的临床问题 | 典型判读 |
|---|---|---|---|
| 疾病预测（发病风险） | **Susceptibility / Risk** | 这个人未来会不会得病？ | 高/低风险分层 |
| 疾病预测（是否已患病） | **Diagnostic** | 现在有没有病、是哪个亚型？ | 阳性/阴性、分子分型 |
| 方便病人入组 / 伴随诊断 | **Predictive** | 这个病人对某疗法会不会有效？ | 富集入组、CDx |
| 疾病治理效果 | **Response（含 PD 药效学）** | 药起没起作用（机制层面）？ | 早期应答信号 |
| 治疗程度分析 | **Monitoring** | 治疗中疾病负荷怎么变？ | 纵向趋势 |
| 进程分析 | **Prognostic** | 不管治不治，病会怎么发展？ | 进展/复发风险 |
| （安全性） | **Safety** | 会不会出现毒性？ | 停药/减量信号 |

**关键洞察：同一个分子实体（一个基因或蛋白），在不同临床语境下可以扮演不同角色。** 例如 PD-L1 既是 predictive（选免疫治疗人群）又能做 monitoring。所以挖掘时不要问"这是不是 biomarker"，而要问"这个候选在哪个 BEST 类别里、和谁比、在什么样本里测得到"。

---

## 2. 能力地图：每个模型/工具回答哪一类 biomarker 问题

把手上的"武器"按**数据模态**分层，每一层天然对应不同的 biomarker 类别。

### 2.1 基因组 / DNA 层 → 主打 Risk、Diagnostic
- **OpenGenome2 / Evo2 风格基因组 FM**（`opengenome2_llama_native_te` recipe）：对**非编码变异**（启动子、增强子、剪接位点）打功能影响分，找遗传易感性位点。
- **CodonFM**（`codonfm_native_te`）：密码子层面的序列建模，评估同义/近同义变异对表达和翻译的影响。
- **ESM-2 / AMPLIFY**（`esm2_native_te`、`amplify`）做**编码变异效应预测（VEP）**：用掩码语言模型的伪对数似然（pseudo-log-likelihood）给错义变异打"有害性"分，无需标签的零样本打分。
- 产出：**susceptibility / risk biomarker**（胚系变异面板）、部分 **diagnostic**（体细胞驱动突变分型）。

### 2.2 转录组 / 单细胞层 → 主打 Diagnostic、Prognostic、Response
- **Geneformer**（`geneformer_native_te_mfsdp_fp8`）：单细胞 RNA 的上下文感知基因网络嵌入。三个杀手锏：
  1. **细胞状态嵌入**：把病人 scRNA-seq 映射到潜空间，找疾病 vs 健康、应答 vs 不应答的可分离信号 → diagnostic / predictive。
  2. **in-silico 扰动（in-silico perturbation）**：在模型里"敲除/激活"某基因，看细胞状态是否被推向健康态，定位**驱动性**基因而非旁观者 → 机制性 biomarker，抗噪。
  3. **网络中心性变化**：找在疾病态里网络地位显著改变的基因 → progression / prognostic。
- 产出：**基因表达 signature**（诊断、预后、疗效应答）。

### 2.3 蛋白层 / PROTO → 主打 Response（PD）、Monitoring、CDx
- **蛋白语言模型（ESM-2 / AMPLIFY）嵌入**：对候选蛋白做功能/家族聚类、把序列变异映射到功能改变。
- **蛋白组差异分析（PROTO）**：血浆/组织蛋白丰度、翻译后修饰（PTM）signature——**这是最贴近"可及性"的一层**，因为血液蛋白最容易做成纵向可测的 monitoring / response 标志物。
- 产出：**可测量的循环蛋白 biomarker**（PD、monitoring、CDx 首选）。

### 2.4 结构 / 机制验证层 → 给候选"背书"
- **Boltz（结构 + 结合亲和力预测）**：候选蛋白/变异是否改变结构或药物结合口袋？为 predictive/CDx 提供机制证据（"为什么这个突变让药无效"）。
- **EDEN（免疫原性预测）**：候选是否涉及免疫识别、新抗原、免疫相关 AE → 免疫治疗场景的 predictive / safety。

### 2.5 药物 / 临床证据层 → 把候选和"可用性"接起来
- **ChEMBL**（target/bioactivity/mechanism/drug/ADMET）：候选靶点有没有已知配体？可成药性如何？→ 判断 biomarker 能否配对一个疗法做 CDx。
- **Inductive Bio**：小分子性质预测，支持 CDx 配对药物的成药性评估。
- **Clinical Trials（ClinicalTrials.gov）**：现有试验用什么 biomarker 做入组/终点？竞争格局？→ 定义"新"的基线、设计入组策略。
- **PubMed / bioRxiv**：候选的先验证据、新颖性核查（是不是真的"新"）、正交实验佐证。

### 2.6 Claude 的角色：编排 + 推理 + 证据三角
Claude 不是又一个预测模型，而是把上面所有层**串成一条可解释的流水线**：
- 生成假设、跨来源三角验证（PubMed + bioRxiv + ChEMBL + Trials）；
- 解读嵌入/打分、对候选排序、做特征归因叙述；
- 核对临床试验终点与入组标准、评估 CDx 可行性；
- 产出可审计的候选清单与证据卡片。

---

## 3. 端到端流水线（6 个阶段）

```
Phase 0  定义问题        →  Phase 1  多组学候选生成
   ↓                              ↓
Phase 5  组合面板 & 决策  ←  Phase 4  临床落地筛
   ↑                              ↑
Phase 3  机制验证         ←  Phase 2  差异/对比发现
```

### Phase 0 — 锁定适应症与临床问题
- 用 **Clinical Trials MCP** 拉该适应症在研试验，`analyze_endpoints` 看主流终点、`search_by_eligibility` 看入组用了哪些 biomarker。
- 用 **PubMed / bioRxiv** 建立已知 biomarker 基线与空白点。
- **明确要打哪个 BEST 类别**（不要一次全做，先挑 1–2 个最有临床价值的缺口）。
- 定义**对比组**：这是全流程的灵魂——responder vs non-responder、progressor vs stable、disease vs healthy。没有清晰对比组，挖出来的都是噪声。

### Phase 1 — 多组学候选生成（并行三条腿）
- 基因组腿：Evo2/OpenGenome2 + CodonFM + ESM-2 VEP → 变异候选。
- 转录组腿：Geneformer 嵌入 + in-silico 扰动 → 基因状态候选。
- 蛋白腿：PROTO 差异蛋白/PTM + ESM-2 嵌入 → 蛋白候选。
- 三条腿各自产出一个排序候选池。

### Phase 2 — 差异 / 对比发现
- 用基础模型嵌入作为**特征**，在带标签队列上训练轻量分类器（linear probe / 梯度提升）区分对比组。
- 用 **SHAP / attention 归因**做特征重要性 → 每个 BEST 类别得到一个候选短名单。
- 关键：**用 Geneformer in-silico 扰动过滤旁观者**——只保留扰动后能改变细胞状态方向的基因，大幅提升可重复性。

### Phase 3 — 机制验证
- 候选蛋白变异 → **Boltz** 看结构/结合改变。
- 免疫相关候选 → **EDEN** 看免疫原性。
- 靶点可成药性 → **ChEMBL** 看已知配体/机制。
- 这一步把"统计相关"升级为"机制可信"，是 predictive/CDx 能不能立住的分水岭。

### Phase 4 — 临床落地筛（三个硬门槛）
Claude 对每个候选打三个分：
1. **新颖性**：PubMed/bioRxiv 检索，是否已知？（决定是不是"新" biomarker）
2. **可及性 / 可测性**：能在血液/常规活检里测吗？有商品化抗体/qPCR/质谱方法吗？（决定能不能落地成检测）
3. **临床对齐**：和现有试验终点、入组标准对得上吗？CDx 能配到哪个在研药？（Clinical Trials + ChEMBL）
- 三门槛过不去的候选，科学上再漂亮也先放一边。

### Phase 5 — 组合面板与决策层
- **不追求单标志物，做多标志物 panel**：单个基因/蛋白很难同时高灵敏高特异，组合面板更稳。
- 每个 BEST 类别产出一个 panel + 一套判读规则（阈值/评分卡）。
- Claude 生成每个候选的**证据卡片**（modality、对比组、效应量、机制证据、可及性、先验文献、拟配对疗法）便于审计与后续实验验证。

---

## 4. 如何"组合起来"：一张矩阵 + 一条病人旅程

### 4.1 用途 × 模态 组合矩阵（核心）
行=临床用途，列=数据/模型层，格子=具体做法。**一个完整方案是把整行都填满、多模态互相印证。**

| 临床用途 (BEST) | 基因组层 (Evo2/CodonFM/ESM-VEP) | 单细胞层 (Geneformer) | 蛋白层 (PROTO/ESM-2) | 机制层 (Boltz/EDEN) | 临床/药物层 (Trials/ChEMBL) |
|---|---|---|---|---|---|
| **风险预测** | 胚系变异功能打分 ★ | 高危个体细胞状态 | 循环蛋白基线 | — | 流行病学佐证 |
| **诊断/分型** | 体细胞驱动突变 ★ | 细胞状态嵌入分型 ★ | 蛋白亚型 signature | 结构确认变异后果 | 现有诊断基线 |
| **入组富集/CDx** | 抗性/敏感突变 ★ | 应答细胞态 signature | 血浆 predictive 蛋白 ★ | Boltz 解释耐药机制 ★ | ChEMBL 配对药物 ★ |
| **疗效 (PD/Response)** | — | 扰动验证通路开关 | PD 血浆蛋白 ★ | 结合占据度 | 终点对齐 |
| **治疗程度 (Monitoring)** | ctDNA 变异负荷 | 残留病细胞态 | 纵向血浆蛋白 ★ | — | 监测终点 |
| **进程 (Prognostic)** | 克隆演化变异 | 网络中心性漂移 ★ | 蛋白轨迹 | — | 预后终点对齐 |

★ = 该用途下最优先、性价比最高的入口。

### 4.2 沿"病人旅程"串成一个系统
把不同 biomarker 按病人经过的时间轴串起来，就是一个完整的伴随分析系统：

```
筛查        诊断         入组/CDx        治疗中           进展
 │           │             │              │               │
Risk      Diagnostic    Predictive    Response/PD      Prognostic
胚系变异    细胞态分型     血浆predictive  血浆PD蛋白       网络漂移
面板       + 蛋白signature 蛋白 + 突变     纵向monitoring   + ctDNA克隆
(基因组)    (单细胞+蛋白)  (蛋白+机制)     (蛋白)          (单细胞+基因组)
```

- **同一批病人的多组学数据贯穿始终**，前一阶段的候选在后一阶段被复用/再验证（例如入组时的 predictive 蛋白，治疗中变成 monitoring 指标）。
- **多模态交叉验证**：一个候选如果在基因组（变异）、转录组（表达）、蛋白（丰度）三层都指向同一方向，可信度远高于单层信号——这是降低假阳性的核心手段。

---

## 5. 一个可落地的样例（免疫治疗，NSCLC 为例）

现状：PD-L1 IHC 和 TMB 是现有 CDx，但预测力有限。缺口清晰，适合挖"超越 PD-L1"的新 predictive biomarker。

1. **Phase 0**：Trials MCP 拉 NSCLC + checkpoint inhibitor 试验，确认入组多用 PD-L1≥1%/50%；对比组 = 用药后 ORR 应答者 vs 原发耐药者。
2. **Phase 1**：
   - Geneformer 嵌入应答者 vs 耐药者的肿瘤/免疫细胞 scRNA → 找 T 细胞耗竭/浸润相关基因态。
   - PROTO 分析基线血浆蛋白 → 找可及的循环 predictive 蛋白。
   - ESM-2 VEP 扫 HLA / 抗原呈递通路变异 → 逃逸相关变异。
3. **Phase 2**：嵌入做特征训练应答分类器，SHAP 归因；Geneformer 扰动验证候选基因确实驱动"冷→热"转化。
4. **Phase 3**：EDEN 评估候选新抗原免疫原性；Boltz 看 HLA-肽结合是否被变异破坏；ChEMBL 看候选靶点有没有在研配对药。
5. **Phase 4/5**：Claude 核查新颖性（是否已发表）、可及性（血浆可测优先）、和在研试验终点对齐 → 输出一个 **3–5 标志物的血浆 predictive panel + 一个 scRNA 诊断 signature**，各附证据卡片。

同一套数据顺带产出：进展期的 prognostic（Geneformer 网络漂移 + ctDNA 克隆演化）、治疗中的 monitoring（纵向血浆蛋白）。**一次数据采集，多个 BEST 类别复用。**

---

## 6. 验证与合规（别跳过）

- **分析验证 vs 临床验证**：模型挖出候选只是发现（discovery）。要落地必须做分析验证（检测方法准不准）+ 临床验证（独立队列、前瞻验证）。基础模型能**大幅缩小 wet-lab 验证的候选空间**，但替代不了验证。
- **避免数据泄漏**：训练/评估队列严格分离，尤其批次效应（不同中心、不同测序平台）——否则"高准确率"是假的。
- **可解释性**：CDx / 入组用途下监管要求机制可解释，所以 Phase 3 的 Boltz/扰动证据不是可选项。
- **新颖性尽调**：Phase 4 的 PubMed/bioRxiv/专利检索决定它到底是不是"新"，也影响 IP。

---

## 7. 先做什么（优先级建议）

1. **先定 1 个适应症 + 1–2 个 BEST 缺口**（建议从 predictive/CDx 或 monitoring 入手——临床价值和可及性最高）。
2. **先跑蛋白层 + 单细胞层**（PROTO + Geneformer）：这两层离"可测量、可落地"最近，出结果快。
3. 基因组层（Evo2/ESM-VEP）作为**风险/耐药的补充证据**，不必一开始全量铺开。
4. 用 Clinical Trials + PubMed **反向锚定**：先看清竞争基线，避免挖出"新但没用"的标志物。
5. 全流程用 Claude 做编排与证据三角，产出**可审计的候选证据卡片**，直接对接 wet-lab 验证。

---

### 一句话总结
**按 BEST 类别定义问题 → 三条组学腿（基因组/单细胞/蛋白）并行生成候选 → 用 Geneformer 扰动和 Boltz/EDEN 做机制过滤 → 用 Trials/ChEMBL/PubMed 做临床可及性与新颖性筛 → 组合成多标志物 panel 并沿病人旅程串成伴随分析系统。多模态交叉印证是降假阳性、挖出真正"新且有用"biomarker 的核心。**
