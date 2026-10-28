# AlphaProtein Novo 复现性审计报告

**审计对象**
- 论文 A：Wu et al. 2026, *Designing enzymes for new-to-nature chemistry and non-natural substrates with AlphaProtein Novo*, bioRxiv [10.64898/2026.10.01.756017](https://www.biorxiv.org/content/10.64898/2026.10.01.756017v1)
- 论文 B：Li et al. 2026, *Design-Evolution Synergy Enables New-to-Nature Enzyme Chemistry*, bioRxiv [10.64898/2026.10.01.754751](https://www.biorxiv.org/content/10.64898/2026.10.01.754751v1)（DeepMind + Arnold 实验室）
- 代码：[google-deepmind/alphaprotein-novo](https://github.com/google-deepmind/alphaprotein-novo) @ `e7ce276`（2026-10-06，仅 1 个 commit 的浅克隆）

**审计日期**：2026-10-28
**方法**：静态分析。**没有下载权重，没有运行流水线**（本机无 GPU）。论文正文与 Methods 已读；**Supplementary Information 未能获取**（>10MB 超出抓取上限，直接下载被 Cloudflare 429 拦截）——这一点在下文相关处标注。

---

## 一句话结论

**论文的核心结论无法复现，而且卡点不止是算力。**

三层封堵，由硬到软：

1. **被实验验证的酶，序列没公开** —— 复现对象根本不存在于公开材料中
2. **湿实验数据一条没发** —— 论文最核心的定量主张（指标能富集出有活性的设计）无从检验
3. **论文自称的两大成功要素之一，开源代码里没有实现**

---

## 1. 数据可得性：复现对象缺失

论文 A 的 Data Availability 自己写明：设计序列与预测结构 "will be made public upon publication"（尚未发表，故尚未公开）。实测核查：

| 项目 | 状态 |
|---|---|
| FASTA / 序列表 / CSV | **仓库内 0 个**（`find` 全盘确认） |
| 实验测试过的设计结构 | **0 个**。`GDM_KE_1483`、`GDM_SE_2937`、`GDM_CT_0103`、`GDM_NT_0151`、`GDM_DEHP_0376` 只在 manifest 的 `description` 字段里被文字提及 |
| 唯一附带的设计结构 | `examples/dehp_esterase_partial_diffusion/GDM_DEHP_0176.cif`（1301 原子）—— 是 90°C 明星设计 **0376 的母本**，不是 0376 本身 |
| 动力学/活性数据（k_obs 等） | **0 条** |
| 模型权重 | 不在仓库内，但 GCS 链接实测可下载（见 §5） |

### 1.1 那个唯一的设计文件，其序列不是蛋白

从 `GDM_DEHP_0176.cif` 提取链 A（191 aa）的组成：

| 残基 | 占比 |
|---|---|
| Ala | **41.9%** |
| Leu | 19.4% |
| Ala+Leu 合计 | **61.3%** |
| Cys / Phe / Lys / Met / Trp / Tyr | **全部为 0** |

一个不含任何赖氨酸、不含任何芳香族残基的 191 aa「酯酶」不可能是被表达纯化过的蛋白。这是扩散模型共生成的原始序列，按设计会被下游 LigandMPNN 整体替换。

**自洽性佐证**：仓库自带的筛选阈值要求 `sequence/max_single_aatype_frac ≤ 0.25`，而这条序列 41.9% —— 会被它自己的过滤器淘汰。

> **因此**：论文 A 的湿实验结论（GDM_DEHP_0378 达到天然酶 QHH 活性的 6%@4h / 13%@24h、天然 Q7SIG1 的 26%/34%；GDM_DEHP_0376 在 90°C 比室温活性高 14 倍、在 75% 乙腈中高 2 倍；新骨架命中率 11%、回收骨架 39%）**无法独立核验**。这些需要基因合成、表达、纯化、LC-MS 测活，是数月的实验室工作，没有任何代码能替代。

---

## 2. 论文自称的两大成功要素，一个缺失

论文 A 摘要把成功归因于两点。逐一核查：

### ✅ (a)「基于 AF3 预测的机制启发指标」—— 基本已实现

定义层面可核验，代码是真的：

| 论文指标 | 代码实现 | 判定 |
|---|---|---|
| 多随机种子 AF3 集合打分 | `evaluate_design.py:506` `_aggregate_seed_metrics`，输出 `/mean` `/median` `/min` `/max` `/std` | 机制吻合 |
| 多催化态打分 | 各 manifest 的 `folding.states`；`multi_state_metrics.compute_multi_state_consistency` | 吻合 |
| Kemp 催化碱 → 6-NBT 的 N3 距离 | `kemp_metrics.py:68,76` `KEMP_ELIMINASE_LIGAND_METADATA` 将 `6NT` 映射到 `acid_atom='N3'` | 吻合 |
| 丝氨酸三联体极性相互作用 | `serine_esterase_metrics.py:285`，几何判据来自 `hbonds.py:142` `LAUKO_HBOND_THRESHOLDS`（2.0–3.6 Å + 角度/二面角，引 Lauko et al. 2025） | 吻合 |
| ES vs TI1 配体构象一致性 | `serine_esterase_metrics.py:955-956` | 吻合 |
| 催化碱疏水埋藏 | `kemp_metrics.py:438` `cat_base_tip_ca_neighbors_mean` | **部分**——是 CA 近邻计数的堆积代理，不区分疏水性 |
| 卡宾/环丙烷化立体选择性 | `carbene_metrics.py` 只输出 `alanine_ratio`、`heme_rmsd`、`transition_state_rmsd` | **缺失**——无立体选择性指标 |

**但阈值无法核对**：论文几乎不公布任何数值 cutoff（它公布的是富集统计量：EF@R50、Mann-Whitney U=10503、Spearman ρ=−0.19）。仓库里的 cutoff 只存在于示例 notebook `examples/filter_designs.ipynb`，而该 notebook 自称是 "benchmark filtering thresholds"，**并未声称等同于论文所用**。

另外：`run_pipeline.py` 只有 4 个阶段（generation / resequence / folding / evaluation，`:68-71`），**流水线里根本没有「筛选」这一阶段**。阈值不在 `src/` 任何地方。

### ❌ (b)「基于序列集合的骨架打分」—— 开源代码里没有

论文把它定义为多种子打分的推广：用一组 LigandMPNN 衍生序列来评价一个骨架，**要求每条重设计序列都通过同一套过滤指标**。具体规则论文有写：Kemp 要求「7 条里至少 6 条通过」，丝氨酸酯酶要求「10 条全部通过」，MPD 要求「10 条里 9 条通过」。

| 组成部分 | 状态 |
|---|---|
| 每个骨架生成 N 条重设计序列 | **有**（`resequence.num_sequences`，每条成为独立设计 `<prefix>_seq<NN>`） |
| 对每条序列折叠 + 打分 | **有**（阶段 3-4 把每条当独立设计处理） |
| **按母本骨架聚合通过/不通过（6-of-7 / 10-of-10 规则）** | **没有** |

我独立核实过：`evaluate_design.py` 从不按母本前缀分组，没有任何 `groupby`，没有剥离 `_seqNN` 后缀的逻辑；`evaluation_summary.csv` 是「一行一个设计」（即一行一条序列）；notebook 的 `get_filter_mask` 严格按行计算，`load_eval_results` 保留 `sample_name`（如 `kemp_0000_seq00`）却从不推导骨架 ID。

更关键的是默认值：`design_manifest.py:356` 写的是 `num_sequences: pydantic.PositiveInt = 1` —— **开箱即用根本没有序列集合**；示例 manifest 也只提到 2。

> 这个缺口比 README 那句轻描淡写的 "reduced settings" 严重得多：被论文列为关键成功要素的机制，在开源版本里被降格成了「生成几条序列，各自独立打分」。

---

## 3. 哌啶那个指标：代码与命名反向

这是你最关心的「合成类药物哌啶」那条线的核心计算指标。`nitrene_metrics.py:233-285`：

```python
piperidine_propensity = d(N1–C5) − d(N1–C4)
```

两个原子**都在底物配体 `pd00` 内部**。不涉及任何蛋白原子，不涉及 Fe，不涉及角度 —— 纯粹是 AF3 预测出的底物自身构象（蛋白的影响只是间接地通过预测姿态体现）。

**化学上对底物是正确的**：示例 manifest 折叠的底物 SMILES 为 `NCCCCCC1=CC=CC=C1`（5-苯基戊胺，氮烯前体），N1–C1…C5(–Ph)。C4 成环 → 五元吡咯烷；C5（苄位）成环 → 六元哌啶。与 manifest 里的产物 `piperidine_R` = 2-苯基哌啶、`pyrrolidine_S` = 2-苄基吡咯烷一致。✅

**但有三个问题：**

1. **命名与数值方向相反。** 哌啶型近攻击构象下 N1–C5 是成键距离（短），N1–C4 是 1,3 环关系（约 2.5 Å），故 `diff ≈ 0` 或为负；吡咯烷型构象下 C5 悬空，`diff` 显著为正。而筛选阈值是 `piperidine_propensity/min ≤ 0.35` —— **选的是数值低的**。过滤方向化学上说得通，但标识符名字是误导性的，且仓库里没有任何地方说明符号约定。

2. **论文从未描述这是「两个距离之差」。** 论文正文能查到的唯一描述是 "near-attack positions of the bond-forming carbon and nitrogen atoms in piperidine synthases"，细节推给了 SI（Fig S22-23, S29 / Supplementary Notes）—— 而 SI 我拿不到。这个描述同样可以指单个 N···C 距离。0.35 Å 这个 cutoff 在论文里完全找不到对应。**此项判定为「无法确认是否吻合」**，需要 SI 才能定论。

3. **静默失效 + 错误状态下输出无意义值。** 原子按名字在 `pd00` 中查找，缺失时三个键全返回 `np.nan`（`:258-267`），所以一旦 AF3 的 SMILES→原子命名有变，指标会无声失效。另外 `evaluate_design.py:371` 把氮烯指标套用到 `('complex', 'heme_substrate', 'substrate')` 三个状态；在 `complex` 状态下 `pd00` 是合并后的 Fe-血红素-哌啶配体，此时 `N1/C4/C5` 是任意的卟啉原子 —— `complex/piperidine_propensity` 会被输出但毫无意义。（notebook 本身正确地只用 `heme_substrate`。）

---

## 4. 示例设置 vs 论文实际用量

README `:330-331` 自己承认是 Google 内部流水线的**移植版**，且 "some settings have been reduced relative to the paper"。量化这个差距：

| 设置 | 论文 | 示例 manifest | 差距 |
|---|---|---|---|
| 每骨架重设计序列数 | 7（Kemp）/ 10（丝氨酸酯酶）/ 10（MPD） | `"num_sequences": 2` | **3.5–5×**；schema 默认是 1 |
| k-of-n 通过规则 | ≥6/7、10/10、9/10 | 未实现 | **缺失** |
| AF3 种子数 | Methods **未给数字**；Fig 3B 图注 7 个，Fig 3D 10 个 | `"seeds": [0,1,2,3,4]` → 5 | 约 1.4–2×；但**论文从未给出流水线实际设置**，无法严格对比 |
| 生成设计数 | **从未给出**（论文只报「测试了」多少：Kemp 1536 条 /7 轮、丝氨酸酯酶 1569 条 /8 轮、氮烯 86+90 /2 轮、卡宾 91 /1 轮、DEHP 酶 81 条首轮 + 两轮） | `"num_designs": 2` | **比率无法计算**（论文给的是测试数而非生成数） |
| 迭代轮数 | 7–8 轮（Kemp、丝氨酸）/ 2 轮（氮烯）/ 3 轮（DEHP） | 每个 manifest 1 轮 | 多轮campaign不可复现 |
| 部分扩散 | 对母本序列同一性 ~88%、RMSD ~0.5 Å；10 条 MPD 样本 | 575/1000 步，2 条样本 | 参数化方式不同（按步数而非同一性/RMSD） |
| `use_side_chain_context` | — | 全部 `false`，README `:404-409` 称这是「论文所用的确切设置」 | 标为忠于论文（schema 默认为 `true`） |

### 一处论文内部不一致（未能解决）

Methods 说 Kemp 的 MR 判据是「7 条中至少 6 条通过」，但 Figure 4 图注描述的选择标准是「7 条重设计中至少 **1** 条通过，通常 7/7」。这两个判据不同，论文没有调和。（经抓取摘要转述，我无法拿到完整图注原文。）

---

## 5. 权重与许可

**实测三个权重链接全部可达，且完全无需申请或认证**：

| 文件 | 大小 | 状态 |
|---|---|---|
| `generator.bin.zst`（AP Novo 生成器） | 547 MB | HTTP 200 |
| `af3_leaving_atom.bin.zst`（AF3-LA） | 1.02 GB | HTTP 200 |
| `af3.bin.zst`（标准 AF3） | 1.02 GB | HTTP 200 |

AF3 权重通常需要走申请表，此处直接匿名可下 —— 实际约束是「名誉制」的使用条款，而非技术门禁。

**法律硬约束：权重仅限非商业用途。** `WEIGHTS_PROHIBITED_USE_POLICY.md` 明确：只有非商业机构（大学等）可用；研究者不得与商业机构共享；禁止用输出去训练其他模型；禁止规避生物安全筛查；分发输出必须显著标注来源。**如果你或你的机构是商业实体，这一条就直接封死了合法复现路径。**

一个容易误读的点：`num_blocklisted_motifs` 过滤器**不是生物安全筛查**，而是过滤 Strep-tag、His-tag、FLAG、c-Myc、HA、V5、TEV/3C/凝血酶切位点这类纯化标签（`sequence_metrics.py:25`），目的是避免设计序列意外包含克隆假象干扰纯化。

---

## 6. 代码质量：是真货

和「无法复现结论」要分开评价 —— 代码本身是认真的工程产物：

- 49 个 Python 文件、18141 行，**全部语法编译通过**，仅 1 处 `NotImplementedError`（`distance_difference.py:241`）、1 处 `TODO`
- 5 套完整酶设计流水线：Kemp 消除酶、4MU-Ac 丝氨酸酯酶、DEHP 酶（新骨架 + 部分扩散两种）、卡宾转移酶、氮烯转移酶/哌啶合成酶
- 依赖声明完整（`jax==0.10.2`、AF3 从 git 装、`rdkit==2026.3.3`），多 GPU / 多 CPU 分片齐备
- 筛选阈值**公开**在 notebook 里（Kemp 18 条 / 丝氨酸 41 条 / DEHP 41-44 条 / 氮烯 23 条），这是实在的可复现性加分

### 已发现的疑似缺陷

1. **氮烯 notebook 很可能跑不起来。** `filter_designs.ipynb` 里有**两份** `NITRENE_TRANSFERASE_THRESHOLDS`：cell 3（**markdown，不执行**）四个配体门控键用 `complex/` 前缀，cell 4（**code，实际执行**）改成了 `heme_substrate/`。我已逐键核实这个差异真实存在。按静态分析，`heme_substrate` 那版会出问题：设计结构的配体是链 `C`/`pd00`，而 `heme_substrate` 状态折叠的是链 `X`=`HEM` + 链 `Y`=`pd00`；`ligand_metrics.has_matching_ligands` 比较 `(chain_id, res_id, res_name, atom_name)` 键集会返回 `False`，于是 `structure_metrics.py:641-643` 的两个 pocket-RMSD 键根本不会被输出 → notebook 取 `df[metric]` 会 `KeyError`；`interface_lddt` 在同条件下被显式置为 `np.nan`（`:581`），那条过滤永远无法通过。两份只差状态前缀，疑为编辑疏漏。*（此项为静态推断，未实际执行验证。）*
2. **卡宾 campaign 没有筛选配方**：`carbene_transfer` 是受支持的 evaluation suite 且有 manifest，但 notebook 只覆盖 Kemp/丝氨酸/DEHP/氮烯 —— 91 条设计选出 48 条做手性分析的标准没有公开，且 `carbene_metrics.py` 不输出任何选择性指标。
3. **16 条 H 键过滤器用的是最宽松判据**：丝氨酸/DEHP 的 H 键阈值全是对 `/mean` 取 `≥ 0.01`。底层是 0/1 指示量在 5 个种子上求均值，所以 `≥0.01` 等于「5 个种子里至少 1 个形成」—— 是对种子集合最宽松的用法，不是共识要求。
4. **README 示例与实际 manifest 不符**（无害）：README `:174` 写 `kemp_eliminase_motif.cif` 和 `"seeds": [230]`，实际文件是 `set1_tsa_top20_E_S_X_tsa-310.cif`、`seeds: [0..4]`。

---

## 7. 无法复现的论文内容清单

1. **Figure 3（全文核心论证）不可复现。** Fig 3B/3D 是「指标 vs 实测 k_obs」的回溯/前瞻富集分析。仓库**零实验数据** —— 没有 k_obs、没有设计-活性对照表、没有任何 csv，也没有那 1536 条 Kemp / 1569 条丝氨酸设计。所以「回溯 2.7 倍 / 前瞻 3 倍命中率富集」和「最佳 k_obs 从 0.49 s⁻¹ 提升到 1.38 s⁻¹」都无从检验。
2. **MR / MPD 分析原理上就不可复现** —— 聚合规则缺失（§2b），且每轮设计集合均未提供。
3. **「这个 manifest 产出了设计 X」无法验证** —— X 的结构和序列都没发。
4. **哌啶选择性指标是否与论文吻合，需 SI 才能定论**（SI 不可获取）。

---

## 8. 如果要真正推进验证，需要什么

| 目标 | 所需 |
|---|---|
| 复现在硅筛选流水线 | 多张 GPU（A100/H100 级），数天到数周；需确认非商业身份 |
| 验证论文的富集主张 | 等作者发表后公布序列 + 活性数据，或向作者索取 |
| 独立核验某个具体设计确实可折叠/几何合理 | 等序列公布；届时可用独立工具（如 Boltz）交叉预测 |
| 核实哌啶指标定义 | 获取 SI（Fig S22-23, S29） |
| **验证酶真的有催化活性 / 90°C 耐热** | **只能做湿实验**：基因合成 → 表达纯化 → LC-MS 测活。没有任何计算替代品 |

---

## 附：论文 B 的定位

论文 B（*Design-Evolution Synergy*，Arnold 实验室主导）是姊妹工作：以 motif scaffolding + 计算再设计 + **定向进化**构建从头血红素蛋白，做区域/对映选择性 C–H 胺化合成手性哌啶。它的关键词是 design-**evolution** synergy —— 即从头设计给出起点，真实活性靠定向进化获得。这与论文 A 的 DEHP 酶数据相互印证：首轮 81 条设计只有 4 条高于背景，经两轮迭代新骨架命中率才到 11%。

> **解读提醒**：「AlphaProtein Novo 生成全新酶」是准确的；但「生成即可用」不是这两篇论文的主张。公开数据显示的是：从头设计提供可进化的起点，绝对活性仍显著低于天然酶（DEHP 酶最佳设计为天然 QHH 的 6%–34%），而 90°C 那条是**相对**主张 —— 设计 0376 在 90°C 比自身室温活性高 14 倍、且此条件下两个天然酯酶均无活性，不等于它在 90°C 下活性很高。

---

*审计者：Claude Code。方法为静态分析，未下载权重、未运行流水线、未获取 SI。标注「静态推断」处未经执行验证。*
