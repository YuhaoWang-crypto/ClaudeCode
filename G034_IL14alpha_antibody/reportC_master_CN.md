# G034 — TXLNA / α-taxilin 多模态药物发现主报告

**靶点:** TXLNA / α-taxilin（别名 IL-14α），UniProt **P40222**，546 aa
**用户构建体:** G034 sequence = P40222 残基 **210–546**（syntaxin 结合 coiled-coil 域）
**适应症背景:** 干燥综合征（Sjögren's / 干眼）、类风湿关节炎（RA）
**日期:** 2026-07-15
**性质:** 纯计算分析（无 wet-lab 实验），三模态（抗体 / 共价 PROTAC / siRNA）

> **重要前提说明（诚实基调）:** TXLNA 的经典身份是**胞内、无信号肽的运输蛋白**（结合 syntaxin，调控 SNARE 介导的胞吐）。数据库将其别名标注为 "IL-14α"（一个历史上描述的分泌型 B 细胞生长因子），但**这一分泌型身份从未在分子层面被确证**（HPA：胞质定位、secretome=null、血浆浓度未检出；且从未鉴定出 IL-14 受体）。因此整个"胞外可中和"治疗前提是一个**待验证的假说**，而非既定事实。本报告对每个模态都标注了相应的机制风险。

---

## 目录
1. 靶点生物学与结构
2. 抗原/表位表征（含结构图）
3. 蛋白互作网络与机制
4. 专利与 FTO
5. 模态一：抑制性抗体（含对接结构）
6. 模态二：共价 PROTAC（含口袋图）
7. 模态三：siRNA 敲低（含结合位点图）
8. 三模态对比与决策框架
9. 结构可视化总览
10. 诚实的局限与后续实验
11. 交付物清单

---

## 1. 靶点生物学与结构

**身份（已核实）:** 用户提供的 G034 DNA 序列翻译为 337 aa，NCBI BLAST 100% 匹配人 α-taxilin（TXLNA），精确对应 canonical P40222 的残基 **210–546**。基因身份：TXLNA（HGNC:30685，Entrez 200081，染色体 1p35.2）。

**结构域架构（AlphaFold P40222）:**
| 区域 | 残基 | 平均 pLDDT | 特征 |
|---|---|---|---|
| N 端 α1 | 1–209 | 中等 | 部分无序 |
| **中央 coiled-coil（= 本抗原核心）** | **210–460** | **94.8（极高）** | syntaxin 结合面 |
| C 端尾 | 461–546 | 49.9（无序） | 低置信度 |

全长 pLDDT 剖面：

![pLDDT profile]({{artifact:art_50d352ee-3369-4893-b2e5-4f802e69e002}})

**定位（Human Protein Atlas）:** 主定位 = 胞质（Cytosol），另见核质、中心体；**secretome location = null**；血浆浓度（免疫分析法 + 质谱法）**均未检出**；组织特异性 = 低（全组织表达）。蛋白类别 = "Predicted intracellular protein"。→ 这是支持"胞内运输蛋白"身份、削弱"分泌型细胞因子"身份的直接证据。

**疾病证据（诚实评估）:**
- **RA（强证据）:** TXLNA 蛋白在 RA 滑液、组织、血浆中经 ELISA 验证升高（n=106 队列，PMID 39634288；滑液成纤维样细胞与糖酵解酶相互作用，PMID 32377534）。
- **Sjögren's（弱/间接证据）:** 依赖 IL-14α 转基因小鼠发展出 Sjögren's 样表型（PMID 19038581, 17015757），而这一关联依赖有争议的 IL-14α=taxilin 身份。TXLNA 蛋白与 Sjögren's 之间**没有直接的原始文献分子链接**。

---

## 2. 抗原/表位表征

本抗原（残基 210–460 coiled-coil 核心）是一根延展的 α 螺旋，syntaxin 通过 coiled-coil 相互作用结合于此。基于结构，我们定义了两个空间上分离的表位面（相距 59.6 Å，确认为不同表面）：

![Target structure]({{artifact:art_8e492920-e4ed-4920-a85a-4c4eca3058a5}})

| 表位 | 残基（近似） | 治疗轴 | 阻断目标 |
|---|---|---|---|
| **Axis-2（胞内运输轴）** | ~295–320（syntaxin 沟槽，近 C245/C373） | 阻断 α-taxilin ↔ syntaxin-1A/3A/4A 结合 | 干扰 SNARE 胞吐调控 |
| **Axis-1（胞外免疫轴，假说）** | ~335–360（推定 IL-14R 面） | 阻断推定的 IL-14α ↔ IL-14R | 需先证实胞外池存在 |

**四个高置信度半胱氨酸（用户构建体内）:** C245、C373（pLDDT 97，coiled-coil 核心，最佳共价锚点）；C471、C523（无序尾区，pLDDT<63，较差锚点）。

---

## 3. 蛋白互作网络与机制

STRING（≥700）+ BioGRID + humanPPI 构建的 TXLNA 网络：8 节点，20 边，密度 0.71，TXLNA 为枢纽（度 7）。

![PPI network]({{artifact:art_d81a5d8a-9c59-4161-8f17-e0971a9e98b3}})

**通路富集（ORA，FDR≤0.05）:** SNARE binding（FDR 6.0e-10）、Taxilin family（1.2e-08）、Other interleukin signaling（2.4e-08，*仅通过 IL-14 别名注释产生，非独立证据*）、Exocytosis（1.0e-06）、Vesicle fusion to plasma membrane（1.2e-06）、Secretion by cell（8.1e-06）。

**机制解读（两轴）:**
- **运输轴（胞内）:** TXLNA 结合游离 syntaxin-1A/3/4，作为 SNARE 介导 Ca²⁺ 依赖胞吐的**负调控因子**——正是控制外分泌腺泪液/唾液分泌的机器。与干眼的联系机制上合理但间接。网络邻居 TSG101（ESCRT）提示：若有胞外 TXLNA，可能经非经典/外泌体途径分泌（与缺乏信号肽一致）。
- **免疫轴（胞外，治疗性）:** IL-14/白介素信号项仅经别名出现；网络邻居 TBK1（先天免疫激酶）。**无独立分子证据**表明 P40222 是作用于 IL-14 受体的分泌型细胞因子。

---

## 4. 专利与 FTO

**US7622574（唯一专用 TXLNA 专利，已过期）:** 优先权 2006-10-26，调整后到期 2028-02-28，**法律状态：过期**（维护费失效）。授权 claim 1 窄：分离的寡核糖核苷酸（SEQ ID NO:1/2/3）+ ≥1 个核糖 2'-O-(2,4-二硝基苯基/DNP) 修饰——非对全部 IL-14α RNA 抑制剂的属类权利要求。

> **注（专利标题分歧，待核实）:** 外部 G034 报告对同一专利号引用了不同标题（"Compositions and methods for detecting and treating autoimmune disorders"）以及抗体表位区（aa 493-523，mAb 1C6/1F2）。两个标题的分歧尚未逐字核对该专利授权文本；但**过期结论一致**，故 FTO 风险实际已消除。

**FTO 结论:** 采用标准 2'-OMe/2'-F 化学、不同靶位点的新型 siRNA，以及 de novo 抗体和共价 PROTAC，均**落在已过期的窄权利要求之外**。

**ChEMBL 化学物质:** CHEMBL6066423（alpha-taxilin）全库仅 1 个配体 CHEMBL5653589，Kd ≈3,539 nM（3.5 µM，kinobead 脱靶命中）。**无可成药小分子起点。**

---

## 5. 模态一：抑制性抗体

**方法:** RFdiffusion-Ab → ProteinMPNN → RF2 de novo 设计流程（Modal A100-80GB GPU），针对两个结构表位各设计抗体，RF2 过滤标准 interaction_pae<10 且 pred_lddt>0.8，然后 Boltz-2/Chai-1 co-fold 验证 + 完整 IgG1 组装。

**当前结果（Axis-2 syntaxin 表位，部分数据）:** 流程完成 40 个 RFdiffusion 骨架 + 160 个 ProteinMPNN 序列；RF2 在作业超时前评分了 17/160 个序列。

| 指标 | 值 |
|---|---|
| 最佳 interaction_pae | 21.6（**>10 = 预测结合较弱**）|
| 中位 interaction_pae | 24.1 |
| 最佳 pred_lddt | 0.89（结构折叠置信度高）|
| 通过过滤（i_pae<10）的设计数 | **0 / 17** |
| VH / VL 长度 | 118 / 110 aa |

![Antibody docking]({{artifact:art_894bebf4-f2b6-409f-affd-d6b2d752f3f2}})

> **诚实评估:** 当前评分的 17 个设计中**没有一个通过 interaction_pae<10 的结合置信度阈值**（最佳 21.6）。结构上抗体折叠正确（Fv β-三明治清晰），但预测的结合界面较弱——这在针对延展 coiled-coil（缺少凹陷口袋供 CDR 环抓握）的 de novo 设计中并不意外。**完整重跑正在进行**（axis-2 对全部 160 序列重新 RF2 评分 + axis-1 全新 24 设计），以确认是否存在更好的结合体。这是一个**真实的负结果信号**，不应粉饰为成功。

---

## 6. 模态二：共价 PROTAC

**策略:** 共价弹头锚定 C245/C373（pLDDT 97 的高置信度半胱氨酸）+ E3 连接酶募集（CRBN/VHL），将结合转化为催化性降解——绕开胞内靶点的抗体区室问题。

![PROTAC pocket]({{artifact:art_d1102cb2-a2ff-4815-b3d3-5b5154d6ef86}})

**6 个化学有效设计（3 弹头 × 2 E3）:**
| ID | 弹头 | E3 | MW (Da) | cLogP | TPSA (Å²) |
|---|---|---|---|---|---|
| PROTAC-01 | acrylamide | CRBN | 530.5 | -0.23 | 169.4 |
| PROTAC-02 | acrylamide | VHL | 588.7 | 2.26 | 135.3 |
| PROTAC-03 | vinylsulfonamide | CRBN | 566.6 | -0.47 | 186.5 |
| PROTAC-04 | vinylsulfonamide | VHL | 624.8 | 2.02 | 152.4 |
| PROTAC-05 | chloroacetamide | CRBN | 553.0 | -0.18 | 169.4 |
| PROTAC-06 | chloroacetamide | VHL | 611.2 | 2.31 | 135.3 |

**关键发现（诚实）:** TXLNA 无深口袋（如图 A/B，黄色 SG 硫醇溶剂暴露于表面而非埋于口袋）。策略依赖**共价捕获反应性表面硫醇**，而非口袋结合。C245 与 C373 在延展 coiled-coil 上相距 ~128 残基，是**两个独立的锚点位点，而非共享口袋**。首选：acrylamide（临床验证弹头，如 afatinib/ibrutinib）+ CRBN。所有设计均为**模块化示意结构 + 计算理化性质，非对接姿态**——需 wet-lab 共价占据（intact-MS）、三元复合物（TR-FRET）、细胞降解（HiBiT）验证。

---

## 7. 模态三：siRNA 敲低

**靶 mRNA（已核实）:** RefSeq **NM_175852**（transcript variant 1，canonical），4,865 nt，CDS 135–1775（1,641 nt = 546 aa）——直接从 NCBI 获取并翻译确认。

**设计:** CDS 全长 19-mer 扫描（Reynolds 2004 评分），GC 36–58% 过滤 + 同聚物/免疫刺激基序排除 → 696 个通过候选 → top-5 → top-3。

![siRNA binding sites]({{artifact:art_da685983-f5b3-48cc-9202-3b3443d9980e}})

**Top-3（最低脱靶风险）:**
| 排名 | 靶位（CDS+）| mRNA 位 | 反义链（guide）5'→3' | GC% | Reynolds | 强脱靶 | Evo2 |
|---|---|---|---|---|---|---|---|
| #1 | 957 | 1091 | UUUGUCGAUAUGCUCCUCGUU | 47.4 | 8/8 | 0 | -0.49 |
| #2 | 512 | 646 | UGCAUCAGCAACGUGAUCUUU | 47.4 | 8/8 | 0 | -0.84 |
| #3 | 848 | 982 | UCCAUCUGCAGCUGAAUGUUU | 47.4 | 8/8 | 0 | -0.49 |

**脱靶筛选:** BLASTn（human RefSeq RNA）——全部 5 候选 100% 匹配所有 TXLNA 转录本，**0 个强脱靶**（最差脱靶 ≤30 bits，E≥2.6）。**Evo2-7b 序列上下文评分**（在 A100 上执行）：全部 5 位点落在正常编码区（mean log-lik −0.43 至 −0.84），无低复杂度/重复异常。

**shRNA 构建体:** U6 启动子 + sense(19) + loop(TTCAAGAGA) + antisense(19) + TTTTTT 终止子；提供 LNA 增强的 guide 末端选项（核酸酶抗性，保留种子区未修饰以维持 RISC 装载）。

**推荐化学:** 2'-OMe/2'-F 交替 + 末端硫代磷酸酯（PS）——与 US7622574 的 2'-O-DNP claim 不同，FTO 清晰。**递送警示:** LNP/GalNAc 偏好肝脏；泪腺/唾液腺递送需组织特异工程或局部给药。

---

## 8. 三模态对比与决策框架

| 维度 | 抗体（Axis 1/2）| 共价 PROTAC | siRNA |
|---|---|---|---|
| 靶区室 | 胞外（需证实池）/ 胞内需内化 | 胞内 ✓ | 胞内 mRNA ✓ |
| 设计成熟度 | de novo 序列，当前 i_pae 弱，重跑中 | 6 化学有效设计，无对接 | top-3 完整，0 强脱靶，Evo2 验证 |
| 关键风险 | 胞外池是否存在；对 coiled-coil 结合弱 | 无口袋，依赖共价；细胞渗透 | 递送到外分泌腺 |
| FTO | 清晰（de novo）| 清晰 | 清晰（新位点+标准化学）|
| 到 POC 估计 | 18–24 mo（若靶点成立）| 36–48 mo | 9–18 mo（最快）|

**决策关口（gating decision）:** **pSS 血清 TXLNA ELISA** 是决定一切下游的 go/no-go 实验：
- 若存在胞外可中和池 → 抗体（Axis-1）成为 lead；LYTAC/AbTAC 可复用抗体。
- 若无胞外池 → 胞内模态（siRNA 首选 / 共价 PROTAC）为唯一可行路径，"清除蛋白"优于"阻断界面"。

**当前推荐优先级:** siRNA（P1，最快、设计最完整、脱靶清洁）> 抗体（P2，待关口实验 + 重跑确认结合体）> 共价 PROTAC（P3，概念验证，需片段筛选）。

---

## 9. 结构可视化总览

本报告已获得并可视化以下结构信息：

| 结构 | 状态 | 图 |
|---|---|---|
| 靶点结构（TXLNA coiled-coil）| ✅ AlphaFold P40222，构建体 210–460 | §2 图 |
| 靶点+受体复合物 | ❌ 无实验/预测的 TXLNA-syntaxin 复合物结构 | — |
| 靶点抗体结构（de novo VH/VL）| ✅ 部分（17 设计，重跑中）| §5 图 |
| 靶点+抗体对接结构 | ✅ RF2 共建模复合物（i_pae 弱）| §5 图 |
| PROTAC + 结合口袋 | ◐ 共价锚点 C245/C373（无深口袋，无对接姿态）| §6 图 |
| siRNA 结合位点 | ✅ NM_175852 上精确 mRNA 位置 | §7 图 |

---

## 10. 诚实的局限与后续实验

1. **IL-14α 分泌身份未确证** —— 整个胞外治疗轴依赖此假说；血清 ELISA 是决定性实验。
2. **抗体当前无通过阈值的结合体** —— 17/160 部分数据 i_pae 全部 >10；重跑进行中，结果可能仍为弱结合（coiled-coil 是困难靶面）。
3. **无 TXLNA-syntaxin 复合物结构** —— 表位定义基于单体 AlphaFold + 文献，非共复合结构。
4. **PROTAC 为示意设计** —— 无三元复合物建模、无对接姿态；MW 530–625 Da 偏小，真实连接子优化会升至 700–1000 Da 并降低渗透性。
5. **siRNA 效力为预测** —— 脱靶为 BLAST + Evo2 上下文，未穷尽种子介导的 miRNA 样脱靶；敲低效率未经实验测定。
6. **递送未解决** —— 所有胞内模态面临外分泌腺递送难题。

**建议实验优先级:** (1) pSS 血清 TXLNA ELISA（n=30+30）定成败；(2) 确认 top-3 siRNA 敲低效率（RT-qPCR/Western）；(3) 若抗体重跑得到 i_pae<10 设计，做 Boltz-2/Chai-1 co-fold + SPR。

---

## 11. 交付物清单

**结构图:** fig_target_structure.png、fig_antibody_docking.png、fig_protac_pocket.png、sirna_binding_map.png
**siRNA:** TXLNA_siRNA_Design_Report_v2.txt、top3_and_shRNA.json、off_target_analysis_v2.json、evo2_scores.json、txlna_mrna_NM_175852.fasta
**PROTAC:** protac_designs_valid.csv、PROTAC_Feasibility_v2.txt、protac_structures_grid.png
**抗体:** axis2_best_complex.pdb、axis2_partial_designs.json（重跑完成后更新）
**对比:** G034_Parallel_Comparison_CN.md（外部报告 vs 本分析）

---

*本报告为纯计算分析，所有 wet-lab 验证步骤均已明确标注为"未执行"。定量结合预测为设计估算，非实验测量。诚实的机制/递送警示贯穿全文。*
