# IL-23 / JAK 应答预测缺口深挖 + IL23R/TYK2 调控变异清单（喂给 PROTO）

> 承接 `ibd-candidate-panel-A-predictive.md`。本文件专攻两个**新颖性最高的空白**：IL-23 类（ustekinumab/risankizumab）与 JAK 类（tofacitinib/upadacitinib）的应答预测，并把 **IL23R / TYK2 的变异清单**整理成 PROTO（Evo2 / AlphaGenome / Enformer·Borzoi / SpliceAI）可直接消费的输入。
>
> **署名**：文献证据来自 PubMed（According to PubMed），每条附 DOI 链接。变异 rsID 来自公开 GWAS 文献知识整理，**运行前须对 GWAS Catalog / Open Targets 再校验**（见 §3 注）。

---

## 1. 缺口证实（检索即证据）

系统检索（PubMed，2026-07）显示：**IL-23 与 JAK 类在 IBD 中几乎没有已验证的"应答预测" biomarker**——这不是检索不足，而是真实空白，恰是最大机会。

| 机制类 | 检索命中情况 | 结论 |
|---|---|---|
| anti-TNF | 大量直接证据（OSM、TREM1、GIMATS 等） | 已拥挤，做差异化难 |
| α4β7 (vedo) | 中等（α4β7/MAdCAM 机制自洽） | 半拥挤 |
| **IL-23 (uste/risa)** | **应答预测证据薄弱**；命中多为银屑病或"定位/换药"（非预测） | **空白 → 高机会** |
| **JAK (tofa/upa)** | **应答预测近乎为零**；命中为安全性/其他适应症 | **空白 → 高机会 + 安全性 biomarker 需求** |

**遗传学线索（支撑机制可预测性）**（据 PubMed）：
- **IL23R R381Q（rs11209026）是功能缺失等位、保护 CD/UC**：Pidasheva 等在原代 T 细胞中证明 Q381 携带者 **IL-23 刺激后 STAT3 磷酸化降低、IL-23 响应 T 细胞减少**（*PLoS One* 2011，[DOI](https://doi.org/10.1371/journal.pone.0025038)）；综述见 Abdollahi 等（*J Immunotoxicol* 2016，[DOI](https://doi.org/10.3109/1547691X.2015.1115448)）。→ **IL-23 通路活性个体差异有遗传基础，可作为 IL-23 类应答预测的生物学锚。**
- **TYK2 蛋白编码变异（P1104A/rs34536443、A928V/rs35018800、I684S/rs12720356）降低激酶活性、保护自身免疫**，对 IBD 有提示性关联：Diogo 等（*PLoS One* 2015，[DOI](https://doi.org/10.1371/journal.pone.0122271)）。→ **TYK2 是 JAK 类的天然"部分抑制"人体模型，其变异可指示 JAK/TYK2 通路依赖度 → 应答预测假设。**

---

## 2. 缺口挖掘策略：从零挖出 IL-23 / JAK "响应细胞态"

因为没有现成标志物，主力是 **Geneformer 从数据里无监督地挖**，PROTO/ESM-2 提供遗传-调控层的正交证据。

### IL-23 类（ustekinumab / risankizumab）
- **Geneformer**：在 ustekinumab/risankizumab 应答分层的黏膜 scRNA 上，
  1. 无监督聚类找 **IL-23 响应细胞态**（Th17、ILC3、IL23R+ 单核吞噬细胞）；
  2. **in-silico 扰动** IL23R / IL12B 看能否把"炎症态"推向"缓解态"，验证驱动性；
  3. 用网络中心性找 IL-23 轴在应答者 vs 无应答者中的地位差异。
- **PROTO（独特贡献）**：IL23R 遗传学最丰富，把其调控/剪接变异（§3）用 Evo2+AlphaGenome 打表达影响、SpliceAI 打剪接影响——**判断个体 IL-23 通路的"基因型激活水平"**，作为应答预测特征。
- **ESM-2 VEP**：R381Q 等编码变异有害性零样本打分。
- **候选假设**：治疗前"IL-23 轴激活 + 低保护性基因型剂量" → ustekinumab/risankizumab 更可能应答。

### JAK 类（tofacitinib / upadacitinib）
- **Geneformer**：找 **JAK/STAT + I 型/II 型干扰素响应细胞态**；扰动 JAK1/TYK2 验证。
- **PROTO**：TYK2 变异（§3）调控+剪接打分，判断 TYK2/JAK 通路依赖度。
- **候选假设**：治疗前黏膜"高 JAK/STAT·IFN 签名" + 非保护性 TYK2 基因型 → JAK 抑制剂更可能应答。同时布局 **安全性 biomarker**（JAK 类有 VTE/肿瘤风险，需风险分层）。

---

## 3. IL23R / TYK2 变异清单 → PROTO 输入

> 下表即 PROTO 的"待打分变异列表"。**分工**：编码错义 → ESM-2 VEP（+可选 Boltz 结构）；非编码/内含子/UTR → PROTO 的 Evo2+AlphaGenome（表达/染色质）与 SpliceAI/Pangolin（剪接）。机器可读版见同目录 `proto_variants_il23r_tyk2.tsv`。

### IL23R（1p31.3）
| rsID | 变异 | 类型 | 已知效应 | PROTO/引擎分析 | 证据 |
|---|---|---|---|---|---|
| rs11209026 | p.Arg381Gln (R381Q) | 编码 missense | **功能缺失、保护 CD/UC** | ESM-2 VEP（有害性）+ Boltz（结构）+ SpliceAI（是否近剪接） | Pidasheva 2011 [DOI](https://doi.org/10.1371/journal.pone.0025038) |
| rs7517847 | 内含子 | 非编码 | IBD 关联（调控假说） | Evo2+AlphaGenome 表达/eQTL、Enformer/Borzoi 染色质 | GWAS（待校验*） |
| rs11465804 | 内含子（近剪接） | 非编码 | IBD 关联 | **SpliceAI/Pangolin 剪接** + AlphaGenome | GWAS（待校验*） |
| rs1004819 | 内含子 | 非编码 | IBD 关联 | Evo2+AlphaGenome 表达 | GWAS（待校验*） |
| rs10889677 | 3′UTR | 非编码 | 影响 mRNA 稳定/miRNA 结合（表达假说） | AlphaGenome 表达 + RNA 结构(ViennaRNA) | GWAS（待校验*） |

### TYK2（19p13.2）
| rsID | 变异 | 类型 | 已知效应 | PROTO/引擎分析 | 证据 |
|---|---|---|---|---|---|
| rs34536443 | p.Pro1104Ala (P1104A) | 编码 missense | **降低激酶活性、保护自身免疫；对 IBD 有提示** | ESM-2 VEP + Boltz（激酶域结构） | Diogo 2015 [DOI](https://doi.org/10.1371/journal.pone.0122271) |
| rs35018800 | p.Ala928Val (A928V) | 编码 missense | 保护、罕见 | ESM-2 VEP + Boltz | Diogo 2015 [DOI](https://doi.org/10.1371/journal.pone.0122271) |
| rs12720356 | p.Ile684Ser (I684S) | 编码 missense | 保护 | ESM-2 VEP | Diogo 2015 [DOI](https://doi.org/10.1371/journal.pone.0122271) |
| rs2304256 | p.Val362Phe (V362F) | 编码 missense | 与自身免疫关联（常与上共分析） | ESM-2 VEP + SpliceAI（外显子剪接增强子） | GWAS（待校验*） |

\* **待校验**：标注"GWAS（待校验）"的 rsID 系从公开文献知识整理，运行前应用 GWAS Catalog / Open Targets Genetics 对齐最新 IBD 关联与精细定位（fine-mapping）后再喂 PROTO，避免过时/错配。R381Q、P1104A、A928V、I684S 已由上引文献直接支持。

### PROTO 运行方案（伪流程）
```
for variant in variants.tsv:
    if variant.type == "coding_missense":
        esm2_vep(seq_ref, seq_alt)        # 有害性/功能改变
        boltz(structure_ref, structure_alt)  # 结构/激酶域(可选)
    else:  # 非编码
        evo2_alphagenome(window±)         # 表达/调控影响
        enformer_borzoi(window)           # 染色质可及性/增强子
    spliceai_pangolin(variant)            # 所有变异查剪接改变
=> 输出每个变异的"功能影响分"，聚合成基因型层"通路激活/依赖度"特征，进入应答模型
```

---

## 4. 三角结论

- **IL-23 与 JAK 类是差异化价值最高的方向**：应答预测近乎空白，但**遗传学（IL23R LoF、TYK2 部分抑制变异）提供了"通路可预测性"的生物学基础**——意味着值得从数据里挖，而非无据空想。
- **PROTO 在 IL23R 上贡献最独特**：该基因非编码/调控变异多，正是 Evo2/AlphaGenome/SpliceAI 的用武之地，可产出"基因型层的 IL-23 通路激活度"这一新特征。
- **JAK 类同时要做 safety biomarker**（VTE/肿瘤风险分层），与 predictive 并行。
- **落地形态**：把 Geneformer 挖出的细胞态签名（组织）+ PROTO 基因型特征（血液可测）+ 血清蛋白，组合成"IL-23 应答概率""JAK 应答概率"两张评分卡，补齐用途 A 缺失的两类机制。
