# IBD 用途 A（Predictive / 选药 CDx）— 初始候选 Panel 与证据卡片

> 对应 `biomarker-discovery-autoimmune-IBD.md` 的用途 A：**预测病人对哪一类机制生物制剂应答**（anti-TNF / IL-23·IL-12 / α4β7 整合素 / JAK）。
>
> 本文件是**起点候选清单**，由三引擎细化：**Geneformer**（黏膜 scRNA 嵌入 + in-silico 扰动）、**PROTO**（用 Evo2 / AlphaGenome / Enformer·Borzoi / SpliceAI 做候选基因的**调控变异功能验证**——注意 PROTO 是核酸/调控层，不是血浆蛋白）、**ESM-2 VEP**（编码变异有害性零样本打分）。
>
> **证据来源与署名**：以下候选的文献证据均来自 PubMed 检索（2026-07）。According to PubMed，每条引用附 DOI 链接。证据强度按"直接 IBD 应答预测证据 / 机制或旁证 / 缺口"分级，矛盾证据显式标注（证据三角）。

---

## 0. 方法（可复现）

1. **系统检索**（PubMed，按机制类别）：anti-TNF（OSM/OSMR/TREM1）、ustekinumab/risankizumab（IL-23/IL23R）、vedolizumab（α4β7/MAdCAM）、JAK（tofacitinib/upadacitinib）；辅以血清蛋白组（Olink）与单细胞检索。
2. **证据三角**：同一候选跨多队列/多模态交叉；显式记录**支持 vs 反驳**证据、样本量、队列类型（成人/儿科、CD/UC）、样本可及性（黏膜 vs 血清）。
3. **分级**：★★★ 有直接 IBD 应答预测证据（多队列）；★★ 单队列或强机制旁证；★ 假设/缺口（新颖性机会）。
4. **交棒计算引擎**：每张卡片末尾写明 Geneformer / PROTO / ESM-2 的具体下一步。

---

## 1. 证据卡片模板（可复用 schema）

> 每个候选一张卡；字段固定，便于后续自动填充与审计。

```yaml
candidate_id:            # 如 A-TNF-01
symbol:                  # 基因/蛋白官方符号 (HGNC)，如 OSM
name:                    # 全称
mechanism_class:         # anti-TNF | IL23/IL12 | alpha4beta7 | JAK | cross
biomarker_role:          # Predictive (对哪类药应答)
direction:               # 高表达/变异 → 应答 or 无应答
sample_type:             # mucosal biopsy | serum/plasma | stool | blood cells
assayability:            # 现成检测? (IHC/qPCR/ELISA/Olink/测序)  + 落地可及性评级
evidence:
  supporting:            # [{claim, cohort, n, disease, PMID, DOI, strength}]
  contradicting:         # 反驳/不一致证据 (证据三角关键字段)
  mechanistic:           # 机制支持 (通路/靶点)
evidence_grade:          # ★ / ★★ / ★★★
novelty:                 # 已知标志物 / 半新 / 新 (决定 IP 与差异化价值)
paired_drug:             # 拟配对药 (CDx 目标)
compute_nextstep:
  geneformer:            # 嵌入分类? 扰动验证驱动性?
  proto:                 # 哪些变异做 Evo2/AlphaGenome 调控打分? 剪接?
  esm2_vep:              # 哪些编码变异做有害性打分?
validation_cohort:       # 拟用哪个真实队列/试验 (NCT)
risk_flags:              # 批次效应/可及性/矛盾证据等
```

---

## 2. 初始候选 Panel

### 机制类 1：anti-TNF 应答预测

#### A-TNF-01 — OSM / OSMR ★★★（旗舰候选）
- **角色/方向**：治疗前**黏膜及血清 OSM 高** → **anti-TNF 无应答**风险高。
- **样本**：黏膜活检（原始证据）+ 血清（Olink 可测，利于落地）。
- **支持证据**：
  - West 等，*Nat Med* 2017：>200 IBD 患者（含 infliximab、golimumab 两个 Ph3 队列），治疗前 OSM 高表达强烈关联 anti-TNF 失败；~40% 原发无应答背景。[DOI](https://doi.org/10.1038/nm.4307)
  - Jongsma 等，*J Crohns Colitis* 2023（TISKids RCT，Olink 92 蛋白）：治疗前 **CD-hi** 聚类（含 **OSM**、TNFSF14、HGF、TGF-α 等高值）预测一线英夫利昔 52 周缓解率显著更低（CD-hi 24% vs CD-lo 58%）。[DOI](https://doi.org/10.1093/ecco-jcc/jjad049)
- **矛盾证据（证据三角）**：
  - Ezirike Ladipo 等，*J Pediatr Gastroenterol Nutr* 2021：**儿科**队列中 OSM **不**预测应答 → 提示年龄/队列依赖，落地需按人群分层验证。[DOI](https://doi.org/10.1097/MPG.0000000000003201)
- **机制支持**：Kokkotis 等，*Inflamm Bowel Dis* 2024：OSM 驱动肠道上皮下肌成纤维细胞促炎表型（IL-6、ICAM1、趋化因子），解释 anti-TNF 抵抗的基质轴。[DOI](https://doi.org/10.1093/ibd/izae098)
- **evidence_grade**：★★★（成人）/ 需人群分层。**novelty**：已知，但可**扩展为多机制血清 panel** 并补齐 IL-23/JAK 缺口。
- **compute_nextstep**：
  - Geneformer：在应答/无应答黏膜 scRNA 上验证 OSM+OSMR+ 基质模块是否驱动"炎症态"，in-silico 敲低 OSMR 看细胞态是否回落。
  - PROTO：对 *OSM/OSMR* 位点的 IBD GWAS/表达 QTL 变异做 Evo2+AlphaGenome 调控打分（是否改变增强子活性/表达）。
  - ESM-2 VEP：OSMR 错义变异有害性打分。
- **validation_cohort**：NCT06227910（Takeda 双靶 Ph3）、NCT07181525（预测模型队列）。

#### A-TNF-02 — TNFSF14 (LIGHT) / HGF / TGF-α ★★（血清可测，成套）
- **角色/方向**：治疗前**血清高值**（与 OSM 同属 CD-hi 聚类）→ 一线英夫利昔难缓解。
- **样本**：血清（Olink 92-plex，已验证可测）→ **落地性强**。
- **支持证据**：Jongsma 等，*J Crohns Colitis* 2023（同上 CD-hi 聚类的核心蛋白）。[DOI](https://doi.org/10.1093/ecco-jcc/jjad049)
- **evidence_grade**：★★（单 RCT，儿科 CD；需成人/UC 外部验证）。**novelty**：半新（作为组合标志物）。
- **compute_nextstep**：PROTO 对 *TNFSF14/HGF/TGFA* 调控变异做 AlphaGenome 表达影响打分；Geneformer 看这三者所在细胞模块在应答分层里的差异。
- **validation_cohort**：NCT07132333（Molecular Inflammation Board 多组学）。

#### A-TNF-03 — TREM1 ★（先验候选，引用待下一轮确认）
- 说明：TREM1（外周血/黏膜）在多个报道中与 anti-TNF 应答相关，但本轮检索未取到干净的 IBD-应答 PMID。**标记为待验证**，下一轮用 `search_articles("TREM1 anti-TNF response inflammatory bowel disease")` + `get_article_metadata` 补齐后再定级。**不在未确认前给出引用。**

---

### 机制类 2：α4β7 整合素（vedolizumab）应答预测

#### A-ITG-01 — α4β7 (ITGA4/ITGB7) + MAdCAM-1 (MADCAM1) 轴 ★★★
- **角色/方向**：治疗前**黏膜 α4β7+ 淋巴细胞 / MAdCAM+ 微静脉多**、外周 **CD4 T 细胞对 MAdCAM-1 黏附强** → vedolizumab **应答好**（机制自洽：药物即阻断 α4β7–MAdCAM）。
- **样本**：黏膜（IHC/流式）+ 外周血功能试验（黏附）→ 血液端可探索落地。
- **支持证据**：
  - Roosenboom 等，*Inflamm Bowel Dis* 2023：黏膜 α4β7+ 淋巴细胞与 MAdCAM+ 微静脉预测 UC 对 vedolizumab 应答。[DOI](https://doi.org/10.1093/ibd/izad123)
  - Allner 等，*BMC Gastroenterol* 2020：基线**动态 CD4 T 细胞–MAdCAM-1 黏附**与 UC 临床应答相关。[DOI](https://doi.org/10.1186/s12876-020-01253-8)
  - Holmer 等，*Ther Adv Gastroenterol* 2020：多个 biomarker 与 CD 对 vedolizumab 的临床/内镜结局相关。[DOI](https://doi.org/10.1177/1756284820971214)
- **evidence_grade**：★★★（机制直接、多队列）。**novelty**：已知机制，可**量化成可及性更好的血液检测**。
- **compute_nextstep**：
  - Geneformer：应答者 vs 无应答者黏膜 scRNA 中 gut-homing 淋巴细胞态（ITGA4/ITGB7/CCR9）丰度与状态；扰动验证。
  - PROTO：*ITGA4/ITGB7/MADCAM1* 调控变异 AlphaGenome/Enformer 打分。
  - ESM-2 VEP：整合素亚基错义变异。
- **validation_cohort**：NCT05428345（vedolizumab 上市后监测）、NCT06227910。

---

### 机制类 3：IL-23 / IL-12（ustekinumab / risankizumab）应答预测 — ★ 高新颖性缺口

- **现状**：本轮检索中，IL-23 类**直接的 IBD 应答预测标志物证据薄弱**——命中多为银屑病（如 risankizumab 的 UltIMMa Ph3，Gordon 等 *Lancet* 2018，[DOI](https://doi.org/10.1016/S0140-6736(18)31713-6)，**属银屑病、非 IBD 应答预测**，仅作机制背景）或定位/换药策略（D'Amico 等，*J Crohns Colitis* 2022，[DOI](https://doi.org/10.1093/ecco-jcc/jjac011)；Latras-Cortés 等，*Dig Dis Sci* 2025，anti-TNF 失败后 ustekinumab 疗效，[DOI](https://doi.org/10.1007/s10620-025-08978-0)）。
- **判断**：**这是最大的新颖性机会**——谁能先做出 IBD 里"IL-23 通路应答"的 predictive 标志物，价值最高。
- **候选假设（待数据验证）**：IL23R、IL12B、IL22、JAK2/TYK2 通路下游黏膜签名；治疗前 IL-23 轴激活程度作为应答预测。
- **compute_nextstep（重点跑）**：
  - Geneformer：在 ustekinumab/risankizumab 应答分层的黏膜 scRNA 上，无监督找 IL-23 响应细胞态；in-silico 扰动 IL23R/IL12B 验证驱动性。
  - PROTO：*IL23R*（含已知 IBD 保护性变异如 R381Q 区域）、*IL12B*、*TYK2* 的调控与剪接变异，用 Evo2/AlphaGenome/SpliceAI 打分——**PROTO 在这一类最有独特贡献**（IL23R 遗传学丰富）。
  - ESM-2 VEP：IL23R/TYK2 编码变异有害性。
- **validation_cohort**：NCT05387031（ustekinumab 狭窄型 CD 前瞻）、NCT07071519（risankizumab 儿科 UC Ph3）。

---

### 机制类 4：JAK（tofacitinib / upadacitinib）应答预测 — ★ 空白缺口

- **现状**：检索未见 UC 中已验证的 JAK **应答预测**标志物；命中为其他适应症（强直性脊柱炎 Ph3，Deodhar 等 *Ann Rheum Dis* 2021，[DOI](https://doi.org/10.1136/annrheumdis-2020-219601)）、个案与**安全性**（肿瘤风险 meta，Bezzio 等 *Cancers* 2023，[DOI](https://doi.org/10.3390/cancers15082197)）。
- **判断**：应答预测**几乎空白**；但 JAK 有明确**安全性 biomarker**需求（VTE/肿瘤风险分层）——可同时布局 predictive + safety。
- **候选假设（待验证）**：治疗前黏膜 JAK/STAT 激活签名、干扰素应答模块；外周免疫细胞 STAT 磷酸化。
- **compute_nextstep**：Geneformer 找 JAK/STAT 与 IFN 响应细胞态；PROTO 对 JAK1/TYK2 调控变异打分；ESM-2 VEP 编码变异。
- **validation_cohort**：NCT06227910（含 upadacitinib 臂）。

---

## 3. 跨候选整合与三角结论

- **可直接进入验证的最强候选**：A-TNF-01（OSM/OSMR，★★★，但需人群分层）、A-ITG-01（α4β7/MAdCAM 轴，★★★，机制自洽）。
- **落地性最好**：A-TNF-02（血清 Olink 三蛋白组合）——血清可重复采样，最接近 CDx 检测形态。
- **新颖性最高（差异化价值）**：IL-23 类与 JAK 类应答预测（当前近空白）——**建议把最多的 Geneformer/PROTO 算力压在这两类**，因为这里最可能挖出"新且有用"的 biomarker。
- **证据三角要点**：OSM 的成人阳性 / 儿科阴性冲突，说明任何候选都必须**按 CD/UC、成人/儿科、黏膜/血清分层验证**，不能一概而论；批次/中心效应是多组学最大陷阱。
- **组合而非单点**：最终 A 的产出应是"每类机制一个应答概率"的**多标志物评分卡**（如：OSM+TNFSF14+HGF → anti-TNF 无应答概率；α4β7/MAdCAM 轴 → vedolizumab 应答概率；IL-23 签名 → ustekinumab 应答概率），而非单基因。

---

## 4. 下一轮工作（明确 TODO）

1. 补齐 **A-TNF-03 TREM1** 与单细胞奠基文献（Crohn GIMATS 病理模块、UC 炎症相关细胞态）的**确证引用**（`search_articles` + `get_article_metadata`），未确认前不写入引用。
2. 对每个 ★★★/★★ 候选，跑 **PROTO** 的调控变异打分（IL23R 优先，遗传学最丰富）。
3. 用 **Geneformer** 在公开 IBD 治疗前黏膜 scRNA（应答分层）上，把 IL-23 类与 JAK 类的"响应细胞态"从零挖出来（这是最大增量）。
4. 把本 panel 落成机器可读的 `candidates.yaml`（按第 1 节 schema），接自动化打分与验证追踪。

---

### 署名
本文件的文献证据均来自 **PubMed**（According to PubMed）。所有被引用文章的 DOI 链接已在正文内逐条给出。临床试验编号来自 ClinicalTrials.gov（2026-07 检索，在招状态）。
