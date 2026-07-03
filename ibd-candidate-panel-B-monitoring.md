# IBD 用途 B（Response / Monitoring）— 初始候选 Panel 与证据卡片

> 对应 `biomarker-discovery-autoimmune-IBD.md` 的用途 B：**疗效（PD/Response）+ 治疗程度（Monitoring）**——药有没有起作用、疾病负荷怎么变、能否减少肠镜。
>
> **与 A 的引擎重心不同**：A 靠"基线预测"；B 靠"纵向趋势"。B 的主力是 **Geneformer（残留病/炎症细胞态）+ 血清蛋白纵向签名**；PROTO（核酸/调控层）与 ESM-2 在 B 中作用较小。
>
> **署名**：文献证据来自 PubMed（According to PubMed），逐条附 DOI。分级同 A：★★★ 多队列直接证据 / ★★ 单队列或强旁证 / ★ 假设或缺口。

---

## 0. 临床缺口

现有监测手段各有短板（据 PubMed，Sakurai 综述 *Digestion* 2022，[DOI](https://doi.org/10.1159/000527846)）：
- **肠镜**：金标准但侵入、不能频繁做；
- **粪钙卫蛋白 (FCP)**：最好的无创指标，但依从性差（粪便采样）、个体波动大；
- **CRP**：不特异，且对"疗效"判读差；
- **缺口**：需要**血液可测、机制特异、纵向优于 FCP/CRP** 的 PD/monitoring 标志物，实现 treat-to-target 的紧密监测。

---

## 1. 证据卡片（沿用 A 的 schema，新增纵向字段）

新增字段：`kinetics`（多久出现变化）、`target_correlation`（与内镜/组织学缓解相关性）。

---

## 2. 初始候选 Panel

### B-01 — 粪钙卫蛋白 (FCP / S100A8·S100A9) ★★★（现有金标准对照）
- **角色**：Monitoring + 复发预测。**方向**：升高→黏膜炎症/复发。
- **样本**：粪便。**可及性**：商品化 ELISA，但采样依从性是短板。
- **支持证据**：Cassinotti 等 *Clin Gastroenterol Hepatol* 2021——硫唑嘌呤停药后，**FCP 阳性是复发的唯一显著预测因子**（UC HR 3.3；CD HR 4.5），[DOI](https://doi.org/10.1016/j.cgh.2021.06.014)。
- **evidence_grade**：★★★。**novelty**：已知——作为**基准/弱标签**，新候选须证明优于它。
- **compute_nextstep**：Geneformer 关联 FCP 高低与残留中性粒细胞/上皮损伤细胞态；作为纵向模型的锚。

### B-02 — 血清内镜愈合指数 EHI（13 蛋白 panel）★★（血清替代内镜）
- **角色**：Monitoring（血清替代黏膜炎症判读）。
- **样本**：血清。**可及性**：已商品化（Monitr/EHI），落地性强。
- **支持证据**：Alsoud 等 *Am J Gastroenterol* 2023——EHI 在 UC 中随内镜评分变化，检出 MES 0-1 的 AUC 77.8%，与 FCP 相当；10 分下降对应 MES 下降 1 分的 OR 提升 89%，[DOI](https://doi.org/10.14309/ajg.0000000000002518)。
- **evidence_grade**：★★（原为 CD 开发，UC 已验证）。**novelty**：已知 panel——可作为**多蛋白监测的模板**，用我们的方法扩展/优化。
- **compute_nextstep**：作为血清 PD panel 的起点，Geneformer 找其组成蛋白对应的组织细胞态来源。

### B-03 — LRG1（亮氨酸富集 α2 糖蛋白）★★
- **角色**：Monitoring（内镜活动度）。**方向**：升高→内镜活动。
- **样本**：血清。**可及性**：ELISA。
- **支持证据**：Sakurai 综述——**LRG 诊断内镜疾病活动度优于 CRP**，[DOI](https://doi.org/10.1159/000527846)。
- **evidence_grade**：★★。**novelty**：半新——比 CRP 好但未广泛落地，值得纳入组合。
- **compute_nextstep**：纳入血清 PD panel；Geneformer 关联 LRG1 来源细胞（中性粒/上皮）。

### B-04 — OSM（血清/黏膜，跨用途复用）★★
- **角色**：Response/Monitoring。**方向**：治疗中下降→应答。
- **样本**：黏膜/血清。**跨用途**：A 里是 predictive（基线高→anti-TNF 无应答），B 里可作 PD（治疗中动态）。
- **支持证据**：Sakurai 综述列为新兴监测候选，[DOI](https://doi.org/10.1159/000527846)；机制见 A-TNF-01（West 2017 [DOI](https://doi.org/10.1038/nm.4307)）。
- **evidence_grade**：★★（作为 PD 需纵向验证）。**novelty**：半新（PD 用途）。
- **compute_nextstep**：Geneformer 扰动验证 OSM 轴随治疗回落；纵向血清测量。

### B-05 — CRP ★★（基准，须知其局限）
- **角色**：Monitoring（急性期）。**样本**：血。
- **支持/局限证据**：Sakurai 综述（不特异）[DOI](https://doi.org/10.1159/000527846)；Liu 等 *J Cancer* 2023 在免疫相关性结肠炎中显示 **CRP 对疗效/复发判读不可靠**（虽反映严重度），[DOI](https://doi.org/10.7150/jca.84261)——提示 CRP 只能做辅助。
- **evidence_grade**：★★（作基准）。**novelty**：已知，低差异化。

### B-06 — 新兴/假设候选（★ 缺口）
- **PBMC CDC42**：Liu 等 *J Clin Lab Anal* 2022——CDC42 随英夫利昔应答上升，可监测 UC 应答，[DOI](https://doi.org/10.1002/jcla.24477)（小样本，需验证）。
- **抗 αvβ6 抗体、PGE-MUM（尿）、循环 miRNA**：Sakurai 综述列为前景候选，[DOI](https://doi.org/10.1159/000527846)。
- **判断**：这些是 B 的**新颖性机会**——用 Geneformer 从纵向 scRNA 挖"残留病细胞态"，再找其血清代理，是最可能超越 FCP/CRP 的路径。

---

## 3. 非分子锚点（监测策略）
- **肠道超声 (IUS, 肠壁厚度 BWT)**：Dolinger 等 *J Crohns Colitis* 2025——IUS（BWT≤2.7mm）+ FCP 早期应答预测 treat-to-target 内镜结局，[DOI](https://doi.org/10.1093/ecco-jcc/jjaf075)。**不是分子 biomarker，但可作为监测队列的客观参照终点**，与血清 PD panel 联合。

---

## 4. 三角结论与产出形态

- **B 的最强现成组合**：FCP（B-01）+ 血清 EHI/LRG1（B-02/03），配 IUS 客观终点。
- **差异化机会（★ 缺口）**：血液可测、机制特异、纵向优于 FCP 的 **PD panel**——靠 Geneformer 挖残留病细胞态 + 找血清代理。
- **跨用途复用**：OSM 在 A（predictive）与 B（PD）双重身份；α4β7 轴基线预测 + 治疗中监测。
- **产出**：一张"**疾病负荷趋势评分**"——FCP + 血清多蛋白（EHI/LRG1/OSM）纵向 + Geneformer 残留细胞态，目标是**减少肠镜频次、早期识别继发失应答**（指导换药）。

---

## 5. 下一轮 TODO
1. 取纵向 IBD 队列（治疗前/中血清 + 配对内镜），把 EHI/LRG1/OSM/FCP 建成纵向模型。
2. Geneformer 在纵向黏膜 scRNA 上定义"残留病细胞态"，映射到血清可测代理。
3. 把 B 候选并入 `candidates.yaml`（新增 use_case: B）与 `score_candidates.py` 追踪。
4. 与真实监测终点（IUS/内镜）对齐验证时点。

---

### 署名
文献证据来自 **PubMed**（According to PubMed），DOI 链接逐条见正文。
