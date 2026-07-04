# 具体案例：A375（BRAF^V600E 黑色素瘤）+ vemurafenib 加药臂

把 `07/08/09` 的通用流程落到一个**完全指定**的体系，剂量与临床信息用
ChEMBL v34 / ClinicalTrials.gov 实查数据填充。你若换成自己的体系，照此结构替换
即可。

> **数据出处：** 下方 IC50/机制来自 **ChEMBL v34**（vemurafenib = CHEMBL1229517，
> 靶点 BRAF = CHEMBL5145，UniProt P15056 V600E）；临床剂量与联用格局来自
> **ClinicalTrials.gov**。均为本会话实查，非记忆值。

---

## 1. 体系规格

| 项 | 设定 |
|----|------|
| 细胞系 | **A375**（人皮肤黑色素瘤），**BRAF^V600E**，倍增时间约 17–24 h |
| 工程改造 | 先建 **A375-Cas9** 稳定株并验活性（见 `04_文库构建.md`） |
| 文库 | 激酶组 ~3,508 引导（518 激酶 × 6 + 400 对照，见 `05_筛选与测序.md`） |
| 加药臂 | **vemurafenib（PLX4032 / Zelboraf）**，BRAF^V600E 抑制剂 |
| 遗传背景（可选替代） | 若做基因型合成致死，可用 A375（BRAF^V600E）vs. BRAF-WT 黑色素瘤系做等基因对比 |

**为什么这个体系是好范例：** A375 对 BRAF 呈**癌基因成瘾**——BRAF 敲除本身即强
丢失，vemurafenib 施加的正是通路（MAPK）压力，能让"敲除后改变 MAPK/药物依赖"
的激酶显形。

---

## 2. 真实活性数据（ChEMBL v34，供设定剂量-反应范围）

| 测定层次 | IC50 | 说明 |
|----------|------|------|
| 生化 BRAF^V600E | **3.2 – 31 nM**（pChEMBL 7.5–8.5，多来源） | 靶点内在效力 |
| A375 细胞内 pERK 抑制 | **33 nM**（WB, 90 min）／**150 nM**（ELISA, 72 h）／190–260 nM（其他） | **靶点接合**层次 |
| 抗增殖 GI50（A375） | 文献多在 **~0.5–1 µM**（本库未含，需自测） | 决定筛选剂量的关键 |

> 关键区分：**靶点接合 IC50（~0.1 µM）≠ 抗增殖 GI50（~0.5–1 µM）**。筛选压力由
> 抗增殖曲线决定，所以第 3 步的剂量-反应必须自测，不能直接搬 ChEMBL 的 pERK 值。

---

## 3. 剂量选择（把 `08` 第 1 步具体化）

1. **剂量-反应范围建议**：8–10 点，从 **10 nM 到 10 µM**（3 倍稀释，跨越细胞
   IC50 ~0.1 µM 上下两个数量级）。
2. **72 h CellTiter-Glo** 拟合 GI50；预期落在 ~0.5–1 µM 附近（以实测为准）。
3. **取 IC20–IC30** 作为筛选终浓度 `C_screen`——本体系多半落在 **~0.2–0.5 µM**。
4. **慢性校验**：以 `C_screen` 持续给药 10–14 天，使加药臂倍增速率约为 DMSO 臂
   的 70–80%，据此微调。

---

## 4. 具体样本与 drugZ 矩阵

样本命名（对应 `templates/sample_sheet.csv`，把 `DRUG` 记为 `VEM`）：

```
plasmid | T0_1..3 | DMSO_1..3 | VEM_1..3        （10 个样本）
```

### A 线 · 必需性（溶剂 vs. T0，BAGEL2）

```bash
BAGEL.py fc -i kinome_screen.count.txt -o kinome_fc -c T0_1,T0_2,T0_3
BAGEL.py bf -i kinome_fc.foldchange -o results/essentiality_bagel.tsv \
  -e CEGv2.txt -n NEGv2.txt -c DMSO_1,DMSO_2,DMSO_3
```

### B 线 · 合成致死（vemurafenib vs. 溶剂，drugZ）

```bash
python drugz.py -i kinome_screen.count.txt -o results/synlethal_drugz.tsv \
  -c DMSO_1,DMSO_2,DMSO_3 -x VEM_1,VEM_2,VEM_3
# normZ < 0：敲除使 A375 对 vemurafenib 更敏感（增敏/合成致死）
# normZ > 0：敲除赋予耐药
```

---

## 5. 本体系的预期命中（兼作阳性对照与判读先验）

### 阳性对照（先确认筛选成立）

- **必需性臂**：泛必需激酶（PLK1、CDK1、AURKB）+ **BRAF 本身**强丢失——A375
  BRAF 成瘾，BRAF 敲除即使无药也显著丢失。若这些不丢，先别看新命中。
- **加药臂**：MAPK 通路增敏方向应见信号（见下）。

### 增敏 / 合成致死（normZ < 0，加药下更丢失）——候选激酶锚点

方向可靠、具体基因以你的数据+文献为准；这些类别正是临床联用的生物学基础：

| 类别 | 候选激酶 | 临床对应 |
|------|----------|----------|
| MAPK 纵向加压 | MAP2K1/MAP2K2（MEK）、MAPK1/MAPK3（ERK） | **BRAFi + MEKi**（vemurafenib+cobimetinib 已获批） |
| 适应性 RTK 反馈（ERK 反弹） | EGFR、ERBB2/3、FGFR、IGF1R、PDGFRB、MET | 旁路重激活是耐药主轴 |
| 平行生存通路 | PIK3CA、AKT1/2/3、PDPK1、MTOR | BRAFi + PI3K/AKT 抑制协同 |
| 细胞周期 | CDK4、CDK6 | **BRAFi + CDK4/6i** 联用 |
| 生存/黏附 | SRC 家族、PTK2(FAK) | 维持存活信号 |

> 注意 MEK/ERK：它们在**溶剂臂也基础必需**，drugZ 用"加药 vs. 溶剂"正好扣掉基础
> 必需、只留药物特异的**增量**——这也是不能用"加药 vs. T0"的原因。

### 耐药 / 富集（normZ > 0，敲除→存活）——预期稀疏，需管理预期

**方法学要点：** 临床 BRAFi 耐药多由**激活型**事件驱动（RTK 上调、MAP3K8/COT
过表达、CRAF 代偿），这些是"获得功能"，**敲除筛选查不到**；而丢失型耐药事件
（NF1、PTEN、CIC、MED12）多为**非激酶**。因此**激酶组敲除筛的耐药方向命中会
偏少**——这是设计的固有限制，不是失败。若要系统查耐药激活型机制，应改用
**过表达/ORF 或 CRISPRa** 筛选，而非敲除。

---

## 6. 临床相关性（ClinicalTrials.gov 实查）

- vemurafenib 联用 MEK 抑制剂 cobimetinib 的注册试验 **≥36 项**（如 NCT01656642、
  NCT04722575 等），BRAFi+MEKi 已成标准联用——正是本筛"增敏方向 = MAPK 反馈
  重激活节点"的临床印证。
- 说明：把筛选命中映射回**已获批/在研联用**，能快速判断哪些增敏命中具转化价值
  （已验证的 MEK/CDK4-6 vs. 新颖节点）。

---

## 7. 一页落地清单（本体系）

1. 建 A375-Cas9，验活性 → 克隆激酶组文库。
2. vemurafenib 剂量-反应（10 nM–10 µM）→ 定 `C_screen`（预期 ~0.2–0.5 µM）→ 慢性校验。
3. 转导（MOI 0.3，≥13.5M 细胞）→ 嘌呤筛选 → T0（3 重复）。
4. 分流 DMSO 臂 / vemurafenib@`C_screen` 臂，各 3 重复，全程 ≥3.5M 细胞，给药至 ~21 天。
5. 收样 → gDNA（~23 µg/样本）→ 两步 PCR → 测序（~44M 读数）。
6. BAGEL2（必需性，含 BRAF 阳性对照）+ drugZ（增敏方向 = MEK/ERK/PI3K/CDK4-6 先验）。
7. 命中映射到 BRAFi 联用格局做转化优先级排序。
