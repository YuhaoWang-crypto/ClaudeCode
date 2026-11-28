# 分子 QTL 因果定位拓展:共定位能把 GWAS 信号落到"正确的因果基因"吗?

这是对第一轮分子 QTL 工作(cis-eQTL 复现 / TWAS / SuSiE·coloc 模拟 / SORT1 单位点)的拓展,直击这一行最核心的痛点:
**在一个 GWAS 位点上,eQTL 共定位到底能不能唯一指认已验证的因果基因?**

用金标准:LDL 胆固醇的 15 个已验证因果基因(多为药靶/孟德尔基因:PCSK9、SORT1、APOB、HMGCR、LDLR、APOE、CETP、LPA…)。
在每个位点,对**区域内所有带 eQTL 的基因**跑共定位(肝 + 全血 eQTL vs LDL GWAS),看能否落到已知因果基因。

结论先行:**eQTL 共定位在 15 个已验证位点里只正确锁定 2 个(SORT1、LIPC)**。失败分两类——① 半数因果基因在可及组织里根本没有可检测的 cis-eQTL;② 剩下有 eQTL 的,大多数因果变异是**蛋白编码/功能变异(错义/失活),不是表达改变**,所以表达 eQTL 抓不到(PCSK9、LDLR、APOB、LPA);而真能共定位的 SORT1 又同时点亮 3 个基因(多基因歧义)。

## 数据(全部公开,远程 tabix)

| 角色 | 来源 |
|---|---|
| GWAS | LDL 胆固醇 GWAS(GCST90239658,GRCh38 harmonised,N≈1.3M) |
| eQTL | eQTL Catalogue:GTEx 肝(QTD000266,n=208)、全血(QTD000356,n=670) |
| 因果基因真值 | 15 个 LDL 已验证因果基因(药靶 + 孟德尔) |
| 共定位 | coloc.abf;每位点区域内每个基因 vs GWAS(共 1,256 次检验) |

## 结果

### 各位点能否落到因果基因(见图 A/C)

| 结果 | 位点数(/15) |
|---|---|
| 共定位且命中因果基因 | **2**(SORT1、LIPC) |
| 因果基因有 eQTL 但不共定位 | 5(PCSK9、LDLR、LPA、LCAT、LDLRAP1) |
| 因果基因在可及组织无可检测 eQTL | 8(HMGCR、APOB、NPC1L1、CETP、APOE、ABCG8、ABCA1、PNPLA3) |

### 为什么失败 —— 两个各自重要的原因

1. **组织/功效不足(8/15)**:HMGCR(eQTL minP 2e-3)、APOB、NPC1L1 等药靶在肝(仅 208 例)/血里的 cis-eQTL 太弱或测不到。相关组织 + 足够样本量是硬前提;可及的 eQTL 资源往往不满足。
2. **因果机制不在"表达"上(5/15)**:PCSK9、LDLR、APOB、LPA 有 cis-eQTL,但 GWAS 因果变异是蛋白编码/失活/剪接变异,通过改变**蛋白功能**而非 mRNA 水平起作用——表达 eQTL 与 GWAS 不共定位是**正确的**,但这意味着 eQTL 共定位会系统性漏掉这类(恰恰是最好的药靶)基因。

### 即便共定位成功,也不唯一(见图 B)

SORT1 位点:SORT1(因果)、PSRC1、CELSR2 三个基因 PP.H4 都 >0.8——它们共享同一个强 eQTL/处在同一 LD 块,共定位无法区分哪个是因果基因。这正是第一轮 SORT1 结论的系统化确认。

## 对这一行的意义与建议

1. 共定位结果应表述为"信号共享的候选基因",**不能**当作"因果基因已锁定";同一位点多基因阳性要如实并列。
2. 交付因果定位必须组合证据:① 匹配的、足够功效的组织 eQTL;② pQTL / sQTL(抓蛋白、剪接层面);③ 对蛋白编码/失活驱动的基因,直接看编码变异与功能(coding-aware),而非只看表达 eQTL;④ 精细定位 + 功能实验兜底。
3. 对"药靶类"基因(PCSK9/LDLR/APOB)尤其警惕:它们最可能通过蛋白功能起效,纯 eQTL 共定位会漏掉——这也解释了为什么领域转向 pQTL 与编码变异 MR。

## 限制

- 只用了两个可及组织(肝 n=208、血),功效有限;更大的肝/肠 eQTL 或 eQTLGen(n≈3万血)会提高"有 eQTL"的比例,但不改变"蛋白机制基因漏检"和"多基因歧义"两个结构性结论。
- coloc 的单因果变异假设在多信号位点会失真(第一轮模拟已量化);此处用标准 coloc.abf。
- GWAS 的 N 用中位近似;eQTL 用 all.tsv 的边际统计,未做条件分析/SuSiE-coloc(后者能部分缓解多基因歧义,是自然的下一步)。
- "因果基因"真值取文献共识的 LDL 药靶/孟德尔基因;个别位点可能有多个真因果基因。
- 单一性状(LDL);其他性状的因果机制构成(更偏调控 vs 编码)会不同。

## 复现

```bash
pip install pandas numpy; Rscript -e 'install.packages(c("coloc","jsonlite"))'  # + tabix
export CW_DIR=...
python 01_harvest.py && Rscript 02_coloc.R && python 03_eval.py && python 04_fig.py
```
数字在 `results/summary.json`、`results/diagnosis.csv`、`results/per_locus.csv`,图 `results/fig1_causal_gene.png`。
