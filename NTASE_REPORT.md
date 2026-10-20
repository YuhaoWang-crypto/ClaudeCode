# 从「植物蛋白从头合成 RNA」提取 AI 方法 → 构建 skill → 跨门类扩展挖掘

**源论文:** Li, Yu, Li, Gu, Zhan & Shi. *De novo* RNA synthesis from free nucleotides
by a family of conserved plant proteins. *Vita* **2026**, 10.0078.

本文档分三部分:(一)论文里到底用了哪些 AI / 计算方法;(二)把这套方法抽象成一个
可复用的 skill;(三)用这个 skill 真的跑了一遍跨门类挖掘,并给出带诚实标注的候选表。

**标注规则**(沿用本仓库 `network-biomarker` 的纪律):

- ✅ **严格** — 确定性计算、真实测得的统计量、或对文献数值的精确复现。
- ⚠️ **假说** — 合理但未经实验验证的推断(候选蛋白、阈值类比、功能推测)。

**最重要的一条:结构相似性永远不是活性的证据。** 下面所有新提名的蛋白都是 ⚠️。

---

## 一、论文中使用的 AI / 计算方法

论文的实验部分很重(NMR、MS、晶体),但**发现**这三个酶的过程完全是计算驱动的。
提取出的方法清单:

| # | 方法 | 在论文中的具体用途 | 关键参数 / 结果 |
|---|---|---|---|
| 1 | **AlphaFold2 蛋白结构预测** | 提供三个植物物种的全蛋白质组预测结构 | 玉米、拟南芥、水稻;AFDB |
| 2 | **AlphaFold 蛋白结构数据库(AFDB)** | 被搜索的目标库 | 覆盖 2.14 亿条序列 |
| 3 | **DALI 结构同源搜索** | 以人 cGAS(PDB 4KM5)为 3D 查询去搜 AFDB | Z≥2.0 → 266;Z≥10.0 → 29 |
| 4 | **pLDDT 置信度** | 评估预测模型可靠性(Fig 1h 配色) | — |
| 5 | **结构叠合 / RMSD** | OsNTP4 模型 vs cGAS 晶体 | 全局 4.623 Å/1182 原子;**催化区 0.735 Å/66 原子** |
| 6 | **全对全结构 + 序列矩阵聚类** | 19 个候选的两两 Z-score 与序列一致性 → 分出 4 个亚群 | 亚群内 40–80%,对 cGAS 仅 ~15% |
| 7 | **催化残基保守性过滤** | 要求三个酸性残基(cGAS E225/D227/D319)保守;Glu↔Asp 可互换 | 29 → **19**(淘汰 10 个) |
| 8 | **最大似然系统发育(MEGA 7.0.26)** | 47 个 NTP 家族蛋白建树 | JTT 模型 + gamma 5 类;95% 位点覆盖;50 次 bootstrap;NNI → **5 个 clade** |
| 9 | **AlphaFold 复合物建模** | 预测 OsNTP4 结合 RNA 五核苷酸 / UTP,分析活性口袋静电 | Supplementary Fig S13 |
| 10 | **纳米孔测序 + 序列统计学** | 分析产物 RNA 的序列特征 | Shannon 熵、语言学复杂度(k=1–6)、嘌呤/嘧啶比、二核苷酸 SD、同聚物生存统计 P(Lmax≥k)、零阶置换零模型(1000 次重采样)、4096 个 6-mer 富集 + BH 校正、MEME |

**真正可迁移的核心是第 1、3、5、7、8 条**,串起来就是:

> 用**催化架构**当搜索键(而不是序列)→ 分层打分 → 用**化学**(催化残基)过滤 →
> 只看**活性位点**的几何 → 聚类建树找 clade → 输出可实验验证的候选。

这正是本次构建的 skill。

---

## 二、构建的 skill:`catalytic-scaffold-mining`

位置:`.claude/skills/catalytic-scaffold-mining/`

```
SKILL.md                    方法总览、六阶段、阈值规则、诚实标注纪律、输出契约
reference/methodology.md    逐阶段细节 + 论文原始数值作为校准参考 + 失败模式
reference/data-access.md    所有 API 配方(Foldseek web / AFDB / EBI FTP / UniProt / RCSB)
reference/validation.md     从候选到结论必须爬的五级对照阶梯
scripts/structure.py        PDB/mmCIF 解析、Kabsch 叠合、离群剔除的 core RMSD
scripts/catalytic.py        催化残基映射、化学类保守性、唯一且保序的分配、活性位点几何
scripts/fsweb.py            Foldseek 网页服务客户端(提交/轮询/取回/重算 TM-score)
scripts/mine.py             主驱动,web 与 local 两种模式
scripts/annotate.py         谱系(UniProt taxonomy)+ 已知酶类注释
scripts/calibrate.py        用已知阳性/阴性反推阈值,而不是照搬文献数字
```

### 六个阶段

| 阶段 | 回答什么 | 产出 |
|---|---|---|
| 0 定义查询 | 什么架构?哪几个残基干化学? | 结构 + 催化残基编号 |
| 1 结构扫描 | 这些蛋白质组里什么长得像? | DALI / Foldseek 打分 |
| 2 分层 | 哪些是信号而不是折叠噪声? | 宽松阈值(清点) + 严格阈值(入围) |
| 3 化学过滤 | 催化残基星座还在吗? | 通过**结构**比对映射,按化学类判定 |
| 4 几何核验 | 活性位点真能叠上吗? | 仅催化壳层的 core RMSD + 羧基间距 |
| 5 家族结构 | 几个亚家族?谁和谁是姐妹? | 全对全矩阵 + ML 建树 → clade |
| 6 分诊与标注 | 什么是证实的,什么是提名的? | 排序表 + ✅/⚠️ + 对照阶梯 |

### 相对论文做的三处方法学改进

这三处都是**实际踩到坑之后**加的,不是设计时想出来的:

**(a) 离群剔除的 core RMSD。** 活性位点旁边的柔性 loop 会把整壳层 RMSD 从 ~1 Å
拉到 >3 Å。实测 cGAS 自身命中(4KM5 晶体 vs 自己的 AlphaFold 模型,362 个位置映射
全部正确):

```
整壳层      3.31 Å / 53 对
core(4 Å 剔除后)  1.33 Å / 44 对(保留 83%)
逐残基偏差: 中位 0.62 Å,p90 2.26 Å,最大 11.63 Å
```

离群点落在残基 212–218 和 364–369,正是 cGAS 的 activation loop。**排序用
`core_rmsd`,并始终同时报告 `n_core`** —— 12 个残基上的 0.6 Å 和 40 个残基上的
0.8 Å 完全不是一回事。

**(b) 配准窗口(registration window)。** 在 ~10–15% 一致性下,结构比对会把催化
motif 整体错位 1–2 个残基。实测 cGAS(4KM5)比对 AtURT1(6L8K)时,E225 落进了
gap,D227 被映射到 D491。生物学上正确的对应是 D491/D493/D547(来自 motif
`...FGFPKS**D491**I**D493**V...` 与 `...GISC**D547**IC...`)。

所以加了 `--window N`,但配了三条保险:**偏移量必须显式报告**(`max_offset`、
残基列里的 `[+2]`)、`n_conserved_exact` 单独计数、**分配必须唯一且保序**。最后
一条是必需的 —— 第一版没有保序约束时,E225 和 D227 双双抢走了同一个 D491,看起来
像保守,其实是假象。

**(c) 参考结构可能是催化失活突变体。** AtURT1 的两个晶体结构 6L3F / 6L8K 都带
`SEQADV ... ALA A 547 ... ASP 547 ENGINEERED MUTATION`,即 D547A 失活突变。
不先查 `SEQADV` 就会"发现"野生型酶缺了自己的催化残基。

---

## 三、方法学验证(✅ 严格)

在拿这套流程去挖新东西之前,先证明它能重现已知的东西。

### 3.1 自身命中对照
以 4KM5 搜 AFDB,人 cGAS(Q8N884)被正确召回,362 个比对位置**零错配**,
TM-score 0.956,三个催化残基 E225→E225、D227→D227、D319→D319 精确对应。✅

### 3.2 独立推导出 AtURT1 的催化三联体
仅凭结构比对(序列一致性 11%),从 cGAS 的 E225/D227/D319 映射到
**AtURT1 D491 / D493 / D547**,与序列 motif 独立推出的位置完全一致;第三位在晶体
里是 Ala,正是上面那个工程突变。✅

### 3.3 精确重现论文的 OsNTP4 催化残基 ← 最强验证
以 AtURT1(6L8K)为查询搜 `afdb-proteome` + `afdb-swissprot`,在 816 个唯一命中里:

| 排名 | 蛋白 | core RMSD | 推导出的催化残基 | 论文报告 |
|---|---|---|---|---|
| 0 | **AtURT1** (O64642,查询蛋白自身) | 0.692 Å / 55 | D491/D493/D547 | — |
| 4 | **OsNTP10** (A0A0P0XSC2,已验证有活性) | 0.815 Å / 55 | D105/D107/D161 | — |
| 5 | **OsNTP4** (A0A0P0VE25,已验证有活性) | 0.894 Å / 55 | **D100 / D102 / D157** | **Asp100, Asp102, Asp157** |

论文正文写的是"conserved acidic residues **Asp100, Asp102 and Asp157**",且这三个
位点的 D→A 突变各自都会废掉活性。流程在不知道答案的情况下把三个残基**逐一精确
推了出来**,并把三个已验证的酶全部排进前 6。✅

### 3.4 阈值校准(而不是照搬)
用 `calibrate.py`,阳性 = 5 个植物 URT1 直系同源(含两个已验证的水稻酶),
阴性 = 古菌 CCA 添加酶(同超家族,但催化的是完全不同的反应):

```
core_rmsd  阳性:  0.692 – 0.894 (中位 0.787)   n_core 全部 = 55
           阴性:  1.202 – 1.267 (中位 1.222)   n_core 46–49
→ 可分;建议阈值 core_rmsd ≤ 1.05(间隔 0.31)
```
✅ 这个阈值是从数据里推出来的,不是从论文搬来的。

### 3.5 local 模式端到端(含 pLDDT 与侧链几何)
古菌 *Methanocaldococcus jannaschii* 全蛋白质组(1773 个模型)本地搜索耗时 **12 秒**。
催化壳层 pLDDT 92–97(很高),羧基间距可测。✅

---

## 四、跨门类扩展挖掘的结果

### 4.1 为什么换查询蛋白

论文用 cGAS 当查询是为了**找到**植物家族。要往外扩展,更合适的查询是**已被证明
有从头合成活性的酶本身** —— AtURT1。实测也支持这个选择:cGAS 搜 `afdb50` 的
1000 个命中几乎全是后生动物的 cGAS/cGLR(光 *Dreissena polymorpha* 一个物种就占
102 个),**一个植物或真菌蛋白都没够到**。

⚠️ **因此:被截断的搜索里某个类群缺席,不能当作该类群没有这个折叠的证据。**

### 4.2 筛选条件与结果

查询 AtURT1(6L8K,催化位点 491/493/547),库 = `afdb-proteome` + `afdb-swissprot`,
`--window 2`。1004 条命中 → 816 个唯一蛋白 → 要求**三个酸性残基全保守 +
core_rmsd ≤ 1.30 Å + n_core ≥ 50**,剩 **153 个**:

| 谱系 | 数量 | 其中已知酶类 | 无实质注释 |
|---|---|---|---|
| Metazoa 后生动物 | 78 | poly(A) polymerase 60、TENT 8、TUTase 6 | 4 |
| Viridiplantae 绿色植物 | 36 | poly(A) pol 20、TUTase 7、其他 8 | 1 |
| Fungi 真菌 | 24 | poly(A) polymerase 22、TUTase 1 | 0 |
| Discoba 眼虫类(锥虫/利什曼) | 8 | TUTase 2、poly(A) pol 1 | 5 |
| Amoebozoa 变形虫 | 5 | poly(A) polymerase 3 | 2 |
| Alveolata 囊泡虫(疟原虫) | 2 | poly(A) polymerase 1 | 1 |

把几何条件放宽到 core ≤ 2.5 Å,还会多出 **细菌 56 个**(主要是 GlnD 类
uridylyltransferase)和 **古菌 20 个**(主要是 CCA 添加酶)。

**结论(✅ 计算层面):支撑三个植物酶从头合成 RNA 的那个催化构型,不是植物特有的 ——
它在后生动物、真菌、眼虫类、变形虫、囊泡虫、细菌和古菌里都保守存在。**

### 4.3 一个事后验证:TENT4B

在非植物命中里,**人 TENT4B(Q8NDF8)core RMSD = 1.055 Å**,排在 117 个非植物
命中的第 13 位左右。

而论文自己在讨论里引用了参考文献 44(Nagpal et al., *Nucleic Acids Res.* 2026):
**人 TENT4B 被独立报道可以进行 primer- 与 template-independent 的 RNA 聚合**。

也就是说,这套从植物酶出发的扩展逻辑,把那个唯一已被实验证实具有同类活性的非植物
蛋白,排进了前列。✅ 这是对扩展方向的一次事后(prospective)验证。

### 4.4 结构域架构的独立佐证(✅)

三个已验证的植物酶与排名最前的跨门类候选,**共享同一对 Pfam 结构域**:

```
PF22600 + PF03828
  AtURT1 (O64642, 764 aa)      OsNTP4 (A0A0P0VE25, 367 aa)   OsNTP10 (A0A0P0XSC2, 372 aa)
  人 TENT4B (Q8NDF8, 572 aa)   人 TENT2/GLD2 (Q6PIY7, 484 aa)
  盘基网柄菌 Q54Y43 (1090 aa)   Q54TZ1 (920 aa)   C7G070 (938 aa)
  恶性疟原虫 A0A5K1K8X2 (858 aa)
```

这是一条和结构比对**彼此独立**的证据线。

### 4.5 ⚠️ A 组候选:已命名、可直接上手测的酶

这些不是"新蛋白"。恰恰相反 —— 论文最重要的隐含启示是:**AtURT1 本来就是已知的
RNA 3′ 尿苷酰转移酶,从头合成只是同一个催化骨架的另一种输出**。所以最省力的
实验预测是:**把这些早就被研究透的末端核苷酸转移酶,在"不加引物、不加模板"
的条件下再测一次**。

| core RMSD | 蛋白 | 物种 | 现有注释 |
|---|---|---|---|
| 0.945 | Gld-2 homolog A (Q9VD44) | 果蝇 | poly(A) RNA polymerase |
| 0.988 | Gld-2 homolog B (Q9VYS4) | 果蝇 | poly(A) RNA polymerase |
| 1.009 / 1.010 | GLD2 (Q5U315 / Q91YI6) | 大鼠 / 小鼠 | poly(A) RNA polymerase |
| **1.022** | **TENT2/GLD2 (Q6PIY7)** | **人** | 胞质 poly(A) 聚合酶 |
| **1.055** | **TENT4B (Q8NDF8)** | **人** | 末端核苷酸转移酶 ← 已被独立证实有该活性 |
| 1.082 | TENT4A (Q5XG87) | 人 | 末端核苷酸转移酶 |
| 1.103 | polynucleotide adenylyltransferase (Q5A2B2) | 白念珠菌 | — |
| 1.115 | **Cid13 (Q9UT49)** | 粟酒裂殖酵母 | poly(A) RNA polymerase |
| 1.162 / 1.177 | Pap1 / Pap2 (P48561 / P53632) | 酿酒酵母 | poly(A) RNA polymerase |
| 1.170 | **Cid11 (O74326)** | 粟酒裂殖酵母 | poly(A) RNA polymerase |
| 1.183 | RNA uridylyltransferase (A4HYQ1) | 婴儿利什曼原虫 | TUTase |
| 1.197 | TUT7 (Q5BLK4) | 小鼠 | 末端尿苷酰转移酶 |
| 1.197 | TUT3 (Q6B9Y8) | 布氏锥虫 | 末端尿苷酰转移酶 |
| 1.203 | TUT4 (Q5TAX3) | 人 | 末端尿苷酰转移酶 |
| 1.263 | **Cid1 (O13833)** | 粟酒裂殖酵母 | 末端尿苷酰转移酶(该家族的原型) |

⚠️ 全部为提名,均未验证。

### 4.6 ⚠️ B 组候选:几乎没有功能注释的蛋白

非植物、三残基全保守、core ≤ 1.30 Å / n_core ≥ 50,且注释里没有任何实质信息:

| core RMSD | 一致性 | 登录号 | 谱系 | 物种 | 现有注释 |
|---|---|---|---|---|---|
| 1.040 | 0.25 | A0A3P7FL62 | 后生动物 | 班氏吴策线虫 | Uncharacterized protein (片段) |
| **1.125** | **0.38** | **Q54Y43** | 变形虫 | 盘基网柄菌 | Uncharacterized protein |
| 1.146 | 0.26 | Q54TZ1 | 变形虫 | 盘基网柄菌 | Uncharacterized protein |
| 1.163 | 0.17 | A4I1Y4 | 眼虫类 | 婴儿利什曼原虫 | DNA polymerase sigma-like |
| 1.179 | 0.28 | A0A8R1U0D9 | 后生动物 | 旋盘尾丝虫 | CCHC-type domain-containing |
| 1.192 | 0.29 | Q09409 | 后生动物 | 秀丽隐杆线虫 | 4Fe-4S ferredoxin-type domain-containing |
| 1.252 | 0.24 | A0A5K1K8X2 | 囊泡虫 | 恶性疟原虫 | Nucleotidyltransferase, putative |
| 1.260 | 0.17 | Q57WQ6 | 眼虫类 | 布氏锥虫 | DNA polymerase sigma-like |
| 1.278 | 0.21 | Q57VP5 | 眼虫类 | 布氏锥虫 | DNA polymerase sigma, putative |
| 1.283 | 0.13 | Q38EI4 | 眼虫类 | 布氏锥虫 | Polymerase nucleotidyl transferase domain |
| 1.288 | 0.20 | Q4DSJ9 | 眼虫类 | 克氏锥虫 | DNA polymerase sigma, putative |

**最值得先看的是 `Q54Y43`(盘基网柄菌)** —— 对 AtURT1 序列一致性 **38%**,是所有
非植物命中里最高的,结构 core RMSD 1.125 Å,结构域架构 PF22600+PF03828 与三个已验证
的酶完全一致,却只被注释为 "Uncharacterized protein"。⚠️

### 4.7 古菌:一条独立的线索(local 模式,含 pLDDT)

*M. jannaschii* 全蛋白质组本地搜索,9 个蛋白三残基全保守且 core ≤ 2.5 Å,
催化壳层 pLDDT **92.3–96.9**:

| core RMSD | 登录号 | 催化残基(本流程推出) | 注释 |
|---|---|---|---|
| 1.197 | Q58612 | D47 / D49 / D79 | 推定蛋白腺苷酸转移酶 MJ1215(VII 型 TA 系统抗毒素) |
| 1.317 | Q58942 | D23 / D25 / D73 | 推定蛋白腺苷酸转移酶 MJ1547 |
| 1.399 | Q58614 | D43 / D45 / D75 | 推定蛋白腺苷酸转移酶 MJ1217 |
| 1.413 | Q58511 | D69 / D71 / D124 | **CCA 添加酶**(已知功能,阳性对照) |
| 1.562 | Q58021 | D41 / D43 / E75 | 推定蛋白腺苷酸转移酶 MJ0604 |
| 1.610 | Q57877 | D38 / D40 / D70 | 推定蛋白腺苷酸转移酶 MJ0435 |
| 1.624 | Q58486 | D157 / D159 / D229 | **Uncharacterized protein MJ1086** |
| 1.629 | Q57592 | D43 / D45 / D75 | 推定蛋白腺苷酸转移酶 MJ0128 |
| 2.183 | Q60312 | D26 / D28 / D62 | **Uncharacterized protein MJ0002** |

⚠️ 有意思的地方:这批古菌蛋白被注释为**蛋白**腺苷酸化酶(AMP 化修饰蛋白质)
和毒素-抗毒素系统组分,而不是核酸聚合酶,但它们保留了完整的三羧基构型。
它们能不能聚合 NTP,**从未被测过**。

### 4.8 真菌:酿酒酵母全蛋白质组(local 模式)

*S. cerevisiae* 全蛋白质组本地搜索,三残基全保守且 n_core ≥ 30 的只有 4 个 ——
而且全部是已知的核苷酸转移酶,没有漏网也没有噪声:

| core RMSD | 登录号 | 基因 | 催化残基(本流程推出) | 壳层 pLDDT | 蛋白 |
|---|---|---|---|---|---|
| 1.162 / 56 | P48561 | **TRF5** | D233 / D235 / D290 | 92.5 | poly(A) RNA polymerase 1(TRAMP 复合物) |
| 1.177 / 56 | P53632 | **PAP2/TRF4** | D236 / D238 / D293 | 90.7 | poly(A) RNA polymerase 2(TRAMP 复合物) |
| 1.250 / 51 | P29468 | **PAP1** | **D100 / D102 / D154** | 95.5 | 经典 poly(A) 聚合酶 |
| 1.655 / 41 | P21269 | CCA1 | D81 / D83 / D142 | 82.4 | 线粒体 CCA tRNA 核苷酸转移酶 |

两点值得注意:(1)这是一次**干净的全蛋白质组筛选** —— 6000 多个蛋白里只有 4 个
通过,且 4 个全是真正的核苷酸转移酶,说明过滤条件既不漏也不滥。(2)酵母 PAP1 的
催化三联体被推为 **D100/D102/D154**,与 OsNTP4 的 **D100/D102/D157** 几乎同位;
真正的共性信号是 `D-x-D` 加下游约 +52~+55 处的第三个 D 这一间距模式。✅

⚠️ TRF4/TRF5/PAP1 都是研究得极透的酶,但据我们所知**都没有在"不加引物、不加模板"
的条件下被测过**,因此它们和 4.5 节的 A 组候选属于同一类最省力的实验预测。

### 4.9 粟酒裂殖酵母:一次完整的家族召回(local 模式)

*S. pombe* 全蛋白质组(5196 个模型)本地搜索,三残基全保守、core ≤ 2.5 Å、
n_core ≥ 30 的共 9 个,**全部配准偏移为 0**,催化壳层 pLDDT 88–96:

| core RMSD | 登录号 | 基因 | 催化残基(本流程推出) | 壳层 pLDDT | 蛋白 |
|---|---|---|---|---|---|
| 1.115 / 55 | Q9UT49 | **cid13** | D110 / D112 / D164 | 94.5 | poly(A) RNA polymerase |
| 1.153 / 37 | Q9UTQ0 | cca1 | D56 / D58 / D114 | 88.2 | tRNA 核苷酸转移酶 |
| 1.165 / 55 | O13798 | **cid16** | D940 / D942 / D992 | 93.7 | caffeine-induced protein 16(1202 aa) |
| 1.170 / 55 | O74326 | **cid11** | D106 / D108 / D160 | 93.2 | poly(A) RNA polymerase |
| 1.186 / 38 | Q9Y7U9 | cca2 | D55 / D57 / D117 | 89.6 | tRNA 核苷酸转移酶 |
| 1.260 / 51 | Q10295 | **pla1** | D99 / D101 / D153 | 94.9 | 经典 poly(A) 聚合酶 |
| 1.263 / 51 | O13833 | **cid1** | D101 / D103 / D160 | 93.6 | 末端尿苷酰转移酶(该家族原型) |
| 1.304 / 57 | Q9UTN3 | **cid14** | D298 / D300 / D356 | 94.2 | poly(A) RNA polymerase |
| 1.393 / 52 | O74518 | **cid12** | D77 / D79 / D129 | 96.4 | poly(A) RNA polymerase |

**这是一次干净的全家族召回:cid1 / cid11 / cid12 / cid13 / cid14 / cid16 六个 Cid
家族成员一个不漏,外加 pla1 和两个 CCA 酶,没有任何可疑命中。** 5196 个蛋白里只有
这 9 个通过,而这 9 个恰好就是裂殖酵母已知的全部末端核苷酸转移酶。✅
这是对筛选条件灵敏度与特异性最直接的证明。

⚠️ Cid1 正是这一家族被研究得最透的原型(poly(U) polymerase),但据我们所知
同样没有在"不加引物、不加模板"条件下测过。

### 4.10 恶性疟原虫(local 模式)

*P. falciparum* 3D7(5168 个模型)筛出 5 个,配准偏移全为 0,壳层 pLDDT 90.9–96.5:

| core RMSD | 登录号 | 催化残基 | 壳层 pLDDT | 注释 |
|---|---|---|---|---|
| 1.252 / 56 | **A0A5K1K8X2** | D484 / D486 / D539 | 90.9 | **Nucleotidyltransferase, putative**(858 aa) |
| 1.273 / 56 | Q8I5E4 | D467 / D469 / D524 | 93.2 | poly(A) RNA polymerase mitochondrial-like(1172 aa) |
| 1.320 / 51 | C6KT92 | D175 / D177 / D229 | 96.4 | poly(A) polymerase |
| 1.334 / 46 | **Q8IJP0** | D336 / D338 / D398 | 96.5 | **Nucleotidyltransferase, putative**(615 aa) |
| 1.419 / 49 | **A0A144A2R0** | D184 / D186 / D241 | 93.6 | **Nucleotidyltransferase, putative**(473 aa) |

⚠️ 其中三个只被注释为 "Nucleotidyltransferase, putative",是囊泡虫类里最值得
跟进的候选;`A0A5K1K8X2` 的结构域架构(PF22600+PF03828)与三个已验证的植物酶一致。

### 4.11 一个说明 `n_core` 为什么必须报告的例子

搜索结果里排第 7–13 位的是假单胞菌的双功能 uridylyltransferase(GlnD),
core RMSD 低至 0.909 Å —— 看着比大多数真核命中还好。但它们的 **n_core 只有 21**,
而植物/动物命中是 50–57。也就是说那只是一个**小结构域层面**的局部吻合。
只看 RMSD 会被这类命中骗过去。✅

---

## 五、边界、局限与必须做的对照

### 这次做到了什么(✅)
- 完整复现了论文的发现链路,并**精确推出了 OsNTP4 的三个催化残基**(D100/D102/D157)。
- 证明该催化构型跨越后生动物、真菌、眼虫类、变形虫、囊泡虫、细菌、古菌。
- 四个物种做了全蛋白质组本地筛选(含 pLDDT 与侧链几何):裂殖酵母**完整召回了
  整个 Cid 家族 + pla1 + 两个 CCA 酶,5196 个蛋白里只有这 9 个通过**,不漏不滥,
  是筛选条件灵敏度与特异性的直接证据。
- 阈值是从已知阳性/阴性校准出来的,不是搬的。
- 跨门类排序把唯一已被独立证实具有该活性的非植物蛋白(TENT4B)排进前列。

### 这次没有做到什么(⚠️,必须明说)
- **没有任何活性证据。** 全部候选都是结构 + 残基保守性提名。
- **没有做系统发育分树。** 论文第 8 步(MEGA 建树 → 5 个 clade → 活性全在 Clade I)
  是先验强度最高的一步,本次只做到聚类前的排序,**clade 归属尚未计算**。这是
  最该补的一步。
- **没有做复合物建模。** 论文第 9 步(AlphaFold 预测结合 RNA/UTP)未复现。
- **截断偏倚。** 搜索有 1000 条上限,因此"某类群没出现"不构成阴性结论。
- **注释分类是关键词匹配**,很粗糙;"unknown" 只表示"注释没给出有用信息",
  不等于"新"。
- **web 模式拿不到 pLDDT**,只有 local 模式有;A/B 组候选的催化壳层置信度尚未逐个核查。

### 任何活性主张之前必须爬完的阶梯
(详见 `reference/validation.md`,全部照搬论文自己做的对照)

1. 完整蛋白质 MS —— 排除共价结合的核苷酸/辅因子(论文:37,681.0 vs 理论 37,681.53 Da)。
2. 蛋白酶 K 全消化 + 离子交换 —— 排除游离 NTP / 寡核苷酸残留,**必须带同样处理的标准品**。
3. 多步正交纯化 —— 亲和 → 离子交换 → 分子筛。
4. **洗过的晶体做活性测定** + 同一液滴里的非晶体对照(最关键的一条)。
5. 催化残基 → Ala 各自失活,**且表面残基突变不失活**(论文:D100A/D102A/D157A 全废,K194A 基本无影响)。

产物化学本身还需要独立证据:RNase 敏感性**不能**确定键型,论文用 ¹H-³¹P HMBC
看桥连磷酸与 H3′/H5′ 的相关来确证 3′–5′ 连接(2′–5′ 会给出 H2′ 相关)。

---

## 复现方式

```bash
# 依赖
pip install numpy scipy biopython
curl -sSL -o foldseek.tar.gz https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xzf foldseek.tar.gz

# 1) 复现论文:cGAS 作查询
python3 .claude/skills/catalytic-scaffold-mining/scripts/mine.py \
  --query 4km5.pdb --chain A --catalytic 225,227,319 \
  --mode web --databases afdb-proteome afdb-swissprot --out out/cgas

# 2) 跨门类扩展:AtURT1 作查询
python3 .claude/skills/catalytic-scaffold-mining/scripts/mine.py \
  --query 6L8K.pdb --chain A --catalytic 491,493,547 --window 2 \
  --mode web --databases afdb-proteome afdb-swissprot --out out/urt1

# 3) 加谱系与已知酶类注释
python3 .claude/skills/catalytic-scaffold-mining/scripts/annotate.py \
  out/urt1.candidates.tsv --out out/urt1.annotated.tsv

# 4) 用已知阳性/阴性校准阈值
python3 .claude/skills/catalytic-scaffold-mining/scripts/calibrate.py \
  out/urt1.annotated.tsv \
  --positives O64642 A0A0P0VE25 A0A0P0XSC2 \
  --negatives Q9HI95 Q46AZ4 Q8TGX2

# 5) 单个蛋白质组本地搜索(有 pLDDT 与侧链几何)
python3 .claude/skills/catalytic-scaffold-mining/scripts/mine.py \
  --query 6L8K.pdb --chain A --catalytic 491,493,547 --window 2 \
  --mode local --target-dir extracted/METJA --out out/metja
```
