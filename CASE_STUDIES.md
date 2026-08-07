# 用真实临床数值检验 Recorder 候选：数据空间、可行性与速度

针对《Disease Recorder Biomarker：系统化重构与新检测内容挖掘》里的 3+2 组合，
把每个候选的**已发表正常值与病理值**取回来，算出效应量、可达 AUC、分析误差预算和
Gate-2 样本量，再回答"哪个能最快做出来"。

配套：`recorder_pipeline/r5_datacases.py`（全部数字由此计算）、
`RECORDER_FRAMEWORK.md`（三条定律）、`figures/recorder_r5_dataspace.png`。

标签：**[OK]** 计算所得 · **[LIT]** 已发表数值 · **[EST]** 只有中位数发表、离散度为假设 ·
**[HYP]** 无直接人体数据。

---

## 0. 三个可以立刻用的结论

**① 最快能做的方案是 SLE 的载体年龄校正，而它不需要任何新 assay。**
EC4d/BC4d 已是商业化流式检测；网织红细胞比例（IRF）与未成熟血小板比例（IPF）是
Sysmex 等血球仪的**常规输出，每一张 CBC 上都有** [LIT]。所以"按载体年龄校正的 C4d"
可以完全作为**算法**建立在既有数据上——零新试剂、零新采样、可回溯既往样本库。而且已有
先例：红细胞结合 C4d / CR1 比值已被用于区分 SLE flare 与感染 [LIT]。
报告把它排在第 4（74 分，新颖性验证篮）。按"速度×确定性"，它应该排第 1。

**② 报告的优先级与数据空间不一致，最严重的两处：**

- **CKD 组合的价值几乎全在 uromodulin，不在 C-Alb。** 血清 uromodulin 从非 CKD 的
  228 ng/mL 跌到 G3 的 52、G5 的 13 [LIT]——**4.4×–17.5×**。C-Alb 在 ESRD vs 非尿毒症
  只有 **2.1×**；而在它真正的 intended use（ESRD 内部的一年死亡预测）只有
  0.77% → 1.01%，即 **1.31×**，算出 AUC ≈ 0.65 [OK]，**够不到 0.70 的临床有用线**。
  报告写的是"C-Alb 为主 + uromodulin 做二维模型"，按数据应当反过来。
- **PRO-C3 是全表效应量最小、assay CV 最大的候选，却拿了最高的非基准分（88）。**
  健康参考区间 6.1–14.7 ng/mL，而 **F3 的中位数 14.5 落在这个区间之内** [LIT]。
  F4/F0-1 = 16.3/9.5 = **1.7×**。同时 PRO-C3 ELISA 的批间 CV 是 **11.03%** [LIT]，
  是全表最高。

**③ 用健康参考区间做样本量估算会系统性 underpower 约 2 倍。**
PRO-C3 由健康参考区间推出的 AUC 是 0.808，而临床队列实测 0.79 [LIT]。反解得真实人群的
对数尺度 SD 是健康参考区间的 **2.1 倍** [OK]。原因是临床对照不是健康献血者，是 F0–F2 的
NAFLD 患者——他们的 PRO-C3 已经抬高了。所有 Gate 2 设计都必须用**疾病对照**的离散度。

---

## 1. 数据空间总表（全部为已发表数值）

| 候选 | 正常 | 病理 | 倍数 | assay CV | 来源 |
|---|---|---|---|---|---|
| **GFAP** (pg/mL) | 健康中位 **8.0**（IQR 3–14） | CT−/MRI+ **414.4**（IQR 139–813）；CT+ 均值 3970 (SD 7820) | **51.8×** | ~10% | [LIT] |
| **TAT** (µg/L) | HISCL 参考 **<4** | DIC 进展阈值 **≥10.8** | **5.4×** | **3.67%** | [LIT] |
| **血清 uromodulin** (ng/mL) | 非 CKD **228**（献血者 207） | G1 153 · G2 107 · **G3 52** · G4 30 · G5 13 | **4.4×**（G3）**17.5×**（G5） | ~8% | [LIT] |
| **C-Alb** (%Lys549) | 非尿毒症 **0.42** | ESRD **0.90** | **2.1×** | ~8% | [LIT] |
| **C-Alb（真正用途）** | 存活>1y **0.77** | 死亡<1y **1.01** | **1.31×** | ~8% | [LIT] |
| **PRO-C3** (ng/mL) | 健康 **6.1–14.7**；F0/F1 9.5 | F3 **14.5** · F4 **16.3** | **1.7×** | **11.03%**（批间） | [LIT] |
| **EC4d** (net MFI) | — | >8.9：70%/93%(健康)；>14：45%/95%(他病) | — | 流式 | [LIT] |
| **BC4d** (net MFI) | — | >48：66%/96%；>60：54%/95% | — | 流式 | [LIT] |

计算结果（`r5_datacases.py`）：

| 候选 | Δln | 可达 AUC | 分析噪声/生物噪声 | Gate-2 n/组（AUC>0.70） |
|---|---|---|---|---|
| GFAP | 3.95 | 0.988 | 0.08 | **5** |
| TAT | 1.69 | 0.960 | 0.05 | **5** |
| uromodulin | 1.48 | 0.895 [EST] | 0.10 | 10 |
| C-Alb（vs 正常） | 0.76 | 0.850 [EST] | 0.15 | 23 |
| PRO-C3 | 0.54 | 0.801 [EST] | **0.25** | **61** |
| C-Alb（死亡预测） | 0.27 | 0.649 [EST] | 0.16 | **达不到** |

> 离散度为 [EST] 的行做了 0.7×–1.5× 敏感性扫描：**排序在整个扫描区间内不变**，
> 绝对 AUC 最多移动 0.10。把这些数字当作排序，不要当作规格。

---

## 2. 六个 Case Study

### Case 1 — SLE：载体年龄校正的 C4d ★ 最快可实现

**现状数值** EC4d >8.9 net MFI → 70% 敏感 / 93% 特异（对健康）；同一指标提高到
>14 → 45% / 95%（对其他自身免疫病）[LIT]。换算到同一尺度 [OK]：

| 操作点 | d′ | AUC |
|---|---|---|
| EC4d >8.9（vs 健康） | 2.00 | 0.921 |
| EC4d >14（vs 他病） | 1.52 | 0.859 |
| BC4d >48（vs 健康） | 2.16 | 0.937 |
| BC4d >60（vs 他病） | 1.75 | 0.891 |

**第一个必须纠正的读法**：同一个标志物对健康人比对其他自身免疫病"好 0.06 AUC"。
**只有后者是 intended use**。报告在这一点上是对的（强调 intended use），但候选评分表里
没有区分这两个数字。

**报告提出的增量**：`cell-bound C4d / carrier count` + 载体年龄代理（IRF、IPF）。

**为什么这是最快的**：
- EC4d/BC4d：**已商业化**（流式，AVISE）
- IRF 参考区间 **1.6–12.1%**，IPF **0.8–5.6%**（女 0.8–4.7 / 男 0.7–6.1）[LIT]——
  **这两个数已经在每一份常规血常规报告里**
- 已有同类先例：EC4d/CR1 比值区分 flare 与感染 [LIT]

→ **可以在既有样本库/既有 LIS 数据上做回顾性验证，不需要新试剂或新采血。**
按 `RECORDER_FRAMEWORK.md` §3，这也正是"每单位载体的印记量"这条强制要求的实现。

**最小证伪实验**：取既有 EC4d + CBC 配对数据，比较三个模型对 flare 的判别：
① EC4d 原始值 ② EC4d/RBC 计数 ③ EC4d ~ f(RBC 计数, IRF, 年龄) 的条件残差。
**Go 条件**：③ 在校正后仍保留独立信息且优于 ①（报告 Gate 2 的停止规则）。

**风险**：溶血、输血、近期失血都会同时改变载体年龄与 C4d，必须作为混杂前置富集。

---

### Case 2 — 脓毒症：TAT/PIC/tPAIC ★ 平台已现成

**现状数值** HISCL 厂商参考区间：TAT **<4 µg/L**、TM 3.8–13.3 TU/mL、
tPAIC **<10.5 µg/L**、PIC **<0.8 mg/L**；ISTH-DIC 评分 <5 的患者中，
TAT **≥10.8 ng/mL** 提示 7 天内进展为不可逆 DIC [LIT]。
批内 CV：TAT **3.67%**、PIC 6.51%、TM 3.64%、tPAIC **2.46%** [LIT]。

**计算** [OK]：Δln = 1.69（5.4×），可达 AUC ≈ 0.96，分析噪声只占生物噪声的 5%，
Gate-2 只需 **每组 5 例**就能证明 AUC>0.70。

**判断**：这是全表**分析性能最好**的一组——CV 低一个数量级，参考区间已建立，仪器已上市。
按 `RECORDER_FRAMEWORK.md` 的带误差判据，TAT/PIC 配对所需的 ρ 只从 0.50 抬到 **0.51**，
即**分析噪声对这一组几乎没有影响**。

**但要注意时间窗**（这是我与报告的一处分歧）：TAT 循环半衰期约 15 分钟，按定律二，
它记录的是**过去约一小时**的凝血酶生成，不是"累积凝血负荷"。所以它适合报告说的
"早期凝血表型分型"，**不适合**做累积负担或疗效积分。报告的 intended use 写对了，
但把它归进"recorder"类时应当标注这个时间尺度。

---

### Case 3 — CKD：应该把 uromodulin 提到主位

**现状数值** [LIT]：

```
血清 uromodulin：非CKD 228 → G1 153 → G2 107 → G3 52 → G4 30 → G5 13 ng/mL
                 预后阈值 55.6 ng/mL（高于此值 ESKD/心血管/死亡风险下降）
C-Alb（%Lys549）：非尿毒症 0.42% → ESRD 0.90%
                 ESRD 内部：一年内死亡 1.01% vs 存活 0.77%（Q4 vs Q1 校正 HR 1.90）
```

**计算** [OK]：uromodulin Δln=1.48（AUC≈0.90）；C-Alb 对正常 Δln=0.76（AUC≈0.85）；
**C-Alb 对死亡 Δln=0.27（AUC≈0.65）**。

**这就是关键问题**：C-Alb 的"疾病 vs 正常"效应尚可，但它在**自己的 intended use**
（ESRD 内部风险分层）上只有 1.31 倍，AUC 0.65。按 Hanley-McNeil，**任何样本量都无法
证明它超过 0.70** [OK]。敏感性扫描的乐观端刚好触到 0.70，所以结论是
**"重新定义 intended use"，不是"直接杀掉"**——例如把它用于 CKD 3–4 期的进展预测
（那里跨度更大），而不是 ESRD 内部的死亡预测。

**二维模型的正确形式**：报告说"累积尿毒负担 × 肾实质储备"，并正确指出不应相乘。
补充一点：两个轴的**效应量差 5 倍**，所以联合模型里 uromodulin 会主导，C-Alb 的增量
必须单独证明（即：`AUC(uromodulin + C-Alb) > AUC(uromodulin)`，而不是
`> AUC(C-Alb)`）。这个比较应当预注册。

**可行性**：uromodulin ELISA 已商品化，效应量大，是 CKD 线**先做的那一半**。
C-Alb 需要 MS（或新抗体），属于第二步。

---

### Case 4 — 纤维化：PRO-C3 的两个硬问题

**现状数值** [LIT]：健康参考 **6.1–14.7 ng/mL**；NAFLD F0/F1 **9.5**、F3 **14.5**、
F4 **16.3**；PRO-C3 单独识别显著/进展期纤维化 AUC **0.83 / 0.79**；ADAPT（含年龄、
糖尿病、血小板）AUC **0.86–0.87**。

**问题一：效应量小到与健康区间重叠。** F3 中位数 14.5 **落在健康参考区间 6.1–14.7 内**。
这不是可以靠 cutoff 调整解决的——它意味着 F3 的一半人在健康区间里。ADAPT 之所以有效，
**是因为它加了年龄/糖尿病/血小板，不是因为 PRO-C3 本身强**。

**问题二：跨研究标定不一致。** 同为"NAFLD 进展期纤维化"，一项研究中位数 16.3，
另一项报告成人 **31.7**（none/mild 17.9）[LIT]——**接近 2 倍的系统差**。叠加
**11.03% 的批间 CV**，意味着**cutoff 不能跨实验室迁移**。报告把关键缺口写为
"器官与病因特异阈值"（生物学问题），但更紧迫的是**标定与精密度问题**。

**带测量误差的比值判据** [OK]：PRO-C3 的健康对数 SD ≈ 0.2245，而批间 CV 11% →
分析噪声约为生物离散度的 **0.49 倍**。代入

```
ρ_bio > (κ² + a²) / (2κ)
```

得：PRO-C3/C3M 比值要跑赢单指标，需要两者的生物学相关 **ρ > 0.62**，而不是教科书的
0.50。也就是说，**分母的不精密度把门槛抬高了 24%**。这个 ρ 必须先在健康/稳定队列里
实测，再决定做不做比值——不能假设。

**建议**：这条线仍值得做（它是平台校准项目），但要把
"**先把 CV 从 11% 压到 5%**"写成 Gate 1 的前置条件，并把 ADAPT 类协变量模型作为
对照臂——否则很难证明 recorder 形式本身的增量。

---

### Case 5 — TBI：GFAP 的数据空间是全表最好的

**现状数值** [LIT]：

```
健康对照        中位 8.0 pg/mL（IQR 3.0–14.0）
CT−/MRI−        中位 74.0 （17.5–214.4）
CT−/MRI+        中位 414.4（139.3–813.4）
CT+ （全谱）    均值 3970.1（SD 7819.6）；CT− 均值 363.8（SD 706.3）
Abbott i-STAT   GFAP cutoff 30 pg/mL（已 FDA 上市）
```

**计算** [OK]：健康 → CT−/MRI+ 为 **51.8×**（Δln 3.95），可达 AUC 0.988；
分析噪声只占生物离散度的 8%；Gate-2 每组 **5 例**。

**这解释了一件事**：报告给 GFAP-BDP 78 分，低于 PRO-C3 的 88 分。但从可测量性看，
GFAP 的效应量比 PRO-C3 **大 30 倍**，平台已 FDA 清关，母体丰度极高（分数形式的分母
不缺）。评分低是因为"证据 C 级 + 端点未确认"——那是**发现风险**，不是**测量风险**。
这两种风险不应该压进同一个分数。

**片段假设的证据状态** [LIT]：2025 年的工作报告 20–26 kDa 产物超过 15–19 kDa 产物，
与 coil1-BDP 优先释放/coil2 保留一致，并且**GFAP 片段谱（而非未切割 GFAP）预测预后**。
需要标注：该结果目前是 **bioRxiv 预印本**，尚未见同行评审定稿——报告给的 C 级是恰当的，
但"预印本"这一点应当明写。钙蛋白酶切位点已定位（N 端 A-56*A-61、C 端 T-383*Q-388，
限制性片段 38K）[LIT]，这为 neo-terminus 抗体提供了明确靶点。

**可行性**：需要新抗体 → 慢。但**先做免疫富集-PRM 确认端点**（报告的路线）是对的，
且因为效应量大，验证队列可以很小。

---

### Case 6 — 动脉粥样硬化：抗体已经存在，但化学计量是拦路虎

**现状** [LIT]：Tyr192 是 MPO 对 ApoA-I 氯化/硝化的**主位点**（脂质游离 ApoA-I 中
3-氯酪氨酸产率约 80% 集中于此）；**只有氯化**显著损害 ABCA1 依赖的胆固醇外排；
CAD 与 ACS 患者 HDL-ApoA-I 的 3-氯酪氨酸含量高于健康对照；氧化 ApoA-I 最高三分位者
CVD 风险为最低三分位的 **6 倍（硝化）与 16 倍（氯化）**；病变中约
**5 个 HDL 分子里有 1 个**带 MPO 修饰印记。

**两条对报告的重要补充**：

1. **抗 chloro-192-Tyr ApoA-I 单克隆抗体已经被开发出来**（2019 年起有多篇）[LIT]。
   报告写"用 PRM 做分析验证，**再决定是否开发 neoepitope 抗体**"——实际上这一步已有
   现成起点，可以直接评估既有抗体，省掉一轮抗体工程。
2. **真正的拦路虎不是灵敏度，是特异性**。血浆 ApoA-I 浓度约 1–1.5 mg/mL，属于高丰度
   蛋白；即使修饰比例只有 ~0.1%，绝对量仍有 ~1 µg/mL，**够测**。问题是要在
   99.9% 未修饰母体的背景里认出 0.1% 的修饰型，**要求抗体的判别比 >1000:1**。
   这才是 Gate 1 的核心指标，而报告的质控表里列的是"人工氧化对照"——需要再加一条
   **未修饰母体的交叉反应上限**。
   注意"1/5 HDL 带印记"是**病变内**的数字，不能外推到循环血浆 [OK 区分]。

---

## 3. 按"多快能做出来"排序

| 排名 | 项目 | 新试剂 | 新采样 | 可回顾既有样本 | 预计到 Gate 2 |
|---|---|---|---|---|---|
| **1** | **SLE C4d + 载体年龄校正** | **无** | **无** | **是（含 LIS 数据）** | **最快：纯算法** |
| 2 | 脓毒症 TAT/PIC/tPAIC | 无（HISCL 已上市） | 需前瞻队列 | 部分 | 快 |
| 3 | CKD uromodulin（先做） | 无（ELISA 已商品化） | 否 | 是 | 快 |
| 4 | 纤维化 PRO-C3/C3M | 无，但需先降 CV | 否 | 是 | 中（受标定拖累） |
| 5 | CKD C-Alb | MS 或新抗体 | 否 | 是 | 中 |
| 6 | TBI GFAP-BDP | **需新抗体** | 否 | 是（CSF/血浆库） | 慢，但验证队列小 |
| 7 | ApoA-I Tyr192 氧化 | 已有候选抗体，需 PRM 定标 | 否 | 是 | 慢 |

**"能快速实现的不错方案"的直接回答**：有，就是第 1 项。它是唯一一个
**不需要新试剂、不需要新采样、可以在既有数据上当天开始**的方案，而且它命中的正是
`RECORDER_FRAMEWORK.md` 里"载体混杂是记录器通用软肋"这条强制要求。第 2、3 项紧随其后，
都是"仪器已上市、只差临床设计"。

---

## 4. 与报告评分不一致之处（供讨论）

| 项目 | 报告分 | 数据空间给出的意见 |
|---|---|---|
| 纤维化 PRO-C3+C3M | 88（最高非基准） | 效应量最小（1.7×）、CV 最大（11%）、跨研究标定差 2×。仍值得做为平台校准，但不应领跑 |
| CKD C-Alb + uromodulin | 84 | 价值主要在 uromodulin（4.4–17.5× vs 2.1×）；C-Alb 在其 intended use 上 AUC≈0.65，够不到线 |
| TBI GFAP-BDP | 78 | 测量学上最优（52×、平台已清关）；78 分反映的是发现风险，不是测量风险——两者应拆开评 |
| SLE platelet-C4d/载体年龄 | 74 | 实现速度最快（零新 assay）；按"速度×确定性"应上调 |

**根因**：报告的七维权重里有"检测可行性 15%"，但那是**定性判断**，没有把效应量与
assay CV 放进去。建议加一维：**可测量性 = f(Δln, CV, 母体丰度)**，可直接由
`r5_datacases.py` 计算，不需要人工打分。

---

## 5. 对 90 天行动的修订建议

报告的 3+2：`PRO-C3/C3M/CTX-III` + `EC4d/BC4d` + `C-Alb+uromodulin`，
新颖项 `GFAP-BDP` + `Tyr192-oxApoA-I`。按数据空间建议调整为：

1. **第 0–2 周（零成本）**：SLE 既有 EC4d/BC4d + CBC 数据的回顾性载体年龄校正分析。
   这一项在 90 天内就能出 Gate-2 级别的结论。
2. **第 0–4 周**：在健康/稳定队列上实测 PRO-C3 与 C3M 的 **ρ 和 log-log 斜率**。
   这是决定"做不做比值"的前置数据，成本极低，而且直接执行报告自己的
   `log(A)=α+βlog(B)` 判据。
3. **调整 CKD 线的顺序**：先 uromodulin（大效应、已商品化），C-Alb 作为增量证明。
4. **PRO-C3 线增加 Gate 1 前置条件**：批间 CV 从 11% 降到 ≤5%，并建立跨批标定品。
5. **ApoA-I 线**：先评估已有的 anti-chloro-192-Tyr 抗体，把"未修饰母体交叉反应
   ≤0.1%"写成 Go/No-Go，再考虑自研。
6. **GFAP 线**：保持免疫富集-PRM 优先；因效应量大，把验证队列规模从默认值下调。

---

## 6. 局限

1. **[EST] 行的离散度是假设的。** 只有 GFAP 的正常与病理离散度都已发表。其余候选只有
   中位数，IQR 为假设。已做 0.7×–1.5× 敏感性扫描，**排序稳健，绝对 AUC 不稳健**。
2. **AUC 用等方差二项正态近似。** 真实分布右偏（GFAP 的 SD 大于均值），实际 AUC 会与
   计算值有偏差；本文用它做排序与量级判断，不作性能承诺。
3. **Gate-2 样本量是统计下限。** 真实研究还要覆盖各案例点名的混杂（eGFR、BMI、溶血、
   血小板计数、感染），这通常才是决定 n 的因素。
4. **跨研究数值合并有风险。** PRO-C3 的 2 倍标定差就是例证；不同平台的 GFAP 值同样不可
   直接比较（ARCHITECT 中位 143.3 vs i-STAT 116.0 pg/mL）[LIT]。
5. **本文只评测量学，不评临床价值。** 一个 AUC 0.99 的标志物如果不改变决策，仍然没有
   产品价值——这一点报告的 Gate 3（decision curve）是对的。

---

## 参考来源

- PRO-C3 健康参考区间与分期中位数 — [Determining a healthy reference range for PRO-C3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8385245/)、[ADAPT 算法](https://journals.lww.com/hep/fulltext/2019/03000/adapt__an_algorithm_incorporating_pro_c3.19.aspx)、[PRO-C3 in NAFLD](https://www.sciencedirect.com/science/article/pii/S2589555919300618)
- PRO-C3 ELISA 精密度（批内 4.11% / 批间 11.03%）— [neo-epitope specific PRO-C3 ELISA](https://pmc.ncbi.nlm.nih.gov/articles/PMC3633973/)
- EC4d/BC4d 操作点与 AVISE 性能 — [ACR abstract: CB-CAPs vs C3/C4](https://acrabstracts.org/abstract/cell-bound-complement-activation-products-have-higher-sensitivity-than-serum-c3-and-c4-levels-in-systemic-lupus-erythematosus/)
- 红细胞 C4d/CR1 区分 flare 与感染 — [J Immunol Res](https://pmc.ncbi.nlm.nih.gov/articles/PMC4529962/)
- IRF / IPF 参考区间 — [Determination of reference ranges for IPF and IRF](https://pubmed.ncbi.nlm.nih.gov/27863758/)
- 血清 uromodulin 分期数值与预后阈值 — [Serum uromodulin and CKD progression](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-018-1693-2)、[uromodulin as marker of parenchymal integrity](https://pubmed.ncbi.nlm.nih.gov/28206617/)
- 碳酰化白蛋白与 ESRD 死亡 — [Sci Transl Med, Berg et al.](https://www.science.org/doi/10.1126/scitranslmed.3005218)、[homocitrulline vs C-Alb 比较](https://pmc.ncbi.nlm.nih.gov/articles/PMC10689527/)
- TAT/PIC/TM/tPAIC 参考区间与 CV — [HISCL 高敏分析仪性能验证](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12704508/)、[TAT 与 DIC 进展](https://ashpublications.org/blood/article/142/Supplement%201/1268/503156/The-Proactive-Diagnostic-Value-of-Thrombin)
- GFAP 健康对照与 TRACK-TBI 分层数值 — [Lancet Neurol CT-negative TBI](https://pubmed.ncbi.nlm.nih.gov/31451409/)、[i-STAT vs ARCHITECT 比较](https://pmc.ncbi.nlm.nih.gov/articles/PMC8086519/)、[i-STAT TBI 综述](https://www.tandfonline.com/doi/full/10.1080/14737159.2024.2306876)
- GFAP 钙蛋白酶切位点与片段谱预后（**预印本**）— [Calpain/caspase-6 GFAP BDP](https://pmc.ncbi.nlm.nih.gov/articles/PMC9409281/)、[GFAP Degradation in TBI, bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.08.01.668181v5)
- ApoA-I Tyr192 氯化与抗体 — [JBC: Tyr192 是主位点](https://www.jbc.org/article/S0021-9258(19)63106-6/fulltext)、[anti-chloro-192-Tyr 抗体开发](https://pubmed.ncbi.nlm.nih.gov/31386835/)、[MPO 与 apoA-I 位点特异氧化](https://pubmed.ncbi.nlm.nih.gov/22219194/)
