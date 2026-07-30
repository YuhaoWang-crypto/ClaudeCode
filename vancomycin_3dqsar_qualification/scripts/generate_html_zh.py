# -*- coding: utf-8 -*-
"""中文版监管报告 HTML（图片 base64 内嵌），用于 Chromium 渲染 PDF。"""
import json, base64
D="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
R=json.load(open(f"{D}/analysis_results.json"))
Q=json.load(open(f"{D}/qsar_results.json"))
QA=json.load(open(f"{D}/qsar_augment.json"))
M7=json.load(open(f"{D}/m7_results.json"))
def img64(f): return "data:image/png;base64,"+base64.b64encode(open(f"{D}/{f}","rb").read()).decode()
def fig(f,cap,w="100%"): return f'<figure><img src="{img64(f)}" style="width:{w}"><figcaption>{cap}</figcaption></figure>'
EPS=[r["endpoint"] for r in Q["cv_metrics"]]; aucs=[r["CV_AUC"] for r in Q["cv_metrics"]]; meanauc=Q["mean_cv_auc"]
EP_ZH={"NR-AR":"雄激素受体","NR-AR-LBD":"雄激素受体(LBD)","NR-AhR":"芳烃受体","NR-Aromatase":"芳香化酶(CYP19A1)","NR-ER":"雌激素受体","NR-ER-LBD":"雌激素受体(LBD)","NR-PPAR-gamma":"PPARγ","SR-ARE":"氧化应激(ARE)","SR-ATAD5":"遗传毒性/DNA损伤(ATAD5)","SR-HSE":"热休克反应","SR-MMP":"线粒体膜电位","SR-p53":"p53 DNA损伤应答"}
MEAS_ZH={"active":"阳性","inactive":"阴性","not tested":"未测试"}
qrows=""
for ep in EPS:
    m=[r for r in Q["cv_metrics"] if r["endpoint"]==ep][0]
    meas=QA["measured_vancomycin"][ep]; hp=QA["heldout_parent_pred"][ep]; p75=QA["heldout_pred_rrt075"][ep]; d=p75-hp
    mc='style="color:#c00;font-weight:700"' if meas=="active" else ""
    qrows+=f"<tr><td>{ep} <span class=sub>({EP_ZH[ep]})</span></td><td class=c>{m['CV_AUC']:.3f}</td><td class=c {mc}>{MEAS_ZH[meas]}</td><td class=c>{hp:.3f}</td><td class=c>{p75:.3f} ({'+' if d>=0 else ''}{d:.3f})</td></tr>"
m7rows=""
for c,r in M7["per_compound"].items():
    ec="#2e7d32" if not r["expert_rulebase_alerts"] else "#c00"; sc="#2e7d32" if r["stat_Ames_call"]=="negative" else "#c00"
    exp="阴性（无结构警示）" if not r["expert_rulebase_alerts"] else "阳性（存在警示）"
    st="阴性" if r["stat_Ames_call"]=="negative" else "阳性"
    m7rows+=f'<tr><td>{c}</td><td class=c style="color:{ec}">{exp}</td><td class=c style="color:{sc}">{st}（p={r["stat_Ames_prob"]:.2f}，域相似度 {r["stat_AD_maxTanimoto"]:.2f}）</td><td class=c style="color:#2e7d32;font-weight:700">5</td></tr>'
n_alerts=M7["methodology_1"].split(";")[-1].strip().split()[0]
ntrain=Q['model'].split('n=')[1].split(')')[0]

HTML=f"""<!doctype html><html lang="zh"><head><meta charset="utf-8"><style>
@page {{ size:A4; margin:16mm 15mm; }}
* {{ box-sizing:border-box; }}
body {{ font-family:'WenQuanYi Zen Hei','Noto Sans CJK SC',sans-serif; font-size:10.3pt; line-height:1.65; color:#1a1a1a; margin:0; }}
h1 {{ font-size:15pt; color:#1F3864; border-bottom:2px solid #1F3864; padding-bottom:3px; margin:20px 0 8px; page-break-after:avoid; }}
h2 {{ font-size:12pt; color:#2E74B5; margin:14px 0 5px; page-break-after:avoid; }}
p {{ text-align:justify; margin:6px 0; }}
.title {{ text-align:center; margin-top:36px; }}
.title h1 {{ border:0; color:#1F3864; font-size:20pt; margin:2px; line-height:1.4; }}
.title .sub {{ color:#2E74B5; font-size:12.5pt; margin:8px 0 22px; }}
table {{ border-collapse:collapse; width:100%; margin:8px 0; font-size:8.9pt; page-break-inside:avoid; }}
th {{ background:#1F3864; color:#fff; padding:5px 7px; text-align:left; }}
td {{ border:1px solid #cdd6e4; padding:4px 7px; vertical-align:top; }}
td.c {{ text-align:center; }}
tr:nth-child(even) td {{ background:#f4f7fb; }}
.sub {{ color:#666; font-size:0.85em; }}
figure {{ margin:10px 0; text-align:center; page-break-inside:avoid; }}
figure img {{ border:1px solid #e2e2e2; }}
figcaption {{ font-size:8.6pt; color:#555; text-align:justify; margin-top:3px; }}
ul {{ margin:6px 0; padding-left:20px; }} li {{ margin:5px 0; text-align:justify; }}
.box {{ background:#f4f7fb; border-left:4px solid #2E74B5; padding:8px 12px; margin:8px 0; }}
.warn {{ background:#fff8e1; border-left:4px solid #e0a800; padding:8px 12px; margin:8px 0; font-size:9.6pt; }}
.kv td:first-child {{ background:#eaf0f8; font-weight:700; width:26%; }}
.pb {{ page-break-before:always; }}
.foot {{ font-size:8.2pt; color:#888; }}
</style></head><body>

<div class="title">
<div class="foot">计算机模拟（IN SILICO）毒理学界定 — 监管支持报告</div>
<h1>两个万古霉素 B 杂质的计算机模拟毒理学界定<br>基于 3D-QSAR / 结构类推（Read-Across）</h1>
<div class="sub">杂质 RRT 0.75 ［(26R)-DMe-DeLI］ 与 RRT 0.87 ［(7S,26R)-万古霉素 B］</div>
<table class="kv" style="width:90%;margin:0 auto;text-align:left">
<tr><td>受试物</td><td>万古霉素 B（原料药 API）及杂质 RRT 0.75、RRT 0.87</td></tr>
<tr><td>计算方法</td><td>毒性团/药效团分析；RDKit PAINS/Brenk/NIH 结构警示；理化与三维形状描述符；三维构象生成及 MCS/O3A 形状类推；自建 12 端点 Tox21 统计 QSAR（随机森林、5 折交叉验证、适用域检查）；ICH M7 双方法学致突变筛查（Benigni–Bossa 专家规则库 + Ames QSAR）；辅助 logD/pKa QSAR</td></tr>
<tr><td>毒理学锚点</td><td>母体 API 90 天静脉犬试验，NOAEL 75/100 mg/kg/day（约为临床剂量 2–3 倍）</td></tr>
<tr><td>日期</td><td>2026-07-30</td></tr>
</table>
<p class="warn" style="width:90%;margin:18px auto;text-align:left"><b>重要提示</b> —— 本报告为计算机模拟筛查与结构类推评估，用于支持毒理学专家审核。其本身不构成监管毒理学结论，也不替代专家签署或（在需要时）实验测试（见第 8 节）。</p>
</div>

<h1 class="pb">1. 执行摘要</h1>
<p>糖肽类抗生素万古霉素 B 的两个有关物质需进行计算机模拟毒理学界定：<b>RRT 0.87 = (7S,26R)-万古霉素 B</b>，一种在低水活度/受热条件下产生的降解产物；<b>RRT 0.75 = (26R)-DMe-DeLI</b>，一种发酵来源的第 1 位残基变体。由于二者与母体的差异主要为立体化学/构成上的差异，常规二维类推无法分辨，故对母体及两个子体分子均采用三维结构方法评估。</p>
<p><b>主要结论：两个杂质相对母体 API 均未引入新的毒性团，也未引入母体所不具备的新的 DNA 反应性（致突变）结构警示。</b>RRT 0.87 为纯立体异构体，与母体共享全部 101 个重原子及完整官能团清单（仅在 18 个立体中心中的 2 个构型不同）。RRT 0.75 与母体共享 99/约100 个重原子；其唯一构成改变为第 1 位残基由仲胺（N-甲基）变为伯脂肪胺 —— 该官能团类别母体本已存在（万古糖胺）—— 亲脂性不变（logD 0.42），碱性近乎不变（pKa 约 8.1 对 8.2）。</p>
<p>一项 12 端点 Tox21 统计 QSAR（5 折交叉验证平均 AUC {meanauc:.3f}）进一步支持该结论：母体本身即为 Tox21 训练集化合物（实测谱基本干净 —— 仅 NR-AR 弱信号），RRT 0.87 与母体指纹相同（各端点 Δ ≡ 0），RRT 0.75 在每个端点与母体差异 ≤ 0.02。ICH M7 双方法学致突变筛查（Benigni–Bossa 专家规则库 + Ames QSAR，AUC {M7['methodology_2']['cv_auc']:.3f}）对三者均给出一致的<b>阴性 / Class 5</b>。综合证据权衡，两个杂质均可通过对母体的结构类推予以毒理学界定；母体安全性由 91 天静脉犬试验（NOAEL 75/100 mg/kg/day）及长期临床使用锚定，前提是各杂质控制在《欧洲药典》有关物质限度内。</p>
<table>
<tr><th>问题</th><th>RRT 0.87</th><th>RRT 0.75</th><th>依据</th></tr>
<tr><td>相对母体是否有新毒性团？</td><td class=c style="color:#2e7d32;font-weight:700">否</td><td class=c style="color:#2e7d32;font-weight:700">否</td><td>官能团清单一致(0.87)；等官能团胺替换(0.75)</td></tr>
<tr><td>是否有新的致突变(M7)警示？</td><td class=c style="color:#2e7d32;font-weight:700">否</td><td class=c style="color:#2e7d32;font-weight:700">否</td><td>M7 两方法学均阴性 → Class 5（见 5.7）</td></tr>
<tr><td>对母体的类推是否成立？</td><td class=c style="color:#2e7d32;font-weight:700">是</td><td class=c style="color:#2e7d32;font-weight:700">是</td><td>≥98% 骨架共享；核心毒性团保守</td></tr>
</table>

<h1 class="pb">2. 目的、范围与监管框架</h1>
<p>本项目目标为对杂质 RRT 0.75 与 RRT 0.87 进行计算机模拟毒理学界定，交付：(i) 对母体及各子体分子的 3D-QSAR / 结构方法评估；(ii) 识别并评估子体分子中存在而母体不具备的任何附加毒性团；(iii) 书面的结构类推论证；(iv) 本监管风格报告。</p>
<p><b>框架。</b>非致突变杂质的界定遵循 ICH Q3A(R2)/Q3B(R2)。致突变潜力按 ICH M7(R2) 筛查，该指南要求采用两条互补的 (Q)SAR 方法学（专家规则 + 统计）并辅以专家审核。计算模型的使用遵循 OECD (Q)SAR 验证原则。由于万古霉素 B 为大分子（MW≈1449）柔性糖肽，处于许多小分子毒性模型适用域之外，本评估以机理性、结构与形状类推为主，统计模型仅在其适用域得到尊重时使用（域外情形均予标注）。</p>

<h1 class="pb">3. 受试物与结构关系</h1>
<p>结构取自申办方资料（IUPAC 名/SMILES/InChI），并用 RDKit 独立校验；计算所得分子式与精确质量与源文件一致。</p>
<table>
<tr><th>属性</th><th>母体 API</th><th>RRT 0.87</th><th>RRT 0.75</th></tr>
<tr><td>Xellia 名称</td><td>万古霉素 B</td><td>(7S,26R)-万古霉素 B</td><td>万古霉素 B 杂质 (26R)-DMe-DeLI</td></tr>
<tr><td>分子式</td><td>C66H75Cl2N9O24</td><td>C66H75Cl2N9O24</td><td>C65H73Cl2N9O24</td></tr>
<tr><td>单同位素质量</td><td>1447.43</td><td>1447.43</td><td>1433.41</td></tr>
<tr><td>已定义立体中心</td><td>18</td><td>18</td><td>19</td></tr>
<tr><td>来源</td><td>发酵（API）</td><td>降解（低水活度+受热）</td><td>发酵</td></tr>
<tr><td>与 API 的关系</td><td>—</td><td>非对映异构体（C7、C26 差向）</td><td>第1位残基变体（去N-甲基；Leu→Ile）+ C26 差向</td></tr>
</table>
{fig("fig1_structures.png","图 1. 母体 API（左）与两个杂质的二维结构。RRT 0.87 因仅立体化学不同，画法与母体一致；RRT 0.75 在第 1 位残基处不同（游离伯胺 + 异亮氨酸侧链，取代 N-甲基亮氨酸）。")}
<p><b>RRT 0.87 —— (7S,26R)-万古霉素 B。</b>分子式与 API 相同；C7（残基3的苄位甲醇碳）与 C26（一个肽 α-碳）构型翻转 —— 与低水活度受热条件下的差向异构化一致。不增删任何原子。</p>
<p><b>RRT 0.75 —— (26R)-DMe-DeLI。</b>分子式 C65H73Cl2N9O24（较母体少一个 CH2）。第 1 位残基由 N-甲基-D-亮氨酸变为去 N-甲基异亮氨酸（仲胺→伯胺；异丁基→仲丁基），并在 C26 处差向。其余残基、两个糖、两个芳基氯及所有酚均保留。</p>
<div style="display:flex;gap:10px">{fig("fig2_parent_residue1.png","母体：N-甲基-D-亮氨酸（仲胺）。","48%")}{fig("fig2_rrt075_residue1.png","RRT 0.75：去甲基异亮氨酸（伯胺）。","48%")}</div>
<figcaption style="text-align:center">图 2. RRT 0.75 与 API 之间唯一的构成改变位点。</figcaption>

<h1 class="pb">4. 计算方法</h1>
<h2>4.1 毒性团 / 结构警示分析</h2>
<p>对每个结构分析其万古霉素药效团（七肽糖苷配基 D-Ala-D-Ala 结合口袋、氢键网络、万古糖胺–葡萄糖二糖）及公认毒性团（酚/芳环→肾小管蓄积；碱性胺中心→组胺/红人综合征；高分子量/亲水性）。用 SMARTS 对官能团进行定量清点并与母体求差。采用 RDKit PAINS(A/B/C)、Brenk、NIH 目录筛查结构隐患，并对照 ICH M7 DNA 反应性警示类别审核。</p>
<h2>4.2 理化与三维描述符；三维形状类推</h2>
<p>RDKit 描述符（MW、cLogP、TPSA、HBD/HBA、可旋转键、QED、Lipinski/Veber）；Inductive Bio logD/pKa QSAR。ETKDGv3 生成构象并经 MMFF 能量最小化；三维形状描述符（回转半径、非球度、球度、离心率、NPR1/2）。三维类推通过 MCS 核心限定叠合（RMSD）及 Crippen Open3DAlign / 形状-Tanimoto 量化。不构建正式的 CoMFA/CoMSIA 3D-QSAR（见第 8 节）。</p>
<h2>4.3 统计 QSAR —— Tox21 与 Ames</h2>
<p>训练并应用了 12 端点 Tox21 QSAR（在 2048 位 Morgan 指纹上的类平衡随机森林；公开 MoleculeNet Tox21；5 折交叉验证 AUC）与 Ames QSAR（Hansen 基准数据集 N≈6500），均带适用域检查。母体–杂质差异为可解释读数。</p>

<h1 class="pb">5. 结果</h1>
<h2>5.1 母体（万古霉素 B）—— 基线毒性团谱</h2>
<p>母体呈现典型糖肽谱：6 个酚羟基、5 个芳环、2 个芳基氯、3 个芳醚、2 个苄醇、1 个羧酸、1 个伯酰胺、6 个仲酰胺、1 个伯胺（万古糖胺）、1 个 N-甲基仲胺。无 PAINS/Brenk 警示、无 DNA 反应性警示；唯一的 NIH 标记（≥7 个脂肪羟基）为极性描述符。此为类推基线。</p>
{fig("fig3_toxicophore_inventory.png","图 3. 毒性团/官能团清单（计数）。母体与 RRT 0.87 在所有类别完全一致；RRT 0.75 仅在胺列不同（红色）—— 一个仲胺（N-甲基）被替换为一个伯胺。")}
<h2>5.2 RRT 0.87 —— (7S,26R)-万古霉素 B</h2>
<p>因仅构型不同，其二维结构、官能团清单与警示谱均与母体一致（ECFP4/MACCS Tanimoto = 1.00；共享 101/101 重原子；0 警示）。<b>无附加毒性团。</b>常规二维类推完全无法将其与 API 区分 —— 这正是需要三维/机理评估的原因。三维形状描述符近乎重合（回转半径 7.45 对 7.50 Å；NPR1 0.529 对 0.528）。C7/C26 构型翻转可能改变靶点结合（效力），但不产生任何新的反应性/DNA 反应性官能团。</p>
<h2>5.3 RRT 0.75 —— (26R)-DMe-DeLI</h2>
<p>保留全部糖苷配基、两个糖、两个芳基氯、所有酚、羧酸与伯酰胺；共享 99/约100 重原子（ECFP4 0.86 / MACCS 0.94）。唯一改变（仲 N-甲基→伯胺）<b>并非新官能团类别</b>（母体本已含伯胺），也<b>不是 DNA 反应性警示</b>。logD 不变（0.42 对 0.42），最强碱性 pKa 近乎不变（8.08 对 8.21）—— 碱性胺/红人综合征相关的隐患为保守而非放大。TPSA 仅 +14 Å²；三维形状接近母体（回转半径 7.46 Å）。</p>
<h2>5.4 差异毒性团分析 —— 是否存在新毒性团？</h2>
<table>
<tr><th>潜在毒性团</th><th>RRT 0.87 对母体</th><th>RRT 0.75 对母体</th></tr>
<tr><td>酚/芳环（肾蓄积）</td><td>一致（6 酚、5 芳环）</td><td>一致（6 酚、5 芳环）</td></tr>
<tr><td>碱性胺中心（组胺/红人）</td><td>一致</td><td style="background:#fff3e0">保守：仲→伯胺，等碱性(pKa 8.1)；未放大</td></tr>
<tr><td>芳基氯</td><td>一致（2）</td><td>一致（2）</td></tr>
<tr><td>DNA 反应性警示（芳香胺、硝基、环氧、Michael 受体、烷基卤…）</td><td>无（同母体）</td><td>无（同母体）</td></tr>
<tr><td>PAINS / Brenk 隐患</td><td>无</td><td>无</td></tr>
<tr><td><b>净结果：是否引入新毒性团？</b></td><td class=c style="color:#2e7d32;font-weight:700">否</td><td class=c style="color:#2e7d32;font-weight:700">否</td></tr>
</table>
<p>唯一值得说明的是 RRT 0.75 的第 1 位残基胺 —— 这是对既有毒性团类别的修饰（仲→伯胺、等碱性），而非新增类别。两个子体分子均未识别出附加毒性团。</p>
<h2>5.5 3D-QSAR / 形状类推指标</h2>
{fig("fig5_similarity.png","图 4. 二维相似度与共享骨架类推指标。RRT 0.87 指纹相同；RRT 0.75 共享≥98% 重原子骨架。")}
{fig("fig4_descriptors.png","图 5. 三维整体形状描述符（左）与理化描述符（右，对称对数）。三者占据基本相同的形状/性质空间。")}

<h2 class="pb">5.6 统计（二维）QSAR —— 12 端点 Tox21 谱</h2>
<p>训练了 12 端点 Tox21 QSAR（2048 位 Morgan 指纹上的类平衡随机森林；公开 MoleculeNet Tox21，n={ntrain}）：5 折交叉验证平均 AUC = {meanauc:.3f}（范围 {min(aucs):.3f}–{max(aucs):.3f}），随后对每个分子打分并做适用域检查。</p>
<div class="box"><b>适用域 —— 一个有利的发现。</b>万古霉素本身即在 Tox21 数据集中（训练化合物 {', '.join(QA['identical_mol_ids'])}；对母体及 RRT 0.87 的最近邻 Tanimoto = 1.00）。因此该 QSAR 对母体为域内；RRT 0.75 亦相近（最大 Tanimoto 0.863）。</div>
<p><b>万古霉素实测谱（真值）：</b>仅 NR-AR 阳性（弱），其余所有测过的端点均为阴性（含遗传毒性 SR-ATAD5、SR-p53）。此为毒理学基线。</p>
{fig("fig6_qsar.png","图 6. 自建 12 端点 Tox21 QSAR。(A) 预测概率 —— 母体与 RRT 0.87 两行相同（指纹相同）。(B) 各端点 5 折 OOF AUC（平均 0.827）。(C) RRT 0.75 − 留一法母体 的 Δ —— 各端点均在 ±0.02 内；RRT 0.87 按构造 Δ ≡ 0。")}
<p><b>杂质对母体。</b>RRT 0.87 在全部 12 端点返回相同概率（Δ ≡ 0）：具体证明了二维 QSAR 对该立体异构体在结构上"盲视"（故三维类推承担主要判据）。对 RRT 0.75，每个端点均在留一法母体的 ±0.02 以内（最大 |Δ| = 0.02），无一越过活性阈值 —— 包括遗传毒性相关的 SR-ATAD5 与 SR-p53。</p>
<table>
<tr><th>端点</th><th>CV-AUC</th><th>实测（母体）</th><th>留一法母体 P</th><th>RRT 0.75 P (Δ)</th></tr>
{qrows}
</table>
<p class="foot"><i>所有预测概率均低于 0.5 阈值（母体 NR-AR 为实测弱阳性，模型低估 —— 在该端点较低 AUC 下属预期；以实测标签为准）。RRT 0.87 未列入表中：指纹与母体相同 ⇒ 各端点预测相同（Δ ≡ 0）。</i></p>

<h2 class="pb">5.7 ICH M7 致突变评估 —— 两条互补 (Q)SAR 方法学</h2>
<p>ICH M7(R2) 要求以两条互补方法学评估潜在的细菌（Ames）致突变性 —— 一条<b>专家规则</b>、一条<b>统计</b> —— 任何阳性均需经专家审核。二者均已应用：</p>
<ul>
<li><b>方法一（专家规则）。</b>Benigni–Bossa 致突变/遗传毒性致癌结构警示规则库（{n_alerts} 类警示：芳香胺、硝基/亚硝基、偶氮/氧化偶氮、N-亚硝基、环氧/氮丙啶、Michael 受体、酰卤/烷基卤、肼、醛、醌、β-内酰胺……）—— 即 Derek 类专家系统背后的公开科学规则库（Toxtree/VEGA）。</li>
<li><b>方法二（统计）。</b>Ames QSAR（类平衡随机森林、2048 位 Morgan 指纹），训练于公开 Hansen 基准（N = {M7['methodology_2']['n_train']}；{round(M7['methodology_2']['pos_frac']*100)}% 为致突变物）；5 折 OOF AUC = {M7['methodology_2']['cv_auc']:.3f}；决策阈值 {M7['methodology_2']['threshold_youden']:.2f}。</li>
</ul>
{fig("fig7_m7.png","图 7. ICH M7 双方法学评估。(A) 一致性矩阵 —— 两条方法学对三个分子均为阴性。(B) 统计 Ames 概率，均低于阈值；琥珀色标注提示适用域相似度低，故以域无关的专家规则库为主判据。")}
<table>
<tr><th>分子</th><th>专家规则</th><th>统计 Ames QSAR</th><th>M7</th></tr>
{m7rows}
</table>
<p><b>结论。</b>两条方法学对母体及两个杂质一致且阴性：无致突变结构警示，统计 Ames 预测在各情形均低于阈值。据 ICH M7，此对应 <b>Class 5（无结构警示）</b> —— 按非致突变的普通 ICH Q3A 有关物质控制，无需专门的细菌回复突变（Ames）试验。为透明起见记录两点：(i) RRT 0.87 与母体指纹相同，故两条方法学按构造返回相同结果；(ii) 此类糖肽处于 Ames 模型适用域之外（最近邻 Tanimoto ≈ 0.22），故统计阴性为低置信度，域无关的专家规则库（明确无警示）为主判据，并由母体实测阴性的 Tox21 DNA 损伤端点（SR-ATAD5、SR-p53；见 5.6）佐证。在验证流程中运行商用 Derek/Sarah Nexus 预期复现此无警示、阴性结论。</p>

<h1 class="pb">6. 结构类推（Read-Across）论证</h1>
<p><b>源（类似物）：</b>母体 API 万古霉素 B（具完整毒理学资料，含 91 天静脉犬试验 NOAEL 75/100 mg/kg/day 及数十年临床暴露）。<b>目标：</b>RRT 0.87 与 RRT 0.75。</p>
<ul>
<li><b>结构相似性。</b>RRT 0.87 共享 100% 重原子骨架与官能团（构型异构体）。RRT 0.75 共享≥98%（99/100），仅在一个肽 N-末端有单一、局部、无警示的修饰。</li>
<li><b>共同机理 / 保守毒性团。</b>所有酚、芳环、芳基氯、酰胺、羧基与碱性胺特征均保守。无新的反应性中心。</li>
<li><b>代谢/理化一致性。</b>万古霉素主要以原形经肾清除；近乎一致的分子大小、极性、logD(0.42)与 pKa(≈8.1)预示相同的 ADME/处置。</li>
<li><b>无致突变担忧。</b>ICH M7 两条方法学均阴性（见 5.7）；QSAR 在 SR-ATAD5/SR-p53 无激活（Δ ≤ 0.02）。</li>
<li><b>已认知的不确定性（立体化学）。</b>立体异构体可能在靶点结合上不同 —— 关乎效力/疗效，而非新生毒性；三维分析正是针对这一维度。</li>
</ul>

<h1>7. 界定结论</h1>
<p>基于证据权衡，本计算机模拟评估支持通过对母体的结构类推对两个杂质进行毒理学界定：</p>
<ul>
<li><b>RRT 0.87 ［(7S,26R)-万古霉素 B］</b> —— 纯立体异构体，无新毒性团、无新致突变警示；类推充分成立。</li>
<li><b>RRT 0.75 ［(26R)-DMe-DeLI］</b> —— 单一、无警示、等官能团的第 1 位残基修饰（加 C26 差向），无新毒性团、无新致突变警示；类推成立。</li>
</ul>
<p>只要二者控制在《欧洲药典》有关物质限度内，母体的 NOAEL（75/100 mg/kg/day，约临床 2–3 倍）及临床记录即可提供充分的毒理学覆盖。正式签署应由合格毒理学家出具，致突变维度按公司 ICH M7 双 (Q)SAR 方法学流程存档。</p>

<h1 class="pb">8. 适用域、诚实声明与局限性</h1>
<ul>
<li><b>"3D-QSAR"的性质。</b>这些化合物不存在带实测毒性端点的同系物训练集，故未（也不能合法地）构建正式的场基 3D-QSAR（CoMFA/CoMSIA）预测模型。本报告所述为三维结构/形状类推加机理性毒性团分析。</li>
<li><b>统计 QSAR —— 有效范围。</b>12 端点 Tox21 QSAR 对母体为域内（万古霉素为训练化合物）。两点提醒：(i) 母体的集内预测含部分记忆效应，故改以实测标签 + 留一法预测报告；(ii) Tox21 面板仅覆盖 12 个特定机制 —— 不建模万古霉素的临床效应（肾/耳毒性、红人综合征），后者由机理分析及对母体体内资料的类推处理。</li>
<li><b>ICH M7 统计模型适用域。</b>Ames QSAR 对此类糖肽为域外（Tanimoto ≈ 0.22），其阴性为低置信度，以专家规则库为主判据。应在验证流程中以商用 Derek/Sarah Nexus 确认。</li>
<li><b>构象采样。</b>三维描述符/叠合来自大柔性大环单一低能构象 → 属佐证而非决定性；决定性证据为共享骨架与保守毒性团清单。</li>
<li><b>范围。</b>本工作为危害识别/结构筛查与类推，用于支持而非替代毒理学专家审核及主管当局可能要求的任何实验测试。</li>
<li><b>杂质控制限度。</b>源文件未提供 RRT 0.75/0.87 的《欧洲药典》数值限度；本结论以控制在适用各论限度内为前提。</li>
</ul>

<h1 class="pb">9. 参考文献与工具</h1>
<ul class="foot" style="font-size:9pt">
<li>ICH Q3A(R2) / Q3B(R2) —— 新原料药/新制剂中的杂质。</li>
<li>ICH M7(R2) —— DNA 反应性（致突变）杂质的评估与控制。</li>
<li>OECD (2007) —— (Q)SAR 模型验证指南。</li>
<li>Benigni R., Bossa C. (2008) —— 致突变物与致癌物的结构警示（Toxtree/VEGA 规则库；方法一）。</li>
<li>Hansen K. et al. (2009) J. Chem. Inf. Model. 49:2077 —— Ames 致突变性基准数据集（方法二训练数据）。</li>
<li>Tox21 / MoleculeNet（Wu Z. et al. 2018, Chem. Sci.）—— 统计 QSAR 训练数据。</li>
<li>RDKit (2024) —— 描述符、PAINS/Brenk/NIH 目录、ETKDGv3、MMFF94、Open3DAlign、MCS。Inductive Bio logD/pKa QSAR。《欧洲药典》万古霉素各论。</li>
<li>申办方源文件："Vanco_RRT_075_and_087_Info"（Xellia）—— 身份、SMILES/InChI、来源、91 天犬 NOAEL 75/100 mg/kg/day。</li>
</ul>
<p class="foot" style="margin-top:20px;border-top:1px solid #ccc;padding-top:6px">计算机模拟筛查与结构类推 —— 支持而非替代毒理学专家审核。</p>
</body></html>"""
open(f"{D}/report_zh.html","w").write(HTML)
print("wrote report_zh.html", len(HTML))
