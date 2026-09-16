"""Build the Chinese-language ADMET + safety-panel PDF report."""

from __future__ import annotations
import os
from reportlab.lib.units import mm
from reportlab.platypus import (Spacer, PageBreak, NextPageTemplate, KeepTogether,
                                Paragraph)

import pdfkit as K
import report_data as R
import safety_panel as SP

OUT = R.OUT
FIG = f"{OUT}/figures"
PDF = f"{OUT}/admet_analysis_report_zh_v2.pdf"
DATE = "2026-09-16"

CW = (210 - 36) * mm          # usable content width


def cover(S, d, v):
    el = []
    el.append(Spacer(1, 42 * mm))
    el.append(K.P("FDA 已批准药物面板的<br/>ADMET 与类药性计算分析", S["title"]))
    el.append(K.P("附：Safety Panel 对比与计算覆盖度分析", S["subtitle"]))
    el.append(Spacer(1, 7 * mm))
    el.append(K.P(
        "本报告对 30 个 FDA 已批准的小分子药物进行系统的吸收、分布、代谢、排泄与毒性 "
        "(ADMET) 计算表征，并以已知药理学阳性对照药物对预测结果做 demo 验证。"
        "在此基础上，将预测终点映射到药物安全药理学中通用的次级药理 (secondary "
        "pharmacology) 靶点面板 Safety-44 / Safety-77 / SafetyScreen87，"
        "量化计算方法能覆盖与不能覆盖的部分。", S["body"]))
    el.append(Spacer(1, 10 * mm))

    meta = [
        ["报告版本", "v2 (含 Safety Panel 对比与计算覆盖度章节)"],
        ["生成日期", DATE],
        ["化合物数", f"{len(d)} 个 FDA 已批准小分子药物"],
        ["预测终点", "41 个 ADMET-AI 终点 + 17 个理化/规则描述符"],
        ["结构警示目录", "PAINS A/B/C、Brenk、NIH (RDKit FilterCatalog)"],
        ["预测模型", "ADMET-AI (Chemprop-RDKit 集成, TDC ADMET Benchmark Group)"],
        ["参考背景", "DrugBank 已批准药物参考集百分位数 (n = 2,579)"],
        ["验证结果", f"DMPK/hERG/BBB 类阳性对照 {v['dmpk_rec']}/{v['dmpk_tot']} 命中；"
                     f"Tox21 核受体类 {v['nr_rec']}/{v['nr_tot']} 命中"],
    ]
    el.append(K.make_table(["项目", "内容"], meta, [38 * mm, CW - 38 * mm],
                           S, font_size=8.4))
    el.append(Spacer(1, 10 * mm))
    el.append(K.callout(
        "<b>用途声明</b>　本报告为计算毒理学筛选与优先级排序工具，用于研究阶段的"
        "风险分诊。所有预测均为基于公开生物活性数据训练的统计估计，存在明确的适用域"
        "边界，<b>不构成</b>法规毒理学评价或临床安全性结论，也<b>不能替代</b> GLP "
        "体外/体内安全药理学试验。", S, colour=K.ACCENT_L, border=K.ACCENT))
    return el


def sec_intro(S):
    el = [K.P("1　引言", S["h1"])]
    el.append(K.P(
        "药物研发中约 <b>20–40%</b> 的临床前与临床期项目终止可归因于安全性问题，"
        "其中相当一部分源于化合物对预期治疗靶点以外的分子发生作用，即<b>脱靶</b>"
        "(off-target) 或<b>次级药理</b>作用。为在早期识别这类风险，制药行业建立了"
        "标准化的体外次级药理靶点面板；与之并行，基于机器学习的 ADMET 预测在近年"
        "已能对数十个吸收、代谢与毒性终点给出可用的早期估计。", S["body"]))
    el.append(K.P(
        "本报告的目的有两层。第一层是<b>方法学演示</b>：在一个成分已知、药理学背景"
        "清楚的 30 药物面板上完整跑通 ADMET 计算流程，并用既定药理学事实检验预测"
        "的可信度——如果模型连已上市药物的公认作用都无法重现，那么它对新化合物的"
        "预测同样不可信。第二层是<b>能力边界的量化</b>：把计算终点逐一对到 "
        "Safety-44 / Safety-77 / SafetyScreen87 面板的靶点上，明确回答"
        "「计算能替代多少体外安全筛选」这一问题。", S["body"]))
    el.append(K.P(
        "选择 FDA 已批准药物作为测试集有两个好处：其一，它们的药理学、代谢与安全性"
        "谱系在文献与监管文件中有充分记载，可作为 ground truth；其二，它们全部通过了"
        "完整的监管审评，因此在类药性指标上应当整体表现良好，任何系统性的高风险预测"
        "都提示模型本身的偏倚而非化合物的问题。", S["body"]))
    return el


def sec_methods(S):
    el = [K.P("2　方法", S["h1"])]

    el.append(K.P("2.1　化合物面板与结构标准化", S["h2"]))
    el.append(K.P(
        "面板包含 30 个 FDA 已批准（或曾获批准后因安全性撤市）的小分子药物，覆盖"
        "非甾体抗炎药、他汀、质子泵抑制剂、抗精神病药、抗心律失常药、抗组胺药、"
        "唑类抗真菌药、氟喹诺酮、激酶抑制剂、核受体配体等类别。面板刻意纳入若干"
        "<b>已知阳性对照</b>：特非那定、西沙必利、胺碘酮、氟哌啶醇、奎尼丁 (hERG/QT)；"
        "酮康唑、维拉帕米 (CYP3A4 与 P-糖蛋白)；他莫昔芬 (ER)、比卡鲁胺 (AR)、"
        "阿那曲唑 (芳香化酶)、罗格列酮 (PPAR-γ)，用于第 3.6 节的验证。", S["body"]))
    el.append(K.P(
        "每个 SMILES 的分子式均与公开数据库参考值逐一核对后方进入流程。结构标准化使用"
        "RDKit <font face='Courier'>MolStandardize</font>：Cleanup → 选取最大片段 → "
        "去电荷 → 生成规范 SMILES；30 个结构全部解析成功，1 个在标准化中发生改变。",
        S["body"]))

    el.append(K.P("2.2　理化描述符与类药性规则", S["h2"]))
    el.append(K.P(
        "由 RDKit 计算分子量、cLogP (Crippen)、拓扑极性表面积 (TPSA)、氢键供体/受体数、"
        "可旋转键数、芳环数与总环数、sp³ 碳比例、摩尔折射率与 QED。在此基础上评估四组"
        "经验规则：<b>Lipinski 五规则</b>（MW≤500、cLogP≤5、HBD≤5、HBA≤10，允许 1 项违反）、"
        "<b>Veber 规则</b>（可旋转键≤10 且 TPSA≤140 Å²）、<b>Ghose 过滤</b>与 "
        "<b>Egan 卵</b>。", S["body"]))

    el.append(K.P("2.3　结构警示", S["h2"]))
    el.append(K.P(
        "使用 RDKit <font face='Courier'>FilterCatalog</font> 的五个目录做子结构匹配："
        "PAINS A/B/C（泛测定干扰化合物）、Brenk（反应性或不宜成药基团）与 NIH 过滤集。"
        "需要强调的是，结构警示标记的是<b>需要复核的子结构</b>，而非毒性判定；本面板中"
        "多个已上市药物同样带有警示。", S["body"]))

    el.append(K.P("2.4　ADMET 终点预测", S["h2"]))
    el.append(K.P(
        "使用 <b>ADMET-AI</b>（Swanson 等, <i>Bioinformatics</i> 2024）。该工具在 "
        "Therapeutics Data Commons (TDC) ADMET Benchmark Group 的 41 个数据集上训练 "
        "Chemprop-RDKit 图神经网络集成模型，输出 <b>41 个终点</b>：其中分类终点包括 "
        "hERG 阻断、药物性肝损伤 (DILI)、Ames 致突变、临床试验毒性 (ClinTox)、致癌性、"
        "五个主要 CYP 亚型的抑制与底物判定、P-糖蛋白抑制、血脑屏障穿透、人肠道吸收、"
        "口服生物利用度，以及 12 个 Tox21 核受体 (NR-) 与应激反应 (SR-) 终点；"
        "回归终点包括水溶解度、亲脂性、Caco-2 渗透性、血浆蛋白结合率、稳态分布容积、"
        "半衰期、肝细胞与微粒体清除率、急性毒性 LD50 等。", S["body"]))
    el.append(K.P(
        "全部计算在 CPU 上离线完成。分类终点输出为 0–1 的概率，本报告统一采用 "
        "<b>0.5</b> 作为判读阈值，并在第 6 节讨论该阈值在本面板上的局限。", S["body"]))

    el.append(K.P("2.5　参考背景百分位数", S["h2"]))
    el.append(K.P(
        "ADMET-AI 内置一个由 <b>2,579 个 DrugBank 已批准药物</b>构成的参考集，"
        "对每个终点给出被测化合物在该参考集中的百分位数。百分位数把绝对预测值转换为"
        "「相对于已上市药物处于什么位置」的可比尺度：例如 hERG 预测 0.95 对应第 93 "
        "百分位，意味着仅约 7% 的已批准药物被预测出更高的 hERG 风险。", S["body"]))
    el.append(K.P(
        "说明：原计划使用 ChEMBL 已批准药物百分位参考集；本次运行中实际使用的是 "
        "ADMET-AI 随包分发的 DrugBank 已批准药物参考集，两者用途相同，本报告中所有"
        "百分位数均按 DrugBank 参考集标注。", S["note"]))

    el.append(K.P("2.6　Demo 验证设计", S["h2"]))
    el.append(K.P(
        "验证不依赖模型自身的训练/测试划分，而是用<b>外部既定事实</b>作为标签："
        "hERG 阳性取自 CredibleMeds 已知 TdP 风险清单及因 QT 延长撤市的药物；"
        "CYP 与 P-糖蛋白抑制剂取自 FDA《药物相互作用》临床指数抑制剂列表；"
        "核受体阳性取该药物本身获批的药理机制。阴性对照选取药理学上明确不作用于"
        "该靶点的药物（如二甲双胍、阿司匹林、咖啡因）。对每个终点报告阳性对照的"
        "预测值、其在 30 药物面板中的排名，以及在阳性/阴性对照数≥3 时的 ROC-AUC。",
        S["body"]))
    return el


def sec_results(S, d, v):
    el = [K.P("3　结果", S["h1"])]

    # 3.1
    el.append(K.P("3.1　理化性质分布", S["h2"]))
    el.append(K.P(
        f"面板的中位分子量为 {d['MW'].median():.0f} Da，中位 cLogP {d['cLogP'].median():.2f}，"
        f"中位 TPSA {d['TPSA'].median():.0f} Å²，整体落在口服小分子药物的典型区间内 "
        f"(表 1)。分布的两端分别由二甲双胍 (MW {d.set_index('name').loc['Metformin','MW']:.0f} Da，"
        f"cLogP {d.set_index('name').loc['Metformin','cLogP']:.2f}，高度亲水) 与阿托伐他汀、"
        f"胺碘酮 (MW > 600 Da 或 cLogP > 6) 界定。", S["body"]))
    hdr, rows = R.table_physchem(d)
    el.append(K.P("<b>表 1</b>　理化性质与类药性描述符汇总 (n = 30)", S["tcap"]))
    el.append(K.make_table(hdr, rows, [40*mm, 22*mm, 34*mm, 36*mm, CW-132*mm], S,
                           font_size=8))
    el.append(Spacer(1, 4 * mm))
    el.append(K.figure(f"{FIG}/fig1_physicochemical_overview.png", CW, S,
                       "<b>图 1</b>　六个核心理化性质的分布。红色虚线为 Lipinski / Veber "
                       "阈值，深色实线为面板中位数。"))

    el.append(K.P("3.2　类药性规则符合度", S["h2"]))
    hdr, rows = R.table_rules(d)
    el.append(K.P("<b>表 2</b>　四组类药性规则的通过情况", S["tcap"]))
    el.append(K.make_table(hdr, rows, [56*mm, 20*mm, 18*mm, CW-94*mm], S,
                           font_size=8))
    el.append(K.P(
        "Lipinski 与 Veber 的通过率分别为 93% 与 90%，符合对已批准口服药物的预期。"
        "Ghose 过滤的通过率仅 57%，但未通过者多为阿司匹林、布洛芬、对乙酰氨基酚、"
        "咖啡因这类<b>低于</b>分子量下限 (160 Da) 的小分子药物，属于规则本身的边界"
        "效应而非真实的成药性缺陷——这正说明经验过滤规则不应被机械套用。", S["body"]))
    el.append(Spacer(1, 2 * mm))
    el.append(K.figure(f"{FIG}/fig2_lipinski_space.png", CW * 0.86, S,
                       "<b>图 2</b>　分子量–cLogP 化学空间。气泡面积正比于 QED，"
                       "颜色表示 Lipinski 违反项数，阴影区为五规则合规象限。"))

    el.append(K.P("3.3　QED 与结构警示", S["h2"]))
    n_alert = int((d["alert_total"] > 0).sum())
    n_pains = int((d["PAINS_total"] > 0).sum())
    el.append(K.P(
        f"面板平均 QED 为 {d['QED'].mean():.2f}。{n_alert}/30 个药物至少命中一条结构警示，"
        f"其中 {n_pains} 个命中 PAINS 目录。这一比例本身就是重要的参照：这些化合物全部"
        f"通过了 FDA 审评并在临床长期使用，说明结构警示是复核提示而非否决条件。"
        f"QED 最低的三个药物 (阿托伐他汀 {d.set_index('name').loc['Atorvastatin','QED']:.2f}、"
        f"胺碘酮 {d.set_index('name').loc['Amiodarone','QED']:.2f}、"
        f"二甲双胍 {d.set_index('name').loc['Metformin','QED']:.2f}) 分别代表分子过大、"
        f"过度亲脂与过度亲水三种不同的偏离方式，却都是各自领域的一线药物。", S["body"]))
    el.append(K.figure(f"{FIG}/fig3_qed_developability.png", CW * 0.68, S,
                       "<b>图 3</b>　按 QED 排序的成药性排名；⚑ 标注结构警示命中数。"))

    el.append(K.P("3.4　关键 ADMET 终点", S["h2"]))
    hdr, rows = R.table_endpoints(d)
    el.append(K.P("<b>表 3</b>　主要分类型 ADMET 终点在面板中的分布 (阈值 0.5)", S["tcap"]))
    el.append(K.make_table(hdr, rows, [36*mm, 26*mm, 20*mm, 28*mm, CW-110*mm], S,
                           font_size=7.8))
    el.append(K.P(
        f"最显著的观察是 hERG：30 个已批准药物中有 "
        f"{int((d['hERG']>0.5).sum())} 个 (中位预测值 {d['hERG'].median():.2f}) 被判为阳性。"
        f"真正的高风险药物确实排在最前 (西沙必利 0.98、特非那定 0.97、维拉帕米 0.97)，"
        f"但伊马替尼 (0.98)、吉非替尼 (0.96) 等临床上并非典型致 TdP 药物同样被标为高风险。"
        f"DILI 终点亦然：咖啡因得到全面板最高的 0.93。这提示分类阈值 0.5 在"
        f"「全部为上市药物」的集合上<b>敏感度高而特异性低</b>，实践中应使用排名或"
        f"百分位而非绝对阈值 (见第 4 节)。", S["body"]))
    el.append(Spacer(1, 3 * mm))
    el.append(K.figure(f"{FIG}/fig4_admet_heatmap.png", CW, S,
                       "<b>图 4</b>　30 个药物 × 23 个代表性 ADMET 终点的预测概率热图，"
                       "按平均预测风险降序排列。三个分区依次为毒性、DMPK 与吸收、"
                       "Tox21 核受体与应激反应终点。注意右侧 Tox21 分区整体接近 0。"))

    el.append(K.P("3.5　连续型终点", S["h2"]))
    hdr, rows = R.table_regression(d)
    el.append(K.P("<b>表 4</b>　回归型 ADMET 终点的分布", S["tcap"]))
    el.append(K.make_table(hdr, rows, [46*mm, 20*mm, 28*mm, 38*mm, CW-132*mm], S,
                           font_size=7.8))
    el.append(K.P(
        "回归终点的排序在化学上是合理的：他莫昔芬被预测为溶解度最低、亲脂性最高、"
        "血浆蛋白结合率最高，二甲双胍在所有这些维度上均处于相反极端。但绝对值需谨慎——"
        "他莫昔芬的血浆蛋白结合率被预测为 105.5%，超出物理可能范围，这是模型在训练"
        "数据分布边缘外推的直接证据，也是适用域评估必要性的实例。", S["body"]))

    el.append(K.P("3.6　高风险化合物", S["h2"]))
    hdr, rows = R.table_flagged(d)
    el.append(K.P("<b>表 5</b>　触发标记最多的 12 个化合物", S["tcap"]))
    el.append(K.make_table(hdr, rows, [26*mm, 38*mm, 15*mm, 15*mm, 15*mm, 17*mm,
                                       CW-126*mm], S, font_size=7.6))
    el.append(K.P(
        "排名靠前的化合物与其已知临床安全性特征有相当程度的吻合：西沙必利与特非那定"
        "均因 QT 延长撤市；酮康唑因肝毒性与 CYP3A4 强抑制而被大幅限制口服适应证；"
        "胺碘酮以多器官毒性著称。但伊马替尼与吉非替尼被列在最前，主要由 hERG 与 DILI "
        "两个终点驱动，而这两者在本面板上都表现出偏高的假阳性倾向。", S["body"]))
    return el


def sec_validation(S, d, v):
    el = [K.P("4　Demo 验证：预测能否重现已知药理学", S["h1"])]
    el.append(K.P(
        "本节用面板内药物<b>已经确立的</b>药理学事实检验预测。每个终点选取作用机制明确的"
        "阳性对照与药理学上无关的阴性对照，检验预测值能否将两者分开。", S["body"]))
    hdr, rows = R.table_validation()
    el.append(K.P("<b>表 6</b>　按终点的阳性/阴性对照验证结果", S["tcap"]))
    el.append(K.make_table(hdr, rows, [34*mm, CW-104*mm, 22*mm, 16*mm, 18*mm], S,
                           font_size=7.6))
    el.append(Spacer(1, 3 * mm))
    el.append(K.figure(f"{FIG}/fig5_demo_validation.png", CW * 0.92, S,
                       "<b>图 5</b>　各终点上阳性对照 (绿色菱形 = 命中，红色叉 = 未命中) "
                       "与阴性对照 (灰点) 的预测值。阴影区为 Tox21 核受体终点。"))

    el.append(K.P("4.1　DMPK 与心脏通道终点：验证通过", S["h2"]))
    el.append(K.P(
        "hERG、CYP3A4、CYP2D6、P-糖蛋白与血脑屏障五个终点在阳性/阴性对照之间实现了"
        "<b>完全分离</b>，ROC-AUC 均为 1.00。具体而言：五个已知致 TdP 药物的 hERG 预测值"
        "全部 ≥ 0.87 而六个阴性对照全部 ≤ 0.13；酮康唑的 CYP3A4 预测值 0.99 为全面板最高，"
        "与其作为 FDA 强效指数抑制剂的地位一致；P-糖蛋白终点上酮康唑 0.98 亦居首。"
        "血脑屏障终点正确区分了中枢活性药物 (利培酮 0.99、地西泮 0.98) 与外周限制性药物 "
        "(阿托伐他汀 0.24、二甲双胍 0.40)。", S["body"]))
    el.append(K.P(
        "唯二的例外具有诊断意义。奎尼丁的 P-糖蛋白预测值仅 0.35，尽管它是 FDA 列出的"
        "临床指数 P-gp 抑制剂；环丙沙星的 CYP1A2 预测值仅 0.015，而它是该酶的强效临床"
        "抑制剂——后者提示训练集所用的 Veith 高通量荧光底物测定与临床相关的 CYP1A2 "
        "抑制之间存在系统性差异。", S["body"]))

    el.append(K.P("4.2　Tox21 核受体终点：验证失败", S["h2"]))
    el.append(K.callout(
        f"<b>核心发现</b>　DMPK / hERG / BBB 类终点共 {v['dmpk_tot']} 个阳性对照，"
        f"命中 <b>{v['dmpk_rec']}</b> 个；Tox21 核受体类终点共 {v['nr_tot']} 个阳性对照，"
        f"仅命中 <b>{v['nr_rec']}</b> 个。四个以该受体为<b>获批作用机制</b>的药物"
        "——他莫昔芬 (ER, 预测 0.18)、比卡鲁胺 (AR, 0.06)、罗格列酮 (PPAR-γ, 0.30)、"
        "阿那曲唑 (芳香化酶, 0.27)——预测值全部低于 0.5，甚至在面板内的排名也未居首。",
        S, colour=K.colors.HexColor("#fdf0ee"), border=K.WARN))
    el.append(K.P(
        "这一结果有明确的机理解释，并非模型实现错误。Tox21 的 NR-ER 与 NR-AR 是"
        "<b>激动模式</b>的报告基因测定，而他莫昔芬与比卡鲁胺在该测定条件下是拮抗剂，"
        "本就不产生激动信号；Tox21 训练集阳性率极低（多数终点 < 8%），类别极度不平衡，"
        "模型倾向于输出接近 0 的概率；此外这些测定筛选的是环境化学物与工业品，其化学"
        "空间与处方药相差较远。", S["body"]))
    el.append(K.P(
        "值得注意的反例是：酮康唑在芳香化酶终点上取得全面板最高的 0.68，与其抑制 "
        "CYP19A1 并因此干扰甾体合成的已知作用一致；奥美拉唑在 AhR 终点上得到 0.83，"
        "与其作为 AhR 激动剂诱导 CYP1A2 的文献记载吻合。因此这些模型并非毫无信号，"
        "而是<b>不能用于判定某化合物是否为特定核受体的药理学配体</b>。", S["body"]))
    el.append(K.P(
        "这一结论直接决定了第 6 章映射分析的判读：即使某个 Tox21 终点在名义上"
        "对应安全药理面板中的一个核受体靶点，也不足以替代该靶点的体外测定。", S["body"]))
    return el


def sec_percentile(S, d):
    el = [K.P("5　参考百分位数背景", S["h1"])]
    el.append(K.P(
        "绝对概率值难以直接判读，而相对于已批准药物的位置更有操作意义。"
        "下表给出代表性药物在六个关键终点上的预测值及其在 2,579 个 DrugBank "
        "已批准药物中的百分位数。", S["body"]))
    hdr, rows = R.table_percentiles(d)
    el.append(K.P("<b>表 7</b>　代表性药物的预测值与 DrugBank 已批准药物百分位数", S["tcap"]))
    el.append(K.make_table(hdr, rows, [26*mm] + [(CW-26*mm)/6]*6, S, font_size=7.4))
    el.append(K.P(
        "百分位数显著改变了判读。以特非那定为例，其 hERG 预测值 0.97 对应第 96 百分位，"
        "确实处于已批准药物的极端；而酮康唑的 Ames 预测值 0.50 虽然「超过阈值」，"
        "却仅处于第 79 百分位，属于常见水平。将阈值判读替换为百分位判读，可以在很大"
        "程度上缓解第 3.4 节指出的高假阳性问题。", S["body"]))
    el.append(K.P(
        "参考集来源：ADMET-AI 随包分发的 DrugBank 已批准药物子集 (n = 2,579)，"
        "由 Swanson 等人整理并在 ADMET-AI 中用于百分位计算。", S["note"]))
    return el


def sec_discussion(S, d, v):
    el = [K.P("7　讨论", S["h1"])]
    el.append(K.P("7.1　计算 ADMET 的可信范围", S["h2"]))
    el.append(K.P(
        "本次验证给出的图景相当清晰，并且在两个方向上都是可操作的。在 DMPK 与心脏离子"
        "通道这一侧，预测质量足以支持实际决策：hERG、CYP3A4、CYP2D6、P-糖蛋白与血脑屏障"
        "五个终点在外部标签上实现了对照的完全分离。这些终点的共同特征是训练数据来自"
        "标准化程度高、化学空间以药物为主、且测定终点定义明确的公开数据集。", S["body"]))
    el.append(K.P(
        "在 Tox21 核受体这一侧，结论同样明确但方向相反：不应将其用于靶点级的药理学"
        "判定。需要强调的是，这不等于说这些模型没有价值——它们在其原本用途（环境"
        "化学物的内分泌干扰初筛）上仍有意义；问题在于把它们当作安全药理面板的"
        "计算替代品。", S["body"]))

    el.append(K.P("7.2　阈值与排序", S["h2"]))
    el.append(K.P(
        f"0.5 阈值在本面板上产生了 {int((d['hERG']>0.5).sum())}/30 的 hERG 阳性率与 "
        f"{int((d['DILI']>0.5).sum())}/30 的 DILI 阳性率。由于面板成员全部为已获批准的"
        f"药物，这一比例本身即构成对阈值的否定证据。更稳健的做法是在项目内部使用"
        f"<b>相对排序</b>或参考集百分位：前者用于同系列化合物的先导优化排序，"
        f"后者用于判断绝对风险是否超出已上市药物的常规范围。", S["body"]))

    el.append(K.P("7.3　适用域", S["h2"]))
    el.append(K.P(
        "他莫昔芬 105.5% 的血浆蛋白结合率预测是一个有用的警示：图神经网络在回归任务上"
        "不受物理边界约束，超出训练分布时会给出不可能的数值。实践中应对每个预测附带"
        "适用域判断（如与训练集的 Tanimoto 最近邻距离），并对超出物理范围的输出直接"
        "作废而非截断。", S["body"]))
    return el


def sec_conclusion(S, d, v):
    el = [K.P("8　结论", S["h1"])]
    items = [
        f"<b>类药性</b>：30 个 FDA 已批准药物中 28 个通过 Lipinski 五规则、27 个通过 "
        f"Veber 规则，平均 QED {d['QED'].mean():.2f}，符合对已上市口服药物的预期；"
        f"12 个药物携带结构警示，印证结构警示为复核提示而非否决标准。",

        f"<b>验证结论（分层）</b>：DMPK、hERG 与血脑屏障类终点的 {v['dmpk_tot']} 个阳性对照"
        f"命中 {v['dmpk_rec']} 个，ROC-AUC 达 1.00，可用于早期风险分诊；Tox21 核受体类终点"
        f"的 {v['nr_tot']} 个阳性对照仅命中 {v['nr_rec']} 个，不可用于靶点级药理学判定。",

        "<b>阈值</b>：固定 0.5 阈值在已批准药物面板上产生 70% 的 hERG 阳性率与 60% 的 "
        "DILI 阳性率，特异性不足；应改用参考集百分位或项目内相对排序。",

        "<b>安全面板覆盖度</b>：在 Bowes-44 面板的 44 个靶点中，现有 ADMET-AI 终点仅能"
        "覆盖约 1–2 个（hERG 可靠，核受体不可靠），覆盖率低于 5%；即便把 DMPK 类的 "
        "CYP 与 P-糖蛋白终点计入（它们严格来说不属于次级药理面板），可计算的靶点也只有"
        "十余个。计算预测与体外面板是<b>互补关系</b>，不存在替代关系（详见第 6 章）。",

        "<b>建议的使用方式</b>：将计算 ADMET 置于化合物合成之前，用于在类似物之间排序"
        "并提前规避明显的 hERG / CYP 抑制风险；体外安全药理面板置于候选化合物提名前后，"
        "用于系统性地排查计算方法完全无法触及的 GPCR、离子通道与转运体脱靶作用。",
    ]
    el += K.bullets(items, S["body"], S["bullet"])
    return el


def sec_limits(S):
    el = [K.P("9　局限性", S["h1"])]
    items = [
        "面板规模为 30 个化合物，用于方法演示与定性结论；表 6 中若干终点的阳性对照"
        "少于 3 个，因而未计算 ROC-AUC，相应结论为定性判断。",
        "阳性/阴性对照标签来自监管文件与文献记载的临床事实，与模型训练所用的体外测定"
        "终点不完全等价（环丙沙星/CYP1A2 的偏差即源于此）。",
        "ADMET-AI 的预测不含立体化学信息，手性药物的对映体无法区分。",
        "本报告第 6 章的面板构成数据来自公开文献与 CRO 产品资料，按<b>靶点家族粒度</b>"
        "呈现；三个面板的完整逐条靶点清单未在本报告中重建。",
        "全部结果为计算筛选，不构成法规毒理学评价，不能替代 GLP 体外/体内试验。",
    ]
    el += K.bullets(items, S["body"], S["bullet"])
    return el


def sec_refs(S):
    el = [K.P("10　参考文献", S["h1"])]
    for r in SP.REFERENCES:
        el.append(K.P(r, S["ref"]))
    return el


def sec_appendix(S, d):
    el = [K.P("附录 A　化合物面板", S["h1"])]
    rows = [[i + 1, r["name"], r["drug_class"], r["formula"],
             f"{r['MW']:.1f}", f"{r['cLogP']:.2f}", f"{r['QED']:.2f}",
             int(r["alert_total"])]
            for i, (_, r) in enumerate(d.iterrows())]
    el.append(K.make_table(
        ["#", "药物", "药理类别", "分子式", "MW", "cLogP", "QED", "警示"],
        rows, [8*mm, 28*mm, 46*mm, 26*mm, 16*mm, 16*mm, 14*mm, CW-154*mm], S,
        font_size=7.2))
    el.append(Spacer(1, 4 * mm))
    el.append(K.P("附录 B　输出文件", S["h1"]))
    files = [
        ("all_properties.csv", "30 药物 × 全部 122 列（理化、规则、警示、41 个 ADMET 终点、百分位）"),
        ("druglikeness_summary.csv", "理化与类药性规则汇总"),
        ("admet_predictions.csv", "41 个 ADMET-AI 终点预测值"),
        ("structural_alerts.csv", "逐条结构警示匹配记录"),
        ("flagged_compounds.csv", "触发风险标记的化合物"),
        ("demo_validation.csv", "阳性/阴性对照验证明细"),
        ("safety_panel_coverage.csv", "Safety panel 计算覆盖度映射（表 9）"),
        ("safety_panel_comparison.csv", "Safety-44 / 77 / 87 面板对比（表 8）"),
        ("analysis_report.md", "Markdown 版分析摘要"),
        ("analysis_object.pkl", "完整结果对象（pickle）"),
        ("figures/*.png / *.svg", "全部 6 张图的位图与矢量版本"),
    ]
    el.append(K.make_table(["文件", "内容"], files, [52*mm, CW-52*mm], S,
                           font_size=7.6))
    return el


# --------------------------------------------------------------------------- #
def build():
    K.register_fonts()
    S = K.styles()
    d = R.load()
    v = R.validation_stats()

    story = []
    story += cover(S, d, v)
    story.append(NextPageTemplate("main"))
    story.append(PageBreak())
    story += sec_intro(S)
    story += sec_methods(S)
    story.append(PageBreak())
    story += sec_results(S, d, v)
    story.append(PageBreak())
    story += sec_validation(S, d, v)
    story += sec_percentile(S, d)
    story.append(PageBreak())
    story += SP.section_safety_panel(S, CW, FIG)     # chapter 6
    story.append(PageBreak())
    story += sec_discussion(S, d, v)
    story += sec_conclusion(S, d, v)
    story.append(KeepTogether(sec_limits(S)))
    story += sec_refs(S)
    story += sec_appendix(S, d)

    doc = K.Doc(PDF, "FDA 药物面板 ADMET 与 Safety Panel 计算覆盖度分析 · v2")
    doc.build(story)
    print(f"✓ wrote {PDF} ({os.path.getsize(PDF)/1024:.0f} KB)")


if __name__ == "__main__":
    build()
