"""Chapter 6 — Safety panel comparison and computational coverage mapping.

Literature facts below are graded:
  [V]      verified against the primary paper, an official regulatory document,
           or the vendor's own published product document
  [V-sec]  stated by a credible secondary source (vendor catalogue, open-access
           commentary) but not read in the primary paper
  [NV]     could not be verified — reported as unknown, never guessed

Both Nature Reviews Drug Discovery papers (Bowes 2012, Brennan 2024) are
paywalled. The 44-target list used here is the Eurofins SafetyScreen44
implementation, which the vendor explicitly attributes to Bowes et al. 2012.
The full Safety-77 target list could not be obtained and is NOT reconstructed.
"""

from __future__ import annotations
from reportlab.lib.units import mm
from reportlab.platypus import Spacer

import pdfkit as K
import mapping_core as MC

# --------------------------------------------------------------------------- #
# Bowes-44 family composition  [V-sec: Eurofins SafetyScreen87 catalogue and
# Reaction Biology InVEST44 independently state 24/8/6/3/2/1]
# "covered" = an ADMET-AI endpoint exists for a target in this family
# "validated" = that endpoint also passed the section-4 positive-control test
# --------------------------------------------------------------------------- #
PANEL_FAMILIES = [
    {"family": "GPCR", "bowes44": 24, "covered": 0, "validated": 0},
    {"family": "离子通道", "bowes44": 8, "covered": 1, "validated": 1},
    {"family": "酶（非激酶）", "bowes44": 6, "covered": 0, "validated": 0},
    {"family": "转运体", "bowes44": 3, "covered": 0, "validated": 0},
    {"family": "核受体", "bowes44": 2, "covered": 1, "validated": 0},
    {"family": "激酶", "bowes44": 1, "covered": 0, "validated": 0},
]

TOTAL_44 = sum(f["bowes44"] for f in PANEL_FAMILIES)
TOTAL_COV = sum(f["covered"] for f in PANEL_FAMILIES)
TOTAL_VAL = sum(f["validated"] for f in PANEL_FAMILIES)


# --------------------------------------------------------------------------- #
# Table 8 — panel comparison
# --------------------------------------------------------------------------- #
COMPARISON_HDR = ["对比维度", "Safety-44 (Bowes 2012)",
                  "Safety-77 (Brennan / IQ DruSafe 2024)",
                  "SafetyScreen87 (Eurofins)"]

COMPARISON_ROWS = [
    ["提出背景",
     "AstraZeneca、GSK、Novartis、Pfizer 四家公司首次公开各自的体外药理谱系策略，"
     "凝练出一个「最小推荐面板」",
     "IQ 联盟 DruSafe 工作组对 <b>18 家公司</b>的问卷调查，总结当代实践并提出"
     "扩展的推荐靶点集",
     "CRO 商业产品：在 Bowes-44 基础上扩充，早于 Safety-77，与之相互独立"],

    ["靶点总数", "44", "77", "87"],

    ["靶点构成",
     "GPCR 24、离子通道 8、酶（非激酶）6、转运体 3、核受体 2、激酶 1",
     "激酶 20（26%）、非激酶酶类约 12%；<b>完整清单未公开获取</b>，"
     "其余家族构成本报告不作推断",
     "Bowes-44 全部 44 个 + 额外 43 个（以 GPCR 与离子通道为主）；"
     "GPCR 约 35–44、离子通道约 18–19、核受体约 6、激酶 3"],

    ["激酶占比",
     "1 / 44 = <b>2.3%</b>",
     "20 / 77 = <b>26%</b>（第二大靶点类别）",
     "3 / 87 = 3.4%"],

    ["核受体",
     "仅 <b>雄激素受体 AR</b> 与<b>糖皮质激素受体 GR</b> 两个",
     "文献描述「增强核受体等代表性不足家族的覆盖」，具体成员未获证实",
     "AR、GR + <b>ERα、PPARγ、孕激素受体 PR、RARα</b>"],

    ["心脏离子通道粒度",
     "通用型结合测定：hERG（膜制备）、L-型钙通道（DHP 位点）、"
     "钠通道（site 2，脑制备）、通用 Kv",
     "扩展离子通道与 GABA-A 功能评估（CRO 实现层面）",
     "分子层面明确：<b>Nav1.5</b>、Cav1.2 三个位点（维拉帕米/DHP/地尔硫䓬）、"
     "Cav2.2、Kv1.1、KATP 等"],

    ["靶器官覆盖",
     "以心血管与中枢神经系统为主",
     "明确扩展至心血管与中枢以外的器官系统",
     "以 GPCR/离子通道扩充为主，器官覆盖随之加宽"],

    ["检测模式",
     "以放射性配体<b>结合</b>测定为主（单浓度，常用 10 µM，n=2），"
     "酶类为酶活测定；原文建议结合测定与功能测定并行以提高灵敏度",
     "推荐功能性测定与剂量-反应；CRO 实现中激酶在<b>生理 1 mM ATP</b> 下检测",
     "以放射性配体结合为主，区分拮抗剂/激动剂放射配体；酶类与激酶为酶活测定"],

    ["剂量-反应",
     "默认单一浓度筛查，命中后再做浓度梯度确证",
     "强调命中后的浓度-反应确证与安全裕度（margin）计算",
     "标准产品为单浓度（10 µM），IC50/Ki 需另行订购"],

    ["法规对齐",
     "非法规强制。ICH S7A 仅将「提示潜在不良作用的配体结合或酶测定数据」"
     "列为设计核心组合试验时<b>应考虑</b>的因素之一",
     "论文摘要指出监管机构「<b>正日益要求</b>」提供已知不良反应关联靶点的活性数据",
     "厂商宣称「ICH S7A/B 对齐」，实指所筛靶点对应 S7A 关注的器官系统，"
     "而非满足任何 S7A/S7B 要求"],

    ["代表 CRO 产品",
     "Eurofins SafetyScreen44™ (P270)、Reaction Biology InVEST44",
     "Reaction Biology InVEST77",
     "Eurofins SafetyScreen87 (P342 / PP223)"],
]


# --------------------------------------------------------------------------- #
# Table 9 — coverage mapping (literature half; drug examples filled from data)
# --------------------------------------------------------------------------- #
# (family, panel target(s), ADMET endpoint key or None, in-44, in-87, note)
COVERAGE_SPEC = [
    ("离子通道", "hERG (Kv11.1)", "hERG", "是", "是",
     "<b>唯一干净映射。</b>面板为膜制备结合测定，TDC 标签多源于膜片钳 IC50，"
     "靶点相同但测定原理不同", "covered"),

    ("核受体", "雄激素受体 (AR)", "NR-AR", "是", "是",
     "靶点名义对应，但面板为放射配体<b>结合</b>测定，Tox21 NR-AR 为"
     "<b>激动模式</b>转录报告测定；本报告第 4.2 节验证显示比卡鲁胺预测值仅 0.06，"
     "<b>不可用</b>", "nominal"),

    ("核受体", "糖皮质激素受体 (GR)", None, "是", "是",
     "ADMET-AI / Tox21 终点集中<b>无对应终点</b>", "none"),

    ("核受体", "雌激素受体 ERα、PPARγ、PR、RARα", "NR-ER", "否", "是",
     "仅存在于 SafetyScreen87，不在 Bowes-44 中。ER 与 PPARγ 有名义对应的 Tox21 终点，"
     "但第 4.2 节验证显示他莫昔芬 0.18、罗格列酮 0.30，<b>均不可用</b>", "nominal87"),

    ("GPCR", "5-HT2B（瓣膜病）、5-HT1A/1B/2A、多巴胺 D1/D2S、"
             "肾上腺素 α1A/α2A/β1/β2、毒蕈碱 M1/M2/M3、组胺 H1/H2、"
             "阿片 µ/δ/κ 等共 24 个", None, "是", "是",
     "<b>全部无计算终点。</b>这是覆盖缺口最大的家族（占 Bowes-44 的 54.5%）。"
     "5-HT2B 的风险来自<b>激动</b>作用，方向性判断本身即超出当前预测能力", "none"),

    ("离子通道", "Cav1.2（L-型钙通道）、Nav1.5、GABA-A(BZD)、NMDA、"
                 "nAChR α4β2、5-HT3、Kv 等其余 7 个", None, "是", "是",
     "无计算终点。CiPA 框架下 Cav1.2 与 Nav1.5 对整合致心律失常风险判断至关重要，"
     "而 Bowes-44 本身仅以通用型结合测定覆盖，SafetyScreen87 才细化到分子亚型", "none"),

    ("酶（非激酶）", "COX-1、COX-2、乙酰胆碱酯酶、MAO-A、PDE3A、PDE4D2", None,
     "是", "是",
     "无计算终点。COX-1（消化道出血）、MAO-A（高血压危象）、PDE3A（心衰死亡率）、"
     "PDE4D2（恶心呕吐）均为有明确临床后果的靶点", "none"),

    ("转运体", "DAT、NET、SERT", None, "是", "是",
     "无计算终点。ADMET-AI 的 P-糖蛋白终点<b>不属于</b>本面板（见下）", "none"),

    ("激酶", "Lck（Bowes-44）；Safety-77 扩展至 20 个激酶", None, "是", "是",
     "无计算终点。这是 Safety-77 相对 Bowes-44 最大的扩充方向，"
     "而计算覆盖度为零", "none"),
]

# endpoints that are NOT part of any of the three panels
NON_PANEL_SPEC = [
    ("CYP1A2 / 2C9 / 2C19 / 2D6 / 3A4 抑制",
     "CYP3A4_Veith",
     "<b>不属于</b>三个面板中的任何一个。CYP 抑制属于 DMPK / 药物相互作用 (DDI) 学科，"
     "依 FDA/EMA 的 DDI 指导原则执行，在 CRO 的业务线上也与安全药理学分属不同部门"),

    ("P-糖蛋白 (ABCB1) 抑制",
     "Pgp_Broccatelli",
     "<b>不属于</b>三个面板。Bowes-44 的 3 个转运体仅为 DAT/NET/SERT；"
     "SafetyScreen87 增加的是 Na⁺/K⁺-ATP 酶、ENT1、GAT-1，均非 ABC 外排转运体。"
     "P-gp 属于 DMPK 吸收分布与 DDI 终点"),

    ("芳香化酶 (CYP19A1)",
     "NR-Aromatase",
     "<b>不在</b>任何一个面板中。属于 Tox21 内分泌干扰筛查终点"),

    ("芳香烃受体 (AhR)",
     "NR-AhR",
     "<b>不在</b>任何一个面板中。AhR 为 bHLH-PAS 类转录因子，并非经典核受体，"
     "属于 Tox21 / CALUX 毒理学终点"),

    ("BSEP (ABCB11) 抑制",
     None,
     "<b>不在</b>任何一个面板中，且 ADMET-AI 亦<b>无对应终点</b>。"
     "胆汁盐外排泵抑制与胆汁淤积型肝损伤相关，需专门的肝转运体测定"
     "（BSEP 转染的倒置膜囊泡或三明治培养肝细胞）"),
]


REFERENCES = [
    "1.　Bowes J, Brown AJ, Hamon J, Jarolimek W, Sridhar A, Waldron G, Whitebread S. "
    "Reducing safety-related drug attrition: the use of in vitro pharmacological profiling. "
    "<i>Nature Reviews Drug Discovery</i> 2012; 11(12): 909–922. DOI: 10.1038/nrd3845.",

    "2.　Brennan RJ, Jenkinson S, Brown AJ, Delaunois A, Dumotier B, Pannirselvam M, "
    "Rao M, Ribeiro LR, Schmidt F, Sibony A, Timsit Y, Sales VT, Armstrong D, Lagrutta A, "
    "Mittlestadt SW, Naven R, Peri R, Roberts S, Vergis JM, Valentin J-P. "
    "The state of the art in secondary pharmacology and its impact on the safety of new medicines. "
    "<i>Nature Reviews Drug Discovery</i> 2024; 23(7): 525–545. DOI: 10.1038/s41573-024-00942-3.",

    "3.　Maciag M, Karamyan VT. Enzymes in secondary pharmacology screening panels: "
    "is there room for improvement? <i>Nature Reviews Drug Discovery</i> 2025; 24(6): 480–481. "
    "DOI: 10.1038/s41573-025-01173-w.",

    "4.　Papoian T, Chiu H-J, Elayan I, et al. Secondary pharmacology data to assess "
    "potential off-target activity of new drugs: a regulatory perspective. "
    "<i>Nature Reviews Drug Discovery</i> 2015; 14: 294. DOI: 10.1038/nrd3845-c1.",

    "5.　Eurofins Pharma Discovery Services. SafetyScreen44™ Panel product flyer. "
    "Ref. P270, Lit. No. EPDSFL420JUNE16, © Eurofins Cerep S.A. 2016. "
    "（本报告表 9 中 Bowes-44 的逐条靶点依据此实现版本）",

    "6.　Eurofins Discovery. SafetyScreen87 Panel 产品目录页 (P342 / PP223)。",

    "7.　ICH Harmonised Tripartite Guideline S7A. Safety Pharmacology Studies for "
    "Human Pharmaceuticals. ICH Step 4, 8 November 2000; EMA CPMP/ICH/539/00.",

    "8.　ICH Harmonised Tripartite Guideline S7B. The Non-Clinical Evaluation of the "
    "Potential for Delayed Ventricular Repolarization (QT Interval Prolongation) by "
    "Human Pharmaceuticals. ICH Step 4, May 2005; EMA CHMP/ICH/423/02.",

    "9.　ICH E14/S7B Implementation Working Group. Clinical and Nonclinical Evaluation "
    "of QT/QTc Interval Prolongation and Proarrhythmic Potential — Questions and Answers. "
    "Step 4, 21 February 2022. (Federal Register 87 FR 52716, 29 August 2022)",

    "10.　Kenna JG, et al. Can Bile Salt Export Pump Inhibition Testing in Drug Discovery "
    "and Development Reduce Liver Injury Risk? An International Transporter Consortium "
    "Perspective. <i>Clinical Pharmacology &amp; Therapeutics</i> 2018. DOI: 10.1002/cpt.1222.",

    "11.　Swanson K, Walther P, Leitz J, Mukherjee S, Wu JC, Shivnaraine RV, Zou J. "
    "ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical "
    "libraries. <i>Bioinformatics</i> 2024; 40(7): btae416. DOI: 10.1093/bioinformatics/btae416.",

    "12.　Huang K, Fu T, Gao W, et al. Therapeutics Data Commons: machine learning "
    "datasets and tasks for drug discovery and development. "
    "<i>NeurIPS Datasets and Benchmarks</i> 2021.",

    "13.　Veith H, Southall N, Huang R, et al. Comprehensive characterization of "
    "cytochrome P450 isozyme selectivity across chemical libraries. "
    "<i>Nature Biotechnology</i> 2009; 27(11): 1050–1055.",

    "14.　Broccatelli F, Carosati E, Neri A, et al. A novel approach for predicting "
    "P-glycoprotein (ABCB1) inhibition using molecular interaction fields. "
    "<i>Journal of Medicinal Chemistry</i> 2011; 54(6): 1740–1751.",

    "15.　Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. "
    "Quantifying the chemical beauty of drugs. <i>Nature Chemistry</i> 2012; 4: 90–98.",

    "16.　Baell JB, Holloway GA. New substructure filters for removal of pan assay "
    "interference compounds (PAINS) from screening libraries. "
    "<i>Journal of Medicinal Chemistry</i> 2010; 53(7): 2719–2740.",

    "17.　Comprehensive In Vitro Proarrhythmia Assay (CiPA) initiative. cipaproject.org.",
]


# --------------------------------------------------------------------------- #
def section_safety_panel(S, CW, FIG):
    el = [K.P("6　Safety Panel 对比与计算覆盖度分析", S["h1"])]

    el.append(K.P(
        "前面各章评估的是<b>计算能做什么</b>。本章回答互补的问题：相对于药物安全性评价"
        "实际使用的体外次级药理靶点面板，这些计算终点覆盖了多少？为此先比较业内三个"
        "参照面板，再把本次运行得到的 41 个 ADMET-AI 终点逐一对到面板靶点上。", S["body"]))

    # ---------------- 6.1 three panels ----------------
    el.append(K.P("6.1　三个参照面板", S["h2"]))
    el.append(K.P(
        "<b>Safety-44（Bowes 2012）</b>　由 AstraZeneca、GSK、Novartis 与 Pfizer 四家公司"
        "在 <i>Nature Reviews Drug Discovery</i> 上首次公开各自的体外药理谱系策略，并凝练出"
        "一个「最小推荐面板」，共 44 个靶点。它此后成为行业事实标准，本报告中称为 Safety-44。",
        S["body"]))
    el.append(K.P(
        "<b>Safety-77（Brennan / IQ DruSafe 2024）</b>　国际制药创新与质量联盟 (IQ) 的 "
        "DruSafe 工作组对 <b>18 家公司</b>做了问卷调查，总结当代实践并给出扩展的推荐靶点集，"
        "共 77 个靶点。其最显著的变化是激酶从 Bowes-44 的 1 个（2.3%）增加到 20 个（26%），"
        "跃居第二大靶点类别，同时明确把靶器官覆盖扩展到心血管与中枢神经系统之外。", S["body"]))
    el.append(K.P(
        "<b>SafetyScreen87（Eurofins）</b>　CRO 商业产品，在 Bowes-44 的全部 44 个靶点"
        "基础上增加 43 个（以 GPCR 与离子通道为主）。需要特别澄清：它是 Bowes-44 的"
        "<b>商业扩展版</b>，<b>早于</b>且独立于 Safety-77，二者并非同一体系，"
        "87 也不是 77 的扩充。", S["body"]))

    el.append(K.P("<b>表 8</b>　Safety-44 / Safety-77 / SafetyScreen87 面板对比", S["tcap"]))
    el.append(K.make_table(COMPARISON_HDR, COMPARISON_ROWS,
                           [22*mm, (CW-22*mm)/3, (CW-22*mm)/3, (CW-22*mm)/3],
                           S, font_size=6.9))
    el.append(K.P(
        "资料来源：Bowes 2012 与 Brennan 2024 原文均为付费墙文献，本表中 Bowes-44 的"
        "逐条靶点构成依据 Eurofins SafetyScreen44 产品资料（厂商明确标注其依据 Bowes "
        "et al. 2012），Safety-77 的激酶计数依据 Maciag &amp; Karamyan 2025 开放获取评论，"
        "SafetyScreen87 构成依据 Eurofins 产品目录页。<b>Safety-77 的完整靶点清单未能获取，"
        "本报告不对其余家族构成作任何推断。</b>", S["note"]))

    # ---------------- 6.2 regulatory ----------------
    el.append(K.P("6.2　与 ICH 法规框架的关系", S["h2"]))
    el.append(K.P(
        "一个常见误解是次级药理面板属于法规强制要求。实际并非如此。"
        "<b>ICH S7A</b> 强制要求的是<b>安全药理学核心组合试验</b>（心血管、呼吸、"
        "中枢神经三个系统的功能学评价），必须在首次人体给药前完成并通常需符合 GLP。"
        "次级药理面板在 S7A 中仅出现在第 2.2 节：「提示潜在不良作用的配体结合或酶测定数据」"
        "被列为设计核心组合试验时<b>应予考虑</b>的四项因素之一；S7A 第 2.11 节并明确指出"
        "次级药效学研究<b>通常不需要</b>符合 GLP。", S["body"]))
    el.append(K.P(
        "<b>ICH S7B</b> 单独处理延迟心室复极（QT 间期延长）风险，要求体外 IKr (hERG) "
        "测定加体内 QT 研究的整合风险评估；2022 年的 <b>ICH E14/S7B 问答文件</b>进一步"
        "提出体外 hERG 测定的最佳实践，使高质量的非临床数据能实质性参与整合风险评估。"
        "CiPA 计划（多离子通道 + 计算机心肌动作电位重建 + iPSC 心肌细胞）推动了这一转变，"
        "但 CiPA 本身<b>不是</b>法规要求。", S["body"]))
    el.append(K.callout(
        "因此准确的表述是：次级药理面板是<b>行业共识的推荐实践</b>，而非法规强制项目。"
        "其约束力来自企业间趋同与审评实践——正如 Brennan 等 2024 年论文摘要所述，"
        "监管机构「正日益要求」提供针对已知不良反应关联靶点的活性数据。"
        "CRO 宣传中的「ICH S7A/B 对齐」指所筛靶点对应 S7A 关注的器官系统，"
        "而非该面板满足了任何 S7A/S7B 要求。", S))

    el.append(K.P("6.3　计算覆盖度映射", S["h2"]))
    el.append(K.P(
        "下表把 Bowes-44 的六个靶点家族逐一与本次运行得到的 ADMET-AI 终点对照。"
        "「可映射」指存在名义对应的计算终点；「已验证」指该终点还通过了第 4 章的"
        "阳性对照检验。表中的阳性药物示例全部取自本次 30 药物面板的实际预测结果"
        "（阈值 0.5），未引用外部数据。", S["body"]))

    stats = MC.endpoint_stats()
    val = MC.validation_lookup()
    hdr = ["靶点家族", "面板靶点", "对应 ADMET-AI 终点", "在 44 / 87 中",
           "本面板阳性药物示例", "判定"]
    rows, hl = [], {}
    for i, (fam, targets, ep, in44, in87, note, kind) in enumerate(COVERAGE_SPEC):
        if ep:
            epname = ep
            pos = stats[ep]["positives"]
            v = val.get(ep)
            vtxt = f"（对照 {v[0]}/{v[1]}）" if v else ""
            pos = pos + vtxt
        else:
            epname = "—（无对应终点）"
            pos = "—"
        verdict = {"covered": "<b>可映射且已验证</b>",
                   "nominal": "仅名义映射，验证失败",
                   "nominal87": "仅名义映射，验证失败",
                   "none": "<b>需体外实验</b>"}[kind]
        rows.append([fam, targets, epname, f"{in44} / {in87}", pos,
                     verdict + "<br/><font size='6'>" + note + "</font>"])
        hl[(i, 5)] = K.GOOD if kind == "covered" else K.WARN
    el.append(K.P("<b>表 9</b>　Safety panel 靶点与 ADMET-AI 计算终点的覆盖度映射", S["tcap"]))
    el.append(K.make_table(hdr, rows,
                           [16*mm, 40*mm, 24*mm, 15*mm, 34*mm, CW-129*mm],
                           S, font_size=6.6, highlight=hl))

    el.append(K.P("6.4　被误认为属于面板的计算终点", S["h2"]))
    el.append(K.P(
        "在做这类映射时最容易出现的错误，是把 DMPK 终点当作安全药理面板靶点。"
        "下表列出四个常被误归入面板的 ADMET-AI 终点，以及一个面板与计算<b>双双缺失</b>"
        "的重要靶点。", S["body"]))
    rows2 = []
    for label, ep, note in NON_PANEL_SPEC:
        if ep:
            pos = stats[ep]["positives"]
        else:
            pos = "ADMET-AI 无此终点"
        rows2.append([label, pos, note])
    el.append(K.P("<b>表 10</b>　不属于次级药理面板的 ADMET 终点", S["tcap"]))
    el.append(K.make_table(["ADMET-AI 终点", "本面板阳性药物示例", "归属说明"],
                           rows2, [34*mm, 44*mm, CW-78*mm], S, font_size=6.9))
    el.append(K.P(
        "需要说明的是，这些终点本身是<b>有价值的</b>——CYP3A4 与 P-糖蛋白终点在第 4.1 节的"
        "对照验证中表现优异，酮康唑与维拉帕米均被正确识别。问题仅在于它们回答的是"
        "药物相互作用问题，而不是脱靶安全性问题，因此不能计入安全面板的覆盖度。", S["body"]))

    # ---------------- 6.5 quantitative conclusion ----------------
    el.append(K.P("6.5　覆盖度的定量结论", S["h2"]))
    el.append(K.figure(f"{FIG}/fig6_panel_coverage.png", CW * 0.88, S,
                       "<b>图 6</b>　Bowes-44 安全药理面板按靶点家族的计算覆盖度。"
                       "灰色为面板靶点数，绿色为存在名义对应 ADMET-AI 终点的靶点数。"))
    el.append(K.callout(
        f"在 Bowes-44 的 <b>{TOTAL_44}</b> 个靶点中，存在名义对应计算终点的仅 "
        f"<b>{TOTAL_COV}</b> 个（hERG 与雄激素受体），占 "
        f"<b>{100*TOTAL_COV/TOTAL_44:.1f}%</b>；而其中还能通过本报告第 4 章阳性对照验证的"
        f"只有 <b>{TOTAL_VAL}</b> 个（hERG），占 <b>{100*TOTAL_VAL/TOTAL_44:.1f}%</b>。"
        f"若以 SafetyScreen87 为基准，名义可映射靶点为 4 个（hERG、AR、ERα、PPARγ），"
        f"约 4.6%，经验证可用的仍只有 hERG 一个。", S,
        colour=K.colors.HexColor("#fdf0ee"), border=K.WARN))
    el.append(K.P(
        "这一数字远低于直觉预期，原因在于两类学科边界。其一，占 Bowes-44 逾半数的 "
        "<b>24 个 GPCR</b>、以及全部 6 个非激酶酶类、3 个单胺转运体与激酶，"
        "在 TDC / ADMET-AI 的终点集中<b>完全没有对应模型</b>——公开的大规模构效数据"
        "集中在 ADME 与通用毒性终点上，而非逐个受体的结合亲和力。其二，几个看似对应的"
        "终点（CYP、P-糖蛋白）实际归属 DMPK 学科，不在面板范围内。", S["body"]))
    el.append(K.P(
        "更值得注意的是<b>方向性问题</b>。即便将来出现覆盖 GPCR 的计算模型，"
        "安全药理最关键的判断往往不是「是否结合」而是「激动还是拮抗」——"
        "5-HT2B 导致心脏瓣膜病的是<b>激动</b>作用，这一点连以结合测定为主的 "
        "Bowes-44 本身都无法区分，这也正是 Bowes 等人建议结合测定与功能测定并行、"
        "以及 Safety-77 强调功能性剂量-反应的原因。", S["body"]))

    el.append(K.P("6.6　互补而非替代", S["h2"]))
    el.append(K.P(
        "综合第 4 章的验证结果与本章的映射结果，计算方法与体外面板的分工是清晰的：",
        S["body"]))
    el += K.bullets([
        "<b>计算方法的优势区间</b>是 DMPK 与心脏钾通道：hERG、CYP 抑制、P-糖蛋白、"
        "血脑屏障、溶解度与渗透性。这些终点在化合物<b>合成之前</b>即可给出排序，"
        "成本近乎为零，适合在先导优化阶段筛除明显有问题的系列。",

        "<b>体外面板不可替代的部分</b>是广谱受体脱靶：24 个 GPCR、其余 7 个离子通道、"
        "6 个酶、3 个转运体与激酶家族，合计约占 Bowes-44 的 95%。这些必须依赖实验，"
        "且对 5-HT2B 等靶点还需功能性测定以判断激动/拮抗方向。",

        "<b>两者的时间顺序也不同</b>：计算置于化合物合成之前用于排序；体外面板置于"
        "候选化合物提名前后用于系统性排查。把面板结果反馈给计算模型（作为新的训练数据）"
        "是缩小这一差距的现实路径，但在公开数据层面尚未实现。",

        "<b>本报告第 4 章提供了一个可复用的检验方法</b>：在把任何计算终点用于安全决策之前，"
        "先用该终点已知的阳性/阴性对照药物检验它能否重现既定事实。本次检验使 hERG "
        "与四个核受体终点得到了截然不同的结论，而这一差别单看模型发表的 AUC 指标是"
        "看不出来的。",
    ], S["body"], S["bullet"])
    return el
