"""Chinese-language campaign report.

Mirrors `report.py` in structure and claims. Two typographic differences:
the CJK face has a single weight, so hierarchy comes from size and colour
rather than bold; and leading is opened up, since CJK text needs more line
spacing than Latin at the same point size.
"""

from __future__ import annotations

import json
import pathlib
from datetime import date

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                PageTemplate, Paragraph, Spacer, Table,
                                TableStyle)

from . import cjk
from .report import (ACCENT, HAIR, INK, INK2, MUTED, RULE, SURFACE, WARN,
                     MARGIN, PAGE_H, PAGE_W)

F = "CJK"


def styles():
    """Paragraph styles. wordWrap="CJK" is essential: without it reportlab
    breaks only at spaces, and a Chinese paragraph becomes a handful of giant
    "words" whose inter-character spacing gets stretched grotesquely by
    justification."""
    ss = getSampleStyleSheet()
    s = {}
    s["title"] = ParagraphStyle("t", parent=ss["Title"], fontName=F, fontSize=20,
                                leading=27, textColor=INK, alignment=TA_LEFT,
                                spaceAfter=2,
                               wordWrap="CJK")
    s["subtitle"] = ParagraphStyle("st", parent=ss["Normal"], fontName=F,
                                   fontSize=10, leading=15.5, textColor=INK2,
                                   spaceAfter=2,
                               wordWrap="CJK")
    s["slot"] = ParagraphStyle("sl", parent=ss["Normal"], fontName=F,
                               fontSize=7.5, leading=11, textColor=MUTED,
                               wordWrap="CJK")
    s["h1"] = ParagraphStyle("h1", parent=ss["Heading1"], fontName=F,
                             fontSize=14.5, leading=20, textColor=ACCENT,
                             spaceBefore=15, spaceAfter=7,
                               wordWrap="CJK")
    s["h2"] = ParagraphStyle("h2", parent=ss["Heading2"], fontName=F,
                             fontSize=11, leading=16, textColor=INK,
                             spaceBefore=11, spaceAfter=5,
                               wordWrap="CJK")
    s["body"] = ParagraphStyle("b", parent=ss["BodyText"], fontName=F,
                               fontSize=9.6, leading=16, textColor=INK,
                               alignment=TA_JUSTIFY, spaceAfter=7,
                               wordWrap="CJK")
    s["lead"] = ParagraphStyle("ld", parent=s["body"], textColor=ACCENT,
                               fontSize=9.8, spaceAfter=2,
                               wordWrap="CJK")
    s["small"] = ParagraphStyle("sm", parent=ss["BodyText"], fontName=F,
                                fontSize=8.4, leading=13.6, textColor=INK2,
                                alignment=TA_JUSTIFY, spaceAfter=6,
                               wordWrap="CJK")
    s["caption"] = ParagraphStyle("cap", parent=ss["BodyText"], fontName=F,
                                  fontSize=8.2, leading=12.5, textColor=INK2,
                                  spaceBefore=4, spaceAfter=11,
                               wordWrap="CJK")
    s["ref"] = ParagraphStyle("rf", parent=ss["BodyText"], fontName=F,
                              fontSize=8.1, leading=12.6, textColor=INK,
                              alignment=TA_LEFT, spaceAfter=5,
                              leftIndent=15, firstLineIndent=-15,
                               wordWrap="CJK")
    s["cell"] = ParagraphStyle("c", parent=ss["BodyText"], fontName=F,
                               fontSize=7.1, leading=9.4, textColor=INK,
                               alignment=TA_LEFT, spaceAfter=0,
                               wordWrap="CJK")
    s["cellh"] = ParagraphStyle("ch", parent=s["cell"], textColor=ACCENT,
                               wordWrap="CJK")
    s["callout"] = ParagraphStyle("co", parent=ss["BodyText"], fontName=F,
                                  fontSize=9, leading=15, textColor=INK,
                                  alignment=TA_JUSTIFY, leftIndent=8,
                                  rightIndent=8, spaceBefore=5, spaceAfter=9,
                                  borderPadding=8,
                                  backColor=colors.HexColor("#fdf2ec"),
                               wordWrap="CJK")
    return s


def _footer(canvas, doc):
    canvas.saveState()
    canvas.setFillColor(SURFACE)
    canvas.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    canvas.setStrokeColor(HAIR)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 13 * mm, PAGE_W - MARGIN, 13 * mm)
    canvas.setFont(F, 7.4)
    canvas.setFillColor(MUTED)
    canvas.drawString(MARGIN, 9 * mm, "GSK3β 从头设计项目 — 全部为计算设计分子，未经实验验证")
    canvas.drawRightString(PAGE_W - MARGIN, 9 * mm, f"第 {doc.page} 页")
    canvas.restoreState()


def stat_row(items, width):
    cells = []
    for value, label in items:
        cells.append(Paragraph(
            f'<font size="15" color="#2a78d6">{value}</font><br/>'
            f'<font size="7.4" color="#52514e">{label}</font>',
            ParagraphStyle("s", fontName=F, leading=15)))
    t = Table([cells], colWidths=[width / len(cells)] * len(cells))
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 9),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 9),
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, RULE),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
    ]))
    return t


def fig(path, width, caption, s):
    from PIL import Image as PILImage
    iw, ih = PILImage.open(path).size
    return KeepTogether([Image(str(path), width=width, height=width * ih / iw),
                         Paragraph(caption, s["caption"])])


def build(result_dir: str, out_pdf: str):
    cjk.register_reportlab(name=F)
    R = pathlib.Path(result_dir)
    figs, figs_zh = R / "figures", R / "figures_zh"
    s = styles()
    D = json.load(open(R / "top10_designs.json"))
    log = json.load(open(R / "cascade_log.json"))
    params = json.load(open(R / "ga_params.json"))
    hist = json.load(open(R / "ga_history.json"))
    retro = {r["design_id"]: r for r in json.load(open(R / "retro_results.json"))}
    meta = json.load(open(R / "run_meta.json"))

    avail = PAGE_W - 2 * MARGIN
    doc = BaseDocTemplate(out_pdf, pagesize=A4, leftMargin=MARGIN,
                          rightMargin=MARGIN, topMargin=MARGIN,
                          bottomMargin=20 * mm,
                          title="GSK3β 从头小分子设计报告",
                          author="计算设计报告")
    doc.addPageTemplates([PageTemplate(id="all", frames=[
        Frame(MARGIN, 20 * mm, avail, PAGE_H - MARGIN - 20 * mm, id="m")],
        onPage=_footer)])

    E = []
    P = lambda t, st="body": Paragraph(t, s[st])

    slot = Table([[Paragraph("［品牌标识位］<br/>logo / 字标", s["slot"])]],
                 colWidths=[38 * mm], rowHeights=[13 * mm])
    slot.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, HAIR),
                              ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                              ("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    head = Table([[Paragraph("GSK3β 从头小分子设计报告", s["title"]), slot]],
                 colWidths=[avail - 40 * mm, 40 * mm])
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                              ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    E += [head]
    E += [P("靶点：糖原合成酶激酶-3β（GSK3B，UniProt P49841，ChEMBL CHEMBL262）", "subtitle")]
    E += [P(f"生成日期 {date.today().isoformat()} ｜ 图遗传算法 ｜ 共评估 "
            f"{meta['n_unique']:,} 个分子 ｜ 最终入选 10 个设计", "subtitle")]
    E += [Spacer(1, 8)]

    E += [stat_row([
        (f"{meta['n_unique']:,}", "生成并打分的<br/>去重分子总数"),
        (f"{log[0]['out']:,}", "合法且通过<br/>性质门"),
        (f"{log[-1]['out']:,}", "通过完整<br/>筛选级联"),
        ("10", "最终入选<br/>（互不相同的骨架类型）"),
        (f"{meta['n_retro_solved']}/10", "获得可行<br/>逆合成路线"),
    ], avail)]
    E += [Spacer(1, 5)]

    E += [Paragraph(
        "<font color='#eb6834'>本报告是什么。</font>"
        "十个由算法从头生成的分子，依据代理活性模型结合类药性与合成可及性启发式指标排序。"
        "<font color='#eb6834'>这些化合物均未经合成，也未做任何实验测定。</font>"
        "报告中出现的每一个活性数值都是模型预测值而非实测值，"
        "而该模型恰恰就是生成算法所优化的对象——在依据任何排序做决策之前，请先阅读「局限性」一节。",
        s["callout"])]

    # 1 引言
    E += [P("1. 引言", "h1")]
    E += [P(
        "GSK3B 是一种组成型激活的丝氨酸/苏氨酸激酶，位于 Wnt/β-catenin 通路、胰岛素信号通路"
        "与微管调控的交汇点。它的失调与阿尔茨海默病的几项核心病理特征直接相关——磷酸化 tau 蛋白、"
        "调节 β-淀粉样蛋白生成、影响神经发生与突触功能——因此长期以来是中枢神经系统领域被反复"
        "开发的重要靶点[6]。这一历史对设计工作而言是双刃剑：靶点验证充分、数据注释丰富，"
        "数据驱动的生成模型有足够的学习素材；但围绕它的可及化学空间也早已被药物化学家充分开垦，"
        "真正新颖的骨架类型相应地更难触及。")]
    E += [P(
        "本项目要回答的是一个范围明确、可操作的问题：以已知的高活性 GSK3B 结合分子为起点，"
        "图遗传算法能否产生同时满足以下三点的分子——（a）在预训练的 GSK3B 活性代理模型上得分良好；"
        "（b）具备类药性且合成上可行；（c）在结构上与起点分子有实质区别？"
        "交付物是一份含性质、新颖性度量、结构警示与逆合成路线的十分子候选清单，"
        "供药物化学家做初步分诊之用，而非一组结论。")]

    # 2 方法
    E += [P("2. 方法", "h1")]
    E += [P("2.1 参考数据与种子分子", "h2")]
    E += [P(
        f"已知活性分子通过 ChEMBL REST API 获取（版本 <font color='#2a78d6'>{meta['chembl_release']}</font>，"
        f"检索日期 {meta['retrieved']}），靶点 CHEMBL262，限定为 IC50/Ki/Kd 类型且 pChEMBL ≥ 7 的记录"
        f"（即活性优于 100 nM）。共返回 {meta['n_activity_records']:,} 条活性记录，"
        f"对应 <font color='#2a78d6'>{meta['n_reference_actives']:,} 个去重化合物</font>，"
        f"每个化合物取其最高活性记录。该集合承担两个彼此独立的角色：遗传算法初始种群的来源，"
        f"以及后续衡量新颖性的参照系。初始种群取其中结构差异最大的 {params['pop_size']} 个分子，"
        f"由 Morgan 指纹上的 MaxMin 算法挑选，以避免算法一开始就落在单一骨架类型内部。")]

    E += [P("2.2 分子生成", "h2")]
    E += [P(
        f"分子由直接作用于分子图的图遗传算法生成，算子族遵循 Jensen [1] 与 GuacaMol 的 "
        f"graph-GA 基线 [2]，本次基于 RDKit 自行实现，未沿用上述任一代码库。"
        f"交叉操作在两个亲本的随机非环单键处断裂，各取一个片段重新连接；"
        f"变异操作每次施加一种局部图编辑：原子追加（取自常见药物化学环系与取代基库）、"
        f"元素替换、端位原子删除、键级改变、或在键中插入原子。"
        f"每个子代都经过净化并通过 SMILES 往返校验，凡 RDKit 无法解析者一律丢弃。"
        f"选择策略为精英截断，亲本按排名加权采样。参数：种群 {params['pop_size']}，"
        f"迭代 {params['n_generations']} 代，每代产生 {params['offspring_per_gen']} 个子代，"
        f"变异率 {params['mutation_rate']}，随机种子 {params['seed']}。"
        f"CPU 实际运行时间 {params['runtime_s']:.0f} 秒。")]

    E += [P("2.3 目标函数", "h2")]
    E += [P(
        "每个分子的得分为三项归一化指标的几何平均：GSK3B 代理活性、QED 类药性 [3]，"
        "以及合成可及性项 (10 − SA_Score)/9 [4]。采用几何平均是刻意的合取式设计："
        "任一项为零则总分为零，算法无法通过牺牲类药性来换取更高的预测活性。"
        "活性项采用 Therapeutics Data Commons [8] 发布的预训练 GSK3B 代理模型——"
        "一个基于 ECFP 特征的随机森林，输出值域 [0,1]。"
        "此外设置硬性门控，将落在以下范围之外的分子直接归零：分子量 150–600、重原子数 10–50、"
        "cLogP ≤ 6、SA_Score ≤ 6、形式电荷为零、单一片段、无大于七元的环，"
        "元素限定为 C/N/O/S/F/Cl/Br。若无此门控，优化代理模型的遗传算法必然会漂移到"
        "得分虚高但根本无法合成的化学空间中去。")]

    E += [P("2.4 筛选与遴选", "h2")]
    E += [P(
        "候选池取自运行过程中出现过的全部分子，而非仅最终种群。级联筛选依次为："
        "合法性、去重与门控存活；新颖性（对全部参考活性分子的最大 Tanimoto < 0.40，"
        "Morgan 半径 2、2048 位）；活性 ≥ 0.50 且 QED ≥ 0.60；不含 PAINS A/B/C 子结构 [5]；"
        "环系合理性；SA_Score ≤ 4.5。"
        "直接取得分最高的十个分子，结果是同一母核的十种修饰，因此最终遴选改为："
        "对存活分子做聚类（Butina，Tanimoto 0.4），从得分最高的十个簇中各取最优代表。"
        "此前曾尝试以 Bemis-Murcko 骨架唯一性作为判据并予以否决——"
        "更换一个侧挂环即可改变 Murcko 骨架，而分子识别基元丝毫未变，该判据过于宽松。")]

    E += [P("2.5 逆合成与结构警示", "h2")]
    E += [P(
        "逆合成路线由 AiZynthFinder [7] 给出，使用公开的 USPTO 模板库、扩展与过滤策略"
        "以及 ZINC 可购库存，每个分子的搜索预算为 120 秒。"
        "此外对入选分子做了 Brenk 活性/不宜基团警示与迈克尔受体 SMARTS 的显式比对；"
        "这些结果以标注形式记录而非用作筛选条件，原因是若干已验证的 GSK3B 骨架类型本身即为马来酰亚胺。")]

    # 3 结果
    E += [P("3. 结果", "h1")]
    E += [P("3.1 收敛情况", "h2")]
    first, last = hist[0], hist[-1]
    E += [P(
        f"算法收敛平稳。{params['n_generations']} 代中，种群平均分从 {first['mean']:.3f} "
        f"升至 {last['mean']:.3f}，最优个体从 {first['best']:.3f} 升至 {last['best']:.3f}；"
        f"在同时满足 Tanimoto < 0.4 新颖性判据的分子中，最优分从 {first['best_novel']:.3f} "
        f"升至 {last['best_novel']:.3f}。最优分曲线呈阶梯状，这是精英截断策略的预期表现。"
        f"平均分曲线在约第 25 代后趋平，说明种群已基本收敛；未继续延长迭代，"
        f"因为后续世代产出的多是已有骨架类型的修饰，而非新的骨架。")]
    E += [fig(figs_zh / "fig1_convergence.png", avail,
              "<font color='#2a78d6'>图 1.</font> 各代目标函数得分。"
              "「新颖分子最优」曲线始终低于不受约束的最优值，这正是本项目的核心张力："
              "代理模型奖励与已知活性分子的相似性，而新颖性判据恰恰惩罚这种相似性。", s)]

    E += [P("3.2 筛选衰减", "h2")]
    E += [P(
        f"在 {log[0]['in']:,} 个去重分子中，{log[0]['out']:,} 个合法、去重且落在性质门之内。"
        f"新颖性是迄今最严苛的一道筛选，剔除了 {log[1]['removed']:,} 个分子"
        f"（占进入该阶段分子的 {100*log[1]['removed']/log[1]['in']:.0f}%）——"
        f"算法的大部分努力都花在了起点附近。活性与 QED 联合门控又剔除 {log[2]['removed']:,} 个。"
        f"PAINS、环系合理性与合成可及性三项合计仅剔除 "
        f"{log[3]['removed']+log[4]['removed']+log[5]['removed']} 个，"
        f"这与其说反映了设计分子质量优良，不如说是因为性质门已经提前排除了这些筛选所针对的大部分情形。"
        f"最终 {log[-1]['out']} 个分子通过完整级联。")]
    E += [fig(figs_zh / "fig3_cascade.png", avail,
              "<font color='#2a78d6'>图 2.</font> 按图示顺序施加各级筛选后的剩余分子数。", s)]

    E += [P("3.3 新颖性", "h2")]
    nd = meta["novelty_bands"]
    E += [P(
        f"此处的新颖性是一个可量化的指标，而测量结果并不乐观。"
        f"在全部生成分子中，与最近已知活性分子的相似度众数约在 0.45–0.50。"
        f"在通过级联的 {log[-1]['out']} 个分子中，分布被紧紧压在 0.40 这条界线之下："
        f"{nd['0.35-0.40']} 个落在 0.35–0.40 之间，{nd['0.30-0.35']} 个落在 0.30–0.35，"
        f"{nd['0.20-0.30']} 个落在 0.20–0.30，"
        f"<font color='#eb6834'>0.20 以下一个也没有</font>。"
        f"最终入选的十个设计相似度区间为 {meta['sim_min']:.2f}–{meta['sim_max']:.2f}。"
        f"因此这些分子的「新颖」仅在方法一节所定义的、相对于特定阈值的意义上成立。"
        f"它们不是新的骨架类型，而是一片已被充分开垦的区域边缘上的新坐标点。"
        f"换一个阈值就会得到另一套答案，而这种敏感性本身就是结论的一部分。")]
    E += [fig(figs_zh / "fig5_novelty.png", avail,
              f"<font color='#2a78d6'>图 3.</font> 每个生成分子与 "
              f"{meta['n_reference_actives']:,} 个参考活性分子中最近者的最大 Tanimoto 相似度。"
              f"位于 1.0 处的尖峰是种子分子本身，它们原本就是已知活性分子。", s)]

    E += [P("3.4 理化性质", "h2")]
    E += [P(
        "存活的设计分子系统性地比其所源自的已知活性分子更小、极性更低——"
        "生成分子的分子量分布集中在 250–350 Da，而 ChEMBL 集合约在 350–450 Da，"
        "TPSA 也相应更低。这是含 QED 项的打分函数的已知偏倚："
        "QED 在各性质区间的中段取得峰值，而非在许多经优化的激酶抑制剂所处的上段。"
        "合成可及性与参考集合重合良好。其实际含义是：这些设计分子在后续优化中还有向上生长的余量，"
        "对一个起点而言，这个位置是合理的。")]
    E += [fig(figs_zh / "fig4_properties.png", avail,
              "<font color='#2a78d6'>图 4.</font> 通过级联的设计分子（填充）"
              "与随机抽取的 400 个 ChEMBL 参考活性分子（描边）的性质分布对比。", s)]

    E += [P("3.5 入选设计", "h2")]
    E += [P(
        f"十个设计分别取自得分最高的十个骨架类型簇（存活分子共聚为 {meta['n_clusters']} 簇）"
        f"的最优代表。相较于直接取前十名，骨架感知的遴选方式在平均目标分上付出了约 0.03 的代价"
        f"（{meta['mean_score_selected']:.3f} 对 {meta['mean_score_naive']:.3f}），"
        f"换来的是一组识别基元明显不同的候选组合，而不是同一个方案的十种写法。")]
    rows = [[Paragraph(h, s["cellh"]) for h in
             ["编号", "总分", "活性", "QED", "SA", "相似度", "分子量",
              "cLogP", "TPSA", "路线", "结构警示"]]]
    zh_alert = {"maleimide": "马来酰亚胺", "michael_acceptor": "迈克尔受体",
                "thiophene_dione": "噻吩二酮", "nitrile": "腈基",
                "sulfonamide": "磺酰胺", "sulfonate_ester": "磺酸酯"}
    for d in D:
        r = retro.get(d["design_id"], {})
        n_steps = r.get("n_steps_best")
        route = f"{n_steps} 步" if r.get("solved") else "未解出"
        alerts = "、".join(zh_alert.get(a, a) for a in d.get("reactive_motifs", [])) or "—"
        rows.append([Paragraph(x, s["cell"]) for x in [
            d["design_id"].replace("GSK3B-DN-", "DN-"),
            f"{d['score']:.3f}", f"{d['activity']:.2f}", f"{d['qed']:.2f}",
            f"{d['sa']:.2f}", f"{d['max_sim_to_known']:.2f}", f"{d['mw']:.0f}",
            f"{d['logp']:.1f}", f"{d['tpsa']:.0f}", route, alerts]])
    t = Table(rows, colWidths=[13 * mm, 11 * mm, 10 * mm, 10 * mm, 9 * mm,
                               12 * mm, 12 * mm, 12 * mm, 11 * mm, 13 * mm,
                               avail - 113 * mm], repeatRows=1)
    t.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.3, HAIR),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, RULE),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
    ]))
    E += [KeepTogether([t, Paragraph(
        "<font color='#2a78d6'>表 1.</font> 入选设计一览。"
        "「活性」为 GSK3B 代理模型输出；「相似度」为对任一参考活性分子的最大 Tanimoto；"
        "「路线」为 120 秒预算内 AiZynthFinder 找到的最短路线。各分子结构见图 5，"
        "完整 SMILES、骨架与最近邻参考分子见随附的 CSV 文件。", s["caption"])])]

    E += [fig(figs / "fig6_top10.png", avail * 0.82,
              "<font color='#2a78d6'>图 5.</font> 十个入选设计，每个骨架类型簇一个代表。"
              "对应的得分、性质与结构警示见表 1；它们都是代理模型输出与启发式指标，并非实测值。", s)]

    E += [P("3.6 逆合成", "h2")]
    solved = [d for d in D if retro.get(d["design_id"], {}).get("solved")]
    E += [P(
        f"AiZynthFinder 为十个设计中的 <font color='#2a78d6'>{len(solved)} 个</font>"
        f"返回了至少一条完全解出的路线（所有叶节点均在 ZINC 库存中），最短路线为 1–6 步。"
        f"其中若干提议在化学上相当常规、表面看来可行——DN-01 归结为芳基硼酸与溴代吲唑的 Suzuki 偶联，"
        f"DN-08 归结为甲硫基嘧啶与苯胺之间的芳香亲核取代。"
        f"未解出的三例（DN-02、DN-06、DN-07）返回的是部分路线树，其叶节点未能在预算内全部落到可购试剂上；"
        f"这说明的是模板库覆盖范围与搜索预算的限制，并不能证明这些分子不可合成。"
        f"所有路线均为基于模板的提议，不包含任何关于收率、选择性或保护基策略的信息。")]

    # 4 结论
    E += [P("4. 结论", "h1")]
    for lead, rest in [
        ("流程已完整跑通，且成本极低。",
         f"共生成 {meta['n_unique']:,} 个分子，完成三目标打分、筛选、聚类与逆合成，"
         f"CPU 耗时不足 {int(params['runtime_s']//60)+3} 分钟，随机种子固定，结果可复现。"),
        ("优化成功了，新颖性在很大程度上没有。",
         "目标函数提升显著，设计分子具备类药性且合成可行。但没有任何一个存活分子与已知活性分子的"
         "相似度低于 0.20，且存活分子紧贴 0.40 这条界线堆积。"
         "诚实的解读是：算法重新发现并重组了已有的 GSK3B 药效团家族——马来酰亚胺、苯胺基嘧啶、"
         "氮杂吲哚、吲唑——而不是创造了新的家族。"),
        ("这种「重新发现」本身是一个弱阳性对照。",
         "一个优化 GSK3B 代理模型的生成器，最终收敛到真实 GSK3B 项目同样收敛到的骨架类型上，"
         "这说明代理模型确实捕捉到了某些真实信息。但它不能作为这些具体分子具有活性的证据。"),
        ("骨架坍缩是默认结果，必须主动设计对策。",
         "朴素的前十名是同一母核的十种修饰。只有显式聚类才产生了真正的候选组合。"
         "任何不带多样性判据就报告 top-N 清单的设计项目，很可能是把同一个赌注报告了 N 遍。"),
    ]:
        E += [P(f"<font color='#2a78d6'>{lead}</font> {rest}")]

    # 5 局限
    E += [P("5. 局限性", "h1")]
    E += [Paragraph("以下局限并非例行免责声明；其中数条会实质性影响上述排序该如何被解读。",
                    s["callout"])]
    for head_txt, body in [
        ("活性评分是代理模型，而生成算法正是直接针对它做优化的。",
         "这是最主要的一条。GSK3B 代理模型是在公开数据上训练的、基于指纹的随机森林；"
         "遗传算法以数千次查询去搜索它，必然会找到部分因模型假象而得分偏高的输入。"
         "代理模型输出高，只能说明该分子与模型训练集中的活性分子相似，不能说明它会结合靶点。"
         "本次未做任何分子对接、自由能计算或结构验证。"),
        ("新颖性是相对于阈值的，且相当勉强。",
         "0.40 这个 Tanimoto 阈值是一个约定，而非化学本身的性质。存活分子集中在 0.35–0.40，"
         "因此阈值、指纹类型或半径的小幅改动，都会显著改变哪些分子算作「新颖」。"),
        ("十个设计中六个带有结构警示，其中五个是马来酰亚胺。",
         "马来酰亚胺以及 DN-07 中的噻吩二酮都是迈克尔受体，存在共价反应性与广谱干扰的风险。"
         "PAINS 过滤并不能识别它们。保留它们是因为已验证的 GSK3B 抑制剂"
         "（SB-216763、双吲哚基马来酰亚胺）本身就含这一基元——"
         "但一份以马来酰亚胺为主的候选清单是实实在在的风险，而非中性的观察结果；"
         "DN-01、DN-06、DN-08、DN-09 是其中不带警示的四个。"),
        ("QED 与 SA_Score 都是启发式指标。",
         "QED 编码的是一种聚合的类药性观念，会惩罚那些在经优化的激酶抑制剂中常见的较大、"
         "极性较强的分子，这一点在图 4 中清晰可见。SA_Score 是基于片段频率的替代指标，"
         "系统性地偏袒由常见片段拼搭而成的分子——而这恰恰就是本生成器所构造的那类分子。"),
        ("逆合成路线是 USPTO 模板给出的提议。",
         "路线「解出」只意味着搜索触及了可购买的叶节点，不蕴含任何关于收率、选择性、规模或保护基的信息，"
         "且模板库反映的是历史上已发表的反应。"),
        ("单次运行，单一随机种子。",
         "仅在种子 42 下运行了一次遗传算法，没有重复实验，因此结果的方差未经测量。"
         "这十个具体分子应被视为某个分布中的一次抽样，而非最优解。"),
        ("完全未评估选择性。",
         "GSK3A 与 GSK3B 的 ATP 位点结构几乎完全一致，而本流程没有涉及任何激酶组范围的选择性、"
         "细胞活性、膜通透性、血脑屏障穿透或毒性评估。"),
    ]:
        E += [Paragraph(f"<font color='#eb6834'>{head_txt}</font> {body}", s["small"])]

    # 6 后续
    E += [P("6. 后续工作建议", "h1")]
    for lead, rest in [
        ("打破代理模型的主导。",
         "用算法从未见过的正交方法重新给 248 个存活分子打分——"
         "对接到 GSK3B 的 ATP 位点结构（可用结构不少，如 1Q5K、6B8J），或采用共折叠模型——"
         "只保留两种方法都认可的分子。在没有实验的条件下，两种方法之间的分歧是信息量最大的信号。"),
        ("做重复实验并报告方差。",
         "跑 5–10 个随机种子，然后统计每个骨架类型复现的频率。"
         "跨种子稳定出现的骨架，价值高于单次运行中的最高分。"),
        ("把新颖性放进目标函数，而不只是放在筛选里。",
         "加入显式的新颖性或多样性奖励，可以让算法主动用预测活性去交换真正的结构突破，"
         "而不是在事后被筛掉——后者正是造成分子在阈值处堆积的原因。"),
        ("就马来酰亚胺做明确的取舍。",
         "在项目层面决定共价/反应性骨架是否在范围之内。若不在，可在门控中加入迈克尔受体排除项重跑；"
         "不带警示的四个设计（DN-01、DN-06、DN-08、DN-09）说明本流程有能力产出这类分子。"),
        ("补充 GSK3A 反向筛选与 ADMET 评估，",
         "再做任何合成决策，并请药物化学家复核所提议的路线——模板上成立不等于一份合成方案。"),
    ]:
        E += [P(f"<font color='#2a78d6'>{lead}</font> {rest}")]

    # 7 参考文献
    E += [P("7. 参考文献", "h1")]
    refs = [
        "Jensen JH. A graph-based genetic algorithm and generative model/Monte "
        "Carlo tree search for the exploration of chemical space. Chem Sci "
        "2019;10(12):3567-3572. doi:10.1039/c8sc05372c",
        "Brown N, Fiscato M, Segler MHS, Vaucher AC. GuacaMol: Benchmarking "
        "Models for de Novo Molecular Design. J Chem Inf Model "
        "2019;59(3):1096-1108. doi:10.1021/acs.jcim.8b00839",
        "Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. "
        "Quantifying the chemical beauty of drugs. Nat Chem 2012;4(2):90-98. "
        "doi:10.1038/nchem.1243",
        "Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of "
        "drug-like molecules based on molecular complexity and fragment "
        "contributions. J Cheminform 2009;1(1):8. doi:10.1186/1758-2946-1-8",
        "Baell JB, Holloway GA. New substructure filters for removal of pan "
        "assay interference compounds (PAINS) from screening libraries and for "
        "their exclusion in bioassays. J Med Chem 2010;53(7):2719-2740. "
        "doi:10.1021/jm901137j",
        "Lauretti E, Dincer O, Praticò D. Glycogen synthase kinase-3 signaling "
        "in Alzheimer's disease. Biochim Biophys Acta Mol Cell Res "
        "2020;1867(5):118664. doi:10.1016/j.bbamcr.2020.118664",
        "Genheden S, Thakkar A, Chadimová V, Reymond JL, Engkvist O, Bjerrum E. "
        "AiZynthFinder: a fast, robust and flexible open-source software for "
        "retrosynthetic planning. J Cheminform 2020;12(1):70. "
        "doi:10.1186/s13321-020-00472-1",
        "Huang K, Fu T, Gao W, et al. Therapeutics Data Commons: machine "
        f"learning datasets and tasks for drug discovery and development. "
        f"tdcommons.ai（软件与预训练 GSK3B 代理模型，访问日期 {meta['retrieved']}）。",
        f"ChEMBL 数据库，版本 {meta['chembl_release']}，EMBL-EBI。"
        f"靶点 CHEMBL262 的活性数据于 {meta['retrieved']} 经 ChEMBL REST API 获取。"
        f"ebi.ac.uk/chembl",
        "RDKit：开源化学信息学工具包。rdkit.org",
    ]
    for i, r in enumerate(refs, 1):
        E += [Paragraph(f"[{i}]&nbsp;&nbsp;{r}", s["ref"])]

    E += [Spacer(1, 10)]
    E += [Paragraph(
        "文献 [1]–[7] 经 PubMed 检索核对，生物活性数据来自 ChEMBL。"
        "全部十个设计分子均为计算提议，未经合成，也未做任何实验测定。", s["small"])]

    doc.build(E)
    return out_pdf
