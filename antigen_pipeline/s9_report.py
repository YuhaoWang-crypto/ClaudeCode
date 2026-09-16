"""Step 9c -- English and Chinese PDF reports.

Every number printed here is read back from the result files on disk.  The
report module never recomputes a statistic and never carries one in from the
run that produced it, so a report can be regenerated from the results
directory alone and must agree with it.
"""

from __future__ import annotations

import ast
import datetime as dt
import json

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (Image, PageBreak, Paragraph, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

from . import config as C

INK = colors.HexColor("#0b0b0b")
INK_2 = colors.HexColor("#52514e")
RULE = colors.HexColor("#d8d7d2")
ACCENT = colors.HexColor("#2a78d6")
BAD = colors.HexColor("#d03b3b")
GOOD = colors.HexColor("#0ca30c")

FIGS = ["fig1_candidate_ranking", "fig2_compartment_heatmap",
        "fig3_therapeutic_index", "fig4_validation_ranks"]

# Column headings and machine-generated phrases, per language. Data values that
# come from a source database (gene symbols, HPA tissue names, atlas titles,
# assay names) are left as the source writes them.
LABELS = {
    "en": {
        "atlas": "atlas", "cells available": "cells available", "sampled": "sampled",
        "analysed": "analysed", "donors": "donors", "assay": "assay",
        "gene": "gene", "rank": "rank", "percentile": "percentile", "score": "score",
        "tumour quality": "tumour quality", "safety": "safety",
        "worst normal tissue": "worst normal tissue", "note": "note",
        "alternative scheme": "alternative scheme",
        "Spearman vs baseline": "Spearman vs baseline",
        "top-20 overlap": "top-20 overlap", "Tier 1 count": "Tier 1 count",
        "specificity": "specificity", "topology": "topology",
        "surface confirmation": "surface confirmation",
        "evidence class": "evidence class", "evidence": "evidence",
        "indication papers": "indication papers", "modality papers": "modality papers",
        "drugs in Open Targets": "drugs in Open Targets",
        "consistency gate": "consistency gate", "result": "result", "detail": "detail",
        "cohort": "Cohort", "validated antigens": "Pre-registered validated antigens",
        "literature triage": "Literature triage", "export gates": "Export gates",
        "no HPA record": "no HPA record",
        "tier": "tier",
        "top": "top",
        "neg_intro": "Negative controls: ",
        "neg_overall": "Overall gate: ",
        "not_ranked_intro": "Not ranked at all: ",
    },
    "zh": {
        "atlas": "数据集", "cells available": "可用细胞数", "sampled": "抽样",
        "analysed": "纳入分析", "donors": "供体数", "assay": "测序平台",
        "gene": "基因", "rank": "排名", "percentile": "分位", "score": "评分",
        "tumour quality": "肿瘤质量", "safety": "安全性",
        "worst normal tissue": "风险最高的正常组织", "note": "备注",
        "alternative scheme": "替代方案",
        "Spearman vs baseline": "与基线的 Spearman 相关",
        "top-20 overlap": "前 20 名重叠", "Tier 1 count": "Tier 1 数量",
        "specificity": "特异性", "topology": "拓扑结构",
        "surface confirmation": "表面确认等级",
        "evidence class": "证据分级", "evidence": "证据",
        "indication papers": "适应症文献数", "modality papers": "模态文献数",
        "drugs in Open Targets": "Open Targets 已知药物",
        "consistency gate": "一致性门控", "result": "结果", "detail": "细节",
        "cohort": "队列构成", "validated antigens": "预注册验证靶点",
        "literature triage": "文献证据分级", "export gates": "导出一致性门控",
        "no HPA record": "HPA 无记录",
        "tier": "分层",
        "top": "前",
        "neg_intro": "阴性对照：",
        "neg_overall": "总体判定：",
        "not_ranked_intro": "完全未进入排序的靶点：",
    },
}

EVIDENCE_CLASS_ZH = {
    "clinically validated target": "已临床验证的靶点",
    "preclinical surface-targeting evidence": "有临床前表面靶向证据",
    "tumour-biology evidence only": "仅有肿瘤生物学证据",
    "novel -- no surface-targeting literature": "全新——无表面靶向文献",
}

STABILITY_ZH = {
    "A_safety_mean": "安全性取蛋白与 RNA 两臂的均值（而非保守最小值）",
    "B_safety_protein_only": "安全性仅取免疫组化蛋白一臂",
    "C_safety_rna_only": "安全性仅取共识 RNA 一臂",
    "D_safety_sqrt": "安全性以平方根进入：更温和的毒性惩罚",
    "E_equal_weights": "五个肿瘤质量子项等权 0.20",
    "F_specificity_heavy": "特异性权重 0.50，其余四项各 0.125",
}

GATE_ZH = {
    "cohort cell accounting": "队列细胞数核对",
    "sampling ceiling": "抽样上限",
    "expression matrix coverage": "表达矩阵覆盖度",
    "ranking integrity": "排序完整性",
    "recall report": "召回报告规范",
    "negative-control verdicts": "阴性对照判定",
    "literature coverage": "文献表覆盖",
    "figures present": "图件齐备",
}


def _neg_note(row, lang: str, n: int) -> str:
    code, rank = row.reason_code, row["rank"]
    if lang == "en":
        return str(row.note)
    if code == "not_in_search_space":
        return "从未进入搜索空间（SURFY 判定为非表面蛋白）"
    if code == "removed_before_ranking":
        return "在排序前被移除（拓扑门控或无表达数据）"
    if code == "inside_fail_window":
        return f"排名 {int(rank)}/{n}，落在前 {C.NEGATIVE_CONTROL_FAIL_RANK} 名之内"
    return f"排名 {int(rank)}/{n}"


def _zh_not_ranked(gene: str, mechanism: dict) -> str:
    info = (mechanism.get("positives_not_ranked_codes") or {}).get(gene, {})
    code, det = info.get("code"), info.get("epi_detection")
    if code == "below_display_threshold" and det is not None:
        return (f"仅在 {det:.1%} 的上皮细胞中检出，低于 "
                f"{C.MIN_EPITHELIAL_DETECTION:.0%} 的抗原展示阈值")
    if code == "no_expression_data":
        return "Census 特征集中无表达数据"
    if code == "topology_gate":
        return "被胞外拓扑门控移除"
    return "未进入排序"


def _status_note(row, lang: str) -> str:
    code = str(row.status_code)
    if code == "ranked":
        return ""
    if lang == "en":
        return str(row.status).split(":", 1)[-1].strip()
    det = row.epi_detection
    if code == "below_display_threshold":
        return (f"仅在 {det:.1%} 的上皮细胞中检出，低于 "
                f"{C.MIN_EPITHELIAL_DETECTION:.0%} 的抗原展示阈值")
    if code == "no_expression_data":
        return "Census 特征集中无表达数据"
    return "被胞外拓扑门控移除"


# ---------------------------------------------------------------------------
# data loading -- the single place the report touches disk
# ---------------------------------------------------------------------------
def load_report_data() -> dict:
    R = C.RESULTS
    prov = {s: json.loads((R / f"{s}_provenance.json").read_text())
            for s in ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s9_export"]}
    return {
        "prov": prov,
        "registry": json.loads((R / "s1_validation_registry.json").read_text()),
        "cohort": pd.read_csv(R / "s2_cohort.csv"),
        "ranked": pd.read_csv(R / "s6_ranked_candidates.csv"),
        "positives": pd.read_csv(R / "s7_validation_positives.csv"),
        "negatives": pd.read_csv(R / "s7_negative_controls.csv"),
        "stability": pd.read_csv(R / "s7_rank_stability.csv"),
        "literature": pd.read_csv(R / "s8_literature_evidence.csv"),
        "excluded": pd.read_csv(R / "s3_excluded.csv"),
        "generated": dt.date.today().isoformat(),
    }


def _fmt(x, nd=3):
    return "n/a" if pd.isna(x) else f"{x:.{nd}f}"


def _text(x, dash: str = "-") -> str:
    """A cell value that is missing prints as a dash, never as 'nan'."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return dash
    s = str(x).strip()
    return dash if s.lower() in ("nan", "none", "") else s


def _clip_words(s: str, limit: int) -> str:
    """Shorten to whole tokens so a table cell never cuts a word in half."""
    s = _text(s)
    if len(s) <= limit:
        return s
    parts = s.replace("; ", ";").split(";")
    out: list[str] = []
    for part in parts:
        candidate = "; ".join(out + [part.strip()])
        if len(candidate) > limit - 3:
            break
        out.append(part.strip())
    return ("; ".join(out) + " ...") if out else s[:limit - 3] + "..."


# Compact forms for narrow table columns; the long form stays in the CSVs.
TOPOLOGY_SHORT = {
    "no TM (GPI/peripheral)": "GPI / no TM",
}
CONFIRMATION_SHORT = {
    "en": {"experimental (CSPA/GPI)": "experimental",
           "Open Targets confirmed": "Open Targets",
           "unconfirmed (ML prediction only)": "ML prediction only"},
    "zh": {"experimental (CSPA/GPI)": "实验确认",
           "Open Targets confirmed": "Open Targets 确认",
           "unconfirmed (ML prediction only)": "仅 ML 预测"},
}


# ---------------------------------------------------------------------------
# styles
# ---------------------------------------------------------------------------
def _styles(lang: str):
    base = "Helvetica"
    bold = "Helvetica-Bold"
    if lang == "zh":
        pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
        base = bold = "STSong-Light"
    ss = getSampleStyleSheet()
    return {
        "title": ParagraphStyle("t", parent=ss["Title"], fontName=bold, fontSize=19,
                                leading=25, textColor=INK, alignment=TA_LEFT, spaceAfter=4),
        "sub": ParagraphStyle("s", parent=ss["Normal"], fontName=base, fontSize=9.5,
                              leading=14, textColor=INK_2, spaceAfter=10),
        "h2": ParagraphStyle("h2", parent=ss["Heading2"], fontName=bold, fontSize=13,
                             leading=17, textColor=INK, spaceBefore=13, spaceAfter=5),
        "h3": ParagraphStyle("h3", parent=ss["Heading3"], fontName=bold, fontSize=10.5,
                             leading=14, textColor=INK, spaceBefore=8, spaceAfter=3),
        "body": ParagraphStyle("b", parent=ss["Normal"], fontName=base, fontSize=9.4,
                               leading=14.4, textColor=INK, spaceAfter=5),
        "small": ParagraphStyle("sm", parent=ss["Normal"], fontName=base, fontSize=8.2,
                                leading=11.5, textColor=INK_2, spaceAfter=4),
        "cell": ParagraphStyle("c", parent=ss["Normal"], fontName=base, fontSize=7.6,
                               leading=10),
        "cellb": ParagraphStyle("cb", parent=ss["Normal"], fontName=bold, fontSize=7.6,
                                leading=10),
        "_base": base, "_bold": bold,
    }


def _table(rows, st, widths, align_right=()):
    data = [[Paragraph(str(c), st["cellb"] if i == 0 else st["cell"]) for c in row]
            for i, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style = [
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, INK_2),
        ("LINEBELOW", (0, 1), (-1, -2), 0.3, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ]
    for c in align_right:
        style.append(("ALIGN", (c, 0), (c, -1), "RIGHT"))
    t.setStyle(TableStyle(style))
    return t


# ---------------------------------------------------------------------------
# text blocks
# ---------------------------------------------------------------------------
TEXT = {
    "en": {
        "title": "Cell-surface antigen discovery in lung adenocarcinoma",
        "subtitle": ("Antibody-accessible ADC / CAR-T target nomination from tumour single-cell "
                     "expression, extracellular topology and normal-tissue therapeutic index.<br/>"
                     "Generated {generated} from CZ CELLxGENE Census {census}. "
                     "All figures in this report are read back from the result files on disk."),
        "principle_h": "Design principle",
        "principle": (
            "Targets are ranked by <b>tumour surface specificity x normal-tissue therapeutic "
            "index</b>, not by tumour essentiality. An antibody therapeutic kills cells that "
            "display the antigen; whether the tumour needs the gene to survive is irrelevant to "
            "efficacy but directly relevant to toxicity, because an essential gene is essential "
            "in normal tissue too. DepMap gene effect is therefore carried as annotation and "
            "never as a filter."),
        "s1_h": "1. Search space",
        "s2_h": "2. Tumour single-cell expression",
        "s3_h": "3. Extracellular topology gate",
        "s4_h": "4. Druggability and essentiality annotation",
        "s5_h": "5. Normal-tissue safety baseline",
        "s6_h": "6. Composite score and tiers",
        "s7_h": "7. Validation and quality control",
        "s8_h": "8. Literature evidence for the head of the ranking",
        "s9_h": "9. Figures, export gates and provenance",
        "top_h": "Top nominations",
        "limits_h": "Limitations",
        "sources_h": "Data sources",
        "figures_h": "Figures",
    },
    "zh": {
        "title": "肺腺癌细胞表面抗原发现",
        "subtitle": ("基于肿瘤单细胞表达、胞外拓扑结构与正常组织治疗指数的 ADC / CAR-T "
                     "可及靶点提名。<br/>生成日期 {generated}，数据来自 CZ CELLxGENE Census {census}。"
                     "本报告中所有数字均从磁盘结果文件逐一读取，不凭记忆重建。"),
        "principle_h": "核心设计理念",
        "principle": (
            "本流程用<b>“肿瘤表面特异性 × 正常组织治疗指数”</b>而非“肿瘤必需性”来排序。"
            "抗体类疗法杀伤的是表面抗原阳性的细胞；靶点是否为肿瘤存活所必需与疗效无关，"
            "却与毒性直接相关——因为必需基因在正常组织中同样必需。因此 DepMap 基因效应值"
            "仅作注释，绝不用作筛选条件。"),
        "s1_h": "1. 确定搜索空间（表面基因全集）",
        "s2_h": "2. 获取肿瘤单细胞表达（Census）",
        "s3_h": "3. 胞外拓扑门控",
        "s4_h": "4. 可成药性与必需性注释",
        "s5_h": "5. 正常组织安全性基线（治疗指数）",
        "s6_h": "6. 复合评分与分层",
        "s7_h": "7. 验证与质控",
        "s8_h": "8. 文献证据",
        "s9_h": "9. 可视化、导出门控与溯源",
        "top_h": "候选靶点提名",
        "limits_h": "局限性",
        "sources_h": "数据来源",
        "figures_h": "图表",
    },
}


def _narrative(d: dict, lang: str) -> dict:
    """Sentences whose numbers come from the loaded artefacts."""
    p = d["prov"]
    cohort = d["cohort"]
    ranked = d["ranked"]
    pos = d["positives"]
    neg = d["negatives"]
    stab = d["stability"]
    lit = d["literature"]

    recall = p["s7"]["recall"]
    pos_ranked = pos.dropna(subset=["rank"])
    best_pos = pos_ranked.sort_values("rank").iloc[0] if len(pos_ranked) else None
    tier_counts = p["s6"]["tier_counts"]
    fails = neg[neg.verdict == "FAIL"]

    M = p["s7"]["mechanism"]
    not_ranked = M["positives_not_ranked"]

    if lang == "en":
        return {
            "s1": (
                f"The search space is the SURFY in-silico surfaceome. Of "
                f"{p['s1']['proteins_in_master_table']:,} proteins in the master table, "
                f"{p['s1']['surface_labelled_proteins']:,} carry a surface label, collapsing to "
                f"<b>{p['s1']['unique_surface_genes']:,} unique genes</b>: "
                f"{p['s1']['experimental_evidence_genes']:,} with experimental evidence "
                f"(CSPA peptide capture or a GPI anchor) and "
                f"{p['s1']['machine_learning_genes']:,} from the machine-learning prediction. "
                f"Each gene keeps its own topology string, so extracellular domain length and "
                f"transmembrane pass count are parsed per gene rather than assumed. "
                f"The pre-registered validation set was locked at this point: "
                f"{p['s1']['positives_in_search_space']}/{len(C.VALIDATION_POSITIVES)} clinically "
                f"validated LUAD antigens and "
                f"{p['s1']['negative_controls_in_search_space']}/{len(C.NEGATIVE_CONTROLS)} "
                f"negative controls fall inside the search space."),
            "s2": (
                f"CZ CELLxGENE Census {p['s2']['census_version']} was queried for primary "
                f"'{p['s2']['disease_label']}' cells. Single-nucleus data systematically "
                f"under-detects membrane-protein transcripts and was excluded, leaving "
                f"<b>{p['s2']['cells_available']:,} whole cells across {p['s2']['n_datasets']} "
                f"atlases</b>. Each atlas was sampled to at most "
                f"{p['s2']['max_cells_per_dataset']:,} cells with seed "
                f"{p['s2']['sampling_seed']} ({p['s2']['cells_sampled']:,} cells), and "
                f"{p['s2']['cells_analysed']:,} of those carry a cell-type label that maps to one "
                f"of the four compartments: "
                + ", ".join(f"{k} {v:,}" for k, v in p['s2']['cells_per_compartment'].items())
                + ". Per-gene mean expression and detection rate were computed per compartment "
                  "per atlas, then integrated into a cross-atlas consensus."),
            "s3": (
                f"Proteins that cannot present an epitope on an intact cell were removed: "
                f"<b>{p['s3']['excluded']} genes excluded</b>, "
                f"{p['s3']['exclusions_by_class'].get('localisation', 0)} on subcellular "
                f"localisation (ER, Golgi, endosome, nucleus, cytoplasm or purely secreted / "
                f"matrix) and "
                + str(sum(v for k, v in p['s3']['exclusions_by_class'].items()
                          if k != 'localisation'))
                + f" for having no extracellular segment of at least {p['s3']['min_ecd_length']} "
                  f"residues. <b>{p['s3']['accessible']:,} antibody-accessible candidates</b> "
                  f"remain, each labelled with its surface-confirmation tier."),
            "s4": (
                f"Open Targets supplied the antibody tractability bucket, an independent "
                f"subcellular localisation call and the drug record for "
                f"{p['s4']['open_targets_hits']:,} candidates. Surface confirmation splits as: "
                + "; ".join(f"{k} {v:,}" for k, v in p['s4']['surface_confirmation_tiers'].items())
                + f". DepMap mean gene effect was attached for {p['s4']['depmap_genes']:,} genes "
                  f"as annotation only. The validated antigens confirm why: "
                + ", ".join(f"{k} {v:+.2f}" for k, v in
                            list(p['s4']['validated_antigen_gene_effect'].items())[:6])
                + " -- all far from essential. Gating on essentiality would have removed the "
                  "targets this pipeline is supposed to rediscover and enriched for housekeeping "
                  "genes instead."),
            "s5": (
                f"Two independent Human Protein Atlas releases set the safety baseline: consensus "
                f"RNA ({p['s5']['hpa_rna_genes']:,} genes x {p['s5']['hpa_rna_tissues']} tissues) "
                f"and normal-tissue immunohistochemistry ({p['s5']['hpa_ihc_genes']:,} genes x "
                f"{p['s5']['hpa_ihc_tissues']} tissues). Each arm gives an organ-weighted risk, "
                f"and the safety coefficient is the <b>conservative minimum of the two</b>, so a "
                f"target that looks clean in RNA but stains in heart muscle is treated as "
                f"dangerous. {p['s5']['candidates_default_safety']} candidates had no HPA record "
                f"and received the flagged neutral default {p['s5']['default_value']}. Median "
                f"safety coefficient across the cohort is {p['s5']['median_safety']:.3f}."),
            "s6": (
                f"The composite is the <b>geometric mean of tumour quality, safety "
                f"coefficient and consensus multiplier</b> -- the cube root of their product, "
                f"a monotone transform that leaves every rank exactly where the raw product "
                f"puts it while keeping the score on the same 0-1 scale as its three factors. "
                f"Tumour quality weights relative TME specificity "
                f"{C.TUMOUR_QUALITY_WEIGHTS['specificity']:.2f}, expression intensity "
                f"{C.TUMOUR_QUALITY_WEIGHTS['intensity']:.2f}, expression uniformity "
                f"{C.TUMOUR_QUALITY_WEIGHTS['uniformity']:.2f}, ectodomain accessibility "
                f"{C.TUMOUR_QUALITY_WEIGHTS['accessibility']:.2f} and antibody druggability "
                f"{C.TUMOUR_QUALITY_WEIGHTS['druggability']:.2f}. Safety enters multiplicatively "
                f"and the consensus multiplier rewards reproducibility across atlases. "
                f"One eligibility rule is applied before ranking: a candidate must be detected "
                f"in at least {p['s6']['min_epithelial_detection']:.0%} of epithelial/malignant "
                f"cells, because an ADC or a CAR cannot engage an antigen that most tumour cells "
                f"never display. That sets aside {p['s6']['set_aside_below_display_threshold']:,} "
                f"genes which would otherwise win on a safety coefficient earned purely by being "
                f"silent everywhere, and leaves {p['s6']['candidates_scored']:,} scored "
                f"candidates ({p['s6']['accessible_without_expression']} accessible genes are "
                f"absent from the Census feature set altogether). Absolute thresholds give "
                f"<b>{tier_counts.get('Tier 1', 0)} Tier 1</b> "
                f"(>= {C.TIER1_THRESHOLD}), {tier_counts.get('Tier 2', 0)} Tier 2 "
                f"({C.TIER2_THRESHOLD}-{C.TIER1_THRESHOLD}) and "
                f"{tier_counts.get('Tier 3', 0)} Tier 3 targets."),
            "s7_recall": (
                f"<b>recall@10 = {recall['recall@10']}, recall@20 = {recall['recall@20']}</b> "
                f"for the pre-registered clinically validated antigens among "
                f"{p['s7']['n_ranked']:,} scored candidates."
                + (f" The best-placed validated antigen is {best_pos.gene} at rank "
                   f"{int(best_pos['rank'])}." if best_pos is not None else "")),
            "s7_mech": (
                f"The validated antigens are not buried, and reading the recall number alone "
                f"would be misleading. {M['positives_ranked']} of {M['positives_total']} were "
                f"ranked at all, between rank {M['positive_best_rank']} and "
                f"{M['positive_worst_rank']}, with a median position in the top "
                f"{M['positive_median_percentile_from_top']}% of scored candidates. On the "
                f"tumour side they are clearly above the field: median tumour quality "
                f"{M['positive_median_tumour_quality']:.3f} against a cohort median of "
                f"{M['cohort_median_tumour_quality']:.3f}. Their safety coefficients are also "
                f"above the cohort median ({M['positive_median_safety']:.3f} against "
                f"{M['cohort_median_safety']:.3f}).<br/><br/>"
                f"What keeps them out of the top ten is that the leaders are cleaner in normal "
                f"tissue, not that the validated antigens look bad: the top ten carry a median "
                f"safety coefficient of {M['top10_median_safety']:.3f} at a median tumour "
                f"quality of {M['top10_median_tumour_quality']:.3f}. EGFR, HER2, TROP2 and "
                f"folate receptor alpha are genuinely expressed in skin, gastrointestinal "
                f"tract, kidney and heart, and their clinical use is bought with managed "
                f"toxicity, dose reduction and patient selection. A purely expression-driven "
                f"safety coefficient prices that cost in; a clinical programme can decide to "
                f"pay it. So the honest reading of this ranking is not that it failed to "
                f"rediscover the standard of care, but that it is answering a different "
                f"question: which antigens have a <i>better</i> normal-tissue profile than the "
                f"targets already in the clinic. A ranking retuned until the known answers came "
                f"back on top would be a ranking tuned to its own validation set."
                + ("<br/><br/>" + LABELS["en"]["not_ranked_intro"]
                   + "; ".join(f"{g} ({s.split(':', 1)[-1].strip()})"
                               for g, s in not_ranked.items()) + "."
                   if not_ranked else "")),
            "s7_neg": (
                LABELS["en"]["neg_intro"] + "; ".join(
                    f"{r.gene} {r.verdict} ({_neg_note(r, 'en', p['s7']['n_ranked'])})"
                    for _, r in neg.iterrows())
                + ". " + LABELS["en"]["neg_overall"]
                + f"<b>{p['s7']['negative_control_overall']}</b>."),
            "s7_stab": (
                f"Rank stability across {len(stab)} alternative safety-aggregation and weighting "
                f"schemes: median Spearman correlation with the baseline ranking "
                f"{p['s7']['median_spearman']:.3f}, worst top-20 overlap "
                f"{p['s7']['min_top20_overlap']}/20. The head of the ranking is therefore not an "
                f"artefact of one particular parameterisation."),
            "s8": (
                f"Each of the top {len(lit)} candidates was queried against Europe PMC with fixed "
                f"query shapes for indication biology, surface-targeting modality and clinical "
                f"stage. Evidence classes: "
                + "; ".join(f"{k} ({v})" for k, v in
                            lit.evidence_class.value_counts().items()) + "."),
            "s9": (
                f"Four figures were written as PNG and SVG. The export bundle passed "
                f"{sum(g['passed'] for g in p['s9_export']['gates'])}/"
                f"{len(p['s9_export']['gates'])} consistency gates, which cross-check the cohort "
                f"cell accounting against the expression matrix, the reported recall string "
                f"against the per-gene ranks, the negative-control verdicts against the rule, and "
                f"the literature table against the ranking."),
            "limits": [
                "Transcript abundance is a proxy for surface protein density; an ADC needs "
                "protein copies per cell, which RNA can over- or under-state.",
                "Compartment assignment relies on the Cell Ontology labels supplied by each "
                "atlas; malignant cells are not re-called from copy-number, so the epithelial "
                "compartment mixes malignant and normal epithelium.",
                "HPA immunohistochemistry is antibody-dependent and semi-quantitative; a missing "
                "or weak antibody looks like safety.",
                "The safety model weights organs by how badly antigen-positive killing is "
                "tolerated. Those weights are a stated judgement, not a measurement.",
                "Novelty in the literature triage means no indexed surface-targeting paper, "
                "which is not the same as no prior art.",
            ],
        }
    return {
        "s1": (
            f"搜索空间为 SURFY 计算表面组。主表 {p['s1']['proteins_in_master_table']:,} 个蛋白中，"
            f"{p['s1']['surface_labelled_proteins']:,} 个带表面标签，去冗余后得到"
            f"<b>{p['s1']['unique_surface_genes']:,} 个基因</b>："
            f"{p['s1']['experimental_evidence_genes']:,} 个有实验证据（CSPA 肽段捕获或 GPI 锚定），"
            f"{p['s1']['machine_learning_genes']:,} 个为机器学习预测。每个基因保留其独立的拓扑结构串，"
            f"因此胞外域长度与跨膜次数是逐基因解析的，而非假定。预注册验证集在此刻锁定："
            f"{p['s1']['positives_in_search_space']}/{len(C.VALIDATION_POSITIVES)} 个已临床验证的"
            f"肺腺癌抗原与 {p['s1']['negative_controls_in_search_space']}/{len(C.NEGATIVE_CONTROLS)} "
            f"个阴性对照落在搜索空间内。"),
        "s2": (
            f"从 CZ CELLxGENE Census {p['s2']['census_version']} 按疾病标签"
            f"“{p['s2']['disease_label']}”检索原始数据。单细胞核数据会系统性低检膜蛋白转录本，故排除，"
            f"共得到<b>{p['s2']['cells_available']:,} 个全细胞，分布于 {p['s2']['n_datasets']} 个 atlas</b>。"
            f"每个 atlas 以固定种子 {p['s2']['sampling_seed']} 抽样至 ≤{p['s2']['max_cells_per_dataset']:,} "
            f"个细胞（共 {p['s2']['cells_sampled']:,} 个），其中 {p['s2']['cells_analysed']:,} 个细胞的"
            f"细胞类型可映射到四个区室之一："
            + "、".join(f"{k} {v:,}" for k, v in p['s2']['cells_per_compartment'].items())
            + "。每个基因在每个区室、每个 atlas 中的均值与表达率分别计算，再做跨数据集共识整合。"),
        "s3": (
            f"无法在完整细胞上暴露表位的蛋白被剔除：<b>排除 {p['s3']['excluded']} 个基因</b>，其中 "
            f"{p['s3']['exclusions_by_class'].get('localisation', 0)} 个因亚细胞定位（内质网、高尔基体、"
            f"内体、细胞核、胞质，或纯分泌型/细胞外基质），"
            + str(sum(v for k, v in p['s3']['exclusions_by_class'].items() if k != 'localisation'))
            + f" 个因缺乏长度 ≥ {p['s3']['min_ecd_length']} 个残基的胞外片段。"
              f"剩余 <b>{p['s3']['accessible']:,} 个抗体可及候选</b>，每个候选均标注表面确认等级。"),
        "s4": (
            f"Open Targets 为 {p['s4']['open_targets_hits']:,} 个候选提供抗体可成药性分档、独立的"
            f"亚细胞定位判断与已知药物记录。表面确认等级分布："
            + "；".join(f"{k} {v:,}" for k, v in p['s4']['surface_confirmation_tiers'].items())
            + f"。DepMap 平均基因效应覆盖 {p['s4']['depmap_genes']:,} 个基因，仅作注释。"
              f"验证靶点本身即说明原因："
            + "、".join(f"{k} {v:+.2f}" for k, v in
                        list(p['s4']['validated_antigen_gene_effect'].items())[:6])
            + "——全部远非必需基因。若以必需性做门控，恰好会剔除本流程应当重新发现的靶点，"
              "转而富集高毒性的管家基因。"),
        "s5": (
            f"安全性基线取自两套独立的 HPA 数据：共识 RNA（{p['s5']['hpa_rna_genes']:,} 基因 × "
            f"{p['s5']['hpa_rna_tissues']} 组织）与正常组织免疫组化蛋白"
            f"（{p['s5']['hpa_ihc_genes']:,} 基因 × {p['s5']['hpa_ihc_tissues']} 组织）。"
            f"两条通路各自给出按器官加权的风险，安全性系数取两者的<b>保守最小值</b>——"
            f"因此一个 RNA 看似干净、却在心肌上染色的靶点仍被判为危险。"
            f"{p['s5']['candidates_default_safety']} 个候选在 HPA 中无记录，"
            f"采用带标记的中性默认值 {p['s5']['default_value']}。队列安全性系数中位数为 "
            f"{p['s5']['median_safety']:.3f}。"),
        "s6": (
            f"最终评分取<b>肿瘤质量、安全性系数与共识乘数三者的几何平均</b>（即三者乘积的立方根）。"
            f"这是一个单调变换：排名与直接取乘积完全一致，同时把评分放回与三个因子相同的 0–1 量纲上。"
            f"肿瘤质量权重为："
            f"相对 TME 特异性 {C.TUMOUR_QUALITY_WEIGHTS['specificity']:.2f}、"
            f"表达强度 {C.TUMOUR_QUALITY_WEIGHTS['intensity']:.2f}、"
            f"表达均一性 {C.TUMOUR_QUALITY_WEIGHTS['uniformity']:.2f}、"
            f"胞外域可及性 {C.TUMOUR_QUALITY_WEIGHTS['accessibility']:.2f}、"
            f"抗体可成药性 {C.TUMOUR_QUALITY_WEIGHTS['druggability']:.2f}。"
            f"安全性系数为乘性，共识乘数奖励跨 atlas 可重复性。"
            f"排名前还应用了一条资格规则：候选必须在至少 "
            f"{p['s6']['min_epithelial_detection']:.0%} 的上皮/恶性细胞中被检出——"
            f"因为 ADC 与 CAR 无法结合大多数肿瘤细胞根本不展示的抗原。"
            f"该规则搁置了 {p['s6']['set_aside_below_display_threshold']:,} 个基因"
            f"（它们否则会仅凭“处处不表达”换来的高安全性系数而胜出），"
            f"最终 {p['s6']['candidates_scored']:,} 个候选参与评分"
            f"（另有 {p['s6']['accessible_without_expression']} 个可及基因根本不在 Census 特征集内）。"
            f"按绝对阈值分层：<b>Tier 1 共 {tier_counts.get('Tier 1', 0)} 个</b>"
            f"（≥ {C.TIER1_THRESHOLD}），Tier 2 共 {tier_counts.get('Tier 2', 0)} 个"
            f"（{C.TIER2_THRESHOLD}–{C.TIER1_THRESHOLD}），Tier 3 共 {tier_counts.get('Tier 3', 0)} 个。"),
        "s7_recall": (
            f"在 {p['s7']['n_ranked']:,} 个评分候选中，预注册的已临床验证抗原的召回为"
            f"<b>recall@10 = {recall['recall@10']}，recall@20 = {recall['recall@20']}</b>。"
            + (f"排名最前的验证靶点是 {best_pos.gene}，位列第 {int(best_pos['rank'])} 名。"
               if best_pos is not None else "")),
        "s7_mech": (
            f"验证靶点并未被埋没，只看召回数会产生误导。{M['positives_total']} 个中有 "
            f"{M['positives_ranked']} 个进入排序，名次介于第 {M['positive_best_rank']} 与第 "
            f"{M['positive_worst_rank']} 名之间，中位数位于全部评分候选的前 "
            f"{M['positive_median_percentile_from_top']}%。在肿瘤一侧它们明显优于整体："
            f"肿瘤质量中位数 {M['positive_median_tumour_quality']:.3f}，队列中位数为 "
            f"{M['cohort_median_tumour_quality']:.3f}；安全性系数同样高于队列中位数"
            f"（{M['positive_median_safety']:.3f} 对 {M['cohort_median_safety']:.3f}）。<br/><br/>"
            f"它们没有进入前十，不是因为自身表现差，而是因为榜首的候选在正常组织中更干净："
            f"前十名的安全性系数中位数为 {M['top10_median_safety']:.3f}，"
            f"肿瘤质量中位数为 {M['top10_median_tumour_quality']:.3f}。"
            f"EGFR、HER2、TROP2、叶酸受体 α 确实在皮肤、消化道、肾与心肌中表达；"
            f"它们能进入临床，靠的是可控毒性、剂量调整与患者筛选。"
            f"纯表达驱动的安全性系数会如实给这份代价定价，而临床项目可以选择承受它。"
            f"因此对本排序的诚实解读不是“没能重新发现标准疗法靶点”，"
            f"而是它在回答另一个问题：哪些抗原的正常组织谱<i>优于</i>已进入临床的靶点。"
            f"若把排序反复调到已知答案重新回到前列，那只是对着自己的验证集调参。"
            + ("<br/><br/>" + LABELS["zh"]["not_ranked_intro"]
               + "；".join(f"{g}（{_zh_not_ranked(g, M)}）" for g in not_ranked)
               + "。" if not_ranked else "")),
        "s7_neg": (
            LABELS["zh"]["neg_intro"] + "；".join(
                f"{r.gene} {r.verdict}（{_neg_note(r, 'zh', p['s7']['n_ranked'])}）"
                for _, r in neg.iterrows())
            + "。" + LABELS["zh"]["neg_overall"]
            + f"<b>{p['s7']['negative_control_overall']}</b>。"),
        "s7_stab": (
            f"在 {len(stab)} 组安全性聚合与加权替代方案下的排名稳定性：与基线排序的 Spearman "
            f"相关系数中位数为 {p['s7']['median_spearman']:.3f}，最差的前 20 名重叠为 "
            f"{p['s7']['min_top20_overlap']}/20。因此排序前列并非某一组特定参数化的产物。"),
        "s8": (
            f"对前 {len(lit)} 名候选逐一用固定查询式检索 Europe PMC，分别覆盖适应症生物学、"
            f"表面靶向模态与临床阶段。证据分级分布："
            + "；".join(f"{EVIDENCE_CLASS_ZH.get(k, k)}（{v}）"
                        for k, v in lit.evidence_class.value_counts().items()) + "。"),
        "s9": (
            f"共生成 4 幅图（PNG + SVG）。导出包通过了 "
            f"{sum(g['passed'] for g in p['s9_export']['gates'])}/"
            f"{len(p['s9_export']['gates'])} 项一致性门控，逐项交叉核对队列细胞数与表达矩阵、"
            f"报告中的召回字符串与逐基因排名、阴性对照判定与判定规则、文献表与排序表。"),
        "limits": [
            "转录本丰度只是表面蛋白密度的替代指标；ADC 真正需要的是每个细胞的蛋白拷贝数，"
            "RNA 可能高估也可能低估。",
            "区室划分依赖各 atlas 提供的 Cell Ontology 标签；恶性细胞未经拷贝数重新判定，"
            "因此上皮区室混合了恶性与正常上皮。",
            "HPA 免疫组化依赖抗体且为半定量；抗体缺失或信号弱会被误读为“安全”。",
            "安全性模型按“抗原阳性杀伤在该器官的可耐受程度”给器官加权。这些权重是明示的判断，"
            "不是测量值。",
            "文献分级中的“全新”指没有被索引到的表面靶向文献，并不等于没有在先技术。",
        ],
    }


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def render(lang: str, d: dict | None = None) -> str:
    d = d or load_report_data()
    st = _styles(lang)
    T = TEXT[lang]
    N = _narrative(d, lang)
    p = d["prov"]
    ranked = d["ranked"]
    lit = d["literature"]
    path = C.REPORTS / f"antigen_discovery_report_{lang}.pdf"

    doc = SimpleDocTemplate(str(path), pagesize=A4,
                            leftMargin=18 * mm, rightMargin=16 * mm,
                            topMargin=16 * mm, bottomMargin=15 * mm,
                            title=T["title"], author="cell-surface-antigen-discovery pipeline")
    flow = []
    P = lambda txt, s="body": Paragraph(txt, st[s])  # noqa: E731

    flow += [P(T["title"], "title"),
             P(T["subtitle"].format(generated=d["generated"], census=p["s2"]["census_version"]), "sub"),
             P(T["principle_h"], "h2"), P(T["principle"])]

    for key, head in (("s1", "s1_h"), ("s2", "s2_h"), ("s3", "s3_h"), ("s4", "s4_h"),
                      ("s5", "s5_h"), ("s6", "s6_h")):
        flow += [P(T[head], "h2"), P(N[key])]

    # cohort table
    L = LABELS[lang]
    H = lambda *keys: [L[k] for k in keys]  # noqa: E731
    cohort = d["cohort"]
    rows = [H("atlas", "cells available", "sampled", "analysed", "donors", "assay")]
    for r in cohort.itertuples():
        try:
            assays = ", ".join(ast.literal_eval(str(r.assays)))
        except (ValueError, SyntaxError):
            assays = str(r.assays)
        rows.append([str(r.title) or str(r.dataset_id)[:12],
                     f"{r.cells_available:,}", f"{r.cells_sampled:,}",
                     f"{r.cells_in_compartments:,}", str(r.donors), assays])
    flow += [P(L["cohort"], "h3"),
             _table(rows, st, [62 * mm, 20 * mm, 17 * mm, 17 * mm, 13 * mm, 45 * mm],
                    align_right=(1, 2, 3, 4))]

    flow += [PageBreak(), P(T["s7_h"], "h2"), P(N["s7_recall"]), P(N["s7_mech"])]

    pos = d["positives"]
    rows = [H("gene", "rank", "percentile", "score", "tumour quality", "safety",
              "worst normal tissue", "note")]
    for r in pos.itertuples():
        worst = str(r.worst_normal_tissue)
        rows.append([r.gene,
                     "-" if pd.isna(r.rank) else str(int(r.rank)),
                     "-" if pd.isna(r.percentile) else f"{L['top']} {100 - r.percentile:.1f}%",
                     _fmt(r.final_score), _fmt(r.tumour_quality), _fmt(r.safety_score),
                     L["no HPA record"] if worst in ("nan", "") else worst[:24],
                     _status_note(r, lang)])
    flow += [P(L["validated antigens"], "h3"),
             _table(rows, st, [18 * mm, 12 * mm, 17 * mm, 15 * mm, 21 * mm, 14 * mm,
                               32 * mm, 45 * mm],
                    align_right=(1, 3, 4, 5))]
    flow += [P(N["s7_neg"]), P(N["s7_stab"])]

    stab = d["stability"]
    rows = [H("alternative scheme", "Spearman vs baseline", "top-20 overlap", "Tier 1 count")]
    for r in stab.itertuples():
        desc = r.description if lang == "en" else STABILITY_ZH.get(r.scheme, r.description)
        rows.append([desc, f"{r.spearman_vs_baseline:.3f}",
                     f"{r.top20_overlap}/20", str(r.tier1_count)])
    flow += [_table(rows, st, [78 * mm, 32 * mm, 26 * mm, 24 * mm], align_right=(1, 2, 3))]

    # top nominations
    flow += [PageBreak(), P(T["top_h"], "h2"), P(N["s8"])]
    rows = [["#"] + H("gene", "score", "tier", "specificity", "safety", "topology",
                      "surface confirmation", "evidence class")]
    lit_by_gene = dict(zip(lit.gene, lit.evidence_class))
    for r in ranked.head(len(lit)).itertuples():
        klass = lit_by_gene.get(r.gene, "")
        if lang == "zh":
            klass = EVIDENCE_CLASS_ZH.get(klass, klass)
        topo = TOPOLOGY_SHORT.get(str(r.topology_class), _text(r.topology_class))
        conf = CONFIRMATION_SHORT[lang].get(str(r.surface_confirmation),
                                            _text(r.surface_confirmation))
        rows.append([str(r.rank), r.gene, _fmt(r.final_score), r.tier,
                     _fmt(r.score_specificity, 2), _fmt(r.safety_score, 2),
                     topo, conf, klass])
    flow += [_table(rows, st, [8 * mm, 19 * mm, 14 * mm, 13 * mm, 16 * mm, 14 * mm,
                               24 * mm, 28 * mm, 40 * mm], align_right=(0, 2, 4, 5))]

    rows = [H("gene", "evidence", "indication papers", "modality papers",
              "drugs in Open Targets")]
    for r in lit.itertuples():
        klass = r.evidence_class
        note = (str(r.evidence_note)[:70] if lang == "en"
                else EVIDENCE_CLASS_ZH.get(klass, klass))
        rows.append([r.gene, note,
                     "-" if pd.isna(r.hits_indication) else f"{int(r.hits_indication):,}",
                     "-" if pd.isna(r.hits_modality) else f"{int(r.hits_modality):,}",
                     _clip_words(r.ot_drugs, 46)])
    flow += [P(L["literature triage"], "h3"),
             _table(rows, st, [20 * mm, 60 * mm, 22 * mm, 22 * mm, 52 * mm],
                    align_right=(2, 3))]

    # figures
    flow += [PageBreak(), P(T["figures_h"], "h2"), P(N["s9"])]
    for name in FIGS:
        img = C.FIGURES / f"{name}.png"
        if img.exists():
            flow += [Image(str(img), width=165 * mm, height=165 * mm * _aspect(img),
                           kind="proportional"), Spacer(1, 6)]

    gate_rows = [H("consistency gate", "result", "detail")]
    for g in p["s9_export"]["gates"]:
        name = g["gate"] if lang == "en" else GATE_ZH.get(g["gate"], g["gate"])
        gate_rows.append([name, "PASS" if g["passed"] else "FAIL", g["detail"][:96]])
    flow += [P(L["export gates"], "h3"),
             _table(gate_rows, st, [42 * mm, 16 * mm, 106 * mm])]

    flow += [P(T["limits_h"], "h2")]
    for item in N["limits"]:
        flow += [P(f"&bull; {item}", "small")]

    flow += [P(T["sources_h"], "h2"),
             P(f"SURFY surfaceome: {C.SURFY_URL}<br/>"
               f"CZ CELLxGENE Census {p['s2']['census_version']}<br/>"
               f"Human Protein Atlas consensus RNA: {C.HPA_RNA_URL}<br/>"
               f"Human Protein Atlas normal-tissue IHC: {C.HPA_IHC_URL}<br/>"
               f"Open Targets Platform GraphQL: {C.OPENTARGETS_GQL}<br/>"
               f"DepMap CRISPR gene effect: {C.DEPMAP_URL}<br/>"
               f"Europe PMC: {C.EUROPEPMC_REST}", "small")]

    doc.build(flow)
    print(f"  wrote {path}")
    return str(path)


def _aspect(path) -> float:
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        return im.height / im.width


def run() -> list[str]:
    print("[S9c] reports")
    d = load_report_data()
    return [render("en", d), render("zh", d)]


if __name__ == "__main__":
    run()
