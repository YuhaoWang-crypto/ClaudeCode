#!/usr/bin/env python3
"""生成案例 2（37HVPX_1_D12）的中文 PDF 分析报告。

排版件来自上层目录的 report_style.py。
用法: cd case2_37HVPX && python3 build_pdf.py <out.pdf>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from reportlab.platypus import (KeepTogether, NextPageTemplate, PageBreak,
                                Paragraph, Spacer, Table, TableStyle)

from report_style import (ACCENT, BAD, CW, GOOD, GOODTINT, MUTED, RULE, ST,
                          TINT, P, callout, codebox, em, img, make_doc,
                          numbered, table)


def build(out):
    doc = make_doc(out, "纳米孔测序数据核查报告 37HVPX_1_D12",
                   "纳米孔测序数据核查报告 · 37HVPX_1_D12")
    E = []
    A = E.append

    # ══ 封面 ═══════════════════════════════════════════════════════════
    A(Spacer(1, 6))
    A(P("测序数据核查报告", "title"))
    A(P("样品 37HVPX_1_D12 &nbsp;·&nbsp; 双价 VHH CAR 慢病毒载体", "sub"))
    A(Spacer(1, 14))

    A(table([
        ["项目", "内容"],
        ["文件名", '<font name="Courier" size="8">37HVPX_1_D12.ab1</font>'],
        ["文件格式", "ABIF（.ab1）—— 但只有 <b>10</b> 个标签，仅够查看器打开的最小集"],
        ["仪器元数据", em("完全没有（无机型、染料、泳道、日期、运行方法）")],
        ["原始道 DATA1-4", em("不存在")],
        ["读长", "<b>4580 bp</b>（毛细管电泳物理上限约 1000–1200 bp）"],
        ["质量值", "Q4–50，均值 10.6，中位数 10，众数 9"],
        ["内容", "HIV-1 慢病毒转移载体表达盒 + 二代双价 VHH CAR"],
    ], [88, CW - 88]))

    A(Spacer(1, 15))
    A(callout([
        "结论一：不是 Sanger，也不是纳米孔的原始读段。",
        "4580 bp 远超毛细管电泳上限；但文件声称 Q9–12（≈10% 错误率），"
        "而序列在有公共参考的地方几乎零错误。<b>质量值和序列互相矛盾。</b>",
    ]))
    A(Spacer(1, 8))
    A(callout([
        "结论二：这份文件比案例 1 干净得多，性质不同。",
        "没有伪造仪器元数据，没有他人样品的原始道，色谱图是透明的方波渲染"
        "（峰高 = 19.7 × QV，r = 0.9965）。这像一次<b>诚实的格式转换</b>"
        "（FASTQ / FASTA → .ab1 供 SnapGene 查看），可疑的只有那条质量串。",
    ], ACCENT, TINT))
    A(Spacer(1, 8))
    A(callout([
        "结论三：序列内容本身合理，但无法自证。",
        "完整的二代 BBζ 双价 VHH CAR 表达盒，1488 nt 阅读框无移码、无提前终止，"
        "元件齐全、边界准确。作为一张<b>图纸</b>没问题；但它证明不了手上的质粒真的长这样。",
    ], GOOD, GOODTINT))

    A(Spacer(1, 15))
    A(P("本报告结构", "h2"))
    A(table([
        ["一", "这是什么文件：与案例 1 对比 + 色谱图渲染证据"],
        ["二", "质量值与序列对不上：滑窗一致性 + 两条判据"],
        ["三", "一处确实像纳米孔：3' 端 36 bp 折返"],
        ["四", "序列内容：载体图谱、CAR 结构域、VHH 特征"],
        ["五", "能不能用 + 该要什么数据"],
    ], [34, CW - 34], head=False, zebra=False))

    A(NextPageTemplate("rest"))
    A(PageBreak())

    # ══ 一 ═════════════════════════════════════════════════════════════
    A(P("一、这是什么文件", "h1"))
    A(P("文件头 4 字节仍是 <font name=\"Courier\" size=\"8.4\">ABIF</font>，但内部与案例 1 的"
        "那份伪造 Sanger 谱图完全不是一回事："))
    A(Spacer(1, 4))
    A(table([
        ["", "案例 1　DAPA21", "案例 2　37HVPX_1_D12"],
        ["ABIF 标签数", "126", em("10（最小集）", ACCENT)],
        ["仪器元数据", em("完整伪造：3730xl / POP7 / BigDye v3 / 泳道 / 日期"),
         em("完全没有", ACCENT)],
        ["原始道 DATA1-4", em("有，但属于另一份 2023 年的哺乳载体样品"),
         em("没有", ACCENT)],
        ["读长", "735 bp", "4580 bp"],
        ["质量值", em("全部恒为 50"), "Q4–50，均值 10.6"],
        ["色谱图", em("等高等宽高斯，伪装成真实电泳"),
         em("方波渲染，4 采样点/碱基，峰高 = 19.7 × QV", ACCENT)],
        ["性质判断", em("伪造的 Sanger 谱图"), em("诚实的格式转换", GOOD)],
    ], [72, (CW - 72) / 2, (CW - 72) / 2]))

    A(Spacer(1, 10))
    A(KeepTogether([
        img("fig1_render_evidence.png"),
        P("图 1　① 每个碱基固定占 4 个采样点，波形是方波而非电泳峰。② 把每个碱基的峰高对它的"
          "质量值作图，落在一条直线上（峰高 ≈ 19.7 × QV + 14，r = 0.9965）—— 色谱图是直接由"
          "质量串画出来的，不含任何独立的测量信息。", "cap")]))

    A(PageBreak())

    # ══ 二 ═════════════════════════════════════════════════════════════
    A(P("二、质量值与序列对不上", "h1"))
    A(P("沿读段每 500 bp 取一个窗口，把「质量值隐含的错误率」与「实测同公共参考的不一致率」"
        "并排比较："))
    A(Spacer(1, 4))
    A(KeepTogether([
        img("fig2_quality_vs_accuracy.png", 468),
        P("图 2　黄色区间是自定义 VHH / CAR 序列，公共库里没有确切参考，"
          "低一致性只说明「查不到」，不代表测序错误。区间之外的每个窗口，"
          "实测一致性都是 99–100%。", "cap")]))

    A(table([
        ["读段窗口", "平均 QV", "QV 隐含错误率", "BLAST 最佳一致性", "最佳命中"],
        ["1–500", "16.4", "5.0%", em("464/464 = 100%", GOOD), "合成构建体"],
        ["501–1000", "10.2", "9.6%", "267/311 = 86%", "CAR 构建体（自定义区）"],
        ["1001–1500", "9.8", "10.6%", "271/345 = 79%", "CAR 构建体（自定义区）"],
        ["1501–2000", "9.1", "12.5%", "180/229 = 79%", "双峰驼免疫球蛋白 mRNA"],
        ["2001–2500", "9.0", "12.5%", "287/291 = 99%", "载体 pGZ-DSB-SDSA"],
        ["2501–3000", "9.0", "12.5%", em("500/500 = 100%", GOOD), "人基因组序列"],
        ["3001–3500", "9.8", "10.6%", em("500/500 = 100%", GOOD), "载体 pGA1backbone"],
        ["3501–4000", "10.5", "9.0%", "439/440 = 99%", "载体 pLVX.TRE3G"],
        ["4001–4500", "11.2", "7.6%", em("500/500 = 100%", GOOD), "合成构建体"],
    ], [62, 48, 66, 96, CW - 62 - 48 - 66 - 96], align_center=(1, 2, 3), font_size=8.0))

    A(P("两条判据", "h2"))
    for f in numbered([
        "三个独立的 500 bp 窗口对公共参考 <b>500/500 全对</b>。若真是 Q10（10% 每碱基错误率），"
        "单个 500 bp 窗口完全无误的概率约 "
        "<font name=\"Courier\" size=\"8.4\">0.9^500 ~ 1e-23</font>，三个窗口同时如此更无可能。",
        "1488 nt 的 CAR 阅读框译出 495 aa，<b>无移码、无提前终止</b>，结构域边界完全正确。"
        "纳米孔的错误以插入缺失为主，Q10 的原始读段会把这么长的阅读框打得千疮百孔。",
    ]):
        A(f)
    A(Spacer(1, 4))
    A(callout([
        "所以：这条序列是打磨过的（consensus / 组装 / 或设计图），"
        "那条 Q9–12 的质量串不是它的质量。",
    ], BAD))

    A(PageBreak())

    # ══ 三 + 四 ════════════════════════════════════════════════════════
    A(P("三、一处确实像纳米孔", "h1"))
    A(P("读段 3' 末端有 <b>36 bp 折返</b>：在回文位点 "
        "<font name=\"Courier\" size=\"8.4\">CAGCTG</font> 处掉头，"
        "接着读回前面序列的反向互补。"))
    A(P("载体正义方向最后 60 bp（竖线为掉头处）：", "small"))
    A(codebox("...GACCAATGACTTACAAGGCAGCTG <b>|</b> CCTTGTAAGTCATTGGTCCCTCATATCTCCTCCTCC<br/>"
              "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
              "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
              "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
              "&nbsp;&nbsp;&nbsp;= revcomp(GGAGGAGGAGATATGAGGGACCAATGACTTACAAGG)"))
    A(Spacer(1, 6))
    A(P("这是纳米孔建库时双链未分开、测完一条链折回去测互补链的典型假象。"
        "所以底层数据很可能<b>确实是纳米孔的</b>—— 但这份 .ab1 是成品，不是证据。"))

    A(P("四、序列内容：二代双价 VHH CAR 慢病毒转移载体", "h1"))
    A(P("以下坐标为载体正义方向（即文件内序列的反向互补），全长 4580 bp。"))
    A(Spacer(1, 4))
    A(KeepTogether([
        img("fig3_vector_map.png", 468),
        P("图 3　上：慢病毒转移载体表达盒，从 RRE 一直读到 3'ΔU3 LTR 起始，"
          "完整覆盖整个「有效载荷」区。下：CAR 的结构域组成。", "cap")]))

    A(PageBreak())
    A(table([
        ["位置", "元件"],
        ["1 .. ~850", "HIV-1 env 片段 / <b>RRE</b>（Rev 应答元件，茎环 136–171）"],
        ["861 .. 899", "<b>cPPT / CTS</b>（中央多聚嘌呤序列）"],
        ["1176 .. ~2350", "<b>人 EF-1α 启动子</b>（完整，含内含子；供体 1400，受体 2341）"],
        ["2358 .. 2401", "多克隆位点 NheI–AgeI–BglII–XhoI–SacI–HindIII–EcoRI"],
        ["2402 .. 2416", "Kozak <font name=\"Courier\" size=\"8\">GCCGCCACC</font> + 起始密码子"],
        ["<b>2411 .. 3898</b>", em("<b>CAR 阅读框，1488 nt / 495 aa</b>", BAD)],
        ["3890 .. 3904", "终止密码子 TAA + XbaI"],
        ["3919 .. ~4500", "<b>WPRE</b>（土拨鼠肝炎病毒转录后调控元件）"],
        ["4510 .. 4544", "<b>3'ΔU3 LTR</b> 起始"],
        ["4545 .. 4580", em("折返假象（36 bp）")],
    ], [88, CW - 88]))

    A(P("CAR 结构域", "h1"))
    A(table([
        ["氨基酸位置", "结构域", "说明"],
        ["1 – 21", "CD8α 信号肽", "MALPVTALLLPLALLLHAARP"],
        ["22 – 142", em("VHH1", BAD), "骆驼科单域抗体"],
        ["143 – 147", "(G₄S) 连接子", "串联双价"],
        ["148 – 260", em("VHH2", BAD), "骆驼科单域抗体"],
        ["273 – 317", "CD8α 铰链", "TTTPAPRPPTPAPTIASQPLSLRPEACRPAAGGAVHTRGLDFACD"],
        ["318 – 341", "CD8α 跨膜", "IYIWAPLAGTCGVLLLSLVITLYC"],
        ["342 – 383", em("4-1BB（CD137）共刺激域", GOOD),
         "KRGRKKLLYIFKQPFMRPVQTTQEEDGCSCRFPEEEEGGCEL"],
        ["384 – 495", em("CD3ζ 信号域", GOOD), "RVKFSRSADAPAYQQGQNQLYNELNLGRR…"],
    ], [66, 108, CW - 66 - 108], font_size=8.0))
    A(P("即标准的 <b>二代 BBζ CAR</b>。共刺激域是 <b>4-1BB</b> 而非 CD28"
        "（CD28 的胞内段 RSKRSRLLHSDYMNMTPRRPGPTRKHYQPYAPPRDFAAYRS 未出现）。", "body"))

    A(P("两个 VHH 的骆驼科特征", "h2"))
    A(table([
        ["", "VHH1（aa 22–142）", "VHH2（aa 148–260）"],
        ["FR2 hallmark", em("WYRQ"), em("WGRQ")],
        ["CDR1", "RYTSSGCNVG　<font color=\"#c0392b\">（含 Cys）</font>", "GFAFSNYAMS"],
        ["CDR3", "NIAPPYTTHCASRRF　<font color=\"#c0392b\">（含 Cys）</font>",
         "AENVDCYGGYCYRANY　<font color=\"#c0392b\">（含 2 个 Cys）</font>"],
        ["额外二硫键", "CDR1 – CDR3", "CDR3 内部"],
    ], [72, (CW - 72) / 2, (CW - 72) / 2], font_size=8.2))
    A(P("常规抗体重链 FR2 此处是 <font name=\"Courier\" size=\"8.4\">WVRQ</font>；"
        "两个结构域都带 VHH 特有的替换，并各有一对额外二硫键 —— 是真正的纳米抗体，不是 scFv。",
        "small"))
    A(P("blastp 查不到确切来源：最好的命中是各种羊驼纳米抗体，只有 70–75% 一致。"
        "属于<b>自有 / 未公开的 VHH</b>，从序列判断不出靶点。", "body"))

    A(P("五、能不能用", "h1"))
    A(P("<b>序列内容本身合理。</b>阅读框干净、元件齐全、结构域边界准确、Kozak 正确 —— "
        "作为一张图纸没问题。"))
    A(P("<b>但它不能自证。</b>序列是打磨过的、质量串又和序列对不上，"
        "这份 .ab1 证明不了手上的质粒真的长这样。"))
    A(Spacer(1, 4))
    A(P("若对方做的是纳米孔全质粒测序（常见服务），完整交付应包含：", "body"))
    for f in numbered([
        "<b>原始 reads</b> —— FASTQ（或 POD5 / FAST5）",
        "<b>比对结果</b> —— BAM + 覆盖深度曲线",
        "<b>组装 / consensus FASTA</b> + 组装工具与版本",
        "<b>CAR 插入片段（2411–3898）上的覆盖深度</b>，以及有没有 reads 完整跨过整个插入片段",
    ]):
        A(f)
    A(P("有了第 2 和第 4 项才谈得上「验证」；现在这份只是把结论渲染成了一张图。", "body"))

    # ══ 附录 ═══════════════════════════════════════════════════════════
    A(P("附录 A　方法", "h1"))
    A(table([
        ["步骤", "做法"],
        ["解析", "Biopython <font name=\"Courier\" size=\"8\">SeqIO(\"abi\")</font> "
                 "读取 ABIF 标签；直接比对 PLOC 间距、DATA9-12 长度与碱基数"],
        ["渲染判定", "取每个碱基在 PLOC 处的峰高，与 PCON 的质量值做线性回归"],
        ["滑窗一致性", "沿读段每 500 bp 取一段，各自提交 NCBI BLAST（blastn 对 nt）"],
        ["折返检验", "读段与自身反向互补做局部比对（Biopython PairwiseAligner）"],
        ["阅读框", "定位 Kozak + ATG 后翻译至终止密码子，比对已知结构域基序"],
        ["VHH 溯源", "两个可变域分别提交 blastp 对 nr"],
    ], [72, CW - 72]))
    A(P("依赖：biopython · numpy · matplotlib · reportlab；BLAST 经 NCBI URL API。"
        "序列同源判定依赖 NCBI 公共数据库 2026-08 版本。", "small"))

    A(P("附录 B　CAR 蛋白序列（495 aa）", "h1"))
    prot = "".join(l.strip() for l in open("CAR_protein_495aa.fasta")
                   if not l.startswith(">"))
    rows = []
    for i in range(0, len(prot), 60):
        rows.append([Paragraph(f'<font color="#{MUTED.hexval()[2:]}">{i + 1:>4}</font>',
                               ST["monos"]),
                     Paragraph(" ".join(prot[i:i + 60][j:j + 10] for j in range(0, 60, 10)),
                               ST["monos"])])
    t = Table(rows, colWidths=[30, CW - 30])
    t.setStyle(TableStyle([
        ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6), ("BACKGROUND", (0, 0), (-1, -1), TINT),
        ("BOX", (0, 0), (-1, -1), 0.6, RULE)]))
    A(t)
    A(Spacer(1, 6))
    A(P("完整核酸序列见同目录 "
        "<font name=\"Courier\" size=\"8\">37HVPX_D12_vector_sense_4580bp.fasta</font>"
        "（载体正义方向）与 <font name=\"Courier\" size=\"8\">CAR_orf_1488nt.fasta</font>。",
        "small"))

    doc.build(E)
    print("wrote", out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "37HVPX_D12_测序数据核查报告.pdf")
