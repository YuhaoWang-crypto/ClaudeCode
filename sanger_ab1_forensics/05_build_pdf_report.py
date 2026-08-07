#!/usr/bin/env python3
"""生成案例 1（DAPA21）的中文 PDF 分析报告。

排版件来自 report_style.py。
用法: python3 05_build_pdf_report.py <out.pdf>
"""
import sys

from reportlab.lib import colors
from reportlab.platypus import (KeepTogether, NextPageTemplate, PageBreak,
                                Paragraph, Spacer, Table, TableStyle)

from report_style import (ACCENT, BAD, CW, GOOD, GOODTINT, INK, MUTED, RULE,
                          ST, TINT, P, callout, em, img, make_doc, table)


def build(out):
    doc = make_doc(out, "Sanger 测序数据质控取证报告 DAPA21",
                   "Sanger 测序数据质控取证报告 · DAPA21")

    E = []
    A = E.append

    # ══ 封面区 ═════════════════════════════════════════════════════════
    A(Spacer(1, 6))
    A(P("Sanger 测序数据质控取证报告", "title"))
    A(P("样品 DAPA21 &nbsp;·&nbsp; dapA 缺失位点验证", "sub"))
    A(Spacer(1, 14))

    A(table([
        ["项目", "内容"],
        ["文件名", '<font name="Courier" size="7.6">DAPA21_TSM20260319057101751_'
                   '202603190571BAN171_E08.1_1.ab1</font>'],
        ["测序服务商", "擎科生物 TSINGKE（头文件 RGNm1 = TSINGKE）"],
        ["仪器 / 耗材", "ABI 3730xl · POP7 聚合物 · BigDye Terminator v3 · KB 1.4.1.8"],
        ["运行信息", "2026-03-20 20:03:14 · 泳道 8 · 96 孔板 E08"],
        ["报告读长", "735 bp"],
    ], [95, CW - 95]))

    A(Spacer(1, 16))
    A(callout([
        "结论一：这份 .ab1 的色谱图是合成的，不是真实电泳信号。",
        "7 项独立检验全部异常；文件内的「原始道」还属于另一份 2023 年的哺乳动物载体样品。"
        "该文件<b>不能</b>作为克隆验证凭证使用。",
    ]))
    A(Spacer(1, 9))
    A(callout([
        "结论二：报告的 735 bp 序列本身是真实生物学序列，内容为 dapA 无痕敲除的连接点。",
        "缺失 878 bp（占 dapA CDS 的 99.9%），邻近的 bamC / gcvR / purC 完好，"
        "无 FRT / loxP / 抗性盒残留。侧翼序列 100% 匹配 BL21 与 K-12 MG1655（与 Nissle 1917 差 7 处），"
        "故背景为 BL21 / K-12 而非 EcN。",
    ], GOOD, GOODTINT))

    A(Spacer(1, 16))
    A(P("结论速览", "h1"))
    A(table([
        ["问题", "结论"],
        ["这份 .ab1 合理吗？", em("不合理。色谱图为合成/重绘，7 项检验全部异常。")],
        ["判读序列可信吗？", "序列内容像是真实的（完美比对到 E. coli），"
                             "但<b>该文件无法为它提供任何证据支持</b>。"],
        ["序列是什么？", "ΔdapA 无痕缺失连接点，跨 gcvR ‖ bamC 断点。"],
        ["菌株背景？", em("侧翼 100% 匹配 BL21 与 K-12 MG1655，与 Nissle 1917 差 7 处"
                          " —— 不是 EcN 背景。", GOOD)],
        ["缺失多大？", em("878 bp，占 dapA CDS 的 99.9%，无标记、无疤痕。", GOOD)],
        ["原始道是什么？", em("一条完全不相关的哺乳动物慢病毒载体"
                              "（pCDH-CMV-MCS-EF1-Puro 类骨架的 EF1α 启动子区）。")],
    ], [110, CW - 110]))

    A(Spacer(1, 18))
    A(P("本报告结构", "h2"))
    A(table([
        ["一", "真实性取证：7 项检验 + 元数据矛盾 + 波形对比图"],
        ["二", "文件内含三条互不相关的序列：比对统计与 BLAST 鉴定"],
        ["三", "序列的具体生物学分析：dapA 断点定位、位点示意图、解读与建议"],
        ["附录", "方法与脚本清单；报告的 735 bp 序列原文"],
    ], [34, CW - 34], head=False, zebra=False))

    A(NextPageTemplate("rest"))
    A(PageBreak())

    # ══ 第一部分 ═══════════════════════════════════════════════════════
    A(P("一、真实性取证：为什么判定色谱图是合成的", "h1"))
    A(P("真实的毛细管电泳数据必然带有几项物理特征：峰高随迁移时间指数衰减、峰形因扩散与"
        "拖尾而不对称、峰间距随读长漂移、基线经扣除后有正负涨落、读段两端质量值下降。"
        "合成或重绘的色谱图不具备这些性质。以下 7 项检验相互独立。"))
    A(Spacer(1, 5))

    A(table([
        ["#", "检验项", "本文件实测", "真实数据应有"],
        ["T1", "峰间距一致性", em("标准差 0.000 scan（恒为 13）"),
         "随读长单调漂移，sd &gt; 1"],
        ["T2", "分析道长度", em("9555 = 735 × 13，整除"), "与碱基数无整除关系"],
        ["T3", "峰高分布", em("CV = 0.003；首 50 峰 / 末 50 峰 = 1.00；"
                              "735 个峰仅 11 种高度"), "CV 0.3–0.8；末端衰减数倍"],
        ["T4", "峰形对称性", em("峰心等距两点强度差 0.0057（近乎完美高斯）"),
         "&gt; 0.03（毛细管拖尾）"],
        ["T5", "基线特征", em("51.2% 的点恰好等于 0，零个负值"),
         "扣基线后必然出现正负涨落"],
        ["T6", "质量值 QV", em("735 个碱基全部 QV = 50，唯一取值 1 个"),
         "读段两端必然衰减"],
        ["T7", "内部数组自洽", em("PBAS/PCON/PLOC = 735，"
                                  "而 P1AM/P2AM/P2BA/P1WD = 252"), "全部相等"],
    ], [24, 74, CW - 24 - 74 - 128, 128], align_center=(0,)))
    A(P("检验代码见 <font name=\"Courier\" size=\"8\">01_ab1_qc_forensics.py</font>，"
        "原始输出见附带的 qc_forensics_report.txt。", "small"))

    A(Spacer(1, 10))
    A(KeepTogether([
        img("results/fig1_trace_comparison.png", 452),
        P("图 1　同一个 .ab1 文件内部的两套数据。①「分析道」DATA9-12 是软件写入的理想化高斯"
          "序列；②「原始道」DATA1-4 保留了真实电泳形态；③ 全长包络对比最直观——真实信号随"
          "片段增长指数衰减，合成信号 735 个碱基全程平坦。", "cap")]))
    A(PageBreak())

    A(P("元数据自相矛盾", "h2"))
    A(P("头部 <font name=\"Courier\" size=\"8.4\">RUND = 2026-03-20</font>（与订单号 "
        "TSM2026031905710175 一致），但同一文件里："))
    A(table([
        ["字段", "取值", "含义"],
        ["RunN1", "Run_ABI-PC_2023-03-13_20-03_0179", "运行名仍是 2023 年"],
        ["BCTS1", "2023-03-13 21:38:20 +08:00", "采集时间戳 2023 年"],
        ["CTID1 / CTNM1", "20230313-FY041-020-020-f", "容器 ID 2023 年"],
        ["SMED1", "Jan 01, 2023", "样品制备日期 2023 年"],
        ["TUBE1", "E2", "与文件名标注的 E08 孔不符"],
    ], [78, 200, CW - 78 - 200]))
    A(P("这些残留说明：该文件是<b>以一份 2023 年的旧 .ab1 为模板重写生成</b>的，"
        "而不是仪器直接输出。"))


    # ══ 第二部分 ═══════════════════════════════════════════════════════
    A(P("二、文件内含三条互不相关的序列", "h1"))
    A(P("原始道（DATA1-4）是文件中唯一未被改写的部分。对它做基线扣除，由数据自估 4×4 染料"
        "串扰矩阵后求逆解混，再逐通道寻峰，可独立反演出一条 898 nt 序列（因缺少 .mob 迁移率"
        "校正，约 10% 误差，足够做同源判定）。把它与文件内各序列做局部比对，并以 30 次序列"
        "随机打乱作对照："))
    A(Spacer(1, 4))
    A(table([
        ["比较对象", "局部比对分", "打乱对照", "Z 值"],
        ["原始道反演 vs 自身前 300 nt（阳性对照）", "600", "20 ± 2", em("+262.8", GOOD)],
        ["原始道反演 vs PBAS1 判读序列", "24", "21 ± 2", em("+1.3")],
        ["原始道反演 vs P2BA1 内嵌 252 nt", "19", "19 ± 1", em("−0.2")],
        ["PBAS1 vs P2BA1", "20", "19 ± 2", em("+0.3")],
    ], [CW - 210, 70, 70, 70], align_center=(1, 2, 3)))
    A(P("即：<b>原始波形与所报告的序列毫无关系</b>，比对分与随机打乱无异。", "body"))

    A(P("BLAST 鉴定（NCBI nt / nr）", "h2"))
    A(table([
        ["序列", "最佳命中", "一致性"],
        ["PBAS1（735 bp，报告序列）", "Escherichia coli Nissle 1917 染色体 CP171242.1",
         "728/735 = 99.05%<br/>E = 0"],
        ["原始道反演（898 nt）",
         "慢病毒/哺乳表达载体 pCDH-CMV-MCS-lox-EF1-Puro-lox、pLV、pFUSE、piggyBac 等同源骨架",
         "505/564 = 90%"],
        ["P2BA1（252 nt）", "无显著同源（符合「次选碱基」本应是噪音的性质）", "—"],
    ], [142, CW - 142 - 118, 118]))
    A(P("原始道对应的是 2023 年那份模板文件里的哺乳动物载体样品；判读序列则另有来源。", "small"))

    A(PageBreak())

    # ══ 第三部分 ═══════════════════════════════════════════════════════
    A(P("三、序列的具体生物学分析：dapA 无痕敲除", "h1"))
    A(P("把 PBAS1 反向互补后比对到 E. coli K-12 MG1655（NC_000913.3）。该位点野生型基因顺序为 "
        "purC – bamC(nlpB/dapX) – dapA – gcvR（均在负链）。读段呈<b>两臂分裂</b>："))
    A(Spacer(1, 4))
    A(table([
        ["读段片段", "对应 MG1655 坐标", "比对分"],
        ["右臂（bamC 侧）　读段 1–365", "2,598,517 – 2,598,881", "730"],
        ["左臂（gcvR 侧）　读段 365–735", "2,599,760 – 2,600,130", "742"],
    ], [CW - 210, 140, 70], align_center=(1, 2)))
    A(P("两臂在读段上重叠 1 bp（微同源 “T”）。由此定出：", "body"))
    A(Spacer(1, 3))
    A(callout([
        "缺失区间 = MG1655 2,598,882 – 2,599,759 = 878 bp",
        "占 dapA CDS（2,598,882–2,599,760，879 bp）的 99.9%；bamC、gcvR、purC 无一碱基落入缺失区。",
    ], GOOD, GOODTINT))
    A(Spacer(1, 10))

    A(table([
        ["基因", "MG1655 坐标", "长度", "落入缺失区"],
        ["purC", "2,596,905 – 2,597,618", "714 bp", "0 bp（0%）"],
        ["bamC (nlpB/dapX)", "2,597,831 – 2,598,865", "1035 bp", "0 bp（0%）"],
        [em("dapA"), em("2,598,882 – 2,599,760"), em("879 bp"), em("878 bp（99.9%）")],
        ["gcvR", "2,599,906 – 2,600,478", "573 bp", "0 bp（0%）"],
    ], [110, 150, 60, CW - 110 - 150 - 60], align_center=(2, 3)))

    A(Spacer(1, 12))
    A(KeepTogether([
        img("results/fig2_dapA_junction.png"),
        P("图 2　位点示意。上行为野生型真实坐标；下行按缺失折叠绘制，即 gcvR 一侧整体左移 878 bp "
          "与 bamC 侧直接相接——这才是样品分子的真实构型。735 bp 读段跨断点连续覆盖。", "cap")]))

    A(P("断点邻近序列", "h2"))
    A(P("小写 = bamC 侧　|　竖线 = 连接点　|　大写 = gcvR 侧", "small"))
    A(Table([[Paragraph(
        "...aggcgcgacttttgaacagagtaagccatcaaatctccctaaact <b>|</b> "
        "GGGCCATCCTCTGTGCAAACAAGTGTCTCAATGGTACGTTTGGTA...", ST["mono"])]],
        colWidths=[CW], style=TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), TINT),
            ("BOX", (0, 0), (-1, -1), 0.6, RULE),
            ("LEFTPADDING", (0, 0), (-1, -1), 9),
            ("RIGHTPADDING", (0, 0), (-1, -1), 9),
            ("TOPPADDING", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 8)])))

    A(PageBreak())
    A(P("解读", "h1"))
    for n, t in enumerate([
        "<b>dapA 被完整敲除。</b>878/879 bp 被删除，只余起始密码子的 1 个碱基，无残留可读框。",
        "<b>无痕（markerless）敲除。</b>连接处只有 1 bp 微同源，"
        "<b>没有</b> FRT(34 bp)、loxP(34 bp) 或任何抗性表达盒——"
        "符合 λ-Red + FLP/I-SceI 无痕或 CRISPR 辅助的无标记缺失策略。",
        "<b>邻近基因完好。</b>读段覆盖的 bamC 前 116 个密码子翻译后与参考 BamC 蛋白 "
        "116/116 = 100% 一致（MAYSVQKSRLAKVAGVSLVLLLAACSSDSRYKRQ…，含典型脂蛋白 lipobox "
        "“LAAC”）；该段对野生型 E. coli 基因组为 365/365 = 100%。gcvR 侧同样完整。"
        "说明缺失精确，未波及 BAM 外膜蛋白装配复合体组分与甘氨酸裂解系统调控子。",
        "<b>表型含义。</b>dapA 编码 4-羟基-四氢吡啶二羧酸合酶（DHDPS），是赖氨酸 / DAP 支路"
        "第一个限速步骤。ΔdapA 无法合成 meso-DAP，肽聚糖无法交联，"
        "<b>离开外源 DAP 即裂解死亡</b>——这正是活菌药物 / 接合供体菌标准的"
        "<b>营养缺陷型生物封闭（biocontainment）</b>设计。",
        "<b>菌株背景是 BL21 / K-12，不是 EcN。</b>把两臂分开比对（避开断点）：bamC 臂对 BL21 "
        "与 MG1655 均 365/365 = 100%，gcvR 臂均 371/371 = 100%；对 Nissle 1917 则分别是 "
        "362/365 与 367/371（共 7 处差异）。该位点 BL21 与 K-12 序列相同，故只能排除 EcN，"
        "不能区分这两者。",
        "<b>唯一含此断点的已发表基因组。</b>GenBank CP171242.1（“E. coli Nissle 1917”，2024 年"
        "提交）在同一位点同样缺失 dapA，读段对它 735/735 全长比对 99.05%——这是全长唯一能比上的"
        "基因组，但那 7 处差异说明它只是碰巧含同一连接点（两边都按“精确删掉 dapA CDS”设计，"
        "自然产生同一断点），并非本样品的背景。",
        "文件名中的 <b>“DAPA21”</b> 与此完全自洽——很可能是 dapA 敲除克隆 #21 的验证测序。",
    ], 1):
        A(table([[f'<font color="#{ACCENT.hexval()[2:]}">{n}</font>', t]],
                [20, CW - 20], head=False, zebra=False))
        A(Spacer(1, 3))

    A(P("建议", "h1"))
    for n, t in enumerate([
        "<b>向擎科索要真实原始 .ab1</b>，并说明本文件色谱图为合成、且原始道属于另一份 2023 年的"
        "哺乳载体样品。在拿到真实数据前，这份文件不能进实验记录、不能作为克隆验证凭证，"
        "更不能用于报批或论文补充材料。",
        "拿到新文件后先跑 <font name=\"Courier\" size=\"8.4\">01_ab1_qc_forensics.py</font> 复核。",
        "dapA 敲除的关键结论（878 bp 无痕缺失、两侧基因完整）本身很可能是对的，但需要"
        "<b>真实色谱图 + 双向测序</b>（正反向引物各一条，完整覆盖断点）来确认。",
        "补做表型验证：无 DAP 的 M9 / LB 平板上不生长、补加 DAP（50–100 µg/mL）恢复生长。",
    ], 1):
        A(table([[f'<font color="#{ACCENT.hexval()[2:]}">{n}</font>', t]],
                [20, CW - 20], head=False, zebra=False))
        A(Spacer(1, 3))

    A(PageBreak())
    A(P("附录 A　方法", "h1"))
    A(table([
        ["脚本", "作用"],
        ['<font name="Courier" size="8">01_ab1_qc_forensics.py</font>',
         "解析 ABIF 头部与数据块；常规 QC（读长 / QV / 信噪比 / 简并碱基）；7 项合成痕迹检验"],
        ['<font name="Courier" size="8">02_raw_trace_basecall.py</font>',
         "对原始道滑动分位数扣基线；按主导通道聚类自估染料串扰矩阵；求逆解混后逐通道寻峰反演序列"],
        ['<font name="Courier" size="8">03_dapA_junction.py</font>',
         "两臂分别做局部比对定位；求断点、微同源与缺失区间；计算各基因落入缺失区的比例"],
        ['<font name="Courier" size="8">04_report_figures.py</font>', "生成图 1、图 2"],
        ['<font name="Courier" size="8">05_build_pdf_report.py</font>', "生成本报告"],
    ], [150, CW - 150]))
    A(P("依赖：biopython · numpy · matplotlib · reportlab；03 需联网访问 NCBI E-utilities。"
        "比对使用 Biopython PairwiseAligner（local，match +2 / mismatch −3 / gap −10,−4）；"
        "同源鉴定使用 NCBI BLAST URL API（megablast 对 nt、blastx 对 nr）。", "small"))

    A(P("附录 B　报告的 735 bp 序列", "h1"))
    A(P("文件内 PBAS1 字段原文（FASTA 见 DAPA21_reported_read.fasta）：", "small"))
    seq = "".join(l.strip() for l in
                  open("results/DAPA21_reported_read.fasta") if not l.startswith(">"))
    rows = []
    for i in range(0, len(seq), 60):
        rows.append([Paragraph(f'<font color="#{MUTED.hexval()[2:]}">{i + 1:>4}</font>',
                               ST["monos"]),
                     Paragraph(" ".join(seq[i:i + 60][j:j + 10] for j in range(0, 60, 10)),
                               ST["monos"])])
    t = Table(rows, colWidths=[30, CW - 30])
    t.setStyle(TableStyle([
        ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6), ("BACKGROUND", (0, 0), (-1, -1), TINT),
        ("BOX", (0, 0), (-1, -1), 0.6, RULE)]))
    A(t)
    A(Spacer(1, 6))
    A(P("说明：本报告所有数值均由附带脚本对该 .ab1 文件直接计算得出，可复现。"
        "序列同源判定依赖 NCBI 公共数据库 2026-08 版本。", "small"))

    doc.build(E)
    print("wrote", out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "results/DAPA21_测序数据分析报告.pdf")
