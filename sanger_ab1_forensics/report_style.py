#!/usr/bin/env python3
"""中文 PDF 报告的共用排版件：字体、配色、段落样式、表格 / 提示框 / 图片助手。

被 05_build_pdf_report.py（案例 1）与 case2_37HVPX/build_pdf.py（案例 2）共用。
"""
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, PageTemplate,
                                Paragraph, Table, TableStyle)

CJK_TTC = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
pdfmetrics.registerFont(TTFont("CJK", CJK_TTC, subfontIndex=0))
pdfmetrics.registerFontFamily("CJK", normal="CJK", bold="CJK",
                              italic="CJK", boldItalic="CJK")
MONO = "Courier"

INK = colors.HexColor("#1f2933")
MUTED = colors.HexColor("#7b8794")
BAD = colors.HexColor("#c0392b")
GOOD = colors.HexColor("#1e7a5a")
ACCENT = colors.HexColor("#41608f")
RULE = colors.HexColor("#d8dee4")
TINT = colors.HexColor("#f4f6f8")
BADTINT = colors.HexColor("#fdf1ef")
GOODTINT = colors.HexColor("#eef6f2")

PW, PH = A4
LM = RM = 45
CW = PW - LM - RM


def S(name, size=9.5, leading=None, color=INK, font="CJK", **kw):
    return ParagraphStyle(name, fontName=font, fontSize=size,
                          leading=leading or size * 1.65, textColor=color,
                          alignment=TA_LEFT, **kw)


ST = {
    "title":   S("title", 19, 26, INK, spaceAfter=2),
    "sub":     S("sub", 10, 15, MUTED, spaceAfter=0),
    "h1":      S("h1", 13.5, 19, ACCENT, spaceBefore=16, spaceAfter=7),
    "h2":      S("h2", 11, 16, INK, spaceBefore=11, spaceAfter=4),
    "body":    S("body", 9.5, 15.5, INK, spaceAfter=5),
    "small":   S("small", 8.3, 13, MUTED, spaceAfter=4),
    "cap":     S("cap", 8.3, 12.5, MUTED, spaceBefore=4, spaceAfter=10),
    "cell":    S("cell", 8.6, 12.5, INK),
    "cellc":   S("cellc", 8.6, 12.5, INK),
    "head":    S("head", 8.6, 12.5, colors.white),
    "mono":    S("mono", 7.2, 10.5, INK, font=MONO),
    "monos":   S("monos", 6.6, 9.2, INK, font=MONO),
    "verdict": S("verdict", 11, 17, BAD),
}


# WenQuanYi Zen Hei 缺这几个字形，ReportLab 会静默丢弃或画成方框
# （"−0.2" 渲染成 "0.2"，"10⁻²³" 渲染成 "10■²³"）。所有文本统一过一遍替换。
_SUB = {"−": "-", "⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
        "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}


def sane(t):
    if not isinstance(t, str):
        return t
    for a, b in _SUB.items():
        t = t.replace(a, b)
    return t


def P(t, s="body"):
    return Paragraph(sane(t), ST[s])


def em(t, c=BAD):
    return f'<font color="#{c.hexval()[2:]}">{t}</font>'


def table(rows, widths, head=True, zebra=True, align_center=(), font_size=None):
    """rows[0] 为表头。单元格内容可为 str（自动包成 Paragraph）或 Flowable。"""
    cell_st = ST["cell"]
    if font_size:
        cell_st = S("cell_s", font_size, font_size * 1.45, INK)
    data = []
    for r, row in enumerate(rows):
        out = []
        for cell in row:
            if isinstance(cell, str):
                out.append(Paragraph(sane(cell),
                                     ST["head"] if (head and r == 0) else cell_st))
            else:
                out.append(cell)
        data.append(out)
    cmds = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
        ("BOX", (0, 0), (-1, -1), 0.6, RULE),
    ]
    if head:
        cmds += [("BACKGROUND", (0, 0), (-1, 0), ACCENT),
                 ("LINEBELOW", (0, 0), (-1, 0), 0, colors.white)]
    if zebra:
        for i in range(1 if head else 0, len(rows)):
            if (i - (1 if head else 0)) % 2 == 1:
                cmds.append(("BACKGROUND", (0, i), (-1, i), TINT))
    for c in align_center:
        cmds.append(("ALIGN", (c, 0), (c, -1), "CENTER"))
    t = Table(data, colWidths=widths, repeatRows=1 if head else 0)
    t.setStyle(TableStyle(cmds))
    return t


def callout(lines, color=BAD, tint=BADTINT):
    head = S("verdict_c", 11, 17, color)          # 标题颜色跟随强调色
    inner = [[Paragraph(sane(l), head if i == 0 else ST["body"])]
             for i, l in enumerate(lines)]
    t = Table(inner, colWidths=[CW - 26])
    t.setStyle(TableStyle([
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 1),
    ]))
    outer = Table([[t]], colWidths=[CW])
    outer.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), tint),
        ("LINEBEFORE", (0, 0), (0, -1), 3, color),
        ("LEFTPADDING", (0, 0), (-1, -1), 13), ("RIGHTPADDING", (0, 0), (-1, -1), 13),
        ("TOPPADDING", (0, 0), (-1, -1), 10), ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
    ]))
    return outer


def img(path, width=CW):
    from PIL import Image as PILImage
    w, h = PILImage.open(path).size
    im = Image(path, width=width, height=width * h / w)
    im.hAlign = "CENTER"
    return im


def codebox(text, style="mono"):
    t = Table([[Paragraph(sane(text), ST[style])]], colWidths=[CW])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), TINT),
        ("BOX", (0, 0), (-1, -1), 0.6, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 9), ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 8), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
    return t


def numbered(items):
    """带编号的说明条目，返回 flowable 列表。"""
    from reportlab.platypus import Spacer
    out = []
    for n, t in enumerate(items, 1):
        out.append(table([[f'<font color="#{ACCENT.hexval()[2:]}">{n}</font>', t]],
                         [20, CW - 20], head=False, zebra=False))
        out.append(Spacer(1, 3))
    return out


def make_doc(path, title, header, footer="分析：Claude Code   ·   生成日期 2026-08-07"):
    """返回配好首页/内页模板的 BaseDocTemplate。首页无页眉。"""
    def decorate(canvas, doc, first=False):
        canvas.saveState()
        if not first:
            canvas.setFont("CJK", 7.5)
            canvas.setFillColor(MUTED)
            canvas.drawString(LM, PH - 32, header)
            canvas.setStrokeColor(RULE)
            canvas.setLineWidth(0.5)
            canvas.line(LM, PH - 38, PW - RM, PH - 38)
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(LM, 40, PW - RM, 40)
        canvas.setFont("CJK", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(LM, 30, footer)
        canvas.drawRightString(PW - RM, 30, f"第 {doc.page} 页")
        canvas.restoreState()

    doc = BaseDocTemplate(path, pagesize=A4, leftMargin=LM, rightMargin=RM,
                          topMargin=52, bottomMargin=52,
                          title=title, author="Claude Code")
    doc.addPageTemplates([
        PageTemplate(id="first", frames=[Frame(LM, 52, CW, PH - 52 - 40, id="first")],
                     onPage=lambda c, d: decorate(c, d, first=True)),
        PageTemplate(id="rest", frames=[Frame(LM, 52, CW, PH - 52 - 62, id="rest")],
                     onPage=decorate),
    ])
    return doc
