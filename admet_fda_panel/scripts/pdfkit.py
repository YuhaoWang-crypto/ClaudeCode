"""Shared ReportLab scaffolding for the Chinese-language ADMET report."""

from __future__ import annotations
import os
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph,
                                Spacer, Table, TableStyle, Image, KeepTogether)

CJK_TTC = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
FONT = "WQY"
FONT_B = "WQY"          # WenQuanYi Zen Hei ships a single weight

INK = colors.HexColor("#1f2933")
MUTED = colors.HexColor("#6b7680")
RULE = colors.HexColor("#c8d1d9")
ACCENT = colors.HexColor("#2f6f9f")
ACCENT_L = colors.HexColor("#eaf1f7")
BAND = colors.HexColor("#f5f7f9")
GOOD = colors.HexColor("#3d8b62")
WARN = colors.HexColor("#c0392b")


def register_fonts():
    pdfmetrics.registerFont(TTFont(FONT, CJK_TTC, subfontIndex=0))
    pdfmetrics.registerFontFamily(FONT, normal=FONT, bold=FONT,
                                  italic=FONT, boldItalic=FONT)


# --------------------------------------------------------------------------- #
def styles():
    base = dict(fontName=FONT, textColor=INK)
    return {
        "title": ParagraphStyle("title", **base, fontSize=22, leading=30,
                                spaceAfter=6),
        "subtitle": ParagraphStyle("subtitle", fontName=FONT, fontSize=12,
                                   leading=18, textColor=MUTED, spaceAfter=4),
        "meta": ParagraphStyle("meta", fontName=FONT, fontSize=9, leading=15,
                               textColor=MUTED),
        "h1": ParagraphStyle("h1", **base, fontSize=15.5, leading=22,
                             spaceBefore=16, spaceAfter=7),
        "h2": ParagraphStyle("h2", **base, fontSize=12, leading=18,
                             spaceBefore=11, spaceAfter=5),
        "h3": ParagraphStyle("h3", **base, fontSize=10.5, leading=16,
                             spaceBefore=8, spaceAfter=3),
        "body": ParagraphStyle("body", **base, fontSize=9.6, leading=15.6,
                               spaceAfter=6, alignment=0),
        "bullet": ParagraphStyle("bullet", **base, fontSize=9.6, leading=15.4,
                                 leftIndent=11, bulletIndent=2, spaceAfter=3),
        "caption": ParagraphStyle("caption", fontName=FONT, fontSize=8.3,
                                  leading=12.6, textColor=MUTED, spaceBefore=3,
                                  spaceAfter=9),
        "tcap": ParagraphStyle("tcap", **base, fontSize=9.2, leading=14,
                               spaceBefore=9, spaceAfter=4),
        "cell": ParagraphStyle("cell", fontName=FONT, fontSize=7.4, leading=10.2,
                               textColor=INK),
        "cellh": ParagraphStyle("cellh", fontName=FONT, fontSize=7.4,
                                leading=10.2, textColor=colors.white),
        "note": ParagraphStyle("note", fontName=FONT, fontSize=8.3, leading=13,
                               textColor=MUTED, spaceAfter=6),
        "ref": ParagraphStyle("ref", fontName=FONT, fontSize=8.5, leading=13.4,
                              textColor=INK, leftIndent=13, firstLineIndent=-13,
                              spaceAfter=4),
    }


# --------------------------------------------------------------------------- #
class Doc(BaseDocTemplate):
    def __init__(self, path, title, **kw):
        super().__init__(path, pagesize=A4, title=title,
                         author="ADMET 计算毒理分析", leftMargin=18 * mm,
                         rightMargin=18 * mm, topMargin=18 * mm,
                         bottomMargin=17 * mm, **kw)
        frame = Frame(self.leftMargin, self.bottomMargin,
                      self.width, self.height, id="body")
        self.addPageTemplates([
            PageTemplate(id="cover", frames=[frame]),
            PageTemplate(id="main", frames=[frame], onPage=self._decorate),
        ])
        self._doc_title = title

    def _decorate(self, canv, doc):
        canv.saveState()
        canv.setFont(FONT, 7.6)
        canv.setFillColor(MUTED)
        canv.drawString(doc.leftMargin, A4[1] - 12 * mm, self._doc_title)
        canv.setStrokeColor(RULE)
        canv.setLineWidth(0.5)
        canv.line(doc.leftMargin, A4[1] - 13.6 * mm,
                  A4[0] - doc.rightMargin, A4[1] - 13.6 * mm)
        canv.line(doc.leftMargin, 13.5 * mm, A4[0] - doc.rightMargin, 13.5 * mm)
        canv.drawRightString(A4[0] - doc.rightMargin, 10 * mm,
                             f"第 {canv.getPageNumber()} 页")
        canv.drawString(doc.leftMargin, 10 * mm,
                        "计算筛选结果 — 不替代 GLP 安全药理学试验")
        canv.restoreState()


def emphasize(text):
    """WenQuanYi ships a single weight, so <b> alone is invisible. Render
    emphasis as a distinct colour instead, keeping <b> for accessibility."""
    t = str(text)
    t = t.replace("<b>", '<font color="#0d4a75"><b>')
    t = t.replace("</b>", "</b></font>")
    return t


# --------------------------------------------------------------------------- #
def P(text, st):
    return Paragraph(emphasize(text), st)


def bullets(items, st, bstyle=None):
    return [Paragraph(emphasize(f"• {t}"), bstyle or st) for t in items]


def make_table(hdr, rows, widths, S, align=None, zebra=True,
               highlight=None, font_size=7.4):
    """highlight: dict {(row_idx, col_idx): color} applied to cell text."""
    cell = ParagraphStyle("c", parent=S["cell"], fontSize=font_size,
                          leading=font_size * 1.38)
    cellh = ParagraphStyle("ch", parent=S["cellh"], fontSize=font_size,
                           leading=font_size * 1.38)
    data = [[Paragraph(emphasize(h), cellh) for h in hdr]]
    for i, r in enumerate(rows):
        line = []
        for j, c in enumerate(r):
            st = cell
            if highlight and (i, j) in highlight:
                st = ParagraphStyle(f"hl{i}{j}", parent=cell,
                                    textColor=highlight[(i, j)])
            line.append(Paragraph(emphasize(c), st))
        data.append(line)

    cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT),
        ("GRID", (0, 0), (-1, -1), 0.25, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 3.6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3.6),
        ("TOPPADDING", (0, 0), (-1, -1), 3.2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.2),
    ]
    if zebra:
        for i in range(1, len(data)):
            if i % 2 == 0:
                cmds.append(("BACKGROUND", (0, i), (-1, i), BAND))
    if align:
        for col, a in align.items():
            cmds.append(("ALIGN", (col, 0), (col, -1), a))
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle(cmds))
    return t


def figure(path, width, S, caption):
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        w, h = im.size
    img = Image(path, width=width, height=width * h / w)
    return KeepTogether([img, P(caption, S["caption"])])


def callout(text, S, colour=ACCENT_L, border=ACCENT):
    p = Paragraph(emphasize(text), ParagraphStyle("co", parent=S["body"], fontSize=9,
                                       leading=14.4, spaceAfter=0))
    t = Table([[p]], colWidths=[None], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colour),
        ("LINEBEFORE", (0, 0), (0, -1), 2.2, border),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]))
    return t
