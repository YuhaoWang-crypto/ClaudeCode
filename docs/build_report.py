import markdown, pathlib, base64, re, mimetypes

FIGDIR = pathlib.Path("docs/figures")
src = pathlib.Path("docs/obesity_program_report_zh.md").read_text(encoding="utf-8")
body = markdown.markdown(src, extensions=["tables","fenced_code","sane_lists","attr_list","md_in_html"], output_format="html5")

# inline <img src="figures/NAME"> as base64 data URIs
def inline(m):
    src_attr = m.group(1)
    name = src_attr.split("/")[-1]
    p = FIGDIR / name
    if not p.exists():
        return m.group(0)
    b64 = base64.b64encode(p.read_bytes()).decode()
    mime = mimetypes.guess_type(name)[0] or "image/png"
    return f'src="data:{mime};base64,{b64}"'
body = re.sub(r'src="(figures/[^"]+)"', inline, body)

CSS = """
:root{ --ink:#1a2230; --muted:#5a6472; --line:#e3e8ef; --accent:#0b6b5e; --accent2:#0d47a1; --bg:#fff; --code:#eef2f6; }
*{ box-sizing:border-box; } html{ -webkit-print-color-adjust:exact; print-color-adjust:exact; }
body{ font-family:"Noto Sans","Helvetica Neue",Arial,"WenQuanYi Zen Hei","Microsoft YaHei",system-ui,sans-serif; color:var(--ink); background:var(--bg); line-height:1.62; margin:0; }
.wrap{ max-width:880px; margin:0 auto; padding:56px 48px 80px; }
h1{ font-size:26px; line-height:1.25; margin:0 0 6px; padding-bottom:14px; border-bottom:3px solid var(--accent); }
h2{ font-size:19px; margin:34px 0 10px; padding-left:11px; border-left:4px solid var(--accent); color:#0e2a45; }
h3{ font-size:15.5px; margin:22px 0 8px; color:#22303f; }
p{ margin:9px 0; } a{ color:var(--accent2); text-decoration:none; }
blockquote{ margin:14px 0; padding:10px 16px; background:#f5f9f8; border-left:4px solid var(--accent); color:var(--muted); font-size:13.5px; border-radius:0 6px 6px 0; }
code{ background:var(--code); padding:1.5px 5px; border-radius:4px; font-size:12px; font-family:"WenQuanYi Zen Hei Mono",ui-monospace,Menlo,Consolas,monospace; word-break:break-all; }
pre{ background:#0e2233; color:#e6edf3; padding:12px 14px; border-radius:8px; overflow:auto; font-size:11.5px; line-height:1.5; }
pre code{ background:none; padding:0; color:inherit; word-break:normal; white-space:pre; }
ul,ol{ margin:8px 0 8px 4px; padding-left:22px; } li{ margin:3px 0; }
hr{ border:none; border-top:1px solid var(--line); margin:30px 0; }
table{ border-collapse:collapse; width:100%; margin:14px 0; font-size:12.5px; }
th,td{ border:1px solid var(--line); padding:7px 10px; text-align:left; vertical-align:top; }
thead th{ background:#0e2a45; color:#fff; font-weight:600; border-color:#0e2a45; }
tbody tr:nth-child(even){ background:#f7fafc; }
strong{ color:#0b3d2e; }
img{ display:block; max-width:100%; height:auto; margin:12px auto 6px; border:1px solid var(--line); border-radius:8px; box-shadow:0 2px 10px rgba(0,0,0,.08); background:#0b0f14; }
.wrap>h1:first-child{ margin-top:0; }
@page{ size:A4; margin:15mm 13mm; }
@media print{ .wrap{ padding:0; max-width:none; } h2,h3{ page-break-after:avoid; } table,pre,blockquote,img{ page-break-inside:avoid; } }
"""
doc = f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>肥胖多靶点多模态药物开发 — 完整报告</title><style>{CSS}</style></head><body><div class="wrap">{body}</div></body></html>'
out = pathlib.Path("docs/obesity_report.html")
out.write_text(doc, encoding="utf-8")
print("HTML bytes:", out.stat().st_size, "| imgs inlined:", doc.count("data:image"))
