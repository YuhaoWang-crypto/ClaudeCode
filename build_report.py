#!/usr/bin/env python3
"""Build the HTML demo report from a finished `demo_out/` run.

    python build_report.py                       # -> dist/neoantigen-selection-report.html
    python build_report.py --demo demo_out --out dist

Two outputs from one template, because they need different wrappers:

  --mode standalone  (default)  a complete HTML document -- doctype, <head>,
                                charset, viewport -- that opens correctly from
                                the filesystem and can be emailed around
  --mode artifact               body content only, for publishing as an
                                Artifact, where the host supplies the skeleton

Figures are inlined as base64 data URIs so a single file is the whole report.
Fonts are linked from Google Fonts; offline the page falls back to the declared
serif / sans / monospace stacks, which is why those stacks are real and not
decoration.
"""

from __future__ import annotations

import argparse
import base64
import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(ROOT, "report_html", "report.template.html")

FIGURES = {
    "__IMG_SUMMARY__": "summary.png",
    "__IMG_JUNCTIONS__": "junctions.png",
    "__IMG_BENCH__": "benchmark.png",
}

SKELETON = """<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="{description}">
<meta name="generator" content="neoantigen-selection {version}">
{head}
</head>
<body>
{body}
</body>
</html>
"""

DESCRIPTION = ("从 388 个体细胞突变到一条 mRNA 里的 34 个新抗原：完整选择流程、"
               "构建体设计、TESLA 真实标签基准，以及与另一个新抗原包的逐项对比。")


def data_uri(path: str) -> str:
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def selected_rows(demo: str) -> str:
    import pandas as pd
    sel = pd.read_csv(os.path.join(demo, "selected.csv")).head(12)
    out = []
    for r in sel.itertuples():
        out.append(
            f"<tr><td>{int(r.slot)}</td><td>{r.gene}</td><td>{r.protein_change}</td>"
            f"<td>{r.allele}</td><td class='pep'>{r.mut_peptide}</td>"
            f"<td class='num'>{r.mut_rank:.2f}</td><td class='num'>{r.wt_rank:.2f}</td>"
            f"<td class='num'>{r.tpm:,.0f}</td><td class='num'>{r.neo_score:.3f}</td></tr>")
    return "\n".join(out)


def build(demo: str, outdir: str, mode: str = "standalone") -> str:
    sys.path.insert(0, ROOT)
    from neoantigen_pipeline import __version__

    tpl = open(TEMPLATE, encoding="utf-8").read()
    tpl = tpl.replace("__ROWS__", selected_rows(demo))
    for token, name in FIGURES.items():
        path = os.path.join(demo, name)
        if not os.path.exists(path):
            raise SystemExit(f"missing figure {path} -- run the demo with figures first")
        tpl = tpl.replace(token, data_uri(path))

    left = re.findall(r"__[A-Z_]+__", tpl)
    if left:
        raise SystemExit(f"unfilled placeholders: {sorted(set(left))}")

    os.makedirs(outdir, exist_ok=True)
    if mode == "artifact":
        path = os.path.join(outdir, "neoantigen-selection-artifact.html")
        open(path, "w", encoding="utf-8").write(tpl)
        return path

    # Split the template's own <title>/<link>/<style> preamble from the body,
    # so the standalone document puts them where a browser expects them.
    m = re.search(r"\n<div class=\"wrap\">", tpl)
    if not m:
        raise SystemExit("template no longer starts its body with <div class=\"wrap\">")
    head, body = tpl[:m.start()].strip(), tpl[m.start():].strip()
    doc = SKELETON.format(head=head, body=body, description=DESCRIPTION,
                          version=__version__)
    path = os.path.join(outdir, "neoantigen-selection-report.html")
    open(path, "w", encoding="utf-8").write(doc)
    return path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", default=os.path.join(ROOT, "demo_out"))
    p.add_argument("--out", default=os.path.join(ROOT, "dist"))
    p.add_argument("--mode", choices=["standalone", "artifact"], default="standalone")
    a = p.parse_args(argv)
    path = build(a.demo, a.out, a.mode)
    print(f"{path}  ({os.path.getsize(path) / 1e6:.2f} MB, {a.mode})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
