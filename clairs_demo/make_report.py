#!/usr/bin/env python3
"""Inline figures/*.png into report_template.html as data URIs.

The published report has to be a single self-contained file (no external
requests), so each figure is base64-embedded rather than linked.

    python3 clairs_demo/make_report.py [OUT_HTML]
"""
import base64
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = {
    "__FIG_DETECTION__": "detection_by_vaf.png",
    "__FIG_AF__": "af_vs_seqc2.png",
    "__FIG_DIST__": "distributions.png",
    "__FIG_SIG__": "spectrum_and_loh.png",
}

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "report.html")
html = open(os.path.join(HERE, "report_template.html")).read()
for token, name in FIGURES.items():
    with open(os.path.join(HERE, "figures", name), "rb") as fh:
        html = html.replace(token, "data:image/png;base64," +
                            base64.b64encode(fh.read()).decode())
if "__FIG" in html:
    raise SystemExit("unreplaced figure placeholder in report_template.html")
with open(out, "w") as fh:
    fh.write(html)
print(f"wrote {out} ({os.path.getsize(out) / 1024:.0f} KB)")
