#!/usr/bin/env python3
"""从 vcff/data 重新生成 demo/bundle.json 与 demo/vcff_demo.html。

前端引擎 (demo/app.js) 与 Python 引擎 (vcff/) 共用同一套组合律, 两者数值一致;
改了任一侧的组合律都要同步另一侧, 并跑 python -m vcff.tests 与 demo 的场景比对。
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "vcff" / "data"


def load(name):
    return json.loads((DATA / f"{name}.json").read_text(encoding="utf-8"))


def build_bundle():
    tox, prt, phe, cmb, twn = (load(n) for n in ("tox", "prt", "phe", "cmb", "twn"))
    return {
        # norm2idx 在前端从 genes 重建; scatter 点在前端用不到
        "tox": {"genes": tox["genes"],
                "hepg2": [round(v, 4) for v in tox["hepg2"]],
                "rpe1": [round(v, 4) for v in tox["rpe1"]],
                "drug2target": tox["drug2target"]},
        "prt": {"genes": prt["genes"], "keys": prt["keys"], "lines": prt["lines"],
                "combos": {k: {"live_r": v["live_r"], "ceiling": v["ceiling"],
                               "pw_names": v["pw_names"], "pw_vals": v["pw_vals"],
                               "up": v["up"][:6], "dn": v["dn"][:6]}
                           for k, v in prt["combos"].items()}},
        "phe": {"drugs": phe["drugs"],
                "per": {k: {"loo_r": v["loo_r"], "known_moa": v["known_moa"],
                            "pred_moa": v["pred_moa"], "pw_names": v["pw_names"],
                            "pw_vals": v["pw_vals"], "neighbors": v["neighbors"],
                            "up": v["up"][:6], "dn": v["dn"][:6]}
                        for k, v in phe["per"].items()}},
        "cmb": {k: cmb[k] for k in ("pairs", "non", "epi", "rank", "ncells")},
        "twn": twn,
    }


def main():
    bundle = json.dumps(build_bundle(), ensure_ascii=False, separators=(",", ":"))
    (ROOT / "bundle.json").write_text(bundle, encoding="utf-8")
    html = (ROOT / "shell.html").read_text(encoding="utf-8")
    html += (ROOT / "page.html").read_text(encoding="utf-8")
    html += '\n<script type="application/json" id="vcff-data">' + bundle + "</script>\n"
    html += "<script>\n" + (ROOT / "app.js").read_text(encoding="utf-8") + "\n</script>\n"
    (ROOT / "vcff_demo.html").write_text(html, encoding="utf-8")
    print(f"bundle {len(bundle):,} B  →  vcff_demo.html {len(html):,} B")


if __name__ == "__main__":
    main()
