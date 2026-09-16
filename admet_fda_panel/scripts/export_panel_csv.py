"""Export tables 8 and 9 as CSV (plain text, markup stripped)."""

from __future__ import annotations
import re
import csv
import pandas as pd

import safety_panel as SP
import mapping_core as MC

OUT = "/home/user/results/admet_fda_panel"


def strip(s):
    s = re.sub(r"<br\s*/?>", " ", str(s))
    s = re.sub(r"<[^>]+>", "", s)
    return s.strip()


def main():
    # table 8
    rows = [[strip(c) for c in r] for r in SP.COMPARISON_ROWS]
    pd.DataFrame(rows, columns=[strip(h) for h in SP.COMPARISON_HDR]).to_csv(
        f"{OUT}/safety_panel_comparison.csv", index=False)

    # table 9 + table 10
    stats = MC.endpoint_stats()
    val = MC.validation_lookup()
    out = []
    for fam, targets, ep, in44, in87, note, kind in SP.COVERAGE_SPEC:
        v = val.get(ep) if ep else None
        out.append({
            "target_family": fam,
            "panel_targets": strip(targets),
            "admet_ai_endpoint": ep or "",
            "in_bowes44": in44,
            "in_safetyscreen87": in87,
            "panel_positives_in_this_study": strip(stats[ep]["positives"]) if ep else "",
            "n_positive_of_30": stats[ep]["n_positive"] if ep else "",
            "positive_control_recovered": f"{v[0]}/{v[1]}" if v else "",
            "verdict": {"covered": "mapped and validated",
                        "nominal": "nominal mapping only - validation failed",
                        "nominal87": "nominal mapping only (87 panel) - validation failed",
                        "none": "requires in vitro assay"}[kind],
            "note": strip(note),
        })
    for label, ep, note in SP.NON_PANEL_SPEC:
        out.append({
            "target_family": "不属于面板 (DMPK/Tox21)",
            "panel_targets": strip(label),
            "admet_ai_endpoint": ep or "(none)",
            "in_bowes44": "否", "in_safetyscreen87": "否",
            "panel_positives_in_this_study": strip(stats[ep]["positives"]) if ep else "",
            "n_positive_of_30": stats[ep]["n_positive"] if ep else "",
            "positive_control_recovered":
                (lambda v: f"{v[0]}/{v[1]}" if v else "")(val.get(ep)) if ep else "",
            "verdict": "not a secondary-pharmacology panel target",
            "note": strip(note),
        })
    pd.DataFrame(out).to_csv(f"{OUT}/safety_panel_coverage.csv", index=False)
    print(f"✓ wrote safety_panel_comparison.csv ({len(rows)} rows)")
    print(f"✓ wrote safety_panel_coverage.csv ({len(out)} rows)")


if __name__ == "__main__":
    main()
