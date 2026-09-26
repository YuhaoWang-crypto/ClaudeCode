"""Stage 6 - the panel report: a self-contained HTML page with inline SVG charts.

Three charts, each chosen for the job its data does:

* **Scaffold Spearman per target** with each target's own permutation-null 95th
  percentile drawn alongside. Two series, so a legend is present; the null is the
  bar the model has to clear, and putting it on the same axis is the only way a
  reader can check that claim without cross-referencing a table.
* **RMSE as a ratio to that target's noise floor**, diverging around 1.0. Below
  1.0 the model predicts held-out scaffolds more tightly than two independent
  assays agree - which is a statement about the CV task being easier than
  reproducing an assay, not about beating the experiment, and the caption says so.
* **Predicted activity profile heatmap**, sequential single hue, for the
  compounds with the strongest in-domain predictions. This is the panel's reason
  to exist: one library against many targets gives a profile per compound.

Every chart ships a table view, because three of the palette's light-mode slots
sit below 3:1 contrast on the light surface and the relief rule applies.
"""

from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# Reference palette (validated default instance), light / dark pairs.
PAL = {
    "surface": ("#fcfcfb", "#1a1a19"),
    "text": ("#0b0b0b", "#ffffff"),
    "text2": ("#52514e", "#c3c2b7"),
    "s1": ("#2a78d6", "#3987e5"),   # blue   - primary series / sequential hue
    "s2": ("#eb6834", "#d95926"),   # orange - second series
    "s3": ("#1baf7a", "#199e70"),   # aqua
    "grid": ("#e4e3df", "#33332f"),
}

BAR_H, BAR_GAP, LABEL_W = 20, 10, 92


def css_vars() -> str:
    """Token block, declared on :root as well as .viz.

    Declaring them only on .viz left body's own `color: var(--text)` resolving
    against an undefined token: it fell back to initial black, which is
    invisible on the dark canvas that `color-scheme: dark` paints. Light mode hid
    the bug because black-on-white happened to be right. Both scopes get the
    tokens, and color-scheme is set on :root so the canvas and the text agree.
    """
    light = "\n".join(f"    --{k}: {v[0]};" for k, v in PAL.items())
    dark = "\n".join(f"    --{k}: {v[1]};" for k, v in PAL.items())
    return f""":root, .viz {{
    color-scheme: light;
{light}
  }}
  @media (prefers-color-scheme: dark) {{
    :root:where(:not([data-theme="light"])),
    :root:where(:not([data-theme="light"])) .viz {{
      color-scheme: dark;
{dark}
    }}
  }}
  :root[data-theme="dark"],
  :root[data-theme="dark"] .viz {{
    color-scheme: dark;
{dark}
  }}"""


BAR_H, BAR_GAP, LABEL_W = 20, 10, 92


def bar_chart(rows: list[dict]) -> str:
    """Scaffold Spearman per target with its permutation-null p95 alongside."""
    n = len(rows)
    h = n * (BAR_H + BAR_GAP) + 46
    w, plot_x, plot_w = 720, LABEL_W, 720 - LABEL_W - 56
    vmax = max(0.9, max(r["rho"] for r in rows) * 1.08)
    sx = lambda v: plot_x + plot_w * max(v, 0) / vmax

    parts = [
        f'<svg viewBox="0 0 {w} {h}" role="img" width="100%" '
        f'aria-label="Scaffold-split Spearman per target against its permutation null">'
    ]
    for gv in (0.2, 0.4, 0.6, 0.8):
        if gv <= vmax:
            x = sx(gv)
            parts.append(
                f'<line x1="{x:.1f}" y1="14" x2="{x:.1f}" y2="{h - 32}" '
                f'stroke="var(--grid)" stroke-width="1"/>'
                f'<text x="{x:.1f}" y="{h - 18}" text-anchor="middle" font-size="10" '
                f'fill="var(--text2)">{gv:.1f}</text>'
            )
    for i, r in enumerate(rows):
        y = 20 + i * (BAR_H + BAR_GAP)
        bw = sx(r["rho"]) - plot_x
        nx = sx(r["null"])
        parts.append(
            f'<text x="{LABEL_W - 8}" y="{y + BAR_H * 0.72:.1f}" text-anchor="end" '
            f'font-size="11" fill="var(--text)">{html.escape(r["target"])}</text>'
            # 4px rounded data-end, anchored to the baseline
            f'<rect x="{plot_x}" y="{y}" width="{max(bw, 0.1):.1f}" height="{BAR_H}" '
            f'rx="4" fill="var(--s1)"><title>{html.escape(r["target"])}: '
            f'scaffold Spearman {r["rho"]:+.3f} +- {r["sd"]:.3f} ({r["model"]}), '
            f'null p95 {r["null"]:+.3f}</title></rect>'
            # null marker: 2px tick, with a 2px surface ring so it reads over the fill
            f'<rect x="{nx - 3:.1f}" y="{y - 2}" width="6" height="{BAR_H + 4}" rx="1" '
            f'fill="var(--surface)"/>'
            f'<rect x="{nx - 1:.1f}" y="{y - 2}" width="2" height="{BAR_H + 4}" '
            f'fill="var(--s2)"><title>{html.escape(r["target"])} permutation-null '
            f'95th percentile: {r["null"]:+.3f}</title></rect>'
            f'<text x="{sx(r["rho"]) + 6:.1f}" y="{y + BAR_H * 0.72:.1f}" font-size="10" '
            f'fill="var(--text2)">{r["rho"]:.3f}</text>'
        )
    parts.append("</svg>")
    legend = (
        '<div class="legend">'
        '<span><i style="background:var(--s1)"></i>scaffold-split Spearman</span>'
        '<span><i style="background:var(--s2);width:3px"></i>permutation-null 95th pct</span>'
        "</div>"
    )
    return legend + '<div class="scroll">' + "".join(parts) + "</div>"


def ratio_chart(rows: list[dict]) -> str:
    """RMSE / noise floor, diverging about 1.0 with a neutral midpoint."""
    n = len(rows)
    h = n * (BAR_H + BAR_GAP) + 46
    w, mid_x, half = 720, LABEL_W + 250, 230
    vmax = max(1.6, max(abs(r["ratio"] - 1) for r in rows) * 1.25 + 1)
    sx = lambda v: mid_x + half * (v - 1) / (vmax - 1)

    parts = [
        f'<svg viewBox="0 0 {w} {h}" role="img" width="100%" '
        f'aria-label="Model RMSE as a ratio to each target assay noise floor">'
    ]
    parts.append(
        f'<line x1="{mid_x}" y1="14" x2="{mid_x}" y2="{h - 32}" stroke="var(--text2)" '
        f'stroke-width="1"/>'
        f'<text x="{mid_x}" y="{h - 18}" text-anchor="middle" font-size="10" '
        f'fill="var(--text2)">1.0 = at the noise floor</text>'
    )
    for i, r in enumerate(rows):
        y = 20 + i * (BAR_H + BAR_GAP)
        x = sx(r["ratio"])
        # Warm pole above the floor, cool pole below it; gray would be the midpoint.
        colour = "var(--s2)" if r["ratio"] > 1 else "var(--s3)"
        x0, x1 = (mid_x, x) if r["ratio"] > 1 else (x, mid_x)
        parts.append(
            f'<text x="{LABEL_W - 8}" y="{y + BAR_H * 0.72:.1f}" text-anchor="end" '
            f'font-size="11" fill="var(--text)">{html.escape(r["target"])}</text>'
            f'<rect x="{x0:.1f}" y="{y}" width="{max(x1 - x0, 0.1):.1f}" height="{BAR_H}" '
            f'rx="4" fill="{colour}"><title>{html.escape(r["target"])}: RMSE '
            f'{r["rmse"]:.3f} / noise floor {r["noise"]:.3f} = {r["ratio"]:.2f}</title></rect>'
            f'<text x="{x + (7 if r["ratio"] > 1 else -7):.1f}" y="{y + BAR_H * 0.72:.1f}" '
            f'font-size="10" fill="var(--text2)" '
            f'text-anchor="{"start" if r["ratio"] > 1 else "end"}">{r["ratio"]:.2f}</text>'
        )
    parts.append("</svg>")
    return (
        '<div class="legend">'
        '<span><i style="background:var(--s3)"></i>RMSE below the noise floor</span>'
        '<span><i style="background:var(--s2)"></i>RMSE above the noise floor</span>'
        "</div>" + '<div class="scroll">' + "".join(parts) + "</div>"
    )


def heatmap(mat: pd.DataFrame, tiers: pd.DataFrame, names: pd.Series) -> str:
    """Sequential single-hue heatmap of predicted pAffinity, compound x target."""
    cell, gap = 26, 2
    left, top = 210, 96
    w = left + mat.shape[1] * (cell + gap) + 16
    h = top + mat.shape[0] * (cell + gap) + 16
    vmin, vmax = float(np.nanmin(mat.values)), float(np.nanmax(mat.values))

    def shade(v: float) -> str:
        if not np.isfinite(v):
            return "var(--grid)"
        t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
        # One hue, light -> dark, as a sequential ramp must be.
        return f"color-mix(in oklab, var(--s1) {12 + 82 * t:.0f}%, var(--surface))"

    parts = [
        f'<svg viewBox="0 0 {w} {h}" role="img" width="100%" '
        f'aria-label="Predicted activity profile across the panel for top in-domain compounds">'
    ]
    for j, t in enumerate(mat.columns):
        x = left + j * (cell + gap) + cell / 2
        parts.append(
            f'<text x="{x}" y="{top - 8}" font-size="10" fill="var(--text2)" '
            f'transform="rotate(-55 {x} {top - 8})" text-anchor="start">{html.escape(str(t))}</text>'
        )
    for i, cid in enumerate(mat.index):
        y = top + i * (cell + gap)
        label = str(names.get(cid, cid))[:30]
        parts.append(
            f'<text x="{left - 8}" y="{y + cell * 0.7:.1f}" text-anchor="end" '
            f'font-size="10" fill="var(--text)">{html.escape(label)}</text>'
        )
        for j, t in enumerate(mat.columns):
            v = mat.iat[i, j]
            x = left + j * (cell + gap)
            tier = tiers.iat[i, j] if tiers is not None else ""
            ring = (
                f'<rect x="{x + 0.5}" y="{y + 0.5}" width="{cell - 1}" height="{cell - 1}" '
                f'rx="3" fill="none" stroke="var(--text)" stroke-width="1.5"/>'
                if tier == "high"
                else ""
            )
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="3" '
                f'fill="{shade(v)}"><title>{html.escape(label)} x {html.escape(str(t))}: '
                f'predicted pAffinity {v:.2f}'
                + (f" ({tier})" if tier else "")
                + "</title></rect>" + ring
            )
    parts.append("</svg>")
    return (
        '<div class="legend"><span>predicted pAffinity '
        f'{vmin:.1f} <i class="ramp"></i> {vmax:.1f}</span>'
        '<span><i class="ringed"></i>outlined = inside the applicability domain</span>'
        "</div>" + '<div class="scroll">' + "".join(parts) + "</div>"
    )


def table(df: pd.DataFrame, caption: str) -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape('' if pd.isna(v) else str(v))}</td>" for v in row) + "</tr>"
        for row in df.itertuples(index=False)
    )
    return (
        f'<details class="tv"><summary>{html.escape(caption)}</summary>'
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{body}</tbody></table></div></details>"
    )


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    bench = json.loads((RESULTS / "benchmark.json").read_text())
    prov = json.loads((RESULTS / "panel_provenance.json").read_text())
    by_name = {t["name"]: t for t in prov["targets"]}
    summary = pd.read_csv(RESULTS / "panel_target_summary.csv")
    profiles = pd.read_csv(RESULTS / "compound_profiles.csv", index_col=0)
    recall = pd.read_csv(RESULTS / "control_recall_summary.csv")

    perf = sorted(
        (
            {
                "target": b["name"],
                "model": b["best_model"],
                "rho": b["models"][b["best_model"]]["scaffold_spearman_mean"],
                "sd": b["models"][b["best_model"]]["scaffold_spearman_sd"],
                "null": b["permutation_null"]["p95"],
                "rmse": b["models"][b["best_model"]]["scaffold_rmse_mean"],
                "noise": b["noise_floor_log"] or float("nan"),
                "ratio": b["rmse_vs_noise_floor"] or float("nan"),
                "leak": b["leakage_factor_random_over_scaffold"],
                "n": b["n_compounds"],
            }
            for b in bench
        ),
        key=lambda r: -r["rho"],
    )

    # Heatmap: strongest in-domain, non-frequent-hitter profiles.
    prof = profiles[(profiles["n_targets_high_ad"] >= 2) & (~profiles["frequent_hitter"])]
    top = prof.nlargest(26, "best_pred_pAffinity")
    targets = [r["target"] for r in sorted(perf, key=lambda r: r["target"])]
    mat = top[[t for t in targets if t in top.columns]]
    tier_wide = pd.read_csv(RESULTS / "panel_matrix_ad_tier.csv", index_col=0)
    tiers = tier_wide.reindex(index=mat.index, columns=mat.columns)

    # --- selectivity validation (p07), if it has been run ------------------
    sel_path = RESULTS / "selectivity_validation.json"
    sel_section = ""
    if sel_path.exists():
        sv = json.loads(sel_path.read_text())
        sel_rows = pd.DataFrame(
            [
                {
                    "pair": x["pair"],
                    "subset": r["subset"],
                    "n": r["n"],
                    "rho(delta)": f"{r['spearman_delta']:+.3f}",
                    "sign agreement": f"{r['sign_agreement']:.0%}",
                    "median |err| of delta": f"{r['median_abs_error_of_delta']:.2f}",
                    "measured delta sd": f"{r['measured_delta_sd']:.2f}",
                }
                for x in sv
                for r in (x["all"], x["gap_above_noise"])
                if "spearman_delta" in r
            ]
        )
        verdicts = "".join(
            f"<li><b>{html.escape(x['pair'])}</b> — {html.escape(x.get('verdict', ''))} "
            f"({x['n_shared_compounds']} compounds measured on both; per-target "
            "Spearman "
            + ", ".join(
                f"{k} {v:+.2f}" for k, v in x["per_target_scaffold_spearman"].items()
            )
            + ")</li>"
            for x in sv
        )
        sel_section = f"""
<h2>Can the panel call selectivity, or only potency?</h2>
<p>A profile's most-used derived quantity is the gap between two targets, and a
gap is a difference of two predictions. Two good models can still produce a
worthless difference, and nothing in a single-target benchmark tests that. The
panel contains two paralog pairs where measured data on both members makes it
directly testable, with scaffold-split out-of-fold predictions on both sides.</p>
<ul>{verdicts}</ul>
<p><strong>The answer is pair-specific, so a blanket claim would have been
false.</strong> HDAC1/HDAC6 has 3,847 compounds measured on both and selectivity
is an explicit design objective across that series, so the models learn it.
PTP1B/PTPN11 has 150, and PTPN11 inhibitors are allosteric while PTP1B inhibitors
target the active site — little shared chemistry from which selectivity could be
learned. Both PTP models are individually strong (Spearman +0.72 and +0.86);
their difference is close to a coin flip.</p>
{table(sel_rows, "Table view: measured vs predicted difference")}
"""

    perf_tbl = pd.DataFrame(
        [
            {
                "target": r["target"],
                "BPS family": by_name[r["target"]]["bps_family"],
                "n train": r["n"],
                "model": r["model"],
                "scaffold rho": f"{r['rho']:+.3f}±{r['sd']:.3f}",
                "null p95": f"{r['null']:+.3f}",
                "RMSE": f"{r['rmse']:.3f}",
                "noise floor": f"{r['noise']:.2f}",
                "RMSE/noise": f"{r['ratio']:.2f}",
                "leakage": r["leak"],
            }
            for r in perf
        ]
    )

    body = f"""<div class="viz">
<h1>Virtual enzyme-activity panel</h1>
<p class="sub">In silico counterparts of {len(perf)} BPS Bioscience biochemical
screening assays. {sum(r['n'] for r in perf):,} training compounds from ChEMBL,
screened against {len(profiles):,} clinical and marketed small molecules.
Every number below is a scaffold-split, out-of-fold estimate.</p>

<div class="tiles">
  <div class="tile"><b>{len(perf)}</b><span>targets modellable</span></div>
  <div class="tile"><b>{sum(r['n'] for r in perf):,}</b><span>training compounds</span></div>
  <div class="tile"><b>{len(profiles):,}</b><span>library compounds profiled</span></div>
  <div class="tile"><b>{int(summary['n_hits'].sum()):,}</b><span>in-domain hits</span></div>
</div>

<h2>Does each virtual assay beat its own null?</h2>
<p>A target enters the screen only if its scaffold-split Spearman clears its own
label-permutation 95th percentile - not zero. The null is drawn on the same axis
so the claim is checkable here rather than in a table.</p>
{bar_chart(perf)}
{table(perf_tbl, "Table view: per-target performance")}

<h2>Is the model near the limit of the data?</h2>
<p>RMSE divided by that target's inter-assay noise floor. A ratio near 1 means
the model is as precise as two independent assays are with each other.
<strong>Below 1 is not "better than the experiment"</strong> - it means the
cross-validation task (predict a median over pooled sources) is easier than
reproducing an independent assay, because folds share source-specific
consistency. NAMPT's 0.17-log floor makes its ratio the least meaningful of the
set.</p>
{ratio_chart(perf)}
{sel_section}
<h2>Activity profile across the panel</h2>
<p>The reason to run a panel rather than a single assay: each compound gets a
profile. Outlined cells are inside that target's applicability domain; unoutlined
cells are extrapolation and should not be read as predictions. Frequent hitters -
top-decile on four or more targets - are excluded here and listed separately,
because a compound that looks active everywhere is more often exploiting a
fingerprint shortcut than genuinely polypharmacological.</p>
{heatmap(mat, tiers, profiles["pref_name"])}

<h2>Do the real assays' reference inhibitors come back?</h2>
<p>Each panel target's actual control compounds, looked up in the screening
output. Controls were never held out, so an in-training control validates the fit
rather than generalisation, and the count is shown separately.</p>
{table(recall, "Table view: reference-inhibitor recall")}

<h2>Per-target hit counts</h2>
{table(summary, "Table view: per-target screen summary")}
</div>"""

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Virtual Enzyme Panel</title>
<style>
  {css_vars()}
  * {{ box-sizing: border-box; }}
  body {{ margin:0; background:var(--surface); color:var(--text);
    font:14px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif; }}
  .viz {{ max-width:900px; margin:0 auto; padding:32px 16px 64px; background:var(--surface); }}
  h1 {{ font-size:26px; margin:0 0 6px; letter-spacing:-.01em; }}
  h2 {{ font-size:17px; margin:40px 0 6px; }}
  p {{ color:var(--text2); margin:0 0 14px; max-width:66ch; }}
  p strong {{ color:var(--text); }}
  .sub {{ font-size:14px; }}
  .tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
    gap:10px; margin:20px 0 8px; }}
  .tile {{ border:1px solid var(--grid); border-radius:10px; padding:12px 14px; }}
  .tile b {{ display:block; font-size:24px; letter-spacing:-.02em; }}
  .tile span {{ color:var(--text2); font-size:12px; }}
  .legend {{ display:flex; gap:18px; flex-wrap:wrap; color:var(--text2);
    font-size:12px; margin:10px 0 4px; align-items:center; }}
  .legend i {{ display:inline-block; width:10px; height:10px; border-radius:2px;
    margin-right:6px; vertical-align:-1px; }}
  .legend i.ramp {{ width:80px; height:10px; border-radius:2px;
    background:linear-gradient(90deg,
      color-mix(in oklab, var(--s1) 12%, var(--surface)),
      color-mix(in oklab, var(--s1) 94%, var(--surface))); }}
  .legend i.ringed {{ border:1.5px solid var(--text); background:transparent; }}
  svg {{ display:block; max-width:100%; margin:4px 0 6px; }}
  .scroll svg {{ min-width:520px; }}
  /* Tables keep nowrap cells for scanability, so they scroll inside their own
     box rather than pushing the page wide - at 390px the page must not scroll
     horizontally. */
  .tv {{ margin:8px 0 4px; }}
  .tv[open] {{ overflow-x:auto; }}
  .scroll {{ overflow-x:auto; -webkit-overflow-scrolling:touch; }}
  .tv summary {{ cursor:pointer; color:var(--text2); font-size:12px; }}
  table {{ border-collapse:collapse; width:100%; margin:10px 0; font-size:12px; }}
  th,td {{ border-bottom:1px solid var(--grid); padding:5px 8px; text-align:left;
    white-space:nowrap; }}
  th {{ color:var(--text2); font-weight:600; }}
  @media (max-width:640px) {{ .viz {{ padding:20px 16px 48px; }} h1 {{ font-size:21px; }} }}
</style></head><body>{body}</body></html>"""

    out = FIGURES / "panel_report.html"
    out.write_text(doc, encoding="utf-8")
    print(f"wrote {out} ({len(doc) / 1024:.0f} KB)")
    print(f"  {len(perf)} targets charted, heatmap {mat.shape}")


if __name__ == "__main__":
    main()
