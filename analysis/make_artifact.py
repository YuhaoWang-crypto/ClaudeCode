#!/usr/bin/env python3
"""Generate a self-contained HTML artifact summarizing the GSE120575 reanalysis
and external validation. Figures are embedded as base64 data URIs so the page
is fully self-contained (Artifact CSP blocks external hosts)."""
import base64
import os
import pandas as pd

OUT = "/home/user/ClaudeCode/results"
FIG = os.path.join(OUT, "figures")


def img(name):
    with open(os.path.join(FIG, name), "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return f"data:image/png;base64,{b64}"


def fig(name, caption):
    return (f'<figure class="fig"><img alt="{caption}" src="{img(name)}"/>'
            f'<figcaption>{caption}</figcaption></figure>')


def chip(label):
    cls = "r" if label == "R" else "nr"
    return f'<span class="chip {cls}">{label}</span>'


# ---- data tables -----------------------------------------------------------
def pop_table():
    df = pd.read_csv(os.path.join(OUT, "population_enrichment_R_vs_NR.csv"))
    rows = ""
    for _, r in df.iterrows():
        direction = "r" if r.log2FC_R_vs_NR > 0 else "nr"
        sig = "sig" if r.pval_mannwhitney < 0.05 else ""
        bar = min(abs(r.log2FC_R_vs_NR) / 2.0, 1.0) * 100
        rows += f"""<tr>
          <td class="mono">{r.cell_type}</td>
          <td class="num">{r['mean_frac_Responder']:.3f}</td>
          <td class="num">{r['mean_frac_Non-responder']:.3f}</td>
          <td class="lfc {direction}"><span class="lfcbar {direction}" style="width:{bar:.0f}%"></span><span class="lfcval">{r.log2FC_R_vs_NR:+.2f}</span></td>
          <td class="num {sig}">{r.pval_mannwhitney:.03g}</td></tr>"""
    return rows


def riaz_table():
    df = pd.read_csv(os.path.join(OUT, "riaz_validation_stats.csv"))
    rows = ""
    for _, r in df.iterrows():
        sig = "sig" if r.pval < 0.05 else ""
        rows += (f'<tr><td class="mono">{r.visit}</td><td class="mono">{r.score}</td>'
                 f'<td class="num">{r.AUC:.2f}</td><td class="num {sig}">{r.pval:.03g}</td></tr>')
    return rows


def cd8_table():
    df = pd.read_csv(os.path.join(OUT, "cd8_substate_summary.csv"))
    df = df.sort_values("log2FC_R_vs_NR", ascending=False)
    rows = ""
    for _, r in df.iterrows():
        d = "r" if r.log2FC_R_vs_NR > 0 else "nr"
        rows += (f'<tr><td class="mono">{r.cd8_state}</td><td class="num">{int(r.n)}</td>'
                 f'<td class="num">{r.CD8_G:.2f}</td><td class="num">{r.CD8_B:.2f}</td>'
                 f'<td class="lfc {d}"><span class="lfcval">{r.log2FC_R_vs_NR:+.2f}</span></td></tr>')
    return rows


HTML = f"""<title>Immunotherapy response signatures — single-cell reanalysis</title>
<style>
:root {{
  --bg:#f7f9fa; --surface:#ffffff; --surface2:#f0f4f6; --ink:#16232b;
  --muted:#5c6f7a; --line:#e2e8ec; --accent:#1f7a9c; --accent-ink:#12586f;
  --danger:#c85a3f; --navy:#0e2b39; --shadow:0 1px 2px rgba(14,43,57,.06),0 8px 24px rgba(14,43,57,.05);
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg:#0c1317; --surface:#131e24; --surface2:#182530; --ink:#e7eef2;
    --muted:#93a4ae; --line:#22323a; --accent:#4fb3d6; --accent-ink:#8fd3ea;
    --danger:#e07b5f; --navy:#081a24; --shadow:0 1px 2px rgba(0,0,0,.4),0 10px 30px rgba(0,0,0,.35);
  }}
}}
:root[data-theme="light"] {{
  --bg:#f7f9fa; --surface:#ffffff; --surface2:#f0f4f6; --ink:#16232b; --muted:#5c6f7a;
  --line:#e2e8ec; --accent:#1f7a9c; --accent-ink:#12586f; --danger:#c85a3f; --navy:#0e2b39;
  --shadow:0 1px 2px rgba(14,43,57,.06),0 8px 24px rgba(14,43,57,.05);
}}
:root[data-theme="dark"] {{
  --bg:#0c1317; --surface:#131e24; --surface2:#182530; --ink:#e7eef2; --muted:#93a4ae;
  --line:#22323a; --accent:#4fb3d6; --accent-ink:#8fd3ea; --danger:#e07b5f; --navy:#081a24;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 10px 30px rgba(0,0,0,.35);
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--bg); color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  line-height:1.6; -webkit-font-smoothing:antialiased;
}}
.serif {{ font-family:ui-serif,Georgia,"Iowan Old Style","Times New Roman",serif; }}
.mono {{ font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace; }}
.num {{ font-variant-numeric:tabular-nums; text-align:right; }}

/* masthead */
header.top {{
  position:sticky; top:0; z-index:10; display:flex; justify-content:space-between;
  align-items:center; gap:1rem; padding:.7rem clamp(1rem,4vw,2.5rem);
  background:color-mix(in srgb,var(--surface) 88%,transparent); backdrop-filter:blur(8px);
  border-bottom:1px solid var(--line);
}}
header.top .brand {{ font-family:ui-monospace,monospace; font-size:.78rem; letter-spacing:.04em;
  color:var(--muted); text-transform:uppercase; }}
header.top .brand b {{ color:var(--ink); }}
header.top .repo {{ font-family:ui-monospace,monospace; font-size:.74rem; color:var(--accent-ink); }}

.wrap {{ max-width:1080px; margin:0 auto; padding:0 clamp(1rem,4vw,2.5rem); }}
.read {{ max-width:66ch; }}

/* hero */
.hero {{ padding:clamp(2.5rem,7vw,5rem) 0 2rem; border-bottom:1px solid var(--line); }}
.eyebrow {{ font-family:ui-monospace,monospace; font-size:.75rem; letter-spacing:.18em;
  text-transform:uppercase; color:var(--accent-ink); margin-bottom:1rem; }}
h1 {{ font-size:clamp(2rem,5.4vw,3.4rem); line-height:1.06; margin:.2rem 0 1rem;
  font-weight:600; letter-spacing:-.01em; text-wrap:balance; }}
h1 em {{ font-style:italic; color:var(--accent); }}
.lede {{ font-size:1.14rem; color:var(--muted); max-width:60ch; }}
.lede b {{ color:var(--ink); font-weight:600; }}

/* stat tiles */
.stats {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:1px;
  background:var(--line); border:1px solid var(--line); border-radius:12px; overflow:hidden;
  margin-top:2.2rem; }}
.stat {{ background:var(--surface); padding:1.1rem 1.2rem; }}
.stat .v {{ font-family:ui-serif,Georgia,serif; font-size:1.9rem; font-weight:600;
  font-variant-numeric:tabular-nums; line-height:1; color:var(--accent); }}
.stat .v.warn {{ color:var(--danger); }}
.stat .k {{ font-size:.76rem; color:var(--muted); margin-top:.45rem; line-height:1.35; }}

/* sections */
section {{ padding:clamp(2.2rem,5vw,3.6rem) 0; border-bottom:1px solid var(--line); }}
.kicker {{ font-family:ui-monospace,monospace; font-size:.74rem; letter-spacing:.14em;
  text-transform:uppercase; color:var(--danger); display:flex; align-items:center; gap:.6rem; }}
.kicker::before {{ content:""; width:1.6rem; height:2px; background:var(--danger); display:inline-block; }}
h2 {{ font-family:ui-serif,Georgia,serif; font-size:clamp(1.5rem,3.4vw,2.1rem); font-weight:600;
  margin:.6rem 0 1rem; letter-spacing:-.01em; text-wrap:balance; }}
p {{ margin:0 0 1rem; }}
section p {{ color:var(--ink); }}
.muted {{ color:var(--muted); }}
gene, .gene {{ font-family:ui-monospace,monospace; font-size:.92em; background:var(--surface2);
  padding:.05em .35em; border-radius:4px; color:var(--accent-ink); }}

/* figure */
.fig {{ margin:1.6rem 0; background:var(--surface); border:1px solid var(--line);
  border-radius:12px; padding:1rem; box-shadow:var(--shadow); }}
.fig img {{ width:100%; height:auto; display:block; border-radius:6px; }}
.fig figcaption {{ font-size:.82rem; color:var(--muted); margin-top:.7rem; padding:0 .2rem; }}
.grid2 {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:1.2rem; }}
.grid2 .fig {{ margin:0; }}

/* tables */
.tablewrap {{ overflow-x:auto; margin:1.4rem 0; border:1px solid var(--line); border-radius:12px; }}
table {{ border-collapse:collapse; width:100%; font-size:.9rem; background:var(--surface); }}
th, td {{ padding:.6rem .9rem; text-align:left; border-bottom:1px solid var(--line); white-space:nowrap; }}
thead th {{ font-family:ui-monospace,monospace; font-size:.72rem; letter-spacing:.05em;
  text-transform:uppercase; color:var(--muted); font-weight:600; background:var(--surface2); }}
tbody tr:last-child td {{ border-bottom:none; }}
td.sig {{ color:var(--accent-ink); font-weight:600; }}
.lfc {{ position:relative; min-width:120px; }}
.lfcbar {{ position:absolute; left:0; top:50%; transform:translateY(-50%); height:60%;
  border-radius:3px; opacity:.22; }}
.lfcbar.r {{ background:var(--accent); }}
.lfcbar.nr {{ background:var(--danger); }}
.lfcval {{ position:relative; font-variant-numeric:tabular-nums; font-family:ui-monospace,monospace; }}
.lfc.r .lfcval {{ color:var(--accent-ink); }}
.lfc.nr .lfcval {{ color:var(--danger); }}
.chip {{ font-family:ui-monospace,monospace; font-size:.72rem; padding:.1em .5em; border-radius:20px;
  font-weight:600; }}
.chip.r {{ background:color-mix(in srgb,var(--accent) 18%,transparent); color:var(--accent-ink); }}
.chip.nr {{ background:color-mix(in srgb,var(--danger) 18%,transparent); color:var(--danger); }}

/* callout / formula */
.callout {{ background:var(--navy); color:#dfeaef; border-radius:14px; padding:1.4rem 1.6rem;
  margin:1.6rem 0; box-shadow:var(--shadow); }}
.callout .lab {{ font-family:ui-monospace,monospace; font-size:.72rem; letter-spacing:.14em;
  text-transform:uppercase; color:#7fc4dd; margin-bottom:.7rem; }}
.callout code {{ font-family:ui-monospace,monospace; font-size:.9rem; line-height:1.9; display:block;
  white-space:pre-wrap; word-break:break-word; }}
.callout .plus {{ color:#7fc4dd; }} .callout .minus {{ color:#e79b86; }}

.two {{ display:grid; grid-template-columns:1.1fr .9fr; gap:2rem; align-items:start; }}
@media (max-width:760px) {{ .two {{ grid-template-columns:1fr; }} }}

.note {{ border-left:3px solid var(--accent); background:var(--surface2); padding:.9rem 1.1rem;
  border-radius:0 10px 10px 0; font-size:.92rem; color:var(--muted); margin:1.4rem 0; }}
.note b {{ color:var(--ink); }}

footer {{ padding:3rem 0 4rem; color:var(--muted); font-size:.86rem; }}
footer .repro {{ background:var(--surface); border:1px solid var(--line); border-radius:12px;
  padding:1.1rem 1.3rem; font-family:ui-monospace,monospace; font-size:.8rem; overflow-x:auto;
  color:var(--ink); white-space:pre; margin:1rem 0; }}
a {{ color:var(--accent-ink); }}
:focus-visible {{ outline:2px solid var(--accent); outline-offset:2px; border-radius:3px; }}
</style>

<header class="top">
  <div class="brand"><b>GSE120575</b> · scRNA-seq reanalysis</div>
  <div class="repo mono">claude/scrna-immunotherapy-signatures</div>
</header>

<div class="wrap">
  <div class="hero">
    <div class="eyebrow">Single-cell immuno-oncology · melanoma · anti–PD-1 / anti–CTLA-4</div>
    <h1 class="serif">It isn't how many CD8 T cells a tumor has — it's <em>what state</em> they're in.</h1>
    <p class="lede">A reproducible reanalysis of <b>16,291 immune cells</b> from 48 melanoma
    biopsies taken before and after checkpoint immunotherapy (Sade-Feldman <i>et al.</i>,
    <i>Cell</i> 2018), extended with external validation on the <b>Riaz</b> anti–PD-1 cohort and
    <b>TCR clonal dynamics</b> from the <b>Yost</b> basal/squamous carcinoma cohort.</p>
    <div class="stats">
      <div class="stat"><div class="v">16,291</div><div class="k">CD45⁺ single cells across 48 biopsies</div></div>
      <div class="stat"><div class="v">0.83</div><div class="k">AUC — CD8 memory−exhaustion signature separates responders</div></div>
      <div class="stat"><div class="v">+1.94</div><div class="k">log₂FC B-cell enrichment in responders (p=0.003)</div></div>
      <div class="stat"><div class="v warn">54%</div><div class="k">of post-treatment expanded T-cell clones are newly recruited</div></div>
    </div>
  </div>

  <section>
    <div class="kicker">The immune atlas</div>
    <h2 class="serif">Eight immune populations, resolved from raw counts</h2>
    <div class="read">
      <p>Starting from the raw genes × cells TPM matrix (GEO series <b>GSE120575</b>), the pipeline
      runs QC → log-normalization → 2,000 highly-variable genes → PCA → Leiden clustering (24 clusters)
      → UMAP, then annotates clusters to lineage by canonical markers. The immune compartment resolves
      cleanly into CD8 T, CD4 T, Treg, B, plasma, NK, myeloid and pDC.</p>
    </div>
    {fig("umap_celltypes.png", "UMAP of 16,291 CD45⁺ cells, colored by annotated lineage.")}
  </section>

  <section>
    <div class="kicker">Who expands, who contracts</div>
    <h2 class="serif">B and CD4 T cells mark responders; myeloid and pDC mark non-responders</h2>
    <div class="read">
      <p>Per-lesion population fractions, contrasted between responding and non-responding tumors.
      Positive log₂ fold-change = expanded in responders. Notably, the <i>fraction</i> of CD8 T cells
      does <b>not</b> separate response (p=0.15) — their <i>state</i> does (next section). Aggregated
      across patients, no population shifts significantly with treatment timepoint alone, reproducing
      the original study.</p>
    </div>
    <div class="tablewrap"><table>
      <thead><tr><th>population</th><th class="num">frac {chip("R")}</th><th class="num">frac {chip("NR")}</th><th>log₂FC (R/NR)</th><th class="num">p</th></tr></thead>
      <tbody>{pop_table()}</tbody>
    </table></div>
  </section>

  <section>
    <div class="kicker">The signature</div>
    <h2 class="serif">A memory-minus-exhaustion score stratifies responders</h2>
    <div class="two">
      <div>
        <p>Within the 6,699 CD8 T cells, two programs anchor the response axis: a
        <span class="gene">memory-like</span> set (CD8_G) and an <span class="gene">exhaustion</span>
        set (CD8_B). The per-lesion difference of their mean scores separates responders with
        <b>AUC 0.83</b> (Mann–Whitney p = 2.3×10⁻⁴, 47 lesions). Responders sit at +0.10, non-responders
        at −0.30. The single marker <span class="gene">TCF7</span> alone reaches AUC 0.73.</p>
        <div class="callout">
          <div class="lab">Proposed responder signature</div>
          <code>score = <span class="plus">mean(</span>TCF7, IL7R, SELL, CCR7,
       CD28, LEF1, CD27, GZMK<span class="plus">)</span>
      <span class="minus">− mean(</span>PDCD1, HAVCR2, LAG3, TIGIT,
        CTLA4, ENTPD1, CD38, TOX<span class="minus">)</span></code>
        </div>
      </div>
      {fig("signature_boxplot.png", "Per-lesion CD8 memory−exhaustion score, responders vs non-responders.")}
    </div>
  </section>

  <section>
    <div class="kicker">External validation · Riaz anti–PD-1</div>
    <h2 class="serif">The state signature is a single-cell property — bulk sees infiltration</h2>
    <div class="read">
      <p>Applied to bulk RNA-seq of the independent Riaz nivolumab cohort (GEO <b>GSE91061</b>,
      109 samples), the <i>difference</i> signature attenuates (AUC ≈ 0.4–0.54). The reason is
      biological: in bulk tissue, exhaustion genes track overall <b>T-cell infiltration</b>, which is
      favorable — the opposite of their meaning <i>within</i> a CD8 cell. What does validate in bulk is
      CD8/immune <b>abundance</b>, strongest on-treatment (<b>AUC 0.78, p=0.006</b>).</p>
    </div>
    <div class="two">
      <div class="tablewrap"><table>
        <thead><tr><th>visit</th><th>score</th><th class="num">AUC</th><th class="num">p</th></tr></thead>
        <tbody>{riaz_table()}</tbody>
      </table></div>
      {fig("riaz_validation.png", "Riaz cohort: AUC of each signature component (CR/PR vs PD), pre and on-treatment.")}
    </div>
    <div class="note"><b>Takeaway:</b> deploying a single-cell state signature on bulk requires
    isolating CD8 state from CD8 quantity — so the right cross-check is other single-cell data.</div>
  </section>

  <section>
    <div class="kicker">CD8 sub-states</div>
    <h2 class="serif">Five CD8 states — proliferating &amp; exhausted flag non-responders</h2>
    <div class="read">
      <p>Sub-clustering the CD8 compartment yields five states. Memory scores (CD8_G) peak in the
      naive/memory state; exhaustion scores (CD8_B) peak in the exhausted and proliferating states.
      Enrichment mirrors outcome: naive/memory and cytotoxic states lean toward responders, while the
      <b>proliferating</b> (log₂FC −1.43) and <b>exhausted</b> (−0.63) states lean strongly toward
      non-responders — matching the original study's cell-cycle⁺exhaustion cluster.</p>
    </div>
    <div class="tablewrap"><table>
      <thead><tr><th>CD8 state</th><th class="num">n</th><th class="num">CD8_G</th><th class="num">CD8_B</th><th>log₂FC (R/NR)</th></tr></thead>
      <tbody>{cd8_table()}</tbody>
    </table></div>
    <div class="grid2">
      {fig("umap_cd8_states.png", "CD8 T cells re-embedded and colored by transcriptional state.")}
      {fig("cd8_state_enrichment.png", "CD8 states enriched in responders (blue) vs non-responders (red).")}
    </div>
  </section>

  <section>
    <div class="kicker">Cross-cancer &amp; TCR clonal dynamics · Yost anti–PD-1</div>
    <h2 class="serif">The gene sets generalize — and exhausted clones expand by replacement</h2>
    <div class="read">
      <p>In an entirely different tumor type (basal/squamous cell carcinoma, GEO <b>GSE123813</b>),
      the same gene sets cleanly separate that study's CD8 states: memory score peaks in
      <span class="gene">CD8_mem</span>, exhaustion score peaks in <span class="gene">CD8_ex</span> and
      <span class="gene">CD8_ex_act</span>. Paired TCR sequencing shows exhausted CD8 cells are the most
      clonally expanded, the activated-exhausted state jumps from 40% → 80% expanded after therapy, and
      on average <b>54% of post-treatment expanded clones are novel</b> — direct evidence of clonal
      replacement rather than reinvigoration of pre-existing clones.</p>
    </div>
    <div class="grid2">
      {fig("yost_geneset_validation.png", "Signature gene sets scored on Yost CD8 states — memory vs exhaustion separate as expected.")}
      {fig("yost_tcr_dynamics.png", "TCR clonal expansion by CD8 state (left) and pre vs post anti–PD-1 (right).")}
    </div>
  </section>

  <footer>
    <div class="read">
      <p class="muted"><b style="color:var(--ink)">Reproduce.</b> All results regenerate from raw
      GEO data. Large matrices are referenced by name, not committed.</p>
    </div>
    <div class="repro">python3 -m venv .venv && .venv/bin/pip install -r analysis/requirements.txt
bash analysis/download_data.sh          # GSE120575 (+ GSE91061, GSE123813)
.venv/bin/python analysis/run_scrna_analysis.py
.venv/bin/python analysis/validate_riaz.py
.venv/bin/python analysis/cd8_states_tcr.py</div>
    <p class="muted">Datasets: <b>GSE120575</b> (Sade-Feldman, Cell 2018) · <b>GSE91061</b>
    (Riaz, Cell 2017) · <b>GSE123813</b> (Yost, Nat Med 2019). Reanalysis for research illustration;
    not a clinical tool.</p>
  </footer>
</div>
"""

with open(os.path.join(OUT, "artifact.html"), "w") as f:
    f.write(HTML)
print("wrote", os.path.join(OUT, "artifact.html"), f"({len(HTML)/1024:.0f} KB)")
