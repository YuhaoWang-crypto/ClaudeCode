#!/usr/bin/env python3
"""
Assemble the HTML report from results/*.tsv and figures/*.png.

Everything numeric on the page is read from the result files - no number is
typed into the template - so a re-run of the pipeline regenerates a report that
matches its own data.
"""
import base64
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, read_metadata, data_path,  # noqa: E402
                    results_path, figures_path, ROOT)


def tsv(name):
    p = results_path(name)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        head = f.readline()
        if head.startswith("(none)"):
            return []
        f.seek(0)
        return list(csv.DictReader(f, delimiter="\t"))


def jsn(name):
    p = results_path(name)
    return json.load(open(p)) if os.path.exists(p) else {}


def img(name):
    p = figures_path(name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" alt="{name}">'


BAND_TAG = {"comparable-to-benchmark": "good", "modestly-elevated": "warn",
            "elevated": "bad", "high": "bad"}
BAND_SHORT = {"comparable-to-benchmark": "comparable", "modestly-elevated": "modest",
              "elevated": "elevated", "high": "high", "n/a (control)": "control"}


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def table(rows, cols, labels=None, cls="", fmt=None):
    if not rows:
        return "<p class='muted'>no rows</p>"
    labels = labels or cols
    fmt = fmt or {}
    # a label prefixed with "!" is raw HTML and is not escaped
    h = "".join(f"<th>{l[1:] if l.startswith('!') else esc(l)}</th>" for l in labels)
    body = []
    for r in rows:
        tds = []
        for c in cols:
            v = r.get(c, "")
            v = fmt[c](r) if c in fmt else esc(v)
            tds.append(f"<td>{v}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return (f"<div class='tw'><table class='{cls}'><thead><tr>{h}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


CSS = """
:root{
  --bg:#ffffff; --fg:#16202b; --muted:#63748a; --line:#e2e8ef; --soft:#f5f8fa;
  --accent:#2f6f9f; --bad:#b23b3b; --warn:#c9752b; --good:#3d8a6b; --violet:#7a5ea8;
  --chip:#eef3f8;
}
:root:not([data-theme="light"]){}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#11161c; --fg:#e6edf5; --muted:#94a6ba; --line:#26313d; --soft:#171e26;
    --accent:#6fa8d4; --bad:#e07b7b; --warn:#e0a165; --good:#6cc0a0; --violet:#b39ddb;
    --chip:#1d2732;
  }
}
:root[data-theme="dark"]{
  --bg:#11161c; --fg:#e6edf5; --muted:#94a6ba; --line:#26313d; --soft:#171e26;
  --accent:#6fa8d4; --bad:#e07b7b; --warn:#e0a165; --good:#6cc0a0; --violet:#b39ddb;
  --chip:#1d2732;
}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--fg);margin:0;
  font:16px/1.62 ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  -webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px 96px}
header.hero{border-bottom:1px solid var(--line);padding:56px 0 34px;margin-bottom:8px}
.kicker{font-size:12px;letter-spacing:.13em;text-transform:uppercase;color:var(--muted);
  font-weight:650}
h1{font-size:38px;line-height:1.15;margin:12px 0 10px;letter-spacing:-.02em}
h2{font-size:25px;margin:56px 0 6px;letter-spacing:-.015em;scroll-margin-top:20px}
h3{font-size:17px;margin:30px 0 6px}
.lede{color:var(--muted);font-size:17px;max-width:74ch;margin:0}
p{max-width:78ch}
.muted{color:var(--muted)}
small{color:var(--muted)}
code{background:var(--soft);border:1px solid var(--line);border-radius:4px;
  padding:1px 5px;font-size:.87em;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.9em}
.rule{height:1px;background:var(--line);border:0;margin:44px 0}
.grid{display:grid;gap:14px;margin:22px 0}
.g4{grid-template-columns:repeat(auto-fit,minmax(168px,1fr))}
.g3{grid-template-columns:repeat(auto-fit,minmax(230px,1fr))}
.g2{grid-template-columns:repeat(auto-fit,minmax(320px,1fr))}
.stat{border:1px solid var(--line);border-radius:10px;padding:16px 16px 14px;background:var(--soft)}
.stat .n{font-size:29px;font-weight:680;letter-spacing:-.02em;line-height:1.1}
.stat .l{font-size:12px;color:var(--muted);margin-top:5px;text-transform:uppercase;
  letter-spacing:.07em;font-weight:620}
.stat .s{font-size:13px;color:var(--muted);margin-top:7px}
.card{border:1px solid var(--line);border-radius:10px;padding:18px 20px;background:var(--bg)}
.card h3{margin-top:0}
.tw{overflow-x:auto;margin:18px 0;border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th{background:var(--soft);text-align:left;padding:9px 12px;font-weight:640;
  border-bottom:1px solid var(--line);white-space:nowrap;font-size:12px;
  letter-spacing:.03em;text-transform:uppercase;color:var(--muted)}
th .unit{text-transform:none;letter-spacing:0}
td{padding:8px 12px;border-bottom:1px solid var(--line);vertical-align:top}
tbody tr:last-child td{border-bottom:0}
.num td:nth-child(n+3){text-align:right;font-variant-numeric:tabular-nums}
.num td:last-child{text-align:left}
td .wrapmono{display:block;max-width:34ch;white-space:normal;word-break:break-word;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;line-height:1.45}
img{max-width:100%;height:auto;display:block;border-radius:10px;border:1px solid var(--line);
  margin:20px 0}
.tag{display:inline-block;font-size:11px;font-weight:660;padding:2px 8px;border-radius:999px;
  background:var(--chip);color:var(--muted);letter-spacing:.03em;white-space:nowrap}
.tag.bad{background:rgba(178,59,59,.14);color:var(--bad)}
.tag.warn{background:rgba(201,117,43,.16);color:var(--warn)}
.tag.good{background:rgba(61,138,107,.16);color:var(--good)}
.tag.acc{background:rgba(47,111,159,.14);color:var(--accent)}
.tag.vio{background:rgba(122,94,168,.16);color:var(--violet)}
.callout{border:1px solid var(--line);border-left:3px solid var(--accent);
  border-radius:0 10px 10px 0;padding:14px 18px;background:var(--soft);margin:22px 0}
.callout.bad{border-left-color:var(--bad)}
.callout.good{border-left-color:var(--good)}
.callout p{margin:0;max-width:none}
ul,ol{max-width:78ch}
li{margin:5px 0}
.pass{color:var(--good);font-weight:700}
.fail{color:var(--bad);font-weight:700}
.toc{display:flex;flex-wrap:wrap;gap:8px;margin:26px 0 4px}
.toc a{font-size:13px;color:var(--muted);text-decoration:none;border:1px solid var(--line);
  border-radius:999px;padding:4px 12px}
.toc a:hover{color:var(--fg);border-color:var(--muted)}
footer{border-top:1px solid var(--line);margin-top:60px;padding-top:22px;font-size:13px;
  color:var(--muted)}
@media (max-width:640px){h1{font-size:29px}h2{font-size:21px}.wrap{padding:0 16px 64px}}
"""


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    meta = read_metadata()
    panel = jsn("m2_panel.json")
    qc = tsv("m1_sequence_qc.tsv")
    rank = tsv("m6_calibrated_ranking.tsv")
    suit = jsn("m6_system_suitability.json")
    val = jsn("m4_filter_validation.json")
    calib = jsn("m6_promiscuity_calibration.json")
    clusters = tsv("m5_clusters.tsv")
    epitopes = tsv("m5_epitopes.tsv")
    summary = {r["id"]: r for r in tsv("m5_ligand_summary.tsv")}
    coincide = tsv("m7_tb_coincidence.tsv")
    grid = tsv("m8_exposure_grid.tsv")
    deimm = tsv("m9_deimmunization_scan.tsv")
    binding = tsv("m3_binding_long.tsv")

    test_id = next(r["id"] for r in tsv("m5_ligand_summary.tsv")
                   if r["role"] == "test_article")
    test = summary[test_id]
    test_rank = next(r for r in rank if r["id"] == test_id)
    anchor = next(r for r in rank if r["id"] == cfg["benchmarks"]["anchor_low"])

    n_el = sum(1 for r in binding if r["call_el"] == "SB")
    n_cons = sum(1 for r in binding if r["call_consensus"] == "SB")
    test_clusters = [c for c in clusters if c["id"] == test_id]
    test_foreign = [c for c in test_clusters if c["tolerance_class"] != "all_tolerised"]
    top = max(test_foreign, key=lambda c: float(c["pop_presenting"])) if test_foreign else None

    H = []
    A = H.append

    # ------------------------------------------------------------------ hero
    A(f"""<header class="hero"><div class="kicker">Process-related impurity · immunogenicity risk assessment · RUO</div>
<h1>HLA-DR immunogenicity of an AAV affinity ligand</h1>
<p class="lede">A {len(seqs[test_id])}-residue camelid VHH affinity-chromatography ligand,
screened across a {panel['panel_size_total']}-molecule HLA-DR panel covering
{panel['weighted_coverage']*100:.1f}% of the weighted US/EU population, scored against
benchmark affinity ligands and assay controls run in the same batch.</p>
<div class="toc">
<a href="#verdict">Verdict</a><a href="#panel">Panel</a><a href="#suitability">Controls</a>
<a href="#promiscuity">Promiscuity scale</a><a href="#binding">Binding</a><a href="#tolerance">Tolerance filter</a>
<a href="#clusters">Epitopes</a><a href="#ada">B-cell / ADA</a>
<a href="#exposure">Exposure</a><a href="#deimm">Deimmunisation</a>
<a href="#limits">Limits &amp; next steps</a></div></header>""")

    A('<div class="wrap">')

    # --------------------------------------------------------------- verdict
    A('<h2 id="verdict">Verdict</h2>')
    band_cls = {"comparable-to-benchmark": "good", "modestly-elevated": "warn",
                "elevated": "bad", "high": "bad"}.get(test_rank["risk_band"], "acc")
    A(f"""<div class="grid g4">
<div class="stat"><div class="n">{float(test_rank['fold_vs_ProteinA_Z']):.1f}&times;</div>
  <div class="l">vs Protein A Z-domain</div>
  <div class="s">the affinity ligand with the longest controlled clinical leachate history</div></div>
<div class="stat"><div class="n">{test['n_foreign_epitopes']}</div>
  <div class="l">non-self epitopes</div>
  <div class="s">of {test['n_epitopes']} total; {test['n_tolerised_epitopes']} matched human self and
      {'was' if int(test['n_tolerised_epitopes']) == 1 else 'were'} down-weighted</div></div>
<div class="stat"><div class="n">{float(test_rank['pop_at_risk_pct']):.0f}%</div>
  <div class="l">US/EU population at risk</div>
  <div class="s">carry a DR molecule predicted to present at least one non-self epitope</div></div>
<div class="stat"><div class="n">{test['max_promiscuity']}<span style="font-size:17px;color:var(--muted)">/{panel['panel_size_total']}</span></div>
  <div class="l">peak promiscuity</div>
  <div class="s">DR molecules presenting the dominant epitope &mdash;
      {calib.get('test_peak_vs_universal_sb','?')}&times; what the tetanus universal epitopes reach on this panel</div></div>
</div>""")

    if top:
        A(f"""<div class="callout bad"><p><strong>The risk is concentrated, not diffuse.</strong>
One region — residues {top['start']}–{top['end']}, binding core
<code>{top['peak_core']}</code>, peak 15-mer <code>{top['peak_peptide']}</code> — is presented by
{top['union_sb_alleles']} of the {panel['panel_size_total']} DR molecules and reaches
{float(top['pop_presenting'])*100:.0f}% of the weighted US/EU population. That single peptide is what a
confirmatory assay should be built around; the rest of the sequence is unremarkable for a
non-human protein of this size.</p></div>""")

    A(f"""<p>On the intrinsic scale the ligand sits at
<strong>{float(test_rank['pIRS']):.2f}</strong> population-weighted non-self epitope units per
100 residues, against <strong>{float(anchor['pIRS']):.2f}</strong> for the Protein A Z-domain —
a <strong class="tag {band_cls}">{test_rank['risk_band']}</strong> result on the calibrated scale. Both numbers come from
the same batch, the same panel and the same thresholds, which is the only reason the comparison
means anything. Read on for what the batch controls say about whether the batch is reportable at
all, and for the exposure arithmetic that turns this into a dose-level call.</p>""")

    # ----------------------------------------------------------------- panel
    A('<h2 id="panel">The panel had to be designed, not inherited</h2>')
    A(f"""<p>The usual starting point is the HLA-DR subset of the IEDB class II reference set —
15 molecules, widely described as broadly representative. Measured against the IEDB
allele-frequency tables it is not: it reaches
<strong>{panel['legacy_per_population']['United States Caucasoid']*100:.1f}%</strong> US-Caucasian
and <strong>{panel['legacy_per_population']['Europe']*100:.1f}%</strong> European DRB1 phenotypic
coverage, <strong>{panel['legacy_weighted_coverage']*100:.1f}%</strong> on the weighted US/EU
composite. A 95–98% requirement is not met by that set at any size, because the missing coverage
sits in alleles it does not contain.</p>
<p>Growing the panel greedily against measured coverage instead reaches the target at
<strong>{panel['panel_size_drb1']} DRB1 molecules</strong>
({panel['weighted_coverage']*100:.2f}% weighted US/EU), plus the four DRB3/4/5 molecules for the
second DR molecule most people express —
<strong>{panel['panel_size_total']} DR molecules in total</strong>.</p>""")
    A(img("fig1_panel_coverage.png"))

    cov_rows = [{"pop": p.replace("United States ", "US "),
                 "new": f"{panel['per_population'][p]*100:.1f}%",
                 "old": f"{panel['legacy_per_population'][p]*100:.1f}%",
                 "d": f"+{(panel['per_population'][p]-panel['legacy_per_population'][p])*100:.1f}"}
                for p in panel["per_population"]]
    A(table(cov_rows, ["pop", "new", "old", "d"],
            ["population", "designed panel", "legacy 15-DR panel", "Δ points"], "num"))
    A(f"""<p class="muted"><small>Coverage is single-locus Hardy–Weinberg phenotypic frequency over
the IEDB tables, including the tool's renormalisation for populations whose allele frequencies sum
above 1 — reproducing <code>calculate_population_coverage.py</code> to two decimals. DRB3/4/5 carry
no frequencies in those tables, so they add presentation breadth without entering this arithmetic
and are never double-counted. US Asian coverage lands at
{panel['per_population']['United States Asian']*100:.0f}% and US Amerindian at
{panel['per_population']['United States Amerindian']*100:.0f}%: a US/EU-weighted objective does not
buy those populations, and reaching them needs DRB1*04:05, *12:02, *15:02 and *14:54 added
explicitly.</small></p>""")

    A('<h3>Panel composition</h3>')
    A('<p class="mono" style="font-size:13px;color:var(--muted)">'
      + " · ".join(a.replace("HLA-", "") for a in
                   panel["drb1_panel"] + panel["drb345_panel"]) + "</p>")

    # ----------------------------------------------------------- suitability
    A('<h2 id="suitability">The batch carries its own controls</h2>')
    A("""<p>A ligand run on its own produces a number nobody can size. This batch runs it beside
two positive controls (tetanus toxin universal T-helper epitope regions, which the panel must
find), two human self proteins (which the tolerance filter must suppress), and four benchmark
affinity ligands. The batch is reportable only if the controls behave.</p>""")
    A(table(suit.get("checks", []), ["check", "detail", "pass"],
            ["system suitability check", "observed", "result"],
            fmt={"pass": lambda r: ('<span class="pass">PASS</span>' if r["pass"]
                                    else '<span class="fail">FAIL</span>')}))
    A(f"""<p>Batch reportable: <strong class="{'pass' if suit.get('batch_valid') else 'fail'}">
{'yes' if suit.get('batch_valid') else 'no'}</strong>.</p>""")
    A(img("fig3_calibrated_ranking.png"))

    role_tag = {"test_article": "bad", "benchmark_ligand": "acc",
                "clinical_anchor": "vio", "class_comparator": "",
                "positive_control": "warn", "negative_control_self": "good"}
    A(table(rank, ["id", "role", "pIRS", "pIRS_raw", "fold_vs_ProteinA_Z",
                   "max_promiscuity", "pop_at_risk_pct", "risk_band"],
            ["sequence", "role", "pIRS", "pIRS unfiltered", "× Protein A Z",
             "peak promiscuity", "% pop at risk", "band"], "num",
            fmt={"role": lambda r: f'<span class="tag {role_tag.get(r["role"],"")}">'
                                   f'{r["role"].replace("_"," ")}</span>',
                 "max_promiscuity": lambda r: f'{r["max_promiscuity"]}/{panel["panel_size_total"]}',
                 "pop_at_risk_pct": lambda r: f'{float(r["pop_at_risk_pct"]):.1f}%',
                 "pIRS": lambda r: f'{float(r["pIRS"]):.2f}',
                 "pIRS_raw": lambda r: f'{float(r["pIRS_raw"]):.2f}',
                 "fold_vs_ProteinA_Z": lambda r: f'{float(r["fold_vs_ProteinA_Z"]):.2f}×',
                 "risk_band": lambda r: (f'<span class="tag {BAND_TAG.get(r["risk_band"], "")}">'
                                         f'{BAND_SHORT.get(r["risk_band"], r["risk_band"])}'
                                         f'</span>')}))

    # --------------------------------------------------- promiscuity anchor
    if calib:
        u = calib["universal_epitopes"]
        n = calib["panel_size"]
        A('<h2 id="promiscuity">How promiscuous is "promiscuous"?</h2>')
        A(f"""<p>Reporting that a peptide binds {test['max_promiscuity']} of
{n} DR molecules invites the question nobody usually answers: compared to what? The batch answers it
by carrying two peptides whose promiscuity is an experimental fact rather than a prediction — the
tetanus toxin universal T-helper epitopes p2 (830–844) and p30 (947–967), which drive CD4 responses
in most donors <em>in vitro</em> and are used as universal helper epitopes in peptide vaccines.</p>""")
        test_best = min((float(e["best_el_rank"]) for e in epitopes
                         if e["id"] == test_id), default=None)
        prom_rows = [{"ep": k.replace("_region", "").replace("TT_", "tetanus toxin "),
                      "sb": f"{v['n_sb']}/{n}", "wb": f"{v['n_wb']}/{n}",
                      "best": v["best_rank"], "hi": False}
                     for k, v in u.items()]
        prom_rows.append({
            "ep": f"{test_id} peak epitope ({top['peak_core'] if top else '-'})",
            "sb": f"{test['max_promiscuity']}/{n}", "wb": "—",
            "best": test_best if test_best is not None else "—", "hi": True})
        strong = lambda k: (lambda r: (f"<strong>{esc(r[k])}</strong>" if r["hi"]
                                       else esc(r[k])))
        A(table(prom_rows, ["ep", "sb", "wb", "best"],
                ["peptide", f"!DR molecules at EL %Rank &lt; {calib['sb_threshold']:g}",
                 f"!at %Rank &lt; {calib['wb_threshold']:g}", "best %Rank"], "num",
                fmt={c: strong(c) for c in ("ep", "sb", "wb", "best")}))
        A(f"""<div class="callout bad"><p><strong>The strong-binder tier does not reproduce the
textbook universal epitopes.</strong> p2 clears EL %Rank &lt; 1 on
{u['TT_p2_region']['n_sb']} of {n} DR molecules and p30 on {u['TT_p30_region']['n_sb']}; both reach
{calib['universal_epitope_ceiling_wb']}/{n} only at the weak-binder tier. So %Rank &lt; 1 is a
high-specificity, low-sensitivity criterion — and a peptide that <em>does</em> clear it across
{test['max_promiscuity']} molecules, as <code>{top['peak_core'] if top else '-'}</code> does, is
<strong>{calib['test_peak_vs_universal_sb']}× more promiscuous by this measure than the universal
epitopes are</strong>. That is a far more useful statement than "six strong binders".</p></div>""")
        A("""<p class="muted"><small>Read in the honest direction, this also bounds the method: a
predictor that ranks known universal epitopes as weak binders is not a sensitive detector of
promiscuity, so peptides below the strong-binder tier in this ligand are not cleared — they are
merely not flagged. That asymmetry is why the confirmatory assays below are scoped on the flagged
peptides <em>plus</em> the intact ligand.</small></p>""")

    # ------------------------------------------------------------- binding
    A('<h2 id="binding">Binding prediction: two heads, not one</h2>')
    A(f"""<p>Every 15-mer of every sequence was scored on both NetMHCIIpan heads across all
{panel['panel_size_total']} DR molecules. Eluted-ligand scoring alone called
<strong>{n_el}</strong> strong binders across the batch; requiring the binding-affinity head to
agree (BA %Rank &lt; {cfg['prediction']['ba_confirm_rank']:.0f}) leaves
<strong>{n_cons}</strong> — <strong>{100*(n_el-n_cons)/max(n_el,1):.0f}%</strong> of EL-only calls
are peptides that look presented but do not measurably bind. Every downstream number uses the
consensus call.</p>""")
    A(img("fig2_binding_landscape.png"))

    # ------------------------------------------------------------ tolerance
    A('<h2 id="tolerance">Framework hits are not risk</h2>')
    A(f"""<p>A VHH framework is close to human germline V. Most of its predicted DR binders have
9-mer cores that a human repertoire has already been tolerised against, and counting them makes
every antibody-derived ligand score the same. Each core here is compared to every 9-mer of the
human Swiss-Prot proteome: {val.get('real_hit_rate_9of9',0)*100:.1f}% are exact human 9-mers,
{(val.get('real_hit_rate_8of9',0)-val.get('real_hit_rate_9of9',0))*100:.1f}% are one substitution
away.</p>
<p>The threshold is validated rather than asserted. Re-running the identical screen on cores drawn
from shuffled versions of the same sequences — a null with the same composition and no real
self-similarity — gives a hit rate of
<strong>{val.get('null_hit_rate_8of9',0)*100:.1f}%</strong> at the 8/9 cut against
<strong>{val.get('real_hit_rate_8of9',0)*100:.1f}%</strong> for real cores, an enrichment of
<strong>{val.get('enrichment_8of9_real_over_null','n/a')}×</strong>. A 5-of-9 "TCR-face" screen,
the obvious shortcut, was tried first and discarded: with only five positions specified it matches
the human proteome by chance several times per query and flags nearly everything.</p>""")
    germ = next((r for r in rank if r["id"] == "HumanVH3_23_germline"), None)
    if germ:
        A(f"""<div class="callout bad"><p><strong>Why this filter is not optional.</strong> Human
germline VH3-23 — a sequence every human is tolerised to by construction — scores a peak promiscuity
of <strong>{germ['max_promiscuity']}/{panel['panel_size_total']}</strong> DR molecules on this panel,
the same as the test article, and an unfiltered pIRS of {float(germ['pIRS_raw']):.2f} against the test
article's {float(test_rank['pIRS_raw']):.2f}. Raw binder counts cannot tell a camelid VHH from the
human framework it resembles. After the filter the germline control falls to
{float(germ['pIRS']):.2f} and the ranking becomes readable.</p></div>""")
    A(f"""<div class="callout"><p>The filter is doing work, not deleting indiscriminately: the two
human self controls lose
{max(float(r['tolerance_drop_pct']) for r in rank if r['role']=='negative_control_self'):.0f}% of
their unfiltered score, while the test article loses
{float(test_rank['tolerance_drop_pct']):.0f}%.</p></div>""")

    # ------------------------------------------------------------- clusters
    A('<h2 id="clusters">Where the risk actually sits</h2>')
    A(table([c for c in test_clusters],
            ["start", "end", "peak_core", "peak_peptide", "peak_el_rank",
             "union_sb_alleles", "pop_presenting", "tolerance_class"],
            ["from", "to", "core", "peak 15-mer", "best EL %Rank", "DR molecules",
             "% US/EU presenting", "class"], "num",
            fmt={"peak_core": lambda r: f'<code>{r["peak_core"]}</code>',
                 "peak_peptide": lambda r: f'<code>{r["peak_peptide"]}</code>',
                 "pop_presenting": lambda r: f'{float(r["pop_presenting"])*100:.1f}%',
                 "tolerance_class": lambda r: f'<span class="tag '
                    f'{"bad" if r["tolerance_class"]=="foreign" else "warn" if r["tolerance_class"]=="mixed" else "good"}">'
                    f'{r["tolerance_class"].replace("_"," ")}</span>'}))

    te = sorted([e for e in epitopes if e["id"] == test_id],
                key=lambda e: -int(e["n_sb_alleles"]))[:10]
    A('<h3>Top epitopes of the test article</h3>')
    A(table(te, ["pos", "core", "peptide", "best_el_rank", "n_sb_alleles",
                 "pop_presenting", "tolerance_class", "sb_alleles"],
            ["pos", "core", "15-mer", "EL %Rank", "DR", "% pop", "class", "presenting molecules"],
            "num",
            fmt={"core": lambda r: f'<code>{r["core"]}</code>',
                 "peptide": lambda r: f'<code>{r["peptide"]}</code>',
                 "pop_presenting": lambda r: f'{float(r["pop_presenting"])*100:.1f}%',
                 "sb_alleles": lambda r: ('<span class="wrapmono muted">'
                                          + esc(r["sb_alleles"].replace(";", "  "))
                                          + "</span>")}))

    # ------------------------------------------------------------------ ADA
    A('<h2 id="ada">The measured endpoint is an antibody, not a T cell</h2>')
    A("""<p>An ADA assay detects antibodies. A B cell only class-switches with help from a CD4
T cell recognising a peptide from the same protein, so the regions worth carrying into wet-lab work
first are those where a non-self DR cluster and a predicted B-cell epitope coincide.</p>""")
    A(img("fig4_tb_coincidence.png"))
    tc = [c for c in coincide if c["id"] == test_id]
    if tc:
        A(table(sorted(tc, key=lambda c: -float(c["t_pop_presenting"])),
                ["t_cluster", "t_peak_core", "t_pop_presenting", "b_region",
                 "overlap_aa", "region_peptide"],
                ["T-cell cluster", "core", "% pop presenting", "B-cell region",
                 "overlap (aa)", "region sequence"], "num",
                fmt={"t_peak_core": lambda r: f'<code>{r["t_peak_core"]}</code>',
                     "t_pop_presenting": lambda r: f'{float(r["t_pop_presenting"])*100:.1f}%',
                     "region_peptide": lambda r: f'<code>{r["region_peptide"]}</code>'}))
    A("""<p class="muted"><small>Linear B-cell prediction (BepiPred-2.0) is the weakest model in
this pipeline — most real ADA epitopes are conformational. This layer prioritises regions; it is
never a standalone claim.</small></p>""")

    # ------------------------------------------------------------- exposure
    A('<h2 id="exposure">Intrinsic risk × dose</h2>')
    A(f"""<p>A leached ligand is an impurity, and impurity risk is intrinsic potential multiplied by
how much a patient receives. At the ligand's {cfg['exposure']['ligand_mw_kda']} kDa, the exposure
grid across plausible leachate levels and doses is:</p>""")
    band_tag = {"negligible": "good", "low": "good", "moderate": "warn", "elevated": "bad"}
    A(table(grid, ["leachate_ng_per_mg", "dose_mg", "ug_ligand_per_dose",
                   "nmol_per_dose", "exposure_band"],
            ["leachate (ng/mg)", "dose (mg)",
             '!ligand per dose (<span class="unit">µg</span>)',
             '!<span class="unit">nmol</span>/dose', "band"], "num",
            fmt={"ug_ligand_per_dose": lambda r: f'{float(r["ug_ligand_per_dose"]):,.3f}',
                 "nmol_per_dose": lambda r: f'{float(r["nmol_per_dose"]):.4f}',
                 "exposure_band": lambda r: f'<span class="tag {band_tag.get(r["exposure_band"],"")}">'
                                            f'{r["exposure_band"]}</span>'}))
    A("""<p class="muted"><small>Bands are an internal triage convention, not a regulatory
threshold. No agency publishes a numeric leachate immunogenicity limit; the ICH Q6B / EMA
expectation is control to a justified, consistently achieved level. Replace the grid in
<code>config/config.yaml</code> with your product's actual leachate data and dose.</small></p>""")

    # ----------------------------------------------------------- deimm
    if deimm:
        wt = next(r for r in deimm if r["variant"] == "WT")
        A('<h2 id="deimm">If the ligand can be re-engineered</h2>')
        muts = [r for r in deimm if r["variant"] != "WT"]
        knock = [r for r in muts if float(r["pop_presenting"]) == 0]
        cons = [r for r in muts if int(r["blosum62"]) >= 0]
        # A free cysteine is a disulfide-scrambling liability in an engineered
        # binder, so it is not offered as a design candidate even when it scores.
        cons_designable = [r for r in cons if not r["variant"].endswith("C")]
        best_cons = sorted(cons_designable, key=lambda r: float(r["pop_presenting"]))[:1]
        best_knock_bl = max((int(r["blosum62"]) for r in knock), default=None)
        A(f"""<p>Every one of the 19 substitutions was placed at each MHC-II anchor pocket (P1, P4,
P6, P9) of the dominant core and the mutated 15-mer re-scored across the panel — {len(muts)}
variants in all. Wild type is presented by {wt['n_sb_alleles']} DR molecules reaching
{float(wt['pop_presenting'])*100:.0f}% of the weighted population.</p>""")
        A(img("fig5_deimmunization.png"))
        if knock and best_cons:
            bc = best_cons[0]
            bcbl = f"{int(bc['blosum62']):+d}" if int(bc["blosum62"]) else "0"
            A(f"""<div class="callout"><p><strong>The knockouts and the designable changes are not
the same substitutions.</strong> {len(knock)} variants abolish predicted presentation entirely, but
every one of them is chemically disruptive — the best BLOSUM62 score among them is
{best_knock_bl:+d}, and they sit at buried framework positions where that is a real stability risk.
The most conservative useful change is <code>{bc['variant']}</code>
({bc['substitution']}, BLOSUM62 {bcbl} — cysteine substitutions are excluded as
disulfide-scrambling liabilities), which cuts presentation from
{float(wt['pop_presenting'])*100:.0f}% to {float(bc['pop_presenting'])*100:.1f}% without a
chemically radical substitution. That is the trade-off to take to a stability screen — not a
knockout that would probably cost the fold.</p></div>""")
        A('<h3>Conservative substitutions (BLOSUM62 ≥ 0)</h3>')
        A(table(sorted(cons, key=lambda r: float(r["pop_presenting"]))[:8],
                ["variant", "substitution", "core", "n_sb_alleles", "pop_presenting",
                 "delta_pop_presenting", "blosum62", "germline_residue"],
                ["variant", "substitution", "core", "DR", "% pop", "Δ % pop",
                 "BLOSUM62", "germline aa"], "num",
                fmt={"core": lambda r: f'<code>{r["core"]}</code>',
                     "pop_presenting": lambda r: f'{float(r["pop_presenting"])*100:.1f}%',
                     "delta_pop_presenting": lambda r: f'{float(r["delta_pop_presenting"])*100:+.1f}',
                     "blosum62": lambda r: (f'{int(r["blosum62"]):+d}'
                                            if int(r["blosum62"]) else "0")}))
        A(f"""<p class="muted"><small>This scan uses the eluted-ligand head alone, so wild-type
breadth reads {wt['n_sb_alleles']}/{panel['panel_size_total']} here against
{test['max_promiscuity']}/{panel['panel_size_total']} on the EL+BA consensus call used everywhere
else — the comparison between variants is what carries the meaning, not the absolute breadth. The
"% pop" column counts DRB1 molecules only, so a variant can read 0.0% while one or two DRB3/4/5
molecules still bind it. And this models DR presentation only: effect on ligand–target affinity,
resin capacity and alkaline stability is not modelled, and any candidate has to go back through
binding and stability screens before it means anything.</small></p>""")

    # ---------------------------------------------------------------- limits
    A('<h2 id="limits">What this is, and what it is not</h2>')
    A("""<p>In-silico DR screening ranks and localises risk. It does not measure it. pIRS is a
relative scale readable only against the benchmarks in this batch; it is not a predicted ADA
incidence, and no in-silico method available today predicts one. The specific limits that matter
for how far this output can be pushed:</p>
<ul>
<li><strong>Prediction, not presentation.</strong> NetMHCIIpan scores peptide–MHC binding. It does
not model antigen uptake, endosomal proteolysis, HLA-DM editing or the ligand's conformational
stability in the endosome — all of which decide what is actually presented.</li>
<li><strong>DR only.</strong> DP and DQ contribute to CD4 responses and are excluded here by the
sponsor's panel specification. DQ in particular is implicated in several biologic ADA
responses.</li>
<li><strong>The tolerance filter is a screen, not JanusMatrix.</strong> It does not require the
human counterpart peptide to bind the same allele, so it errs toward calling more peptides
tolerised. Every flagged core is written out with the human protein it matched, so each call is
checkable.</li>
<li><strong>No aggregation, no adjuvant effect, no repeat-dose modelling.</strong> Aggregated or
particulate impurity is substantially more immunogenic than the monomer, and none of that is in
scope for a sequence-based method.</li>
<li><strong>The benchmark anchor is an argument, not a measurement.</strong> Protein A Z-domain
leachate has a long controlled clinical history, but "comparable predicted epitope content to
Protein A" is an inference about risk, not evidence of it.</li>
</ul>""")

    A('<h3>Recommended confirmatory work (all available RUO)</h3>')
    A(f"""<ol>
<li><strong>HLA-DR competitive binding</strong> on the {len(test_foreign)} flagged peptides against
the panel's dominant molecules — directly confirms the predicted binding, days not months, lowest
cost.</li>
<li><strong>MAPPs</strong> (MHC-associated peptide proteomics) on monocyte-derived dendritic cells
from HLA-typed donors — shows what is actually processed and presented from the intact ligand,
which no predictor can.</li>
<li><strong>Ex-vivo PBMC / CD4 T-cell proliferation</strong> across ~50 HLA-typed donors matched to
this panel, using the flagged peptides plus whole ligand — the closest available surrogate for
clinical ADA risk.</li>
</ol>
<p>All three are scoped by the peptide list this pipeline produces, which is the practical point of
running it: it turns "assess the ligand" into a defined, costed experiment on
{len(test_foreign)} peptides.</p>""")

    # ---------------------------------------------------------------- inputs
    A('<hr class="rule">')
    A('<h2 id="inputs">Inputs and provenance</h2>')
    A(f"""<p>The test article is a <strong>public stand-in</strong> — the AAVX affinity ligand from
PDB 9DC3 — so this demonstration is reproducible end to end. Proprietary ligands drop into
<code>data/sequences.fasta</code> as de-identified FASTA and every module runs unchanged; nothing
in the pipeline requires a database identifier for the test article.</p>""")
    rows = []
    for r in qc:
        rows.append({**r, "note": meta[r["id"]]["note"]})
    A(table(rows, ["id", "role", "length", "mw_kda", "pI",
                   "pct_id_human_IGHV3_23", "vhh_hallmarks", "source"],
            ["sequence", "role", "aa", "kDa", "pI", "% id human IGHV3-23",
             "VHH hallmark tetrad", "source"], "num",
            fmt={"role": lambda r: f'<span class="tag {role_tag.get(r["role"],"")}">'
                                   f'{r["role"].replace("_"," ")}</span>',
                 "pct_id_human_IGHV3_23": lambda r: f'{float(r["pct_id_human_IGHV3_23"]):.1f}%',
                 "source": lambda r: f'<small>{esc(r["source"])}</small>'}))
    A("""<p class="muted"><small>The hallmark tetrad is FR2 Kabat 37/44/45/47: camelid VHHs carry
F/Y–E–R–G/L where a human VH carries V–G–L–W. It separates a humanised VHH from a native camelid
one at a glance, and that distinction is a better predictor of clinical ADA than epitope count
alone.</small></p>""")

    A(f"""<footer><p>Generated by the <code>hla_ligand_immunogenicity</code> pipeline.
Predictions: NetMHCIIpan (EL + BA heads) and BepiPred-2.0 via the IEDB REST API.
Allele frequencies: IEDB population-coverage tables v3.0.2.
Human proteome: UniProt Swiss-Prot <em>Homo sapiens</em>.
Sequences: RCSB PDB and UniProt, accessions in the table above.
Research use only — not for regulatory submission without confirmatory wet-lab data.</p></footer>""")
    A("</div>")

    html = (f"<title>AAV Ligand Immunogenicity</title>\n<style>{CSS}</style>\n"
            + "\n".join(H))
    out = os.path.join(ROOT, "report.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"wrote {out} ({len(html)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
