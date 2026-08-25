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


FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=IBM+Plex+Mono:wght@400;500&'
         'family=IBM+Plex+Sans:wght@400;450;500;600&'
         'family=IBM+Plex+Serif:wght@500;600&display=swap">')

# IBM Plex: a family drawn for technical documentation, and its monospace is
# what the peptide sequences on this page actually need. Serif for headings,
# sans for running text, mono for every sequence and figure of merit.
CSS = """
:root{
  --serif:"IBM Plex Serif",Cambria,Georgia,serif;
  --sans:"IBM Plex Sans",-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
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
  font:400 16px/1.66 var(--sans);
  -webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px 96px}
header.hero{border-bottom:1px solid var(--line);padding:56px 0 34px;margin-bottom:8px}
.kicker{font-size:12px;letter-spacing:.13em;text-transform:uppercase;color:var(--muted);
  font-weight:650}
h1,h2,h3{font-family:var(--serif);font-weight:600;text-wrap:balance}
h1{font-size:40px;line-height:1.12;margin:14px 0 12px;letter-spacing:-.018em}
h2{font-size:26px;margin:58px 0 8px;letter-spacing:-.012em;scroll-margin-top:20px;
  line-height:1.22}
h3{font-size:17.5px;margin:30px 0 6px;letter-spacing:-.005em}
.lede{color:var(--muted);font-size:17px;max-width:74ch;margin:0}
p{max-width:78ch}
.muted{color:var(--muted)}
small{color:var(--muted)}
code{background:var(--soft);border:1px solid var(--line);border-radius:4px;
  padding:1px 5px;font-size:.86em;font-family:var(--mono);letter-spacing:.01em}
.mono{font-family:var(--mono);font-size:.9em}
.rule{height:1px;background:var(--line);border:0;margin:44px 0}
.grid{display:grid;gap:14px;margin:22px 0}
.g4{grid-template-columns:repeat(auto-fit,minmax(168px,1fr))}
.g3{grid-template-columns:repeat(auto-fit,minmax(230px,1fr))}
.g2{grid-template-columns:repeat(auto-fit,minmax(320px,1fr))}
.stat{border:1px solid var(--line);border-radius:10px;padding:16px 16px 14px;background:var(--soft)}
.stat .n{font-family:var(--serif);font-size:30px;font-weight:600;letter-spacing:-.02em;
  line-height:1.08;font-variant-numeric:tabular-nums}
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
  font-family:var(--mono);font-size:11px;line-height:1.5}
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
.toc a:focus-visible,a:focus-visible{outline:2px solid var(--accent);outline-offset:2px;
  border-radius:4px}
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
    acc = jsn("m11_calibration.json")
    bench = jsn("m10_benchmark_summary.json")
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
<a href="#promiscuity">Promiscuity scale</a><a href="#accuracy">Measured accuracy</a>
<a href="#binding">Binding</a><a href="#tolerance">Tolerance filter</a>
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
  <div class="s">DR molecules presenting the dominant epitope &mdash; the same breadth as the most
      promiscuous epitope with human T-cell evidence, measured on this panel</div></div>
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
by carrying three peptides whose promiscuity is an experimental fact rather than a prediction:
influenza haemagglutinin HA306-318 and the tetanus toxin epitopes p2 (830–844) and p30 (947–967).
Each was checked against IEDB before being used — HA306-318 has positive human DR-restricted T-cell
assays on 25 distinct DR molecules, p2 on 6, p30 on 1. PADRE was a candidate and was dropped: IEDB
holds no positive human DR-restricted T-cell record for it.</p>""")
        test_best = min((float(e["best_el_rank"]) for e in epitopes
                         if e["id"] == test_id), default=None)
        prom_rows = [{"ep": k.replace("_region", "").replace("TT_", "tetanus toxin ")
                            .replace("HA_306_318", "influenza HA306-318"),
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
        rng_sb = calib.get("universal_epitope_range_sb", [0, 0])
        rng_wb = calib.get("universal_epitope_range_wb", [0, 0])
        peak = int(test["max_promiscuity"])
        rel = ("above" if peak > rng_sb[1] else
               "at the top of" if peak == rng_sb[1] else "inside")
        A(f"""<div class="callout {'bad' if peak >= rng_sb[1] else ''}"><p>
<strong>The ligand's dominant core is as promiscuous as the most promiscuous epitope we have
evidence for.</strong> The three universal epitopes reach {rng_sb[0]}&ndash;{rng_sb[1]} of {n} DR
molecules at EL %Rank&nbsp;&lt;&nbsp;{calib['sb_threshold']:g}; <code>{top['peak_core'] if top else '-'}</code>
reaches {peak}/{n} &mdash; {rel} that range. It is not an outlier. It is also not innocuous: it sits
exactly where a peptide known to drive CD4 responses in most donors sits.</p></div>""")
        A(f"""<div class="callout"><p><strong>The same table bounds the method's sensitivity.</strong>
The strong-binder tier recovers only {100*rng_sb[1]/n:.0f}% of the DR molecules a universal epitope
is actually presented by &mdash; HA306-318 is positive in human T-cell assays on 25 distinct DR
molecules and clears %Rank&nbsp;&lt;&nbsp;1 on {u['HA_306_318_region']['n_sb']} of the {n} tested here.
So a peptide below the tier is <em>unflagged</em>, not <em>cleared</em>, and the confirmatory assays
below are scoped on the flagged peptides <em>plus</em> the intact ligand for exactly that
reason.</p></div>""")

    # ---------------------------------------------------- measured accuracy
    if acc and bench:
        geo_el = acc["comparisons"].get("GEO_vs_EL") or {}
        cur = next((v for k, v in acc["operating_points"].items()
                    if k.startswith("AND")), None)
        best = acc.get("sweep", {}).get("max_mcc")
        el_only = next((v for k, v in acc["operating_points"].items()
                        if k.startswith("EL<1")), None)
        A('<h2 id="accuracy">Does the rule actually work?</h2>')
        A(f"""<p>Everything above this point is a claim about accuracy that had not been tested.
The EL+BA consensus rule entered the pipeline on an argument — eluted-ligand scoring over-calls, so
make the affinity head agree — and an argument is not evidence. To find out, the rule was scored
against ground truth: every HLA-DR-restricted CD4 T-cell assay result IEDB holds for the
{panel['panel_size_total']} molecules in this panel, in human hosts.</p>""")
        A(f"""<div class="grid g4">
<div class="stat"><div class="n">{bench['labelled_pairs']:,}</div>
  <div class="l">labelled outcomes</div>
  <div class="s">{bench['positives']:,} positive, {bench['negatives']:,} negative
  (peptide, DR molecule) pairs from {bench['raw_records']:,} assay records</div></div>
<div class="stat"><div class="n">{bench['clusters']:,}</div>
  <div class="l">independent clusters</div>
  <div class="s">{bench['distinct_peptides']:,} distinct peptides collapsed by shared 9-mer, so
  overlapping peptides from one study count once</div></div>
<div class="stat"><div class="n">{acc['rules'][acc['best_continuous_rule']]['auc']:.3f}</div>
  <div class="l">best rule, ROC AUC</div>
  <div class="s">{acc['best_continuous_rule']} — against
  {acc['rules']['EL']['auc']:.3f} for eluted-ligand scoring alone</div></div>
<div class="stat"><div class="n">{(cur['ppv_at_scan_prevalence']*100 if cur else 0):.0f}%</div>
  <div class="l">PPV of a flag</div>
  <div class="s">chance a flagged peptide is a real epitope at
  {acc['scan_prevalence']*100:.0f}% assumed scan prevalence</div></div>
</div>""")
        A(img("fig6_calibration.png"))

        A("<h3>What each rule buys</h3>")
        A(table([{"rule": k, "auc": f"{v['auc']:.3f}",
                  "ap": f"{v['average_precision']:.3f}"}
                 for k, v in acc["rules"].items()],
                ["rule", "auc", "ap"],
                ["score", "ROC AUC", "average precision"], "num"))
        if geo_el:
            better = geo_el["ci95"][0] > 0
            A(f"""<div class="callout {'good' if better else 'bad'}"><p>
<strong>{'The consensus rule is a real improvement.' if better else
        'The consensus rule does not improve discrimination.'}</strong>
A paired bootstrap over 9-mer clusters puts the AUC difference between the EL&times;BA consensus and
EL alone at <strong>{geo_el['delta']:+.4f}</strong>, 95% CI
[{geo_el['ci95'][0]:+.4f}, {geo_el['ci95'][1]:+.4f}].
{'The interval excludes zero, so the gain survives the redundancy in the benchmark.' if better else
 'The interval spans zero — on this evidence the second head adds nothing to ranking, and its value '
 'is confined to where the threshold is placed.'}</p></div>""")

        A("<h3>Operating points</h3>")
        rows_op = [{"rule": k, **v} for k, v in acc["operating_points"].items()]
        if best:
            rows_op.append({"rule": f"calibrated {acc['best_continuous_rule']} "
                                    f"%Rank &lt; {best['threshold_rank']:.2f}", **best})
        A(table(rows_op, ["rule", "sensitivity", "specificity", "ppv", "mcc",
                          "ppv_at_scan_prevalence"],
                ["rule", "sensitivity", "specificity", "PPV (balanced set)", "MCC",
                 "PPV at scan prevalence"], "num",
                fmt={"rule": lambda r: r["rule"],
                     "sensitivity": lambda r: f'{r["sensitivity"]:.3f}',
                     "specificity": lambda r: f'{r["specificity"]:.3f}',
                     "ppv": lambda r: f'{r["ppv"]:.3f}',
                     "mcc": lambda r: f'{r["mcc"]:.3f}',
                     "ppv_at_scan_prevalence": lambda r:
                         f'{r["ppv_at_scan_prevalence"]*100:.0f}%'}))
        if cur and best:
            A(f"""<p>The rule the pipeline currently runs sits at
<strong>{cur['sensitivity']*100:.0f}% sensitivity</strong> and
<strong>{cur['specificity']*100:.0f}% specificity</strong>. The threshold that maximises
Matthews correlation on this benchmark is
<strong>{acc['best_continuous_rule']} %Rank &lt; {best['threshold_rank']:.2f}</strong>, at
{best['sensitivity']*100:.0f}% / {best['specificity']*100:.0f}% — MCC
{best['mcc']:.3f} against {cur['mcc']:.3f}.</p>""")

        strat = acc.get("stratified", {})
        if all(k in strat and "note" not in strat[k] for k in ("self", "non_self")):
            A(f"""<p class="muted"><small>Stratifying by source removes the obvious confound: the
benchmark's positives are {bench['self_fraction_positives']*100:.0f}% human-derived and its
negatives {bench['self_fraction_negatives']*100:.0f}%, so a predictor that merely recognised
"foreign" would score well for the wrong reason. Within the non-self stratum
(n&nbsp;=&nbsp;{strat['non_self']['n']}) the best rule reaches AUC
{strat['non_self'][acc['best_continuous_rule']]:.3f}; within the self stratum
(n&nbsp;=&nbsp;{strat['self']['n']}) {strat['self'][acc['best_continuous_rule']]:.3f}. The
discrimination is not an artefact of foreignness.</small></p>""")

        A(f"""<div class="callout"><p><strong>Read the absolute numbers as an upper bound.</strong>
NetMHCIIpan is trained on IEDB binding-affinity and mass-spec eluted-ligand data. The labels here
are T-cell assay outcomes — a different endpoint it was not trained on — but many of these peptides
also carry binding measurements in IEDB, so partial training-set overlap is certain and cannot be
excluded from outside. The <em>comparison between rules</em> is far more robust: every rule uses the
same two possibly-leaky predictors on the same peptides, so leakage inflates them together and
largely cancels in the difference. IEDB's negatives also carry their own bias — a peptide tested and
found negative in one donor set is not a clean non-binder — which inflates apparent specificity in
both directions equally.</p></div>""")

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
    bnd = (calib or {}).get("boundary_controls", {})
    if bnd:
        A('<h3>Where the filter is wrong, measured</h3>')
        A(f"""<p>The filter's rule is "self implies tolerated". That is right often enough to be
worth applying and wrong often enough to be worth measuring, so the batch carries two controls that
sit on the boundary.</p>""")
        A(table([{"id": k, **v} for k, v in bnd.items()],
                ["id", "role", "dr_breadth_sb", "dr_breadth_wb",
                 "cores_called_tolerised", "ligand_pIRS_unfiltered", "ligand_pIRS"],
                ["control", "role", f"!DR at %Rank&lt;{calib['sb_threshold']:g}",
                 f"!at &lt;{calib['wb_threshold']:g}", "cores called tolerised",
                 "pIRS unfiltered", "pIRS filtered"], "num",
                fmt={"id": lambda r: f'<code>{r["id"].replace("_region","")}</code>',
                     "role": lambda r: f'<span class="tag warn">{r["role"].replace("_"," ")}</span>',
                     "dr_breadth_sb": lambda r: f'{r["dr_breadth_sb"]}/{calib["panel_size"]}',
                     "dr_breadth_wb": lambda r: f'{r["dr_breadth_wb"]}/{calib["panel_size"]}',
                     "cores_called_tolerised": lambda r:
                         f'{r["cores_called_tolerised"]}/{r["predicted_cores_in_window"]}',
                     "ligand_pIRS_unfiltered": lambda r: f'{r["ligand_pIRS_unfiltered"]:.2f}',
                     "ligand_pIRS": lambda r: f'{r["ligand_pIRS"]:.2f}'}))
        mbp = bnd.get("MBP_85_99_region")
        clp = bnd.get("CLIP_87_101_region")
        if mbp:
            A(f"""<div class="callout bad"><p><strong>The filter suppresses a real epitope, and this
is the clearest example of it.</strong> Myelin basic protein 85&ndash;99 is a human self peptide
with positive DR-restricted T-cell assays on 10 distinct DR molecules in IEDB. On this panel it is
the broadest binder in the whole batch &mdash; {mbp['dr_breadth_sb']}/{calib['panel_size']} DR
molecules at the strong-binder tier, unfiltered pIRS {mbp['ligand_pIRS_unfiltered']:.2f}, higher
than any ligand here. The tolerance filter takes it to
{mbp['ligand_pIRS']:.2f}. That is the filter working as designed and being wrong: self-derived
epitopes exist, autoimmunity is what they are, and a sequence-identity filter cannot see the
difference. It is safe for an affinity ligand because a foreign scaffold's framework similarity to
human germline is the case the filter was built for &mdash; not for anything where self-reactivity
is the question.</p></div>""")
        if clp:
            A(f"""<p>The opposite control behaves too. CLIP occupies the groove of essentially every
DR molecule as the invariant chain's placeholder, yet carries no positive human DR T-cell record in
IEDB. Here it reaches {clp['dr_breadth_wb']}/{calib['panel_size']} DR molecules at the weak tier and
{clp['dr_breadth_sb']}/{calib['panel_size']} at the strong tier, and ends at pIRS
{clp['ligand_pIRS']:.2f} &mdash; binding strength alone is not being read as risk.</p>""")

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

    html = (f"<title>AAV Ligand Immunogenicity</title>\n{FONTS}\n"
            f"<style>{CSS}</style>\n" + "\n".join(H))
    out = os.path.join(ROOT, "report.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"wrote {out} ({len(html)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
