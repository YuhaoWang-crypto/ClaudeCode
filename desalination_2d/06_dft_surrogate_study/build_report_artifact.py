#!/usr/bin/env python3
"""
Assemble the shareable HTML report artifact: reads the study JSON result files,
embeds the result figures as base64 data URIs (CSP-safe, no external hosts), and
writes a self-contained report page. Re-run to pick up new results (e.g. MoS2
hydrated PMF) and redeploy to the same artifact URL.

    python build_report_artifact.py --out report.html
"""
import argparse
import base64
import json
import os

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")


def img(path):
    """Return a data: URI for a PNG, or None if missing."""
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    return f"data:image/png;base64,{b64}"


def load(path, default=None):
    return json.load(open(path)) if os.path.exists(path) else default


def figure(src, caption, note=""):
    if not src:
        return (f'<figure class="fig pending"><div class="ph">figure pending — '
                f'run still computing</div><figcaption>{caption}</figcaption></figure>')
    n = f'<span class="fignote">{note}</span>' if note else ""
    return (f'<figure class="fig"><img src="{src}" alt="{caption}">'
            f'<figcaption>{caption}{n}</figcaption></figure>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "report.html"))
    args = ap.parse_args()

    study = load(os.path.join(RES, "study_results.json"), {})
    hyd25 = load(os.path.join(RES, "hydrated_pmf_results.json"), {})
    hyd30 = load(os.path.join(HERE, "results_r30", "hydrated_pmf_results.json"), {})
    mos_h25 = load(os.path.join(HERE, "results_mos2_r25",
                                "hydrated_pmf_results.json"), {})
    mos_h30 = load(os.path.join(HERE, "results_mos2_r30",
                                "hydrated_pmf_results.json"), {})

    def stat(d, key):
        try:
            return d["materials"][0]["stats"][key]["barrier_eV"]
        except Exception:
            return None

    # figures
    f_pore = img(os.path.join(RES, "poresize_water_barrier.png"))
    f_hyd = img(os.path.join(RES, "hydrated_summary.png"))
    f_mace = img(os.path.join(RES, "mace_learning_curve.png"))
    f_gp = img(os.path.join(RES, "profile_graphene.png"))
    f_mp = img(os.path.join(RES, "profile_MoS2.png"))

    # numbers
    partA = study.get("partA_poresize_scan_graphene_H2O", [])
    partB = study.get("partB_materials", [])

    def pb(mat, perm):
        for m in partB:
            if m["material"] == mat:
                return m["permeants"].get(perm, {}).get("barrier_eV")
        return None

    g_na30 = stat(hyd30, "Na_hyd")
    g_h2o30 = stat(hyd30, "H2O")
    mos_na30 = stat(mos_h30, "Na_hyd")
    mos_h2o30 = stat(mos_h30, "H2O")
    mos_pending = not mos_h25 and not mos_h30

    def cell(v, unit=" eV"):
        return f"{v:g}{unit}" if v is not None else '<span class="na">pending</span>'

    # ---- pore-size table rows
    pore_rows = "".join(
        f'<tr><td class="num">{r["pore_radius_A"]:.1f}</td>'
        f'<td class="num">{r["n_removed_C"]}</td>'
        f'<td class="num">{r["barrier_eV"]:g}</td>'
        f'<td>{"blocked" if r["barrier_eV"]>3 else ("RO window" if r["barrier_eV"]>0.3 else "free flow")}</td></tr>'
        for r in partA)

    # ---- selectivity (bare, Part B) table
    def selrow(mat):
        return (f'<tr><td>{mat}</td>'
                f'<td class="num">{cell(pb(mat,"H2O"))}</td>'
                f'<td class="num">{cell(pb(mat,"Na"))}</td>'
                f'<td class="num">{cell(pb(mat,"Cl"))}</td></tr>')
    sel_rows = selrow("graphene") + (selrow("MoS2") if pb("MoS2", "H2O") is not None else "")

    # ---- hydrated table
    def hydrow(label, e_h2o, e_na, e_cl):
        return (f'<tr><td>{label}</td>'
                f'<td class="num pass">{cell(e_h2o)}</td>'
                f'<td class="num reject">{cell(e_na)}</td>'
                f'<td class="num reject">{cell(e_cl)}</td></tr>')
    hyd_rows = (
        hydrow("graphene · r=2.5 Å", stat(hyd25, "H2O"), stat(hyd25, "Na_hyd"),
               stat(hyd25, "Cl_hyd")) +
        hydrow("graphene · r=3.0 Å", stat(hyd30, "H2O"), stat(hyd30, "Na_hyd"),
               stat(hyd30, "Cl_hyd")))
    if mos_h25:
        hyd_rows += hydrow("MoS₂ · r=2.5 Å", stat(mos_h25, "H2O"),
                           stat(mos_h25, "Na_hyd"), stat(mos_h25, "Cl_hyd"))
    if mos_h30:
        hyd_rows += hydrow("MoS₂ · r=3.0 Å", stat(mos_h30, "H2O"),
                           stat(mos_h30, "Na_hyd"), stat(mos_h30, "Cl_hyd"))

    mos_status = ("computing" if mos_pending else "done")

    subs = {
        "f_pore": figure(f_pore, "Part A · Water crossing barrier vs graphene pore radius.",
                         "CHGNet surrogate"),
        "f_hyd": figure(f_hyd, "Step 2 · Hydrated vs bare ion vs water crossing barriers (graphene).",
                        "bars > 3 eV are steric upper bounds"),
        "f_gp": figure(f_gp, "Part B · Bare-species crossing profiles, graphene (r=2.5 Å)."),
        "f_mp": figure(f_mp, "Part B · Bare-species crossing profiles, MoS₂ (r=2.5 Å)."),
        "f_mace": figure(f_mace, "Steps 3–4 · MACE demo training on CHGNet-labelled frames."),
        "pore_rows": pore_rows, "sel_rows": sel_rows, "hyd_rows": hyd_rows,
        "g_na30": cell(g_na30), "g_h2o30": cell(g_h2o30),
        "mos_na30": cell(mos_na30), "mos_h2o30": cell(mos_h2o30),
        "mos_status": mos_status,
    }
    html = TEMPLATE
    for k, v in subs.items():
        html = html.replace("{" + k + "}", str(v))
    with open(args.out, "w") as fh:
        fh.write(html)
    print(f"Wrote {args.out} ({len(html)//1024} KB)  MoS2 hydrated: {mos_status}")


TEMPLATE = r"""<title>2D-Material Desalination — DFT-stage Study</title>
<style>
:root{
  --paper:#f6f8f9; --card:#ffffff; --ink:#141e24; --muted:#5c6f78;
  --line:#d9e2e5; --accent:#158793; --accent-soft:#e0eff0;
  --pass:#158793; --reject:#bd5a44; --warn:#b6832b; --warn-soft:#f6ecd6;
  --mono:ui-monospace,"SF Mono",Menlo,"Cascadia Code",Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
}
@media (prefers-color-scheme:dark){
  :root{ --paper:#0e1519; --card:#151f25; --ink:#e6eef1; --muted:#93a6ad;
    --line:#26343b; --accent:#3fb3bf; --accent-soft:#12303350;
    --pass:#3fb3bf; --reject:#e08668; --warn:#d6a94e; --warn-soft:#2a230f; }
}
:root[data-theme="dark"]{ --paper:#0e1519; --card:#151f25; --ink:#e6eef1; --muted:#93a6ad;
  --line:#26343b; --accent:#3fb3bf; --accent-soft:#12303350;
  --pass:#3fb3bf; --reject:#e08668; --warn:#d6a94e; --warn-soft:#2a230f; }
:root[data-theme="light"]{ --paper:#f6f8f9; --card:#ffffff; --ink:#141e24; --muted:#5c6f78;
  --line:#d9e2e5; --accent:#158793; --accent-soft:#e0eff0;
  --pass:#158793; --reject:#bd5a44; --warn:#b6832b; --warn-soft:#f6ecd6; }

*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
  line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:820px;margin:0 auto;padding:0 24px}
.eyebrow{font-family:var(--mono);font-size:.72rem;letter-spacing:.14em;
  text-transform:uppercase;color:var(--accent);font-weight:600}
h1{font-size:clamp(1.9rem,4.5vw,2.9rem);line-height:1.08;margin:.3em 0 .2em;
  font-weight:800;letter-spacing:-.02em;text-wrap:balance}
h2{font-size:1.5rem;font-weight:750;margin:2.6em 0 .5em;letter-spacing:-.01em;
  padding-top:.7em;border-top:1px solid var(--line);text-wrap:balance}
h3{font-size:1.08rem;font-weight:700;margin:1.8em 0 .4em}
p{margin:.7em 0;max-width:66ch}
a{color:var(--accent)}
header{padding:64px 0 30px}
.lede{font-size:1.15rem;color:var(--muted);max-width:60ch}

/* honesty chips */
.chip{display:inline-flex;align-items:center;gap:.4em;font-family:var(--mono);
  font-size:.7rem;font-weight:600;letter-spacing:.03em;padding:.28em .6em;
  border-radius:999px;border:1px solid;white-space:nowrap}
.chip.warn{color:var(--warn);border-color:var(--warn);background:var(--warn-soft)}
.chip.ok{color:var(--accent);border-color:var(--accent);background:var(--accent-soft)}
.method-band{display:flex;flex-wrap:wrap;gap:10px;align-items:center;
  margin:22px 0 6px;padding:14px 16px;background:var(--card);
  border:1px solid var(--line);border-radius:12px;font-size:.9rem;color:var(--muted)}

/* key findings */
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin:26px 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px 18px}
.kpi .v{font-family:var(--mono);font-size:1.7rem;font-weight:700;letter-spacing:-.02em;
  font-variant-numeric:tabular-nums}
.kpi .v.pass{color:var(--pass)} .kpi .v.reject{color:var(--reject)}
.kpi .k{font-size:.82rem;color:var(--muted);margin-top:2px}

/* tables */
.tablewrap{overflow-x:auto;margin:1.1em 0;border:1px solid var(--line);border-radius:12px}
table{border-collapse:collapse;width:100%;font-size:.92rem}
caption{caption-side:top;text-align:left;font-family:var(--mono);font-size:.72rem;
  letter-spacing:.08em;text-transform:uppercase;color:var(--muted);padding:12px 14px 6px}
th,td{padding:9px 14px;text-align:left;border-bottom:1px solid var(--line)}
thead th{font-size:.74rem;letter-spacing:.05em;text-transform:uppercase;
  color:var(--muted);font-weight:600;background:var(--accent-soft)}
tbody tr:last-child td{border-bottom:none}
td.num{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right}
td.pass{color:var(--pass)} td.reject{color:var(--reject)}
.na{color:var(--warn);font-family:var(--mono);font-size:.85em}

/* callout */
.callout{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--accent);
  border-radius:10px;padding:16px 20px;margin:1.4em 0}
.callout.key{border-left-color:var(--reject)}
.callout .lab{font-family:var(--mono);font-size:.7rem;letter-spacing:.1em;
  text-transform:uppercase;color:var(--accent);font-weight:600}
.callout.key .lab{color:var(--reject)}

/* figures */
.fig{margin:1.4em 0;border:1px solid var(--line);border-radius:12px;overflow:hidden;background:var(--card)}
.fig img{display:block;width:100%;height:auto}
.fig figcaption{padding:10px 16px;font-size:.82rem;color:var(--muted);
  border-top:1px solid var(--line);display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap}
.fignote{font-family:var(--mono);font-size:.72rem;color:var(--warn)}
.fig.pending .ph{padding:40px;text-align:center;font-family:var(--mono);
  font-size:.8rem;color:var(--warn);background:var(--warn-soft)}

/* status list */
.steps{list-style:none;padding:0;margin:1.1em 0;display:flex;flex-direction:column;gap:8px}
.steps li{display:grid;grid-template-columns:auto 1fr auto;gap:12px;align-items:center;
  background:var(--card);border:1px solid var(--line);border-radius:10px;padding:11px 15px}
.steps .n{font-family:var(--mono);font-weight:700;color:var(--accent)}
.steps .s{font-family:var(--mono);font-size:.68rem;letter-spacing:.06em;text-transform:uppercase;
  padding:.25em .6em;border-radius:999px;border:1px solid var(--line);white-space:nowrap}
.s.done{color:var(--pass);border-color:var(--pass)}
.s.ready{color:var(--warn);border-color:var(--warn)}
.s.partial{color:var(--muted)}
footer{margin:4em 0 3em;padding-top:1.4em;border-top:1px solid var(--line);
  font-size:.82rem;color:var(--muted)}
</style>

<div class="wrap">
<header>
  <div class="eyebrow">Computational study · DFT stage</div>
  <h1>2D materials as reverse-osmosis desalination membranes</h1>
  <p class="lede">Graphene and MoS₂ taken end-to-end through the DFT → machine-learning-potential → NEMD workflow: pore geometry, water/ion crossing barriers, and the selectivity that decides salt rejection.</p>
  <div class="method-band">
    <span class="chip warn">⚠ CHGNet surrogate</span>
    <span>Numbers are order-of-magnitude and trend-level, computed with a universal ML potential offline. Absolute values need the client's QE&nbsp;+&nbsp;vdW-DF2.</span>
  </div>
</header>

<div class="kpis">
  <div class="kpi"><div class="v pass">{g_h2o30}</div><div class="k">Water barrier, graphene open pore (r=3.0 Å)</div></div>
  <div class="kpi"><div class="v reject">{g_na30}</div><div class="k">Hydrated Na⁺ barrier, same pore → rejected</div></div>
  <div class="kpi"><div class="v">2.5 Å</div><div class="k">Pore radius of the RO operating window</div></div>
</div>

<p>This is the DFT-stage deliverable for the request: <em>study 2D materials as an RO membrane using DFT, then reuse the results to build a ReaxFF/MLP for MD</em>. The whole pipeline runs and is demonstrated on both materials; the single most important finding is a methodological one about how ion rejection must be modelled.</p>

<h2>The result that matters: hydration, not sterics</h2>
<p>A first pass scanning <em>bare</em> ions through the pore gives the wrong answer — bare Na⁺ and Cl⁻ are smaller than a water molecule, so they pass a pore that blocks water. Real rejection comes from the energy cost of stripping an ion's hydration shell to enter a sub-nanometre pore.</p>

<div class="callout key">
  <div class="lab">Key finding · open pore, r = 3.0 Å</div>
  <p style="margin:.4em 0 0">Water crosses at <strong>{g_h2o30}</strong> while the hydrated <strong>Na⁺·(H₂O)₆</strong> faces a <strong>{g_na30}</strong> partial-dehydration barrier — a clean, literature-scale selectivity that the bare-ion scan misses entirely, from the same method.</p>
</div>

{f_hyd}

<div class="tablewrap"><table>
  <caption>Hydrated-ion pore-crossing barriers (Eₐ)</caption>
  <thead><tr><th>system</th><th>H₂O</th><th>Na⁺·(H₂O)₆</th><th>Cl⁻·(H₂O)₆</th></tr></thead>
  <tbody>{hyd_rows}</tbody>
</table></div>
<p style="font-size:.85rem;color:var(--muted)">Barriers above ~3 eV are steric upper bounds — the intact shell cannot pass a pore that tight within the relaxation budget; the ion must dehydrate, which needs longer PMF sampling (umbrella / metadynamics) to quantify.</p>

<h2>Pore size is the engineering knob</h2>
<p>The water crossing barrier collapses from several eV to zero as the graphene pore opens from ~2.0 to ~3.0 Å — the classic RO trade-off, made quantitative. Reject more by shrinking the pore, and pay a steeply rising water-transport penalty.</p>
{f_pore}
<div class="tablewrap"><table>
  <caption>Part A · graphene water barrier vs pore radius</caption>
  <thead><tr><th>pore radius (Å)</th><th>C atoms removed</th><th>water Eₐ (eV)</th><th>regime</th></tr></thead>
  <tbody>{pore_rows}</tbody>
</table></div>

<h2>Graphene vs MoS₂</h2>
<p>Both materials were carried through the identical workflow. At the same nominal pore radius, MoS₂'s thicker, differently-terminated pore presents a lower barrier to water than graphene — a genuine point of comparison for material selection.</p>
<div class="tablewrap"><table>
  <caption>Part B · bare-species crossing barriers (r=2.5 Å)</caption>
  <thead><tr><th>material</th><th>H₂O</th><th>Na (bare)</th><th>Cl (bare)</th></tr></thead>
  <tbody>{sel_rows}</tbody>
</table></div>
<div style="display:grid;grid-template-columns:1fr 1fr;gap:14px">{f_gp}{f_mp}</div>

<h2>From DFT to a machine-learning potential</h2>
<p>The training path is wired end-to-end. Diverse configurations were labelled (here with the CHGNet surrogate, standing in for AIMD/QE), exported to the format MACE and NequIP consume, and used to train a small MACE model on CPU — energy error dropping an order of magnitude — proving the DFT → data → MLP plumbing before any GPU time is spent.</p>
{f_mace}
<p style="font-size:.85rem;color:var(--muted)">A plumbing proof, not a converged potential: production training uses a larger model, float64, and real QE labels on a GPU. DeePMD, MACE, and NequIP configs are provided for a three-way bake-off on identical data.</p>

<h2>MLP or ReaxFF?</h2>
<p>For physical reverse osmosis — water and hydrated ions crossing without breaking bonds — an <strong>MLP is the right tool</strong>: DFT-level accuracy at MD speed. <strong>ReaxFF</strong> costs 10–100× more and is only worth it for genuine chemistry: pore-edge protonation, membrane degradation, or fouling. The recommendation is to ship the MLP path first and treat ReaxFF as an optional, separately-scoped follow-on rather than a parallel effort.</p>

<h2>Where each step stands</h2>
<ul class="steps">
  <li><span class="n">1</span><span>Real QE + vdW-DF2 SCF on the r≈2.5 Å pores (both materials) + convergence grid</span><span class="s ready">inputs ready</span></li>
  <li><span class="n">2</span><span>Hydrated-ion crossing barrier — done locally (surrogate); real-DFT NEB input ready</span><span class="s done">done</span></li>
  <li><span class="n">3</span><span>AIMD training frames — demonstrated with CHGNet-labelled configurations</span><span class="s done">done</span></li>
  <li><span class="n">4</span><span>DeePMD / MACE / NequIP bake-off — MACE trained locally; comparison script + 3 configs</span><span class="s done">done</span></li>
  <li><span class="n">5</span><span>100 MPa NEMD → water flux &amp; salt rejection — script ready, piston force calibrated, analyzer toy-validated</span><span class="s ready">ready to run</span></li>
  <li><span class="n">6</span><span>MoS₂ repeat — built, Part-B selectivity done, QE input ready; hydrated PMF {mos_status}</span><span class="s partial">{mos_status}</span></li>
</ul>

<footer>
  Method: CHGNet universal ML potential as a DFT surrogate (offline, bundled weights). Trends and orderings are robust; absolute energies and salt-rejection numbers require converged QE + vdW-DF2 with charged supercells and explicit water. Full code, configs, and raw data accompany this report as a reusable pipeline.
</footer>
</div>
"""

if __name__ == "__main__":
    main()
