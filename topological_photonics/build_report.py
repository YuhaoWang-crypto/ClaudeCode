"""
Build a single self-contained HTML report (all figures embedded as data URIs).
Output: ../topological_photonics/Topological_Photonics_Report.html
Open in any browser; print-to-PDF for a PDF copy.
"""
import base64
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
OUT = os.path.join(HERE, "Topological_Photonics_Report.html")


def data_uri(fname):
    path = os.path.join(FIG, fname)
    mime = "image/gif" if fname.endswith(".gif") else "image/png"
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return f"data:{mime};base64,{b64}"


def fig(fname, num, title, caption, badge=None):
    b = f'<span class="badge {badge[0]}">{badge[1]}</span>' if badge else ""
    return f"""
    <figure>
      <img src="{data_uri(fname)}" alt="{title}">
      <figcaption><b>Figure {num}. {title}.</b> {caption} {b}</figcaption>
    </figure>"""


RIG = ("rig", "✅ rigorous / ab-initio")
EFF = ("eff", "⚠️ effective model")
ILL = ("eff", "⚠️ illustrative constants")

BODY = f"""
<h1>Topological Photonics &amp; Phononics on a Silica/Silicon Platform</h1>
<p class="sub">Computational study for the proposal <i>&ldquo;Fabrication and Study of
Topological Insulators with Silica-Coated Bi<sub>2</sub>Se<sub>3</sub> for Ultralow-loss
Phononic Waveguides and Photonic Chips&rdquo;</i> &mdash; geometry-induced topology on the
Wu&ndash;Hu C<sub>6</sub> lattice, computed end-to-end in open-source NumPy/SciPy.</p>

<div class="toc">
<b>Contents</b>
<ol>
  <li><a href="#scope">Scientific scope &amp; honest caveats</a></li>
  <li><a href="#system">The model system</a></li>
  <li><a href="#basic">Basic package: bands, invariant, edges, scan</a></li>
  <li><a href="#complete">Complete package: real-rod edges, FDTD, phononic twin</a></li>
  <li><a href="#rigor">Disorder invariant &amp; loss budget</a></li>
  <li><a href="#ledger">Honesty ledger</a></li>
  <li><a href="#run">How to reproduce</a></li>
  <li><a href="#refs">References</a></li>
</ol>
</div>

<h2 id="scope">1. Scientific scope &amp; honest caveats</h2>
<p>These bound what the study asserts and correct several framings in the original
proposal.</p>
<ol class="scope">
  <li><b>Silica coating does not create topology.</b> SiO<sub>2</sub> is a trivial
  wide-gap dielectric. Ref [1] (Bi<sub>2</sub>Se<sub>3</sub>@SiO<sub>2</sub>) shows
  colloidal stability, biocompatibility and <i>retained optical properties</i> &mdash;
  not a topological band structure. Here topology is <b>geometry-induced</b>;
  SiO<sub>2</sub>/Si serve only as low-loss host media.</li>
  <li><b>The cited chip [2] is a different platform</b> (programmable Si
  waveguide/microring/MZI near 1525&nbsp;nm) &mdash; a methodology reference, not a
  Bi<sub>2</sub>Se<sub>3</sub>@SiO<sub>2</sub> material model.</li>
  <li><b>&ldquo;silicon&rdquo; &ne; &ldquo;silica&rdquo;.</b> The six-petal holey lattice [4]
  is <i>silicon</i> (a stiff elastic solid). Elastic modelling must use the correct
  density/elastic constants and fabrication; the phononic module here is a
  <b>scalar-acoustic</b> analogue, not the full-vector elastic plate of [4].</li>
  <li><b>A Dirac gap is not proof of non-trivial topology.</b> We therefore also provide
  the C<sub>6</sub> band-inversion indicator, a topological invariant, interface states
  between two <i>distinct</i> phases, edge-state dispersion with field localisation, and
  direct transport.</li>
  <li><b>The invariant is symmetry-dependent &mdash; not hard-coded &ldquo;Chern&rdquo;.</b>
  These crystals are time-reversal symmetric, so the global Chern number is
  <b>0</b> (verified). We report the <b>spin-Chern / spin-Bott</b> index.</li>
  <li><b>Under disorder, k-space Berry curvature is invalid.</b> We use the
  <b>real-space spin Bott index</b> instead.</li>
  <li><b>&ldquo;Ultralow-loss&rdquo; is not implied by topology.</b> Topology removes
  <i>backscattering</i>; absorption/radiation/interface losses set the floor, which we
  quantify from the real mode.</li>
</ol>

<h2 id="system">2. The model system</h2>
<p>A triangular Bravais lattice carries a hexagonal cluster of six dielectric rods per
cell (the Wu&ndash;Hu crystal, Wu &amp; Hu, <i>PRL</i> <b>114</b>, 223901, 2015). The single
geometric knob is the cluster radius <i>R/a</i>. At <i>R&nbsp;=&nbsp;a/3</i> the rods complete a
honeycomb and the two graphene Dirac cones fold onto &Gamma; into a
<b>double Dirac cone</b>. Shrinking (<i>R&lt;a/3</i>) opens a <b>trivial</b> gap; expanding
(<i>R&gt;a/3</i>) opens a <b>band-inverted (topological)</b> gap in which the p-like
(l=&plusmn;1) and d-like (l=&plusmn;2) doublets swap order &mdash; the photonic analogue of the
quantum spin Hall effect, with the pseudospin supplied by C<sub>6</sub> symmetry.</p>

<h2 id="basic">3. Basic package</h2>
{fig("fig1_bands_dirac_gap.png", 1, "Photonic band structure (plane-wave expansion)",
     "TM bands of the Si-rod crystal along &Gamma;&ndash;M&ndash;K&ndash;&Gamma;. The critical radius shows the "
     "four-fold double Dirac cone at &Gamma;; shrinking and expanding open a gap. The C6 "
     "angular-momentum of the modes shows the p/d doublets invert between the two phases.", RIG)}
{fig("fig2_berry_spin_chern.png", 2, "Berry curvature &amp; spin-Chern number",
     "Fukui&ndash;Hatsugai&ndash;Suzuki lattice computation on the BHZ blocks calibrated from the "
     "PWE gap: C&#8593;=&minus;1, C&#8595;=+1 &rarr; total Chern 0 (time-reversal symmetric) but "
     "spin-Chern &plusmn;1; the invariant jumps at the gap-closing.", EFF)}
{fig("fig3_edge_states_robustness.png", 3, "Helical edge states &amp; disorder (effective model)",
     "BHZ ribbon: a helical Kramers pair crosses the gap only in the topological ribbon, "
     "and the mid-gap crossing survives on-site disorder (20/20 realisations).", EFF)}
{fig("fig4_param_scan_phase.png", 4, "Geometry scan / phase diagram",
     "&Gamma; gap vs cluster radius from the PWE: the gap closes at R=a/3 (Dirac) and reopens "
     "as a topological gap for R&gt;a/3, locating the phase boundary and the operating point "
     "of largest protected bandwidth.", RIG)}

<h2 id="complete">4. Complete package</h2>
{fig("fig5_fdfd_real_edge.png", 5, "Real-rod edge states (FDFD supercell)",
     "Ab-initio finite-difference frequency-domain supercell of the actual Si-rod crystal "
     "with a topological|trivial domain wall. Left: projected bands coloured by wall "
     "localisation &mdash; edge branches traverse the bulk gap. Right: the edge mode hugs the "
     "domain wall. This removes the effective-model caveat of Figure 3.", RIG)}
{fig("fig6_fdtd_snapshots.png", 6, "FDTD wave propagation",
     "Pure-NumPy 2D TM FDTD of the real crystal: an edge wave guides along a straight wall, "
     "turns a sharp double-bend, and passes a point defect (green &times;) &mdash; backscattering-"
     "immune transport.", RIG)}
{fig("fdtd_bend.gif", 7, "FDTD animation of the sharp bend",
     "Time evolution of the edge wave negotiating the two 90&deg; corners (animated GIF).", RIG)}
{fig("fig7_energy_delivery.png", 8, "Backscattering-immunity metric",
     "Steady-state energy delivered to the channel output at f=0.46 c/a: a point defect "
     "barely dents the topological channel (0.61 vs 0.70 straight) while a trivial guide "
     "with no domain wall delivers almost nothing (0.08).", RIG)}
{fig("fig8_photonic_vs_phononic.png", 9, "Two-system comparison: photonic vs phononic",
     "The identical Wu&ndash;Hu lattice solved with two wave kernels (Maxwell TM light; scalar "
     "acoustic sound). Both host a double Dirac cone at R=a/3 and a C6-inverted topological "
     "gap; the phononic gap sits on higher bands at higher frequency, and its topological "
     "side is R&lt;a/3 (opposite to photonic) for the chosen material contrast.", RIG)}

<h2 id="rigor">5. Disorder invariant &amp; loss budget</h2>
{fig("fig9_spin_bott_disorder.png", 10, "Real-space spin Bott index under disorder",
     "Because disorder breaks periodicity, k-space Berry curvature no longer applies. The "
     "Loring&ndash;Hastings spin Bott index (L=14 torus, no Brillouin zone) stays quantised then "
     "collapses at strong disorder (clean topological, blue), and jumps 0&rarr;1 at intermediate "
     "disorder (clean trivial, red) &mdash; a Topological Anderson Insulator, the mechanism of "
     "ref [3].", RIG)}
{fig("fig10_loss_budget.png", 11, "Loss budget: absorption sets the floor, not topology",
     "From the real FDFD edge mode (confinement factor &Gamma;<sub>rod</sub>=0.93) via the "
     "confinement-factor modal-loss formula: a low-loss dielectric gives L<sub>p</sub>&asymp;1.8&nbsp;mm, "
     "but bare Bi<sub>2</sub>Se<sub>3</sub> (&epsilon;&Prime;&asymp;8) gives L<sub>p</sub>&asymp;0.2&nbsp;&mu;m. "
     "&ldquo;Ultralow-loss&rdquo; and &ldquo;Bi<sub>2</sub>Se<sub>3</sub> carries the mode&rdquo; are in "
     "direct tension. Optical constants are representative.", ILL)}

<h2 id="ledger">6. Honesty ledger</h2>
<table>
<tr><th>Claim</th><th>Basis</th><th>Label</th></tr>
<tr><td>Double Dirac cone at R=a/3</td><td>PWE eigenvalues of the real crystal</td><td class="rigc">✅ rigorous</td></tr>
<tr><td>Gap opens on symmetry breaking</td><td>PWE</td><td class="rigc">✅ rigorous</td></tr>
<tr><td>Band inversion (p&harr;d)</td><td>C6 indicator on PWE modes</td><td class="rigc">✅ rigorous</td></tr>
<tr><td>Gap-vs-R phase boundary</td><td>PWE scan</td><td class="rigc">✅ rigorous</td></tr>
<tr><td>Spin-Chern = &plusmn;1 (Chern total 0)</td><td>BHZ + FHS (sign from PWE)</td><td class="effc">⚠️ effective model</td></tr>
<tr><td>Helical edge states in the gap</td><td>FDFD supercell (real rods)</td><td class="rigc">✅ ab-initio</td></tr>
<tr><td>Guiding / bend / defect bypass</td><td>FDTD</td><td class="rigc">✅ (qualitative transport)</td></tr>
<tr><td>Phononic double Dirac + gap + inversion</td><td>acoustic PWE</td><td class="rigc">✅ ab-initio</td></tr>
<tr><td>Disorder invariant + TAI (no BZ)</td><td>spin Bott index</td><td class="rigc">✅ rigorous</td></tr>
<tr><td>Absorption-limited L<sub>p</sub> / dB-cm</td><td>&Gamma; (FDFD) + perturbation</td><td class="effc">✅ method; ⚠️ illustrative &epsilon;&Prime;</td></tr>
</table>

<h2 id="run">7. How to reproduce</h2>
<pre>pip install -r requirements.txt
cd src
python run_all.py          # basic package  (~90 s)  -> fig1..fig4
python run_all.py --full   # + complete pkg (~4 min) -> fig5..fig10 + GIF</pre>
<p>All results are deterministic (fixed disorder seeds). No proprietary software
(no Lumerical/COMSOL/MEEP); pure NumPy/SciPy/Matplotlib.</p>
<p><b>Remaining upgrades</b> (documented in <code>REPORT.md</code> &sect;5): full-vector
elastodynamics for the true ref-[4] system; UPML + eigenmode-source FDTD for calibrated
dB transmission; spin-Chern from real modes via Wilson loops; local Chern marker maps;
and measured Bi<sub>2</sub>Se<sub>3</sub>/SiO<sub>2</sub> optical constants for the loss budget.</p>

<h2 id="refs">8. References</h2>
<ol class="refs">
<li>Belec B. <i>et&nbsp;al.</i> Silica-coated Bi<sub>2</sub>Se<sub>3</sub> topological-insulator
nanoparticles. <i>Nanomaterials</i> <b>13</b>, 809 (2023).</li>
<li>Dai T. <i>et&nbsp;al.</i> A programmable topological photonic chip. <i>Nat. Mater.</i>
<b>23</b>, 928 (2024).</li>
<li>Li Z. <i>et&nbsp;al.</i> Time-reversal-symmetric 2D photonic topological Anderson
insulator. arXiv:2501.11251.</li>
<li>Xu X.-B. <i>et&nbsp;al.</i> Gigahertz topological phononic circuits on unsuspended
Si waveguide arrays. <i>Nat. Electron.</i> <b>8</b>, 689 (2025).</li>
<li>Wu L.-H. &amp; Hu X. Scheme for achieving a topological photonic crystal by using
dielectric material. <i>PRL</i> <b>114</b>, 223901 (2015).</li>
<li>Loring T. A. &amp; Hastings M. B. Disordered topological insulators via C*-algebras
(Bott index). <i>EPL</i> <b>92</b>, 67004 (2010).</li>
</ol>
<p class="foot">Generated from the <code>topological_photonics</code> package. Figures are
computed, not drawn. Print this page to PDF for a portable copy.</p>
"""

CSS = """
:root{--bg:#ffffff;--fg:#1a1a2e;--mut:#556;--card:#f6f7fb;--bd:#dde;--acc:#1b3b6f;--rig:#1a7a3a;--eff:#a86400}
@media(prefers-color-scheme:dark){:root{--bg:#12131a;--fg:#e9ebf3;--mut:#9aa;--card:#1b1d28;--bd:#333a4d;--acc:#7fa8e0;--rig:#5fce8a;--eff:#e0b45f}}
:root[data-theme=dark]{--bg:#12131a;--fg:#e9ebf3;--mut:#9aa;--card:#1b1d28;--bd:#333a4d;--acc:#7fa8e0;--rig:#5fce8a;--eff:#e0b45f}
:root[data-theme=light]{--bg:#fff;--fg:#1a1a2e;--mut:#556;--card:#f6f7fb;--bd:#dde;--acc:#1b3b6f;--rig:#1a7a3a;--eff:#a86400}
*{box-sizing:border-box}body{background:var(--bg);color:var(--fg);font:16px/1.6 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;padding:2.2rem 1rem}
.wrap{max-width:900px;margin:0 auto}
h1{font-size:1.85rem;line-height:1.25;color:var(--acc);margin:.2em 0 .3em}
h2{font-size:1.3rem;margin-top:2.2em;padding-bottom:.25em;border-bottom:2px solid var(--bd)}
.sub{color:var(--mut);font-size:1.02rem}
.toc{background:var(--card);border:1px solid var(--bd);border-radius:10px;padding:.6rem 1.2rem;margin:1.6rem 0}
.toc ol{margin:.3em 0}.toc a{color:var(--acc);text-decoration:none}.toc a:hover{text-decoration:underline}
figure{margin:1.6rem 0;text-align:center}
img{max-width:100%;height:auto;border:1px solid var(--bd);border-radius:8px;background:#fff}
figcaption{color:var(--mut);font-size:.9rem;margin-top:.5rem;text-align:left}
.badge{display:inline-block;font-size:.75rem;padding:.1em .5em;border-radius:6px;font-weight:600;white-space:nowrap}
.badge.rig{background:rgba(26,122,58,.15);color:var(--rig)}.badge.eff{background:rgba(168,100,0,.15);color:var(--eff)}
table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.92rem;overflow-x:auto;display:block}
th,td{border:1px solid var(--bd);padding:.45em .6em;text-align:left}th{background:var(--card)}
.rigc{color:var(--rig);font-weight:600}.effc{color:var(--eff);font-weight:600}
ol.scope>li{margin:.5em 0}
pre{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:.9rem 1rem;overflow-x:auto;font-size:.9rem}
code{background:var(--card);padding:.1em .35em;border-radius:4px}
.refs li,.foot{color:var(--mut)}.foot{margin-top:2.5em;font-size:.85rem}
a{color:var(--acc)}
"""

HTML = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Topological Photonics &amp; Phononics — Full Report</title>
<style>{CSS}</style></head>
<body><div class="wrap">{BODY}</div></body></html>"""

with open(OUT, "w") as f:
    f.write(HTML)
print("wrote", OUT, f"({os.path.getsize(OUT)//1024} KB)")
