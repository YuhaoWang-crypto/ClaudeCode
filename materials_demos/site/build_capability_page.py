"""Build a self-contained capability-showcase HTML page with the real demo
figures embedded as data URIs. Downscales each figure for the web, base64-
encodes it, and injects it into the HTML template."""
import base64
import io
import os
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.dirname(HERE)
OUT = os.path.join(HERE, "materials_capability.html")

FIGS = {
    "FIG_QD": "qd_fss_demo.png",
    "FIG_METAL": "metal_center_screen.png",
    "FIG_MS": "moltensalt_demo.png",
    "FIG_CD": "cd_binding_sites.png",
    "FIG_SERS": "sers_chirality_demo.png",
}


def data_uri(path, maxw=1000, quality=82):
    im = Image.open(path).convert("RGB")
    if im.width > maxw:
        h = round(im.height * maxw / im.width)
        im = im.resize((maxw, h), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, format="JPEG", quality=quality, optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode()
    return f"data:image/jpeg;base64,{b64}", len(b64)


def main():
    html = TEMPLATE
    total = 0
    for token, fname in FIGS.items():
        uri, n = data_uri(os.path.join(DEMO, fname))
        total += n
        html = html.replace("{{" + token + "}}", uri)
    with open(OUT, "w") as fh:
        fh.write(html)
    print(f"wrote {OUT}  ({os.path.getsize(OUT)/1024:.0f} KB, images {total/1024:.0f} KB b64)")


TEMPLATE = r"""<style>
  :root{
    --ink:#0e1119; --panel:#151a25; --line:#e2e6ec;
    --paper:#eceff4; --card:#ffffff;
    --fg:#161b24; --muted:#5a6474; --faint:#8a94a4;
    --signal:#2f74ff; --signal-soft:#dbe6ff;
    --local:#1f9d57; --local-soft:#dcf2e6;
    --hpc:#c9791a; --hpc-soft:#f8ebd6;
    --sans:system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",sans-serif;
    --mono:ui-monospace,"SF Mono","JetBrains Mono","Cascadia Code",Menlo,monospace;
    --maxw:1080px;
  }
  @media (prefers-color-scheme:dark){
    :root{
      --paper:#0e1119; --card:#161c28; --line:#252c3a;
      --fg:#e7ebf2; --muted:#9aa5b6; --faint:#6b7688;
      --signal:#5b93ff; --signal-soft:#1b2740;
      --local:#43c47e; --local-soft:#123024;
      --hpc:#e2a04a; --hpc-soft:#332413;
    }
  }
  :root[data-theme="dark"]{
    --paper:#0e1119; --card:#161c28; --line:#252c3a;
    --fg:#e7ebf2; --muted:#9aa5b6; --faint:#6b7688;
    --signal:#5b93ff; --signal-soft:#1b2740;
    --local:#43c47e; --local-soft:#123024;
    --hpc:#e2a04a; --hpc-soft:#332413;
  }
  :root[data-theme="light"]{
    --paper:#eceff4; --card:#ffffff; --line:#e2e6ec;
    --fg:#161b24; --muted:#5a6474; --faint:#8a94a4;
    --signal:#2f74ff; --signal-soft:#dbe6ff;
    --local:#1f9d57; --local-soft:#dcf2e6;
    --hpc:#c9791a; --hpc-soft:#f8ebd6;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--paper);color:var(--fg);font-family:var(--sans);
    line-height:1.6;-webkit-font-smoothing:antialiased}
  .wrap{max-width:var(--maxw);margin:0 auto;padding:0 24px}
  a{color:var(--signal)}
  .eyebrow{font-family:var(--mono);font-size:12px;letter-spacing:.18em;
    text-transform:uppercase;color:var(--muted)}

  /* hero */
  header{border-bottom:1px solid var(--line);padding:64px 0 52px;position:relative;overflow:hidden}
  .prompt{font-family:var(--mono);font-size:13.5px;color:var(--muted);
    border:1px solid var(--line);border-radius:10px;padding:14px 16px;background:var(--card);
    max-width:640px}
  .prompt .c{color:var(--signal)}
  h1{font-size:clamp(30px,5vw,50px);line-height:1.05;letter-spacing:-.02em;
    text-wrap:balance;margin:22px 0 14px;font-weight:680}
  h1 .hl{color:var(--signal)}
  .lede{font-size:18px;color:var(--muted);max-width:60ch;margin:0}
  .metrics{display:flex;flex-wrap:wrap;gap:28px;margin-top:34px}
  .metric .n{font-family:var(--mono);font-size:26px;font-weight:600;color:var(--fg)}
  .metric .l{font-family:var(--mono);font-size:11.5px;letter-spacing:.1em;
    text-transform:uppercase;color:var(--faint)}

  /* two-tier band */
  .tiers{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:44px 0}
  @media(max-width:720px){.tiers{grid-template-columns:1fr}}
  .tier{border:1px solid var(--line);border-radius:14px;padding:22px;background:var(--card)}
  .tier h3{margin:10px 0 8px;font-size:18px}
  .tier p{margin:0;color:var(--muted);font-size:14.5px}
  .tier.local{border-top:3px solid var(--local)}
  .tier.hpc{border-top:3px solid var(--hpc)}
  .badge{font-family:var(--mono);font-size:11px;letter-spacing:.08em;text-transform:uppercase;
    padding:3px 9px;border-radius:999px;font-weight:600;display:inline-block}
  .badge.local{background:var(--local-soft);color:var(--local)}
  .badge.hpc{background:var(--hpc-soft);color:var(--hpc)}

  section{padding:46px 0;border-top:1px solid var(--line)}
  .sec-head{margin-bottom:28px}
  .sec-head h2{font-size:24px;margin:8px 0 6px;letter-spacing:-.01em}
  .sec-head p{margin:0;color:var(--muted);max-width:64ch}

  /* demo cards */
  .grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}
  @media(max-width:760px){.grid{grid-template-columns:1fr}}
  .card{border:1px solid var(--line);border-radius:16px;background:var(--card);
    overflow:hidden;display:flex;flex-direction:column}
  .card .fig{background:#fff;border-bottom:1px solid var(--line);aspect-ratio:16/7;
    overflow:hidden;display:flex;align-items:center;justify-content:center}
  .card .fig img{width:100%;height:100%;object-fit:cover;object-position:center}
  .card .body{padding:18px 20px 20px;display:flex;flex-direction:column;gap:10px;flex:1}
  .card .num{font-family:var(--mono);font-size:12px;color:var(--faint);letter-spacing:.12em}
  .card h3{margin:0;font-size:18px;letter-spacing:-.01em}
  .ask{font-size:13.5px;color:var(--muted);border-left:2px solid var(--signal);
    padding:2px 0 2px 12px;font-style:italic}
  .rows{display:flex;flex-direction:column;gap:7px;margin-top:2px}
  .row{display:flex;gap:9px;align-items:flex-start;font-size:13.5px}
  .row .k{font-family:var(--mono);font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;
    padding:2px 7px;border-radius:6px;white-space:nowrap;margin-top:1px;font-weight:600}
  .row .k.here{background:var(--local-soft);color:var(--local)}
  .row .k.hpc{background:var(--hpc-soft);color:var(--hpc)}
  .row .v{color:var(--fg)}
  .row .v .file{font-family:var(--mono);font-size:12px;color:var(--muted)}

  /* matrix */
  .tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:14px}
  table{border-collapse:collapse;width:100%;font-size:14px;min-width:640px}
  th,td{text-align:left;padding:12px 16px;border-bottom:1px solid var(--line)}
  thead th{font-family:var(--mono);font-size:11px;letter-spacing:.1em;text-transform:uppercase;
    color:var(--muted);font-weight:600;background:var(--card)}
  tbody tr:last-child td{border-bottom:none}
  td .mono{font-family:var(--mono);font-size:12.5px;color:var(--muted)}
  .chk{color:var(--local);font-weight:700}
  .warn{color:var(--hpc);font-weight:700}

  /* honesty */
  .honesty{background:var(--card);border:1px solid var(--line);border-radius:16px;padding:26px 28px}
  .honesty h2{margin:0 0 14px;font-size:20px}
  .honesty ul{margin:0;padding-left:20px;color:var(--muted)}
  .honesty li{margin:8px 0}
  .honesty strong{color:var(--fg)}

  footer{padding:40px 0 64px;color:var(--faint);font-family:var(--mono);font-size:12px}

  .reveal{opacity:0;transform:translateY(14px);transition:opacity .6s ease,transform .6s ease}
  .reveal.in{opacity:1;transform:none}
  @media(prefers-reduced-motion:reduce){.reveal{opacity:1;transform:none;transition:none}}
</style>

<header>
  <div class="wrap">
    <div class="prompt">&gt; <span class="c">client:</span> "can you compute this for my material?"<br>&gt; <span class="c">reply:</span> a runnable demo today &nbsp;+&nbsp; the HPC input files for production</div>
    <h1>Theoretical calculation, <span class="hl">de-risked before the cluster runs.</span></h1>
    <p class="lede">We turn a computational-materials request into two things: a small
      calculation that runs immediately and produces real numbers &amp; a figure, and a
      validated scaffold — workflow plus HPC input files — for the full periodic job.
      Six worked examples below, each drawn from a real request.</p>
    <div class="metrics">
      <div class="metric"><div class="n">6</div><div class="l">worked demos</div></div>
      <div class="metric"><div class="n">DFT · MD · docking</div><div class="l">methods covered</div></div>
      <div class="metric"><div class="n">VASP·QE·CP2K·DeePMD</div><div class="l">production targets scaffolded</div></div>
    </div>
  </div>
</header>

<div class="wrap">
  <div class="tiers">
    <div class="tier local">
      <span class="badge local">runs immediately</span>
      <h3>Local demo / proxy</h3>
      <p>Analytical &amp; tight-binding models, small-molecule DFT (PySCF), classical MD,
        structure builders, geometric predictors, training-set packaging. Real output,
        validated against something known — proof the workflow and physics hold.</p>
    </div>
    <div class="tier hpc">
      <span class="badge hpc">needs HPC</span>
      <h3>Production run</h3>
      <p>Periodic DFT (VASP / Quantum ESPRESSO / CP2K), AIMD at scale, electron-phonon
        (EPW), GW-BSE, phonons, reaction barriers, ML-potential training. We hand over
        submittable input files; the numbers come from the client's cluster.</p>
    </div>
  </div>
</div>

<section>
  <div class="wrap">
    <div class="sec-head">
      <div class="eyebrow">Portfolio</div>
      <h2>Six calculations, six real requests</h2>
      <p>Each card: what the client asked, what we ran locally (with the figure it
        produced), and what the production version needs. Proxies are labeled as proxies.</p>
    </div>

    <div class="grid">

      <article class="card reveal">
        <div class="fig"><img src="{{FIG_QD}}" alt="Quantum-dot fine-structure-splitting figure"></div>
        <div class="body">
          <div class="num">DEMO 01 · QUANTUM DOTS</div>
          <h3>Exciton fine-structure splitting under stress</h3>
          <div class="ask">"Can you do this theoretical calculation — Figs 2 &amp; 3?" (exciton polarization / FSS vs. uniaxial stress)</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">Analytical two-level exciton model in numpy — reproduces the paper's limits exactly. <span class="file">qd_fss_model.py</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">Atomistic empirical-pseudopotential supercell for a specific dot.</span></div>
          </div>
        </div>
      </article>

      <article class="card reveal">
        <div class="fig"><img src="{{FIG_METAL}}" alt="Metal-center DFT screening figure"></div>
        <div class="body">
          <div class="num">DEMO 02 · SINGLE-ATOM CATALYSIS</div>
          <h3>Metal-center electronic pre-screen</h3>
          <div class="ask">"Screen metal centers by stability / electronic structure / charge for a biomimetic (peroxidase-like) catalyst."</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">PySCF UKS DFT on an [M–OH] site across Cr/Mn/Fe/Co/Ni/Cu: gap, charge, spin. <span class="file">metal_center_screen.py</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">Periodic Cr-Co₃O₄ slab + NEB reaction barriers.</span></div>
          </div>
        </div>
      </article>

      <article class="card reveal">
        <div class="fig"><img src="{{FIG_MS}}" alt="Molten-salt MD and DeePMD figure"></div>
        <div class="body">
          <div class="num">DEMO 03 · MOLTEN-SALT ML POTENTIAL</div>
          <h3>MD → properties + DeePMD packaging</h3>
          <div class="ask">"AIMD at 6 temperatures → a DeePMD training set → density &amp; viscosity."</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">6-temperature MD → g(r), diffusion Arrhenius; frames packaged in DeePMD/npy (client's format, validated) + real FLiNaK/LiCl-KCl builders. <span class="file">moltensalt_mlp_pipeline.py</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">CP2K/VASP AIMD → dp train → LAMMPS NPT/Green-Kubo.</span></div>
          </div>
        </div>
      </article>

      <article class="card reveal">
        <div class="fig"><img src="{{FIG_CD}}" alt="Cd binding-site prediction figure"></div>
        <div class="body">
          <div class="num">DEMO 04 · METALLOPROTEIN</div>
          <h3>Cd²⁺ binding-site prediction</h3>
          <div class="ask">"Find which sites bind Cd²⁺, so I can design mutants for MST."</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">HSAB donor scoring + clique clustering → ranked sites &amp; mutation targets. Validated: top site on 1CA2 = the real His94/96/119 metal site. <span class="file">cd_binding_site_predictor.py</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">Metal-aware docking + QM/MM coordination refinement.</span></div>
          </div>
        </div>
      </article>

      <article class="card reveal">
        <div class="fig"><img src="{{FIG_SERS}}" alt="SERS chirality DFT figure"></div>
        <div class="body">
          <div class="num">DEMO 05 · CHIRAL SERS</div>
          <h3>Why L/D differ on gold</h3>
          <div class="ask">"Compare L- vs D-arginine on gold: adsorption energy &amp; Raman / SERS difference."</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">DFT proves free-molecule L≡D to ~0.002 cm⁻¹ → the SERS difference must live in the surface adsorption geometry. <span class="file">sers_chirality_dft.py</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">Au(111) slab + D3, both enantiomers' poses, surface Raman.</span></div>
          </div>
        </div>
      </article>

      <article class="card reveal">
        <div class="fig" style="aspect-ratio:16/7;background:linear-gradient(135deg,var(--signal-soft),var(--hpc-soft));display:flex;align-items:center;justify-content:center">
          <div style="font-family:var(--mono);font-size:13px;color:var(--muted);text-align:center;padding:20px;line-height:1.9">
            phonon DOS &nbsp;→&nbsp; S = ½(ΔQ)²ω/ħ<br>ΔSCF STE geometry &nbsp;→&nbsp; CC diagram<br><span style="color:var(--hpc)">reviewer rebuttal, costed</span>
          </div>
        </div>
        <div class="body">
          <div class="num">DEMO 06 · SCINTILLATOR / STE</div>
          <h3>Reviewer-rebuttal calculation plan</h3>
          <div class="ask">"Reviewers: no phonon / electron-phonon work; STE mechanism unsupported; only ground state."</div>
          <div class="rows">
            <div class="row"><span class="k here">ran here</span><span class="v">A costed method map: each reviewer point → method → software → effort. <span class="file">demo3_method_checklist.md</span></span></div>
            <div class="row"><span class="k hpc">production</span><span class="v">DFPT phonons + Huang-Rhys factor + ΔSCF STE geometry (highest ROI first).</span></div>
          </div>
        </div>
      </article>

    </div>
  </div>
</section>

<section>
  <div class="wrap">
    <div class="sec-head">
      <div class="eyebrow">Coverage</div>
      <h2>Capability matrix</h2>
      <p>What runs in a lightweight environment today vs. what we scaffold for the cluster.</p>
    </div>
    <div class="tablewrap">
      <table>
        <thead><tr><th>Direction</th><th>Local demo</th><th>Production target</th><th>Runs here</th></tr></thead>
        <tbody>
          <tr><td>Quantum-dot FSS</td><td class="mono">analytical exciton model</td><td class="mono">atomistic pseudopotential</td><td class="chk">✓</td></tr>
          <tr><td>Catalyst metal screen</td><td class="mono">[M–OH] PySCF DFT</td><td class="mono">Cr-Co₃O₄ slab + NEB</td><td class="chk">✓</td></tr>
          <tr><td>Molten-salt ML potential</td><td class="mono">MD + DeePMD packaging</td><td class="mono">AIMD + dp train + LAMMPS</td><td class="chk">✓</td></tr>
          <tr><td>Protein Cd²⁺ sites</td><td class="mono">HSAB site predictor</td><td class="mono">docking + QM/MM</td><td class="chk">✓</td></tr>
          <tr><td>Chiral SERS on Au</td><td class="mono">molecular DFT frequencies</td><td class="mono">Au(111) slab + Raman</td><td class="chk">✓</td></tr>
          <tr><td>Scintillator / STE</td><td class="mono">costed method map</td><td class="mono">DFPT + Huang-Rhys + ΔSCF</td><td class="warn">plan</td></tr>
        </tbody>
      </table>
    </div>
  </div>
</section>

<section>
  <div class="wrap">
    <div class="honesty">
      <h2>What we don't do</h2>
      <ul>
        <li>A molecular cluster is <strong>not</strong> the periodic material; a classical potential is <strong>not</strong> DFT; a geometric predictor is <strong>not</strong> docking. Every proxy is labeled as one.</li>
        <li>Template defaults — cutoffs, k-meshes, timesteps, functionals — are <strong>reasonable starting points, not converged production settings</strong>. We flag what still needs convergence on your system.</li>
        <li>We report what validated and what didn't. A packaging demo built on cheap forces proves the <strong>plumbing</strong>; the physics still needs the real AIMD.</li>
      </ul>
    </div>
  </div>
</section>

<footer>
  <div class="wrap">computational-materials capability demos · figures generated by the runnable scripts in materials_demos/ · proxy vs. production labeled throughout</div>
</footer>

<script>
  (function(){
    var els=document.querySelectorAll('.reveal');
    if(!('IntersectionObserver' in window)){els.forEach(function(e){e.classList.add('in')});return;}
    var io=new IntersectionObserver(function(en){en.forEach(function(x){if(x.isIntersecting){x.target.classList.add('in');io.unobserve(x.target);}})},{threshold:.14});
    els.forEach(function(e){io.observe(e);});
  })();
</script>
"""

if __name__ == "__main__":
    main()
