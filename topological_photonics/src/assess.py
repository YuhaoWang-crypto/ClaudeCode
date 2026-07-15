"""
Pre-fabrication assessment engine — a SERVICE-DELIVERABLE wrapper.

Input: a geometry/material spec (a dict). Output: a structured go / conditional /
no-go assessment with computed numbers for the four decision axes a fab customer
cares about:

    1. TOPOLOGY   — is the gap actually band-inverted? which invariant?
    2. BANDWIDTH  — gap size, centre frequency, fractional bandwidth, λ mapping.
    3. YIELD      — fabrication tolerance: geometric ΔR that halves the gap, and
                    the disorder strength W_c at which topology is lost.
    4. LOSS       — propagation length / dB·cm⁻¹ from the real mode + material ε''.

Each axis gets a PASS / WATCH / FAIL flag and the overall verdict is the worst of
them. `assess(spec)` returns the dict; `render_html(result)` writes a one-page
assessment card. This turns the physics pipeline into a repeatable report.
"""
from __future__ import annotations
import os
import numpy as np

import geometry as geo
import symmetry as sym
from bands import gamma_gap
import loss_budget

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")

DEFAULT_SPEC = {
    "name": "Wu-Hu photonic TI (Si rods, telecom)",
    "system": "photonic",
    "R_over_a": 0.36,
    "eps_imag": 1e-3,       # rod absorption (ε''); Si/SiO2 ~ 1e-3, Bi2Se3 ~ 8
    "lambda_um": 1.55,
    "device_length_um": 100.0,   # required functional length
    "min_fractional_bw": 0.05,   # bandwidth spec
}


def _flag(ok, watch):
    return "PASS" if ok else ("WATCH" if watch else "FAIL")


def assess(spec=None):
    s = {**DEFAULT_SPEC, **(spec or {})}
    R = s["R_over_a"]
    pc = geo.wu_hu_crystal(R, nmax=7)

    # --- 1. topology ---
    gap = gamma_gap(pc)
    diag = sym.diagnose(pc) if gap > 0.008 else {"phase": "gapless", "inverted": False}
    topo = diag["phase"] == "topological"

    # --- 2. bandwidth ---
    f = pc.solve_k([0, 0], 6)
    fc = 0.5 * (f[2] + f[3])
    frac_bw = gap / fc if fc > 0 else 0.0
    a_um = s["lambda_um"] * fc          # a so that lambda = a/fc
    bw_ok = frac_bw >= s["min_fractional_bw"]

    # --- 3. yield: geometric tolerance (ΔR halving the gap) + disorder W_c ---
    dR = 0.0
    for d in np.linspace(0, 0.06, 25):
        g2 = gamma_gap(geo.wu_hu_crystal(R + d, nmax=7))
        g3 = gamma_gap(geo.wu_hu_crystal(max(R - d, 0.2), nmax=7))
        if min(g2, g3) < 0.5 * gap:
            dR = d
            break
    else:
        dR = 0.06
    dR_nm = dR * a_um * 1000            # tolerance in nm
    Wc = _disorder_Wc() if topo else 0.0

    # --- 4. loss ---
    gamma = loss_budget.edge_mode_confinement()
    _, db_cm, Lp_um = loss_budget.modal_loss(gamma, s["eps_imag"])
    loss_ok = Lp_um >= s["device_length_um"]
    loss_watch = Lp_um >= 0.3 * s["device_length_um"]

    axes = {
        "Topology": {
            "flag": _flag(topo, diag["phase"] == "gapless"),
            "detail": f"phase = {diag['phase']}; invariant = spin-Chern/valley "
                      f"(global Chern 0, T-symmetric)",
            "value": diag["phase"],
        },
        "Bandwidth": {
            "flag": _flag(bw_ok, frac_bw >= 0.5 * s["min_fractional_bw"]),
            "detail": f"gap Δf={gap:.3f}, f_c={fc:.3f} (ωa/2πc), fractional "
                      f"BW={frac_bw*100:.1f}%  →  a≈{a_um:.2f} µm @ λ={s['lambda_um']} µm",
            "value": f"{frac_bw*100:.1f}%",
        },
        "Yield / tolerance": {
            "flag": _flag(dR_nm > 15 and (Wc > 4 or not topo), dR_nm > 8),
            "detail": f"geometric tolerance ±{dR_nm:.0f} nm (halves gap); "
                      f"disorder W_c≈{Wc:.1f} (spin-Bott collapse)"
                      + ("" if topo else "  [n/a: trivial]"),
            "value": f"±{dR_nm:.0f} nm",
        },
        "Loss": {
            "flag": _flag(loss_ok, loss_watch),
            "detail": f"Γ_rod={gamma:.2f}, ε''={s['eps_imag']}  →  α={db_cm:.0f} dB/cm, "
                      f"L_p={Lp_um:.1f} µm (need ≥{s['device_length_um']:.0f} µm)",
            "value": f"L_p={Lp_um:.0f} µm",
        },
    }
    order = {"FAIL": 0, "WATCH": 1, "PASS": 2}
    worst = min(axes.values(), key=lambda a: order[a["flag"]])["flag"]
    verdict = {"FAIL": "NO-GO", "WATCH": "CONDITIONAL", "PASS": "GO"}[worst]

    return {"spec": s, "axes": axes, "verdict": verdict,
            "summary": _summary(verdict, axes)}


def _disorder_Wc(L=10, M=1.0):
    """Quick estimate of the disorder strength where the spin-Bott index drops."""
    import bott
    rng = np.random.default_rng(0)
    for W in np.linspace(2, 12, 6):
        vals = [bott.spin_bott(L, M, W, rng)[2] for _ in range(3)]
        if np.mean(vals) < 0.5:
            return W
    return 12.0


def _summary(verdict, axes):
    fails = [k for k, v in axes.items() if v["flag"] == "FAIL"]
    watch = [k for k, v in axes.items() if v["flag"] == "WATCH"]
    if verdict == "GO":
        return "All four axes pass — proceed to fabrication with the stated tolerances."
    if verdict == "NO-GO":
        return f"Blocked by: {', '.join(fails)}. Revise before fabrication."
    return f"Conditional — watch: {', '.join(watch)}. Fabricate a de-risking test die first."


CSS = """body{font:15px/1.55 -apple-system,Segoe UI,Roboto,Arial,sans-serif;max-width:820px;
margin:2rem auto;padding:0 1rem;color:#1a1a2e}h1{color:#1b3b6f;font-size:1.4rem}
table{border-collapse:collapse;width:100%;margin:1rem 0}th,td{border:1px solid #dde;padding:.5em .7em;
text-align:left;vertical-align:top}th{background:#f6f7fb}.PASS{color:#1a7a3a;font-weight:700}
.WATCH{color:#a86400;font-weight:700}.FAIL{color:#c0392b;font-weight:700}
.verdict{font-size:1.6rem;font-weight:800;padding:.4rem 1rem;border-radius:10px;display:inline-block}
.GO{background:#e6f6ec;color:#1a7a3a}.CONDITIONAL{background:#fdf3e0;color:#a86400}
.NOGO{background:#fbeaea;color:#c0392b}.muted{color:#667;font-size:.9rem}"""


def render_html(result, path=None):
    s = result["spec"]
    rows = "".join(
        f"<tr><td><b>{k}</b></td><td class='{v['flag']}'>{v['flag']}</td>"
        f"<td>{v['detail']}</td></tr>"
        for k, v in result["axes"].items())
    vclass = result["verdict"].replace("-", "")
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Pre-fab assessment — {s['name']}</title><style>{CSS}</style></head><body>
<h1>Pre-fabrication assessment</h1>
<p class="muted">Device: <b>{s['name']}</b> · system: {s['system']} · R/a={s['R_over_a']} ·
ε''={s['eps_imag']} · λ={s['lambda_um']} µm · required length {s['device_length_um']} µm</p>
<p>Verdict: <span class="verdict {vclass}">{result['verdict']}</span></p>
<p>{result['summary']}</p>
<table><tr><th>Axis</th><th>Flag</th><th>Computed basis</th></tr>{rows}</table>
<p class="muted">Generated by the topological-photonics <code>assess.py</code> engine.
Numbers are computed (PWE bands, C6 indicator, spin-Bott disorder, mode-based loss).
Absorption uses representative ε''; drop in measured constants for a final number.</p>
</body></html>"""
    path = path or os.path.join(FIGDIR, "assessment_card.html")
    with open(path, "w") as f:
        f.write(html)
    return path


if __name__ == "__main__":
    for name, spec in [
        ("A: Si rods, telecom (good)", {}),
        ("B: Bi2Se3 core rods (lossy)", {"name": "Wu-Hu with Bi2Se3 rods",
                                         "eps_imag": 8.0}),
        ("C: trivial geometry", {"name": "Shrunken (trivial)", "R_over_a": 0.30}),
    ]:
        r = assess(spec)
        print(f"\n=== {name} ===  VERDICT: {r['verdict']}")
        for k, v in r["axes"].items():
            print(f"  [{v['flag']:5s}] {k}: {v['detail']}")
    render_html(assess({}))
    print("\nwrote assessment_card.html")
