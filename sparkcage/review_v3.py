"""
Review of the third Biomni report, whose new content is section M6.

    python3 sparkcage/review_v3.py

WHAT IS NEW. Sections M1 to M5 are unchanged from the version reviewed in
reproduce_biomni.py. M6 adds two things:

  Part 1  the surface-confined reporter and the diffusion field are now solved
          as one coupled system, so a dye that lets go of the electrode and is
          recaptured can be modelled. This is new physics, not a refinement.
  Part 2  the first ABSOLUTE current predictions in the whole series: seven
          candidate sites on helix alpha6, each with a closed-state and an
          open-state peak, and a recommendation of which site to use at which
          frequency.

Part 1 is the more interesting calculation. Part 2 is the one that will be
acted on, because it names sites and frequencies, so it gets checked first and
hardest here.

The short version: the two implementations agree to a fraction of a percent
where they compute the same quantity, which is the strongest cross-validation
in this whole exchange -- and Part 2's site ranking rests on a tunnelling decay
constant that is ten times too small, which inverts its conclusion.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from echem import swv                      # noqa: E402
from structure import analyze              # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"

# The report's own M6 Part 2 table, transcribed from the HTML design report.
SITES = [
    # name, d_closed (nm), k0_closed (1/s), k0_open (1/s), closed pk, open pk, polarity
    ("K250", 1.83, 717, 51.0,   3.8, 525.9, "inverted"),
    ("R254", 2.41, 402, 29.0,  53.4, 371.3, "inverted"),
    ("R257", 2.87, 253, 18.0, 173.2, 261.4, "inverted (weak)"),
    ("R261", 3.49, 137,  9.8, 492.6, 163.7, "normal"),
    ("E264", 3.91,  90,  6.4, 595.5, 122.7, "normal"),
    ("D265", 4.06,  77,  5.5, 594.9, 112.4, "normal"),
    ("R268", 4.53,  48,  3.5, 511.6,  90.3, "normal"),
]

# Their stated model: k0(d) = 1e3 * exp(-beta (d - 1.5 nm)) with beta in 1/nm,
# and an opening displacement of 2.64 nm applied identically to every site.
K0_REF, D_REF_NM, DISPLACEMENT_NM = 1.0e3, 1.5, 2.64
BETA_THEIRS_PER_NM = 1.0
BETA_LITERATURE_PER_NM = 10.0          # = 1.0 per Angstrom

# Their cell, as stated: 3 mm disk, 3e11 molecules/cm2, endpoint sampling.
CELL = dict(gamma_mol_cm2=3.0e11 / 6.02214076e23, area_cm2=0.071,
            n_electrons=2, sample_window=0.02)


def peak_na(k0, freq_hz=100.0, **kw):
    if k0 <= 0.0:
        return 0.0
    return abs(swv.swv_scan(k0=float(k0), freq_hz=freq_hz,
                            **{**CELL, **kw})["peak_current"]) * 1e9


# ---------------------------------------------------------------------------
def section_agreement():
    print("=" * 78)
    print("1. WHERE THE TWO IMPLEMENTATIONS AGREE, AND IT IS WORTH SAYING FIRST")
    print("=" * 78)

    ks = np.logspace(-1, 5, 160)
    p = np.array([peak_na(k) for k in ks])
    qrm = ks[int(np.argmax(p))] / 100.0
    ref = peak_na(100.0)

    print(f"""
  The report states two numbers that this package computes independently, from
  a different integrator written before the report existed.

    quantity                          report      here      agreement
    ---------------------------------------------------------------------
    surface-confined peak, k0 = 100/s   584 nA   {ref:6.1f} nA    {abs(ref-584)/584*100:5.2f}%
    quasi-reversible maximum, k0/f      0.9      {qrm:5.2f}       {abs(qrm-0.9)/0.9*100:5.1f}%

  Those are not soft agreements. The first is two independent numerical schemes
  -- their finite-element surface solver and this closed-form half-cycle
  integrator -- landing within a quarter of a percent on the same voltammogram.
  Their own third implementation (584.15 nA from the coupled FEM) is inside the
  same bracket. Three codes, one answer.

  The second matters more for design. The quasi-reversible maximum is the
  reason a sensor cannot be improved by making electron transfer faster, and
  both models put it at k0 slightly BELOW the interrogation frequency. This
  package had to have the electron count corrected in the Butler-Volmer
  exponent before it agreed; the report arrived there independently. That is
  the single most trustworthy number in either document.
""")

    print("""  ONE DISCREPANCY WITH A CLEAN CAUSE.

  The report's closed/open contrast is 7.55x under endpoint sampling. This
  package gives a larger number, and the whole difference is the transfer
  coefficient:
""")
    print("     alpha    closed peak    closed/open contrast")
    print("     " + "-" * 50)
    for a in (0.37, 0.44, 0.50, 0.60):
        c, o = peak_na(100.0, alpha=a), peak_na(2.0, alpha=a)
        tag = "   <- the report's 7.55x sits here" if 7.0 <= c / o <= 8.0 else ""
        print(f"     {a:4.2f}   {c:8.1f} nA      {c/o:6.2f}x{tag}")
    c37, o37 = peak_na(100.0, alpha=0.37), peak_na(2.0, alpha=0.37)
    print(f"""
  The report is using the textbook alpha of about 0.5. Methylene blue's
  transfer coefficient has been measured, and it is 0.37 +/- 0.02
  (Dauphin-Ducharme et al., Langmuir 2017, 33, 4407) -- the same paper the
  COMSOL reference model comes from. At the measured value the contrast is
  {c37/o37:.1f}x rather than 7.55x.

  This correction runs in the report's FAVOUR: its sensor is predicted to work
  better than it claims, by a factor of {(c37/o37)/7.55:.1f}. Worth fixing anyway,
  because alpha also sets the peak position and the two peaks' separation, and
  the read-at-the-released-potential trick depends on that separation.
""")


# ---------------------------------------------------------------------------
def section_beta():
    print("=" * 78)
    print("2. THE TUNNELLING DECAY CONSTANT IS TEN TIMES TOO SMALL")
    print("=" * 78)
    print("""
  M6 Part 2 states its model explicitly, which makes it checkable:

      k0(d) = 1e3 * exp( -beta * (d - 1.5 nm) ) per second,  beta = 1.0 /nm

  First, confirm that is really the formula behind the table, rather than a
  simplified caption:
""")
    print("     site   d_c (nm)   k0_c reported   k0_c from the formula   k0_o/k0_c")
    print("     " + "-" * 68)
    for name, dc, k0c, k0o, _, _, _ in SITES:
        pred = K0_REF * np.exp(-BETA_THEIRS_PER_NM * (dc - D_REF_NM))
        print(f"     {name}   {dc:6.2f}   {k0c:11.0f}    {pred:18.1f}   "
              f"{k0o/k0c:9.4f}")
    print(f"""
     every ratio k0_o/k0_c equals {np.exp(-BETA_THEIRS_PER_NM*DISPLACEMENT_NM):.4f}
     = exp(-1.0 * 2.64), so beta = 1.0 per NANOMETRE throughout. Confirmed.

  Now the value. Electron tunnelling through a saturated organic bridge or a
  protein decays with beta of about 1.0 to 1.4 per ANGSTROM, which is 10 to 14
  per nanometre. That is the value in the E-AB literature this project is built
  on -- Dauphin-Ducharme and Plaxco use beta = 1 per Angstrom to convert
  methylene blue rate constants between tether lengths, and it is the value
  this package's own swv module has carried from the start.

  beta = 1 per nanometre is not a typo for nothing: it is the right order for a
  fully conjugated molecular wire. It is wrong by a factor of ten for a dye on
  an aliphatic linker attached to a helical protein, which is what is being
  modelled here.

  The report's sensitivity scan covers beta = 1.0 +/- 0.2 per nanometre. The
  correct value is outside that window by a factor of eight, so the scan cannot
  catch this.
""")

    print("=" * 78)
    print("3. WHAT THE CORRECTION DOES TO THE SITE RANKING")
    print("=" * 78)
    print("""
  Same geometry, same reference rate, same integrator, same frequency. Only
  beta changes, from 1.0 per nm to 1.0 per Angstrom:
""")
    print("     site   d_c(nm)     REPORT (beta = 1/nm)          CORRECTED (beta = 1/A)")
    print("                      k0_c    closed   open        k0_c      closed     open")
    print("     " + "-" * 74)
    for name, dc, k0c, k0o, pc, po, _ in SITES:
        kc = K0_REF * np.exp(-BETA_LITERATURE_PER_NM * (dc - D_REF_NM))
        ko = K0_REF * np.exp(-BETA_LITERATURE_PER_NM * (dc + DISPLACEMENT_NM - D_REF_NM))
        print(f"     {name}   {dc:5.2f}   {k0c:6.0f}  {pc:7.1f}  {po:6.1f}"
              f"     {kc:9.2e}  {peak_na(kc):8.3g}  {peak_na(ko):8.2g}")

    spread_theirs = SITES[0][2] / SITES[-1][2]
    kc_near = K0_REF * np.exp(-BETA_LITERATURE_PER_NM * (SITES[0][1] - D_REF_NM))
    kc_far = K0_REF * np.exp(-BETA_LITERATURE_PER_NM * (SITES[-1][1] - D_REF_NM))
    print(f"""
  Three things change, and they change the conclusion rather than the decimals.

  1. THE SPREAD. Across the seven sites the report's rate constants vary by
     {spread_theirs:.0f}x. With the correct beta they vary by {kc_near/kc_far:.0e}x. Seven sites
     cannot all sit within a factor of a few of the quasi-reversible maximum
     when they span eleven orders of magnitude. The frequency-site matching map
     -- a smooth gradient of polarity across alpha6 -- cannot exist.

  2. THE POLARITY INVERSION DISAPPEARS. It requires the open state to stay
     bright enough to overtake the closed one. At the correct beta a 2.64 nm
     displacement costs exp(-26.4) = {np.exp(-26.4):.1e} in rate. The open state is
     dark by any measure. Inversion through the tunnelling channel is not
     possible at any of these sites.

  3. THE RECOMMENDATION REVERSES. The report advises R268C or D265C at 100 Hz
     and demotes K250C to a high-frequency control. At the correct beta R268C
     and D265C are the SILENT end of their own table -- sub-picoamp, below any
     instrument -- and K250C, the closest site, is the only one with a usable
     current.

  This is not a matter of calibration. The reference rate 1e3 /s is
  acknowledged as uncalibrated, so it is fair to ask whether a larger value
  rescues the far sites. It does not: putting R268C at the quasi-reversible
  maximum would need a contact rate of about {86*np.exp(BETA_LITERATURE_PER_NM*(SITES[-1][1]-D_REF_NM)):.0e} per second, which is
  far above a molecular attempt frequency. No calibration reaches it.

  The earlier version of this report recommended K250C / R254C / R257C. The new
  version overturns that in favour of R261C to R268C. With the literature beta,
  the OLD recommendation was right and the correction is an artefact.
""")


# ---------------------------------------------------------------------------
def section_geometry():
    print("=" * 78)
    print("4. THE GEOMETRY, BY CONTRAST, CHECKS OUT")
    print("=" * 78)
    path = analyze.fetch("7CBC", PDB_DIR)
    ca = analyze.load_ca(path, "A")
    cage = np.array([ca[r] for r in sorted(ca) if 2 <= r <= 240])
    centre = cage.mean(axis=0)
    axis = np.linalg.svd(cage - centre)[2][0]
    height = {r: float((ca[r] - centre) @ axis) for r in ca}

    print("""
  It is worth separating the error from the parts that are sound, because the
  distances feeding that formula were computed independently here.

  The report does not say how it measures d, and the choice matters: a straight
  centre-to-centre distance and a projection onto the electrode normal differ by
  a lot for a bundle lying at an angle. Both are computed below. With the anchor
  at residue 142 sitting ON the electrode plane, the projection is the right
  measure, and it is the one that reproduces their table:
""")
    print("     site   CA-CA to 142   projected on normal   + 1.10 nm tether"
          "   report's d_c")
    print("     " + "-" * 78)
    diffs = []
    for name, dc, *_ in SITES:
        r = int(name[1:])
        if r not in ca:
            continue
        straight = float(np.linalg.norm(ca[r] - ca[142])) / 10.0
        proj = abs(height[142] - height[r]) / 10.0
        diffs.append(proj + 1.10 - dc)
        print(f"     {name}   {straight:9.2f} nm   {proj:15.2f} nm"
              f"   {proj+1.10:15.2f} nm   {dc:11.2f} nm")
    print(f"""
     largest disagreement on any site: {max(abs(d) for d in diffs)*10:.1f} Angstrom

  That is a real reproduction, not a coincidence: seven sites, computed from
  the deposited coordinates by a different route, landing within an Angstrom of
  the report's table. The structural half of M6 Part 2 is sound. It is one
  physical constant applied to it that is off by a factor of ten -- which is
  worth stressing, because a reader who finds the beta error might otherwise
  distrust the geometry too, and the geometry is reusable.

  ONE MODELLING SIMPLIFICATION WORTH NAMING. Every site is given the same
  2.64 nm displacement on opening. A latch that swings out about a hinge does
  not displace its N-terminal turn and its C-terminal turn equally; positions
  near the hinge move least. Since the contrast in the report's table is
  identical across sites by construction, that uniform displacement is doing
  real work in the conclusion and it is not justified anywhere.
""")


# ---------------------------------------------------------------------------
def section_part1():
    print("=" * 78)
    print("5. PART 1: THE COUPLED MODEL IS THE BEST WORK IN THE REPORT")
    print("=" * 78)
    print("""
  Solving the tethered reporter and the diffusion field together is a genuine
  addition, and the validation is done properly: the k_des -> 0 limit recovers
  the static answer to 0.06%, mass conservation drifts by 1e-10, and the pure
  surface limit matches a closed form. Those are the right gates and they pass.

  The result -- that exchange between an adsorbed and a released dye can pump
  the open state's current above the closed state's, inverting the contrast at
  k_des/f between 1 and 10 -- is a mechanism neither document had before, and
  it is the kind of failure mode that would otherwise be discovered on a bench.

  Three reservations, in order of how much they matter.

  A. THE TETHER LAYER IS TEN TIMES THICKER THAN THE SAME SECTION'S GEOMETRY.
     Part 1 confines the released dye to a 50 nm layer, which sets the exchange
     channel's effective rate at k0s/H = 2000 /s and is what makes the pumping
     resonance land where it does. Part 2, two paragraphs later, gives the
     tether as 1.10 nm and the furthest site as 4.53 nm. A dye on a 1.1 nm
     tether attached to a 5 nm protein does not explore 50 nm. With H = 5 nm
     the effective rate is 20000 /s, which is 200x the interrogation frequency
     and far past the quasi-reversible maximum, so the pumped current is
     smaller. The resonance may survive; its magnitude and its position in
     k_des/f will not.

  B. THE MAIN SCENARIO IS NOT THE DESIGN. The whole project specifies a
     COVALENTLY attached dye -- the previous report established that
     non-covalent methylene blue migrates too far to read out. A covalent dye
     cannot desorb; what it can do is stick to the carbon and let go while
     staying tethered, which is a real effect and is presumably what is meant.
     But then the right zero-desorption baseline is not a static mixture, it is
     a linker conformational ensemble, which is a different calculation.
     Labelling this the main scenario, rather than a failure-mode analysis of
     dye adsorption, overstates what it covers.

  C. IT IS A DIAGNOSTIC, NOT A DESIGN RULE. The irreversible-escape control --
     contrast collapsing from 7.55x to 1.0x when the linker breaks -- is the
     most immediately useful thing in M6, because it tells an experimentalist
     what a cleaved linker looks like on a voltammogram before they spend weeks
     on it. That deserves to be a headline conclusion; it is currently a
     parenthetical control.
""")


# ---------------------------------------------------------------------------
def section_charging():
    print("=" * 78)
    print("6. THE CHARGING BUDGET SURVIVES EVERYTHING, AND IT IS THE ACTIONABLE ONE")
    print("=" * 78)
    print("""
  Bare carbon at 20 uF/cm2 and 500 ohm gives tau_dl = 0.71 ms. Checked here
  against the report's numbers:
""")
    print("     f (Hz)   charging here   report   faradaic (site-dependent)   fraction")
    print("     " + "-" * 72)
    for f, claim, their_fara in ((100.0, "77%", 41.7), (400.0, "94%", 2227.0)):
        cap = swv.capacitive_background(freq_hz=f, area_cm2=0.071,
                                        c_dl_uf_cm2=20.0, r_solution_ohm=500.0,
                                        sample_window=0.02) * 1e9
        frac = cap / (cap + their_fara) * 100.0
        print(f"     {f:5.0f}   {cap:11.1f} nA   {claim:>6s}   {their_fara:19.1f} nA"
              f"   {frac:7.1f}%")
    passiv = swv.capacitive_background(freq_hz=1000.0, area_cm2=0.071,
                                       c_dl_uf_cm2=5.0, r_solution_ohm=100.0,
                                       sample_window=0.02) * 1e9
    print(f"""
     passivated interface (5 uF/cm2, 100 ohm) at 1 kHz: {passiv:.2g} nA

  The charging currents reproduce. The percentages reproduce once you read them
  as charging over charging-plus-faradaic rather than charging over faradaic;
  the report does not say which, and the two differ by a lot, so the caption
  needs the formula.

  More important: this conclusion does not depend on beta at all. Charging is
  set by the electrode and the electrolyte, not by where the dye sits, and it
  swamps a few-hundred-nanoamp faradaic signal on a 3 mm disk at any frequency
  worth using. PASSIVATE THE ELECTRODE is the one recommendation in M6 that is
  safe to act on today, and it happens to be the cheapest to implement.
""")


def main():
    section_agreement()
    section_beta()
    section_geometry()
    section_part1()
    section_charging()
    print("=" * 78)
    print("""
  SUMMARY

  KEEP
    the coupled desorption model and its validation gates -- new physics,
      properly checked
    the linker-breakage signature (contrast collapsing to 1.0x) as a bench
      diagnostic, promoted to a headline
    passivate the electrode: charging is 70-95% of the signal on bare carbon
    the quasi-reversible maximum at k0 slightly below f, which this package
      confirms independently at 0.88 against the report's 0.9

  FIX BEFORE ANYONE ORDERS A CONSTRUCT
    beta = 1.0 /nm should be 1.0 /Angstrom. This is the whole of M6 Part 2:
      it creates the polarity inversion, the frequency-site map, and the
      recommendation of R268C/D265C -- all three of which reverse on correction
    alpha = 0.5 should be the measured 0.37 for methylene blue; the contrast is
      better than reported, not worse
    the 50 nm tether layer in Part 1 contradicts the 1.1 nm tether in Part 2
    the uniform 2.64 nm displacement across all seven sites needs a hinge model
      or an explicit caveat

  NET
    M6 Part 1 advanced the project. M6 Part 2 should not be acted on until beta
    is corrected, and when it is, the previous report's recommendation comes
    back: the dye goes at the N-terminal end of alpha6, as close to the anchor
    as the chemistry allows.
""")
    print("=" * 78)


if __name__ == "__main__":
    main()
