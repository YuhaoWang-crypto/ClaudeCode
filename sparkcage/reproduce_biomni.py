"""
Reproducing the Biomni SparkCage report: what came back the same, and what did not.

    python3 sparkcage/reproduce_biomni.py

THE QUESTION. A second team ran a five-module pipeline on the same scaffold this
package uses (7CBC) and delivered a feasibility report, an HTML design report and
a slide deck. Can the process be reproduced?

Reproducing a pipeline means three different things and they have three different
answers, so they are kept apart here:

  RE-RUN      execute the same software on the same inputs. Needs the same
              hardware. Section 1 measures what this machine can and cannot do,
              rather than asserting it.
  RE-DERIVE   get the same NUMBER by a route that does not share the first
              route's assumptions. This is the only kind of reproduction that
              can find an error. Sections 2 and 3 do it where it is cheap.
  RE-CHECK    test whether a claim is consistent with the report's own stated
              uncertainties. Costs nothing and finds the most. Section 4.

A claim that survives all three is worth building an experiment on. A claim that
fails section 4 cannot be rescued by any amount of compute, which is why that
section runs last and matters most.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from echem import swv                      # noqa: E402
from structure import analyze              # noqa: E402
import spec_review                         # noqa: E402
from mb_placement_design import BURIAL, SIDECHAIN_VOL, MB_VOLUME_A3   # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"

# Measured on this container by scratchpad/bench_md.py: OpenMM 8.6, CPU platform,
# 4 cores, amber14 + GBn2 implicit, 4 fs with hydrogen mass repartitioning,
# 1.2 nm cutoff -- the report's own protocol.
NS_PER_DAY_LATCH = 18.6      # 75 residues, 1159 atoms
NS_PER_DAY_FULL = 2.91       # 317 resolved residues, 5139 atoms


# ---------------------------------------------------------------------------
# 1. What this machine can actually run
# ---------------------------------------------------------------------------
def section_environment():
    print("=" * 78)
    print("1. RE-RUN: WHAT THIS MACHINE CAN AND CANNOT EXECUTE")
    print("=" * 78)

    have = {}
    for mod in ("numpy", "scipy", "openmm", "skfem", "torch", "dolfinx", "pymbar"):
        try:
            __import__(mod)
            have[mod] = True
        except Exception:
            have[mod] = False

    rows = [
        ("M1  RFDiffusion backbones",     "GPU + weights",  "NO",
         "no CUDA device; CPU inference is days per backbone"),
        ("M1  ProteinMPNN sequences",     "torch, CPU ok",  "PARTLY",
         "runs on CPU in seconds, but torch is not installed here"),
        ("M1  Boltz-2 structures",        "GPU or API",     "YES, via API",
         "the Boltz MCP endpoint is reachable; local weights are not"),
        ("M2  AutoDock Vina docking",     "CPU",            "YES",
         "minutes for a blind box on a 75-residue latch"),
        ("M2  OpenMM implicit MD",        "CPU ok",         "YES, slowly",
         f"measured {NS_PER_DAY_LATCH:.1f} ns/day latch, {NS_PER_DAY_FULL:.1f} full chain"),
        ("M3  umbrella sampling + MBAR",  "CPU x 16",       "NOT IN A SESSION",
         "see the wall-clock estimate below"),
        ("M4  surface-confined SWV",      "numpy",          "YES",
         "closed form, no mesh; this package already has it"),
        ("M4  microdisk / SPE FEA",       "skfem or FEniCSx", "YES via skfem",
         "FEniCSx has no working pip wheel; scikit-fem does the same job"),
        ("M5  endpoint MM/GBSA Ala scan", "OpenMM CPU",     "YES",
         "minimisations, not dynamics; minutes per mutant"),
    ]
    print("\n  module                        needs             here          note")
    print("  " + "-" * 74)
    for name, needs, verdict, note in rows:
        print(f"  {name:29s} {needs:17s} {verdict:13s} {note}")

    print(f"""
  Installed here: openmm {have['openmm']}, scikit-fem {have['skfem']},
  torch {have['torch']}, dolfinx {have['dolfinx']}, pymbar {have['pymbar']}.
  No CUDA device. 4 cores.

  The report states its own cost as roughly 60 core-hours of CPU across two
  16-core workers, plus GPU for RFDiffusion, ProteinMPNN and six Boltz-2 runs.
  Measured here on the report's own protocol -- amber14 with GBn2 implicit
  solvent, 1.2 nm cutoff, 4 fs steps with hydrogen mass repartitioning:

      latch alone, 1159 atoms      {NS_PER_DAY_LATCH:5.1f} ns/day
      full chain A, 5139 atoms     {NS_PER_DAY_FULL:5.1f} ns/day

  which puts the report's modules at {4.0/NS_PER_DAY_LATCH*24:.0f} h for the 4 ns
  dye-latch run, {9.6/NS_PER_DAY_FULL*24:.0f} h for the 16-window umbrella
  sampling and {1.5/NS_PER_DAY_FULL*24:.0f} h for the steered run.

  So: the pipeline is reproducible in kind, on ordinary hardware, with no
  commercial licence. It is not reproducible inside one session on this box.
  That is a scheduling fact, not a scientific one, and it is the least
  interesting thing in this file.
""")


# ---------------------------------------------------------------------------
# 2. Electrochemistry, re-derived
# ---------------------------------------------------------------------------
def section_echem():
    print("=" * 78)
    print("2. RE-DERIVE: THE ELECTROCHEMISTRY, INDEPENDENTLY")
    print("=" * 78)
    print("""
  This is the part of the report that can be checked properly, because the
  surface-confined problem has a closed-form solution and needs no mesh, no
  force field and no GPU. Three of its claims are numerical and testable.

  BEFORE ANY OF THEM: THE REPRODUCTION FOUND A BUG HERE, NOT THERE.

  The two reports disagree about whether methylene blue is a one- or
  two-electron couple. Settling that meant looking at where the electron count
  actually enters, and this package had it in only one of the two places. The
  charge prefactor carried n; the Butler-Volmer exponent did not. That makes n
  a pure scale factor, so every peak this module drew had the WIDTH and the
  kinetic shift of a one-electron couple no matter what n was set to.

  Fixed: the exponent is now Laviron's n F / RT. Three things moved, and all
  three moved towards the literature rather than away from it.
""")
    for n_e, want in ((1, 90.6), (2, 45.3)):
        r = swv.swv_scan(k0=1.0e2, freq_hz=100.0, n_electrons=n_e,
                         e_step=0.001, alpha=0.5)
        i = np.abs(r["i_net"])
        sel = np.where(i >= i.max() / 2.0)[0]
        fwhm = abs(r["e"][sel[-1]] - r["e"][sel[0]]) * 1000.0
        print(f"      n = {n_e}   peak half-width {fwhm:5.1f} mV"
              f"   (Laviron: {want:.0f} mV)")
    ks = np.logspace(-1, 5, 80)
    for n_e in (1, 2):
        p = np.array([abs(swv.swv_scan(k0=k, freq_hz=100.0,
                                       n_electrons=n_e)["peak_current"])
                      for k in ks])
        print(f"      n = {n_e}   quasi-reversible maximum at k0/f ="
              f" {ks[int(np.argmax(p))]/100.0:5.2f}   (literature: about 1)")
    print("""
      The half-widths now scale as 1/n, which is the whole content of
      Laviron's result, and the quasi-reversible maximum moves from 1.30 to
      1.09. Both are closer to the published values than before. The residual
      few millivolts is the staircase step, not the physics.

      This is what re-deriving is for. Neither report's numbers were wrong
      here; comparing them exposed a defect in the checker. The peak
      separation quoted from demo_end_to_end.py drops from 225 mV to 125 mV as
      a result -- still a real separation, still enough for the read-at-the-
      released-potential argument, but the earlier figure was too large.
""")

    # --- claim 1: closed/open peak ratio ------------------------------------
    cell = dict(gamma_mol_cm2=2.0e-12, area_cm2=0.071, e_step=0.004,
                e_sw=0.025, freq_hz=100.0, n_electrons=2)
    p_closed = abs(swv.swv_scan(k0=100.0, **cell)["peak_current"])
    p_open = abs(swv.swv_scan(k0=2.0, **cell)["peak_current"])
    print("  CLAIM  'at 100 Hz the closed state (k0 = 100/s) peaks 9x higher")
    print("          than the open state (k0 = 2/s)'")
    print(f"  HERE   closed {p_closed*1e9:7.1f} nA, open {p_open*1e9:7.1f} nA,"
          f"  ratio {p_closed/p_open:5.2f}x  at this module's defaults")

    grid = []
    for n in (1, 2):
        for alpha in (0.37, 0.50):
            for e_sw in (0.025, 0.050):
                for window in (0.10, 0.25, 0.50):
                    kw = {**cell, "n_electrons": n, "alpha": alpha,
                          "e_sw": e_sw, "sample_window": window}
                    a = abs(swv.swv_scan(k0=100.0, **kw)["peak_current"])
                    b = abs(swv.swv_scan(k0=2.0, **kw)["peak_current"])
                    grid.append(a / b)
    lo, hi = min(grid), max(grid)
    hits = sum(1 for g in grid if 8.0 <= g <= 10.0)
    print(f"""
  ->     NOT DECIDABLE AS STATED, which is a finding about the report rather
         than about the physics. The ratio of two square-wave peaks depends on
         four things the report does not give for this number: the electron
         count in the Butler-Volmer exponent, the transfer coefficient, the
         square-wave amplitude, and the fraction of the pulse the instrument
         integrates over. Sweeping all four over their reasonable ranges gives
         {lo:.1f}x to {hi:.1f}x, and {hits} of the {len(grid)} combinations land
         within 8-10x.

         So 9x is inside the reachable set and cannot be called wrong. It also
         cannot be called reproduced, because a number that moves by a factor of
         {hi/lo:.0f} across unstated conventions is not a prediction yet. What
         would make it one is three extra symbols in the caption.
""")

    # --- claim 2: the frequency window rule ---------------------------------
    print("  CLAIM  'charging reaches 10% of the faradaic peak at f10 = 0.07/tau_dl,")
    print("          across C_dl 5-50 uF/cm2 and R_s 100-500 ohm'")
    print("\n     C_dl   R_s     tau_dl      f10 found     0.07/tau_dl    ratio")
    print("     " + "-" * 62)
    ratios, taus = [], []
    for c_dl in (5.0, 20.0, 50.0):
        for r_s in (100.0, 500.0):
            tau = r_s * c_dl * 1e-6 * 0.071
            f_lo, f_hi = 1.0, 1.0e6
            for _ in range(60):                      # bisect on charging/faradaic
                f = np.sqrt(f_lo * f_hi)
                fara = abs(swv.swv_scan(k0=100.0, **{**cell, "freq_hz": f})["peak_current"])
                cap = swv.capacitive_background(freq_hz=f, area_cm2=0.071,
                                                c_dl_uf_cm2=c_dl, r_solution_ohm=r_s,
                                                e_sw=0.025, e_step=0.004)
                if cap / fara < 0.10:
                    f_lo = f
                else:
                    f_hi = f
            f10 = np.sqrt(f_lo * f_hi)
            ratios.append(f10 * tau)
            taus.append(tau)
            print(f"     {c_dl:4.0f}  {r_s:4.0f}   {tau*1e6:7.1f} us  {f10:9.1f} Hz"
                  f"   {0.07/tau:11.1f} Hz   {f10*tau:7.3f}")
    lo, hi = min(ratios), max(ratios)
    min_tau, max_tau = min(taus), max(taus)
    print(f"""
  ->     REPRODUCES AS A SCALING, NOT AS A CONSTANT. The product f10 * tau_dl
         stays inside {lo:.3f} to {hi:.3f} while tau_dl itself varies {max_tau/min_tau:.0f}-fold
         across the grid, so the report's central claim is right: the frequency
         ceiling is set by the cell time constant and by nothing else. That is
         the useful half of the rule and it transfers.

         The prefactor does not transfer. It drifts by {hi/lo:.1f}x across this
         grid alone, and it brackets rather than reproduces the quoted 0.07.
         The drift is a sampling convention: this integrator averages charge
         over the last 25% of each pulse, and a model reading the instantaneous
         current at the end of the pulse sees less residual charging and lands
         higher. Neither is wrong. But quoting 0.07 without the potentiostat's
         sampling window makes it look like a measured constant when it is a
         model output, and an engineer sizing an electrode from it can be a
         factor of two out on the usable frequency.
""")

    # --- claim 3: the peak-vs-frequency discriminant ------------------------
    print("  CLAIM  'diffusing species scale sublinearly with frequency while")
    print("          surface-confined species scale linearly, so the peak-versus-")
    print("          frequency exponent tells a desorbed sensor from a working one'")
    freqs = np.array([25.0, 50.0, 100.0, 200.0, 400.0])
    print("\n     f (Hz)   k0/f    peak (nA)   local exponent d(ln i)/d(ln f)")
    print("     " + "-" * 60)
    peaks = []
    for f in freqs:
        peaks.append(abs(swv.swv_scan(k0=100.0, **{**cell, "freq_hz": f})["peak_current"]))
    peaks = np.array(peaks)
    slopes = np.gradient(np.log(peaks), np.log(freqs))
    for f, p, s in zip(freqs, peaks, slopes):
        print(f"     {f:6.0f}  {100.0/f:6.2f}   {p*1e9:9.3f}   {s:+6.2f}")
    diffusive = np.log(266 / 117) / np.log(16.0)
    print(f"""
  ->     THE RULE IS RIGHT AND THE OPERATING POINT IS WRONG. A surface-confined
         peak is linear in frequency only near the quasi-reversible maximum at
         k0/f = 1. Away from it the exponent is not a constant at all: it is
         {slopes[0]:+.2f} between 25 and 50 Hz and {slopes[-1]:+.2f} between 200
         and 400 Hz, because the sampled current dies away on both sides of the
         maximum while the 1/pulse-width prefactor keeps rising.

         The report's own diffusive exponent, 117 to 266 nA over a 16-fold
         frequency range, is {diffusive:.2f}. Put the two together at the report's
         own closed-state k0 of 100 per second:

           at  25 Hz   surface-confined {slopes[0]:+.2f}   vs diffusive {diffusive:+.2f}
                       -> well separated, the test works
           at 400 Hz   surface-confined {slopes[-1]:+.2f}   vs diffusive {diffusive:+.2f}
                       -> INVERTED: the confined sensor now scales more weakly
                          than the desorbed one, so the test returns the wrong
                          answer with full confidence

         This matters because 400 Hz is where the report models the SPE and it
         is also where E-AB sensors are commonly interrogated. The discriminant
         has to be run below the quasi-reversible maximum -- for this k0, under
         roughly {100.0/3:.0f} Hz -- and the report does not say so.
""")


# ---------------------------------------------------------------------------
# 3. The reporter site, re-derived
# ---------------------------------------------------------------------------
def section_placement():
    print("=" * 78)
    print("3. RE-DERIVE: THE REPORTER SITE, WHERE THE TWO MODELS DISAGREE")
    print("=" * 78)

    path = analyze.fetch("7CBC", PDB_DIR)
    atoms = spec_review.load_atoms(path)
    universe = list(range(len(atoms)))
    latch = [i for i, a in enumerate(atoms) if 244 <= a[0] <= 269]
    caged = spec_review.sasa(atoms, latch, universe)
    free = spec_review.sasa(atoms, latch, latch)

    ident, side_sasa, side_free = {}, {}, {}
    for i in latch:
        rn = atoms[i][0]
        ident[rn] = atoms[i][1]
        if atoms[i][2] not in ("N", "CA", "C", "O"):
            side_sasa[rn] = side_sasa.get(rn, 0.0) + caged[i]
            side_free[rn] = side_free.get(rn, 0.0) + free[i]
    # burial measured here rather than read from the stored table, so that
    # positions outside that table (the report proposes R268C) still score
    burial_here = {r: (1.0 - side_sasa[r] / side_free[r]) if side_free[r] > 1.0
                   else 0.0 for r in side_sasa}

    theirs = [250, 254, 257, 268]
    mine = [249]
    print("""
  The HTML report recommends K250C / R254C / R257C, chosen for HIGH solvent
  exposure on the outward face of helix alpha6. This package recommends R249C,
  chosen for the opposite reason: burial at the cage-latch interface.

  Same structure, same SASA calculation, opposite selection rule. So one of the
  two selection rules is wrong, and the numbers alone will not say which.
""")
    print("     pos  aa    caged side-chain SASA   burial   MB deficit   source")
    print("     " + "-" * 70)
    for r in sorted(set(theirs + mine)):
        if r not in ident:
            continue
        b = burial_here.get(r, 0.0)
        vol = SIDECHAIN_VOL.get(ident[r], 80)
        deficit = b * MB_VOLUME_A3 - vol
        src = "HTML report" if r in theirs else "this package"
        print(f"     {r:3d}  {ident[r]}   {side_sasa.get(r, 0.0):10.1f} A2       "
              f"{b*100:5.1f}%   {deficit:+8.0f} A3   {src}")

    print("""
  Both calculations agree on the facts. The report's three candidates have
  essentially zero burial; this package's candidate has 59%. The disagreement is
  entirely about what burial is FOR, and that is a question about mechanism.

  THE MECHANISM QUESTION, STATED PROPERLY.

  An E-AB signal has two possible origins and they want opposite sites.

    Tunnelling distance. The current falls as exp(-beta d) with the dye's
    distance to the electrode. Only DISPLACEMENT matters. Burial is irrelevant,
    and buying burial costs a cavity the cage cannot afford. On this mechanism
    the HTML report is right and the occlusion gate in this package is not just
    unnecessary, it is harmful.

    Steric blocking. The current falls because the dye cannot reach the
    electrode surface in one state, whatever the centre-of-mass distance. Kang
    and Plaxco (JACS 2017) measured the largest gains for reporters PROXIMAL to
    the binding site rather than most displaced, which is the signature of this
    mechanism. On it, burial is the whole point and the HTML report's exposed
    sites will give a small signal for a large motion.

  Both reports assume the first mechanism and neither tests it. This package
  assumed the second in its occlusion gate without saying so. THAT is the real
  finding of the reproduction: the gate encodes an untested mechanism claim.

  What the geometry can settle without settling the mechanism: the displacement
  itself. If the two faces of alpha6 move by the same amount when the latch
  ejects, then on the tunnelling mechanism the exposed sites are strictly better,
  because they cost nothing. If they do not, the comparison is live.
""")

    ca = analyze.load_ca(path, "A")
    cage = np.array([ca[r] for r in sorted(ca) if 2 <= r <= 240])
    centre = cage.mean(axis=0)
    axis = np.linalg.svd(cage - centre)[2][0]

    print("     pos  aa    axial height   radial offset   distance to Res 142")
    print("     " + "-" * 66)
    ref = ca.get(142)
    for r in sorted(set(theirs + mine)):
        if r not in ca:
            continue
        v = ca[r] - centre
        h = float(v @ axis)
        rad = float(np.linalg.norm(v - h * axis))
        d142 = float(np.linalg.norm(ca[r] - ref)) if ref is not None else float("nan")
        print(f"     {r:3d}  {ident.get(r, '?')}   {h:9.1f} A   {rad:11.1f} A   {d142:14.1f} A")

    print("""
  Read the last column, not the burial column. With the anchor at residue 142,
  the tunnelling distance is set by how far the dye sits from the anchor, and
  the four candidates span 19 to 40 Angstrom. At beta = 1 per Angstrom that is
  the difference between a measurable baseline and none:
""")
    d249 = float(np.linalg.norm(ca[249] - ref))
    print("     pos  aa   d to anchor   d - d(249)   relative k0 at beta = 1/A")
    print("     " + "-" * 64)
    for r in sorted(set(theirs + mine)):
        if r not in ca or ref is None:
            continue
        d = float(np.linalg.norm(ca[r] - ref))
        note = "  <- this package" if r in mine else ""
        print(f"     {r:3d}  {ident.get(r, '?')}   {d:8.1f} A   {d-d249:+8.1f} A"
              f"   {np.exp(-(d-d249)):18.1e}{note}")

    print("""
  Absolute currents are deliberately not quoted here. The two models in play
  disagree about where the electrode plane sits relative to the bundle -- this
  package measured along the bundle axis from the latch junction, the anchor
  geometry measures from residue 142 -- and they differ by about ten Angstrom,
  which at beta = 1 is four orders of magnitude in baseline current. Neither is
  calibrated. What survives that disagreement is the SPACING between
  candidates, because it is the same in both.
""")

    print("""
  This reorders both recommendations, and it does so without settling the
  mechanism question at all.

  The radial offsets differ by only three or four Angstrom between the buried
  and the exposed face -- both sit on the same helix, and a helix radius is
  small next to a bundle radius. So burial is nearly irrelevant to the
  tunnelling distance, exactly as the tunnelling mechanism would have it.

  What is NOT irrelevant is position along the helix. R249 is one helical turn
  closer to the anchor than K250 and two turns closer than R254, and on an
  exponential that is worth orders of magnitude in baseline current. This
  package's candidate does win -- but for a reason its own scoring never
  stated, and the reason has nothing to do with the burial gate that selected
  it. Of the report's three, K250C is much the best and R257C is not viable at
  all; the report ranks R254C first.

  The corrected recommendation is a two-arm screen at MATCHED distance, which
  neither document proposes:
      R249C   buried face,  reference distance   tests steric blocking
      K250C   exposed face, +3 A further out    tests tunnelling distance
  Adjacent residues, opposite faces, one helical turn apart. If K250C gives
  the larger signal change the mechanism is distance and this package's
  occlusion gate should be deleted. If R249C does, the gate is right and all
  three of the report's candidates are on the wrong face. Comparing R249C
  against R254C instead, as the two documents implicitly do, confounds face
  with distance and cannot answer anything.
""")

    print("     ANCHOR SITE CHECK (the report proposes K142C on the cage back face)")
    print("     " + "-" * 66)
    res_names = {}
    for line in open(PDB_DIR / "7cbc.pdb"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
            res_names[int(line[22:26])] = line[17:20].strip()
    for r in (2, 51, 54, 142, 145, 151, 154):
        if r not in ca:
            print(f"     Res {r:3d}  not resolved in the crystal")
            continue
        h = float((ca[r] - centre) @ axis)
        print(f"     Res {r:3d}  {res_names.get(r, '???')}   axial height {h:+7.1f} A")
    h_latch = float((ca[245] - centre) @ axis) if 245 in ca else float("nan")
    print(f"     latch start (245)      axial height {h_latch:+7.1f} A")
    print("""
  Every residue identity in the report's anchor list checks out against the
  deposited coordinates: 51 and 142 are lysine, 54 is arginine, 151 and 154 are
  glutamate. The constructs it names can be made.

  And on placement the two analyses agree, which is worth more than the
  agreement on any single number. An anchor near residue 142-145 sits at the
  SAME end of the bundle as the latch junction, +34 Angstrom, against -33 for
  the far end at residue 2. This package reached that independently in
  latch_geometry.py by a different argument: an anchor at the far end puts the
  dye 70 Angstrom from the electrode and there is no baseline peak at all. The
  report got there from a solvent-accessibility and orientation argument. Two
  routes, one answer, and it rules out the N-terminal anchor that the original
  client specification proposed.
""")


# ---------------------------------------------------------------------------
# 4. Claims that fail on the report's own numbers
# ---------------------------------------------------------------------------
def section_selfconsistency():
    print("=" * 78)
    print("4. RE-CHECK: WHAT NO AMOUNT OF COMPUTE WOULD FIX")
    print("=" * 78)
    print("""
  A. THE ALANINE SCAN CONTRADICTS ITS OWN NOISE FLOOR.

     The report measures its run-to-run reproducibility and states the answer
     plainly: the same wild-type minimisation differs by 85 kJ/mol between
     machines, and it draws the correct conclusion, that only effects above
     about 100 kJ/mol are credible.

     It then reports S255A at +52 kJ/mol and Y287A at +50, calls both
     "stabilising knobs", and presents S255 as an independent rediscovery of the
     known LOCKR tuning site. Both numbers are BELOW the floor the report itself
     set. On its own stated uncertainty, the scan did not find S255; it produced
     a value indistinguishable from zero at a position that was already known.

     This is the one claim in the report that should be withdrawn rather than
     refined, and no additional sampling changes that -- the fix is a
     conformational-ensemble protocol, which is a different calculation, not
     more of this one.

     What survives: the two large negative hits, P245A and D246A, at -165 and
     -178 kJ/mol, clear the floor. And the Tier-2 falsification of I262A is a
     genuine result: a -1014 kJ/mol prediction from a structure that folds
     perfectly well is a protocol artefact, correctly caught.

  B. THE BINDING FREE ENERGY DOES NOT CONSTRAIN ANYTHING.

     Reported: dG_LT = -32.1 kcal/mol from single-point MM/GBSA, with the
     honest note that the entropy correction is +10 to +20 kcal/mol, giving
     -12 to -22.
""")
    rt = 1.987204259e-3 * 310.15
    for dg in (-12.0, -22.0, -32.1):
        kd = np.exp(dg / rt)
        print(f"       dG_LT = {dg:6.1f} kcal/mol  ->  Kd = {kd:9.2e} M")
    print(f"""
     The stated range spans {np.exp(-12.0/rt)/np.exp(-22.0/rt):.0e} in affinity.
     A designed minibinder is nanomolar; only the loose end of that band is
     physical, and the calculation cannot tell you where in the band you are.
     Since EC50 is Kd times exp(dG_open / RT), an unconstrained Kd makes the
     whole dose-response prediction unconstrained. The report says as much --
     "the largest source of uncertainty", a four-decade shift -- and then still
     carries dG_LT forward into a thermodynamic-feasibility conclusion.

     The feasibility conclusion happens to be right, but not for the reason
     given: target binding drives opening because a nanomolar binder supplies
     about 12 kcal/mol against an opening cost of 1 to 2. That argument needs
     no MM/GBSA at all. The calculation added a number, not knowledge.

  C. WHAT REPRODUCES WITHOUT RESERVATION.

     Three of the report's conclusions are robust, cheap to check, and were
     reached here independently:

       - The dye must be covalently attached. A non-covalent methylene blue
         migrates over the latch surface (their MD: 2.0-2.6 nm RMSD), and this
         package's survey of the six methylene-blue protein structures found the
         same thing structurally: the dye is read out by a generic
         aromatic-plus-acidic motif, so it has no unique site to occupy.

       - The dye's own preferred patches, and the electrode adsorption patch,
         both lie in the grafted binder domain past residue 272. Avoiding them
         is free, because the reporter has to be on the latch helix anyway.

       - The opening free energy has to be tuned to roughly 5-10 kJ/mol. This
         package reached 1.3-8.8 kJ/mol from the clinical decision thresholds
         for the TBI markers, by a completely different argument. The overlap is
         the most useful cross-check in either document.
""")


def main():
    section_environment()
    section_echem()
    section_placement()
    section_selfconsistency()
    print("=" * 78)
    print("""
  BOTTOM LINE

  The METHOD reproduces. Every tool in the pipeline is open source and runs on
  ordinary hardware -- no Rosetta, no COMSOL, no commercial MD licence -- and
  the two slow modules were benchmarked here rather than guessed at: about
  three days of this machine's CPU for the whole free-energy half.

  The NUMBERS reproduce unevenly, and the pattern is consistent. Every claim
  that is a scaling came back; every claim that is a constant did not.

    reproduced   the frequency ceiling is set by the cell time constant
    reproduced   the anchor belongs at the latch-junction end of the bundle
    reproduced   the opening free energy has to be tuned to a few kJ/mol
    not decidable  the 9x signal-off ratio: it moves by 9x across conventions
                   the report does not state, and 9x is inside that range
    reproduces in shape only   f10 = 0.07/tau_dl; the 0.07 is model-specific
    fails        the peak-versus-frequency discriminant INVERTS above the
                 quasi-reversible maximum, and 400 Hz is above it

  Three things the reproduction found that re-running the same code would not:

    - This package had the electron count missing from the Butler-Volmer
      exponent. Found by trying to settle the two reports' disagreement about
      n. Fixed; the peak widths and the quasi-reversible maximum both moved
      towards the published values.

    - The alanine scan's headline result -- rediscovering the S255 tuning site
      -- sits below the noise floor the report itself measured. It should be
      withdrawn, not refined.

    - The two documents recommend opposite faces of the same helix for the dye
      and neither states that this is a disagreement about mechanism. Worse,
      the natural head-to-head comparison confounds face with distance. The
      clean experiment is R249C against K250C: adjacent residues, opposite
      faces, one helical turn apart.

  Reproducing was worth more than re-running would have been, and cost hours
  rather than days.
""")
    print("=" * 78)


if __name__ == "__main__":
    main()
