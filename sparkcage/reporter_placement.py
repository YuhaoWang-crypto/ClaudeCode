"""
Where to attach the methylene blue, and why the measured signal is small.

    python3 sparkcage/reporter_placement.py

Answers three questions that come up as soon as someone asks how hard this
project actually is:

1. Methylene blue is attached COVALENTLY, not bound in a designed pocket.
   So the design question is not "what sequence binds MB" but "which residue
   carries the cysteine". This script ranks candidate residues by how far they
   move when the protein switches.

2. The known MB-protein co-crystal structures are surveyed, to show what an MB
   binding site looks like and why none of them is reusable as a module.

3. Real protein E-AB sensors deliver only 4-30% signal change against domain
   motions of 11-21 A. THREE different mechanisms explain that, they need
   different fixes, and the literature favours the third one.

IMPORTANT CAVEAT ON SECTION 2: the displacement ranking assumes the signal
comes from rigid-body motion of the attachment residue.  The only published
protein E-AB study (Kang/Plaxco, JACS 2017) reports the opposite placement
rule -- best gain with the reporter PROXIMAL to the binding site, via steric
blocking.  Section 4 lays this out.  Do not treat section 2 as the answer;
treat it as one of two competing hypotheses to test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from structure import analyze              # noqa: E402
from echem import swv                      # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"

# The six PDB entries whose chemical component MBT is methylene blue, with the
# composition of the pocket around one MB copy (residues within 4.5 A).
MB_STRUCTURES = [
    ("2YVE", "multidrug-binding transcriptional regulator", 1.40, 12, 3, 4),
    ("3ZYX", "human monoamine oxidase B", 2.20, 14, 4, 0),
    ("4XYP", "piscine viral fusion protein", 2.10, None, None, None),
    ("5ACM", "Mcg immunoglobulin variable domain", 1.05, 10, 8, 2),
    ("5DLP", "acetylcholinesterase (no PEG)", 2.70, 8, 7, 1),
    ("5E4T", "acetylcholinesterase (with PEG)", 2.43, None, None, None),
]


def survey_mb_structures():
    print("=" * 74)
    print("1. KNOWN METHYLENE-BLUE PROTEIN STRUCTURES (PDB ligand code MBT)")
    print("=" * 74)
    print("\n  PDB    resolution  pocket  aromatic  acidic   what it is")
    print("  " + "-" * 70)
    for pid, what, res, n, arom, acid in MB_STRUCTURES:
        n_s = f"{n:5d}" if n else "    -"
        a_s = f"{arom:6d}" if arom is not None else "     -"
        c_s = f"{acid:5d}" if acid is not None else "    -"
        print(f"  {pid}   {res:5.2f} A   {n_s}   {a_s}   {c_s}   {what}")
    print("\n  Every one of these binds MB adventitiously or promiscuously:")
    print("  a multidrug regulator whose job is to bind many ligands, an enzyme")
    print("  MB happens to inhibit, a famously promiscuous antibody cavity, and")
    print("  a crystallisation additive. None is a purpose-built MB binder and")
    print("  none has the affinity a sensor would need.")
    print("\n  What they do agree on is the recognition motif: MB is a flat")
    print("  aromatic CATION, so it is read out by aromatic stacking plus an")
    print("  acidic counter-charge. That is also the motif that recognises most")
    print("  other flat drug-like cations, which is exactly why every protein in")
    print("  this list is promiscuous. Designing a SPECIFIC MB pocket means")
    print("  fighting that, for no benefit over a covalent bond.")


def rank_attachment_sites(top_n=15):
    print("\n" + "=" * 74)
    print("2. WHERE TO PUT THE CYSTEINE (MBP clamshell, N-lobe on the electrode)")
    print("=" * 74)
    o = analyze.load_ca(analyze.fetch("1OMP", PDB_DIR))
    c = analyze.load_ca(analyze.fetch("1ANF", PDB_DIR))
    common = sorted(set(o) & set(c))
    anchor = analyze._expand(analyze.MBP_N_DOMAIN, set(common))
    rot, pc, qc, _ = analyze.kabsch(np.array([o[r] for r in anchor]),
                                    np.array([c[r] for r in anchor]))
    aligned = {r: (rot @ (o[r] - pc)) + qc for r in common}
    disp = {r: float(np.linalg.norm(aligned[r] - c[r])) for r in common}
    centroid = np.mean([c[r] for r in common], axis=0)

    order = sorted(disp.items(), key=lambda kv: -kv[1])
    print("\n  residue   displacement   distance from centroid (exposure proxy)")
    print("  " + "-" * 66)
    for r, d in order[:top_n]:
        print(f"   {r:4d}       {d:5.1f} A            {np.linalg.norm(c[r]-centroid):5.1f} A")

    vals = np.array(list(disp.values()))
    print(f"\n  median displacement over all residues : {np.median(vals):5.1f} A")
    print(f"  residues moving more than 15 A        : {np.mean(vals > 15)*100:4.0f}%")
    print(f"  residues moving less than  3 A        : {np.mean(vals < 3)*100:4.0f}%")
    print("\n  Read-out: nearly half of all positions barely move. Picking an")
    print("  attachment site at random is close to a coin flip on whether the")
    print("  sensor works at all. Under a DISPLACEMENT mechanism the good sites")
    print("  cluster in three loops around residues 134-143, 200-202 and")
    print("  352-356, all well exposed. Under a STERIC-BLOCKING mechanism the")
    print("  ranking is different and probably inverted -- see section 4 before")
    print("  committing to any of these positions.")
    print("\n  CAVEAT: this ranks TOTAL displacement. What actually sets the")
    print("  signal is the component NORMAL to the electrode surface, which")
    print("  depends on how the protein is oriented once immobilised. Use this")
    print("  to shortlist, then decide with molecular dynamics on the real")
    print("  surface -- which is what the quote budgets 6-8 weeks for.")
    return disp


def conjugation_chemistry():
    print("\n" + "=" * 74)
    print("3. CONJUGATION CHEMISTRY: THE TWO SCAFFOLDS ARE NOT EQUAL")
    print("=" * 74)
    aa3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
           "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
           "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
           "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
           "MSE": "M"}
    for pid, label in [("1ANF", "MBP clamshell"), ("7CBC", "sCageHA (LucCage family)")]:
        path = analyze.fetch(pid, PDB_DIR)
        res = []
        for line in open(path):
            if line.startswith("SEQRES") and line[11] == "A":
                res.extend(line[19:].split())
        seq = "".join(aa3.get(x, "X") for x in res)
        cys = [i + 1 for i, ch in enumerate(seq) if ch == "C"]
        print(f"\n  {label} ({pid}): {len(seq)} aa")
        print(f"    cysteines : {len(cys)} {cys if cys else '(none)'}")
        print(f"    lysines   : {seq.count('K')}")
    print("\n  Read-out: MBP is CYSTEINE-FREE. Introduce one cysteine anywhere and")
    print("  maleimide-MB labels it, site-specifically, at 1:1 stoichiometry, with")
    print("  nothing to compete. The LucCage-family construct carries two native")
    print("  cysteines inside the grafted binder, so it needs them mutated out")
    print("  first, and any mislabelling puts a reporter somewhere that does not")
    print("  move. On conjugation chemistry alone the clamshell is the easier")
    print("  starting point.")


def discriminate_signal_loss():
    print("\n" + "=" * 74)
    print("4. WHY MEASURED SIGNAL IS ONLY 4-30%: THREE HYPOTHESES, DIFFERENT FIXES")
    print("=" * 74)
    swv.k0_from_distance.__defaults__ = (3.0e2, 1.0, 5.0)
    cell = dict(area_cm2=0.03, freq_hz=87.0, e_step=0.002)
    gamma = 2.0e-12

    print("\n  Hypothesis A: the reporter genuinely barely moves.")
    print("    displacement   signal change")
    print("    " + "-" * 34)
    base = abs(swv.sensor_response_distributed(
        0.0, 6.0, 6.0, 2.5, gamma_total_mol_cm2=gamma, **cell)["peak_current"])
    for dd in [0.5, 1.0, 2.0, 3.0, 4.0, 8.0]:
        b = abs(swv.sensor_response_distributed(
            1.0, 6.0, 6.0 + dd, 2.5, gamma_total_mol_cm2=gamma,
            **cell)["peak_current"])
        ch = (b - base) / base * 100.0
        flag = "   <- in the measured band" if 4 <= abs(ch) <= 30 else ""
        print(f"      {dd:4.1f} A       {ch:7.1f}%{flag}")

    print("\n  Hypothesis B: the reporter moves the full 8 A, but most probes")
    print("  are non-responsive (wrong orientation, denatured, mislabelled).")
    print("    dead fraction   signal change")
    print("    " + "-" * 34)
    for dead in [0.0, 0.25, 0.5, 0.75, 0.9, 0.95]:
        live_lo = abs(swv.sensor_response_distributed(
            0.0, 6.0, 14.0, 2.5, gamma_total_mol_cm2=gamma * (1 - dead),
            **cell)["peak_current"]) if dead < 1 else 0.0
        live_hi = abs(swv.sensor_response_distributed(
            1.0, 6.0, 14.0, 2.5, gamma_total_mol_cm2=gamma * (1 - dead),
            **cell)["peak_current"]) if dead < 1 else 0.0
        stuck = abs(swv.sensor_response_distributed(
            0.0, 6.0, 14.0, 2.5, gamma_total_mol_cm2=gamma * dead,
            **cell)["peak_current"]) if dead > 0 else 0.0
        ch = ((live_hi + stuck) - (live_lo + stuck)) / (live_lo + stuck) * 100.0
        flag = "   <- in the measured band" if 4 <= abs(ch) <= 30 else ""
        print(f"      {dead*100:5.1f}%        {ch:7.1f}%{flag}")

    print("\n  Hypothesis C: STERIC BLOCKING, and this is the one the only")
    print("  protein precedent actually supports.")
    print("    Kang/Plaxco (JACS 2017, 139, 12113) scanned eight single-cysteine")
    print("    positions on CheY (M17C, E37C, T71C, A80C, G89C, K91C, K97C,")
    print("    E117C) and found gain was LARGEST WHEN MB SITS PROXIMAL TO THE")
    print("    BINDING SITE, attributing it to the bound analyte sterically")
    print("    blocking MB's approach to the electrode -- not to rigid-body")
    print("    displacement of the attachment residue.")
    print("    Consequence: the displacement ranking above optimises a DIFFERENT")
    print("    mechanism from the one that demonstrably works. Under C the right")
    print("    site is adjacent to the binding interface, which for MBP means the")
    print("    maltose-site residues, NOT the 134-143 / 200-202 / 352-356 loops.")
    print("    Treat the two rankings as competing hypotheses and test both.")
    print("    The model in echem/ is agnostic: blocking and displacement both")
    print("    act by lowering the effective electron-transfer rate, so the")
    print("    voltammetry is described identically. What differs is WHERE to")
    print("    put the cysteine, and that is the expensive decision.")
    print("\n  A, B and C all reproduce the measured 4-30%, and they demand")
    print("  different fixes:")
    print("  A is fixed by moving the label to a higher-displacement residue;")
    print("  B by changing immobilisation chemistry and surface density; C by")
    print("  moving the label TOWARD the binding site instead of away from it.")
    print("  Guessing wrong wastes a design cycle.")
    print("\n  They are distinguishable in ONE experiment. Under A the whole")
    print("  monolayer shifts a little, so the voltammogram stays single-peaked")
    print("  and simply changes height. Under B the monolayer splits into a large")
    print("  non-switching population and a small switching one, so the")
    print("  voltammogram is a SUM OF TWO components with different k0, and it")
    print("  broadens rather than just shrinking. Fit the measured voltammogram")
    print("  with a k0 distribution (RedoxPySolid does exactly this) and the")
    print("  shape tells you which mechanism you have before you redesign")
    print("  anything.")

    print("\n  Distance dispersion alone is NOT the explanation, for the record:")
    print("    ensemble width   signal change (true 8 A displacement)")
    print("    " + "-" * 50)
    for s in [0.5, 2.5, 6.0, 9.0]:
        a = abs(swv.sensor_response_distributed(
            0.0, 6.0, 14.0, s, gamma_total_mol_cm2=gamma, **cell)["peak_current"])
        b = abs(swv.sensor_response_distributed(
            1.0, 6.0, 14.0, s, gamma_total_mol_cm2=gamma, **cell)["peak_current"])
        print(f"      {s:4.1f} A          {(b-a)/a*100:7.1f}%")
    print("    Even a 9 A wide ensemble still leaves about half the signal, so")
    print("    linker floppiness by itself cannot account for the measured range.")


if __name__ == "__main__":
    survey_mb_structures()
    rank_attachment_sites()
    conjugation_chemistry()
    discriminate_signal_loss()
    print("\n" + "=" * 74)
