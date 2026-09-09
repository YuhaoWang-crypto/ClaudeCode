"""
Does the placement rule generalise, or is it fitted to 7CBC?

    python3 sparkcage/blind_test_7jh5.py           # both scaffolds
    python3 sparkcage/blind_test_7jh5.py 9DCX      # one

WHAT WAS ASKED, AND WHAT ACTUALLY HAPPENED
------------------------------------------
The gates in mb_placement_design.py were all devised while looking at 7CBC, so
they could be convincing and still be fitted to one example. The plan was to
re-run them blind on 7JH5 (Co-LOCKR). Two things went wrong with that plan, and
both are more informative than the test would have been.

1. 7JH5 IS NOT AN INDEPENDENT TEST. Its cage is 97% identical to 7CBC's, the
   first eleven residues of its latch are identical, and the position that
   survives is literally the same arginine, at the same offset of five from the
   latch start. Both derive from the same designed parent scaffold. Re-running
   the rule there is a consistency check, not a test of generalisation, and the
   script now measures and states that rather than letting the verdict
   overstate itself.

2. A GENUINELY INDEPENDENT SCAFFOLD PUT THE RULE OUT OF DOMAIN. 9DCX is a
   designed allosteric switch from a different paper and a different lineage,
   23% identical. The rule scored it and returned an answer, which it should
   not have: 9DCX's latch lies ACROSS the end of its bundle at 86 degrees,
   spanning six Angstrom axially, so every position sits at the same height and
   the distance gate separates nothing. On 7CBC and 7JH5 the same tilt is one
   to three degrees.

So the rule has an unstated precondition, and the useful outcome of this run is
that it is now stated and enforced: the rule declines unless the latch runs
along the bundle. Two bugs were fixed on the way, both of which had produced
confident wrong answers rather than errors: the latch identifier picked a
three-residue helical turn on 9DCX because a fragment has no packing partners,
and the verdict did not check whether the test structure was a homologue.

BOTTOM LINE. The generalisation claim is NOT established. One test was a
homologue and the other was out of domain. What is established is the rule's
domain of validity, and a sequence-sensitivity number that matters more for the
project than either test would have: an arbitrary latch has roughly a
one-in-seven chance of arriving with a usable reporter site already in place.

PRE-REGISTERED PREDICTION, stated before any numbers were computed:

    If the rule generalises, the surviving positions will fall in the FIRST TWO
    TURNS of the latch helix, within about eight residues of the anchored end.

WHAT IS AND IS NOT CARRIED OVER
-------------------------------
Carried over unchanged, same numeric thresholds:
    occlusion   side-chain burial at the cage-latch interface >= 0.35
    volume      buried fraction of the dye fits the native side chain
    baseline    caged peak current >= 0.1 nA
Carried over structurally:
    binder      exclude anything past the latch helix
NOT carried over, and said so:
    tuning site 7CBC position 255 is the engineered V255S handle on the opening
                free energy; the analogue cannot be identified blind, so this
                gate is dropped rather than guessed.

Topology is derived from coordinates. Nothing is assumed from 7CBC.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from structure import analyze                # noqa: E402
from echem import swv                        # noqa: E402
from mb_placement_design import SIDECHAIN_VOL, MB_VOLUME_A3  # noqa: E402

PDB_DIR = HERE / "structure" / "pdb"

# Identical to the 7CBC run
CELL = dict(area_cm2=0.03, gamma_mol_cm2=2.0e-12, e_step=0.002, freq_hz=87.0)
K0_CONTACT, BETA, D_CONTACT = 3.0e2, 1.0, 5.0
BURIAL_MIN, BASELINE_FLOOR_NA, STANDOFF = 0.35, 0.1, 5.0
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90}
PROBE = 1.4


def fibonacci_sphere(n):
    i = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi),
                 np.sin(theta) * np.sin(phi), np.cos(phi)]


def load_atoms(path, chain="A"):
    out = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != chain or line[16] not in (" ", "A"):
            continue
        resn = line[17:20].strip()
        if resn in ("HOH", "EOH"):
            continue
        el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        out.append((int(line[22:26]), resn, line[12:16].strip(), el,
                    np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])])))
    return out


def sasa(atoms, subset, universe, n_points=120):
    sp = fibonacci_sphere(n_points)
    xyz = np.array([atoms[i][4] for i in universe])
    rad = np.array([VDW.get(atoms[i][3], 1.7) for i in universe]) + PROBE
    where = {g: k for k, g in enumerate(universe)}
    out = {}
    for gi in subset:
        k = where[gi]
        pts = xyz[k] + rad[k] * sp
        d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2)
        d[:, k] = 1e9
        out[gi] = float((d >= rad[None, :]).all(axis=1).mean()
                        * 4.0 * np.pi * rad[k] ** 2)
    return out


def helices_from_header(path, chain="A"):
    out = []
    for line in open(path):
        if line.startswith("HELIX") and line[19] == chain:
            out.append((int(line[21:25]), int(line[33:37])))
    return sorted(out)


def identify_latch(ca, helices, min_length=12):
    """Find the latch from topology alone.

    The latch is the helix that (a) is long enough to be a helix rather than a
    turn, (b) terminates the chain, so releasing it does not break the bundle,
    and (c) packs against the fewest DISTINCT partner helices, so it sits in a
    groove rather than inside the core. All three are statements about this
    structure; none is carried over from 7CBC.

    The length floor is not cosmetic. Without it a three-residue helical turn
    scores zero packing partners and wins, which is what happened on the first
    run against 9DCX before this was fixed.
    """
    segs = [np.array([ca[r] for r in range(a, b + 1) if r in ca]) for a, b in helices]
    n = len(helices)
    partners = []
    for i in range(n):
        p = 0
        for j in range(n):
            if i == j or len(segs[i]) == 0 or len(segs[j]) == 0:
                continue
            d = np.linalg.norm(segs[i][:, None, :] - segs[j][None, :, :], axis=-1)
            if (d < 10.0).sum() > 20:
                p += 1
        partners.append(p)
    long_enough = [i for i, (a, b) in enumerate(helices)
                   if b - a + 1 >= min_length]
    if not long_enough:
        raise ValueError("no helix long enough to be a latch")
    fewest = min(partners[i] for i in long_enough)
    candidates = [i for i in long_enough if partners[i] == fewest]
    latch = max(candidates)          # the terminal one among them
    return latch, partners


AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}


def chain_sequence(path, chain="A"):
    d = {}
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA" \
                and line[21] == chain and line[16] in (" ", "A"):
            d[int(line[22:26])] = AA3.get(line[17:20].strip(), "X")
    return d


def best_ungapped_identity(sa, sb, min_overlap=50):
    """Highest identity over any ungapped offset. Crude, but enough to tell a
    homologue from an unrelated design, which is the only question here."""
    best = 0.0
    for off in range(-len(sb) + min_overlap, len(sa) - min_overlap):
        n = m = 0
        for i, ch in enumerate(sa):
            j = i - off
            if 0 <= j < len(sb):
                m += 1
                n += ch == sb[j]
        if m >= min_overlap and n / m > best:
            best = n / m
    return best


def independence_check(path, chain="A"):
    """How much of this test is actually new information.

    A rule devised on 7CBC and re-tested on a close homologue of 7CBC is not
    being tested; it is being re-run. This prints the identity so the reader
    can discount the verdict accordingly, rather than the verdict quietly
    overstating itself.
    """
    ref = chain_sequence(PDB_DIR / "7cbc.pdb")
    sa = "".join(chain_sequence(path, chain)[r]
                 for r in sorted(chain_sequence(path, chain)))
    sb = "".join(ref[r] for r in sorted(ref))
    ident = best_ungapped_identity(sa, sb)
    if ident > 0.60:
        verdict = ("HOMOLOGUE. This is a consistency check, NOT a test of "
                   "generalisation.")
    elif ident > 0.35:
        verdict = "PARTIAL homology. Treat the result as weak evidence."
    else:
        verdict = "INDEPENDENT design lineage. A genuine test."
    return ident, verdict


def run(pdb_id="7JH5", chain="A"):
    path = analyze.fetch(pdb_id, PDB_DIR)
    ca = analyze.load_ca(path, chain, hetatm=True)
    atoms = load_atoms(path, chain)
    helices = helices_from_header(path, chain)

    print("=" * 78)
    print(f"BLIND TEST OF THE PLACEMENT RULE ON {pdb_id}")
    print("=" * 78)
    ident, indep = independence_check(path, chain)
    print(f"\n  Sequence identity to 7CBC, the structure the rule was built on:"
          f" {ident*100:.0f}%")
    print(f"  {indep}")
    print(f"\n  {len(ca)} resolved residues, {len(helices)} helices from the header:")
    for a, b in helices:
        print(f"     {a:4d} - {b:4d}   length {b - a + 1}")

    idx, partners = identify_latch(ca, helices)
    la, lb = helices[idx]
    print(f"\n  distinct packing partners per helix: {partners}")
    print(f"  latch identified from topology alone: helix {idx + 1}, residues {la}-{lb}")

    # bundle axis from everything except the latch, heights zeroed at the low end
    cage = np.array([ca[r] for r in sorted(ca) if r < la])
    centre = cage.mean(axis=0)
    axis = np.linalg.svd(cage - centre)[2][0]
    height = {r: float((ca[r] - centre) @ axis) for r in ca}
    junction = height[la]
    span = max(height.values()) - min(height.values())
    # APPLICABILITY PRECONDITION, added after 9DCX exposed its absence.
    # The distance gate presumes a latch that runs ALONG the bundle away from
    # the anchored end. If the latch lies across the end instead, every position
    # sits at the same height and the gate carries no information, so the rule
    # must decline rather than return an answer it cannot support.
    lat_xyz = np.array([ca[r] for r in range(la, lb + 1) if r in ca])
    lat_axis = np.linalg.svd(lat_xyz - lat_xyz.mean(axis=0))[2][0]
    tilt = float(np.degrees(np.arccos(abs(float(lat_axis @ axis)))))
    hs = [height[r] for r in range(la, lb + 1) if r in ca]
    axial_span = max(hs) - min(hs)
    print(f"  bundle length along its axis: {span:.0f} A")
    print(f"  latch tilt relative to the bundle axis: {tilt:.0f} deg")
    print(f"  latch axial span: {axial_span:.0f} A")
    if tilt > 45.0 or axial_span < 20.0:
        print(f"""
  OUT OF DOMAIN. The rule DECLINES to score this structure.

  The distance gate presumes a latch running along the bundle away from the
  anchored end, so that different positions sit at different heights above the
  electrode. Here the latch lies across the end of the bundle at {tilt:.0f} degrees and
  spans only {axial_span:.0f} A axially, so every position is at effectively the same
  height and the gate separates nothing.

  This is not the rule failing. It is the rule's domain of validity, which
  until this run was implicit and therefore untested. On 7CBC and 7JH5 the
  latch tilt is a few degrees; here it is nearly perpendicular, which is a
  different architecture wearing the same functional name.
""")
        return
    print(f"  latch runs from {abs(height[la]-min(height.values())):.0f} to "
          f"{abs(height[lb]-min(height.values())):.0f} A, so it hangs off one end")
    print(f"  the electrode is bonded at the latch-junction end, residue {la}")

    # burial of every latch residue: caged (whole chain) vs released (latch alone)
    latch_idx = [i for i, a in enumerate(atoms) if la <= a[0] <= lb]
    universe = list(range(len(atoms)))
    caged = sasa(atoms, latch_idx, universe)
    released = sasa(atoms, latch_idx, latch_idx)
    per = {}
    for i in latch_idx:
        r = atoms[i][0]
        per.setdefault(r, [atoms[i][1], 0.0, 0.0])
        if atoms[i][2] not in ("N", "CA", "C", "O"):
            per[r][1] += caged[i]
            per[r][2] += released[i]

    print("\n" + "=" * 78)
    print("THE SAME FOUR GATES, SAME THRESHOLDS, NOTHING RETUNED")
    print("=" * 78)
    print("\n  res  aa   burial  MB deficit  caged d  baseline   occl vol base")
    print("  " + "-" * 66)
    survivors = []
    for r in sorted(per):
        name, f, a_free = per[r]
        if a_free < 1.0:
            continue
        burial = (a_free - f) / a_free
        if burial < 0.30:
            continue
        native = SIDECHAIN_VOL.get(name, 80)
        deficit = burial * MB_VOLUME_A3 - native
        d = abs(junction - height[r]) + STANDOFF
        k0 = K0_CONTACT * np.exp(-BETA * (d - D_CONTACT))
        base = abs(swv.swv_scan(k0, **CELL)["peak_current"]) * 1e9
        g = (burial >= BURIAL_MIN, deficit <= 0.0, base >= BASELINE_FLOOR_NA)
        marks = "".join("  Y " if x else "  . " for x in g)
        star = "  <== PASSES" if all(g) else ""
        if all(g):
            survivors.append((r, name, burial, deficit, d, base))
        print(f"  {r:4d} {name} {burial*100:6.1f}%  {deficit:+8.0f} A3 "
              f"{d:7.1f} A {base:8.2f} nA {marks}{star}")

    print("\n" + "=" * 78)
    print("VERDICT AGAINST THE PRE-REGISTERED PREDICTION")
    print("=" * 78)
    if not survivors:
        print("\n  No position survives. The rule does not transfer.")
        return
    offsets = [r - la for r, *_ in survivors]
    names = ", ".join(f"{n}{r}" for r, n, *_ in survivors)
    print(f"\n  Survivors: {names}")
    print(f"  Offsets from the latch start: {offsets}")
    print(f"  Prediction was: within about 8 residues, the first two turns.")
    ok = all(o <= 8 for o in offsets)
    print(f"\n  RESULT: {'CONFIRMED' if ok else 'NOT CONFIRMED'}")
    if ok:
        print("""
  Every surviving position on a structure the rule had never seen falls in the
  first two turns of the latch helix, which is where it fell on 7CBC. The rule
  is a statement about cage-and-latch architecture, not a fit to one example.

  This does not make the thresholds right. It makes their SHAPE right: burial
  against volume against distance, resolved at the anchored end of the latch.
  The numeric cut-offs still rest on uncalibrated electrochemistry.""")
    else:
        print("""
  The survivors do not sit where the rule predicted. Either the rule was fitted
  to 7CBC, or this construct differs in a way the rule does not capture. Either
  way the R249 result on 7CBC should be treated as a single observation rather
  than an instance of a general rule.""")

    tail = max(ca) - lb
    print(f"""
  Two gates could not be applied and are reported as such rather than fudged.
  The binder-domain exclusion has {'nothing to exclude here: this construct ends only ' + str(tail) + os.linesep + '  residues after the latch helix and carries no grafted binder' if tail < 20 else 'a ' + str(tail) + '-residue tail past the latch to exclude'}. The
  opening-free-energy tuning site cannot be identified without knowing which
  interface mutation the authors used, so it was dropped rather than guessed.
""")


def sequence_sensitivity():
    """How much of the answer is geometry, and how much is the sequence?

    The backbone fixes which positions have the right burial and height. The
    volume gate then asks whether the residue that happens to sit there is
    large enough to vacate room for the dye. That second question is pure
    sequence, so the answer moves when the latch does.
    """
    from mb_placement_design import (BURIAL, load_scaffold, caged_distance,
                                     peak_na)
    ident, height = load_scaffold()
    freq = {"A": .082, "R": .055, "N": .041, "D": .054, "C": .014, "Q": .039,
            "E": .067, "G": .071, "H": .023, "I": .059, "L": .097, "K": .058,
            "M": .024, "F": .039, "P": .047, "S": .066, "T": .054, "W": .011,
            "Y": .029, "V": .069}
    one = {v: k for k, v in AA3.items() if v != "M" or k == "MET"}
    print("=" * 78)
    print("HOW MUCH OF THE ANSWER IS GEOMETRY, AND HOW MUCH IS SEQUENCE")
    print("=" * 78)
    print("\n  pos  burial  baseline   volume needed   residues large enough")
    print("  " + "-" * 74)
    rows = []
    for r, b in sorted(BURIAL.items()):
        if b < 0.35:
            continue
        d = caged_distance(height, r)
        base = peak_na(d)
        if base < 0.1:
            continue
        need = b * MB_VOLUME_A3
        ok = [k for k, v in sorted(SIDECHAIN_VOL.items(), key=lambda kv: -kv[1])
              if v >= need]
        p_ok = sum(freq[AA3[k]] for k in ok)
        rows.append((r, need, ok, p_ok))
        print(f"  {r:4d} {b*100:6.1f}% {base:8.2f} nA  {need:8.0f} A3   "
              f"{', '.join(ok) if ok else 'NONE, no natural residue is big enough'}")
    print("\n  chance a random residue at each position is large enough:")
    for r, need, ok, p_ok in rows:
        print(f"    position {r}: {p_ok*100:5.1f}%")
    all_p = 1.0 - np.prod([1.0 - p for _, _, _, p in rows])
    no_pro = [row for row in rows if ident.get(row[0]) != "PRO"]
    p_no_pro = 1.0 - np.prod([1.0 - p for _, _, _, p in no_pro]) if no_pro else 0.0
    print(f"""
  Chance that at least one usable position qualifies: {all_p*100:.0f}%
  Excluding the helix-capping proline, which cannot be mutated: {p_no_pro*100:.0f}%

  So the R249 result is not a geometric inevitability. The backbone picks the
  window; the sequence decides whether anything in that window can physically
  host a 225 cubic Angstrom dye, and only the four largest residues can at 59%
  burial.

  The number to give Monod Bio is the second one. An arbitrary latch has
  roughly a one-in-seven chance of arriving with a usable site already in
  place. If it does not, the fix is a conservative point substitution to put a
  large residue at the right position in the first two turns, not a search
  elsewhere on the latch, because everywhere else fails on burial, baseline or
  volume regardless of sequence.
""")


if __name__ == "__main__":
    targets = sys.argv[1:] or ["7JH5", "9DCX"]
    for t in targets:
        run(t)
        print("=" * 78 + "\n")
    sequence_sensitivity()
    print("=" * 78)
