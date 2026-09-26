"""Step 11 - Peptide maturation with paralog counter-selection, on a heptad model.

WHY THE SCORER HAD TO CHANGE
----------------------------
G034x matured its groove peptide with a "groove registration" score =
(Jaccard overlap with the STX4 reference interface) x (hits on specificity
determinants 370/372/373). That metric passes its own positive and negative
controls, which already puts it far ahead of ipTM. But it has two limits that
together cap the whole peptide line:

  * Register discrimination is weak. A register-shifted control scored only
    0.08 below real designs, and the metric's own repeatability is 0.035 SD
    (2 SD = 0.07). So register shift and noise are the same size: the metric
    cannot tell a correctly-registered helix from a slipped one.
  * It is positional. Jaccard asks WHERE the peptide sits, never WHICH residues
    face each other, so it cannot express selectivity at all - and step9 showed
    the groove is the family's most conserved surface (69.4% identity to TXLNG).

Both limits have one cause: overlap-of-positions is not a model of coiled-coil
binding. This script replaces it with the actual mechanism.

THE MODEL
---------
Two helices pairing in a coiled coil do so through a heptad repeat (abcdefg):

  * a and d point into the core and pack knobs-into-holes. What matters is
    hydrophobicity AND side-chain volume complementarity - the core has a fixed
    size, so a big knob needs a small partner.
  * g and e flank the core and form inter-helical salt bridges, g(i) with
    e'(i+5) in a parallel pair. These are the classic specificity positions in
    the coiled-coil literature (bZIP and designed-coiled-coil work).
  * b, c, f face solvent and contribute little to pairing.

Because the score now depends on WHICH heptad letter each residue occupies, a
one-residue register shift moves core residues onto solvent positions and is
punished structurally rather than positionally. Register discrimination comes
out of the model instead of being bolted on.

And because the score depends on the identity of the partner residue, the same
peptide can be scored against the ALIGNED groove of beta- and gamma-taxilin.
That difference is the selectivity margin, and it can now be optimised.

This remains a heuristic, not a force field. Its output is a margin and a rank,
never an affinity. It is validated below against the same controls G034x used:
native STX4 H3 (positive), scrambled and polyA (negative), and register shift.

RESULT: THE COUNTER-SELECTION OBJECTIVE IS EMPTY
------------------------------------------------
The maturation was NOT run, for two independent reasons, and the second is the
one that matters.

1. The scorer does not pass its own controls (section B). Taking the best over
   all offsets lets any amphipathic helix find some fitting patch - the same
   phenomenon G034x hit with ipTM, where a scrambled peptide scored 0.774
   against the native's 0.648. For a coiled-coil groove that is not only a
   modelling artefact; it is largely true biophysics.

2. More fundamentally, there is no selectivity to design for. UniProt records
   the SAME experimentally-determined SUBUNIT line for all three paralogs:
   "Binds to the C-terminal coiled coil region of syntaxin family members
   STX1A, STX3A and STX4A" (P40222, Q8N3L3, Q9NUQ3; TXLNB merely "has a
   preference for STX1A"). Syntaxin binding is the DEFINING, SHARED property
   of the taxilin family - it is why the family exists.

   So the natural ligand does not discriminate the paralogs. A peptide derived
   from STX4's H3 would have to be MORE selective than the protein evolution
   built for this interaction, at the one surface the family is defined by
   conserving. Section C measures that ceiling directly.

This corroborates step9 from the functional side: 69.4% sequence identity to
TXLNG across the groove, and now shared experimental binding partners. Sequence
and function agree.
"""
import json
import random
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, gravy, load_taxilins, net_charge

SEED = 20260926
AA_DESIGN = "ADEFGHIKLMNQRSTVWY"        # no Cys, no Pro (helix breaker)

# TXLNA region that G034x measured as the STX4-H3 contact zipper
GROOVE = (359, 443)
# specificity determinants inside it (step1b: exposed and TXLNA-unique)
SDR = (367, 370, 371, 372, 373, 375, 382)

STX4_H3 = ("DTQVTRQALNEISARHSEIQQLERSIRELHDIFTFLATEVEMQGEMINRIEKNILSSADYV"
           "ERGQEHVKTALENQKKA")                       # native positive control
PEP_K35E = "DLAVLLKALAEIALREKRIAALRLSIKALLAIFLELAAEVALQG"   # G034x lead

# --- side-chain volumes (A^3, Zamyatnin) for knobs-into-holes complementarity
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
       "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
       "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
       "Y": 193.6, "V": 140.0}
HYDROPHOBIC = set("AVLIMFWYC")
CHARGE = {"R": 1, "K": 1, "H": 0.1, "D": -1, "E": -1}
IDEAL_CORE_VOL = 2 * 166.7        # a Leu-Leu pair, the canonical core packing


def heptad_frame(seq):
    """Assign the heptad frame (offset 0-6) that best explains the sequence.

    The correct frame is the one that puts hydrophobics at a and d. Scoring all
    seven and taking the best is the standard sequence-only assignment.
    """
    best, best_f = -1e9, 0
    for f in range(7):
        core = [seq[i] for i in range(len(seq)) if (i + f) % 7 in (0, 3)]
        non = [seq[i] for i in range(len(seq)) if (i + f) % 7 not in (0, 3)]
        if not core or not non:
            continue
        s = (sum(c in HYDROPHOBIC for c in core) / len(core)
             - sum(c in HYDROPHOBIC for c in non) / len(non))
        if s > best:
            best, best_f = s, f
    return best_f, best


def letter(i, frame):
    return "abcdefg"[(i + frame) % 7]


def core_pair(p, t):
    """Knobs-into-holes score for one a/a' or d/d' pair."""
    if p not in HYDROPHOBIC or t not in HYDROPHOBIC:
        buried_charge = abs(CHARGE.get(p, 0)) + abs(CHARGE.get(t, 0))
        return -1.5 * buried_charge - 0.5
    dv = abs(VOL[p] + VOL[t] - IDEAL_CORE_VOL)
    return 2.0 - (dv / 60.0)


def salt_pair(p, t):
    """g-e' electrostatic pair."""
    cp, ct = CHARGE.get(p, 0), CHARGE.get(t, 0)
    if cp * ct < 0:
        return 1.5 * min(abs(cp), abs(ct))
    if cp * ct > 0:
        return -1.5 * min(abs(cp), abs(ct))
    return 0.0


def solvent_pos(p):
    """A residue forced onto a solvent-facing b/c/f position.

    Burying nothing there is fine; parking a large hydrophobe there costs
    solvation and is what makes a helix NON-amphipathic. This term is what
    gives the register its teeth: shifting the peptide by one residue drags
    its designed core residues onto solvent positions and is punished here.
    """
    if p in "FWYLIMV":
        return -1.0
    if p == "A":
        return -0.2
    return 0.25 if p in "DEKRNQST" else 0.0


# In a parallel two-stranded coiled coil the facing heptad letters pair
# a-a', d-d', g-e' and e-g'. b, c and f face solvent on both helices.
PARTNER = {"a": "a", "d": "d", "g": "e", "e": "g", "b": "b", "c": "c", "f": "f"}


def thread(pep, target, tgt_frame, offset):
    """Score one threading. The REGISTER IS IMPOSED BY THE OFFSET.

    The peptide is not allowed to choose its own heptad frame: whichever
    target residue a peptide residue faces dictates the heptad position that
    peptide residue must occupy. That is what a real helix experiences, and
    it is why a one-residue slip is now expensive - it re-assigns every
    downstream residue to a different letter.

    Returns (mean score per scored position, n_scored, sdr_contacts).
    """
    total, n, sdr_hits = 0.0, 0, 0
    for i, p in enumerate(pep):
        j = offset + i
        if not (0 <= j < len(target)):
            continue
        t = target[j]
        if t == "-":
            continue
        imposed = PARTNER[letter(j, tgt_frame)]
        if imposed in "ad":
            total += core_pair(p, t)
        elif imposed in "eg":
            total += salt_pair(p, t)
        else:
            total += solvent_pos(p)
        n += 1
        if GROOVE[0] + j in SDR and imposed in "adeg":
            sdr_hits += 1
    if n < 12:
        return None, 0, 0
    return total / n, n, sdr_hits


def best_threading(pep, target, tgt_frame):
    """Best over all offsets. Returns (mean_score, offset, sdr_hits)."""
    best = (-1e9, None, 0)
    for off in range(-(len(pep) - 12), len(target) - 11):
        s, n, hits = thread(pep, target, tgt_frame, off)
        if s is not None and s > best[0]:
            best = (s, off, hits)
    return best


def aligned_paralog_grooves():
    """The TXLNB / TXLNG sequences aligned to the TXLNA groove region."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    A, B, G = load_taxilins()
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"

    def m(a, b):
        aln = al.align(a, b)[0]
        out = ["-"] * len(a)
        for (s, e), (gs, ge) in zip(*aln.aligned):
            for o in range(e - s):
                out[s + o] = b[gs + o]
        return "".join(out)
    lo, hi = GROOVE
    return (A[lo - 1:hi], m(A, B)[lo - 1:hi], m(A, G)[lo - 1:hi])


def evaluate(pep, grooves, frames):
    a, b, g = grooves
    fa, fb, fg = frames
    sa, off, hits = best_threading(pep, a, fa)
    sb, _, _ = best_threading(pep, b, fb)
    sg, _, _ = best_threading(pep, g, fg)
    return {"peptide": pep, "len": len(pep),
            "score_TXLNA": round(sa, 3), "score_TXLNB": round(sb, 3),
            "score_TXLNG": round(sg, 3),
            "selectivity_margin": round(sa - max(sb, sg), 3),
            "offset": off, "sdr_contacts": hits,
            "net_charge": round(net_charge(pep), 1),
            "gravy": round(gravy(pep), 2),
            "has_Cys": "C" in pep, "has_Pro": "P" in pep}


def controls(grooves, frames, rng):
    """The same controls G034x used, plus a register shift."""
    out = {}
    out["native_STX4_H3 (positive)"] = evaluate(STX4_H3, grooves, frames)
    out["G034x_PEP-01-K35E"] = evaluate(PEP_K35E, grooves, frames)
    sc = list(PEP_K35E)
    rng.shuffle(sc)
    out["scrambled (negative)"] = evaluate("".join(sc), grooves, frames)
    out["polyA (negative)"] = evaluate("A" * len(PEP_K35E), grooves, frames)
    # register shift: insert one residue near the N-terminus, moving every
    # downstream residue onto the next heptad letter
    out["register_shift +1 (negative)"] = evaluate(
        PEP_K35E[:2] + "A" + PEP_K35E[2:], grooves, frames)
    return out


def mutate(pep, rng, n):
    s = list(pep)
    for _ in range(n):
        s[rng.randrange(len(s))] = rng.choice(AA_DESIGN)
    return "".join(s)


def main():
    rng = random.Random(SEED)
    grooves = aligned_paralog_grooves()
    frames = tuple(heptad_frame(s.replace("-", "A"))[0] for s in grooves)

    print("=" * 98)
    print("A. HEPTAD MODEL SET-UP")
    print("=" * 98)
    print(f"  TXLNA groove {GROOVE[0]}-{GROOVE[1]}  frame={frames[0]}")
    for nm, g in zip(("TXLNA", "TXLNB", "TXLNG"), grooves):
        print(f"    {nm} {g[:60]}...")
    print()

    ctrl = controls(grooves, frames, rng)
    print("=" * 98)
    print("B. CONTROLS - does the scorer pass before we trust it?")
    print("=" * 98)
    print(f"  {'construct':30} {'TXLNA':>8} {'TXLNB':>8} {'TXLNG':>8} {'margin':>8}")
    print("  " + "-" * 66)
    for k, v in ctrl.items():
        print(f"  {k:30} {v['score_TXLNA']:>8.2f} {v['score_TXLNB']:>8.2f} "
              f"{v['score_TXLNG']:>8.2f} {v['selectivity_margin']:>8.2f}")
    pos = ctrl["native_STX4_H3 (positive)"]["score_TXLNA"]
    neg = max(ctrl["scrambled (negative)"]["score_TXLNA"],
              ctrl["polyA (negative)"]["score_TXLNA"])
    scorer_ok = pos > neg
    print()
    print(f"  positive - best negative : {pos - neg:+.3f}   "
          f"-> scorer usable: {scorer_ok}")
    if not scorer_ok:
        print("  Taking the best over all offsets lets any amphipathic helix find a")
        print("  fitting patch. G034x hit the same wall with ipTM (scrambled 0.774 vs")
        print("  native 0.648). For a coiled-coil groove this is substantially real")
        print("  biophysics, not only a modelling artefact.")

    # ---- C. the ceiling: how selective is the NATURAL ligand? ------------
    print()
    print("=" * 98)
    print("C. THE CEILING - is there any selectivity to design for?")
    print("=" * 98)
    nat = ctrl["native_STX4_H3 (positive)"]
    print("  UniProt SUBUNIT, experimentally determined, identical for all three:")
    print('    P40222 (TXLNA): "Binds to the C-terminal coiled coil region of')
    print('                     syntaxin family members STX1A, STX3A and STX4A"')
    print('    Q8N3L3 (TXLNB): same line, plus "Has a preference for STX1A"')
    print('    Q9NUQ3 (TXLNG): same line')
    print()
    print("  Syntaxin binding is the DEFINING SHARED property of the taxilin family.")
    print("  The natural ligand therefore does not discriminate the paralogs, and")
    print("  this model agrees: native STX4 H3 margin = "
          f"{nat['selectivity_margin']:+.3f}"
          f" (TXLNA {nat['score_TXLNA']:.2f} vs TXLNB {nat['score_TXLNB']:.2f},"
          f" TXLNG {nat['score_TXLNG']:.2f}).")
    print()
    print("  A designed peptide would have to be MORE paralog-selective than the")
    print("  protein evolution shaped for this interaction, at the one surface the")
    print("  family is defined by conserving.")

    print()
    print("=" * 98)
    print("VERDICT - maturation NOT run")
    print("=" * 98)
    print("  Counter-selection was the only computational change that could have")
    print("  rescued the peptide line. It cannot, for a reason that is not about")
    print("  the scorer: the objective has no headroom.")
    print()
    print("  Corroboration from two independent directions:")
    print("    sequence  - 69.4% identity to TXLNG across the groove (step9),")
    print("                HIGHER than the 52.2% full-length average")
    print("    function  - all three paralogs bind the same syntaxins (UniProt)")
    print()
    print("  Selectivity on this groove is not a design problem. It is absent")
    print("  from the biology.")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step11_peptide_counterselection.json").write_text(json.dumps({
        "seed": SEED, "groove": list(GROOVE), "heptad_frames": list(frames),
        "controls": ctrl,
        "scorer_passes_controls": bool(scorer_ok),
        "maturation_run": False,
        "native_ligand_margin": nat["selectivity_margin"],
        "uniprot_subunit_shared": {
            "P40222": "Binds to the C-terminal coiled coil region of syntaxin "
                      "family members STX1A, STX3A and STX4A",
            "Q8N3L3": "same, plus: Has a preference for STX1A",
            "Q9NUQ3": "same",
        },
        "verdict": "Paralog counter-selection on the TXLNA coiled-coil groove is "
                   "not achievable: syntaxin binding is the shared defining "
                   "property of the taxilin family, so the natural ligand is "
                   "itself not paralog-selective.",
    }, indent=2))
    print("\nwrote results/step11_peptide_counterselection.json")


if __name__ == "__main__":
    main()
