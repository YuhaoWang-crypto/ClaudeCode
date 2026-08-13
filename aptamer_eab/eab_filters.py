"""EAB-specific structural filters on top of the generic aptamer-design kernel.

The generic kernel scores "is this a well-folded aptamer?".  An electrochemical
aptamer-based (EAB / E-AB) sensor needs something different: a sequence that is
*conformationally switchable*, because the entire signal comes from a
binding-induced change in how fast the methylene-blue (MB) tag exchanges
electrons with the gold surface.

=== HOW THE SCORE WAS ARRIVED AT, AND WHY YOU SHOULD DISTRUST IT ===

A first version of this score used the raw probability of the MFE structure and
the raw MFE.  Both are strongly length-dependent, so on an 80-nt construct they
collapse - and that version ranked the published, working NPY aptamer 4.31
(230 % signal change) BELOW its own scramble control.  A score that is
anti-correlated with the only ground truth available is worse than no score.

The current version therefore uses only LENGTH-NORMALISED descriptors, and is
checked against the two published E-AB outcomes for this exact target:

    4.31 full 80-nt construct  -> works, 230 % signal change   (must rank HIGH)
    4.31 bare 40-nt core       -> "no appreciable signal"      (must rank LOW)

`anchor_check()` runs that test.  Be explicit about what this means: the
windows below were adjusted AFTER looking at those two sequences, so this is a
2-point ANCHOR, not a validation.  Two points cannot validate a five-term
heuristic.  Treat the score as a way to order a synthesis list, never as a
prediction of affinity or of signal change.
"""

import sys

sys.path.insert(0, "/root/.claude/skills/synced/aptamer-design")
from kernel import (fold_one, gc_frac, max_homopolymer_run)  # noqa: E402

# --------------------------------------------------------------- windows ----
# [DESIGN] (lo, hi, ideal) windows on length-normalised descriptors.
SWITCH_SPEC = {
    # Base pairs joining the 5' anchor end to the 3' MB end in the MFE fold -
    # the element that couples a binding event at the core to MB motion.
    #
    # SCOPE, corrected after the orthogonal-validation round: this term
    # separates FRAMED from UNFRAMED architecture (4.31 full = 2, its inactive
    # bare core = 0), which is what the anchor check needs. It does NOT
    # discriminate within the framed libraries: the arms are identical in every
    # framed candidate, so the composition-matched scramble also scores 2, and
    # the term is near-constant across libraries A-D. The within-set
    # discrimination actually comes from mfe_per_nt and norm_ens_div, which do
    # separate 4.31 (-0.150 / 0.169) from its scramble (-0.076 / 0.344).
    # Read the 0.30 weight as "enforce the architecture", not "rank the cores".
    "terminal_bp": (0.001, 10.0, 3.0),
    # MFE per nucleotide (kcal/mol/nt, DNA Mathews2004, 37 C).  Length-robust
    # stand-in for "folded but not welded shut".  4.31 sits at -0.150; its
    # scramble at -0.076 (too floppy).
    "mfe_per_nt": (-0.32, -0.06, -0.16),
    # Ensemble diversity / length.  LOW = one dominant fold, HIGH = the
    # molecule is a disordered blur.  4.31 = 0.169, its scramble = 0.344.
    "norm_ens_div": (0.02, 0.42, 0.15),
    # Fraction of nucleotides paired - enough structure to present a surface.
    "paired_frac": (0.25, 0.75, 0.52),
}

# [DESIGN] How many of the last 5 nt should be unpaired.  A dangling 3' tail
# lets MB sample the surface; 4.31 has all 5 free.  Scored separately because
# the direction is monotonic, not a window.
TAIL_FREE_IDEAL = 4          # >= this many of the last 5 nt unpaired = full marks
GC_WINDOW = (0.35, 0.65)
MAX_HOMOPOLYMER = 4

WEIGHTS = {"term": 0.30, "stab": 0.25, "order": 0.20, "pair": 0.15,
           "tail": 0.10}


def _bell(value, lo, hi, ideal):
    """[DESIGN] 1.0 at `ideal`, tapering linearly to 0.0 at `lo`/`hi`, 0.0
    outside.  A metric is rewarded for being in a WINDOW, not for being
    extreme - 'most stable fold wins' ranking selects exactly the rigid
    sequences that cannot switch."""
    if value <= lo or value >= hi:
        return 0.0
    if value == ideal:
        return 1.0
    if value < ideal:
        return (value - lo) / (ideal - lo)
    return (hi - value) / (hi - ideal)


def pair_table(struct):
    """[DERIVED] Dot-bracket -> 0-based partner index list (-1 = unpaired)."""
    partner = [-1] * len(struct)
    stack = []
    for i, ch in enumerate(struct):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            partner[i], partner[j] = j, i
    return partner


def terminal_stem_bp(struct, k=10):
    """[PREDICTED] Base pairs in the MFE fold joining the first `k` nucleotides
    to the last `k` - the 5'-thiol anchor end to the 3'-MB end."""
    n = len(struct)
    partner = pair_table(struct)
    k = min(k, n // 2)
    return sum(1 for i in range(k) if partner[i] >= n - k)


def three_prime_free(struct, tail=5):
    """[PREDICTED] How many of the last `tail` nucleotides are UNPAIRED."""
    return sum(1 for ch in struct[-tail:] if ch == ".")


def arm_pairing(struct, arm=20):
    """[PREDICTED] For a FRAMED construct (20-nt arm + core + 20-nt arm), split
    the base pairs by segment: (arm_arm_bp, arm_core_bp).

    Returns (None, None) for anything too short to be framed.

    Why these are reported but NOT scored: the two fixed arms are IDENTICAL in
    every framed candidate, so both numbers are near-constant across the framed
    libraries and neither can discriminate a good candidate from a bad one -
    the composition-matched scramble of 4.31 has arm_arm_bp = 9 and
    arm_core_bp = 4, essentially the same as real 4.31 (9 and 5).

    They earn their place as a DESIGN CONSTRAINT instead: doped mutagenesis of
    the core can destroy the 5'arm<->core duplex, and that duplex is where the
    orthogonal contact-probability run puts roughly a third of the predicted
    protein-DNA interface (strongest 5'-arm contacts T14/C15/A16 - two of which
    are exactly the arm positions this fold pairs into the core). Candidates
    that break it are throwing away the part of 4.31 two independent methods
    agree is load-bearing.
    """
    n = len(struct)
    if n < 3 * arm:
        return None, None
    pt = pair_table(struct)
    aa = sum(1 for i, j in enumerate(pt) if j > i and i < arm and j >= n - arm)
    ac = sum(1 for i, j in enumerate(pt)
             if j > i and i < arm and arm <= j < n - arm)
    return aa, ac


def self_dimer_dg(seq):
    """[PREDICTED, informational only] Best intermolecular duplex energy with a
    copy of itself (ViennaRNA duplexfold, DNA parameters).

    NOT used in the score.  The working 4.31 construct self-dimerises harder
    (-25.8 kcal/mol) than its scramble (-12.2), so penalising it would have
    demoted the one sequence known to work.  On a diluted, mercaptohexanol-
    backfilled monolayer, intermolecular duplex formation is largely suppressed
    anyway.  Reported so a very extreme value can be spotted by eye.
    """
    import RNA
    RNA.params_load_DNA_Mathews2004()
    return round(RNA.duplexfold(seq, seq).energy, 2)


def eab_metrics(seq, temperature=37.0):
    """[PREDICTED] Full EAB metric row for one ssDNA candidate."""
    f = fold_one(seq, is_rna=False, temperature=temperature)
    struct = f["struct"]
    m = {
        "seq": seq,
        "len": len(seq),
        "gc": round(gc_frac(seq), 3),
        "mfe": f["mfe"],
        "p_mfe": f["p_mfe"],
        "ens_div": f["ens_div"],
        "paired_frac": f["paired_frac"],
        "struct": struct,
        "terminal_bp": terminal_stem_bp(struct),
        "tail3_free": three_prime_free(struct),
        "self_dimer_dg": self_dimer_dg(seq),
        "max_run": max_homopolymer_run(seq),
    }
    m["mfe_per_nt"] = round(m["mfe"] / m["len"], 4)
    m["norm_ens_div"] = round(m["ens_div"] / m["len"], 4)
    aa, ac = arm_pairing(struct)
    m["arm_arm_bp"] = aa
    m["arm_core_bp"] = ac
    # [DESIGN] 4.31 pairs 5 arm positions (14,15,18,19,20) into its core.
    # Keep >=4 to retain that duplex through core mutagenesis.
    m["keeps_arm_core_duplex"] = None if ac is None else bool(ac >= 4)
    return m


def eab_score(m):
    """[DESIGN] Composite 0-100 switchability score.  See the module docstring
    for why this is an anchored heuristic and not a validated predictor."""
    s = {
        "term": _bell(m["terminal_bp"], *SWITCH_SPEC["terminal_bp"]),
        "stab": _bell(m["mfe_per_nt"], *SWITCH_SPEC["mfe_per_nt"]),
        "order": _bell(m["norm_ens_div"], *SWITCH_SPEC["norm_ens_div"]),
        "pair": _bell(m["paired_frac"], *SWITCH_SPEC["paired_frac"]),
        "tail": min(1.0, m["tail3_free"] / TAIL_FREE_IDEAL),
    }
    raw = sum(WEIGHTS[k] * v for k, v in s.items())
    if not (GC_WINDOW[0] <= m["gc"] <= GC_WINDOW[1]):
        raw -= 0.15
    if m["max_run"] > MAX_HOMOPOLYMER:
        raw -= min(0.10, 0.033 * (m["max_run"] - MAX_HOMOPOLYMER))
    return round(100 * max(0.0, min(1.0, raw)), 2)


def score_pool(seqs, temperature=37.0, extra=None):
    """[PREDICTED] Metric + score every sequence; returns a DataFrame sorted by
    eab_score.  `extra` is an optional dict-of-dicts keyed by sequence, merged
    into each row (carries `parent`/`arch`/`n_mut` provenance through)."""
    import pandas as pd
    rows = []
    for s in seqs:
        m = eab_metrics(s, temperature=temperature)
        m["eab_score"] = eab_score(m)
        if extra and s in extra:
            m.update(extra[s])
        rows.append(m)
    return (pd.DataFrame(rows)
            .sort_values("eab_score", ascending=False)
            .reset_index(drop=True))


# ------------------------------------------------------- G-quadruplexes -----
# ViennaRNA has NO G-quadruplex term in the DNA parameter set: a real G4 is
# scored as if it were an unstructured coil, so every G4 candidate is pushed to
# the bottom of `eab_score`.  That is a modelling blind spot, not a biological
# result - G4 is the scaffold behind the best-characterised E-AB sensor of all
# (the thrombin-binding aptamer), it is compact (15-30 nt, fast electron
# transfer, large relative signal change) and intrinsically nuclease-resistant.
# G4 candidates are therefore ranked by a sequence-propensity score instead.

def g4hunter_score(seq):
    """[PREDICTED] G4Hunter score (Bedrat, Lacroix & Mergny, NAR 2016).

    Each base scores by the length of the G-run (or C-run) it belongs to,
    capped at 4: G-runs positive, C-runs negative, A/T zero; the score is the
    mean over the sequence.  |score| > ~1.2 is the usual "likely G4" threshold
    for the G-rich strand.  Purely sequence-based - it says a G4 is plausible,
    not that one forms.  Confirm experimentally by CD (positive ~264 nm /
    negative ~240 nm for parallel G4) and a K+ - vs Li+ - dependent thermal melt.
    """
    import itertools
    vals = []
    for base, grp in itertools.groupby(seq.upper()):
        run = len(list(grp))
        n = min(run, 4)
        if base == "G":
            vals.extend([n] * run)
        elif base == "C":
            vals.extend([-n] * run)
        else:
            vals.extend([0] * run)
    return round(sum(vals) / len(vals), 3)


# [MEASURED] The thrombin-binding aptamer, the archetypal E-AB recognition
# element (Xiao, Lubin, Heeger & Plaxco, Angew Chem 2005) and a textbook
# antiparallel chair G4.  It anchors the G4 track: TBA's G4Hunter score is only
# ~1.13, BELOW the usual 1.2 "likely G4" threshold.  Maximising G4Hunter is
# therefore the WRONG objective for an aptamer - it drives the pool to
# near-homopolymeric G runs, which form intermolecular G-wires and aggregates,
# synthesise badly, and make hopeless monolayers.  The track targets the
# TBA-like regime instead.
TBA = "GGTTGGTGTGGTTGG"
G4HUNTER_WINDOW = (0.9, 2.4, 1.4)      # (lo, hi, ideal)
G4_MAX_GC = 0.72
G4_MAX_G_RUN = 4


def g4_tract_pattern(seq):
    """[DERIVED] (number of G-runs of length 2-4, longest G-run).  A canonical
    intramolecular G4 needs exactly four such tracts separated by real loops."""
    import itertools
    runs = [(b, len(list(g))) for b, g in itertools.groupby(seq.upper())]
    g_runs = [n for b, n in runs if b == "G"]
    return sum(1 for n in g_runs if 2 <= n <= 4), (max(g_runs) if g_runs else 0)


def g4_quality(seq):
    """[DESIGN] 0-100 G4 plausibility for an APTAMER (not for a genome scan).

    Rewards being in the TBA-like G4Hunter window and having exactly four
    well-formed G-tracts; hard-penalises the G-homopolymer corner that raw
    G4Hunter maximisation runs into.
    """
    n_tracts, longest = g4_tract_pattern(seq)
    s_hunter = _bell(g4hunter_score(seq), *G4HUNTER_WINDOW)
    s_tracts = 1.0 if n_tracts == 4 else (0.4 if n_tracts in (3, 5) else 0.0)
    raw = 0.55 * s_hunter + 0.45 * s_tracts
    if gc_frac(seq) > G4_MAX_GC:
        raw -= 0.35
    if longest > G4_MAX_G_RUN:
        raw -= 0.35
    return round(100 * max(0.0, min(1.0, raw)), 2)


def g4_metrics(seq, temperature=37.0):
    """[PREDICTED] Metrics for a G4 candidate: the same structural columns (so
    the tables line up) plus G4Hunter and the aptamer-oriented g4_quality,
    which is what they are RANKED by.  `eab_score` is reported for these but is
    NOT meaningful - ViennaRNA cannot see a G-quadruplex at all."""
    m = eab_metrics(seq, temperature=temperature)
    m["g4hunter"] = g4hunter_score(seq)
    n_tracts, longest = g4_tract_pattern(seq)
    m["g4_tracts"] = n_tracts
    m["longest_g_run"] = longest
    m["g4_quality"] = g4_quality(seq)
    m["eab_score"] = eab_score(m)
    m["eab_score_meaningful"] = False
    return m


# ----------------------------------------------------------- anchor check ---
def anchor_check(verbose=True):
    """[DERIVED] Does the score reproduce the two published E-AB outcomes?

    Returns (passed: bool, rows: list).  Called by run_all so a broken scoring
    function fails loudly instead of silently ordering the wrong oligos.
    """
    from . import known
    from .build_library import frame, scramble

    cases = [
        ("4.31 full 80-nt  [MEASURED: works, 230% SC]", known.APT431, "HIGH"),
        ("4.31 bare core   [MEASURED: no signal]", known.APT431_CORE, "LOW"),
        ("scramble control [expected: LOW]", frame(scramble(known.APT431_CORE)),
         "LOW"),
    ]
    rows = []
    for label, seq, want in cases:
        m = eab_metrics(seq)
        m["eab_score"] = eab_score(m)
        rows.append({"case": label, "want": want, "score": m["eab_score"],
                     "len": m["len"], "terminal_bp": m["terminal_bp"],
                     "mfe_per_nt": m["mfe_per_nt"],
                     "norm_ens_div": m["norm_ens_div"]})
    active = rows[0]["score"]
    passed = all(active > r["score"] + 10 for r in rows[1:])
    if verbose:
        for r in rows:
            print(f"    {r['case']:<45} score={r['score']:>6.1f}  "
                  f"termBP={r['terminal_bp']}  mfe/nt={r['mfe_per_nt']:+.3f}  "
                  f"ndiv={r['norm_ens_div']:.3f}")
        print(f"    ANCHOR CHECK: {'PASS' if passed else 'FAIL'} - the working "
              f"construct must out-score both negatives by >10 points.")
        print("    n=2 known outcomes and the windows were set after seeing "
              "them:\n    this is an anchor, NOT a validation.")
    return passed, rows


# ------------------------------------------------------------- constructs ---
def eab_construct(seq, label, mb_position="3prime"):
    """[DESIGN] The orderable EAB construct for one sequence.

    Defaults follow the published NPY E-AB sensor (5'-thiol anchor, 3'-MB,
    terminally labelled), which reached 230 % signal change at 100 Hz SWV.
    """
    if mb_position == "3prime":
        s = f"5'-HS-(CH2)6-{seq}-MB-3'"
    else:
        s = f"5'-HS-(CH2)6-{seq}  [MB internal, see mb_position]-3'"
    return {
        "id": label,
        "construct": s,
        "length_nt": len(seq),
        "anchor": "5' hexanethiol, gold SAM, 6-mercapto-1-hexanol backfill",
        "reporter": "methylene blue (C7 amino linker) at the 3' terminus",
        "readout": "square-wave voltammetry; scan frequency 60-500 Hz, "
                   "pick the gain-optimal frequency empirically",
    }
