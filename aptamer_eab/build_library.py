"""Build the candidate ssDNA libraries for the NPY / PP EAB-sensor campaign.

Five sub-libraries, deliberately different in how much prior knowledge they
carry, because the two targets are in very different starting positions:

  A  NPY, knowledge-seeded   - doped maturation around aptamer 4.31, the ssDNA
                               binder that already drives a published NPY E-AB
                               sensor.  Highest prior, lowest risk.
  B  NPY, de novo structured - naive structured pool, for diversity / to escape
                               the 4.31 epitope (which may be the family-
                               conserved face).
  C  PP,  cross-seeded       - doped maturation around 4.31 at HIGHER doping.
                               Rationale: 4.31 is only 42-fold selective for
                               NPY over hPP, i.e. it measurably binds PP.  That
                               residual affinity is a real binding prior for a
                               PP binder; counter-selection against NPY is what
                               inverts the specificity.
  D  PP,  de novo structured - naive structured pool, independent RNG stream.
  E  degenerate oligo specs  - what you actually order for wet-lab SELEX.

Plus controls: a length- and composition-matched scramble of 4.31, and the
4.31 core alone (published to be INACTIVE in the E-AB format - a negative
control with a citation is worth more than a random one).
"""

import sys
import random

sys.path.insert(0, "/root/.claude/skills/synced/aptamer-design")
from kernel import (generate_aptamer_pool, doped_mutant_library,  # noqa: E402
                    levenshtein, gc_frac)

from . import known
from .eab_filters import score_pool


def frame(core, fwd=known.LIB_FWD, rev=known.LIB_REV):
    """[DESIGN] Put a selected core into the proven 80-nt E-AB architecture.

    The published NPY sensor works with the fixed primer arms in place and
    explicitly FAILS when they are removed, so the arms are treated as part of
    the construct, not as disposable amplification handles.
    """
    return fwd + core + rev


def scramble(seq, seed=1234):
    """[DESIGN] Composition-matched scramble - the standard SELEX / E-AB
    negative control (same length, same base counts, no selected structure)."""
    rng = random.Random(seed)
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


# ---------------------------------------------------------------- A: NPY ----
def library_A_npy_seeded(n=400, doping=0.15, seed=11):
    """Doped affinity-maturation pool around the 4.31 core, primers retained.

    doping=0.15 is the low end of the usual 15-30 % SELEX doping range: 4.31 is
    already a validated binder for this exact target, so the job is fine-tuning
    affinity and switch amplitude, not searching.  Expected ~6 mutations per
    40-nt core.
    """
    muts = doped_mutant_library([known.APT431_CORE], n_per_seed=n,
                                doping=doping, seed=seed, include_parent=True)
    extra = {frame(m["seq"]): {"core": m["seq"], "n_mut": m["n_mut"],
                               "parent": "APT431_core", "sublib": "A_NPY_seeded",
                               "arch": "doped_4.31"}
             for m in muts}
    return score_pool(list(extra.keys()), extra=extra)


# ---------------------------------------------------------------- B: NPY ----
def library_B_npy_denovo(n_per_arch=200, seed=42):
    """Naive structured pool (stem-loop / dual-hairpin / G-quadruplex /
    three-way-junction / random), framed into the 80-nt architecture."""
    pool = generate_aptamer_pool(is_rna=False, n_per_arch=n_per_arch,
                                 len_min=36, len_max=44, seed=seed)
    extra = {frame(c["seq"]): {"core": c["seq"], "n_mut": None,
                               "parent": None, "sublib": "B_NPY_denovo",
                               "arch": c["arch"]}
             for c in pool}
    return score_pool(list(extra.keys()), extra=extra)


# ----------------------------------------------------------------- C: PP ----
def library_C_pp_crossseeded(n=400, doping=0.30, seed=23):
    """Doped pool around 4.31 at HIGH doping, for re-targeting to PP.

    doping=0.30 (~12 mutations per 40-nt core) is deliberately at the top of
    the SELEX doping range: we want to keep 4.31's folding backbone - which is
    known to present a binding surface to this peptide family - while moving
    far enough in sequence space to flip the 42-fold preference from NPY to PP.
    """
    muts = doped_mutant_library([known.APT431_CORE], n_per_seed=n,
                                doping=doping, seed=seed, include_parent=False)
    extra = {frame(m["seq"]): {"core": m["seq"], "n_mut": m["n_mut"],
                               "parent": "APT431_core", "sublib": "C_PP_crossseeded",
                               "arch": "doped_4.31_high"}
             for m in muts}
    return score_pool(list(extra.keys()), extra=extra)


# ----------------------------------------------------------------- D: PP ----
def library_D_pp_denovo(n_per_arch=200, seed=7):
    pool = generate_aptamer_pool(is_rna=False, n_per_arch=n_per_arch,
                                 len_min=36, len_max=44, seed=seed)
    extra = {frame(c["seq"]): {"core": c["seq"], "n_mut": None,
                              "parent": None, "sublib": "D_PP_denovo",
                              "arch": c["arch"]}
             for c in pool}
    return score_pool(list(extra.keys()), extra=extra)


# ------------------------------------------------- short-construct track ----
def library_short(n_per_arch=200, seed=99, len_min=28, len_max=40):
    """Free-standing short switches (no primer arms).

    Higher risk for THIS target family - the published NPY sensor lost all
    signal when 4.31 was cut to its bare 40-nt core - but short constructs give
    faster electron transfer and larger relative signal change when they do
    work, so the track is kept as a clearly-labelled alternative.
    """
    pool = generate_aptamer_pool(is_rna=False, n_per_arch=n_per_arch,
                                 len_min=len_min, len_max=len_max, seed=seed)
    extra = {c["seq"]: {"core": c["seq"], "n_mut": None, "parent": None,
                        "sublib": "S_short_unframed", "arch": c["arch"]}
             for c in pool}
    return score_pool(list(extra.keys()), extra=extra)


# ------------------------------------------------------- G-quadruplex -------
def library_G_g4(n=700, seed=5, len_min=15, len_max=32):
    """Compact G-quadruplex switches, ranked by `g4_quality`, NOT by eab_score.

    Built to the canonical intramolecular grammar: exactly four G-tracts of 2-4
    G separated by G-free loops of 1-5 nt, with short A/C/T flanks.  Keeping
    the loops G-free guarantees four resolvable tracts instead of one long G
    run - the failure mode that raw G4Hunter maximisation selects for.

    Kept as its own track because ViennaRNA cannot see a G4 at all, so these
    sequences can never compete on the secondary-structure score no matter how
    good they are.  Worth the trouble here: the thrombin E-AB sensor - the
    archetype of the whole method - is a 15-nt G4; G4 DNA is compact (fast
    electron transfer, large relative signal change) and intrinsically
    nuclease-resistant, which matters for serum or in vivo operation.
    """
    import pandas as pd
    from .eab_filters import g4_metrics, TBA
    rng = random.Random(seed)
    seen, rows = set(), []
    for _ in range(n * 25):
        if len(rows) >= n:
            break
        tracts = ["G" * rng.choice([2, 2, 3, 3, 3, 4]) for _ in range(4)]
        loops = ["".join(rng.choice("ACT") for _ in range(rng.randint(1, 5)))
                 for _ in range(3)]
        core = "".join(t + l for t, l in zip(tracts, loops + [""]))
        s = ("".join(rng.choice("ACT") for _ in range(rng.randint(0, 3)))
             + core
             + "".join(rng.choice("ACT") for _ in range(rng.randint(0, 3))))
        if not (len_min <= len(s) <= len_max) or s in seen:
            continue
        seen.add(s)
        m = g4_metrics(s)
        if m["g4_quality"] < 50:
            continue
        m.update({"core": s, "n_mut": None, "parent": None,
                  "sublib": "G_quadruplex", "arch": "gquad"})
        rows.append(m)

    # TBA as the track's reference point, flagged so it is never mistaken for
    # a designed candidate.
    ref = g4_metrics(TBA)
    ref.update({"core": TBA, "n_mut": None, "parent": None,
                "sublib": "G_quadruplex_REFERENCE", "arch": "gquad_TBA"})
    df = pd.DataFrame(rows + [ref])
    # Tie-break on length: g4_quality plateaus at 100 for anything sitting in
    # the TBA-like window with four clean tracts, and among equals the shorter
    # construct wins - shorter means faster electron transfer to the electrode
    # and a larger relative signal change.
    return (df.sort_values(["g4_quality", "len"], ascending=[False, True])
              .reset_index(drop=True))


# ------------------------------------------------------------- controls -----
def controls():
    scr_core = scramble(known.APT431_CORE)
    seqs = {
        "CTRL_scramble_4.31_core": frame(scr_core),
        "CTRL_4.31_core_only_published_inactive": known.APT431_CORE,
        "REF_4.31_full80_published_binder": known.APT431,
    }
    extra = {s: {"core": None, "n_mut": None, "parent": None,
                 "sublib": "CONTROL", "arch": k}
             for k, s in seqs.items()}
    df = score_pool(list(seqs.values()), extra=extra)
    inv = {v: k for k, v in seqs.items()}
    df["id"] = df["seq"].map(inv)
    return df


# ------------------------------------------------------------ selection -----
def pick_panel(df, k=16, min_edit=6, exclude=None):
    """Greedy diversity-aware top-k: take the best-scoring candidates but never
    two whose CORES are within `min_edit` edits of each other, so the ordered
    panel samples the space instead of ordering one sequence 16 times."""
    exclude = exclude or []
    picked = []
    for _, r in df.iterrows():
        core = r.get("core") or r["seq"]
        if any(levenshtein(core, p) < min_edit for p in exclude):
            continue
        if any(levenshtein(core, (p.get("core") or p["seq"])) < min_edit
               for p in picked):
            continue
        picked.append(r.to_dict())
        if len(picked) >= k:
            break
    import pandas as pd
    return pd.DataFrame(picked)


def order_specs():
    """[DESIGN] Sub-library E - what to send to the oligo synthesiser."""
    return [
        {"name": "E1_naive_N40_SELEX_library",
         "spec": known.selex_library_template(40),
         "use": "naive CE-SELEX / bead-SELEX round 1 for either target; "
                "same arms as the published NPY selection so the workflow "
                "and the resulting sensor stay directly comparable",
         "scale": "1 umol, desalted; ~10^14-10^15 unique molecules in a "
                  "typical 1 nmol input aliquot"},
        {"name": "E2_doped_4.31_NPY_maturation",
         "spec": known.doped_oligo_spec(doping=0.15)["construct"],
         "use": "affinity/switch maturation of the existing NPY binder; "
                "run with increasing stringency + PP and PYY counter-selection",
         "scale": "0.2 umol, desalted"},
        {"name": "E3_doped_4.31_PP_retargeting",
         "spec": known.doped_oligo_spec(doping=0.30)["construct"],
         "use": "re-target the 4.31 backbone to PP; counter-select against NPY "
                "and PYY in EVERY round from round 1 - this is the step that "
                "inverts the 42-fold NPY preference",
         "scale": "0.2 umol, desalted"},
        {"name": "E4_primers",
         "spec": f"FWD 5'-{known.LIB_FWD}-3'  /  REV 5'-{known.LIB_REV}-3' "
                 f"(order the reverse primer 5'-phosphorylated or "
                 f"poly-A-tailed for ssDNA regeneration)",
         "use": "amplification of every round",
         "scale": "100 nmol each"},
    ]
