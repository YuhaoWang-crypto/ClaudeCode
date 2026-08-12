"""Literature-anchored aptamer sequences for the NPY / PP family.

Everything in this file is [MEASURED] - transcribed from a publication or a
curated aptamer database - EXCEPT where a comment says otherwise.  These are
the binding priors the designed library is seeded from; do not edit a sequence
without updating its citation.
"""

# ---------------------------------------------------------------------------
# 1. NPY ssDNA aptamer "4.31"  (the only ssDNA NPY aptamer with a published
#    EAB sensor built on it - this is the single most valuable prior we have)
# ---------------------------------------------------------------------------
# [MEASURED] Selected by CE-SELEX in 4 rounds from an 80-nt ssDNA library.
#   Mendonsa SD, Bowser MT. J Am Chem Soc 2005;127(26):9382-3.
#   doi:10.1021/ja052406n  (PMID 15984861)
#   Kd(NPY) = 0.3 +/- 0.2 uM; pool Kd range 300-1000 nM.
#   Specificity vs human pancreatic polypeptide (hPP, ~50% identity): 42-fold.
#     -> NOTE: 42-fold means 4.31 DOES bind hPP, roughly 12 uM.  That residual
#        affinity is exactly why 4.31 is also a legitimate SEED for a PP
#        binder (see build_library.py, sub-library C).
#
# [MEASURED] Re-used as the recognition element of a working E-AB sensor in
#   Seibold JM, Abeykoon SW, Ross AE, White RJ. ACS Sens 2023;8(12).
#   doi:10.1021/acssensors.3c00855 (PMID 38033269, PMC11214579):
#     architecture   5'-HS-(CH2)6 - aptamer - MB-3'   ("terminally labelled", TL)
#     signal         signal-ON, square-wave voltammetry at 100 Hz
#     Kd             385 nM in buffer / 56 nM in serum
#     dynamic range  20-600 nM;  LOD 162 pM (buffer), 390 pM (serum)
#     max signal     230 +/- 15 % signal change
#     selectivity    PYY 4 +/- 1 % SC; IFN-gamma 3 +/- 2 % SC; CRP negligible
#     CRITICAL       a 40-nt truncation stripped of the primer regions ("STLA")
#                    gave NO appreciable signal.  The fixed primer arms are
#                    part of the binding-competent / switching structure.
#                    => do NOT hand the synthesiser a bare 40-mer core.
APT431_FWD_PRIMER = "AGCAGCACAGAGGTCAGATG"          # 20 nt, 5' fixed arm
APT431_CORE = "CAAACCACAGCCTGAGTGGTTAGCGTATGTCATTTACGGA"   # 40 nt, selected core
APT431_REV_PRIMER = "CCTATGCGTGCTACCGTGAA"          # 20 nt, 3' fixed arm
APT431 = APT431_FWD_PRIMER + APT431_CORE + APT431_REV_PRIMER   # 80 nt

# [MEASURED] The internally-labelled variant from the same paper: MB moved to
# the 33rd position (after ...AGCCT).  Reported as a design variable, not a
# different sequence.
APT431_IL_MB_AFTER = 33   # 1-based nucleotide after which MB was placed

# [MEASURED] Reported by the same authors as inactive in the EAB format -
# retained here as an explicit NEGATIVE control, not as a candidate.
APT431_TRUNC_INACTIVE = APT431_CORE

# ---------------------------------------------------------------------------
# 2. NPY RNA aptamer "DP3"  (different chemistry, still a useful prior)
# ---------------------------------------------------------------------------
# [MEASURED] Proske D, Hofliger M, Soll RM, Beck-Sickinger AG, Famulok M.
#   "A Y2 receptor mimetic aptamer directed against neuropeptide Y."
#   J Biol Chem 2002;277(13):11416-22.
#   Chemistry: 2'-amino-2'-deoxy pyrimidine RNA (nuclease-stabilised), 95 nt,
#   12 SELEX rounds against N-terminally biotinylated NPY.  Kd = 370 nM.
#   Functionally a Y2-receptor mimic -> it reads the C-TERMINAL face of NPY,
#   which is the face conserved across NPY/PYY/PP.  Expect poor family
#   discrimination from anything derived from this epitope.
DP3_RNA = ("UCGGAGAAAGGGAAGCUUGAGCAGCAGGAGGGCCGGCGUUAGGGUUAGCGAGCCGAUUGA"
           "AAGAAGAAGGAACGAGCGUACGGAUCCGAUC")

# [DESIGN, low prior] Plain-DNA transliteration of the DP3 core.  A U->T swap
# does NOT preserve an RNA fold (no 2'-OH, different helical geometry), so this
# is included only as a diversity seed, clearly flagged.  Do not present it as
# an NPY binder.
DP3_AS_DNA = DP3_RNA.replace("U", "T")

# ---------------------------------------------------------------------------
# 3. Pancreatic polypeptide (PP / PPY)
# ---------------------------------------------------------------------------
# [MEASURED, sequences NOT retrievable here] PP-specific DNA aptamers were
#   selected by SELEX and used in a colorimetric AuNP aptasensor:
#     Ali ASM, El-Halawany MS, Ibrahim SA, Pluckthun O, Khalil ASG, Mayer G.
#     "Aptasensor for Quantifying Pancreatic Polypeptide."
#     ACS Omega 2019;4(2):2948-56.  doi:10.1021/acsomega.8b03131
#     Reported performance: linear range 1-60 nM, LOD 521 pM, serum recovery
#     92.4-93.9 %.
#   The nucleotide sequences live in that paper's text/SI.  pubs.acs.org
#   returns HTTP 403 from this sandbox and the article is not mirrored in
#   Europe PMC, so the sequences are NOT included here.  ACTION FOR THE USER:
#   pull the sequences from the PDF/SI and add them to PP_KNOWN below - they
#   are the highest-value seeds available for the PP sub-library.
PP_KNOWN = []          # <- populate from ACS Omega 2019 SI

# ---------------------------------------------------------------------------
# 4. Fixed arms for the new library
# ---------------------------------------------------------------------------
# [DESIGN] Reusing the Bowser/White primer arms verbatim is the lowest-risk
# choice: they are already proven to (a) amplify in the CE-SELEX workflow and
# (b) support a 230 % signal-ON EAB response with 4.31.  Keeping them means the
# new library is directly comparable to the published sensor.
LIB_FWD = APT431_FWD_PRIMER
LIB_REV = APT431_REV_PRIMER


def selex_library_template(n_random=40):
    """[DESIGN] The degenerate oligo to order for a naive SELEX round."""
    return f"5'-{LIB_FWD}-(N){n_random}-{LIB_REV}-3'"


def doped_oligo_spec(core=APT431_CORE, doping=0.22):
    """[DESIGN] Synthesiser spec for a doped affinity-maturation library.

    A doped library is ordered as ONE degenerate oligo, not as a list of
    sequences: at each core position the synthesiser mixes `1-doping` of the
    wild-type phosphoramidite with `doping/3` of each of the other three.
    """
    wt = round(100 * (1 - doping))
    other = round(100 * doping / 3, 1)
    return {
        "construct": f"5'-{LIB_FWD}-[{core} @ {doping:.0%} doping]-{LIB_REV}-3'",
        "length_nt": len(LIB_FWD) + len(core) + len(LIB_REV),
        "per_position_mix": f"{wt}% wild-type / {other}% each of the other 3 bases",
        "expected_mutations_per_molecule": round(doping * len(core), 1),
        "note": "order as a single doped oligo; the sequence lists in "
                "output/*.csv are the in-silico sample of this same space, "
                "for anyone who prefers a defined array instead.",
    }


KNOWN_NPY = {
    "APT431_full80_ssDNA": APT431,
    "APT431_core40_ssDNA_INACTIVE_in_EAB": APT431_CORE,
    "DP3_2NH2_RNA": DP3_RNA,
}
