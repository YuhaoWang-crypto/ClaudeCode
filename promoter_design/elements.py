"""
Stimulus-responsive DNA modules (response elements) and minimal promoters.

Each entry pairs a stimulus with the transcription factor that terminates its
signalling pathway and the cis-regulatory element that TF binds. The `core`
field is the well-established consensus motif; `seed` is a representative
functional monomer (canonical binding site + minimal flanks) used as the
STARTING point for in-silico optimisation. These are seeds, not final answers:
the design pipeline mutates/multimerises them and ranks variants with Evo2 and
AlphaGenome, so absolute optimality of the seed is not assumed.

Sequences are forward strand, 5'->3'. Uppercase = TF binding site,
lowercase = neutral spacer. IUPAC codes: R=A/G, Y=C/T, W=A/T, N=any.

References (representative): HRE/EPO enhancer (Semenza & Wang 1992);
kB/Ig-HIV LTR (Sen & Baltimore 1986); CRE/somatostatin (Montminy 1986);
ISRE & GAS (Darnell, Kerr & Stark 1994); ERE/vitellogenin A2 (Klein-Hitpass
1986); GRE consensus (Beato 1989); ARE/NQO1 (Rushmore & Pickett 1990);
UPRE/ERSE/AARE (Yoshida, Kaufman, Ron & Walter, 2000s).
"""

# --- Stimulus -> pathway TF -> DNA module -------------------------------------
# Each module: consensus core, a functional monomer seed, and the recommended
# spacer used when tandemly multimerising the monomer.

ELEMENTS = {
    "hypoxia": {
        "stimulus": "Hypoxia / low O2 (also mimicked by CoCl2, DMOG)",
        "tf": "HIF1A-ARNT (HIF-1)",
        "element": "HRE (Hypoxia Response Element)",
        "core": "RCGTG",                       # HIF binding site (HBS)
        "seed": "gtcTACGTGCTgtc",              # HBS from EPO 3' enhancer
        "spacer": "gcgtc",
        "notes": "Multimerise 4-6x. Pair HBS with a downstream HIF ancillary "
                 "sequence (HAS, ~CACAG) for full induction.",
    },
    "inflammation": {
        "stimulus": "TNF-alpha, IL-1, LPS (NF-kB inflammatory input)",
        "tf": "RELA / NF-kB (p65:p50)",
        "element": "kB response element",
        "core": "GGGRNWYYCC",
        "seed": "gcGGGACTTTCCgc",              # Ig-kappa / HIV-1 LTR canonical kB
        "spacer": "gggac",
        "notes": "Canonical strong kB site. 3-5x tandem gives clean TNF/LPS "
                 "induction; keep spacers short to preserve dimer spacing.",
    },
    "camp": {
        "stimulus": "cAMP / forskolin / beta-adrenergic (PKA -> CREB)",
        "tf": "CREB1 (and CREM/ATF1)",
        "element": "CRE (cAMP Response Element)",
        "core": "TGACGTCA",                    # perfect palindrome
        "seed": "ctTGACGTCAct",
        "spacer": "tcacc",
        "notes": "Somatostatin-type CRE. 4-6x tandem; half-site CGTCA also "
                 "responsive but weaker.",
    },
    "interferon_typeI": {
        "stimulus": "Type I interferon (IFN-alpha/beta) -> ISGF3",
        "tf": "STAT1/STAT2/IRF9 (ISGF3)",
        "element": "ISRE (Interferon-Stimulated Response Element)",
        "core": "RNGAAANNGAAACT",
        "seed": "gGGGAAAGTGAAACTg",            # ISG15-type ISRE
        "spacer": "gaacc",
        "notes": "Type-I-IFN selective. Distinct from GAS (below); use ISRE "
                 "for IFN-alpha/beta, GAS for IFN-gamma.",
    },
    "interferon_typeII": {
        "stimulus": "Type II interferon (IFN-gamma) -> STAT1 homodimer",
        "tf": "STAT1 (GAF)",
        "element": "GAS (Gamma-Activated Sequence)",
        "core": "TTCNNNGAA",                   # palindromic
        "seed": "ctTTCCCGGAAct",               # GRR / FcgRI-type GAS
        "spacer": "tgcaa",
        "notes": "IFN-gamma selective. Combine an ISRE + GAS block for a "
                 "pan-interferon sensor.",
    },
    "estrogen": {
        "stimulus": "Estrogen / 17beta-estradiol",
        "tf": "ESR1 (ER-alpha)",
        "element": "ERE (Estrogen Response Element)",
        "core": "GGTCANNNTGACC",               # 13 bp inverted palindrome
        "seed": "caAGGTCACAGTGACCTca",         # vitellogenin A2 ERE
        "spacer": "gatcg",
        "notes": "Perfect palindrome = strongest. Only active in ER+ cells "
                 "(e.g. MCF-7); silent in ER- lines -> intrinsic cell selectivity.",
    },
    "glucocorticoid": {
        "stimulus": "Glucocorticoid (dexamethasone, cortisol)",
        "tf": "NR3C1 / GR",
        "element": "GRE (Glucocorticoid Response Element)",
        "core": "GGTACANNNTGTTCT",             # 15 bp imperfect palindrome
        "seed": "gaGGTACAGGATGTTCTga",
        "spacer": "gatct",
        "notes": "GR and MR/PR/AR share this; combine with cell-restricted GR "
                 "expression for selectivity. 2-4x tandem.",
    },
    "oxidative_stress": {
        "stimulus": "Oxidative / electrophilic stress (sulforaphane, tBHQ, H2O2)",
        "tf": "NFE2L2 / NRF2 (with sMAF)",
        "element": "ARE / EpRE (Antioxidant Response Element)",
        "core": "RTGACNNNGCR",
        "seed": "gtGTGACTCAGCAgaatctg",        # human NQO1 ARE
        "spacer": "gcaat",
        "notes": "The 3' GC and flanks give NRF2 selectivity vs AP-1/TRE. "
                 "Keep the GAATC-type flank of the NQO1 ARE.",
    },
    # ER stress groups three effector arms; provide all three sub-elements.
    "er_stress_xbp1": {
        "stimulus": "ER stress (UPR, IRE1 arm) -> spliced XBP1",
        "tf": "XBP1s",
        "element": "UPRE (Unfolded Protein Response Element)",
        "core": "TGACGTGG",
        "seed": "ctTGACGTGGCAGAct",            # UPRE-1
        "spacer": "gatca",
        "notes": "IRE1/XBP1s-selective arm of the UPR.",
    },
    "er_stress_atf6": {
        "stimulus": "ER stress (UPR, ATF6 arm)",
        "tf": "ATF6 (with NF-Y)",
        "element": "ERSE-I (ER Stress Response Element)",
        "core": "CCAATN9CCACG",
        "seed": "CCAATgacgctgacCCACG",         # CCAAT-N9-CCACG
        "spacer": "gatcg",
        "notes": "Requires the CCAAT (NF-Y) + N9 + CCACG geometry; do NOT "
                 "collapse the 9-nt spacer.",
    },
    "er_stress_atf4": {
        "stimulus": "ER stress / amino-acid starvation (PERK -> ATF4)",
        "tf": "ATF4 (with CHOP/DDIT3)",
        "element": "AARE / CARE (Amino Acid Response Element)",
        "core": "ATTGCATCA",                   # CHOP AARE
        "seed": "gaATTGCATCATga",
        "spacer": "tgatg",
        "notes": "PERK/ATF4 arm; also induced by amino-acid limitation.",
    },
}

# --- Cell-lineage-restricted modules (for stimulus-AND-cell-type gating) ------
# These are bound by lineage-defining TFs that are expressed only in certain
# cell types. Combined with a stimulus module upstream of a minimal promoter,
# they give analog AND behaviour: strong output needs BOTH the signal (stimulus
# RE occupied) AND the right cell (lineage RE occupied). Cores are well
# established; validate lineage selectivity with AlphaGenome contrastive scoring.

LINEAGE_ELEMENTS = {
    "myeloid": {          # monocyte / macrophage (THP-1)
        "tf": "SPI1 / PU.1 (ETS)",
        "core": "GAGGAA",
        "seed": "aaAGAGGAAGTGaa",
        "spacer": "gccta",
        "ontology": "CL:0000576",   # monocyte
    },
    "hepatocyte": {       # liver (HepG2)
        "tf": "HNF4A (DR1)",
        "core": "AGGTCAAAGGTCA",
        "seed": "cAGGTCAAAGGTCAg",
        "spacer": "gttaa",
        "ontology": "CL:0000182",   # hepatocyte
    },
    "hepatocyte_hnf1": {  # liver, orthogonal HNF1 site
        "tf": "HNF1A",
        "core": "GTTAATNATTAAC",
        "seed": "gGTTAATGATTAACg",
        "spacer": "gatca",
        "ontology": "CL:0000182",
    },
    "hematopoietic": {    # erythroid / blood (K562, Jurkat)
        "tf": "GATA1/2 (WGATAR)",
        "core": "AGATAA",
        "seed": "agAGATAAGCag",
        "spacer": "gcaat",
        "ontology": "CL:0000094",   # granulocyte-ish blood lineage placeholder
    },
    "neuronal": {         # neuron (SH-SY5Y)
        "tf": "NeuroD / E-box",
        "core": "CAGCTG",
        "seed": "gcCAGCTGgc",
        "spacer": "gttca",
        "ontology": "CL:0000540",   # neuron
    },
}

# --- Composite sensors: multi-element single-cassette designs -----------------
# pan_interferon: ISRE + GAS in one enhancer -> fires on type I OR type II IFN.
COMPOSITE = {
    "pan_interferon": {
        "description": "Pan-interferon sensor: type I (ISRE) + type II (GAS)",
        "modules": ["interferon_typeI", "interferon_typeII"],  # keys in ELEMENTS
        "layout": "interleave",   # ISRE-GAS-ISRE-GAS...
    },
}

# --- Minimal (core) promoters -------------------------------------------------
# The enhancer module sets WHAT signal turns the promoter on; the minimal
# promoter sets baseline strength and contributes cell-type tuning. Swap these
# and re-score in AlphaGenome per target cell type.

MINIMAL_PROMOTERS = {
    "minCMV": (
        "taggcgtgtacggtgggaggcctatataagcagagctcgtttagtgaaccgtcagatc"
    ),  # widely used minimal CMV (minP); strong, broadly active
    "miniTK": (
        "ccccgcccagcgtcttgtcattggcgaattcgaacacgcagatgcagtcggggcggcg"
        "cggtccc"
    ),  # HSV-TK minimal promoter; moderate, broad
    "E1b_TATA": (
        "gggtatataatggggggcc"
    ),  # adenovirus E1b minimal TATA+Inr; very low baseline, clean induction
    "YB_TATA": (
        "ggggactttccacggggggctataaaagggggtgggg"
    ),  # synthetic super-core (TATA + Inr); tunable
}

# --- Cell-context hints -------------------------------------------------------
# Which pathways are wired/active in common lines, and which AlphaGenome/Enformer
# track families to score against. Use to (a) pick a stimulus that CAN fire in
# the cell and (b) filter for ON-in-target / OFF-elsewhere designs.

# `ontology` = AlphaGenome ontology term for that cell/tissue, used as the
# target (ontology_terms) or off-target (contrastive_ontology_terms) in scoring.
# Tissue-level UBERON/CL terms are given (more robustly supported); confirm the
# exact term against your AlphaGenome build's output_metadata before running
# (cell-line-specific EFO terms like EFO:xxxx for HepG2/THP-1 are often available
# and give sharper contrast).
CELL_CONTEXTS = {
    "HEK293":  {"tissue": "embryonic kidney", "ontology": "UBERON:0002113",  # kidney
                "good_for": ["camp", "er_stress_xbp1", "er_stress_atf6",
                "er_stress_atf4", "hypoxia"],
                "weak_for": ["estrogen", "interferon_typeII"]},
    "HepG2":   {"tissue": "hepatocyte", "ontology": "UBERON:0002107",  # liver
                "good_for": ["glucocorticoid", "oxidative_stress", "hypoxia",
                "er_stress_atf6"], "weak_for": ["estrogen"]},
    "THP1":    {"tissue": "monocyte/macrophage", "ontology": "CL:0000576",  # monocyte
                "good_for": ["inflammation", "interferon_typeI",
                "interferon_typeII"], "weak_for": ["estrogen", "glucocorticoid"]},
    "MCF7":    {"tissue": "ER+ breast", "ontology": "UBERON:0000310",  # breast
                "good_for": ["estrogen", "glucocorticoid", "hypoxia"],
                "weak_for": ["interferon_typeII"]},
    "Jurkat":  {"tissue": "T lymphocyte", "ontology": "CL:0000084",  # T cell
                "good_for": ["inflammation", "interferon_typeI"],
                "weak_for": ["estrogen"]},
    "SHSY5Y":  {"tissue": "neuronal", "ontology": "CL:0000540",  # neuron
                "good_for": ["camp", "oxidative_stress"],
                "weak_for": ["estrogen", "inflammation"]},
}
