"""
pathways.py — pathway component definitions and candidate model-organism panels.

Design principle (matches the repo's ethos): every *number* downstream is
computed, and every *assignment* here that is curated rather than computed is
labelled as such.

  - COMPONENT FAMILIES: which gene families make up the pathway core, and how
    many human paralogues each has. Paralogue counts are used for the
    "low-redundancy / deconvolution" axis (a single-copy pathway means one
    knockout gives a clean phenotype).
  - GATE FAMILIES: the components whose *absence* means "pathway absent" — the
    hard gate (e.g. no Notch receptor or no CSL transcription factor ⇒ the
    organism cannot run canonical Notch signalling ⇒ disqualified as a model).
  - ORTHOLOG SEED TABLE: literature-curated per-species gene symbols for the
    gate families, so we can pull the *real* ortholog sequence from UniProt and
    hand it to Evo2. The ortholog *assignment* is curated; the sequence it
    resolves to and the Evo2 constraint score computed from it are not.

Candidate species mirror the Notch dashboard panel (10 with the pathway + 2
negative controls where canonical Notch is genuinely absent).
"""

from dataclasses import dataclass, field


# --------------------------------------------------------------------------- #
# Candidate species panel: display name -> NCBI taxon id + tractability priors
# --------------------------------------------------------------------------- #
# `tract` (experimental tractability) and `thru` (throughput) are curated
# priors from the model-organism literature (CRISPR/transgenesis maturity,
# generation time, imaging access, cost). They feed the X axis. Labelled
# CURATED — they are not sequence-derived and Evo2 does not touch them.
@dataclass(frozen=True)
class Species:
    name: str
    taxon: int
    tract: float   # 0..1 experimental tractability (CURATED prior)
    thru: float    # 0..1 throughput / cost-speed (CURATED prior)
    clade: str


SPECIES: dict[str, Species] = {
    "human":       Species("Homo sapiens (human, ref)",            9606,  0.55, 0.10, "mammal"),
    "mouse":       Species("Mus musculus (mouse)",                10090,  0.80, 0.30, "mammal"),
    "rat":         Species("Rattus norvegicus (rat)",             10116,  0.60, 0.25, "mammal"),
    "chicken":     Species("Gallus gallus (chicken)",              9031,  0.50, 0.35, "bird"),
    "frog":        Species("Xenopus tropicalis (frog)",            8364,  0.65, 0.55, "amphibian"),
    "zebrafish":   Species("Danio rerio (zebrafish)",              7955,  0.90, 0.85, "fish"),
    "tunicate":    Species("Ciona intestinalis (tunicate)",        7719,  0.50, 0.55, "tunicate"),
    "fly":         Species("Drosophila melanogaster (fly)",        7227,  0.95, 0.90, "insect"),
    "worm":        Species("Caenorhabditis elegans (worm)",        6239,  0.95, 0.95, "nematode"),
    "sea_anemone": Species("Nematostella vectensis (sea anemone)", 45351, 0.40, 0.45, "cnidarian"),
    "yeast":       Species("Saccharomyces cerevisiae (yeast)",   559292,  0.95, 0.98, "fungus"),
    "plant":       Species("Arabidopsis thaliana (plant)",         3702,  0.85, 0.80, "plant"),
    # Regeneration-focused species (Hippo–YAP panel). CURATED tractability priors.
    "axolotl":     Species("Ambystoma mexicanum (axolotl)",        8296,  0.50, 0.30, "amphibian"),
    "naked_mole_rat": Species("Heterocephalus glaber (naked mole-rat)", 10181, 0.31, 0.20, "mammal"),
    "planaria":    Species("Schmidtea mediterranea (planaria)",   79327,  0.76, 0.80, "flatworm"),
}

REFERENCE = "human"


# --------------------------------------------------------------------------- #
# Pathway component families
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Family:
    key: str
    label: str
    human_genes: list[str]     # human paralogues (drives redundancy axis)
    is_gate: bool = False      # absence ⇒ pathway absent ⇒ disqualified


@dataclass(frozen=True)
class Pathway:
    key: str
    label: str
    reference_species: str
    families: list[Family]
    # per-species gene symbol for each GATE family (curated ortholog assignment)
    # {family_key: {species_key: [symbols]}}
    ortholog_seed: dict[str, dict[str, list[str]]] = field(default_factory=dict)

    @property
    def gate_families(self) -> list[Family]:
        return [f for f in self.families if f.is_gate]

    @property
    def human_mean_paralogs(self) -> float:
        return sum(len(f.human_genes) for f in self.families) / len(self.families)


# ---- Notch --------------------------------------------------------------- #
# 7 component families. Gate = receptor + CSL transcription factor: with no
# receptor and no CSL there is no canonical Notch signalling at all.
NOTCH = Pathway(
    key="Notch",
    label="Notch signalling",
    reference_species="human",
    families=[
        Family("receptor", "Notch receptor",
               ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"], is_gate=True),
        Family("dsl_ligand", "DSL ligand",
               ["DLL1", "DLL3", "DLL4", "JAG1", "JAG2"]),
        Family("csl", "CSL / RBPJ transcription factor",
               ["RBPJ", "RBPJL"], is_gate=True),
        Family("maml", "Mastermind co-activator",
               ["MAML1", "MAML2", "MAML3"]),
        Family("fringe", "Fringe glycosyltransferase",
               ["LFNG", "MFNG", "RFNG"]),
        Family("gamma_secretase", "gamma-secretase / S3 processing",
               ["PSEN1", "PSEN2", "NCSTN", "APH1A", "PSENEN"]),
        Family("hes_hey", "HES/HEY bHLH effector",
               ["HES1", "HES5", "HEY1", "HEY2"]),
    ],
    # Curated per-species gene symbols for the two GATE families.
    # Distant species use their own historical gene names (fly N/Su(H),
    # worm lin-12+glp-1 / lag-1). Yeast & plant intentionally have NO entry:
    # canonical Notch is genuinely absent — the pipeline must find nothing and
    # the hard gate must disqualify them.
    ortholog_seed={
        "receptor": {
            "human":       ["NOTCH1"],
            "mouse":       ["Notch1"],
            "rat":         ["Notch1"],
            "chicken":     ["NOTCH1"],
            "frog":        ["notch1"],
            "zebrafish":   ["notch1a"],
            "tunicate":    ["notch"],
            "fly":         ["N"],          # Drosophila Notch
            "worm":        ["lin-12"],     # C. elegans Notch receptor
            "sea_anemone": ["notch"],
            # Negative controls: query the human symbol so the search really
            # runs and really finds nothing (empirical gate, not asserted).
            "yeast":       ["NOTCH1"],
            "plant":       ["NOTCH1"],
        },
        "csl": {
            "human":       ["RBPJ"],
            "mouse":       ["Rbpj"],
            "rat":         ["Rbpj"],
            "chicken":     ["RBPJ"],
            "frog":        ["rbpj"],
            "zebrafish":   ["rbpja", "rbpjb"],   # fallback order
            "tunicate":    ["su(h)"],
            "fly":         ["Su(H)"],      # Suppressor of Hairless
            "worm":        ["lag-1"],      # C. elegans CSL
            "sea_anemone": ["su(h)"],
            # Negative controls (see receptor family).
            "yeast":       ["RBPJ"],
            "plant":       ["RBPJ"],
        },
    },
)


# ---- Hippo–YAP/TAZ–TEAD (regeneration) ----------------------------------- #
# 7 component families. Gate = the transcriptional output module: YAP/TAZ
# effector + TEAD DNA-binding partner. With no YAP/TAZ and no TEAD there is no
# canonical Hippo transcriptional output (analogous to Notch receptor + CSL).
HIPPO = Pathway(
    key="Hippo",
    label="Hippo–YAP/TAZ–TEAD signalling",
    reference_species="human",
    families=[
        Family("mst_kinase", "MST1/2 (Hpo)", ["STK3", "STK4"]),
        Family("sav", "SAV1 (Sav)", ["SAV1"]),
        Family("lats_kinase", "LATS1/2 (Wts)", ["LATS1", "LATS2"]),
        Family("mob", "MOB1 (Mats)", ["MOB1A", "MOB1B"]),
        Family("yap_taz", "YAP/TAZ (Yki) effector", ["YAP1", "WWTR1"], is_gate=True),
        Family("tead", "TEAD1-4 (Sd) transcription factor",
               ["TEAD1", "TEAD2", "TEAD3", "TEAD4"], is_gate=True),
        Family("upstream", "NF2 / KIBRA / FAT4 mechanical input",
               ["NF2", "WWC1", "FAT4"]),
    ],
    ortholog_seed={
        "yap_taz": {
            "human":          ["YAP1", "WWTR1"],
            "mouse":          ["Yap1", "Wwtr1"],
            "zebrafish":      ["yap1", "wwtr1"],
            "fly":            ["yki"],                 # yorkie
            "naked_mole_rat": ["YAP1", "WWTR1"],
            "axolotl":        ["yap1", "wwtr1"],
            "planaria":       ["yki", "yap1"],
        },
        "tead": {
            "human":          ["TEAD1"],
            "mouse":          ["Tead1"],
            "zebrafish":      ["tead1a", "tead1b"],
            "fly":            ["sd"],                  # scalloped
            "naked_mole_rat": ["TEAD1"],
            "axolotl":        ["tead1"],
            "planaria":       ["tead", "tead1"],
        },
    },
)


PATHWAYS = {p.key: p for p in [NOTCH, HIPPO]}


def get_pathway(key: str) -> Pathway:
    return PATHWAYS[key]
