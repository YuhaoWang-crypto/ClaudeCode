"""Shared configuration for the cell-surface antigen discovery pipeline.

Every constant that changes a number in the final report lives here, so the
report generator can print the exact settings a run used.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures" / "antigen"
REPORTS = ROOT / "reports"

for _p in (DATA_RAW, RESULTS, FIGURES, REPORTS):
    _p.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Data source pins (all public, all fetched at run time)
# --------------------------------------------------------------------------
SURFY_URL = "https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx"
SURFY_FILE = DATA_RAW / "table_S3_surfaceome.xlsx"

HPA_RNA_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
HPA_RNA_FILE = DATA_RAW / "rna_tissue_consensus.tsv.zip"

HPA_IHC_URL = "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
HPA_IHC_FILE = DATA_RAW / "normal_ihc_data.tsv.zip"

# DepMap 23Q4 CRISPRGeneEffect (Chronos). Annotation only -- never a filter.
DEPMAP_URL = "https://ndownloader.figshare.com/files/43346616"
DEPMAP_FILE = DATA_RAW / "CRISPRGeneEffect.csv"

OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
EUROPEPMC_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# --------------------------------------------------------------------------
# Step 2 -- tumour single-cell cohort
# --------------------------------------------------------------------------
CENSUS_VERSION = "2025-11-08"          # LTS release
DISEASE_LABEL = "lung adenocarcinoma"
# Single-nucleus data systematically under-detects membrane-protein transcripts,
# so nucleus suspensions are excluded from the discovery cohort.
SUSPENSION_TYPE = "cell"
MAX_CELLS_PER_DATASET = 20_000
SAMPLING_SEED = 0
CENSUS_LAYER = "normalized"

COMPARTMENTS = ["epithelial", "caf", "immune", "endothelial"]

# Census stores library-size-normalised values (per-cell fractions, ~1e-5 for a
# typical gene). Everything downstream works in counts per 10,000 so that the
# pseudocount in the specificity ratio stays small relative to real expression.
EXPRESSION_SCALE = 1e4
EXPRESSION_PSEUDOCOUNT = 0.05        # CP10K
# An atlas counts towards consensus when the target is both enriched in the
# epithelial compartment and actually detected there.
CONSENSUS_MIN_LOG2FC = 0.5
CONSENSUS_MIN_DETECTION = 0.10

# Substring rules applied to the Cell Ontology label, first match wins.
# Ordering matters: "malignant" and epithelial rules are tested before stroma.
COMPARTMENT_RULES: list[tuple[str, str]] = [
    ("malignant", "epithelial"),
    ("epithelial", "epithelial"),
    ("alveolar type", "epithelial"),
    ("pulmonary alveolar type", "epithelial"),
    ("club cell", "epithelial"),
    ("goblet", "epithelial"),
    ("basal cell", "epithelial"),
    ("brush cell", "epithelial"),
    ("ionocyte", "epithelial"),
    ("hillock", "epithelial"),
    ("multiciliated", "epithelial"),
    ("endothelial", "endothelial"),
    ("fibroblast", "caf"),
    ("myofibroblast", "caf"),
    ("stromal cell", "caf"),
    ("t cell", "immune"),
    ("t-regulatory", "immune"),
    ("b cell", "immune"),
    ("plasma cell", "immune"),
    ("natural killer", "immune"),
    ("macrophage", "immune"),
    ("monocyte", "immune"),
    ("dendritic cell", "immune"),
    ("mast cell", "immune"),
    ("neutrophil", "immune"),
    ("granulocyte", "immune"),
    ("myeloid cell", "immune"),
    ("erythrocyte", "immune"),
]
# Cell types that match none of the rules (unknown, smooth muscle, pericyte,
# mesothelial) are carried in the cohort table but excluded from compartment
# statistics -- they belong to no compartment of the four-way model.

# --------------------------------------------------------------------------
# Step 3 -- extracellular topology gate
# --------------------------------------------------------------------------
# UniProt writes subcellular location as a semicolon-separated list in which the
# organelle comes first, so each entry is judged by its prefix. A bare
# "Membrane (Single-pass type I membrane protein)" is UniProt's way of saying
# plasma membrane; "Endoplasmic reticulum membrane (Single-pass ...)" is not.
# A candidate passes when at least one entry starts with one of these.
SURFACE_LOCATION_PREFIXES = [
    "membrane",                 # bare "Membrane ..." = plasma membrane
    "cell membrane",
    "plasma membrane",
    "cell surface",
    "apical",
    "basolateral",
    "basal cell membrane",
    "lateral cell membrane",
    "microvillus",
    "cell projection",
    "cell junction",
    "synaptic",
    "postsynaptic",
    "presynaptic",
    "photoreceptor",
    "cilium membrane",
    "flagellum membrane",
    "acrosome membrane",
]
MIN_ECD_LENGTH = 25       # residues of contiguous non-cytoplasmic sequence

# --------------------------------------------------------------------------
# Step 5 -- normal-tissue safety
# --------------------------------------------------------------------------
# On-target/off-tumour toxicity is not symmetric across organs. Weight = how
# badly antigen-positive killing is tolerated in that tissue.
VITAL_ORGAN_WEIGHTS: dict[str, float] = {
    "heart muscle": 1.00,
    "brain": 1.00,
    "cerebral cortex": 1.00,
    "cerebellum": 1.00,
    "hippocampal formation": 1.00,
    "basal ganglia": 1.00,
    "hypothalamus": 1.00,
    "midbrain": 1.00,
    "medulla oblongata": 1.00,
    "pons": 1.00,
    "spinal cord": 1.00,
    "amygdala": 1.00,
    "thalamus": 1.00,
    "substantia nigra": 1.00,
    "white matter": 1.00,
    "liver": 0.95,
    "kidney": 0.95,
    "lung": 0.90,
    "bone marrow": 0.90,
    "pancreas": 0.85,
    "adrenal gland": 0.75,
    "thyroid gland": 0.70,
    "stomach": 0.70,
    "duodenum": 0.70,
    "small intestine": 0.70,
    "colon": 0.70,
    "rectum": 0.65,
    "esophagus": 0.65,
    "urinary bladder": 0.60,
    "skin": 0.60,
    "spleen": 0.55,
    "thymus": 0.50,
}
DEFAULT_ORGAN_WEIGHT = 0.40      # any normal tissue not named above
# No single organ may zero a target outright: a coefficient of exactly 0 would
# make the rest of the evidence unreadable, and toxicity that is managed in the
# clinic is a cost, not an absolute veto.
MAX_ORGAN_WEIGHT = 0.95
SAFETY_MISSING_DEFAULT = 0.70    # flagged neutral prior when HPA has no record
# nTPM at which a normal tissue counts as fully "on"; risk grows with
# log(nTPM), because the difference between 1 and 10 nTPM matters more than the
# difference between 200 and 400.
RNA_SATURATION_NTPM = 100.0
IHC_LEVEL_SCORE = {
    "not detected": 0.0,
    "low": 0.33,
    "medium": 0.66,
    "high": 1.0,
}

# --------------------------------------------------------------------------
# Step 6 -- composite score
# --------------------------------------------------------------------------
TUMOUR_QUALITY_WEIGHTS = {
    "specificity": 0.30,     # epithelial/malignant vs rest of the TME
    "intensity": 0.20,       # absolute expression level
    "uniformity": 0.20,      # fraction of epithelial cells positive
    "accessibility": 0.15,   # extracellular-domain accessibility
    "druggability": 0.15,    # Open Targets antibody tractability
}
TIER1_THRESHOLD = 0.55
TIER2_THRESHOLD = 0.35

# Antigen display requirement: an ADC or a CAR cannot engage a target that most
# tumour cells never put on their surface. Candidates below this epithelial
# detection rate are annotated and set aside before ranking rather than being
# allowed to win on a perfect safety score earned by not being expressed.
MIN_EPITHELIAL_DETECTION = 0.10

# --------------------------------------------------------------------------
# Step 1/7 -- pre-registered validation set (locked before any ranking)
# --------------------------------------------------------------------------
VALIDATION_POSITIVES = {
    "EGFR": "EGFR-directed ADC/bispecifics in LUAD (amivantamab; patritumab deruxtecan partner axis)",
    "ERBB2": "HER2 ADC, trastuzumab deruxtecan approved for HER2-mutant NSCLC",
    "MET": "telisotuzumab vedotin, c-Met ADC, approved/advanced in c-Met-high NSCLC",
    "TACSTD2": "TROP2 ADC (sacituzumab govitecan / datopotamab deruxtecan) in NSCLC",
    "CEACAM5": "CEACAM5 ADC tusamitamab ravtansine, phase 3 non-squamous NSCLC",
    "MSLN": "mesothelin ADC and CAR-T programmes in lung adenocarcinoma",
    "FOLR1": "folate receptor alpha ADC mirvetuximab class, NSCLC expansion cohorts",
    "CD276": "B7-H3 ADC (ifinatamab deruxtecan) and CAR-T in lung cancer",
    "ROR1": "ROR1 ADC zilovertamab vedotin, NSCLC cohorts",
}
NEGATIVE_CONTROLS = {
    "EPCAM": "pan-epithelial adhesion molecule; high in normal epithelium, catumaxomab-class toxicity",
    "CDH1": "E-cadherin, ubiquitous epithelial junction protein",
    "ATP1A1": "Na+/K+ ATPase alpha-1, housekeeping pump expressed in every tissue",
    "CGN": "cingulin, cytoplasmic tight-junction plaque protein -- must fail the topology gate",
}
RECALL_K = (10, 20)
# A negative control ranking inside this many positions means the safety or
# topology layer failed to price it.
NEGATIVE_CONTROL_FAIL_RANK = 50

# --------------------------------------------------------------------------
# Step 8 -- literature evidence
# --------------------------------------------------------------------------
N_LITERATURE_CANDIDATES = 12

# --------------------------------------------------------------------------
# Provenance helper
# --------------------------------------------------------------------------
def write_provenance(step: str, payload: dict) -> Path:
    """Persist a small JSON sidecar so the report never reconstructs numbers."""
    path = RESULTS / f"{step}_provenance.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
    return path


def env_offline() -> bool:
    return os.environ.get("ANTIGEN_PIPELINE_OFFLINE", "0") == "1"
