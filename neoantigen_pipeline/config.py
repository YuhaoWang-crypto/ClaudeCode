"""Configuration for the personalized neoantigen selection pipeline.

Every threshold and weight in here is a *published, editable default*, not a
reproduction of any proprietary vendor score.  The public Moderna/Merck
description of intismeran autogene (mRNA-4157/V940) fixes the *workflow*
(tumor/normal DNA + tumor RNA -> somatic variants -> HLA type -> presentation
prediction -> <=34 neoantigens -> one concatemeric mRNA in an LNP) but not the
scoring formula, training set, thresholds or construct rules.  Those are what
this module makes explicit and tunable.

Provenance of the defaults is documented in
`.claude/skills/neoantigen-selection/reference/scoring.md`.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import os

CACHE_DIR = os.environ.get(
    "NEOANTIGEN_CACHE", os.path.join(os.path.dirname(__file__), "data", "cache")
)

# --------------------------------------------------------------------------
# Hard gates: a candidate that fails any of these is dropped before scoring.
# Mirrors the public description "expressed, presented, tumor-specific".
# --------------------------------------------------------------------------


@dataclass
class Gates:
    min_tpm: float = 1.0                # gene must be expressed in the tumor
    max_rank_mhc1: float = 2.0          # NetMHCpan-EL %rank weak-binder cutoff
    max_rank_mhc2: float = 10.0         # NetMHCIIpan weak-binder cutoff
    min_dna_vaf: float = 0.05           # somatic variant support in tumor DNA
    min_ccf: float = 0.0                # cancer-cell fraction (0 = no clonality gate)
    require_novel_vs_self: bool = True  # mutant 9mer must not exist in the proteome
    drop_anchor_only: bool = False      # drop peptides whose mutation is at P2/POmega only
    exclude_genes: Tuple[str, ...] = ()  # e.g. HLA genes, IG/TR loci


# --------------------------------------------------------------------------
# Composite score weights. Sum is normalized internally, so relative size is
# what matters. Rationale per feature in reference/scoring.md.
# --------------------------------------------------------------------------


@dataclass
class Weights:
    """Defaults are presentation-dominant *because two benchmarks said so*.

    The original literature-balanced defaults gave presentation 0.30 and spread
    0.40 across agretopicity / dissimilarity / TCR prior / hydrophobicity. Both
    the presentation-controlled IEDB benchmark and the TESLA mirror (real T-cell
    assay labels, real negatives) found those four features at or below the
    random baseline, and the diluted composite scoring *below* the binding
    predictor alone. So the weight moved to where the evidence is.

    Expression, clonality and class-II support keep their weight: neither
    benchmark can evaluate them (an IEDB or TESLA peptide has no tumor RNA
    value and no CCF), and they encode facts about the tumor rather than
    guesses about the T-cell repertoire.

    See reference/benchmark.md for the numbers behind this.
    """

    presentation: float = 0.45    # MHC-I EL %rank (NetMHCpan-4.1)
    expression: float = 0.18      # tumor RNA TPM
    clonality: float = 0.12       # CCF from DNA VAF + purity + CN
    agretopicity: float = 0.08    # WT rank / MUT rank (Duan 2014)
    tcr_prior: float = 0.07       # Luksza-2017-style alignment to known immunogenic epitopes
    mhc2_support: float = 0.05    # a CD4 helper epitope in the same 25mer
    dissimilarity: float = 0.03   # BLOSUM-weighted distance from the self peptide
    hydrophobicity: float = 0.02  # Chowell 2015 / Wells 2020 TCR-contact hydrophobicity


# Named presets. `PRESENTATION_ONLY` is the single best-scoring setting on the
# TESLA mirror (AP 0.207 pooled, 0.266 mean per patient, 31/35 positives inside
# a 34-slot budget); use it when hit rate is the only objective and you have no
# expression or clonality data. `LITERATURE_BALANCED` is what this package
# shipped before the benchmarks were run -- kept so the change is reproducible,
# not because it is recommended.
PRESENTATION_ONLY = Weights(presentation=1.0, expression=0.0, clonality=0.0,
                            agretopicity=0.0, tcr_prior=0.0, mhc2_support=0.0,
                            dissimilarity=0.0, hydrophobicity=0.0)
LITERATURE_BALANCED = Weights(presentation=0.30, agretopicity=0.15, expression=0.15,
                              clonality=0.10, dissimilarity=0.10, tcr_prior=0.10,
                              hydrophobicity=0.05, mhc2_support=0.05)


@dataclass
class SelectionRules:
    """Constraints on the final set that goes into the construct."""

    max_neoantigens: int = 34       # public upper bound for mRNA-4157/V940
    min_neoantigens: int = 20       # public lower bound
    max_per_gene: int = 2           # avoid spending the payload on one gene
    min_alleles_covered: int = 4    # spread across the patient's class-I alleles
    per_allele_cap: Optional[int] = None  # optional hard cap per allele
    force_include_genes: Tuple[str, ...] = ()  # e.g. known drivers (BRAF, NRAS)
    prefer_clonal: bool = True


@dataclass
class ConstructRules:
    """mRNA concatemer design rules."""

    epitope_length: int = 25        # 25mer minigene, mutation centered at position 13
    linker: str = ""                # "" = direct fusion; e.g. "GPGPG" or "AAY"
    signal_peptide: str = ""        # optional secretion / MITD trafficking tag
    max_cds_nt: int = 4000          # payload budget for one LNP-formulated mRNA
    junction_scan_rank: float = 0.5  # flag junction peptides stronger than this %rank
    # Ordering is optimized over these lengths (each extra length multiplies the
    # number of predictions by ~1x, so the default is the dominant class-I length)...
    junction_cost_lengths: Tuple[int, ...] = (9,)
    # ...but the FINAL order is rescanned over all of these, so a binder created
    # at a length the optimizer did not see is still reported rather than missed.
    junction_scan_lengths: Tuple[int, ...] = (8, 9, 10, 11)
    host: str = "human"             # codon table for optimization
    avoid_sites: Tuple[str, ...] = ("BsaI", "BsmBI", "EcoRI", "BamHI", "NotI")
    max_homopolymer: int = 6
    gc_target: Tuple[float, float] = (40.0, 70.0)


@dataclass
class PatientConfig:
    """Everything patient-specific."""

    patient_id: str
    hla_class1: List[str] = field(default_factory=list)   # ["HLA-A*02:01", ...]
    hla_class2: List[str] = field(default_factory=list)   # ["HLA-DRB1*04:01", ...]
    tumor_purity: float = 0.7
    mhc1_lengths: Tuple[int, ...] = (8, 9, 10, 11)
    mhc2_lengths: Tuple[int, ...] = (15,)


@dataclass
class PipelineConfig:
    patient: PatientConfig
    gates: Gates = field(default_factory=Gates)
    weights: Weights = field(default_factory=Weights)
    selection: SelectionRules = field(default_factory=SelectionRules)
    construct: ConstructRules = field(default_factory=ConstructRules)
    cache_dir: str = CACHE_DIR
    mhc1_method: str = "netmhcpan_el"
    mhc2_method: str = "netmhciipan_el"

    def weight_dict(self) -> Dict[str, float]:
        w = self.weights.__dict__.copy()
        total = sum(w.values()) or 1.0
        return {k: v / total for k, v in w.items()}


# Common class-I supertype panel used when a patient HLA type is unavailable.
# NOT a substitute for real typing -- flagged as an assumption in the report.
DEMO_HLA_CLASS1 = [
    "HLA-A*01:01", "HLA-A*02:01",
    "HLA-B*07:02", "HLA-B*08:01",
    "HLA-C*07:01", "HLA-C*07:02",
]
DEMO_HLA_CLASS2 = ["HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01"]
