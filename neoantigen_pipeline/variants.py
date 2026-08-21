"""Step 1-3: somatic variants -> annotated, tumor-specific, expressed variant table.

Input is a list of cBioPortal DETAILED mutation records (or any dict list with
the same keys); a MAF reader is provided for real patient data.

Output columns
--------------
sample, gene, entrez, protein_change, variant_class, ref_aa, alt_aa, aa_pos,
dna_vaf, tumor_depth, tpm, ccf, is_clonal
"""

from __future__ import annotations

import math
import re
from typing import Dict, Iterable, List, Optional

import pandas as pd

PROTEIN_CHANGE_RE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")

# Variant classes that can create a novel peptide sequence.
CODING_CLASSES = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}


def parse_protein_change(pc: str):
    """'R24C' -> ('R', 24, 'C'); returns (None, None, None) if not a simple SNV."""
    if not isinstance(pc, str):
        return None, None, None
    pc = pc.lstrip("p.")
    m = PROTEIN_CHANGE_RE.match(pc)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), m.group(3)


def from_cbioportal(records: Iterable[dict], sample_id: str) -> pd.DataFrame:
    rows = []
    for r in records:
        gene = (r.get("gene") or {}).get("hugoGeneSymbol") or r.get("hugoGeneSymbol")
        alt = r.get("tumorAltCount")
        ref = r.get("tumorRefCount")
        depth = (alt or 0) + (ref or 0)
        vaf = (alt / depth) if depth else float("nan")
        ref_aa, pos, alt_aa = parse_protein_change(r.get("proteinChange", ""))
        rows.append({
            "sample": r.get("sampleId", sample_id),
            "gene": gene,
            "entrez": r.get("entrezGeneId"),
            "protein_change": r.get("proteinChange"),
            "variant_class": r.get("mutationType"),
            "ref_aa": ref_aa,
            "aa_pos": pos,
            "alt_aa": alt_aa,
            "transcript": r.get("refseqMrnaId"),
            "chrom": r.get("chr"),
            "start": r.get("startPosition"),
            "ref_allele": r.get("referenceAllele"),
            "alt_allele": r.get("variantAllele"),
            "tumor_alt": alt,
            "tumor_depth": depth,
            "dna_vaf": vaf,
            "normal_alt": r.get("normalAltCount"),
        })
    return pd.DataFrame(rows)


def from_maf(path: str, sample_id: Optional[str] = None) -> pd.DataFrame:
    """Read a standard somatic MAF (the format a real tumor/normal caller emits)."""
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    cols = {c.lower(): c for c in df.columns}

    def col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None

    hgvsp = col("HGVSp_Short", "Protein_Change", "amino_acid_change")
    alt_c = col("t_alt_count")
    ref_c = col("t_ref_count")
    out = pd.DataFrame({
        "sample": df[col("Tumor_Sample_Barcode")] if col("Tumor_Sample_Barcode") else sample_id,
        "gene": df[col("Hugo_Symbol")],
        "entrez": df[col("Entrez_Gene_Id")] if col("Entrez_Gene_Id") else None,
        "protein_change": df[hgvsp].astype(str).str.lstrip("p.") if hgvsp else None,
        "variant_class": df[col("Variant_Classification")],
        "transcript": df[col("Transcript_ID")] if col("Transcript_ID") else None,
        "chrom": df[col("Chromosome")] if col("Chromosome") else None,
        "start": df[col("Start_Position", "Start_position")] if col("Start_Position", "Start_position") else None,
        "ref_allele": df[col("Reference_Allele")] if col("Reference_Allele") else None,
        "alt_allele": df[col("Tumor_Seq_Allele2")] if col("Tumor_Seq_Allele2") else None,
        "tumor_alt": df[alt_c] if alt_c else None,
        "tumor_depth": (df[alt_c] + df[ref_c]) if (alt_c and ref_c) else None,
    })
    out["dna_vaf"] = out["tumor_alt"] / out["tumor_depth"] if alt_c and ref_c else float("nan")
    parsed = out["protein_change"].apply(parse_protein_change)
    out["ref_aa"] = [p[0] for p in parsed]
    out["aa_pos"] = [p[1] for p in parsed]
    out["alt_aa"] = [p[2] for p in parsed]
    return out


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------

def add_expression(df: pd.DataFrame, expr: Dict[int, float],
                   rescale_to_tpm: bool = False) -> pd.DataFrame:
    """Join per-gene tumor RNA abundance.

    TCGA / cBioPortal `rna_seq_v2_mrna` values are RSEM normalized counts, which
    already sit on a TPM-like scale, so the default is to use them as-is. Set
    `rescale_to_tpm=True` only when the supplied dict really is the *whole*
    transcriptome and you want it renormalized to sum to 1e6; rescaling a
    gene subset silently inflates every value.
    """
    df = df.copy()
    if rescale_to_tpm and expr:
        total = sum(v for v in expr.values() if v and v > 0) or 1.0
        expr = {k: v * 1e6 / total for k, v in expr.items()}
    df["tpm"] = df["entrez"].map(lambda e: expr.get(int(e)) if pd.notna(e) else None)
    return df


def add_clonality(df: pd.DataFrame, purity: float = 0.7,
                  copy_number: Optional[Dict[str, float]] = None,
                  clonal_ccf: float = 0.8) -> pd.DataFrame:
    """CCF = VAF * (purity*CN_tumor + 2*(1-purity)) / (purity * multiplicity).

    With multiplicity=1 and CN=2 this reduces to the usual CCF ~ 2*VAF/purity.
    Supplying a per-gene copy number dict refines it.
    """
    df = df.copy()
    cn = copy_number or {}

    def _ccf(row):
        v = row.get("dna_vaf")
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return float("nan")
        cn_t = cn.get(row.get("gene"), 2.0)
        denom = purity * 1.0
        if denom <= 0:
            return float("nan")
        return min(1.5, v * (purity * cn_t + 2 * (1 - purity)) / denom)

    df["ccf"] = df.apply(_ccf, axis=1)
    df["is_clonal"] = df["ccf"] >= clonal_ccf
    return df


def apply_variant_gates(df: pd.DataFrame, gates) -> pd.DataFrame:
    """Hard gates that need no peptide-level information yet."""
    d = df.copy()
    d["gate_coding"] = d["variant_class"].isin(CODING_CLASSES)
    d["gate_parsed"] = d["aa_pos"].notna() & d["alt_aa"].notna() & (d["alt_aa"] != "*")
    d["gate_vaf"] = d["dna_vaf"].fillna(0) >= gates.min_dna_vaf
    d["gate_expr"] = d["tpm"].fillna(0) >= gates.min_tpm
    d["gate_ccf"] = d["ccf"].fillna(0) >= gates.min_ccf if "ccf" in d else True
    d["gate_gene"] = ~d["gene"].isin(gates.exclude_genes)
    gate_cols = [c for c in d.columns if c.startswith("gate_")]
    d["passes_variant_gates"] = d[gate_cols].all(axis=1)
    return d


def gate_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Waterfall of how many variants survive each gate (for the report)."""
    order = ["gate_coding", "gate_parsed", "gate_vaf", "gate_expr", "gate_ccf", "gate_gene"]
    rows, keep = [], pd.Series(True, index=df.index)
    rows.append({"step": "all somatic variants", "n": int(len(df))})
    for c in order:
        if c not in df:
            continue
        keep = keep & df[c].fillna(False)
        rows.append({"step": c, "n": int(keep.sum())})
    return pd.DataFrame(rows)
