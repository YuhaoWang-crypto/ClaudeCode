"""
T1 - human tissue expression: organ-level matrix, tissue-specificity tau,
tissue-enriched gene sets.

Two preprocessing choices here change the answer more than anything downstream,
so they are made explicitly rather than by default:

1. GTEx v10 ships 68 columns, but 14 of them are the v10 cell-fraction /
   sub-region pilots (Liver_Hepatocyte, Pancreas_Islets, Stomach_Muscularis,
   ... ) and 2 are cultured cell lines. Keeping them would put four near-
   duplicate "liver" columns in the matrix. We keep the 52 classic bulk
   tissues.

2. GTEx has 13 brain sub-regions and only one testis column. tau (Yanai 2005)
   counts columns, so an identically brain-restricted gene looks 13x more
   "broadly expressed" than an identically testis-restricted one. Left
   uncorrected this ALONE would make brain look constrained and testis look
   specific - which is the headline result the pipeline is trying to test.
   We therefore collapse sub-regions to ORGANS (median across sub-regions)
   before computing tau, and keep cerebellum separate from the rest of the
   brain because its transcriptome genuinely is an outlier.

tau = sum_i (1 - x_i / x_max) / (n - 1), on log2(TPM+1), in [0,1];
0 = perfectly uniform, 1 = expressed in exactly one organ.
"""
import os
import gzip

import numpy as np
import pandas as pd

from .common import RAW, DERIVED, FIGDIR, apply_figure_style, log, to_tsv_gz

GTEX_GCT = os.path.join(RAW, "gtex_v10_median_tpm.gct.gz")
GTEX_URL = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
            "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")

# v10-only cell-fraction / sub-region pilots and cell lines - excluded.
DROP_COLS = [
    "Colon_Transverse_Mixed_Cell", "Colon_Transverse_Mucosa",
    "Colon_Transverse_Muscularis", "Liver_Hepatocyte", "Liver_Mixed_Cell",
    "Liver_Portal_Tract", "Pancreas_Acini", "Pancreas_Islets",
    "Pancreas_Mixed_Cell", "Small_Intestine_Terminal_Ileum_Mixed_Cell",
    "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate", "Stomach_Mixed_Cell",
    "Stomach_Mucosa", "Stomach_Muscularis",
    "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes",
]

# GTEx tissue -> organ. Anything unlisted maps to itself.
ORGAN_MAP = {
    "Adipose_Subcutaneous": "Adipose", "Adipose_Visceral_Omentum": "Adipose",
    "Artery_Aorta": "Artery", "Artery_Coronary": "Artery",
    "Artery_Tibial": "Artery",
    "Brain_Amygdala": "Brain", "Brain_Anterior_cingulate_cortex_BA24": "Brain",
    "Brain_Caudate_basal_ganglia": "Brain", "Brain_Cortex": "Brain",
    "Brain_Frontal_Cortex_BA9": "Brain", "Brain_Hippocampus": "Brain",
    "Brain_Hypothalamus": "Brain",
    "Brain_Nucleus_accumbens_basal_ganglia": "Brain",
    "Brain_Putamen_basal_ganglia": "Brain", "Brain_Substantia_nigra": "Brain",
    "Brain_Cerebellar_Hemisphere": "Cerebellum", "Brain_Cerebellum": "Cerebellum",
    "Brain_Spinal_cord_cervical_c-1": "Spinal_cord",
    "Cervix_Ectocervix": "Cervix", "Cervix_Endocervix": "Cervix",
    "Colon_Sigmoid": "Colon", "Colon_Transverse": "Colon",
    "Esophagus_Gastroesophageal_Junction": "Esophagus",
    "Esophagus_Mucosa": "Esophagus", "Esophagus_Muscularis": "Esophagus",
    "Heart_Atrial_Appendage": "Heart", "Heart_Left_Ventricle": "Heart",
    "Kidney_Cortex": "Kidney", "Kidney_Medulla": "Kidney",
    "Skin_Not_Sun_Exposed_Suprapubic": "Skin",
    "Skin_Sun_Exposed_Lower_leg": "Skin",
    "Small_Intestine_Terminal_Ileum": "Small_intestine",
    "Breast_Mammary_Tissue": "Breast", "Muscle_Skeletal": "Skeletal_muscle",
    "Nerve_Tibial": "Nerve_tibial", "Whole_Blood": "Blood",
    "Minor_Salivary_Gland": "Salivary_gland", "Adrenal_Gland": "Adrenal",
    "Fallopian_Tube": "Fallopian_tube",
}


def coding_universe():
    """Nuclear protein-coding genes - the only universe a standard-code dN/dS
    analysis is entitled to talk about.

    The exclusion that matters is the MITOCHONDRIAL genome. MT-encoded genes
    carry 28-65% of a tissue's TPM mass (65% in brain, 28% in testis), so
    leaving them in makes any expression-weighted statistic mostly a readout of
    how oxidative the tissue is. They also cannot be scored by T3 at all: the
    vertebrate mitochondrial code reassigns AGA/AGG to stop, ATA to Met and TGA
    to Trp, so applying the standard-code NG86 tables to them would be wrong
    rather than merely noisy. And mtDNA is uniparentally inherited, non-
    recombining and mutates roughly an order of magnitude faster - it is a
    different evolutionary regime, not a fast corner of the same one.
    """
    from .common import RAW, ENSEMBL_FTP, download, http_get, read_fasta
    path = os.path.join(RAW, "pep_homo_sapiens.fa.gz")
    if not os.path.exists(path):
        base = f"{ENSEMBL_FTP}/homo_sapiens/pep/"
        fn = [t.split("/")[-1] for t in http_get(base, timeout=300).split('"')
              if t.endswith(".pep.all.fa.gz")][0]
        download(base + fn, path)
    nuclear, mito = set(), set()
    for header, _ in read_fasta(path):
        fields = header.split()
        gene = biotype = loc = None
        for f in fields[1:]:
            if f.startswith("gene:"):
                gene = f.split(":", 1)[1].split(".")[0]
            elif f.startswith("gene_biotype:"):
                biotype = f.split(":", 1)[1]
            elif f.startswith(("chromosome:", "scaffold:")):
                loc = f
        if gene is None or biotype != "protein_coding":
            continue
        if loc and loc.split(":")[2] == "MT":
            mito.add(gene)
        else:
            nuclear.add(gene)
    return nuclear, mito


def _ensure_gtex():
    if not os.path.exists(GTEX_GCT):
        from .common import download
        log("downloading GTEx v10 median TPM")
        download(GTEX_URL, GTEX_GCT)
    return GTEX_GCT


def load_organ_matrix():
    """Return (organ TPM DataFrame indexed by unversioned ENSG, gene symbols)."""
    path = _ensure_gtex()
    with gzip.open(path, "rt") as fh:
        fh.readline()                       # '#1.2'
        fh.readline()                       # dims
        df = pd.read_csv(fh, sep="\t")
    df["Name"] = df["Name"].str.split(".").str[0]
    df = df[~df["Name"].duplicated(keep="first")]
    symbols = df.set_index("Name")["Description"]
    mat = df.set_index("Name").drop(columns=["Description"])
    mat = mat.drop(columns=[c for c in DROP_COLS if c in mat.columns])

    organs = pd.Series({c: ORGAN_MAP.get(c, c) for c in mat.columns})
    # median across sub-regions: robust to one odd sub-region (e.g. substantia
    # nigra has visibly lower library complexity than the other brain regions)
    organ_mat = mat.T.groupby(organs).median().T
    organ_mat = organ_mat.loc[organ_mat.max(axis=1) > 0]
    return organ_mat, symbols


def tau(mat, pseudocount=1.0):
    """Yanai tissue-specificity index on log2(TPM+1). NaN if a gene is silent."""
    x = np.log2(mat.to_numpy(dtype=float) + pseudocount)
    xmax = x.max(axis=1, keepdims=True)
    n = x.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (1.0 - x / xmax).sum(axis=1) / (n - 1)
    t[xmax.ravel() <= 0] = np.nan
    return pd.Series(t, index=mat.index, name="tau")


def enriched_sets(mat, fold=4.0, min_tpm=1.0):
    """HPA-style 'tissue enriched': top organ >= `fold` x second organ and
    >= `min_tpm`. Returns organ -> list of ENSG."""
    vals = mat.to_numpy(dtype=float)
    order = np.argsort(-vals, axis=1)
    top = vals[np.arange(len(vals)), order[:, 0]]
    second = vals[np.arange(len(vals)), order[:, 1]]
    keep = (top >= min_tpm) & (top >= fold * np.maximum(second, 1e-9))
    organ_of_top = np.array(mat.columns)[order[:, 0]]
    out = {}
    for organ in mat.columns:
        out[organ] = list(mat.index[keep & (organ_of_top == organ)])
    return out


def group_enriched_sets(mat, fold=4.0, min_tpm=1.0, max_group=5):
    """HPA-style 'group enriched': expression at least `fold` higher in a group
    of 2-`max_group` organs than in every other organ.

    This exists because the strict single-organ rule badly under-counts organs
    that have close transcriptional neighbours. Brain, cerebellum and spinal
    cord are three separate columns here, so a gene expressed throughout the
    CNS is never 4x above its runner-up and counts as enriched NOWHERE - which
    is why the single-organ rule finds only ~25 brain-enriched genes and would
    make every brain conclusion rest on that handful. The same applies to
    artery/heart and colon/small intestine.

    Returns gene -> tuple(organs), keeping the SMALLEST qualifying group.
    """
    vals = mat.to_numpy(dtype=float)
    order = np.argsort(-vals, axis=1)
    srt = np.take_along_axis(vals, order, axis=1)
    cols = np.array(mat.columns)
    out = {}
    for k in range(2, max_group + 1):
        if k >= srt.shape[1]:
            break
        lo_in_group = srt[:, k - 1]
        hi_outside = srt[:, k]
        ok = (lo_in_group >= min_tpm) & (lo_in_group >= fold * np.maximum(hi_outside, 1e-9))
        for i in np.flatnonzero(ok):
            gene = mat.index[i]
            if gene not in out:                      # smallest group wins
                out[gene] = tuple(sorted(cols[order[i, :k]]))
    return out


def build(force=False):
    """Compute and cache the per-gene expression feature table."""
    out = os.path.join(DERIVED, "t1_gene_expression_features.tsv.gz")
    organ_path = os.path.join(DERIVED, "t1_organ_median_tpm.tsv.gz")
    if os.path.exists(out) and not force:
        return pd.read_csv(out, sep="\t", index_col=0), pd.read_csv(
            organ_path, sep="\t", index_col=0)

    mat, symbols = load_organ_matrix()
    t = tau(mat)
    feat = pd.DataFrame({
        "symbol": symbols.reindex(mat.index),
        "tau": t,
        "max_tpm": mat.max(axis=1),
        "top_organ": mat.idxmax(axis=1),
        # expression LEVEL is the dominant known predictor of protein
        # evolutionary rate, so it is carried through as a covariate, not
        # as an afterthought
        "mean_tpm": mat.mean(axis=1),
        "n_organs_over_1tpm": (mat > 1).sum(axis=1),
    })
    enr = enriched_sets(mat)
    enriched_of = {}
    for organ, genes in enr.items():
        for g in genes:
            enriched_of[g] = organ
    feat["enriched_in"] = pd.Series(enriched_of).reindex(feat.index)

    grp = group_enriched_sets(mat)
    feat["group_enriched_in"] = pd.Series(
        {g: ",".join(o) for g, o in grp.items()}).reindex(feat.index)

    to_tsv_gz(mat, organ_path)
    to_tsv_gz(feat, out)
    log(f"T1: {mat.shape[0]} genes x {mat.shape[1]} organs -> {os.path.basename(out)}")
    return feat, mat


def make_figure(feat, mat):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))

    axes[0].hist(feat["tau"].dropna(), bins=60, color="#0072B2", alpha=0.85)
    axes[0].set_xlabel(r"tissue-specificity $\tau$")
    axes[0].set_ylabel("genes")
    axes[0].set_title(f"$\\tau$ over {mat.shape[1]} organs (GTEx v10)")

    counts = feat["enriched_in"].value_counts().head(15)
    axes[1].barh(counts.index[::-1], counts.values[::-1], color="#D55E00")
    axes[1].set_xlabel("tissue-enriched genes ($\\geq$4x next organ)")
    axes[1].set_title("organs by enriched-gene count")
    axes[1].grid(axis="y", alpha=0)

    p = os.path.join(FIGDIR, "t1_expression_overview.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report():
    print("=" * 72)
    print("T1 - human tissue expression, tau, tissue-enriched sets (GTEx v10)")
    print("=" * 72)
    feat, mat = build()
    print(f"organs kept                  : {mat.shape[1]}")
    print(f"genes with any expression    : {mat.shape[0]}")
    print(f"median tau                   : {feat['tau'].median():.3f}")
    print(f"genes with tau > 0.8         : {(feat['tau'] > 0.8).sum()}")
    print(f"tissue-enriched genes total  : {feat['enriched_in'].notna().sum()}")
    print("\ntop 10 organs by tissue-enriched gene count:")
    for organ, n in feat["enriched_in"].value_counts().head(10).items():
        sub = feat[feat["enriched_in"] == organ]
        print(f"    {organ:<18} {n:>5}   median tau {sub['tau'].median():.3f}")
    p = make_figure(feat, mat)
    print(f"\nfigure written to: {p}")
    return feat, mat


if __name__ == "__main__":
    report()
