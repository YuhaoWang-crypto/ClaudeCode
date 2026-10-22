"""Compute GISMO and GISMO-mis from the authors' released inputs.

Faithful Python port of GenerateScore/01_GISMOscores.R (Liao et al., Nature
Genetics 2026, doi:10.1038/s41588-026-02783-9; code doi:10.5281/zenodo.20547722).

GISMO (gene loss, lower = more constrained)
-------------------------------------------
  I_ij = 1 if TOGA calls gene i "one2zero" (ORF lost) in species j, else 0
  w_j  = d_j / D_max, d_j = human divergence time (Myr, TimeTree)
  point statistic   s_i = (1/N) * sum_j I_ij * w_j          (N = 462 species)
  The paper writes s_i = sum I w / sum w; the code divides by N. The two
  differ by the constant sum(w)/N, so ranks/deciles are identical.

  Uncertainty: p_hat is estimated per gene *per mammalian order* (lin.num),
  I*_ijk ~ Bernoulli(p_hat_i,order(j)), k = 1..K (paper K = 10,000).
  The released GISMO score is the 97.5th percentile ("upper.ci") of the
  simulated statistic; deciles are ntile(upper.ci, 10). Using the upper bound
  (like LOEUF) penalises genes whose low loss count is uncertain.

GISMO-mis (fixed missense, lower = more constrained)
----------------------------------------------------
  For each gene, per species counts m_ij (missense) and s_ij (synonymous)
  against a mammalian *consensus* codon sequence (not human).
  GISMO-mis = mean over species of (m_ij / s_ij)   (R: `mean.comb`; species
  with s_ij = 0 dropped). This is what the authors' released deciles and all
  downstream analyses use (verified: 99.98% decile identity). The R script
  also writes `mis.cal` = mean(m)/mean(s) - kept here as gismo_mis_ratio_of_means.
"""
import gzip

import numpy as np
import pandas as pd

from gismo import paths


# ----------------------------------------------------------------- loading
def load_species(root):
    """Species table: dedupe to the highest contig-N50 assembly per species."""
    gs = pd.read_csv(root / "GenerateScore/inputs/genestat_26-07-2023.txt", sep="\t")
    gs = gs.sort_values(["species", "contig_N50_bp"], ascending=[True, False])
    gs = gs.drop_duplicates("species")
    return gs


def load_loss_matrix(root, species):
    """Binary gene x species loss matrix (one2zero -> 1), species ordered by lin.num."""
    f = root / "GenerateScore/inputs/unmerged-species_combined_matrix_2023-07-26.tsv"
    raw = pd.read_csv(f, sep="\t", index_col=0, dtype=str)
    raw.columns = [c.rsplit("_", 1)[-1] for c in raw.columns]   # -> assembly name
    use = species[species.assembly_name.isin(raw.columns)].sort_values("lin.num", kind="stable")
    loss = (raw[use.assembly_name.tolist()] == "one2zero").to_numpy(np.uint8)
    # keep the full orthology class too (for descriptive stats)
    return pd.DataFrame(loss, index=raw.index, columns=use.assembly_name), use, raw


def load_mis_syn(root, species):
    def rd(name):
        df = pd.read_csv(root / f"GenerateScore/inputs/{name}", sep="\t",
                         index_col=0, na_values=["None"], low_memory=False)
        return df[df.index.isin(species.assembly_name)].sort_index()
    mis = rd("consensus_v2_missense-counts.gz")
    syn = rd("consensus_v2_synonymous-counts.gz")
    common = mis.columns.intersection(syn.columns)
    return mis[common].astype(float), syn.loc[mis.index, common].astype(float)


def load_gene_annot(root):
    g = pd.read_csv(root / "partitioned-heritability/inputs/"
                    "gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3.genes", sep="\t")
    g = g.drop_duplicates("gene_id")
    return g.set_index("gene_id")[["gene_name", "seqnames"]]


# ----------------------------------------------------------------- GISMO
def gismo(loss, species, n_iter=2000, seed=1337, chunk=40_000):
    """Return DataFrame with point stat, simulated CI and deciles.

    n_iter=2000 (paper: 10,000) keeps runtime ~1 min; rank agreement with
    the authors' deciles is checked in validate_against_author().
    """
    I = loss.to_numpy(np.float32)                       # G x N
    G, N = I.shape
    w = (species["div.time"].to_numpy(np.float64) / species["div.time"].max()).astype(np.float32)
    grp = species["lin.num"].to_numpy()

    point = (I * w).mean(1)

    # p_hat per gene per order, broadcast back to species
    P = np.zeros_like(I)
    for g in np.unique(grp):
        m = grp == g
        P[:, m] = I[:, m].mean(1, keepdims=True)

    # analytic normal approximation (cheap cross-check of the simulation)
    mu = (P * w).sum(1) / N
    sd = np.sqrt((P * (1 - P) * w ** 2).sum(1)) / N

    # sparse Monte-Carlo: only species with 0 < p < 1 are random
    rng = np.random.default_rng(seed)
    fixed = ((P >= 1) * w).sum(1) / N                  # p==1 always lost
    gi, sj = np.nonzero((P > 0) & (P < 1))
    p_nz, w_nz = P[gi, sj], w[sj] / N
    sims = np.zeros((G, n_iter), np.float32)
    starts = np.r_[0, np.cumsum(np.bincount(gi, minlength=G))]
    # process whole genes per chunk
    g0 = 0
    while g0 < G:
        g1 = g0
        while g1 < G and starts[g1 + 1] - starts[g0] <= chunk:
            g1 += 1
        g1 = max(g1, g0 + 1)
        a, b = starts[g0], starts[g1]
        if b > a:
            draws = rng.random((n_iter, b - a), dtype=np.float32) < p_nz[a:b]
            contrib = draws * w_nz[a:b]                  # n_iter x nnz
            # segment-sum per gene via cumulative sums (gi is sorted)
            seg = gi[a:b] - g0
            cs = np.cumsum(contrib, axis=1, dtype=np.float64)
            ends = np.searchsorted(seg, np.arange(g1 - g0), side="right") - 1
            begs = np.searchsorted(seg, np.arange(g1 - g0), side="left") - 1
            hi = np.where(ends >= 0, cs[:, np.clip(ends, 0, None)], 0.0)
            lo = np.where(begs >= 0, cs[:, np.clip(begs, 0, None)], 0.0)
            empty = ends == begs
            out = np.where(empty, 0.0, hi - lo)
            sims[g0:g1] = out.T
        g0 = g1
    sims += fixed[:, None]
    q = np.quantile(sims, [0.025, 0.5, 0.975], axis=1).T

    df = pd.DataFrame({
        "gene_id": loss.index,
        "n_loss": I.sum(1).astype(int),
        "gismo_point": point,
        "gismo_lower": q[:, 0], "gismo_median": q[:, 1], "gismo": q[:, 2],
        "gismo_normal_upper": mu + 1.959964 * sd,
    })
    df["gismo_decile"] = ntile(df["gismo"], 10)
    return df


def ntile(x, n):
    """dplyr::ntile: rank (ties broken by order), then ceil(rank*n/len)."""
    r = pd.Series(x).rank(method="first").to_numpy()
    return np.floor((r - 1) * n / len(r)).astype(int) + 1


# ----------------------------------------------------------------- GISMO-mis
def gismo_mis(mis, syn):
    mean_mis, mean_syn = mis.mean(axis=0), syn.mean(axis=0)
    ratio = (mis / syn).replace(np.inf, np.nan)
    df = pd.DataFrame({
        "transcript_id": mis.columns,
        "mean_mis": mean_mis.values, "mean_syn": mean_syn.values,
        "gismo_mis": ratio.mean(axis=0).values,                  # released metric (mean.comb)
        "gismo_mis_ratio_of_means": (mean_mis / mean_syn).values,  # R mis.cal
        "n_species_mis": mis.notna().sum(axis=0).values,
    })
    df["gismo_mis_decile"] = ntile(df["gismo_mis"].fillna(np.inf), 10)
    df.loc[df.gismo_mis.isna(), "gismo_mis_decile"] = np.nan
    return df


# ----------------------------------------------------------------- author truth
# NOTE: partitioned-heritability/inputs/decile/d*/gene-set-file look like GISMO
# deciles but are GISMO-mis deciles (README; 100% identical to `dec` below).
# The archive does not ship final GISMO scores (those are Supplementary Table 2),
# so GISMO is validated by reproducing the paper's reported numbers instead.
def author_gismo_mis(root):
    f = root / "partitioned-heritability/inputs/gismo-mis_consensus_genenames_01-08-2023.txt"
    a = pd.read_csv(f, sep="\t")
    return a.rename(columns={"genes": "transcript_id", "dec": "author_gismo_mis_decile",
                             "genename": "symbol_mis", "mean.comb": "author_gismo_mis"})


# ----------------------------------------------------------------- build
def build(n_iter=2000, out=paths.SCORES_TSV):
    root = paths.zenodo_root()
    species = load_species(root)
    loss, sp_used, _ = load_loss_matrix(root, species)
    print(f"[gismo] loss matrix {loss.shape[0]} genes x {loss.shape[1]} species")
    g = gismo(loss, sp_used, n_iter=n_iter)
    ann = load_gene_annot(root)
    g["symbol"] = g.gene_id.map(ann.gene_name)
    g["chrom"] = g.gene_id.map(ann.seqnames).astype(str)

    mis, syn = load_mis_syn(root, sp_used)
    m = gismo_mis(mis, syn)
    a = author_gismo_mis(root)[["transcript_id", "symbol_mis", "author_gismo_mis",
                                "author_gismo_mis_decile"]]
    m = m.merge(a, on="transcript_id", how="left")
    print(f"[gismo] GISMO-mis {len(m)} transcripts x {mis.shape[0]} species")

    # join GISMO (ENSG) and GISMO-mis (ENST) on gene symbol
    m_by_sym = m.dropna(subset=["symbol_mis"]).drop_duplicates("symbol_mis")
    full = g.merge(m_by_sym, left_on="symbol", right_on="symbol_mis", how="outer")
    full["symbol"] = full["symbol"].fillna(full["symbol_mis"])
    full = full.drop(columns=["symbol_mis"])
    cols = ["symbol", "gene_id", "chrom", "transcript_id",
            "n_loss", "gismo", "gismo_decile", "gismo_point", "gismo_lower",
            "gismo_median", "gismo_normal_upper",
            "gismo_mis", "gismo_mis_decile", "mean_mis", "mean_syn",
            "gismo_mis_ratio_of_means", "n_species_mis", "author_gismo_mis",
            "author_gismo_mis_decile"]
    full = full[cols].sort_values(["gismo", "gismo_mis"], na_position="last")
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        full.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    print(f"[gismo] wrote {out} ({len(full)} rows)")
    return full


def load_scores(path=paths.SCORES_TSV):
    return pd.read_csv(path, sep="\t", dtype={"chrom": str})


if __name__ == "__main__":
    import sys
    build(n_iter=int(sys.argv[1]) if len(sys.argv) > 1 else 2000)
