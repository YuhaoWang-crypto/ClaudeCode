import os, numpy as np, pandas as pd
from scipy import stats
DATA = os.environ.get("MQ_DATA", "data"); WORK = os.environ.get("MQ_WORK", "work"); OUT = "results"
os.makedirs(WORK, exist_ok=True)

def int_norm(v):                      # rank-based inverse normal transform
    r = stats.rankdata(v); return stats.norm.ppf((r - 0.5) / len(v))

def residualize(M, C):                # remove covariates (rows = samples)
    C = np.column_stack([np.ones(len(C)), C]); return M - C @ np.linalg.lstsq(C, M, rcond=None)[0]

def pcs(G, k=3, seed=0):              # genotype PCs from a random SNP subset
    rng = np.random.default_rng(seed); X = G[:, rng.choice(G.shape[1], min(20000, G.shape[1]), replace=False)].astype(float)
    X = X - X.mean(0); X = X[:, X.std(0) > 0]; X /= X.std(0)
    return np.linalg.svd(X, full_matrices=False)[0][:, :k]

def load():
    """Genotypes (samples x SNPs, 0/1/2 dosage of bim allele 1), bim, expression (samples x genes), gene table, sample sets."""
    from pandas_plink import read_plink1_bin
    g = read_plink1_bin(f"{WORK}/geno.bed", verbose=False)
    G = g.values.astype(np.float32); bim = pd.DataFrame({"snp": g.snp.values, "chr": g.chrom.values.astype(str), "pos": g.pos.values.astype(int),
                                                         "a_counted": g.a1.values, "a_other": g.a0.values})
    ids = g.iid.values.astype(str); col = np.nanmean(G, 0); ix = np.where(np.isnan(G)); G[ix] = col[ix[1]]
    E = pd.read_csv(f"{WORK}/expr.tsv", sep="\t", index_col=0); genes = pd.read_csv(f"{WORK}/genes.tsv", sep="\t", index_col=0)
    sets = pd.read_csv(f"{WORK}/sets.tsv", sep="\t", index_col=0).loc[ids]
    return G, bim, E.loc[ids], genes, sets
