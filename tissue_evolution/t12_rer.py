"""
T12 - relative evolutionary rates (RERconverge) across the Zoonomia mammals.

Everything before this module compares rates between PAIRS of species, which
cannot say where on the tree a rate change happened. This one works branch by
branch across 400+ mammals, which is what "did this gene speed up in the
lineages that evolved big brains" actually requires.

Data
----
* Codon alignments: TOGA (Kirilenko et al., Science 2023), human hg38 reference,
  `MultipleCodonAlignments` - per-gene MACSE2 alignments over ~540 mammalian
  assemblies with inactivating mutations masked, built for selection screens.
* Topology: UCSC hg38 470-way (`scientificNames.nh`). Only the TOPOLOGY is used;
  every branch length is re-estimated per gene from the alignment.
* Trait: adult brain and body mass for 1349 mammals (Gonzalez-Suarez et al.,
  figshare 12867305). Relative brain size = residual of log brain on log body
  mass, so the scan is about ENCEPHALISATION, not about being a large animal.

Method (a reimplementation, not the R package - differences noted)
-----------------------------------------------------------------
1. Per gene, translate the codon alignment and compute all pairwise
   Poisson-corrected amino-acid distances, vectorised as a Gram matrix.
2. Fit branch lengths on the FIXED topology by least squares. The design matrix
   maps each branch (= each split of the tip set) onto the pairs whose path
   crosses it; it is identical for every gene, so its pseudo-inverse is computed
   ONCE and each gene costs a single matrix-vector product. This replaces
   per-gene ML tree inference, which is what makes 17k genes x 400 taxa
   tractable here; it is less precise than RAxML/PAML branch lengths, and that
   is the main methodological gap versus a published RERconverge analysis.
3. RER: sqrt-transform the gene x branch matrix, then for each BRANCH regress
   its values across genes on each gene's overall rate, and keep the residuals.
   This removes both "this branch is long because it spans a lot of time" and
   "this gene is fast everywhere", leaving rate relative to expectation.
4. Trait per branch: continuous tip values, ancestral states reconstructed under
   Brownian motion by the standard two-pass algorithm, branch value = mean of
   its two endpoints.
5. Per gene, correlate RER across branches with the branch trait
   (Spearman), then Benjamini-Hochberg across genes.

Result, stated up front
-----------------------
The scan returns 230 genes at BH q<0.05 out of 4195, with a coherent-looking
top list (RNA-binding and splicing factors). NONE OF IT SURVIVES THE NULL.
Traits simulated under Brownian motion on the same tree - pure noise carrying
only the tree's phylogenetic autocorrelation - yield a MEDIAN of 194 hits at
the same threshold, 90th percentile 430. The real trait's 230 has an empirical
p of 0.395.

The reason is that BH corrects across GENES while the dependence that inflates
the p-values is across BRANCHES: encephalisation is concentrated in primates and
cetaceans, so any gene with a rate shift anywhere in those clades correlates
with it, and the per-gene t-test treats 239 phylogenetically nested branches as
239 independent observations. This is the same failure mode the RERconverge
authors addressed with "permulations" (Saputra et al., Genetics 2021); the BM
simulation here is that idea. A convergence scan reported without a null of this
kind is reporting its own tree.

What this can and cannot conclude
---------------------------------
A significant correlation says a gene's rate co-varies with the trait across
lineages. It does NOT say the gene caused the trait, and it does not separate
positive selection from relaxed constraint - trait LOSS relaxes constraint and
accelerates genes exactly as adaptation does. Direction of the correlation is
reported so that acceleration and deceleration can be told apart, which is the
minimum needed before anyone reads a hit as adaptive.

And a null result here is not proof of absence. This run is under-powered
relative to a published RERconverge analysis in three specific ways, all
measurable: 4195 genes rather than ~19000 (the complete-core-species
requirement, see CORE_PRESENCE); 240 branches from 121 species rather than
~480 from 241; and least-squares rather than maximum-likelihood branch lengths.
Each of those costs power, so "no signal above a phylogenetic null" here means
exactly that, not "no signal exists".
"""
import os
import re
import gzip
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE, log,
                     to_tsv_gz)

ALN_DIR = os.path.join(RAW, "toga_aln")
TREE = os.path.join(RAW, "tree_470way.nh")
OVERVIEW = os.path.join(RAW, "toga_overview.tsv")
BRAIN = os.path.join(RAW, "brain_complete.csv")

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
CODON_TABLE = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
BASES = "TCAG"
_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA_OF = dict(zip(_CODONS, CODON_TABLE))

# Core-species threshold. There is no free choice here: NO species is present
# in 100% of TOGA alignments, so a fixed complete core trades species against
# genes. Measured on a 300-gene probe:
#     presence >=0.99 ->   5 species, 97% of genes retained
#     presence >=0.97 ->  84 species, 70% of genes retained
#     presence >=0.95 -> 206 species, 33% of genes retained
#     presence >=0.90 -> 359 species,  6% of genes retained
# RER power per gene comes from the NUMBER OF BRANCHES, so species are worth
# more than genes here, up to the point where gene retention collapses. 0.95 is
# that knee. The retained genes are NOT a random third - genes missing species
# are the lineage-specific and fast-diverging ones, the same bias T4 found for
# 1:1 orthologues, and it is reported rather than assumed away.
CORE_PRESENCE = 0.95
CORE_PROBE = 300          # alignments sampled to measure presence
MIN_SITES = 100           # aligned, ungapped codons required for a pairwise distance


# --------------------------------------------------------------------------
# species / tree
# --------------------------------------------------------------------------
def assembly_to_species():
    t = pd.read_csv(OVERVIEW, sep="\t")
    return dict(zip(t["Assembly name"], t["Species"].str.replace(" ", "_")))


def tree_tips():
    nh = open(TREE).read()
    return set(re.findall(r"[\(,]\s*([A-Za-z][A-Za-z0-9_]+):", nh))


def load_tree(keep):
    """Newick -> Bio.Phylo tree pruned to `keep`, with polytomies left alone."""
    from Bio import Phylo
    import io
    tree = Phylo.read(io.StringIO(open(TREE).read()), "newick")
    for tip in [t for t in tree.get_terminals() if t.name not in keep]:
        tree.prune(tip)
    return tree


def tree_splits(tree, tips):
    """Every branch as a split of the tip set, canonicalised and deduplicated.

    Splits, not edges: the two edges either side of the root induce the same
    split and only their sum is identifiable from distances, so treating them
    as separate branches would make the design matrix rank-deficient.
    """
    idx = {t: i for i, t in enumerate(tips)}
    n = len(tips)
    seen, splits = set(), []
    for clade in tree.find_clades():
        below = [idx[t.name] for t in clade.get_terminals() if t.name in idx]
        if not below or len(below) == n:
            continue
        mask = np.zeros(n, dtype=bool)
        mask[below] = True
        if mask[0]:                      # canonical side excludes tip 0
            mask = ~mask
        key = mask.tobytes()
        if key in seen:
            continue
        seen.add(key)
        splits.append(mask)
    return np.array(splits)              # (n_branches, n_tips)


def design_matrix(splits):
    """A[pair, branch] = 1 if the branch lies on the path between the pair."""
    n_tips = splits.shape[1]
    iu = np.triu_indices(n_tips, k=1)
    # a branch separates i and j iff exactly one of them is on its far side
    A = (splits[:, iu[0]] ^ splits[:, iu[1]]).T.astype(np.float32)
    return A, iu


# --------------------------------------------------------------------------
# per-gene distances and branch lengths
# --------------------------------------------------------------------------
def read_alignment(path, asm2sp):
    """gz codon FASTA -> {species: amino-acid string}. REFERENCE is human."""
    out = {}
    with gzip.open(path, "rt") as fh:
        name, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                head = line[1:].split()[0]
                if head == "REFERENCE":
                    name = "Homo_sapiens"
                elif head.startswith("vs_"):
                    name = asm2sp.get(head[3:])
                else:
                    name = None
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            out[name] = "".join(chunks)
    return {k: v for k, v in out.items() if k}


def translate_matrix(seqs, order):
    """-> uint8 matrix (species x codons); 255 = gap/ambiguous/stop."""
    L = min(len(s) for s in seqs.values()) // 3
    M = np.full((len(order), L), 255, dtype=np.uint8)
    for r, sp in enumerate(order):
        s = seqs[sp].upper().replace("U", "T")
        for c in range(L):
            aa = AA_OF.get(s[3 * c:3 * c + 3])
            if aa is not None and aa != "*":
                M[r, c] = AA_IDX[aa]
    return M


def pairwise_distances(M, iu):
    """Poisson-corrected amino-acid distance for every pair, vectorised."""
    n, L = M.shape
    valid = (M != 255)
    onehot = np.zeros((n, L * 20), dtype=np.float32)
    rows, cols = np.nonzero(valid)
    onehot[rows, cols * 20 + M[rows, cols]] = 1.0
    matches = onehot @ onehot.T
    covered = valid.astype(np.float32) @ valid.astype(np.float32).T
    with np.errstate(invalid="ignore", divide="ignore"):
        p = 1.0 - matches / covered
    p = p[iu]
    cov = covered[iu]
    # Poisson correction; p >= 0.95 is treated as saturated rather than clipped
    # silently to a finite value that would look like a measurement
    d = np.full(p.shape, np.nan, dtype=np.float64)
    ok = (cov >= MIN_SITES) & (p < 0.95)
    d[ok] = -np.log(1.0 - p[ok])
    return d


_GLOBAL = {}


def _init_worker(order, pinv, iu, asm2sp):
    _GLOBAL.update(order=order, pinv=pinv, iu=iu, asm2sp=asm2sp)


def _gene_branches(fname):
    """Least-squares branch lengths for one gene, or None if unusable."""
    order, pinv, iu, asm2sp = (_GLOBAL["order"], _GLOBAL["pinv"],
                               _GLOBAL["iu"], _GLOBAL["asm2sp"])
    try:
        seqs = read_alignment(os.path.join(ALN_DIR, fname), asm2sp)
        if not all(sp in seqs for sp in order):
            return None
        M = translate_matrix({sp: seqs[sp] for sp in order}, order)
        if M.shape[1] < MIN_SITES:
            return None
        d = pairwise_distances(M, iu)
        if not np.isfinite(d).all():
            return None
        b = pinv @ d
        return fname, np.clip(b, 0.0, None).astype(np.float32), M.shape[1]
    except Exception:                                     # noqa: BLE001
        return None


def choose_core(asm2sp, tips, sample=CORE_PROBE, seed=0):
    """Species present in >= CORE_PRESENCE of a random sample of alignments."""
    files = sorted(f for f in os.listdir(ALN_DIR) if f.endswith(".fasta.gz"))
    rng = np.random.default_rng(seed)
    probe = rng.choice(files, size=min(sample, len(files)), replace=False)
    counts = {}
    for f in probe:
        for sp in read_alignment(os.path.join(ALN_DIR, f), asm2sp):
            counts[sp] = counts.get(sp, 0) + 1
    keep = sorted(sp for sp, c in counts.items()
                  if sp in tips and c >= CORE_PRESENCE * len(probe))
    # measure, on the same probe, what fraction of genes will survive the core
    kept = set(keep)
    retention = float(np.mean([kept <= set(read_alignment(
        os.path.join(ALN_DIR, f), asm2sp)) for f in probe]))
    log(f"core={len(keep)} species at presence>={CORE_PRESENCE:.2f}; "
        f"expected gene retention {retention:.0%}")
    if "Homo_sapiens" in keep:
        keep.remove("Homo_sapiens")
        keep.insert(0, "Homo_sapiens")          # tip 0 anchors split canonicalisation
    return keep, files


def build_matrix(processes=4, limit=None, force=False):
    """Gene x branch matrix of least-squares branch lengths."""
    out = os.path.join(DERIVED, "t12_branch_lengths.npz")
    if os.path.exists(out) and not force:
        z = np.load(out, allow_pickle=True)
        return (z["B"], list(z["genes"]), list(z["order"]),
                z["splits"], list(z["n_codons"]))

    asm2sp = assembly_to_species()
    tips = tree_tips()
    order, files = choose_core(asm2sp, tips)
    log(f"core species: {len(order)} (present in >={CORE_PRESENCE:.0%} of probed genes)")

    tree = load_tree(set(order))
    present = [t.name for t in tree.get_terminals()]
    order = [s for s in order if s in set(present)]
    splits = tree_splits(tree, order)
    A, iu = design_matrix(splits)
    rank = int(np.linalg.matrix_rank(A))
    log(f"branches (identifiable splits): {splits.shape[0]}, pairs: {A.shape[0]}, "
        f"design rank: {rank}")
    if rank < splits.shape[0]:
        log(f"  {splits.shape[0] - rank} branch(es) not separately identifiable "
            f"(root edge / polytomy); pinv gives the minimum-norm split")
    pinv = np.linalg.pinv(A).astype(np.float32)

    if limit:
        files = files[:limit]
    log(f"estimating branch lengths for {len(files)} genes on {processes} processes")
    genes, rows, ncod = [], [], []
    with Pool(processes, initializer=_init_worker,
              initargs=(order, pinv, iu, asm2sp)) as pool:
        for i, r in enumerate(pool.imap_unordered(_gene_branches, files,
                                                  chunksize=8), 1):
            if r is not None:
                genes.append(r[0])
                rows.append(r[1])
                ncod.append(r[2])
            if i % 2000 == 0:
                log(f"  {i}/{len(files)} processed, {len(genes)} usable")
    B = np.vstack(rows)
    log(f"branch-length matrix: {B.shape[0]} genes x {B.shape[1]} branches "
        f"({len(genes)}/{len(files)} = {len(genes)/len(files):.0%} of genes retained)")
    np.savez_compressed(out, B=B, genes=np.array(genes), order=np.array(order),
                        splits=splits, n_codons=np.array(ncod))
    return B, genes, order, splits, ncod


# --------------------------------------------------------------------------
# RER
# --------------------------------------------------------------------------
def compute_rer(B, min_branch=1e-6):
    """RERconverge-style residuals: sqrt transform, then per-branch regression
    on each gene's overall rate."""
    X = np.sqrt(np.clip(B, 0, None))
    gene_rate = X.mean(axis=1)
    gr = (gene_rate - gene_rate.mean()) / (gene_rate.std() + 1e-12)
    RER = np.empty_like(X)
    for j in range(X.shape[1]):
        y = X[:, j]
        slope = float((gr * (y - y.mean())).mean())      # gr is standardised
        RER[:, j] = y - (y.mean() + slope * gr)
    sd = RER.std(axis=0, keepdims=True)
    sd[sd < min_branch] = 1.0
    return RER / sd


# --------------------------------------------------------------------------
# trait
# --------------------------------------------------------------------------
def relative_brain_size(order):
    """Residual of log10 brain mass on log10 body mass, over the core species."""
    b = pd.read_csv(BRAIN).dropna(subset=["Adult_brain_mass_g", "Adult_body_mass_g"])
    b["sp"] = b["binomial_Hedges_final"].str.replace(" ", "_")
    b = b.drop_duplicates("sp").set_index("sp")
    b = b[b.index.isin(order)]
    x = np.log10(b["Adult_body_mass_g"].to_numpy())
    y = np.log10(b["Adult_brain_mass_g"].to_numpy())
    slope, intercept = np.polyfit(x, y, 1)
    res = y - (slope * x + intercept)
    return pd.Series(res, index=b.index), float(slope)


def bm_ancestral(tree, tips, trait):
    """Brownian-motion ancestral states, standard two-pass algorithm.

    Tips without a trait value contribute no information (infinite variance)
    rather than being dropped, so the tree topology is preserved.
    """
    down_mu, down_v = {}, {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            v = trait.get(clade.name, np.nan)
            down_mu[id(clade)] = v if np.isfinite(v) else 0.0
            down_v[id(clade)] = 0.0 if np.isfinite(v) else np.inf
        else:
            num = den = 0.0
            for ch in clade.clades:
                t = ch.branch_length or 1e-6
                vv = down_v[id(ch)] + t
                if np.isfinite(vv) and vv > 0:
                    num += down_mu[id(ch)] / vv
                    den += 1.0 / vv
            down_mu[id(clade)] = num / den if den > 0 else 0.0
            down_v[id(clade)] = 1.0 / den if den > 0 else np.inf

    state = {id(tree.root): down_mu[id(tree.root)]}
    parent = {}
    for clade in tree.find_clades():
        for ch in clade.clades:
            parent[id(ch)] = clade
    for clade in tree.find_clades(order="preorder"):
        if id(clade) in state:
            continue
        p = parent[id(clade)]
        t = clade.branch_length or 1e-6
        # combine the parent's state with this subtree's own evidence
        vv = down_v[id(clade)]
        if not np.isfinite(vv):
            state[id(clade)] = state[id(p)]
        else:
            w_p, w_d = 1.0 / t, 1.0 / max(vv, 1e-9)
            state[id(clade)] = (w_p * state[id(p)] + w_d * down_mu[id(clade)]) / (w_p + w_d)
    return state, parent


def branch_trait(tree, order, splits, trait):
    """Trait value per identifiable split = mean of the branch's two endpoints."""
    state, parent = bm_ancestral(tree, order, trait)
    idx = {t: i for i, t in enumerate(order)}
    n = len(order)
    by_split = {}
    for clade in tree.find_clades():
        if id(clade) not in parent:
            continue
        below = [idx[t.name] for t in clade.get_terminals() if t.name in idx]
        if not below or len(below) == n:
            continue
        mask = np.zeros(n, dtype=bool)
        mask[below] = True
        if mask[0]:
            mask = ~mask
        val = 0.5 * (state[id(clade)] + state[id(parent[id(clade)])])
        by_split.setdefault(mask.tobytes(), []).append(val)
    return np.array([np.mean(by_split.get(s.tobytes(), [np.nan]))
                     for s in splits])


# --------------------------------------------------------------------------
# scan
# --------------------------------------------------------------------------
def _bh(p):
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    prev = 1.0
    qq = np.empty(n)
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, pv[i] * n / k)
        qq[i] = prev
    q[ok] = qq
    return q


def rer_ranks(RER, keep):
    """Standardised rank of each gene's RER across the usable branches.

    Precomputed once because it does not depend on the trait - which turns each
    null simulation into a single matrix-vector product instead of a rescan.
    """
    R = RER[:, keep]
    Rr = np.apply_along_axis(stats.rankdata, 1, R)
    Rr = (Rr - Rr.mean(axis=1, keepdims=True)) / Rr.std(axis=1, keepdims=True)
    return Rr


def _rho_from_ranks(Rr, t):
    tr = stats.rankdata(t)
    tr = (tr - tr.mean()) / tr.std()
    return Rr @ tr / len(tr)


def simulate_bm_trait(tree, order, rng):
    """Tip values under Brownian motion on the same tree.

    This is the null that matters. Encephalisation is strongly phylogenetically
    clustered (primates and cetaceans at one end), so ANY gene with a lineage-
    specific rate shift in those clades correlates with it. Permuting the trait
    across tips would destroy that clustering and give a null that is far too
    permissive. Simulating under BM keeps the phylogenetic covariance and asks
    the honest question: how many hits does a trait with this much phylogenetic
    structure produce when it has nothing to do with the genes?
    """
    vals = {id(tree.root): 0.0}
    parent = {}
    for clade in tree.find_clades():
        for ch in clade.clades:
            parent[id(ch)] = clade
    for clade in tree.find_clades(order="preorder"):
        if id(clade) in vals:
            continue
        t = clade.branch_length or 1e-6
        vals[id(clade)] = vals[id(parent[id(clade)])] + rng.normal(0.0, np.sqrt(t))
    return {tip.name: vals[id(tip)] for tip in tree.get_terminals()
            if tip.name in set(order)}


def null_calibration(tree, order, splits, Rr, keep, n_sim=200, seed=0,
                     alpha=0.05):
    """How many genes pass BH q<alpha under phylogenetically-structured nulls?"""
    rng = np.random.default_rng(seed)
    n = int(keep.sum())
    counts = []
    for _ in range(n_sim):
        sim = simulate_bm_trait(tree, order, rng)
        bt = branch_trait(tree, order, splits, sim)
        t = bt[keep]
        if not np.isfinite(t).all():
            continue
        rho = _rho_from_ranks(Rr, t)
        with np.errstate(divide="ignore", invalid="ignore"):
            tstat = rho * np.sqrt((n - 2) / (1 - rho ** 2))
        p = 2 * stats.t.sf(np.abs(tstat), df=n - 2)
        counts.append(int((_bh(p) < alpha).sum()))
    return np.array(counts)


def scan(RER, genes, bt, min_branches=50):
    """Spearman between each gene's RER and the branch trait."""
    keep = np.isfinite(bt)
    log(f"branches usable for the trait scan: {keep.sum()}/{len(bt)}")
    if keep.sum() < min_branches:
        raise RuntimeError("too few branches with a trait value")
    t = bt[keep]
    R = RER[:, keep]
    tr = stats.rankdata(t)
    tr = (tr - tr.mean()) / tr.std()
    rho = np.empty(R.shape[0])
    for i in range(R.shape[0]):
        rr = stats.rankdata(R[i])
        rho[i] = float((tr * ((rr - rr.mean()) / rr.std())).mean())
    n = keep.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    p = 2 * stats.t.sf(np.abs(tstat), df=n - 2)
    sym = [g.split(".", 2)[1] if g.count(".") >= 2 else g for g in genes]
    return pd.DataFrame(dict(gene=genes, symbol=sym, rho=rho, p=p,
                             q=_bh(p))).set_index("gene")


# --------------------------------------------------------------------------
# tie the scan back to the tissue question
# --------------------------------------------------------------------------
def organ_enrichment(res, min_genes=25):
    """Do an organ's genes track the trait more than the rest of the genome?

    This is the join between T12 and T1-T11: the scan is per gene, the question
    is per organ. Mann-Whitney on rho (not on q) so that a coherent shift of a
    whole gene set is detectable even when no single gene clears BH on its own -
    which is the expected situation for a diffuse effect (see T11).
    """
    from . import t1_expression, t5_confounders
    feat, _ = t1_expression.build()
    feat = feat[feat["symbol"].notna()]
    sym2organ = {}
    d = feat.copy()
    if "group_enriched_in" not in d.columns:
        d["group_enriched_in"] = np.nan
    for organ in t5_confounders.organ_list(d):
        members = d[t5_confounders.organ_membership(d, organ)]["symbol"]
        sym2organ[organ] = set(members.dropna())

    rows = []
    for organ, syms in sym2organ.items():
        inset = res["symbol"].isin(syms)
        a = res.loc[inset, "rho"].dropna()
        b = res.loc[~inset, "rho"].dropna()
        if len(a) < min_genes:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(organ=organ, n=len(a),
                         median_rho=float(a.median()),
                         median_rho_rest=float(b.median()),
                         rank_biserial=float(2 * u / (len(a) * len(b)) - 1),
                         p=float(p)))
    out = pd.DataFrame(rows).set_index("organ")
    out["q"] = _bh(out["p"].to_numpy())
    return out.sort_values("median_rho")


def rate_heterogeneity(RER, genes, res):
    """How variable is each gene's rate ACROSS lineages?

    The testis question ("why does it change so much?") does not need a trait at
    all: a gene set whose rate swings from lineage to lineage is under a
    different regime from one that is uniformly fast. This is the standard
    deviation of a gene's RER across branches - lineage-to-lineage rate
    heterogeneity, independent of any phenotype.
    """
    from . import t1_expression, t5_confounders
    het = RER.std(axis=1)
    sym = res["symbol"].to_numpy()
    feat, _ = t1_expression.build()
    d = feat[feat["symbol"].notna()].copy()
    if "group_enriched_in" not in d.columns:
        d["group_enriched_in"] = np.nan
    rows = []
    allhet = pd.Series(het, index=sym)
    for organ in t5_confounders.organ_list(d):
        syms = set(d[t5_confounders.organ_membership(d, organ)]["symbol"].dropna())
        a = allhet[allhet.index.isin(syms)]
        b = allhet[~allhet.index.isin(syms)]
        if len(a) < 25:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(organ=organ, n=len(a), median_het=float(a.median()),
                         median_het_rest=float(b.median()),
                         rank_biserial=float(2 * u / (len(a) * len(b)) - 1),
                         p=float(p)))
    out = pd.DataFrame(rows).set_index("organ")
    out["q"] = _bh(out["p"].to_numpy())
    return out.sort_values("median_het", ascending=False), allhet


def make_figure(res, orgres, hetres, bt, RER):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    axes[0].hist(res["rho"].dropna(), bins=70, color=PALETTE[0], alpha=0.85)
    axes[0].axvline(0, color="k", lw=0.9)
    sig = (res["q"] < 0.05).sum()
    axes[0].set_xlabel(r"Spearman $\rho$ (gene RER vs branch relative brain size)")
    axes[0].set_ylabel("genes")
    axes[0].set_title(f"{len(res)} genes scanned, {sig} at BH q<0.05")

    o = orgres.copy()
    colours = [PALETTE[1] if q < 0.05 else "#BBBBBB" for q in o["q"]]
    axes[1].barh(o.index, o["median_rho"], color=colours)
    axes[1].axvline(float(res["rho"].median()), color="k", ls="--", lw=0.9,
                    label="genome median")
    axes[1].set_xlabel(r"organ median $\rho$")
    axes[1].set_title("which organs' genes track\nencephalisation (orange: q<0.05)")
    axes[1].grid(axis="y", alpha=0)
    axes[1].legend(fontsize=7)

    h = hetres.head(14)
    colours = [PALETTE[2] if q < 0.05 else "#BBBBBB" for q in h["q"]]
    axes[2].barh(h.index[::-1], h["median_het"][::-1], color=colours[::-1])
    axes[2].axvline(float(np.median(hetres["median_het_rest"])), color="k",
                    ls="--", lw=0.9, label="rest of genome")
    axes[2].set_xlabel("median lineage-to-lineage\nrate heterogeneity (SD of RER)")
    axes[2].set_title("rate variability across lineages")
    axes[2].grid(axis="y", alpha=0)
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t12_rer_zoonomia.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(processes=4, limit=None, force=False, top=25):
    print("=" * 72)
    print("T12 - RERconverge across the Zoonomia mammals (TOGA codon alignments)")
    print("=" * 72)
    B, genes, order, splits, ncod = build_matrix(processes=processes,
                                                 limit=limit, force=force)
    print(f"genes with a complete core-species alignment : {len(genes)}")
    print(f"core mammal species                          : {len(order)}")
    print(f"identifiable branches                        : {B.shape[1]}")
    print(f"median aligned codons per gene               : {int(np.median(ncod))}")

    RER = compute_rer(B)
    trait, slope = relative_brain_size(order)
    print(f"\ntrait: log10 brain ~ log10 body slope = {slope:.3f} "
          f"(allometric exponent; residual = encephalisation)")
    print(f"core species with brain+body mass            : {len(trait)}/{len(order)}")

    tree = load_tree(set(order))
    bt = branch_trait(tree, order, splits, trait.to_dict())
    res = scan(RER, genes, bt)
    res = res.sort_values("rho")

    n_sig = int((res["q"] < 0.05).sum())
    keep = np.isfinite(bt)
    Rr = rer_ranks(RER, keep)
    null = null_calibration(tree, order, splits, Rr, keep, n_sim=200)
    print(f"\ngenes at BH q<0.05                           : {n_sig} / {len(res)}")
    if len(null):
        exceed = float((null >= n_sig).mean())
        print(f"NULL CALIBRATION - {len(null)} traits simulated under Brownian motion")
        print(f"  on the same tree (same phylogenetic clustering, no real biology):")
        print(f"  hits at q<0.05 under the null: median {int(np.median(null))}, "
              f"90th pct {int(np.percentile(null, 90))}, max {int(null.max())}")
        print(f"  empirical p for observing >= {n_sig} hits: {exceed:.3f}")
        if exceed > 0.05:
            print("  -> the real trait does NOT beat a structureless trait with the")
            print("     same phylogenetic autocorrelation. The gene list below is")
            print("     NOT distinguishable from noise and must not be read as biology.")
        np.save(os.path.join(DERIVED, "t12_null_hits.npy"), null)
    print(f"  accelerated in encephalised lineages (rho>0): "
          f"{int(((res['q'] < 0.05) & (res['rho'] > 0)).sum())}")
    print(f"  decelerated in encephalised lineages (rho<0): "
          f"{int(((res['q'] < 0.05) & (res['rho'] < 0)).sum())}")

    print(f"\ntop {top} genes ACCELERATED with encephalisation:")
    for _, r in res.sort_values("rho", ascending=False).head(top).iterrows():
        print(f"   {r['symbol']:<14} rho={r['rho']:+.3f}  q={r['q']:.2e}")
    print(f"\ntop {top} genes DECELERATED (more constrained) with encephalisation:")
    for _, r in res.head(top).iterrows():
        print(f"   {r['symbol']:<14} rho={r['rho']:+.3f}  q={r['q']:.2e}")

    orgres = organ_enrichment(res)
    print("\nper-organ shift in rho (Mann-Whitney vs rest of genome):")
    print(f"   {'organ':<18}{'n':>6}{'med rho':>10}{'rest':>9}{'rb':>7}{'q':>10}")
    for organ, r in orgres.iterrows():
        star = "*" if r["q"] < 0.05 else " "
        print(f"   {organ:<18}{int(r['n']):>6}{r['median_rho']:>+10.4f}"
              f"{r['median_rho_rest']:>+9.4f}{r['rank_biserial']:>+7.2f}"
              f"{r['q']:>10.2e}{star}")

    hetres, allhet = rate_heterogeneity(RER, genes, res)
    print("\nlineage-to-lineage rate heterogeneity (SD of RER across branches):")
    print(f"   {'organ':<18}{'n':>6}{'median':>9}{'rest':>9}{'rb':>7}{'q':>10}")
    for organ, r in hetres.iterrows():
        star = "*" if r["q"] < 0.05 else " "
        print(f"   {organ:<18}{int(r['n']):>6}{r['median_het']:>9.4f}"
              f"{r['median_het_rest']:>9.4f}{r['rank_biserial']:>+7.2f}"
              f"{r['q']:>10.2e}{star}")

    to_tsv_gz(res, os.path.join(DERIVED, "t12_rer_brain_scan.tsv.gz"))
    orgres.to_csv(os.path.join(DERIVED, "t12_organ_rho.tsv"), sep="\t")
    hetres.to_csv(os.path.join(DERIVED, "t12_organ_heterogeneity.tsv"), sep="\t")
    p = make_figure(res, orgres, hetres, bt, RER)
    print(f"\nfigure written to: {p}")
    return res, orgres, hetres


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    report(processes=a.processes, limit=a.limit, force=a.force)
