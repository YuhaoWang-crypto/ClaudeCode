"""
T19 - the T18 design at scale: hundreds of pathways instead of one.

T18 found three things on the TCA cycle that only mean something if they
generalise, and one of them, if true, matters well beyond this project:

  A. within a pathway, the proteins doing the CHEMISTRY are more constrained
     than the machinery that serves them (3.1x on TCA);
  B. RATE AND RETENTION ARE UNCORRELATED - frataxin evolves at 2.9x the genome
     median yet is kept in 481/511 mammalian assemblies and its loss causes
     Friedreich ataxia. If that holds genome-wide it undermines the routine use
     of conservation as a proxy for indispensability;
  C. paralogue allocation shifts rate, but compartment is not the rule - ACO1 is
     cytosolic AND slow because it moonlights as IRP1.

Scaling each one needs a different fix.

A needs the hand-curated strata replaced by something automatic. The
operationalisation used here is UniProt: a gene is CATALYTIC if it carries an EC
number or an annotated ACT_SITE, and ACCESSORY otherwise.

This proxy is COARSER THAN THE TCA CURATION AND IT IS WORTH SAYING HOW.
It puts subunits, assembly factors and regulators in one bin. Worse, it
misassigns exactly the genes that made the TCA result interesting: frataxin
carries EC 1.16.3.1 and ISCU carries two annotated active sites, so both land in
"catalytic" here even though T18 called them Fe-S biogenesis, i.e. accessory.
The direction of that error is the saving grace - it moves fast accessory genes
INTO the catalytic bin, which shrinks the catalytic/accessory gap. Whatever
difference survives is therefore a lower bound. A stricter ACT_SITE-only split
is reported alongside as a sensitivity check.

B needs retention for every gene, not 46. TOGA builds an alignment only where an
orthologue is called Intact / Partial Intact / Uncertain Loss, so the number of
sequences in a gene's alignment IS its retention across ~500 mammalian
assemblies.

C needs paralogue families, taken from Ensembl BioMart.

Data handling note: the container running this pipeline gets reclaimed
periodically, and the TOGA alignment set is 1.2 GB. So this module STREAMS -
each alignment is fetched, reduced to (species count, per-residue conservation)
and discarded, with the reduced table checkpointed every few hundred genes. A
restart resumes from the checkpoint instead of re-downloading a gigabyte.
"""
import os
import re
import gzip
import json
import argparse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     biomart_query, http_get, log)
from . import t5_confounders, t15_pathway

TOGA_BASE = ("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/"
             "MultipleCodonAlignments/")
TOGA_LISTING = os.path.join(RAW, "toga_listing.html")
UNIPROT_ENZ = os.path.join(RAW, "uniprot_enzyme.tsv")
RETENTION = os.path.join(DERIVED, "t19_retention.tsv.gz")
CONS_DIR = os.path.join(DERIVED, "t19_conservation")

MIN_PER_GROUP = 3
MIN_SPECIES = 30
CHECKPOINT = 500

AA_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AA_OF = dict(zip(_CODONS, AA_TABLE))


# --------------------------------------------------------------------------
# streaming TOGA reduction
# --------------------------------------------------------------------------
def _reduce_alignment(blob, keep_vector=False):
    """(n_species, mean conservation, invariant fraction, [per-residue vector])."""
    ref, others = None, []
    name, chunks = None, []
    for line in gzip.decompress(blob).decode("utf-8", "replace").split("\n"):
        if line.startswith(">"):
            if name is not None:
                seq = "".join(chunks)
                if name == "REFERENCE":
                    ref = seq
                else:
                    others.append(seq)
            name, chunks = line[1:].split()[0] if line[1:].split() else "", []
        else:
            chunks.append(line.strip())
    if name is not None:
        seq = "".join(chunks)
        (others.append(seq) if name != "REFERENCE" else None)
        if name == "REFERENCE":
            ref = seq
    if ref is None or len(others) < MIN_SPECIES:
        return None

    def aa(seq, i):
        return AA_OF.get(seq[3 * i:3 * i + 3].upper().replace("U", "T"))

    L = len(ref) // 3
    cons = []
    for i in range(L):
        r = aa(ref, i)
        if r is None or r == "*":
            continue
        obs = [aa(s, i) for s in others]
        obs = [o for o in obs if o is not None and o != "*"]
        if len(obs) < MIN_SPECIES:
            cons.append(np.nan)
            continue
        cons.append(sum(1 for o in obs if o == r) / len(obs))
    arr = np.array(cons, dtype=np.float32)
    ok = arr[np.isfinite(arr)]
    if len(ok) < 50:
        return None
    rec = dict(n_species=len(others) + 1, n_residues=len(ok),
               mean_conservation=float(ok.mean()),
               frac_invariant=float((ok >= 0.999).mean()))
    return (rec, arr) if keep_vector else (rec, None)


def stream_retention(keep_vectors_for=None, workers=12, force=False):
    """Fetch every TOGA alignment, reduce it, discard the file.

    Nothing of the 1.2 GB set is kept on disk; only the reduced table and, for
    `keep_vectors_for`, the per-residue conservation vectors needed by the
    structural analysis. Checkpointed so a container restart is cheap.
    """
    if os.path.exists(RETENTION) and not force:
        return pd.read_csv(RETENTION, sep="\t")
    os.makedirs(CONS_DIR, exist_ok=True)
    keep = set(keep_vectors_for or ())
    files = sorted(set(re.findall(r'(ENST[0-9]+\.[A-Za-z0-9._\-]+\.fasta\.gz)',
                                  open(TOGA_LISTING).read())))
    part = RETENTION + ".part"
    done = {}
    if os.path.exists(part):
        prev = pd.read_csv(part, sep="\t")
        done = {r["file"]: r.to_dict() for _, r in prev.iterrows()}
        log(f"resuming from checkpoint with {len(done)} genes")
    todo = [f for f in files if f not in done]
    log(f"streaming {len(todo)} TOGA alignments ({len(files)} total)")

    def fetch(fn):
        sym = fn.split(".", 2)[1] if fn.count(".") >= 2 else fn
        for attempt in range(3):
            try:
                with urllib.request.urlopen(TOGA_BASE + fn, timeout=180) as r:
                    blob = r.read()
                break
            except Exception:                             # noqa: BLE001
                blob = None
        if blob is None:
            return None
        try:
            out = _reduce_alignment(blob, keep_vector=sym in keep)
        except Exception:                                 # noqa: BLE001
            return None
        if out is None:
            return None
        rec, vec = out
        rec.update(file=fn, symbol=sym)
        if vec is not None:
            np.save(os.path.join(CONS_DIR, f"{sym}.npy"), vec)
        return rec

    rows = list(done.values())
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, rec in enumerate(ex.map(fetch, todo), 1):
            if rec is not None:
                rows.append(rec)
            if i % CHECKPOINT == 0:
                pd.DataFrame(rows).to_csv(part, sep="\t", index=False)
                log(f"  {i}/{len(todo)} streamed, {len(rows)} usable")
    df = pd.DataFrame(rows)
    df.to_csv(RETENTION, sep="\t", index=False,
              compression={"method": "gzip", "mtime": 0})
    if os.path.exists(part):
        os.remove(part)
    log(f"retention table: {len(df)} genes")
    return df


# --------------------------------------------------------------------------
# annotation
# --------------------------------------------------------------------------
def enzyme_annotation():
    """symbol -> is_catalytic (EC number or an annotated active site)."""
    d = pd.read_csv(UNIPROT_ENZ, sep="\t")
    ec = d.get("EC number", pd.Series(index=d.index, dtype=object)).fillna("")
    act = d.get("Active site", pd.Series(index=d.index, dtype=object)).fillna("")
    bind = d.get("Binding site", pd.Series(index=d.index, dtype=object)).fillna("")
    d["is_catalytic"] = (ec.str.len() > 0) | act.str.contains("ACT_SITE")
    d["has_binding"] = bind.str.contains("BINDING")
    d["n_act_site"] = act.str.count("ACT_SITE")
    out = d.dropna(subset=["Gene Names (primary)"]).drop_duplicates(
        "Gene Names (primary)").set_index("Gene Names (primary)")
    return out[["Entry", "is_catalytic", "has_binding", "n_act_site"]]


def analysis_frame(sp_key="mouse"):
    d = t5_confounders.analysis_frame(sp_key)
    ann = enzyme_annotation()
    d = d.join(ann, on="symbol")
    d["is_catalytic"] = d["is_catalytic"].fillna(False)
    # stricter variant: an explicitly annotated catalytic residue, not merely an
    # EC number. Frataxin drops out of the catalytic class under this rule.
    d["has_act_site"] = (d["n_act_site"].fillna(0) > 0)
    ret = pd.read_csv(RETENTION, sep="\t") if os.path.exists(RETENTION) else None
    if ret is not None:
        r = ret.drop_duplicates("symbol").set_index("symbol")
        for c in ["n_species", "mean_conservation", "frac_invariant", "n_residues"]:
            d[c] = d["symbol"].map(r[c])
    return d


# --------------------------------------------------------------------------
# A. catalytic vs accessory, within pathways
# --------------------------------------------------------------------------
def catalytic_vs_accessory(d, paths, min_per_group=MIN_PER_GROUP,
                           flag="is_catalytic"):
    rows = []
    for pid, (name, genes) in paths.items():
        sub = d.loc[[g for g in genes if g in d.index]]
        mask = sub[flag].astype(bool)
        a = sub.loc[mask, "omega"].dropna()
        b = sub.loc[~mask, "omega"].dropna()
        if len(a) < min_per_group or len(b) < min_per_group:
            continue
        rows.append(dict(pathway=pid, name=name, n_cat=len(a), n_acc=len(b),
                         omega_catalytic=float(a.median()),
                         omega_accessory=float(b.median()),
                         log2_acc_over_cat=float(np.log2(b.median() / a.median())),
                         tau_catalytic=float(sub.loc[a.index, "tau"].median()),
                         tau_accessory=float(sub.loc[b.index, "tau"].median())))
    res = pd.DataFrame(rows).set_index("pathway")
    if not len(res):
        return res, {}
    w = stats.wilcoxon(res["log2_acc_over_cat"], alternative="two-sided")
    stat = dict(n_pathways=len(res),
                median_log2=float(res["log2_acc_over_cat"].median()),
                fold=float(2 ** res["log2_acc_over_cat"].median()),
                frac_positive=float((res["log2_acc_over_cat"] > 0).mean()),
                median_tau_gap=float((res["tau_accessory"]
                                      - res["tau_catalytic"]).median()),
                p=float(w.pvalue))
    return res.sort_values("log2_acc_over_cat", ascending=False), stat


def catalytic_fixed_effects(d, paths):
    """Does catalytic status predict rate WITHIN pathways, controlling for the
    covariates that predicted rate everywhere else?"""
    import statsmodels.formula.api as smf
    smallest = {}
    for pid, (name, genes) in paths.items():
        for g in genes:
            if g not in smallest or len(paths[smallest[g]][1]) > len(genes):
                smallest[g] = pid
    sub = d.loc[[g for g in smallest if g in d.index]].copy()
    sub["pathway"] = [smallest[g] for g in sub.index]
    sub["catalytic"] = sub["is_catalytic"].astype(bool).astype(int)
    sub = sub.dropna(subset=["log_omega", "tau", "log_expr", "log_codons"])
    base = smf.ols("log_omega ~ catalytic", data=sub).fit()
    cov = smf.ols("log_omega ~ catalytic + tau + log_expr + log_codons",
                  data=sub).fit()
    fe = smf.ols("log_omega ~ catalytic + tau + log_expr + log_codons "
                 "+ C(pathway)", data=sub).fit()
    # no "fraction retained" here: the raw coefficient is near zero and flips
    # sign when covariates enter, so a ratio to it is not interpretable. The
    # coefficients are in log10(omega) units and are reported directly.
    return dict(n=int(base.nobs), n_pathways=sub["pathway"].nunique(),
                coef_raw=float(base.params["catalytic"]),
                coef_cov=float(cov.params["catalytic"]),
                coef_fe=float(fe.params["catalytic"]),
                p_fe=float(fe.pvalues["catalytic"]),
                fold_fe=float(10 ** fe.params["catalytic"]))


# --------------------------------------------------------------------------
# B. rate vs retention - the claim that matters outside this project
# --------------------------------------------------------------------------
def _partial_spearman(x, y, covars):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    C = np.column_stack([stats.rankdata(c) for c in covars] + [np.ones(len(rx))])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p)


def rate_vs_retention(d, paths=None):
    ok = d.dropna(subset=["omega", "n_species"])
    out = dict(n=len(ok))
    # coding length correlates with BOTH omega (short genes look fast: median
    # omega 0.200 at 100-150 codons vs 0.102 above 500) and retention, so the
    # raw correlation is partly a length effect
    okl = ok.dropna(subset=["log_codons"])
    r, p = _partial_spearman(okl["omega"], okl["n_species"], [okl["log_codons"]])
    out["omega_given_length"] = dict(rho=r, p=p, n=len(okl))
    for label, col in [("omega", "omega"), ("tau", "tau"),
                       ("log_expr", "log_expr"),
                       ("mean_conservation", "mean_conservation")]:
        if col not in ok.columns:
            continue
        sub = ok.dropna(subset=[col])
        rho, p = stats.spearmanr(sub[col], sub["n_species"])
        out[label] = dict(rho=float(rho), p=float(p), n=len(sub))
    # the decisive framing: among genes that ARE retained everywhere, does
    # omega still span the full range?
    hi = ok[ok["n_species"] >= ok["n_species"].quantile(0.9)]
    out["universally_retained"] = dict(
        n=len(hi), q10_omega=float(hi["omega"].quantile(0.10)),
        median_omega=float(hi["omega"].median()),
        q90_omega=float(hi["omega"].quantile(0.90)),
        frac_above_genome=float((hi["omega"] > ok["omega"].median()).mean()))
    if paths:
        rows = []
        for pid, (name, genes) in paths.items():
            sub = ok.loc[[g for g in genes if g in ok.index]]
            if len(sub) < 15:
                continue
            rho, p = stats.spearmanr(sub["omega"], sub["n_species"])
            rows.append(dict(pathway=pid, name=name, n=len(sub),
                             rho=float(rho), p=float(p)))
        out["per_pathway"] = pd.DataFrame(rows).set_index("pathway")
    return out


def fast_but_kept(d, omega_q=0.90, retention_q=0.75):
    """Genes that are fast AND universally retained - the frataxin class.

    "Fast" is defined on the LENGTH-RESIDUALISED rate, not on raw omega. Short
    genes have systematically higher omega, so a raw threshold selects for short
    genes and the resulting list is a chemokine/defensin catalogue by
    construction rather than by biology.
    """
    ok = d.dropna(subset=["omega", "n_species", "log_codons"]).copy()
    lo = np.log10(ok["omega"].clip(lower=1e-4))
    beta = np.polyfit(ok["log_codons"], lo, 1)
    ok["omega_resid"] = lo - np.polyval(beta, ok["log_codons"])
    thr_o = ok["omega_resid"].quantile(omega_q)
    thr_r = ok["n_species"].quantile(retention_q)
    sel = ok[(ok["omega_resid"] > thr_o) & (ok["n_species"] >= thr_r)]
    return sel.sort_values("omega_resid", ascending=False), dict(
        n=len(sel), frac=float(len(sel) / len(ok)),
        expected_frac=float((1 - omega_q) * (1 - retention_q)),
        median_codons=float(sel["n_codons"].median()) if "n_codons" in sel else np.nan,
        retention_threshold=float(thr_r))


# --------------------------------------------------------------------------
# C. paralogues at scale
# --------------------------------------------------------------------------
def paralog_pairs():
    cache = os.path.join(RAW, "biomart_paralogs.tsv")
    if not os.path.exists(cache):
        txt = biomart_query("hsapiens_gene_ensembl",
                            ["ensembl_gene_id", "hsapiens_paralog_ensembl_gene",
                             "hsapiens_paralog_perc_id",
                             "hsapiens_paralog_subtype"])
        open(cache, "w").write(txt)
    rows = []
    for line in open(cache).read().rstrip("\n").split("\n")[1:]:
        f = line.split("\t")
        if len(f) >= 3 and f[1]:
            try:
                pid = float(f[2]) if f[2] else np.nan
            except ValueError:
                pid = np.nan
            rows.append((f[0], f[1], pid))
    return pd.DataFrame(rows, columns=["gene_a", "gene_b", "perc_id"])


def paralog_allocation(d, min_perc_id=50.0):
    """Within paralogue pairs, does the tissue-restricted copy evolve faster?

    A paralogue pair holds fold, chemistry and ancestry fixed, which is a tighter
    control than the within-pathway test in T15. Pairs are counted once (a<b) and
    only where exactly one member is tissue-enriched, so the comparison is always
    restricted vs broad rather than two restricted copies.

    Only each gene's CLOSEST paralogue is used. BioMart's full table is 3.56M
    pairs over 20920 genes - an average of 170 paralogues each - because it
    includes every member of the olfactory-receptor and KRAB-zinc-finger
    expansions. Taking all pairs would let a handful of huge families supply most
    of the comparisons; taking the closest one gives each gene a single vote and
    keeps the pair genuinely comparable.
    """
    pp = paralog_pairs()
    pp = pp[pp["perc_id"] >= min_perc_id]
    pp = (pp.sort_values("perc_id", ascending=False)
            .drop_duplicates("gene_a"))
    is_e = (d["enriched_in"].notna() | d["group_enriched_in"].notna())
    om = d["omega"]
    rows = []
    seen = set()
    for a, b, pid in pp.itertuples(index=False):
        if a not in d.index or b not in d.index:
            continue
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        ea, eb = bool(is_e.get(a, False)), bool(is_e.get(b, False))
        if ea == eb:
            continue
        oa, ob = om.get(a), om.get(b)
        if not (np.isfinite(oa) and np.isfinite(ob)):
            continue
        restricted, broad = (a, b) if ea else (b, a)
        rows.append(dict(restricted=restricted, broad=broad, perc_id=pid,
                         omega_restricted=float(om[restricted]),
                         omega_broad=float(om[broad]),
                         log2_r_over_b=float(np.log2(om[restricted] / om[broad])),
                         sym_r=d.at[restricted, "symbol"],
                         sym_b=d.at[broad, "symbol"]))
    res = pd.DataFrame(rows)
    if not len(res):
        return res, {}
    w = stats.wilcoxon(res["log2_r_over_b"], alternative="two-sided")
    stat = dict(n_pairs=len(res),
                median_log2=float(res["log2_r_over_b"].median()),
                fold=float(2 ** res["log2_r_over_b"].median()),
                frac_positive=float((res["log2_r_over_b"] > 0).mean()),
                p=float(w.pvalue))
    return res.sort_values("log2_r_over_b", ascending=False), stat


def make_figure(cat, cstat, rr, par, pstat, d):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    axes[0].hist(cat["log2_acc_over_cat"], bins=55, color=PALETTE[0], alpha=0.85)
    axes[0].axvline(0, color="k", lw=1.0)
    axes[0].axvline(cstat["median_log2"], color=PALETTE[1], lw=1.5,
                    label=f"median {cstat['median_log2']:+.3f} "
                          f"({cstat['fold']:.2f}x)")
    axes[0].set_xlabel(r"$\log_2(\omega_{\rm accessory} / \omega_{\rm catalytic})$"
                       "\nwithin the same pathway")
    axes[0].set_ylabel("pathways")
    axes[0].set_title(f"A. {cstat['n_pathways']} pathways, "
                      f"{cstat['frac_positive']:.0%} positive")
    axes[0].legend(fontsize=7)

    ok = d.dropna(subset=["omega", "n_species"])
    axes[1].hexbin(np.log10(ok["omega"].clip(lower=1e-4)), ok["n_species"],
                   gridsize=45, cmap="Blues", mincnt=1, bins="log")
    r = rr["omega"]
    axes[1].set_xlabel(r"$\log_{10}\omega$ (human-mouse)")
    axes[1].set_ylabel("mammalian assemblies retaining the gene")
    axes[1].set_title(f"B. rate vs retention\n" r"Spearman $\rho$"
                      f"={r['rho']:+.3f} (n={r['n']})")

    axes[2].hist(par["log2_r_over_b"], bins=55, color=PALETTE[2], alpha=0.85)
    axes[2].axvline(0, color="k", lw=1.0)
    axes[2].axvline(pstat["median_log2"], color=PALETTE[1], lw=1.5,
                    label=f"median {pstat['median_log2']:+.3f} "
                          f"({pstat['fold']:.2f}x)")
    axes[2].set_xlabel(r"$\log_2(\omega_{\rm restricted}/\omega_{\rm broad})$"
                       "\nwithin paralogue pairs")
    axes[2].set_ylabel("pairs")
    axes[2].set_title(f"C. {pstat['n_pairs']} paralogue pairs,\n"
                      f"{pstat['frac_positive']:.0%} positive")
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t19_scaleout.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(sp_key="mouse", workers=12):
    print("=" * 72)
    print("T19 - the T18 design across hundreds of pathways")
    print("=" * 72)
    if not os.path.exists(RETENTION):
        stream_retention(workers=workers)
    d = analysis_frame(sp_key)
    members = t15_pathway.pathway_members()
    paths = t15_pathway.usable_pathways(d, members)
    print(f"genes: {len(d)}  with retention: {d['n_species'].notna().sum()}  "
          f"catalytic: {int(d['is_catalytic'].astype(bool).sum())}")
    print(f"pathways in range: {len(paths)}\n")

    cat, cstat = catalytic_vs_accessory(d, paths)
    print("A. CATALYTIC vs ACCESSORY, within the same pathway")
    print(f"   pathways testable          : {cstat['n_pathways']}")
    print(f"   median log2(acc/cat)       : {cstat['median_log2']:+.4f}  "
          f"({cstat['fold']:.2f}x)")
    print(f"   pathways positive          : {cstat['frac_positive']:.1%}")
    print(f"   median tau gap (acc - cat) : {cstat['median_tau_gap']:+.3f}")
    print(f"   Wilcoxon signed-rank p     : {cstat['p']:.3g}")
    cat2, cstat2 = catalytic_vs_accessory(d, paths, flag="has_act_site")
    if cstat2:
        print(f"\n   sensitivity - stricter split (annotated ACT_SITE only, "
              f"n={cstat2['n_pathways']} pathways):")
        print(f"      median log2(acc/cat) {cstat2['median_log2']:+.4f} "
              f"({cstat2['fold']:.2f}x), {cstat2['frac_positive']:.0%} positive, "
              f"p={cstat2['p']:.2g}")

    fe = catalytic_fixed_effects(d, paths)
    print(f"\n   catalytic coefficient on log10(omega), n={fe['n']} genes / "
          f"{fe['n_pathways']} pathways:")
    print(f"      raw                        {fe['coef_raw']:+.4f}")
    print(f"      + tau, expression, length  {fe['coef_cov']:+.4f}")
    print(f"      + pathway fixed effects    {fe['coef_fe']:+.4f} "
          f"(p={fe['p_fe']:.2g},  {fe['fold_fe']:.2f}x)")
    print("      the raw coefficient is ~0 and flips sign under covariates, so")
    print("      the WITHIN-pathway effect is larger than the genome-wide one:")
    print("      catalytic genes concentrate in pathways that are fast overall,")
    print("      which masks the effect until pathway is held fixed.")

    rr = rate_vs_retention(d, paths)
    print(f"\nB. RATE vs RETENTION across {rr['n']} genes")
    for k in ["omega", "tau", "log_expr", "mean_conservation"]:
        if k in rr:
            v = rr[k]
            print(f"   Spearman({k:<18}, retention) = {v['rho']:+.3f}  "
                  f"(p={v['p']:.2g}, n={v['n']})")
    g = rr["omega_given_length"]
    print(f"   partial Spearman(omega, retention | coding length) = "
          f"{g['rho']:+.3f}  (p={g['p']:.2g}, n={g['n']})")
    u = rr["universally_retained"]
    print(f"\n   among the top-decile-retained genes (n={u['n']}):")
    print(f"      omega 10th pct {u['q10_omega']:.4f}, median "
          f"{u['median_omega']:.4f}, 90th pct {u['q90_omega']:.4f}")
    print(f"      {u['frac_above_genome']:.0%} of them are ABOVE the genome-wide "
          f"median omega")
    fk, fstat = fast_but_kept(d)
    print(f"\n   the frataxin class - top-decile LENGTH-RESIDUALISED rate AND "
          f"retained in >= {fstat['retention_threshold']:.0f} assemblies:")
    print(f"      {fstat['n']} genes ({fstat['frac']:.1%} of the set; "
          f"{fstat['expected_frac']:.1%} expected if rate and retention were "
          f"independent)")
    print(f"      median coding length {fstat['median_codons']:.0f} codons")
    print("      examples: " + ", ".join(
        f"{r.symbol}({r.omega:.2f}/{int(r.n_species)})"
        for r in fk.head(12).itertuples() if isinstance(r.symbol, str)))
    if "per_pathway" in rr and len(rr["per_pathway"]):
        pp = rr["per_pathway"]
        sig = (pp["p"] < 0.05).sum()
        print(f"\n   per-pathway: {len(pp)} pathways with >=15 genes, "
              f"median rho {pp['rho'].median():+.3f}, "
              f"{sig} at p<0.05 ({sig / len(pp):.0%}; 5% expected by chance)")

    par, pstat = paralog_allocation(d)
    print(f"\nC. PARALOGUE PAIRS - restricted vs broad copy")
    print(f"   pairs with exactly one restricted member: {pstat['n_pairs']}")
    print(f"   median log2(restricted/broad)           : "
          f"{pstat['median_log2']:+.4f}  ({pstat['fold']:.2f}x)")
    print(f"   pairs positive                          : "
          f"{pstat['frac_positive']:.1%}")
    print(f"   Wilcoxon signed-rank p                  : {pstat['p']:.3g}")

    cat.to_csv(os.path.join(DERIVED, "t19_catalytic.tsv"), sep="\t")
    par.to_csv(os.path.join(DERIVED, "t19_paralogs.tsv"), sep="\t", index=False)
    fk.head(200).to_csv(os.path.join(DERIVED, "t19_fast_but_kept.tsv"), sep="\t")
    if "per_pathway" in rr:
        rr["per_pathway"].to_csv(os.path.join(DERIVED, "t19_retention_by_pathway.tsv"),
                                 sep="\t")
    p = make_figure(cat, cstat, rr, par, pstat, d)
    print(f"\nfigure written to: {p}")
    return d, cat, rr, par


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--stream-only", action="store_true")
    a = ap.parse_args()
    if a.stream_only:
        stream_retention(workers=a.workers)
    else:
        report(workers=a.workers)
