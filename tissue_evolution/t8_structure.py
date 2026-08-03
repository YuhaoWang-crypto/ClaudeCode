"""
T8 - structural divergence (Foldseek/TM-align) next to sequence divergence.

Sequence and structure do not diverge at the same speed. Fold is far more
conserved than the residues that build it, which is why Foldseek finds homologues
past the point where alignment-based methods stop working. That gives two things
this pipeline needs:

  * a divergence measure that has NOT saturated at the depths where dS has (T6),
  * a way to separate "this protein's sequence is running fast" from "this
    protein's shape is changing" - a tissue whose proteins accumulate
    substitutions while holding their fold is under a very different regime from
    one whose proteins are being structurally remodelled.

Method: AlphaFold DB models for both members of a 1:1 human/mouse orthologue
pair, then `foldseek easy-search --alignment-type 1` (TM-align) for the TM-score,
and Foldseek's 3Di structural alphabet for a sequence-like structural identity.

Two limits worth stating plainly. First, AlphaFold models are single static
predicted conformers, so a TM-score between two of them measures agreement
between two PREDICTIONS, and predicted models of close homologues are biased
toward looking alike - this compresses the structural-divergence axis and makes
it conservative. Second, TM-score saturates near 1.0 for orthologues this close;
human/mouse is deliberately the wrong depth to see structural change, and is used
here to establish the baseline against which the deep comparisons are read.
"""
import os
import json
import subprocess
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     biomart_query, http_get, download, log)
from . import t4_tissue_rates

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/"
STRUCT_DIR = os.path.join(RAW, "structures")
FOCUS_ORGANS = ["Testis", "Brain", "Cerebellum", "Liver", "Blood",
                "Skeletal_muscle", "Skin", "Pituitary", "Kidney", "Spleen"]


def uniprot_map(dataset, cache_name):
    """ensembl gene -> reviewed SwissProt accession (AFDB covers SwissProt)."""
    path = os.path.join(RAW, cache_name)
    if not os.path.exists(path):
        txt = biomart_query(dataset, ["ensembl_gene_id", "uniprotswissprot"])
        open(path, "w").write(txt)
    txt = open(path).read()
    out = {}
    for line in txt.rstrip("\n").split("\n")[1:]:
        f = line.split("\t")
        if len(f) >= 2 and f[1] and f[0] not in out:
            out[f[0]] = f[1]
    return out


def fetch_structure(acc, subdir):
    """Download the AFDB model for one accession. Returns path or None."""
    d = os.path.join(STRUCT_DIR, subdir)
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, f"{acc}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    try:
        meta = json.loads(http_get(AFDB_API + acc, timeout=120, retries=2))
        if not meta:
            return None
        url = meta[0].get("pdbUrl")
        if not url:
            return None
        download(url, dest, timeout=300, retries=3)
        return dest
    except Exception:                                # noqa: BLE001
        return None


def fetch_all(accs, subdir, workers=8):
    with ThreadPoolExecutor(max_workers=workers) as ex:
        paths = list(ex.map(lambda a: fetch_structure(a, subdir), accs))
    ok = {a: p for a, p in zip(accs, paths) if p}
    log(f"  {subdir}: {len(ok)}/{len(accs)} structures retrieved")
    return ok


def build_sample(n_per_organ=30, seed=0, sp_key="mouse"):
    df, _, _ = t4_tissue_rates.build(sp_key)
    hm = pd.read_csv(os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz"),
                     sep="\t").set_index("human_gene")
    df = df.join(hm[["ortholog_gene"]], how="inner")

    h_map = uniprot_map("hsapiens_gene_ensembl", "uniprot_human.tsv")
    m_map = uniprot_map("mmusculus_gene_ensembl", "uniprot_mouse.tsv")

    df["acc_human"] = pd.Series(h_map).reindex(df.index)
    df["acc_ortholog"] = pd.Series(m_map).reindex(df["ortholog_gene"]).to_numpy()
    df = df.dropna(subset=["acc_human", "acc_ortholog", "dN"])

    rng = np.random.default_rng(seed)
    picks = []
    for organ in FOCUS_ORGANS:
        sub = df[df["enriched_in"] == organ]
        if len(sub) == 0:
            continue
        take = sub.index.to_numpy()
        if len(take) > n_per_organ:
            take = rng.choice(take, n_per_organ, replace=False)
        picks.extend(take)
    bg = df[df["enriched_in"].isna() & (df["tau"] < 0.5)]
    picks.extend(rng.choice(bg.index.to_numpy(), min(n_per_organ * 2, len(bg)),
                            replace=False))
    return df.loc[pd.Index(picks).drop_duplicates()]


def run_foldseek(sample, workdir=None, threads=2):
    workdir = workdir or os.path.join(RAW, "foldseek_run")
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir, exist_ok=True)

    h_ok = fetch_all(sorted(set(sample["acc_human"])), "human")
    m_ok = fetch_all(sorted(set(sample["acc_ortholog"])), "mouse")
    usable = sample[sample["acc_human"].isin(h_ok)
                    & sample["acc_ortholog"].isin(m_ok)]
    log(f"  pairs with both structures: {len(usable)}")
    if len(usable) < 20:
        return None, usable

    qdir, tdir = os.path.join(workdir, "q"), os.path.join(workdir, "t")
    os.makedirs(qdir), os.makedirs(tdir)
    for a in set(usable["acc_human"]):
        shutil.copy(h_ok[a], qdir)
    for a in set(usable["acc_ortholog"]):
        shutil.copy(m_ok[a], tdir)

    aln = os.path.join(workdir, "aln.tsv")
    # `-e inf --max-seqs 5000` makes this an exhaustive all-vs-all TM-align
    # (347x347 here) when only the 347 designated orthologue pairs are ever
    # read back - roughly a 200x waste, and slow enough to starve anything else
    # on the box. The 3Di prefilter ranks a true orthologue in the top handful
    # of hits, so a small --max-seqs costs no recall on the pairs we want.
    cmd = ["foldseek", "easy-search", qdir, tdir, aln,
           os.path.join(workdir, "tmp"),
           "--alignment-type", "1",          # TM-align
           "-e", "10", "--max-seqs", "25", "--threads", str(threads),
           "--format-output", "query,target,fident,alnlen,evalue,alntmscore,lddt"]
    log("  running foldseek easy-search (TM-align mode)")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        log(f"  foldseek failed: {res.stderr[-600:]}")
        return None, usable
    cols = ["query", "target", "fident", "alnlen", "evalue", "alntmscore", "lddt"]
    hits = pd.read_csv(aln, sep="\t", names=cols)
    hits["query"] = hits["query"].str.replace(".pdb", "", regex=False)
    hits["target"] = hits["target"].str.replace(".pdb", "", regex=False)
    return hits, usable


def mean_plddt(pdb_path):
    """AlphaFold writes per-residue pLDDT into the B-factor column."""
    vals = []
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    vals.append(float(line[60:66]))
    except Exception:                                    # noqa: BLE001
        return np.nan
    return float(np.mean(vals)) if vals else np.nan


def add_plddt(res, sample):
    """Attach mean pLDDT for both models.

    This is not a decoration. A TM-score between two AlphaFold models is
    depressed wherever either model is UNCONFIDENT: low-pLDDT regions are
    intrinsically disordered or simply unpredicted, and their coordinates are
    essentially arbitrary, so two models of the same conserved protein can
    disagree there for no evolutionary reason. Proteins of the CNS are strongly
    enriched in intrinsically disordered regions, which is exactly the direction
    that would manufacture a fake "brain proteins are structurally divergent"
    result. Any organ-level TM-score claim has to be read against this column.
    """
    acc_h = sample["acc_human"].reindex(res.index)
    acc_m = sample["acc_ortholog"].reindex(res.index)
    res = res.copy()
    res["plddt_human"] = [mean_plddt(os.path.join(STRUCT_DIR, "human", f"{a}.pdb"))
                          if isinstance(a, str) else np.nan for a in acc_h]
    res["plddt_mouse"] = [mean_plddt(os.path.join(STRUCT_DIR, "mouse", f"{a}.pdb"))
                          if isinstance(a, str) else np.nan for a in acc_m]
    res["plddt_min"] = res[["plddt_human", "plddt_mouse"]].min(axis=1)
    return res


def join_pairs(hits, usable):
    """Keep only the designated orthologue pairs, best hit per pair."""
    key = hits.set_index(["query", "target"])
    rows = []
    for gene, r in usable.iterrows():
        k = (r["acc_human"], r["acc_ortholog"])
        if k not in key.index:
            continue
        h = key.loc[k]
        if isinstance(h, pd.DataFrame):
            h = h.sort_values("alntmscore", ascending=False).iloc[0]
        rows.append(dict(human_gene=gene, organ=r["enriched_in"] or "_background",
                         dN=r["dN"], dS=r["dS"], omega=r["omega"],
                         tau=r["tau"], identity_cds=r["identity"],
                         fident=float(h["fident"]),
                         tmscore=float(h["alntmscore"]),
                         lddt=float(h["lddt"])))
    return pd.DataFrame(rows).set_index("human_gene")


def make_figure(res):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))

    ok = res.dropna(subset=["dN", "tmscore"])
    axes[0].scatter(ok["dN"], ok["tmscore"], s=14, alpha=0.6, color=PALETTE[0])
    rho, p = stats.spearmanr(ok["dN"], ok["tmscore"])
    axes[0].set_xlabel(r"$d_N$ (human vs mouse)")
    axes[0].set_ylabel("TM-score (AlphaFold models)")
    axes[0].set_title(f"sequence vs structure divergence\nSpearman "
                      r"$\rho$" f"={rho:+.2f} (p={p:.1g})")
    axes[0].set_xscale("log")

    organs = [o for o in FOCUS_ORGANS if (res["organ"] == o).sum() >= 8]
    vals = [res.loc[res["organ"] == o, "tmscore"].dropna() for o in organs]
    if vals:
        bp = axes[1].boxplot(vals, tick_labels=organs, patch_artist=True,
                             showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(PALETTE[2])
            patch.set_alpha(0.75)
        for med in bp["medians"]:
            med.set_color("k")
    axes[1].set_ylabel("TM-score")
    axes[1].set_title("structural conservation by organ")
    axes[1].tick_params(axis="x", rotation=40)
    for lab in axes[1].get_xticklabels():
        lab.set_ha("right")

    # residual: structural conservation beyond what sequence divergence predicts
    if len(ok) > 20 and "plddt_min" in ok.columns and ok["plddt_min"].notna().any():
        ok = ok.dropna(subset=["plddt_min"])
        X = np.column_stack([np.log10(ok["dN"].clip(lower=1e-4)),
                             ok["plddt_min"], np.ones(len(ok))])
        y = ok["tmscore"].to_numpy()
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        ok = ok.assign(resid=y - X @ beta)
    elif len(ok) > 20:
        x = np.log10(ok["dN"].clip(lower=1e-4))
        y = ok["tmscore"]
        b = np.polyfit(x, y, 1)
        ok = ok.assign(resid=y - np.polyval(b, x))
        med = ok.groupby("organ")["resid"].median().sort_values()
        axes[2].barh(med.index, med.values, color=PALETTE[1])
        axes[2].axvline(0, color="k", lw=0.9)
        axes[2].set_xlabel("median TM-score residual\n"
                           r"(after $\log_{10} d_N$ AND mean pLDDT)")
        axes[2].set_title("structurally buffered (right)\nvs remodelled (left)")
        axes[2].grid(axis="y", alpha=0)

    p = os.path.join(FIGDIR, "t8_structure_vs_sequence.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(n_per_organ=30, force=False, threads=2):
    print("=" * 72)
    print("T8 - structural vs sequence divergence (Foldseek TM-align + AFDB)")
    print("=" * 72)
    out = os.path.join(DERIVED, "t8_structure.tsv")
    if os.path.exists(out) and not force:
        res = pd.read_csv(out, sep="\t", index_col=0)
        log(f"loaded cached structural comparison for {len(res)} pairs")
    else:
        sample = build_sample(n_per_organ=n_per_organ)
        log(f"sampled {len(sample)} orthologue pairs with SwissProt accessions")
        hits, usable = run_foldseek(sample, threads=threads)
        if hits is None:
            print("foldseek produced no usable output")
            return None
        res = add_plddt(join_pairs(hits, usable), usable)
        res.to_csv(out, sep="\t")

    print(f"orthologue pairs with both AlphaFold models: {len(res)}")
    print(f"median TM-score                             : {res['tmscore'].median():.4f}")
    print(f"median 3Di/structural identity (fident)     : {res['fident'].median():.4f}")
    print(f"median CDS-level codon identity             : {res['identity_cds'].median():.4f}")
    ok = res.dropna(subset=["dN", "tmscore"])
    rho, p = stats.spearmanr(ok["dN"], ok["tmscore"])
    print(f"Spearman(dN, TM-score)                      : {rho:+.3f} (p={p:.2g}, n={len(ok)})")
    if "plddt_min" in res.columns and res["plddt_min"].notna().any():
        okp = res.dropna(subset=["plddt_min", "tmscore"])
        rp, pp = stats.spearmanr(okp["plddt_min"], okp["tmscore"])
        print(f"Spearman(min mean pLDDT, TM-score)          : {rp:+.3f} "
              f"(p={pp:.2g}, n={len(okp)})")
        print("  -> a strong positive value means the TM-score axis is partly")
        print("     AlphaFold CONFIDENCE, not evolutionary structural change")

    print("\nper-organ structural conservation:")
    print(f"   {'organ':<18}{'n':>5}{'TM':>9}{'fident':>9}{'dN':>9}{'pLDDT':>9}")
    g = res.groupby("organ")
    for organ, sub in sorted(g, key=lambda kv: -kv[1]["tmscore"].median()):
        if len(sub) < 5:
            continue
        pl = sub["plddt_min"].median() if "plddt_min" in sub else float("nan")
        print(f"   {organ:<18}{len(sub):>5}{sub['tmscore'].median():>9.4f}"
              f"{sub['fident'].median():>9.4f}{sub['dN'].median():>9.4f}{pl:>9.1f}")
    p = make_figure(res)
    print(f"\nfigure written to: {p}")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-organ", type=int, default=30)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    report(n_per_organ=a.n_per_organ, force=a.force, threads=a.threads)
