"""
T7 - ESM-2 constraint, used for the one thing dN/dS structurally cannot do.

What ESM is NOT good for here
-----------------------------
ESM-2 is trained on UniRef50 - the same set of natural sequences whose
divergence dN/dS measures. Its per-residue likelihood is a learned summary of
family-level conservation, so treating it as an INDEPENDENT estimate of
evolutionary rate is circular: agreement between ESM and omega is close to a
consistency check on ESM, not corroborating evidence about a tissue. This
module reports that correlation only as a calibration.

What ESM IS good for here
-------------------------
Every dN/dS number in T3-T6 requires a 1:1 orthologue. T4 measures how much of
each organ's transcriptome that excludes, and the excluded genes are not a
random sample - they are the lineage-specific, recently duplicated, fastest-
diverging ones, concentrated in exactly the organs the analysis is most
interested in. ESM scores a SINGLE sequence, so it reaches those genes.

Score: mean per-residue log P(observed residue | full sequence) from one forward
pass ("wt-marginal"). Higher = the residue the model expected = more constrained
family. Masked-marginals would be cleaner but need L forward passes per protein.

The bias that matters
---------------------
ESM assigns high likelihood to sequences from families it saw often. A primate-
specific orphan protein scores low partly because its family is small, not
necessarily because it is unconstrained - the same property that makes it
orthologue-less makes it ESM-hard. `family_breadth_bias()` quantifies this
directly by regressing the ESM score on the number of ladder species retaining
a 1:1 orthologue (0-4), and the residual effect is what the organ comparison
should be read against.
"""
import os
import gzip
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, SPECIES, ENSEMBL_FTP,
                     apply_figure_style, PALETTE, download, http_get,
                     read_fasta, log)
from . import t1_expression, t4_tissue_rates

MODEL = "facebook/esm2_t33_650M_UR50D"
MAX_LEN = 1022                      # ESM-2 context is 1024 incl. <cls>/<eos>
FOCUS_ORGANS = ["Testis", "Brain", "Cerebellum", "Liver", "Blood",
                "Skeletal_muscle", "Skin", "Pituitary", "Kidney", "Spleen"]


def human_proteins():
    """gene id -> longest protein sequence."""
    path = os.path.join(RAW, "pep_homo_sapiens.fa.gz")
    if not os.path.exists(path):
        base = f"{ENSEMBL_FTP}/homo_sapiens/pep/"
        fn = [t.split("/")[-1] for t in http_get(base, timeout=300).split('"')
              if t.endswith(".pep.all.fa.gz")][0]
        download(base + fn, path)
    best = {}
    for header, seq in read_fasta(path):
        gene, biotype = None, None
        for f in header.split()[1:]:
            if f.startswith("gene:"):
                gene = f.split(":", 1)[1].split(".")[0]
            elif f.startswith("transcript_biotype:"):
                biotype = f.split(":", 1)[1]
        if gene is None or biotype != "protein_coding":
            continue
        seq = seq.rstrip("*")
        if "*" in seq or len(seq) < 50:
            continue
        if gene not in best or len(seq) > len(best[gene]):
            best[gene] = seq
    return best


def ortholog_breadth():
    """gene -> number of ladder species retaining a usable 1:1 orthologue (0-4)."""
    counts = {}
    n_species = 0
    for sp in SPECIES:
        f = os.path.join(DERIVED, f"t3_dnds_human_{sp['key']}.tsv.gz")
        if not os.path.exists(f):
            continue
        n_species += 1
        d = pd.read_csv(f, sep="\t", usecols=["human_gene", "n_codons"])
        for g in d.loc[d["n_codons"] >= 100, "human_gene"]:
            counts[g] = counts.get(g, 0) + 1
    return pd.Series(counts, dtype="float64"), n_species


def build_sample(n_per_group=50, seed=0):
    """Stratified sample: per focus organ, tissue-enriched genes WITH and
    WITHOUT a 1:1 mouse orthologue."""
    feat, _ = t1_expression.build()
    prot = human_proteins()
    breadth, n_species = ortholog_breadth()

    feat = feat[feat.index.isin(prot)].copy()
    feat["n_ortholog_species"] = breadth.reindex(feat.index).fillna(0.0)
    feat["has_ortholog"] = feat["n_ortholog_species"] > 0
    feat["prot_len"] = pd.Series({g: len(prot[g]) for g in feat.index})

    rng = np.random.default_rng(seed)
    picks = []
    for organ in FOCUS_ORGANS:
        sub = feat[feat["enriched_in"] == organ]
        for has in (True, False):
            pool = sub[sub["has_ortholog"] == has]
            if len(pool) == 0:
                continue
            take = pool.index.to_numpy()
            if len(take) > n_per_group:
                take = rng.choice(take, n_per_group, replace=False)
            for g in take:
                picks.append(dict(gene=g, organ=organ, has_ortholog=has))
    # a genome-wide background of broadly expressed genes for reference
    bg = feat[(feat["enriched_in"].isna()) & (feat["tau"] < 0.5)]
    take = rng.choice(bg.index.to_numpy(), min(n_per_group * 3, len(bg)),
                      replace=False)
    for g in take:
        picks.append(dict(gene=g, organ="_background_broad",
                          has_ortholog=bool(feat.loc[g, "has_ortholog"])))

    s = pd.DataFrame(picks).drop_duplicates("gene").set_index("gene")
    s = s.join(feat[["symbol", "tau", "mean_tpm", "n_ortholog_species", "prot_len"]])
    s["seq"] = [prot[g] for g in s.index]
    return s, n_species


def esm_scores(seqs, model_name=MODEL, log_every=25):
    """Mean per-residue wt-marginal log-probability, plus mean entropy."""
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    torch.set_num_threads(os.cpu_count() or 4)
    tok = AutoTokenizer.from_pretrained(model_name)
    mod = AutoModelForMaskedLM.from_pretrained(model_name).eval()

    out = []
    for i, (gene, seq) in enumerate(seqs.items()):
        s = seq[:MAX_LEN]
        with torch.no_grad():
            enc = tok(s, return_tensors="pt")
            logits = mod(**enc).logits[0]
        logp = torch.log_softmax(logits.float(), dim=-1)
        ids = enc["input_ids"][0]
        # drop <cls> and <eos>
        sel = slice(1, 1 + len(s))
        wt = logp[sel, :].gather(1, ids[sel].unsqueeze(1)).squeeze(1)
        p = logp[sel, :].exp()
        ent = -(p * logp[sel, :]).sum(dim=1)
        out.append(dict(gene=gene, esm_logp=float(wt.mean()),
                        esm_entropy=float(ent.mean()),
                        esm_len_scored=len(s), truncated=len(seq) > MAX_LEN))
        if (i + 1) % log_every == 0:
            log(f"  ESM {i + 1}/{len(seqs)}")
    return pd.DataFrame(out).set_index("gene")


def family_breadth_bias(df):
    """Regress the ESM score on orthologue breadth - the confound that makes
    orphan genes look unconstrained for the wrong reason."""
    d = df.dropna(subset=["esm_logp", "n_ortholog_species"])
    rho, p = stats.spearmanr(d["esm_logp"], d["n_ortholog_species"])
    grouped = d.groupby("n_ortholog_species")["esm_logp"].agg(["median", "count"])
    return dict(rho=float(rho), p=float(p)), grouped


def calibration_vs_omega(df, sp_key=t4_tissue_rates.PRIMARY):
    f = os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz")
    if not os.path.exists(f):
        return None
    r = pd.read_csv(f, sep="\t").set_index("human_gene")
    d = df.join(r[["omega", "dN"]], how="inner").dropna(subset=["omega", "esm_logp"])
    if len(d) < 30:
        return None
    rho_o, p_o = stats.spearmanr(d["esm_logp"], d["omega"])
    rho_d, p_d = stats.spearmanr(d["esm_logp"], d["dN"])
    return dict(n=len(d), rho_omega=float(rho_o), p_omega=float(p_o),
                rho_dN=float(rho_d), p_dN=float(p_d))


def make_figure(df, grouped):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    d = df[df["organ"] != "_background_broad"]
    organs = [o for o in FOCUS_ORGANS if o in set(d["organ"])]
    pos = np.arange(len(organs))
    for k, (has, colour, off) in enumerate([(True, PALETTE[0], -0.18),
                                            (False, PALETTE[1], 0.18)]):
        vals = [d[(d["organ"] == o) & (d["has_ortholog"] == has)]["esm_logp"].dropna()
                for o in organs]
        bp = axes[0].boxplot(vals, positions=pos + off, widths=0.3,
                             patch_artist=True, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(colour)
            patch.set_alpha(0.75)
        for med in bp["medians"]:
            med.set_color("k")
    axes[0].set_xticks(pos)
    axes[0].set_xticklabels(organs, rotation=40, ha="right")
    axes[0].set_ylabel("mean ESM-2 wt log-probability")
    axes[0].set_title("blue = has 1:1 orthologue, orange = none")

    axes[1].plot(grouped.index, grouped["median"], "o-", color=PALETTE[2])
    axes[1].set_xlabel("ladder species retaining a 1:1 orthologue")
    axes[1].set_ylabel("median ESM-2 log-probability")
    axes[1].set_title("the family-breadth bias\n(orphans score low by construction)")

    bg = df[df["organ"] == "_background_broad"]["esm_logp"].dropna()
    tes = df[(df["organ"] == "Testis")]["esm_logp"].dropna()
    axes[2].hist(bg, bins=25, alpha=0.7, color=PALETTE[7],
                 density=True, label="broadly expressed")
    axes[2].hist(tes, bins=25, alpha=0.7, color=PALETTE[1],
                 density=True, label="testis-enriched")
    axes[2].set_xlabel("mean ESM-2 wt log-probability")
    axes[2].set_ylabel("density")
    axes[2].set_title("constraint distributions")
    axes[2].legend(fontsize=7)

    p = os.path.join(FIGDIR, "t7_esm_constraint.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def report(n_per_group=50, force=False):
    print("=" * 72)
    print("T7 - ESM-2 constraint, including genes dN/dS cannot reach")
    print("=" * 72)
    out = os.path.join(DERIVED, "t7_esm_scores.tsv.gz")
    if os.path.exists(out) and not force:
        df = pd.read_csv(out, sep="\t", index_col=0)
        log(f"loaded cached ESM scores for {len(df)} proteins")
    else:
        sample, n_species = build_sample(n_per_group=n_per_group)
        log(f"scoring {len(sample)} proteins with {MODEL} "
            f"(ladder species available: {n_species})")
        scores = esm_scores(sample["seq"])
        df = sample.drop(columns=["seq"]).join(scores)
        df.to_csv(out, sep="\t")

    print(f"proteins scored     : {len(df)}")
    print(f"truncated at {MAX_LEN} aa : {int(df['truncated'].sum())}")

    cal = calibration_vs_omega(df)
    if cal:
        print(f"\ncalibration (NOT independent evidence - ESM trains on the same"
              f"\nsequences dN/dS is measured from):")
        print(f"   Spearman(ESM logp, omega) = {cal['rho_omega']:+.3f} "
              f"(p={cal['p_omega']:.2g}, n={cal['n']})")
        print(f"   Spearman(ESM logp, dN)    = {cal['rho_dN']:+.3f} "
              f"(p={cal['p_dN']:.2g})")

    bias, grouped = family_breadth_bias(df)
    print(f"\nfamily-breadth bias: Spearman(ESM logp, n orthologue species) = "
          f"{bias['rho']:+.3f} (p={bias['p']:.2g})")
    for k, r in grouped.iterrows():
        print(f"   {int(k)} species: median ESM logp {r['median']:+.3f}  "
              f"(n={int(r['count'])})")

    print("\nper-organ ESM constraint, split by orthologue availability:")
    print(f"   {'organ':<20}{'n_orth':>7}{'logp_orth':>11}{'n_none':>8}"
          f"{'logp_none':>11}{'delta':>8}")
    rows = []
    for organ in FOCUS_ORGANS + ["_background_broad"]:
        sub = df[df["organ"] == organ]
        a = sub[sub["has_ortholog"]]["esm_logp"].dropna()
        b = sub[~sub["has_ortholog"]]["esm_logp"].dropna()
        if len(a) < 5 and len(b) < 5:
            continue
        ma = a.median() if len(a) else np.nan
        mb = b.median() if len(b) else np.nan
        print(f"   {organ:<20}{len(a):>7}{ma:>11.3f}{len(b):>8}"
              f"{mb:>11.3f}{(mb - ma):>8.3f}")
        rows.append(dict(organ=organ, n_with=len(a), logp_with=ma,
                         n_without=len(b), logp_without=mb, delta=mb - ma))
    pd.DataFrame(rows).to_csv(os.path.join(DERIVED, "t7_organ_esm.tsv"),
                              sep="\t", index=False)
    p = make_figure(df, grouped)
    print(f"\nfigure written to: {p}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-group", type=int, default=50)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    report(n_per_group=a.n_per_group, force=a.force)
