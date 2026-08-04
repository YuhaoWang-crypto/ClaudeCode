"""
T18 - one pathway, in full: the TCA cycle across strata, paralogues and structure.

Everything upstream aggregated over thousands of genes. This module does the
opposite: it takes ONE pathway whose chemistry is textbook and fixed, and asks
what the evolutionary record says about how that chemistry is implemented.

Three questions, each with a design that answers it.

1. WHY SO MANY PROTEINS FOR EIGHT REACTIONS?
   Reactome's TCA set has 35 members, not 8. They fall into strata that have
   nothing to do with each other functionally: the catalytic enzymes; the extra
   subunits of the multi-subunit steps; the Fe-S cluster biogenesis machinery
   that ACO2 and SDHB cannot work without; the assembly factors that build
   complex II; and post-translational regulators. If those strata are under
   different selective regimes, they should separate by rate - and that is a
   prediction, not a description.

2. PARALOGUES: THE TIGHTEST CONTROL AVAILABLE.
   Several steps exist twice - mitochondrial and cytosolic (ACO2/ACO1,
   MDH2/MDH1, IDH2/IDH1), or as tissue-specific subunit choices sharing a
   partner (SUCLA2 vs SUCLG2, both binding SUCLG1). A paralogue pair holds fold,
   chemistry and ancestry fixed, so comparing the two copies isolates what
   compartment and tissue allocation do to the rate. This is the same allocation
   question as T15 at the finest resolution the data allows.

3. WHERE ON THE PROTEIN DOES THE VARIATION SIT?
   Per-residue conservation across mammals from the TOGA alignments, mapped onto
   the AlphaFold model, against (a) relative solvent accessibility and (b)
   UniProt-annotated catalytic and ligand-binding residues. The prediction that
   would be boring if true and interesting if false: buried and catalytic
   residues are invariant, surface residues carry the change.

Data all re-fetched at small scale (~45 genes), so this module runs without the
1.2 GB TOGA download the genome-wide modules need.
"""
import os
import re
import gzip
import json
import argparse
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from .common import (RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE,
                     download, http_get, log)

TOGA_BASE = ("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/"
             "MultipleCodonAlignments/")
TOGA_LISTING = os.path.join(RAW, "toga_listing.html")
ALN_DIR = os.path.join(RAW, "tca_aln")
STRUCT_DIR = os.path.join(RAW, "tca_struct")
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/"
GNOMAD = os.path.join(RAW, "gnomad_constraint.txt.bgz")

# --------------------------------------------------------------------------
# Curated functional strata. Reactome gives membership, not role, and the whole
# first question is about role - so this mapping is hand-made from the
# biochemistry and is the one hand-curated object in the pipeline.
# --------------------------------------------------------------------------
STRATUM = {
    # the eight chemical steps
    "CS": "catalytic", "ACO2": "catalytic", "IDH2": "catalytic",
    "IDH3A": "catalytic", "OGDH": "catalytic", "SUCLG1": "catalytic",
    "SUCLA2": "catalytic", "SUCLG2": "catalytic", "SDHA": "catalytic",
    "SDHB": "catalytic", "FH": "catalytic", "MDH2": "catalytic",
    # extra subunits of the multi-subunit steps
    "IDH3B": "subunit", "IDH3G": "subunit", "DLST": "subunit",
    "DLD": "subunit", "SDHC": "subunit", "SDHD": "subunit",
    "MRPS36": "subunit",
    # Fe-S cluster biogenesis - ACO2 and SDHB are dead without it
    "FXN": "FeS_biogenesis", "NFS1": "FeS_biogenesis", "ISCU": "FeS_biogenesis",
    "ISCA1": "FeS_biogenesis", "ISCA2": "FeS_biogenesis",
    "LYRM4": "FeS_biogenesis",
    # complex II assembly
    "SDHAF1": "assembly", "SDHAF2": "assembly", "SDHAF3": "assembly",
    "SDHAF4": "assembly",
    # post-translational regulation / redox balancing
    "SIRT3": "regulator", "TRAP1": "regulator", "CSKMT": "regulator",
    "NNT": "regulator", "ACAT1": "regulator",
    # cytosolic / adjacent paralogues, added because the paralogue test needs them
    "ACO1": "cytosolic_paralog", "IDH1": "cytosolic_paralog",
    "MDH1": "cytosolic_paralog", "OGDHL": "cytosolic_paralog",
    "ACLY": "cytosolic_paralog", "GOT1": "cytosolic_paralog",
    "GOT2": "catalytic_adjacent", "PDHA1": "catalytic_adjacent",
    "PDHB": "catalytic_adjacent", "DLAT": "catalytic_adjacent",
    "SLC25A1": "transport", "SLC25A11": "transport",
}
CORE_STRATA = ["catalytic", "subunit", "FeS_biogenesis", "assembly", "regulator"]

# (mitochondrial / shared member, other member, what differs)
PARALOG_PAIRS = [
    ("ACO2", "ACO1", "mitochondrial vs cytosolic (ACO1 moonlights as IRP1)"),
    ("MDH2", "MDH1", "mitochondrial vs cytosolic"),
    ("IDH2", "IDH1", "mitochondrial vs cytosolic, both NADP"),
    ("IDH2", "IDH3A", "NADP- vs NAD-dependent, both mitochondrial"),
    ("OGDH", "OGDHL", "ubiquitous vs cerebellum-enriched paralogue"),
    ("SUCLG1", "SUCLA2", "shared alpha vs ATP-specific beta"),
    ("SUCLG1", "SUCLG2", "shared alpha vs GTP-specific beta"),
    ("SUCLA2", "SUCLG2", "ATP-specific (muscle) vs GTP-specific (liver) beta"),
    ("GOT2", "GOT1", "mitochondrial vs cytosolic"),
]

SPECIES = ["macaque", "mouse", "opossum", "chicken"]
STRUCT_FOCUS = ["CS", "ACO2", "IDH2", "FH", "MDH2", "SDHA", "SDHB", "OGDH",
                "SUCLG1", "SUCLA2", "IDH3A", "ACO1", "MDH1", "IDH1"]


# --------------------------------------------------------------------------
# gene table
# --------------------------------------------------------------------------
def gene_table():
    feat = pd.read_csv(os.path.join(DERIVED, "t1_gene_expression_features.tsv.gz"),
                       sep="\t", index_col=0)
    sym2g = {}
    for g, s in feat["symbol"].dropna().items():
        sym2g.setdefault(s, g)
    rows = []
    for sym, stratum in STRATUM.items():
        g = sym2g.get(sym)
        if g is None:
            continue
        rows.append(dict(gene=g, symbol=sym, stratum=stratum))
    d = pd.DataFrame(rows).set_index("gene")
    d = d.join(feat[["tau", "max_tpm", "top_organ", "enriched_in",
                     "group_enriched_in"]])
    for sp in SPECIES:
        p = os.path.join(DERIVED, f"t3_dnds_human_{sp}.tsv.gz")
        if not os.path.exists(p):
            continue
        r = pd.read_csv(p, sep="\t").set_index("human_gene")
        r = r[r["n_codons"] >= 100]
        d[f"omega_{sp}"] = r["omega"].reindex(d.index)
        d[f"dN_{sp}"] = r["dN"].reindex(d.index)
        if sp == "macaque" and "Nd" in r.columns:
            d["Nd"] = r["Nd"].reindex(d.index)
            d["Sd"] = r["Sd"].reindex(d.index)
    return d


def genome_baseline(sp="mouse"):
    r = pd.read_csv(os.path.join(DERIVED, f"t3_dnds_human_{sp}.tsv.gz"), sep="\t")
    return float(r[r["n_codons"] >= 100]["omega"].median())


# --------------------------------------------------------------------------
# 1. strata
# --------------------------------------------------------------------------
def stratum_test(d, sp="mouse"):
    col = f"omega_{sp}"
    rows = []
    for s in CORE_STRATA + ["cytosolic_paralog", "catalytic_adjacent", "transport"]:
        v = d.loc[d["stratum"] == s, col].dropna()
        if len(v) == 0:
            continue
        rows.append(dict(stratum=s, n=len(v), median_omega=float(v.median()),
                         median_tau=float(d.loc[d["stratum"] == s, "tau"].median())))
    res = pd.DataFrame(rows).set_index("stratum")
    cat = d.loc[d["stratum"] == "catalytic", col].dropna()
    acc = d.loc[d["stratum"].isin(["FeS_biogenesis", "assembly"]), col].dropna()
    if len(cat) >= 4 and len(acc) >= 4:
        u, p = stats.mannwhitneyu(cat, acc, alternative="two-sided")
        cmp_ = dict(n_catalytic=len(cat), n_accessory=len(acc),
                    median_catalytic=float(cat.median()),
                    median_accessory=float(acc.median()),
                    fold=float(acc.median() / cat.median()),
                    rank_biserial=float(2 * u / (len(cat) * len(acc)) - 1),
                    p=float(p))
    else:
        cmp_ = {}
    return res, cmp_


# --------------------------------------------------------------------------
# 2. paralogues
# --------------------------------------------------------------------------
def paralog_table(d, sp="mouse"):
    col = f"omega_{sp}"
    by_sym = d.reset_index().set_index("symbol")
    rows = []
    for a, b, note in PARALOG_PAIRS:
        if a not in by_sym.index or b not in by_sym.index:
            continue
        ra, rb = by_sym.loc[a], by_sym.loc[b]
        if not (np.isfinite(ra[col]) and np.isfinite(rb[col])):
            continue
        rows.append(dict(pair=f"{a} / {b}", note=note,
                         omega_a=float(ra[col]), omega_b=float(rb[col]),
                         log2_b_over_a=float(np.log2(rb[col] / ra[col])),
                         tau_a=float(ra["tau"]), tau_b=float(rb["tau"]),
                         organ_a=ra["top_organ"], organ_b=rb["top_organ"]))
    return pd.DataFrame(rows).set_index("pair")


# --------------------------------------------------------------------------
# 3. alignments, conservation, structure
# --------------------------------------------------------------------------
def toga_files_for(symbols):
    html = open(TOGA_LISTING).read()
    files = re.findall(r'(ENST[0-9]+\.([A-Za-z0-9._\-]+)\.fasta\.gz)', html)
    want = {}
    for fn, sym in files:
        if sym in symbols and sym not in want:
            want[sym] = fn
    return want


def fetch_alignments(symbols):
    os.makedirs(ALN_DIR, exist_ok=True)
    want = toga_files_for(symbols)

    def get(item):
        sym, fn = item
        dest = os.path.join(ALN_DIR, fn)
        if not os.path.exists(dest) or os.path.getsize(dest) == 0:
            try:
                download(TOGA_BASE + fn, dest, timeout=300, retries=3)
            except Exception:                             # noqa: BLE001
                return sym, None
        return sym, dest
    with ThreadPoolExecutor(max_workers=8) as ex:
        got = dict(ex.map(get, want.items()))
    got = {k: v for k, v in got.items() if v}
    log(f"TOGA alignments fetched for {len(got)}/{len(symbols)} symbols")
    return got


AA_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AA_OF = dict(zip(_CODONS, AA_TABLE))


def per_residue_conservation(path):
    """Fraction of species matching the human residue, per human position."""
    seqs = []
    ref = None
    with gzip.open(path, "rt") as fh:
        name, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    (seqs if name != "REFERENCE" else []).append("".join(chunks))
                    if name == "REFERENCE":
                        ref = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            if name == "REFERENCE":
                ref = "".join(chunks)
            else:
                seqs.append("".join(chunks))
    if ref is None or not seqs:
        return None

    def aa(seq, i):
        c = seq[3 * i:3 * i + 3].upper().replace("U", "T")
        return AA_OF.get(c)

    L = len(ref) // 3
    cons, ref_aa, ref_pos = [], [], []
    pos = 0
    for i in range(L):
        r = aa(ref, i)
        if r is None or r == "*":
            continue
        pos += 1                                  # 1-based ungapped human index
        obs = [aa(s, i) for s in seqs]
        obs = [o for o in obs if o is not None and o != "*"]
        if len(obs) < 30:
            continue
        cons.append(sum(1 for o in obs if o == r) / len(obs))
        ref_aa.append(r)
        ref_pos.append(pos)
    if not cons:
        return None
    return pd.DataFrame(dict(pos=ref_pos, aa=ref_aa, conservation=cons))


def uniprot_sites(symbols):
    """symbol -> (accession, set of catalytic/binding residue positions)."""
    out = {}
    q = "+OR+".join(f"gene_exact:{s}" for s in symbols)
    url = ("https://rest.uniprot.org/uniprotkb/stream?query=%28" + q +
           "%29+AND+organism_id:9606+AND+reviewed:true"
           "&fields=accession,gene_primary,ft_act_site,ft_binding,sequence"
           "&format=tsv")
    txt = http_get(url, timeout=300)
    for line in txt.rstrip("\n").split("\n")[1:]:
        f = line.split("\t")
        if len(f) < 5:
            continue
        acc, sym, act, bind, seq = f[0], f[1], f[2], f[3], f[4]
        if sym not in symbols or sym in out:
            continue
        pos = set()
        for m in re.finditer(r"(?:ACT_SITE|BINDING)\s+(\d+)(?:\.\.(\d+))?", act + " " + bind):
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) else a
            pos.update(range(a, b + 1))
        out[sym] = dict(accession=acc, sites=pos, length=len(seq))
    return out


def fetch_structure(acc):
    os.makedirs(STRUCT_DIR, exist_ok=True)
    dest = os.path.join(STRUCT_DIR, f"{acc}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    try:
        meta = json.loads(http_get(AFDB_API + acc, timeout=120, retries=2))
        url = meta[0].get("pdbUrl") if meta else None
        if not url:
            return None
        download(url, dest, timeout=300, retries=3)
        return dest
    except Exception:                                     # noqa: BLE001
        return None


def relative_accessibility(pdb_path):
    """Per-residue relative solvent accessibility via Shrake-Rupley."""
    from Bio.PDB import PDBParser, ShrakeRupley
    # Tien et al. 2013 theoretical maxima
    MAXASA = dict(A=129, R=274, N=195, D=193, C=167, E=223, Q=225, G=104,
                  H=224, I=197, L=201, K=236, M=224, F=240, P=159, S=155,
                  T=172, W=285, Y=263, V=174)
    THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E",
                 GLN="Q", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K",
                 MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W",
                 TYR="Y", VAL="V")
    st = PDBParser(QUIET=True).get_structure("x", pdb_path)
    ShrakeRupley().compute(st[0], level="R")
    rows = []
    for res in st[0].get_residues():
        one = THREE.get(res.get_resname())
        if one is None:
            continue
        rows.append(dict(pos=res.id[1], aa=one,
                         rsa=res.sasa / MAXASA[one],
                         plddt=float(np.mean([a.bfactor for a in res]))))
    return pd.DataFrame(rows)


def structure_analysis(symbols=STRUCT_FOCUS):
    alns = fetch_alignments(set(symbols))
    sites = uniprot_sites(set(symbols))
    rows = []
    for sym in symbols:
        if sym not in alns or sym not in sites:
            continue
        cons = per_residue_conservation(alns[sym])
        if cons is None:
            continue
        acc = sites[sym]["accession"]
        pdb = fetch_structure(acc)
        if pdb is None:
            continue
        rsa = relative_accessibility(pdb)
        m = cons.merge(rsa, on="pos", suffixes=("_aln", "_pdb"))
        # the TOGA transcript and the UniProt/AFDB sequence need not be the same
        # isoform; keep only positions where the residue identity agrees
        m = m[m["aa_aln"] == m["aa_pdb"]]
        if len(m) < 100:
            continue
        m["symbol"] = sym
        m["is_site"] = m["pos"].isin(sites[sym]["sites"])
        rows.append(m)
    if not rows:
        return None
    return pd.concat(rows, ignore_index=True)


def structure_stats(res):
    out = {}
    ok = res[res["plddt"] >= 70]
    rho, p = stats.spearmanr(ok["rsa"], ok["conservation"])
    out["rsa_vs_conservation"] = dict(rho=float(rho), p=float(p), n=len(ok))
    buried = ok[ok["rsa"] < 0.20]["conservation"]
    exposed = ok[ok["rsa"] > 0.50]["conservation"]
    u, p2 = stats.mannwhitneyu(buried, exposed, alternative="two-sided")
    out["buried_vs_exposed"] = dict(
        n_buried=len(buried), n_exposed=len(exposed),
        median_buried=float(buried.median()), median_exposed=float(exposed.median()),
        rank_biserial=float(2 * u / (len(buried) * len(exposed)) - 1), p=float(p2))
    site = ok[ok["is_site"]]["conservation"]
    non = ok[~ok["is_site"]]["conservation"]
    if len(site) >= 20:
        u3, p3 = stats.mannwhitneyu(site, non, alternative="two-sided")
        out["site_vs_rest"] = dict(
            n_sites=len(site), median_site=float(site.median()),
            median_rest=float(non.median()),
            frac_site_invariant=float((site >= 0.999).mean()),
            frac_rest_invariant=float((non >= 0.999).mean()),
            rank_biserial=float(2 * u3 / (len(site) * len(non)) - 1), p=float(p3))
        # sites are buried; the honest test compares them to RSA-matched controls
        matched = []
        rng = np.random.default_rng(0)
        s_rsa = ok[ok["is_site"]]["rsa"].to_numpy()
        pool = ok[~ok["is_site"]]
        for r in s_rsa:
            cand = pool[(pool["rsa"] - r).abs() < 0.05]
            if len(cand):
                matched.append(float(cand["conservation"].iloc[
                    rng.integers(len(cand))]))
        if len(matched) >= 20:
            u4, p4 = stats.mannwhitneyu(site, matched, alternative="two-sided")
            out["site_vs_rsa_matched"] = dict(
                n=len(matched), median_site=float(site.median()),
                median_matched=float(np.median(matched)),
                rank_biserial=float(2 * u4 / (len(site) * len(matched)) - 1),
                p=float(p4))
    return out


def species_retention(symbols):
    """How many mammalian assemblies retain each gene, from the TOGA alignments.

    "Is the cycle the same in every species" is a retention question before it is
    a rate question: TOGA only builds an alignment where an orthologue is called
    Intact / Partial Intact / Uncertain Loss, so the sequence count is a direct
    readout of how universally the gene survives across mammals.
    """
    alns = fetch_alignments(set(symbols))
    rows = []
    for sym, path in alns.items():
        n = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
        rows.append(dict(symbol=sym, n_species=n))
    return pd.DataFrame(rows).set_index("symbol").sort_values("n_species")


def depth_stability(d):
    """Does the stratum ordering hold at every divergence depth?"""
    rows = []
    for sp in SPECIES:
        col = f"omega_{sp}"
        if col not in d.columns:
            continue
        rec = dict(species=sp)
        for s in CORE_STRATA:
            v = d.loc[d["stratum"] == s, col].dropna()
            rec[s] = float(v.median()) if len(v) >= 3 else np.nan
        rows.append(rec)
    res = pd.DataFrame(rows).set_index("species")
    ranks = res.rank(axis=1)
    rho = {}
    keys = list(ranks.index)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ok = ranks.loc[a].notna() & ranks.loc[b].notna()
            if ok.sum() >= 3:
                rho[f"{a}~{b}"] = float(stats.spearmanr(ranks.loc[a][ok],
                                                        ranks.loc[b][ok]).statistic)
    return res, rho


def make_figure(d, strat, par, res, stats_out, sp="mouse"):
    plt = apply_figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    col = f"omega_{sp}"
    base = genome_baseline(sp)

    order = [s for s in CORE_STRATA if s in set(d["stratum"])]
    vals = [d.loc[d["stratum"] == s, col].dropna() for s in order]
    bp = axes[0].boxplot(vals, tick_labels=order, patch_artist=True,
                         showfliers=False)
    for patch, c in zip(bp["boxes"], [PALETTE[0], PALETTE[5], PALETTE[1],
                                      PALETTE[4], PALETTE[2]]):
        patch.set_facecolor(c)
        patch.set_alpha(0.8)
    for m in bp["medians"]:
        m.set_color("k")
    for i, s in enumerate(order, start=1):
        v = d.loc[d["stratum"] == s, col].dropna()
        axes[0].scatter(np.full(len(v), i) + np.random.uniform(-.09, .09, len(v)),
                        v, s=12, color="k", alpha=0.55, zorder=3)
    axes[0].axhline(base, color="k", ls="--", lw=0.9,
                    label=f"genome median {base:.3f}")
    axes[0].set_ylabel(r"$\omega$ (human-mouse)")
    axes[0].set_title("one pathway, five strata")
    axes[0].tick_params(axis="x", rotation=30)
    for lab in axes[0].get_xticklabels():
        lab.set_ha("right")
    axes[0].legend(fontsize=7)

    p = par.sort_values("log2_b_over_a")
    colours = [PALETTE[1] if v > 0 else PALETTE[0] for v in p["log2_b_over_a"]]
    axes[1].barh(p.index, p["log2_b_over_a"], color=colours)
    axes[1].axvline(0, color="k", lw=0.9)
    axes[1].set_xlabel(r"$\log_2(\omega_B / \omega_A)$")
    axes[1].set_title("paralogue pairs: same fold,\nsame chemistry, different post")
    axes[1].tick_params(axis="y", labelsize=6.5)
    axes[1].grid(axis="y", alpha=0)

    if res is not None:
        ok = res[res["plddt"] >= 70]
        bins = np.linspace(0, 1.0, 11)
        idx = np.digitize(ok["rsa"], bins[1:-1])
        med = [ok["conservation"][idx == i].median() for i in range(len(bins) - 1)]
        axes[2].plot(0.5 * (bins[:-1] + bins[1:]), med, "o-", color=PALETTE[0],
                     label="all residues")
        s = ok[ok["is_site"]]
        axes[2].axhline(float(s["conservation"].median()), color=PALETTE[1],
                        ls="--", lw=1.4,
                        label=f"catalytic/binding (n={len(s)})")
        axes[2].set_xlabel("relative solvent accessibility")
        axes[2].set_ylabel("fraction of mammals matching human")
        axes[2].set_title("where the variation sits")
        axes[2].legend(fontsize=7)

    path = os.path.join(FIGDIR, "t18_tca.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def report(sp="mouse", do_structure=True):
    print("=" * 72)
    print("T18 - the TCA cycle: strata, paralogues, structure")
    print("=" * 72)
    d = gene_table()
    base = genome_baseline(sp)
    print(f"genes in the curated TCA set: {len(d)}  "
          f"(genome median omega, human-{sp} = {base:.4f})\n")

    strat, cmp_ = stratum_test(d, sp)
    print("1. STRATA - Reactome lists 35 proteins for 8 reactions. Why?")
    print(f"   {'stratum':<20}{'n':>4}{'median omega':>14}{'vs genome':>11}{'median tau':>12}")
    for s, r in strat.iterrows():
        print(f"   {s:<20}{int(r['n']):>4}{r['median_omega']:>14.4f}"
              f"{r['median_omega'] / base:>11.2f}{r['median_tau']:>12.3f}")
    if cmp_:
        print(f"\n   catalytic enzymes vs Fe-S/assembly machinery: "
              f"{cmp_['median_catalytic']:.4f} vs {cmp_['median_accessory']:.4f}"
              f"  ({cmp_['fold']:.1f}x, rb={cmp_['rank_biserial']:+.2f}, "
              f"p={cmp_['p']:.2g})")

    print("\n   slowest and fastest members:")
    ok = d.dropna(subset=[f"omega_{sp}"]).sort_values(f"omega_{sp}")
    for _, r in ok.head(6).iterrows():
        print(f"      {r['symbol']:<9} {r[f'omega_{sp}']:.4f}  {r['stratum']}")
    print("      ...")
    for _, r in ok.tail(6).iterrows():
        print(f"      {r['symbol']:<9} {r[f'omega_{sp}']:.4f}  {r['stratum']}")

    par = paralog_table(d, sp)
    print("\n2. PARALOGUES - fold, chemistry and ancestry held fixed")
    print(f"   {'pair':<22}{'wA':>8}{'wB':>8}{'log2 B/A':>10}{'tauA':>7}{'tauB':>7}  note")
    for pair, r in par.iterrows():
        print(f"   {pair:<22}{r['omega_a']:>8.4f}{r['omega_b']:>8.4f}"
              f"{r['log2_b_over_a']:>+10.2f}{r['tau_a']:>7.2f}{r['tau_b']:>7.2f}"
              f"  {r['note'][:44]}")

    dep, rho = depth_stability(d)
    print("\n   stratum median omega at each divergence depth:")
    print("   " + dep.round(4).to_string().replace("\n", "\n   "))
    if rho:
        vals = list(rho.values())
        print(f"   pairwise Spearman of the stratum ordering across depths: "
              f"min {min(vals):+.2f}, median {np.median(vals):+.2f}")

    res, sout = None, {}
    if do_structure:
        print("\n3. STRUCTURE - per-residue conservation across mammals vs geometry")
        res = structure_analysis()
        if res is None:
            print("   no usable alignment/structure pairs")
        else:
            sout = structure_stats(res)
            n_sym = res["symbol"].nunique()
            print(f"   {len(res)} residues over {n_sym} enzymes "
                  f"(pLDDT>=70 kept for the tests)")
            a = sout["rsa_vs_conservation"]
            print(f"   Spearman(RSA, conservation) = {a['rho']:+.3f} "
                  f"(p={a['p']:.2g}, n={a['n']})")
            b = sout["buried_vs_exposed"]
            print(f"   buried (RSA<0.2) {b['median_buried']:.3f} vs "
                  f"exposed (RSA>0.5) {b['median_exposed']:.3f}  "
                  f"(rb={b['rank_biserial']:+.2f}, p={b['p']:.2g})")
            c = sout.get("site_vs_rest")
            if c:
                print(f"   catalytic/binding residues (n={c['n_sites']}): "
                      f"{c['median_site']:.3f} vs {c['median_rest']:.3f} elsewhere")
                print(f"      invariant across mammals: {c['frac_site_invariant']:.0%} "
                      f"of annotated sites vs {c['frac_rest_invariant']:.0%} elsewhere")
            e = sout.get("site_vs_rsa_matched")
            if e:
                print(f"   vs RSA-MATCHED controls (n={e['n']}): "
                      f"{e['median_site']:.3f} vs {e['median_matched']:.3f} "
                      f"(rb={e['rank_biserial']:+.2f}, p={e['p']:.2g})")
                print("      this is the test that matters - annotated sites are")
                print("      buried, so an unmatched comparison mostly measures burial")

    print("\n4. CROSS-SPECIES RETENTION - mammalian assemblies with an orthologue")
    ret = species_retention(set(d["symbol"]))
    d = d.join(ret, on="symbol")
    print(f"   {len(ret)} genes; median {int(ret['n_species'].median())} assemblies, "
          f"range {int(ret['n_species'].min())}-{int(ret['n_species'].max())}")
    print("   least universally retained:")
    for sym, r in ret.head(8).iterrows():
        st = d.loc[d["symbol"] == sym, "stratum"].iloc[0]
        print(f"      {sym:<9} {int(r['n_species']):>4} assemblies   {st}")
    by_str = d.groupby("stratum")["n_species"].median().sort_values()
    print("   median retention by stratum: "
          + ", ".join(f"{k} {int(v)}" for k, v in by_str.items()))

    d.to_csv(os.path.join(DERIVED, "t18_tca_genes.tsv"), sep="\t")
    ret.to_csv(os.path.join(DERIVED, "t18_retention.tsv"), sep="\t")
    dep.to_csv(os.path.join(DERIVED, "t18_depth_stability.tsv"), sep="\t")
    strat.to_csv(os.path.join(DERIVED, "t18_strata.tsv"), sep="\t")
    par.to_csv(os.path.join(DERIVED, "t18_paralogs.tsv"), sep="\t")
    if res is not None:
        res.to_csv(os.path.join(DERIVED, "t18_residues.tsv.gz"), sep="\t",
                   index=False, compression={"method": "gzip", "mtime": 0})
    p = make_figure(d, strat, par, res, sout, sp)
    print(f"\nfigure written to: {p}")
    return d, strat, par, res, sout


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="mouse")
    ap.add_argument("--no-structure", action="store_true")
    a = ap.parse_args()
    report(sp=a.species, do_structure=not a.no_structure)
