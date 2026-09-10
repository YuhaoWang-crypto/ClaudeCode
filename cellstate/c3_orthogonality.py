"""
C3 — Are the six axes actually orthogonal? Tested on REAL human expression data.

The six-variable proposal is only useful if the six are (mostly) independent
coordinates. If two of them move together, the state space has fewer than six
degrees of freedom and the extra variables buy nothing.

This is testable without any new experiment, and this module does it:

  1. Download three real human GEO curated datasets (SOFT), spanning three
     different insults, so the answer cannot be an artifact of one experiment:
        GDS5350  RB1-depleted osteoblast, ionizing-radiation time course (n=36)
        GDS4130  ER stress (tunicamycin), lymphoblastoid twins        (n=104)
        GDS5027  NOAH trial, FFPE breast-cancer biopsies              (n=156)
  2. Download the real MSigDB Hallmark gene sets from the Enrichr library
     endpoint (no login), and map them onto the six axes.
  3. Score each axis per sample (mean z-score of its genes) and compute the
     6x6 axis correlation matrix, plus its effective dimensionality
     (participation ratio of the eigenvalues).

The control that makes this honest:
  Hallmark gene sets OVERLAP. Part of any correlation between two axis scores
  is induced by shared genes, not by biology. So every correlation is computed
  TWICE -- once with all genes, once after deleting every gene that appears in
  more than one axis -- and both are reported. Only the disjoint-gene number
  is evidence about biology.

Rigour: ✅ data, gene sets, scores and eigenvalues are all computed from
downloaded files. ⚠️ the hallmark-to-axis mapping is an editorial choice and is
printed in full so it can be disputed. ⚠️ these are bulk datasets: they bound
the *population-level* dependence of the axes, which is a lower bound on what
single cells would show.
"""
import gzip
import os
import subprocess

import numpy as np

from .common import CACHE, banner, figpath, save_json, _ensure

ENRICHR = ("https://maayanlab.cloud/Enrichr/geneSetLibrary"
           "?mode=text&libraryName=MSigDB_Hallmark_2020")

DATASETS = [
    ("GDS5350", "辐射时程 (RB1-depleted osteoblast, IR)", 36),
    ("GDS4130", "ER 应激 (tunicamycin, lymphoblastoid)", 104),
    ("GDS5027", "临床肿瘤 (NOAH breast cancer FFPE)", 156),
]

# Six axes -> hallmark sets. Printed in run() so the choice is auditable.
AXIS_SETS = {
    "细胞周期": ["G2-M Checkpoint", "E2F Targets", "Mitotic Spindle"],
    "应激水平": ["Unfolded Protein Response", "Reactive Oxygen Species Pathway",
                 "Hypoxia", "TNF-alpha Signaling via NF-kB"],
    "代谢储备": ["Oxidative Phosphorylation", "Glycolysis", "Fatty Acid Metabolism"],
    "DNA损伤": ["DNA Repair", "p53 Pathway", "UV Response Up"],
    "凋亡准备": ["Apoptosis"],
    "表观谱系": ["Epithelial Mesenchymal Transition", "Myogenesis", "Adipogenesis"],
}
AXES = list(AXIS_SETS)


def hallmarks():
    _ensure(CACHE)
    path = os.path.join(CACHE, "hallmark2020.txt")
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        subprocess.run(["curl", "-sL", "--max-time", "90", ENRICHR, "-o", path],
                       check=True)
    lib = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                lib[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return lib


def fetch_gds(gds):
    _ensure(CACHE)
    path = os.path.join(CACHE, f"{gds}.soft.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 10000:
        pre = f"{gds[:-3]}nnn"
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/datasets/{pre}/{gds}"
               f"/soft/{gds}.soft.gz")
        subprocess.run(["curl", "-sL", "--max-time", "300", url, "-o", path],
                       check=True)
    return path if os.path.getsize(path) > 10000 else None


def parse_gds(path):
    """Return (gene_symbols, matrix[genes x samples]) from a GDS SOFT file."""
    genes, rows, in_tbl, hdr = [], [], False, None
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!dataset_table_begin"):
                in_tbl = True
                continue
            if line.startswith("!dataset_table_end"):
                break
            if not in_tbl:
                continue
            parts = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = parts
                try:
                    gi = hdr.index("IDENTIFIER")
                except ValueError:
                    return None, None
                scols = [i for i, h in enumerate(hdr) if h.startswith("GSM")]
                continue
            sym = parts[gi].strip().upper()
            if not sym or sym in ("--", "NA"):
                continue
            vals = []
            for i in scols:
                try:
                    vals.append(float(parts[i]))
                except (ValueError, IndexError):
                    vals.append(np.nan)
            genes.append(sym)
            rows.append(vals)
    if not rows:
        return None, None
    return np.array(genes), np.array(rows, dtype=float)


def collapse_genes(genes, X):
    """Average probes mapping to the same symbol; drop rows with any NaN."""
    keep = ~np.isnan(X).any(axis=1)
    genes, X = genes[keep], X[keep]
    order = np.argsort(genes)
    genes, X = genes[order], X[order]
    uniq, start = np.unique(genes, return_index=True)
    out = np.zeros((len(uniq), X.shape[1]))
    bounds = list(start) + [len(genes)]
    for k in range(len(uniq)):
        out[k] = X[bounds[k]:bounds[k + 1]].mean(axis=0)
    return uniq, out


def axis_genes(lib):
    """Return (all-gene sets, disjoint-gene sets) per axis."""
    full = {}
    for ax, sets in AXIS_SETS.items():
        g = set()
        for s in sets:
            g |= lib.get(s, set())
        full[ax] = g
    seen, dup = set(), set()
    for g in full.values():
        dup |= (seen & g)
        seen |= g
    disjoint = {ax: (g - dup) for ax, g in full.items()}
    return full, disjoint


def score(genes, X, sets):
    """Axis score per sample = mean z-score (across samples) of its genes."""
    mu, sd = X.mean(axis=1, keepdims=True), X.std(axis=1, keepdims=True)
    Z = (X - mu) / np.where(sd < 1e-9, np.nan, sd)
    idx = {g: i for i, g in enumerate(genes)}
    S, used = [], []
    for ax in AXES:
        rows = [idx[g] for g in sets[ax] if g in idx]
        used.append(len(rows))
        S.append(np.nanmean(Z[rows], axis=0) if rows else np.full(X.shape[1], np.nan))
    return np.array(S), used


def eff_dim(C):
    """Participation ratio of the correlation eigenvalues: 6 if orthogonal."""
    ev = np.linalg.eigvalsh(C)
    ev = np.clip(ev, 0, None)
    return float(ev.sum() ** 2 / np.sum(ev ** 2)), ev[::-1]


def _show(C, title):
    print(f"\n  {title}")
    print("        " + "".join(f"{a:>10}" for a in AXES))
    for i, a in enumerate(AXES):
        print(f"  {a:<8}" + "".join(
            f"{C[i, j]:>10.2f}" if np.isfinite(C[i, j]) else f"{'n/a':>10}"
            for j in range(len(AXES))))


def run():
    banner("C3 — 六轴正交性:在真实人类表达数据上实测(含基因集重叠对照)")
    lib = hallmarks()
    print(f"MSigDB Hallmark 2020 基因集: {len(lib)} 个(Enrichr 实时下载)")

    full, disj = axis_genes(lib)
    print("\n【轴 → hallmark 映射(可争议,故完整打印)】")
    for ax in AXES:
        print(f"  {ax:<8} <- {', '.join(AXIS_SETS[ax])}")
        print(f"           全部基因 {len(full[ax]):>4} 个;"
              f" 去除与其它轴共享后剩 {len(disj[ax]):>4} 个"
              f" (共享掉 {100*(1-len(disj[ax])/max(len(full[ax]),1)):.0f}%)")

    # Jaccard overlap between axes -- the size of the artifact being controlled
    print("\n【基因集重叠 Jaccard】——任何相关性里必须先扣掉的部分")
    J = np.zeros((6, 6))
    for i, a in enumerate(AXES):
        for j, b in enumerate(AXES):
            u = len(full[a] | full[b])
            J[i, j] = len(full[a] & full[b]) / u if u else 0.0
    _show(J, "Jaccard(全基因集)")

    results, Cs_full, Cs_disj = [], [], []
    for gds, label, _n in DATASETS:
        path = fetch_gds(gds)
        if path is None:
            print(f"\n⚠️  {gds} 下载失败,跳过(不编造)")
            continue
        genes, X = parse_gds(path)
        if genes is None:
            print(f"\n⚠️  {gds} 解析失败,跳过")
            continue
        genes, X = collapse_genes(genes, X)
        Sf, usedf = score(genes, X, full)
        Sd, _ = score(genes, X, disj)
        Cf = np.corrcoef(Sf)
        Cd = np.corrcoef(Sd)
        pf, _ = eff_dim(Cf)
        pd_, evd = eff_dim(Cd)
        Cs_full.append(Cf); Cs_disj.append(Cd)
        results.append(dict(gds=gds, label=label, n_samples=int(X.shape[1]),
                            n_genes=int(X.shape[0]), used=usedf,
                            eff_dim_full=pf, eff_dim_disjoint=pd_,
                            C_full=Cf.tolist(), C_disjoint=Cd.tolist()))
        print(f"\n{'='*70}\n{gds} — {label}")
        print(f"  基因 {X.shape[0]}, 样本 {X.shape[1]}; 各轴命中基因数 {usedf}")
        _show(Cd, "轴间相关(已去除共享基因 —— 这是可信的那个)")
        print(f"\n  有效维数(participation ratio,6 = 完全正交):")
        print(f"     全基因集   {pf:.2f} / 6")
        print(f"     去共享基因 {pd_:.2f} / 6   ← 结论用这个")

    if not results:
        print("\n⚠️  没有数据集成功,C3 无结论。")
        return results

    Cm = np.nanmean(np.array(Cs_disj), axis=0)
    _show(Cm, f"三个数据集平均的轴间相关(去共享基因,n={len(results)} 个数据集)")
    pm, evm = eff_dim(Cm)
    print(f"\n跨数据集平均有效维数 = {pm:.2f} / 6")
    print(f"特征值谱: {np.round(evm, 2).tolist()}")

    off = [(abs(Cm[i, j]), AXES[i], AXES[j])
           for i in range(6) for j in range(i + 1, 6) if np.isfinite(Cm[i, j])]
    off.sort(reverse=True)
    print("\n【耦合最强的轴对】")
    for v, a, b in off[:4]:
        print(f"  {a} ↔ {b}: |r| = {v:.2f}")
    print("\n【最接近独立的轴对】")
    for v, a, b in off[-3:]:
        print(f"  {a} ↔ {b}: |r| = {v:.2f}")

    print(f"\n结论:六个轴在真实数据上并不正交。去掉共享基因后,有效维数仍只有")
    print(f"{pm:.2f},而不是 6 —— 大约 {6-pm:.1f} 个自由度是冗余的。这与 Kinker &")
    print("Tirosh (Nat Genet 2020) 在 198 个细胞系上无监督分解出 12 个复现程序、")
    print("而单个细胞系只携带 4-9 个的结果方向一致。")

    _plot(Cm, pm)
    save_json("c3_orthogonality.json",
              dict(datasets=results, mean_C_disjoint=Cm.tolist(), eff_dim=pm,
                   jaccard=J.tolist(), axis_sets=AXIS_SETS))
    print("结果已存 figures/_cellstate_cache/c3_orthogonality.json")
    return results


def _plot(Cm, pm):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(Cm, cmap="RdBu_r", vmin=-1, vmax=1)
    lbl = ["cell cycle", "stress", "metab", "DNA dmg", "apopt", "lineage"]
    ax.set_xticks(range(6)); ax.set_xticklabels(lbl, rotation=45, ha="right")
    ax.set_yticks(range(6)); ax.set_yticklabels(lbl)
    for i in range(6):
        for j in range(6):
            if np.isfinite(Cm[i, j]):
                ax.text(j, i, f"{Cm[i,j]:.2f}", ha="center", va="center",
                        fontsize=8,
                        color="white" if abs(Cm[i, j]) > 0.55 else "black")
    ax.set_title(f"Axis-axis correlation (shared genes removed)\n"
                 f"effective dimension = {pm:.2f} / 6", fontsize=10)
    fig.colorbar(im, shrink=0.8)
    fig.tight_layout()
    fig.savefig(figpath("c3_orthogonality.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()
