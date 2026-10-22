"""
Module 23 — Which paralogue in a fiber is evolutionarily indispensable? (GISMO)

M1/M11/M14 collapse paralogue families (HRAS/KRAS/NRAS, MAP2K1/2, MAPK1/3 ...)
into one fiber representative: topologically they are interchangeable. Biology
breaks that symmetry - usually one member is essential, the others are not.
This module asks whether *deep-mammalian conservation* (GISMO / GISMO-mis,
Liao et al. Nat Genet 2026) recovers the essential member, and how that
compares with the human-population metric LOEUF (gnomAD v2.1.1).

Truth: mouse knockout lethality (gnomAD benchmark list; mouse symbols
upper-cased, ⚠️ approximate ortholog mapping).
Test: for every multi-member layer of the M14 atlas whose members DISAGREE on
lethality, is the top-ranked (most constrained) member a lethal one?
Baseline = expected hit rate when picking a member at random.

Rigor: ✅ scores and the hit counts are computed; ⚠️ n is small (tens of
families), so treat the comparison as indicative, with a binomial p-value.
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from grn_pipeline.m14_atlas import PATHWAYS


def families():
    rows = []
    for pw, spec in PATHWAYS.items():
        for layer, members in spec["layers"]:
            if len(members) >= 2:
                rows.append((pw, layer, list(members)))
    return rows


def analyse():
    from gismo import paths
    from gismo.benchmark import _gene_list
    from gismo.core import load_scores
    d = load_scores().drop_duplicates("symbol").set_index("symbol")
    lethal = _gene_list("mouse_lethal.tsv", upper=True)
    lo = None
    p = paths.REF / "gnomad_v2.1.1_lof.txt.bgz"
    if p.exists():
        v2 = pd.read_csv(p, sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"])
        lo = v2.dropna().sort_values("oe_lof_upper").drop_duplicates("gene").set_index("gene").oe_lof_upper

    fam_rows, tests = [], {"GISMO": [], "GISMO-mis": [], "LOEUF": []}
    for pw, layer, mem in families():
        rec = []
        for g in mem:
            rec.append({"gene": g, "lethal": g in lethal,
                        "GISMO": d.gismo.get(g, np.nan),
                        "GISMO-mis": d.gismo_mis.get(g, np.nan),
                        "LOEUF": np.nan if lo is None else lo.get(g, np.nan)})
        t = pd.DataFrame(rec)
        disc = t.lethal.any() and not t.lethal.all()
        row = {"pathway": pw, "layer": layer, "members": ",".join(mem),
               "lethal": ",".join(t.gene[t.lethal]) or "-", "discordant": disc}
        for m in tests:
            s = t.dropna(subset=[m])
            if len(s) >= 2:
                top = s.sort_values(m).iloc[0]
                row[f"top_{m}"] = top.gene
                if disc and s.lethal.any() and not s.lethal.all():
                    tests[m].append((bool(top.lethal), float(s.lethal.mean())))
        fam_rows.append(row)

    summary = {}
    for m, v in tests.items():
        if not v:
            continue
        hits = sum(h for h, _ in v); n = len(v); base = float(np.mean([b for _, b in v]))
        summary[m] = {"n_families": n, "hits": int(hits), "hit_rate": round(hits / n, 3),
                      "random_baseline": round(base, 3),
                      "p_binom": float(binomtest(hits, n, base, alternative="greater").pvalue)}
    return pd.DataFrame(fam_rows), summary


def report(fig_path="figures/m23_gismo_paralogs.png"):
    print("M23  GISMO paralogue essentiality — which fiber member is indispensable?")
    fam, summ = analyse()
    disc = fam[fam.discordant]
    cols = ["pathway", "layer", "members", "lethal", "top_GISMO", "top_GISMO-mis", "top_LOEUF"]
    cols = [c for c in cols if c in disc.columns]
    print(disc[cols].to_string(index=False))
    for m, s in summ.items():
        print(f"  {m:9s}: top member is mouse-lethal in {s['hits']}/{s['n_families']} "
              f"families ({s['hit_rate']:.0%}; random {s['random_baseline']:.0%}; "
              f"p={s['p_binom']:.3g})")
    _figure(summ, fig_path)
    return {"families": fam.to_dict("records"), "summary": summ}


def _figure(summ, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ms = list(summ)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    x = np.arange(len(ms))
    ax.bar(x - 0.18, [summ[m]["hit_rate"] for m in ms], 0.36, label="most-constrained member is lethal",
           color=["#b5838d", "#6a9fb5", "#c8b273"][:len(ms)])
    ax.bar(x + 0.18, [summ[m]["random_baseline"] for m in ms], 0.36, color="#cccccc", label="random member")
    for i, m in enumerate(ms):
        ax.text(i - 0.18, summ[m]["hit_rate"] + 0.02, f"{summ[m]['hits']}/{summ[m]['n_families']}",
                ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(ms); ax.set_ylim(0, 1.1)
    ax.set_ylabel("fraction of paralogue families")
    ax.set_title("M23: picking the essential paralogue in M14 fibers", fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(path, dpi=120)


if __name__ == "__main__":
    report()
