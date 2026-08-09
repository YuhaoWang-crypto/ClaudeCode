#!/usr/bin/env python3
"""Cross-platform analysis of the ClairS HCC1395 demo calls.

Reads the three committed ClairS VCFs plus the SEQC2 truth subset, and answers:
  - which truth sites each platform recovers, and how they overlap
  - how well ClairS's allele fraction tracks the SEQC2 consensus VAF
  - whether the matched normal is clean (NAF) and reads are strand-balanced
  - what separates the one missed variant from the 28 recovered ones

Pure stdlib except matplotlib for the figures. Writes analysis.json and
figures/*.png next to this script.

    python3 clairs_demo/analyze.py
"""
import gzip
import json
import os
import statistics as st
from collections import Counter, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
TRUTH = os.path.join(HERE, "data", "truth_chr17_80.0-80.1Mb.vcf")
FIGDIR = os.path.join(HERE, "figures")

PLATFORMS = OrderedDict([
    ("ont",         {"label": "ONT R10.4.1",     "flag": "ont_r10_guppy"}),
    ("ilmn",        {"label": "Illumina NovaSeq", "flag": "ilmn"}),
    ("pacbio_hifi", {"label": "PacBio Revio",     "flag": "hifi_revio"}),
])
# Colour-blind-safe; kept consistent across every figure.
COLOURS = {"ont": "#4C72B0", "ilmn": "#DD8452", "pacbio_hifi": "#55A868"}
REGION = ("chr17", 80_000_000, 80_100_000)


# --------------------------------------------------------------------- parsing
def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(field):
    out = {}
    for item in field.split(";"):
        k, _, v = item.partition("=")
        out[k] = v if _ else True
    return out


def read_vcf(path):
    """-> list of dicts, one per record, with INFO and FORMAT/sample flattened."""
    records = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rec = {"chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": f[4],
                   "qual": f[5], "filter": f[6], "info": parse_info(f[7])}
            if len(f) > 9:
                rec["fmt"] = dict(zip(f[8].split(":"), f[9].split(":")))
            records.append(rec)
    return records


def in_region(rec):
    c, s, e = REGION
    return rec["chrom"] == c and s <= rec["pos"] <= e


def fnum(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ----------------------------------------------------------------- load inputs
truth = [r for r in read_vcf(TRUTH) if in_region(r)]
truth_by_pos = {r["pos"]: r for r in truth}

calls = {}
for key in PLATFORMS:
    recs = [r for r in read_vcf(os.path.join(RESULTS, key, "output.vcf.gz"))
            if in_region(r) and r["filter"] == "PASS"]
    calls[key] = {r["pos"]: r for r in recs}

called_sets = {k: set(v) for k, v in calls.items()}
truth_pos = set(truth_by_pos)

report = {
    "region": f"{REGION[0]}:{REGION[1]:,}-{REGION[2]:,}",
    "n_truth_snv": len(truth),
    "platforms": {},
}

# -------------------------------------------------------- 1. detection & overlap
for key, meta in PLATFORMS.items():
    got = called_sets[key]
    tp = sorted(got & truth_pos)
    report["platforms"][key] = {
        "label": meta["label"], "platform_flag": meta["flag"],
        "n_called": len(got), "tp": len(tp),
        "fp": len(got - truth_pos), "fn": len(truth_pos - got),
        "recall": round(len(tp) / len(truth_pos), 4),
        "precision": round(len(tp) / len(got), 4) if got else None,
    }

all_three = set.intersection(*called_sets.values())
any_one = set.union(*called_sets.values())
report["concordance"] = {
    "called_by_all_three": len(all_three),
    "called_by_at_least_one": len(any_one),
    "truth_missed_by_all": sorted(truth_pos - any_one),
    "per_site_platform_count": Counter(
        sum(p in called_sets[k] for k in PLATFORMS) for p in sorted(truth_pos)
    ),
}
report["concordance"]["per_site_platform_count"] = {
    str(k): v for k, v in sorted(report["concordance"]["per_site_platform_count"].items())
}

# ---------------------------------------------- 2. AF accuracy vs SEQC2 consensus
# SEQC2 TVAF is the consensus tumour VAF from ~60 short-read replicate call sets;
# treat it as the orthogonal reference, not as ground truth to the third decimal.
af_table = []
for pos in sorted(truth_pos):
    t = truth_by_pos[pos]
    row = {"pos": pos, "sub": f'{t["ref"]}>{t["alt"]}',
           "seqc2_tvaf": fnum(t["info"].get("TVAF")),
           "seqc2_nvaf": fnum(t["info"].get("NVAF")),
           "pacbio_tvaf": fnum(t["info"].get("PACB_TVAF")),
           "mq0": fnum(t["info"].get("MQ0")),
           "local_complexity": fnum(t["info"].get("LC"))}
    for key in PLATFORMS:
        rec = calls[key].get(pos)
        if rec:
            row[f"{key}_af"] = fnum(rec["fmt"].get("AF"))
            row[f"{key}_dp"] = fnum(rec["fmt"].get("DP"))
            row[f"{key}_naf"] = fnum(rec["fmt"].get("NAF"))
            row[f"{key}_ndp"] = fnum(rec["fmt"].get("NDP"))
            row[f"{key}_qual"] = fnum(rec["qual"])
        else:
            row[f"{key}_af"] = None
    af_table.append(row)
report["per_site"] = af_table


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = st.fmean(xs), st.fmean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
    return round(num / den, 4) if den else None


def _ranks(vals):
    """Average ranks, ties shared."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    """Rank correlation — the honest read when the VAF range is this narrow."""
    return pearson(_ranks(xs), _ranks(ys)) if len(xs) >= 3 else None


for key in PLATFORMS:
    pairs = [(r["seqc2_tvaf"], r[f"{key}_af"]) for r in af_table
             if r["seqc2_tvaf"] is not None and r.get(f"{key}_af") is not None]
    xs, ys = [p[0] for p in pairs], [p[1] for p in pairs]
    errs = [y - x for x, y in pairs]
    ratios = [y / x for x, y in pairs if x > 0]
    report["platforms"][key].update({
        "af_vs_seqc2_r": pearson(xs, ys),
        "af_vs_seqc2_spearman": spearman(xs, ys),
        "af_mean_signed_error": round(st.fmean(errs), 4),
        "af_median_ratio_to_seqc2": round(st.median(ratios), 3),
        "af_mean_abs_error": round(st.fmean([abs(e) for e in errs]), 4),
        "af_min": round(min(ys), 4), "af_median": round(st.median(ys), 4),
        "af_max": round(max(ys), 4),
    })

# --------------------------------- 3. normal cleanliness, depth, QUAL, strand, H
for key in PLATFORMS:
    recs = list(calls[key].values())
    nafs = [fnum(r["fmt"].get("NAF"), 0.0) for r in recs]
    dps = [fnum(r["fmt"].get("DP"), 0.0) for r in recs]
    ndps = [fnum(r["fmt"].get("NDP"), 0.0) for r in recs]
    quals = [fnum(r["qual"], 0.0) for r in recs]

    # alt-supporting reads on the forward strand / all alt reads
    strand = []
    for r in recs:
        alt = r["alt"]
        fwd = fnum(r["info"].get(f"F{alt}U"))
        rev = fnum(r["info"].get(f"R{alt}U"))
        if fwd is not None and rev is not None and (fwd + rev) > 0:
            strand.append(fwd / (fwd + rev))

    report["platforms"][key].update({
        "normal_af_max": round(max(nafs), 4),
        "normal_af_nonzero_sites": sum(1 for v in nafs if v > 0),
        "tumour_depth_median": round(st.median(dps), 1),
        "normal_depth_median": round(st.median(ndps), 1),
        "qual_min": round(min(quals), 2), "qual_median": round(st.median(quals), 2),
        "qual_max": round(max(quals), 2),
        "alt_forward_fraction_median": round(st.median(strand), 3) if strand else None,
        "alt_forward_fraction_range": [round(min(strand), 3), round(max(strand), 3)] if strand else None,
        "haplotype_tagged_calls": sum(1 for r in recs if "H" in r["info"]),
    })

# ------------------------------------------- 4. substitution spectrum + context
COMP = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


context = {}
ctx_path = os.path.join(HERE, "data", "truth_context.tsv")
if os.path.exists(ctx_path):
    with open(ctx_path) as fh:
        next(fh)
        for line in fh:
            pos, ref, alt, ctx = line.rstrip("\n").split("\t")
            context[int(pos)] = (ref, alt, ctx)

spectrum, tcw = Counter(), Counter()
for r in truth:
    ref, alt = r["ref"], r["alt"]
    ctx = context.get(r["pos"], (None, None, None))[2]
    if ref in ("A", "G"):          # collapse to pyrimidine reference
        ref, alt = COMP[ref], COMP[alt]
        ctx = revcomp(ctx) if ctx else None
    sub = f"{ref}>{alt}"
    spectrum[sub] += 1
    # APOBEC (SBS2/SBS13) acts on TCW = TCA/TCT, producing C>T and C>G
    if ctx and ref == "C":
        tcw["TCW" if (ctx[0] == "T" and ctx[2] in "AT") else "other C"] += 1
report["substitution_spectrum"] = dict(sorted(spectrum.items()))
report["c_site_context"] = dict(sorted(tcw.items())) if tcw else None

# ------------------------------------------ 5. run intermediates (see collect_run_stats.py)
derived_path = os.path.join(HERE, "data", "run_derived.json")
derived = json.load(open(derived_path)) if os.path.exists(derived_path) else {}
report["run_derived"] = derived
for key, d in derived.items():
    gn, gt = d.get("germline_normal"), d.get("germline_tumour")
    if gn and gt and gn["het"]:
        report["platforms"][key]["germline_het_normal"] = gn["het"]
        report["platforms"][key]["germline_het_tumour"] = gt["het"]
        report["platforms"][key]["germline_het_retained_in_tumour"] = round(
            gt["het"] / gn["het"], 3)

# ------------------------------------------------------------ 6. the missed site
missed = sorted(truth_pos - called_sets["ont"])
report["missed_by_long_reads"] = []
for pos in missed:
    t = truth_by_pos[pos]
    report["missed_by_long_reads"].append({
        "pos": pos, "sub": f'{t["ref"]}>{t["alt"]}',
        "seqc2_tvaf": fnum(t["info"].get("TVAF")),
        "seqc2_nvaf": fnum(t["info"].get("NVAF")),
        "seqc2_tvaf95": t["info"].get("TVAF95"),
        "pacbio_tvaf": fnum(t["info"].get("PACB_TVAF")),
        "n_passes": fnum(t["info"].get("nPASSES")),
        "called_by": [k for k in PLATFORMS if pos in called_sets[k]],
        "illumina_af": (calls["ilmn"][pos]["fmt"]["AF"] if pos in calls["ilmn"] else None),
        "illumina_ad": (calls["ilmn"][pos]["fmt"]["AD"] if pos in calls["ilmn"] else None),
        "illumina_nad": (calls["ilmn"][pos]["fmt"]["NAD"] if pos in calls["ilmn"] else None),
    })
# rank of the missed site's VAF among all truth VAFs
tvafs = sorted(r["seqc2_tvaf"] for r in af_table if r["seqc2_tvaf"] is not None)
for m in report["missed_by_long_reads"]:
    report["concordance"]["missed_site_vaf_rank"] = (
        f'{tvafs.index(m["seqc2_tvaf"]) + 1} of {len(tvafs)} (lowest = 1)')

# --------------------------------------------------------------------- figures
def make_figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIGDIR, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 160})

    # Fig 1 — ClairS AF vs SEQC2 consensus VAF
    fig, ax = plt.subplots(figsize=(5.0, 4.4))
    lim = 0.5
    ax.plot([0, lim], [0, lim], color="#999", lw=1, ls="--", zorder=1,
            label="AF = SEQC2 TVAF")
    # 27 of 29 truth VAFs sit inside this band — the r values are driven by the
    # two points outside it, which is the whole point of the figure.
    ax.axvspan(0.178, 0.214, color="#F0F0F0", zorder=0)
    ax.text(0.196, 0.335, "27/29 truth SNVs\nVAF 0.178–0.214", ha="center",
            fontsize=7, color="#777")
    for key, meta in PLATFORMS.items():
        pts = [(r["seqc2_tvaf"], r[f"{key}_af"]) for r in af_table
               if r["seqc2_tvaf"] is not None and r.get(f"{key}_af") is not None]
        d = report["platforms"][key]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, alpha=.85,
                   color=COLOURS[key], edgecolor="white", linewidth=.5,
                   label=f'{meta["label"]}  r={d["af_vs_seqc2_r"]}, ρ={d["af_vs_seqc2_spearman"]}',
                   zorder=3)
    for m in report["missed_by_long_reads"]:
        ax.axvline(m["seqc2_tvaf"], color="#C44E52", lw=1, ls=":", zorder=2)
        ax.annotate(f'missed by long reads\nTVAF={m["seqc2_tvaf"]:.3f}',
                    xy=(m["seqc2_tvaf"], lim * .93), xytext=(m["seqc2_tvaf"] + .03, lim * .93),
                    color="#C44E52", fontsize=7.5, va="top")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("SEQC2 consensus tumour VAF")
    ax.set_ylabel("ClairS reported AF")
    ax.set_title("Per-site AF agreement is unmeasurable in this window", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "af_vs_seqc2.png")); plt.close(fig)

    # Fig 2 — per-site detection matrix, sites ordered by SEQC2 VAF
    order = sorted([r for r in af_table if r["seqc2_tvaf"] is not None],
                   key=lambda r: r["seqc2_tvaf"])
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.2, 3.9), sharex=True, height_ratios=[2, 1.5],
        gridspec_kw={"hspace": .12}, layout="constrained")
    xs = range(len(order))
    ax0.bar(xs, [r["seqc2_tvaf"] for r in order], color="#BBB", width=.72)
    ax0.set_ylabel("SEQC2 VAF")
    ax0.set_title("All 29 truth SNVs, ordered by allele fraction", fontsize=10)
    for i, key in enumerate(PLATFORMS):
        y = len(PLATFORMS) - 1 - i
        hit = [j for j, r in enumerate(order) if r.get(f"{key}_af") is not None]
        miss = [j for j, r in enumerate(order) if r.get(f"{key}_af") is None]
        ax1.scatter(hit, [y] * len(hit), marker="s", s=34, color=COLOURS[key])
        ax1.scatter(miss, [y] * len(miss), marker="x", s=34, color="#C44E52")
    ax1.set_yticks(range(len(PLATFORMS)))
    ax1.set_yticklabels([PLATFORMS[k]["label"] for k in reversed(PLATFORMS)], fontsize=8)
    ax1.set_ylim(-.6, len(PLATFORMS) - .4)
    ax1.set_xlabel("truth SNV (lowest VAF → highest)")
    ax1.grid(axis="x", color="#EEE", lw=.5)
    fig.savefig(os.path.join(FIGDIR, "detection_by_vaf.png")); plt.close(fig)

    # Fig 3 — AF / depth / QUAL spread per platform
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9))
    series = [("AF", "af"), ("tumour depth", "dp"), ("QUAL", "qual")]
    for ax, (title, field) in zip(axes, series):
        for i, key in enumerate(PLATFORMS):
            vals = [r[f"{key}_{field}"] for r in af_table
                    if r.get(f"{key}_{field}") is not None]
            ax.scatter([i + (j % 5 - 2) * .045 for j in range(len(vals))], vals,
                       s=16, alpha=.75, color=COLOURS[key], edgecolor="none")
            ax.plot([i - .22, i + .22], [st.median(vals)] * 2, color="#333", lw=1.6)
        ax.set_xticks(range(len(PLATFORMS)))
        ax.set_xticklabels(["ONT", "ILMN", "HiFi"], fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel("value (bar = median)", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "distributions.png")); plt.close(fig)

    # Fig 4 — mutational spectrum (APOBEC) and the LOH signature
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.4, 2.9))
    order6 = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    sig = ["#3BC8F0", "#000000", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5"]  # COSMIC
    axL.bar(order6, [spectrum.get(s, 0) for s in order6], color=sig,
            edgecolor="#444", linewidth=.4)
    axL.set_ylabel("truth SNVs")
    axL.set_title("Substitution spectrum (n=29)", fontsize=9)
    if report.get("c_site_context"):
        tot = sum(report["c_site_context"].values())
        pct = 100 * report["c_site_context"].get("TCW", 0) / tot
        axL.text(.97, .95, f'{pct:.0f}% of C-site\nmutations at TCW\n(background 13.9%)',
                 transform=axL.transAxes, ha="right", va="top", fontsize=7.5, color="#333")

    labels, normal, tumour = [], [], []
    for key in PLATFORMS:
        d = derived.get(key, {})
        if d.get("germline_normal") and d.get("germline_tumour"):
            labels.append(PLATFORMS[key]["label"])
            normal.append(d["germline_normal"]["het"])
            tumour.append(d["germline_tumour"]["het"])
    if labels:
        x = range(len(labels))
        axR.bar([i - .19 for i in x], normal, width=.38, label="normal (HCC1395BL)",
                color="#8C8C8C")
        axR.bar([i + .19 for i in x], tumour, width=.38, label="tumour (HCC1395)",
                color="#C44E52")
        axR.set_xticks(list(x)); axR.set_xticklabels(labels, fontsize=8)
        axR.set_ylabel("germline het SNPs in window")
        axR.set_title("Loss of heterozygosity", fontsize=9)
        axR.legend(frameon=False, fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "spectrum_and_loh.png")); plt.close(fig)

    return ["af_vs_seqc2.png", "detection_by_vaf.png", "distributions.png",
            "spectrum_and_loh.png"]


if __name__ == "__main__":
    try:
        report["figures"] = make_figures()
    except ImportError:
        report["figures"] = []
        print("[warn] matplotlib missing — skipping figures")

    with open(os.path.join(HERE, "analysis.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    p = report["platforms"]
    print(f'Region {report["region"]}  |  {report["n_truth_snv"]} truth SNVs\n')
    hdr = f'{"platform":<18}{"TP":>4}{"FP":>4}{"FN":>4}{"recall":>9}{"AF r":>8}{"AF bias":>10}{"medDP":>8}'
    print(hdr); print("-" * len(hdr))
    for k, m in PLATFORMS.items():
        d = p[k]
        print(f'{m["label"]:<18}{d["tp"]:>4}{d["fp"]:>4}{d["fn"]:>4}{d["recall"]:>9.4f}'
              f'{d["af_vs_seqc2_r"]:>8}{d["af_mean_signed_error"]:>+10.4f}{d["tumour_depth_median"]:>8}')
    print(f'\nCalled by all three: {report["concordance"]["called_by_all_three"]}'
          f' / {report["n_truth_snv"]}')
    for m in report["missed_by_long_reads"]:
        print(f'Missed by long reads: chr17:{m["pos"]:,} {m["sub"]}  '
              f'SEQC2 TVAF={m["seqc2_tvaf"]}  Illumina AF={m["illumina_af"]}')
    print(f'Substitution spectrum: {report["substitution_spectrum"]}')
    print(f'\nWrote analysis.json and figures/ in {HERE}')
