"""FineDNA demo: stream a sample of each training stage and look inside it.

    pip install -r finedna_demo/requirements.txt
    python finedna_demo/finedna_demo.py            # ~2-4 min, ~150 MB of network reads
    python finedna_demo/finedna_demo.py --n-dna 100 --n-rl 2000   # quicker

Nothing is downloaded in full (the dataset is ~1 TB). Rows are streamed and the
sample is cached under finedna_demo/cache/ so a re-run only redoes the analysis.
Outputs: finedna_demo/figures/*.png and finedna_demo/results/*.
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verifiers import check  # noqa: E402

REPO = "HuggingFaceBio/FineDNA"
DNA_CONFIGS = {  # label -> config
    "Mammal": "mid_dna_eukaryote_vertebrate_mammalian",
    "Invertebrate": "mid_dna_eukaryote_invertebrate",
    "Plant": "mid_dna_eukaryote_plant",
    "Prokaryote": "mid_dna_prokaryote",
    "mRNA": "mid_dna_mrna",
}
CACHE, FIG, RES = HERE / "cache", HERE / "figures", HERE / "results"

# Reference palette (dataviz skill, light mode); categorical slots in fixed order.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


# --- data access ------------------------------------------------------------

def sample(config, n, shuffle_buffer=0):
    """Stream n rows of one config, caching them as JSONL."""
    path = CACHE / f"{config}_{n}_{shuffle_buffer}.jsonl"
    if path.exists():
        return [json.loads(line) for line in path.open()]
    from datasets import load_dataset
    ds = load_dataset(REPO, config, split="train", streaming=True)
    if shuffle_buffer:
        ds = ds.shuffle(seed=0, buffer_size=shuffle_buffer)
    rows = list(ds.take(n))
    CACHE.mkdir(exist_ok=True)
    with path.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"  streamed {len(rows):>5} rows of {config}", flush=True)
    return rows


# --- stage 1: mid-training DNA ----------------------------------------------

DNA_RE = re.compile(r"((?:<[a-z_]+>)*)<dna>([^<]*)</dna>")


def dna_stats(rows):
    out = []
    for r in rows:
        for tags, s in DNA_RE.findall(r["text"]):
            if len(s) < 60:
                continue
            u = s.upper()
            c, g = u.count("C"), u.count("G")
            out.append({
                "id": r["id"], "len": len(s), "gc": 100 * (c + g) / len(s),
                "cpg_oe": u.count("CG") * len(s) / (c * g) if c and g else float("nan"),
                "mult6": len(s) % 6 == 0, "lower_frac": sum(ch.islower() for ch in s) / len(s),
                "n_frac": u.count("N") / len(s), "tags": re.findall(r"<([a-z_]+)>", tags),
            })
    return out


def kmer6_tokens(seq, k=6):
    """Carbon-X reads DNA as non-overlapping 6-mers."""
    return [seq[i:i + k] for i in range(0, len(seq) - len(seq) % k, k)]


# --- stage 2/3: SFT and RL --------------------------------------------------

def sft_stats(rows):
    dom = collections.Counter(r["domain"] for r in rows)
    turns = collections.defaultdict(list)
    examples = {}
    for r in rows:
        turns[r["domain"]].append(len(r["messages"]))
        examples.setdefault(r["domain"], r["messages"])
    return dom, {d: st.median(v) for d, v in turns.items()}, examples


def rl_verify(rows):
    details, per_type = [], collections.defaultdict(collections.Counter)
    for r in rows:
        if r["task_group"] != "basic_dna":
            continue
        prompt = r["prompt"][-1]["content"]
        ours, status = check(r["task_type"], prompt, r["answer_text"])
        per_type[r["task_type"]][status] += 1
        if status in ("agree", "disagree"):
            details.append({"task_type": r["task_type"], "status": status,
                            "gold": r["answer_text"], "ours": ours, "prompt": prompt})
    return per_type, details


# --- figures ----------------------------------------------------------------

def _style(ax, title):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _strip(ax, groups, key, title, xlabel, ref=None, log=False):
    import random
    rnd = random.Random(0)
    labels = list(groups)
    for i, lab in enumerate(labels):
        vals = [d[key] for d in groups[lab] if d[key] == d[key]]
        y = [i + rnd.uniform(-0.28, 0.28) for _ in vals]
        ax.scatter(vals, y, s=9, color=SERIES[i], alpha=0.45, linewidths=0)
        med = st.median(vals)
        ax.plot([med, med], [i - 0.36, i + 0.36], color=INK, linewidth=2)
        ax.annotate(f"{med:.2f}" if key == "cpg_oe" else f"{med:,.0f}", (med, i - 0.36),
                    color=INK, fontsize=8, ha="center", va="bottom", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.9))
    if ref is not None:
        ax.axvline(ref, color=INK2, linestyle="--", linewidth=1)
    if log:
        ax.set_xscale("log")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylim(len(labels) - 0.5, -0.85)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    _style(ax, title)


def plot_dna(groups):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), facecolor=SURFACE)
    _strip(axes[0], groups, "len", "Sequence length", "bp (log scale) · bar = median", log=True)
    _strip(axes[1], groups, "gc", "GC content", "% G+C · bar = median")
    _strip(axes[2], groups, "cpg_oe", "CpG observed / expected",
           "ratio · dashed = 1 (no CpG depletion)", ref=1.0)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    fig.suptitle("FineDNA mid-training DNA — one dot per sampled sequence",
                 x=0.01, ha="left", color=INK, fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "dna_composition.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def _hbar(ax, counter, title, color):
    items = counter.most_common()[::-1]
    labels, vals = [k for k, _ in items], [v for _, v in items]
    ax.barh(labels, vals, color=color, height=0.7)
    total = sum(vals)
    for i, v in enumerate(vals):
        ax.text(v, i, f"  {v} ({100 * v / total:.0f}%)", va="center", color=INK2, fontsize=8)
    ax.set_xlim(0, max(vals) * 1.25)
    _style(ax, title)


def plot_mix(sft_dom, rl_groups):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), facecolor=SURFACE)
    _hbar(axes[0], sft_dom, f"SFT domains (sft_mix, n={sum(sft_dom.values())})", SERIES[0])
    _hbar(axes[1], rl_groups, f"RL task groups (rl_mix, n={sum(rl_groups.values())})", SERIES[0])
    fig.tight_layout()
    fig.savefig(FIG / "sft_rl_mix.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


# --- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-dna", type=int, default=300, help="documents per DNA config")
    ap.add_argument("--n-sft", type=int, default=2000)
    ap.add_argument("--n-rl", type=int, default=5000)
    args = ap.parse_args()
    FIG.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)
    summary = {}

    print("[1/3] Mid-training DNA")
    groups = {}
    for label, cfg in DNA_CONFIGS.items():
        rows = sample(cfg, args.n_dna, shuffle_buffer=1000)
        groups[label] = dna_stats(rows)
    for label, recs in groups.items():
        tagc = collections.Counter(t for r in recs for t in r["tags"])
        summary.setdefault("dna", {})[label] = {
            "config": DNA_CONFIGS[label], "sequences": len(recs),
            "distinct_source_ids": len({r["id"] for r in recs}),
            "median_len_bp": st.median(r["len"] for r in recs),
            "median_gc_pct": round(st.median(r["gc"] for r in recs), 1),
            "median_cpg_oe": round(st.median(r["cpg_oe"] for r in recs if r["cpg_oe"] == r["cpg_oe"]), 3),
            "pct_len_multiple_of_6": round(100 * sum(r["mult6"] for r in recs) / len(recs), 1),
            "pct_with_tags": round(100 * sum(bool(r["tags"]) for r in recs) / len(recs), 1),
            "mean_lowercase_frac": round(st.mean(r["lower_frac"] for r in recs), 3),
            "top_tags": tagc.most_common(6),
        }
        s = summary["dna"][label]
        print(f"  {label:<12} n={s['sequences']:<4} len~{s['median_len_bp']:>7,.0f}  GC={s['median_gc_pct']:>4}%  "
              f"CpG o/e={s['median_cpg_oe']:.2f}  6-multiple={s['pct_len_multiple_of_6']}%  tagged={s['pct_with_tags']}%")
    first = sample(DNA_CONFIGS["Mammal"], args.n_dna, 1000)[0]["text"]
    seq = DNA_RE.search(first).group(2)
    summary["tokenization_example"] = {"raw_prefix": first[:80], "first_6mer_tokens": kmer6_tokens(seq[:60])}
    plot_dna(groups)

    print("[2/3] SFT")
    sft_rows = sample("sft_mix", args.n_sft)
    dom, med_turns, examples = sft_stats(sft_rows)
    summary["sft"] = {"rows": len(sft_rows), "domains": dom.most_common(), "median_messages": med_turns}
    with (RES / "sft_examples.md").open("w") as f:
        for d, msgs in examples.items():
            f.write(f"## {d}\n\n")
            for m in msgs[:4]:
                body = m["content"] if len(m["content"]) < 1500 else m["content"][:1500] + " …[truncated]"
                f.write(f"**{m['role']}**\n\n```\n{body}\n```\n\n")
    for d, c in dom.most_common():
        print(f"  {d:<30} {c:>5}  median messages={med_turns[d]}")

    print("[3/3] RL + independent answer check")
    rl_rows = sample("rl_mix", args.n_rl)
    rl_groups = collections.Counter(r["task_group"] for r in rl_rows)
    per_type, details = rl_verify(rl_rows)
    tot = collections.Counter()
    for c in per_type.values():
        tot.update(c)
    summary["rl"] = {"rows": len(rl_rows), "task_groups": rl_groups.most_common(),
                     "basic_dna_check": {"totals": dict(tot),
                                         "per_task_type": {k: dict(v) for k, v in sorted(per_type.items())
                                                           if v["agree"] + v["disagree"]}}}
    with (RES / "rl_check_details.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["task_type", "status", "gold", "ours", "prompt"])
        w.writeheader(); w.writerows(details)
    checked = tot["agree"] + tot["disagree"]
    print(f"  basic_dna rows: {sum(tot.values())}; checked by our code: {checked}; "
          f"agree: {tot['agree']} ({100 * tot['agree'] / max(checked, 1):.1f}%); "
          f"disagree: {tot['disagree']}; unparsed: {tot['unparsed']}; no checker: {tot['no_checker']}")
    plot_mix(dom, rl_groups)

    (RES / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Done. Figures in {FIG}, results in {RES}")


if __name__ == "__main__":
    main()
    sys.stdout.flush()
    os._exit(0)  # datasets' streaming threads can keep the interpreter alive after main()
