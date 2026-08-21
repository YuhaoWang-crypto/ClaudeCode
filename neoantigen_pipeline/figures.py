"""Figures for a pipeline run. Four panels, one per question a reviewer asks.

    from neoantigen_pipeline import figures
    figures.summary_figure(res, cfg, "demo_out/summary.png", bench=scored)
"""

from __future__ import annotations

import os
from typing import Dict, Optional

import numpy as np
import pandas as pd


def summary_figure(res: Dict[str, object], cfg, path: str,
                   bench: Optional[pd.DataFrame] = None) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Neoantigen selection - {cfg.patient.patient_id}", fontsize=14)

    # (a) how many candidates survive each gate
    wf = res.get("gate_waterfall")
    a = ax[0, 0]
    if isinstance(wf, pd.DataFrame) and not wf.empty:
        labels = [s.replace("gate_", "") for s in wf["step"]]
        a.barh(range(len(wf))[::-1], wf["n"], color="#4C72B0")
        a.set_yticks(range(len(wf))[::-1])
        a.set_yticklabels(labels, fontsize=9)
        for i, n in enumerate(wf["n"]):
            a.text(n, len(wf) - 1 - i, f" {n}", va="center", fontsize=8)
        a.set_xlabel("variants remaining")
    a.set_title("(a) variant gates: what survives, and where it dies", fontsize=11)

    # (b) the score landscape: every candidate, the selected ones marked
    b = ax[0, 1]
    scored = res.get("scored")
    sel = res.get("selected")
    if isinstance(scored, pd.DataFrame) and not scored.empty:
        x = scored["mut_rank"].clip(lower=1e-3)
        b.scatter(x, scored["neo_score"], s=12, alpha=0.35, color="#999999",
                  label=f"candidates (n={len(scored)})")
        if isinstance(sel, pd.DataFrame) and not sel.empty:
            b.scatter(sel["mut_rank"].clip(lower=1e-3), sel["neo_score"], s=42,
                      color="#C44E52", edgecolor="k", linewidth=0.4,
                      label=f"selected (n={len(sel)})", zorder=3)
        b.set_xscale("log")
        b.set_xlabel("NetMHCpan-4.1 EL %rank (lower = better presented)")
        b.set_ylabel("composite neoantigen score")
        b.legend(fontsize=8, loc="lower left")
    b.set_title("(b) presentation alone does not decide the payload", fontsize=11)

    # (c) per-feature contribution of the selected set
    c = ax[1, 0]
    if isinstance(sel, pd.DataFrame) and not sel.empty:
        feats = [k for k in cfg.weight_dict() if f"feat_{k}" in sel.columns]
        w = cfg.weight_dict()
        vals = [(sel[f"feat_{k}"].fillna(0) * w[k]).mean() for k in feats]
        base = [(res["scored"][f"feat_{k}"].fillna(0) * w[k]).mean() for k in feats]
        y = np.arange(len(feats))
        c.barh(y + 0.2, vals, height=0.4, color="#C44E52", label="selected")
        c.barh(y - 0.2, base, height=0.4, color="#999999", label="all candidates")
        c.set_yticks(y)
        c.set_yticklabels(feats, fontsize=9)
        c.set_xlabel("mean weighted contribution to the score")
        c.legend(fontsize=8)
    c.set_title("(c) which features the selected set actually wins on", fontsize=11)

    # (d) benchmark, or the junction-cost result when no benchmark was run
    d = ax[1, 1]
    if bench is not None and not bench.empty and "label" in bench:
        from .benchmark import auc
        for col, colour, name in (("score_netmhcpan_only", "#999999", "NetMHCpan %rank alone"),
                                  ("score_composite_no_tcr", "#4C72B0", "composite (no TCR prior)"),
                                  ("score_composite", "#C44E52", "composite (all features)")):
            if col not in bench:
                continue
            s = bench[col].fillna(0).to_numpy()
            y = bench["label"].to_numpy()
            thr = np.unique(np.concatenate([[-np.inf], np.sort(s), [np.inf]]))
            tpr = [(s[y == 1] >= t).mean() for t in thr]
            fpr = [(s[y == 0] >= t).mean() for t in thr]
            d.plot(fpr, tpr, color=colour, lw=1.8,
                   label=f"{name} (AUC>={auc(y, s):.3f})")
        d.plot([0, 1], [0, 1], "k--", lw=0.8)
        d.set_xlabel("false positive rate (decoys are UNLABELLED, not verified negative)")
        d.set_ylabel("true positive rate")
        d.legend(fontsize=8, loc="lower right")
        d.set_title("(d) validated neoepitopes vs matched decoys", fontsize=11)
    else:
        con = res.get("construct")
        if con:
            d.bar(["input order", "junction-optimized"],
                  [con["junction_cost_naive"], con["junction_cost_optimized"]],
                  color=["#999999", "#C44E52"])
            d.set_ylabel("worst predicted junction presentation (sum over junctions)")
        d.set_title("(d) minigene ordering removes junction epitopes", fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _roc(ax, bench: pd.DataFrame, title: str):
    from .benchmark import auc
    for col, colour, name in (("score_netmhcpan_only", "#999999", "NetMHCpan %rank alone"),
                              ("score_composite_no_tcr", "#4C72B0", "composite (no TCR prior)"),
                              ("score_composite", "#C44E52", "composite (all features)")):
        if col not in bench:
            continue
        s = bench[col].fillna(0).to_numpy()
        y = bench["label"].to_numpy()
        thr = np.unique(np.concatenate([[-np.inf], np.sort(s), [np.inf]]))
        tpr = [(s[y == 1] >= t).mean() for t in thr]
        fpr = [(s[y == 0] >= t).mean() for t in thr]
        ax.plot(fpr, tpr, color=colour, lw=1.8, label=f"{name} (AUC>={auc(y, s):.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(title, fontsize=11)


def benchmark_figure(scored_a: pd.DataFrame, scored_b: pd.DataFrame,
                     path: str) -> str:
    """Side by side: the easy benchmark and the one that means something."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    _roc(ax[0], scored_a, "(a) decoys matched on allele + length only\n"
                          "(mostly a binder / non-binder test)")
    _roc(ax[1], scored_b, "(b) decoys ALSO matched on predicted binding\n"
                          "(the immunogenicity question)")
    for i, (d, lab, colour) in enumerate(((scored_a, "allele-matched", "#999999"),
                                          (scored_b, "binding-matched", "#C44E52"))):
        for lbl, ls in ((1, "-"), (0, "--")):
            v = np.log10(d[d["label"] == lbl]["mut_rank"].clip(lower=1e-3).dropna())
            ax[2].hist(v, bins=25, histtype="step", ls=ls, color=colour, lw=1.5,
                       label=f"{lab}, {'positives' if lbl else 'decoys'}")
    ax[2].set_xlabel("log10 predicted %rank")
    ax[2].set_ylabel("count")
    ax[2].legend(fontsize=7)
    ax[2].set_title("(c) why (a) is easy: the decoy\nbinding distribution moves", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def junction_figure(res: Dict[str, object], path: str) -> Optional[str]:
    """Before/after junction cost plus the final junction %rank distribution."""
    con = res.get("construct")
    if not con:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    jn = con.get("junction_scan")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].bar(["input order", "optimized"],
              [con["junction_cost_naive"], con["junction_cost_optimized"]],
              color=["#999999", "#C44E52"])
    ax[0].set_ylabel("summed junction presentation cost")
    ax[0].set_title("concatemer ordering")
    if isinstance(jn, pd.DataFrame) and not jn.empty:
        ax[1].hist(np.log10(jn["percentile_rank"].clip(lower=1e-3)), bins=30,
                   color="#4C72B0")
        ax[1].axvline(np.log10(0.5), color="#C44E52", ls="--",
                      label="strong-binder threshold (0.5%)")
        ax[1].set_xlabel("log10 %rank of junction peptides")
        ax[1].set_ylabel("count")
        ax[1].legend(fontsize=8)
        ax[1].set_title(f"final junctions: {int(jn['flagged'].sum())} flagged "
                        f"of {len(jn)} peptides")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
