"""
run_audit.py -- reproduce the Boltz-2 triage calibration against real wet-lab labels,
from the Anthropic claude-protein-binder-design release (CC BY 4.0).

Self-contained: fetches only the design_summary parquet (1,440 x 65) through the agent
proxy, no dependency on the binder-benchmark-audit skill path.

    python3 binder_audit/run_audit.py

Answers, for the co-folding interface score this project uses (Boltz-2 ipSAE_min):
how much does triaging by it enrich the binder hit rate, and does the score track
affinity? (No: it predicts WHETHER a design binds, not how tightly.)
"""
import os, json, urllib.request
import numpy as np, pandas as pd

HF = ("https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/"
      "resolve/main/data/tables/design_summary.parquet")
EXCLUDED = "Mature GDF-8"                       # aggregated in both assay formats -> excluded
SCORE3 = ["ipsae_min_ef2full", "ipsae_min_ef2fast", "ipsae_min_ptxv2"]   # campaign selectors


def load(cache="binder_audit/design_summary.parquet"):
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if not os.path.exists(cache):
        urllib.request.urlretrieve(HF, cache)
    return pd.read_parquet(cache)


def within_target_map(d, score_col):
    from sklearn.metrics import average_precision_score
    aps, beats = [], 0
    for _, g in d.groupby("target"):
        if g.binder.sum() < 3 or (~g.binder).sum() < 3:
            continue
        s = (g[score_col] - g[score_col].mean()) / g[score_col].std(ddof=0)
        ap = average_precision_score(g.binder, s)
        aps.append(ap); beats += ap > g.binder.mean()
    return float(np.mean(aps)), f"{beats}/{len(aps)}"


def main():
    raw = load()
    d = raw[raw.target != EXCLUDED].copy()
    d["binder"] = d.binder_final.astype("boolean").fillna(False).astype(bool)
    d["score3"] = d[SCORE3].mean(axis=1)
    d["boltz2"] = d["ipsae_min_boltz2"]

    m3, b3 = within_target_map(d, "score3")
    mb, bb = within_target_map(d, "boltz2")

    enrich = []
    for frac in (1.0, 0.5, 0.3, 0.1, 0.05):
        keep = d.groupby("target", group_keys=False).apply(
            lambda g: g.nlargest(max(1, int(round(len(g) * frac))), "boltz2"),
            include_groups=False)
        enrich.append({"keep_pct": 100 * frac, "n": int(len(keep)),
                       "hit_rate_pct": round(100 * keep.binder.mean(), 1)})

    bnd = d[d.binder & d.kd_nM_final.notna()].copy()
    bnd["logkd"] = np.log10(bnd.kd_nM_final)
    rho = float(bnd[["boltz2", "logkd"]].corr(method="spearman").iloc[0, 1])

    out = {"analysis_set": {"n": int(len(d)), "targets": int(d.target.nunique()),
                            "overall_hit_rate_pct": round(100 * d.binder.mean(), 1)},
           "campaign_score3_within_target_mAP": round(m3, 3),
           "campaign_score3_beats_chance_targets": b3,
           "boltz2_within_target_mAP": round(mb, 3),
           "boltz2_beats_chance_targets": bb,
           "boltz2_triage_enrichment": enrich,
           "boltz2_vs_logKD_spearman_amongBinders": round(rho, 2)}
    json.dump(out, open("binder_audit/audit_boltz2_calibration.json", "w"), indent=2)
    print(json.dumps(out, indent=2))
    print("\noverall %.1f%% (paper 26.8); score3 mAP %.3f (paper 0.51); "
          "Boltz-2 mAP %.3f beats-chance %s; Boltz-2 vs logKD rho %.2f (feasibility, not affinity)"
          % (100 * d.binder.mean(), m3, mb, bb, rho))


if __name__ == "__main__":
    main()
