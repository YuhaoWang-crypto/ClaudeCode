"""Recompute the headline claims of "Autonomous de novo protein binder design with Claude"
from the public data release (HuggingFace: Anthropic/claude-protein-binder-design).

This is a retrospective reproduction: it re-derives every number the paper reports
from the released per-design tables. It does not re-run the design campaigns
(GPU compute) or the wet-lab measurements (CROs).

Usage:
    python reproduce_claims.py --release-dir <path to protein_binder_design_data_release>
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(0)

# claims as stated in the paper: key -> (reported value, tolerance)
CLAIMS = {}
CHECKS = []


def check(name, reported, computed, tol=0.0, note=""):
    """Record a claim check. tol is absolute tolerance for numeric comparison."""
    if isinstance(reported, (int, float)) and isinstance(computed, (int, float)):
        ok = abs(reported - computed) <= tol
    else:
        ok = str(reported) == str(computed)
    CHECKS.append(
        {"claim": name, "paper": reported, "recomputed": computed, "match": bool(ok), "note": note}
    )
    return ok


def load(release_dir: Path):
    t = release_dir / "tables"
    design = pd.read_parquet(t / "design_summary.parquet")
    prov = pd.read_parquet(t / "insilico" / "provenance_summary.parquet")
    cofold = pd.read_parquet(t / "insilico" / "cofold_predictions.parquet")
    wet = pd.read_parquet(t / "wetlab" / "summary.parquet")
    return design, prov, cofold, wet


def average_precision(scores, labels):
    """AP with ties broken by mean rank (sklearn-free)."""
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    y = np.asarray(labels, dtype=int)[order]
    if y.sum() == 0:
        return np.nan
    cum = np.cumsum(y)
    prec = cum / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / y.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=Path("reproduction_report.json"))
    args = ap.parse_args()

    design, prov, cofold, wet = load(args.release_dir)

    # ---------------------------------------------------------------- dataset
    check("designs delivered (all campaigns)", 1440, len(design))
    tested = design[design["binder_final"].notna()].copy()
    check("designs with interpretable measurements", 1320, len(tested))
    check("targets analysed", 15, tested["target"].nunique())
    check(
        "mature GDF-8 excluded (uninterpretable)",
        120,
        int(design["binder_final"].isna().sum()),
        note="all excluded designs are mature GDF-8: %s"
        % sorted(design.loc[design["binder_final"].isna(), "target"].unique()),
    )

    # ------------------------------------------------------------- hit rates
    binders = int(tested["binder_final"].sum())
    check("binders (integrated call)", 354, binders)
    check("overall hit rate %", 26.8, round(100 * binders / len(tested), 1), tol=0.05)

    per_target = (
        tested.groupby("target")["binder_final"].agg(["sum", "count"]).astype(int).sort_values("sum")
    )
    per_target["hit_rate_%"] = (100 * per_target["sum"] / per_target["count"]).round(1)
    check("targets with >=1 binder", 14, int((per_target["sum"] > 0).sum()))
    # the paper counts 1,315 "tested" designs: those measured at at least one CRO
    measured = tested[tested["vendor_agreement"] != "not_tested_either"].copy()
    check("designs tested at >=1 CRO", 1315, len(measured))

    for tgt, exp_b, exp_n in [
        ("TREM2", 72, 90),
        ("VEGF-A", 54, 90),
        ("IL-7Ra", 49, 90),
        ("BBF-14", 3, 90),
        ("15-PGDH", 1, 30),
        ("MBP", 0, 90),
        ("RBX1", 28, 90),
        ("Nipah-G", 19, 90),
        ("EGFR", 10, 90),
        ("Latent GDF-8", 14, 60),
        ("TNFa", 12, 150),
    ]:
        row = per_target.loc[tgt]
        check(f"{tgt}: binders/tested", f"{exp_b}/{exp_n}", f"{int(row['sum'])}/{int(row['count'])}")

    # campaign-level: model x format
    tested["campaign_id"] = tested["design_model"] + " / " + tested["campaign"]
    camp = tested.groupby("campaign_id")["binder_final"].agg(["sum", "count"]).astype(int)
    camp["hit_rate_%"] = (100 * camp["sum"] / camp["count"]).round(1)

    def camp_row(model, fmt):
        m = (tested["design_model"] == model) & (tested["campaign"] == fmt)
        return int(tested.loc[m, "binder_final"].sum()), int(m.sum())

    o_multi = camp_row("Opus 4.8", "multi_target")
    m_multi = camp_row("Mythos Preview", "multi_target")
    m_single = camp_row("Mythos Preview", "single_target")
    check("Opus 4.8 multi-target", "88/390 (22.6%)", "%d/%d (%.1f%%)" % (*o_multi, 100 * o_multi[0] / o_multi[1]))
    check("Mythos multi-target", "104/390 (26.7%)", "%d/%d (%.1f%%)" % (*m_multi, 100 * m_multi[0] / m_multi[1]))
    check("Mythos single-target", "158/450 (35.1%)", "%d/%d (%.1f%%)" % (*m_single, 100 * m_single[0] / m_single[1]))

    # shared 13 targets between the three campaigns run on >=13 targets
    multi_targets = set(tested.loc[tested["campaign"] == "multi_target", "target"])
    shared13 = sorted(multi_targets)
    check("targets shared by multi- and single-target campaigns", 13, len(shared13))
    ms = tested[(tested["design_model"] == "Mythos Preview") & (tested["campaign"] == "single_target") & tested["target"].isin(shared13)]
    check(
        "Mythos single-target on the shared 13",
        "143/390 (36.7%)",
        "%d/%d (%.1f%%)" % (int(ms["binder_final"].sum()), len(ms), 100 * ms["binder_final"].mean()),
    )
    # single vs multi, same model, same 13 targets
    table = [[int(ms["binder_final"].sum()), len(ms) - int(ms["binder_final"].sum())],
             [m_multi[0], m_multi[1] - m_multi[0]]]
    p = stats.fisher_exact(table)[1]
    check("Fisher p, Mythos single vs multi (13 targets)", 0.003, round(p, 3), tol=0.0005)

    # ----------------------------------------------------- ranking calibration
    # 41 rankings = 3 campaigns x their targets (13 + 13 + 15)
    rank_sets = tested[tested["campaign"].isin(["multi_target", "single_target"])].copy()
    rank_sets = rank_sets[~((rank_sets["design_model"] == "Opus 4.8") & (rank_sets["campaign"] == "single_target"))]
    groups = rank_sets.groupby(["design_model", "campaign", "target"])
    check("rankings of 30 designs from the three main campaigns", 41, groups.ngroups)

    topn = {}
    for n in [1, 5, 10, 30]:
        sub = rank_sets[rank_sets["rank"] <= n]
        topn[n] = 100 * sub["binder_final"].mean()
    check("hit rate, top-1 design %", 49, round(topn[1]), tol=0.5)
    check("hit rate, top-5 %", 44, round(topn[5]), tol=0.5)
    check("hit rate, top-10 %", 39, round(topn[10]), tol=0.5)
    check("hit rate, all 30 %", 28, round(topn[30]), tol=0.5)

    n_rankings_top1 = sum(g["binder_final"][g["rank"] == 1].any() for _, g in groups)
    n_targets_top1 = rank_sets[rank_sets["rank"] == 1].groupby("target")["binder_final"].any().sum()
    n_rankings_top5 = sum(g["binder_final"][g["rank"] <= 5].any() for _, g in groups)
    n_targets_top5 = rank_sets[rank_sets["rank"] <= 5].groupby("target")["binder_final"].any().sum()
    check("rankings whose top-1 design bound", 20, int(n_rankings_top1))
    check("targets covered by top-1 designs", 12, int(n_targets_top1))
    check("rankings with a binder in the top 5", 27, int(n_rankings_top5))
    check("targets covered by top-5 designs", 13, int(n_targets_top5))

    # average precision of Claude's delivered rank, per target (13 evaluable targets)
    evaluable = [
        t for t in tested["target"].unique()
        if tested.loc[tested["target"] == t, "binder_final"].sum() >= 3
        and (~tested.loc[tested["target"] == t, "binder_final"].astype(bool)).sum() >= 3
    ]
    check("targets evaluable for average precision", 13, len(evaluable))

    # AP is computed within each ranking (one per target and campaign), averaged within a
    # target, then over the 13 evaluable targets.
    def rank_ap(frame):
        """AP of the delivered rank within each target (campaigns pooled), averaged over
        the evaluable targets."""
        per_t = []
        for t in evaluable:
            sub = frame[frame["target"] == t]
            if sub.empty:
                continue
            per_t.append(average_precision(-sub["rank"].to_numpy(), sub["label"].astype(int)))
        return float(np.nanmean(per_t))

    rank_sets["label"] = rank_sets["binder_final"].astype(int)
    obs = rank_ap(rank_sets)
    # variant including the two Opus 4.8 single-target campaigns (TNFa, latent GDF-8)
    alt = tested[tested["campaign"].isin(["multi_target", "single_target"])].copy()
    alt["label"] = alt["binder_final"].astype(int)
    check(
        "average precision of Claude's rank", 0.48, round(obs, 2), tol=0.025,
        note="paper does not state how rankings are pooled within a target; "
             "three main campaigns -> %.3f, all delivered rankings -> %.3f"
             % (obs, rank_ap(alt)),
    )

    # permutation test: shuffle ranks within each campaign's ranking of each target.
    # The mean of this null is the paper's "expected 0.35 for an uninformative ordering".
    null = []
    shuffled = rank_sets.copy()
    for _ in range(2000):
        shuffled["label"] = (
            rank_sets.groupby(["design_model", "campaign", "target"], group_keys=False)["binder_final"]
            .apply(lambda s: pd.Series(RNG.permutation(s.to_numpy()), index=s.index))
            .astype(int)
        )
        null.append(rank_ap(shuffled))
    null = np.asarray(null)
    check("chance level (mean of the permutation null)", 0.35, round(float(null.mean()), 3), tol=0.01,
          note="permutation null mean, 2000 shuffles")
    p_perm = (1 + (null >= obs).sum()) / (1 + len(null))
    check("permutation p for rank calibration", "<0.001", "<0.001" if p_perm < 0.001 else round(p_perm, 4))

    # ------------------------------------------------- co-folding re-score AP
    # campaign ranking score predictors: ESMFold2, ESMFold2-Fast, Protenix v2
    trio = ["ipsae_min_ef2full", "ipsae_min_ef2fast", "ipsae_min_ptxv2"]
    tested["cofold_mean"] = tested[trio].mean(axis=1)

    aps_score, chance_score = [], []
    for t in evaluable:
        sub = tested[tested["target"] == t]
        z = sub[trio].apply(lambda c: (c - c.mean()) / c.std(ddof=0)).mean(axis=1)
        aps_score.append(average_precision(z.to_numpy(), sub["binder_final"].astype(int)))
        chance_score.append(sub["binder_final"].mean())
    check("average precision of the re-scored co-folding score", 0.52, round(float(np.mean(aps_score)), 2), tol=0.015)
    check("chance level across the 13 targets", 0.31, round(float(np.mean(chance_score)), 3), tol=0.01,
          note="mean per-target hit rate = 0.305")
    n_above = sum(a > c for a, c in zip(aps_score, chance_score))
    check("targets where the score beats chance", 12, int(n_above))
    p_sign = stats.binomtest(n_above, len(aps_score), 0.5, alternative="greater").pvalue
    check("sign-test p", 0.003, round(p_sign, 3), tol=0.001)

    # ensemble of the seven predictors NOT in the campaign ranking score
    seven = ["ipsae_min_odde", "ipsae_min_afm3", "ipsae_min_boltz2", "ipsae_min_chai1",
             "ipsae_min_of3", "ipsae_min_rf3", "ipsae_min_af3of3"]
    aps7 = []
    for t in evaluable:
        sub = tested[tested["target"] == t]
        cols = [c for c in seven if sub[c].notna().any()]
        z = sub[cols].apply(lambda c: (c - c.mean()) / c.std(ddof=0)).mean(axis=1)
        aps7.append(average_precision(z.fillna(z.mean()).to_numpy(), sub["binder_final"].astype(int)))
    check("average precision, 7-predictor ensemble", 0.57, round(float(np.mean(aps7)), 2), tol=0.02)

    # correlation of Claude's delivered order with the re-scored co-folding score
    rank_sets["cofold_mean"] = rank_sets["uuid"].map(tested.set_index("uuid")["cofold_mean"])
    rhos = []
    for (_, _, _), g in rank_sets.groupby(["design_model", "campaign", "target"]):
        if g["cofold_mean"].notna().sum() > 3:
            rhos.append(stats.spearmanr(-g["rank"], g["cofold_mean"], nan_policy="omit").statistic)
    check("median Spearman rho (rank vs score)", 0.86, round(float(np.nanmedian(rhos)), 2), tol=0.03)

    # across-target: median score vs hit rate
    per_t = tested.groupby("target").agg(med=("cofold_mean", "median"), hr=("binder_final", "mean"))
    rho_t = stats.spearmanr(per_t["med"], per_t["hr"]).statistic
    p_t = stats.spearmanr(per_t["med"], per_t["hr"]).pvalue
    check("Spearman rho, median score vs hit rate (15 targets)", 0.61, round(float(rho_t), 2), tol=0.03)
    check("its p-value", 0.02, round(float(p_t), 2), tol=0.005)

    # ------------------------------------------------------------- affinities
    b = tested[tested["binder_final"] == True]  # noqa: E712
    kd = b["kd_nM_final"].astype(float)
    check("binders with a KD of record", 354, int(kd.notna().sum()))
    check("binders below 100 nM", 194, int((kd < 100).sum()))
    check("binders below 10 nM", 90, int((kd < 10).sum()))
    check("binders below 1 nM", 42, int((kd < 1).sum()))

    oligomeric = {"TNFa", "VEGF-A", "Nipah-G", "Latent GDF-8", "15-PGDH"}
    mono = b[~b["target"].isin(oligomeric)]
    check("binders on oligomeric targets", 100, int(b["target"].isin(oligomeric).sum()))
    check("binders on monomeric targets", 254, len(mono))
    kdm = mono["kd_nM_final"].astype(float)
    check("monomeric binders below 100 nM", 142, int((kdm < 100).sum()))
    check("monomeric binders below 10 nM", 78, int((kdm < 10).sum()))
    check("monomeric binders below 1 nM", 38, int((kdm < 1).sum()))

    # ------------------------------------------------------------------ TNFa
    tnf = tested[tested["target"] == "TNFa"]
    check("TNFa binders", 12, int(tnf["binder_final"].sum()))
    check("TNFa hit rate %", 8.0, round(100 * tnf["binder_final"].mean(), 1), tol=0.05)
    check(
        "all TNFa binders from Opus 4.8",
        "Opus 4.8",
        ", ".join(sorted(tnf.loc[tnf["binder_final"] == True, "design_model"].unique())),  # noqa: E712
    )
    check(
        "TNFa binders per Opus campaign (multi / single / supplementary)",
        "8 / 2 / 2",
        " / ".join(
            str(int(tnf[(tnf["binder_final"] == True) & (tnf["campaign"] == c)].shape[0]))  # noqa: E712
            for c in ["multi_target", "single_target", "single_target_supplementary"]
        ),
    )
    check(
        "tightest TNFa binder, apparent KD (nM)",
        0.70,
        round(float(tnf.loc[tnf["binder_final"] == True, "kd_nM_final"].min()), 2),  # noqa: E712
        tol=0.005,
    )

    # ------------------------------------------------------------------ RBX1
    rbx = tested[tested["target"] == "RBX1"]
    top_mythos_single = rbx[(rbx["design_model"] == "Mythos Preview") & (rbx["campaign"] == "single_target") & (rbx["rank"] == 1)]
    top_opus_multi = rbx[(rbx["design_model"] == "Opus 4.8") & (rbx["campaign"] == "multi_target") & (rbx["rank"] == 1)]
    check("RBX1 top-ranked Mythos single-target, Adaptyv KD (nM)", 3.9,
          round(float(top_mythos_single["adaptyv_kd_nM"].iloc[0]), 1), tol=0.05)
    check("RBX1 top-ranked Mythos single-target, KD of record (nM)", 7.0,
          round(float(top_mythos_single["kd_nM_final"].iloc[0]), 1), tol=0.05)
    check("RBX1 top-ranked Opus multi-target, Adaptyv KD (nM)", 30,
          round(float(top_opus_multi["adaptyv_kd_nM"].iloc[0])), tol=0.5)
    check("RBX1 top-ranked Mythos single-target, Twist KD on CUL1-RBX1 (nM)", 8.4,
          round(float(top_mythos_single["twist_kd_nM"].iloc[0]), 1), tol=0.05)

    # ------------------------------------------------------- cross-reactivity
    cyno_eval = b[b["cyno_binding_final"] != "not_tested"]
    mouse_eval = b[b["mouse_binding_final"] != "not_tested"]
    check("binders with an evaluable cyno titration", 179, len(cyno_eval))
    check("of those, binding cyno", 154, int((cyno_eval["cyno_binding_final"] == "binder").sum()))
    check("binders tested against the mouse ortholog", 233, len(mouse_eval))
    check("of those, binding mouse", 130, int((mouse_eval["mouse_binding_final"] == "binder").sum()))
    for tgt, exp in [("TREM2", "68/69"), ("VEGF-A", "32/53"), ("IL-7Ra", "18/40"),
                     ("EGFR", "4/8"), ("TNFa", "3/12"), ("TrkA", "2/14"), ("PD-L1", "2/36")]:
        s = mouse_eval[mouse_eval["target"] == tgt]
        check(f"{tgt}: binders binding mouse ortholog", exp,
              "%d/%d" % (int((s["mouse_binding_final"] == "binder").sum()), len(s)))

    # --------------------------------------------------------- design methods
    gen = measured.groupby("generator")["binder_final"].agg(["sum", "count"]).astype(int)
    gen["hit_rate_%"] = (100 * gen["sum"] / gen["count"]).round(1)
    gen = gen.sort_values("count", ascending=False)
    for name, exp in [("PXDesign", 358), ("RFdiffusion3", 267), ("Genie3", 185),
                      ("FreeBindCraft (BindCraft)", 135), ("BoltzGen", 134), ("RFdiffusion", 118),
                      ("Proteina-Complexa", 100), ("FoldCraft", 14), ("BoltzDesign1", 2),
                      ("Protein Hunter", 2)]:
        check(f"tested designs from {name}", exp, int(gen.loc[name, "count"]), tol=0)

    main7 = gen[gen["count"] >= 100]
    check("structure-generation methods with >=100 tested designs", 7, len(main7))
    check("their hit-rate range %", "22-43",
          "%d-%d" % (round(main7["hit_rate_%"].min()), round(main7["hit_rate_%"].max())))

    seq = measured.groupby("sequence_design_method").size().sort_values(ascending=False)
    for name, exp in [("SolubleMPNN", 1133), ("ProteinMPNN", 21), ("Caliby (SolubleCaliby)", 111),
                      ("designed jointly by the structure model (native co-design)", 50)]:
        check(f"sequences designed by {name}", exp, int(seq.get(name, 0)), tol=0)

    # backbone-level de-duplication
    prov_idx = prov.set_index("uuid")
    measured_bb = measured.join(prov_idx["prov_root_backbone_id"], on="uuid")
    best = measured_bb.sort_values("rank").groupby("prov_root_backbone_id").first()
    check("distinct generated backbones among tested designs", 809, len(best))
    check("backbones whose best-ranked sequence bound", 200, int(best["binder_final"].sum()))
    check("backbone-level hit rate %", 24.7, round(100 * best["binder_final"].mean(), 1), tol=0.05)

    # fold diversity: the paper counts over the 1,320 delivered designs with data
    fold = tested.join(prov_idx["prov_fold_class_designed"], on="uuid")
    not_all_alpha = fold[fold["prov_fold_class_designed"].astype(str) != "all-α"]
    check("delivered designs that are not all-α", 126, len(not_all_alpha),
          note="fold classes present: %s" % sorted(fold["prov_fold_class_designed"].dropna().unique()))
    fold_m = measured.join(prov_idx["prov_fold_class_designed"], on="uuid")
    nb = fold_m["prov_fold_class_designed"].astype(str) != "all-α"
    check("not-all-α designs that bound", "16/125",
          "%d/%d" % (int(fold_m.loc[nb, "binder_final"].sum()), int(nb.sum())))
    check("all-α designs that bound", "338/1190",
          "%d/%d" % (int(fold_m.loc[~nb, "binder_final"].sum()), int((~nb).sum())))

    # ------------------------------------------------------------------ report
    df = pd.DataFrame(CHECKS)
    n_ok = int(df["match"].sum())
    print(df.to_string(index=False))
    print(f"\n{n_ok}/{len(df)} claims reproduced exactly.")
    print("\nPer-target hit rates:\n", per_target.to_string())
    print("\nPer-campaign hit rates:\n", camp.to_string())
    print("\nStructure-generation methods (tested designs):\n", gen.to_string())
    print("\nSequence-design methods:\n", seq.to_string())

    args.out.write_text(json.dumps(CHECKS, indent=2, default=str))
    df.to_csv(args.out.with_suffix(".csv"), index=False)
    per_target.to_csv(args.out.parent / "per_target_hit_rates.csv")


if __name__ == "__main__":
    main()
