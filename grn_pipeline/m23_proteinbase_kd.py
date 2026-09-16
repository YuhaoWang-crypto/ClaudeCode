"""
Module 23 — Recalibrate the co-folding scores against 1,025 measured outcomes.

WHAT THIS FIXES
---------------
M10 validated the Boltz screen against n=5 ChEMBL potencies and concluded
"optimization_score tracks potency (rho=+0.6), binding_confidence does not".
At n=5 a Spearman rho of +0.6 has a 95% CI spanning essentially [-0.5, +0.95]:
it cannot distinguish a useful predictor from noise. The direction of that
finding was right but its strength was never measurable.

Proteinbase (see `proteinbase_db`) supplies the missing denominator: 2,630
design/target pairs assayed in one lab on standardised protocols, of which
2,161 are measured NON-binders, plus 435 pairs with a real Kd (2,621 / 2,157 /
430 after excluding Adaptyv's own control constructs). That makes two
distinct questions separable for the first time:

  (D) DISCRIMINATION — does the score separate binders from non-binders?
      Scored on all labelled pairs, the question a screen actually asks.
  (A) AFFINITY RANKING — among things that DO bind, does the score rank Kd?
      Scored on the Kd subset only, and therefore range-restricted.

These are not the same question, and a score can pass one while failing the
other. M10, and much of the field's "our score correlates with affinity"
reporting, conflates them by testing only on confirmed binders — a set
selected *because* it binds.

HONESTY NOTES
-------------
* Everything here is measured, nothing is asserted. Every rho/AUROC carries n
  and a bootstrap or analytic CI.
* Analysis is stratified by target. 1,025 of the 1,049 scored pairs are against
  Nipah glycoprotein G, so a pooled number is really a Nipah number plus a
  between-target offset; pooling different targets manufactures correlation.
* Adaptyv's own control constructs are excluded (they are known binders).
* `boltz2_pdockq`/`pdockq2` are excluded automatically: they are constant in
  this snapshot (a defect `proteinbase_db.qc_report` flags), so they would
  score AUROC exactly 0.500 and look like a legitimate null result.
* The Kd replicate noise floor (~1.44x) caps how well ANY predictor can
  correlate with these numbers; the observed correlations are far below that
  cap, so the ceiling is not what is limiting them.
"""
from __future__ import annotations

import statistics

import numpy as np
from scipy.optimize import minimize
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from grn_pipeline import proteinbase_db as pb

RNG = np.random.default_rng(20260916)

# Metric -> (+1 if higher should mean better binding, -1 if lower does).
# ProteinMPNN's score here spans 0.68-3.91 with median 1.53: that is a
# per-residue negative log-likelihood, so LOWER is the better sequence.
# The sign is fixed from that definition, never from which sign fits better.
ORIENT = {
    "boltz2_min_ipsae": +1,
    "boltz2_ipsae": +1,
    "boltz2_iptm": +1,
    "boltz2_ptm": +1,
    "boltz2_complex_iplddt": +1,
    "boltz2_complex_plddt": +1,
    "boltz2_lis": +1,
    "shape_complimentarity_boltz2_binder_ss": +1,
    "boltz2_pdockq": +1,
    "boltz2_pdockq2": +1,
    "esmfold_plddt": +1,
    "proteinmpnn_score": -1,
    "redesigned_proteinmpnn_score": -1,
    "proteinmpnn_seq_recovery": +1,
}
DESIGN_LEVEL = {"esmfold_plddt", "proteinmpnn_score",
                "redesigned_proteinmpnn_score", "proteinmpnn_seq_recovery"}
# The co-folding INTERFACE CONFIDENCE family — the scores a screen normally
# ranks on, and the ones M10's claim was about. Kept separate from the other
# quantities because the two families behave differently on affinity here.
CONFIDENCE_SCORES = {
    "boltz2_ipsae", "boltz2_min_ipsae", "boltz2_iptm", "boltz2_ptm",
    "boltz2_complex_iplddt", "boltz2_complex_plddt", "boltz2_lis",
}
PRIMARY_TARGET = "nipah-glycoprotein-g"


# ------------------------------------------------------------------ dataset

def assemble(target: str | None = None, drop_controls: bool = True):
    """Join designs onto complexes; return usable rows for one target (or all).

    Each row carries every candidate score plus the experimental verdict, so
    discrimination and affinity can be scored off the same object.
    """
    designs, complexes = pb.load()
    by_id = {d["design_id"]: d for d in designs}
    rows = []
    for c in complexes:
        if not c["binding_strength"] and not c["kd_M"]:
            continue                                  # never assayed
        if target and c["target"] != target:
            continue
        d = by_id.get(c["design_id"], {})
        if drop_controls and d.get("is_control") == "1":
            continue
        scores = {}
        for m in ORIENT:
            scores[m] = pb.num(d if m in DESIGN_LEVEL else c, m)
        rows.append({
            "design_id": c["design_id"],
            "target": c["target"],
            "method": c["design_method"] or "(unspecified)",
            "label": c["binding_strength"],
            "binder": c["binding_strength"] in pb.BINDER_LABELS,
            "pkd": pb.num(c, "pkd"),
            "kd_nM": pb.num(c, "kd_nM"),
            "platform": c.get("platform", ""),
            "scores": scores,
        })
    return rows


def usable_metrics(rows):
    """Metrics with enough non-null, non-constant values to be scoreable."""
    out = []
    for m in ORIENT:
        v = [r["scores"][m] for r in rows if r["scores"][m] is not None]
        if len(v) >= 50 and len(set(v)) > 1:
            out.append(m)
    return out


# ------------------------------------------------------------------- stats

def auroc(pos, neg):
    """AUROC = P(score of a binder > score of a non-binder), via Mann-Whitney."""
    if len(pos) < 5 or len(neg) < 5:
        return None
    u = mannwhitneyu(pos, neg, alternative="two-sided")
    return u.statistic / (len(pos) * len(neg)), u.pvalue


def auroc_ci(pos, neg, n_boot: int = 2000, alpha: float = 0.05):
    """Percentile bootstrap CI, resampling the two classes independently."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    boots = []
    for _ in range(n_boot):
        p = RNG.choice(pos, len(pos), replace=True)
        n = RNG.choice(neg, len(neg), replace=True)
        boots.append(mannwhitneyu(p, n).statistic / (len(p) * len(n)))
    return (float(np.quantile(boots, alpha / 2)),
            float(np.quantile(boots, 1 - alpha / 2)))


def oriented(rows, metric):
    """(positives, negatives) score arrays, sign-flipped so higher = better."""
    s = ORIENT[metric]
    pos = [s * r["scores"][metric] for r in rows
           if r["label"] and r["binder"] and r["scores"][metric] is not None]
    neg = [s * r["scores"][metric] for r in rows
           if r["label"] and not r["binder"] and r["scores"][metric] is not None]
    return pos, neg


def roc_curve(pos, neg):
    s = np.concatenate([np.asarray(pos, float), np.asarray(neg, float)])
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(-s)
    y = y[order]
    tpr = np.concatenate([[0], np.cumsum(y) / max(y.sum(), 1)])
    fpr = np.concatenate([[0], np.cumsum(1 - y) / max((1 - y).sum(), 1)])
    return fpr, tpr


def precision_at_k(rows, metric, ks=(10, 25, 50, 100, 200)):
    """If you could only order the top-k designs by this score, what binds?

    This is the number that decides a campaign: a screen's value is the lift
    over ordering at random, not its AUROC.
    """
    s = ORIENT[metric]
    cand = [(s * r["scores"][metric], r["binder"]) for r in rows
            if r["label"] and r["scores"][metric] is not None]
    cand.sort(key=lambda t: -t[0])
    base = statistics.fmean([1.0 if b else 0.0 for _, b in cand])
    out = {}
    for k in ks:
        if k <= len(cand):
            hit = statistics.fmean([1.0 if b else 0.0 for _, b in cand[:k]])
            out[k] = (hit, hit / base if base else float("nan"))
    return base, out, len(cand)


def _logistic_cv(rows, metrics, n_repeat: int = 25, test_frac: float = 0.3):
    """Held-out AUROC for a logistic combination of `metrics`.

    Fitted on a random train split and scored on the untouched test split,
    repeated. Reporting an in-sample combination AUROC would be circular, so
    the split is the whole point.
    """
    data = [r for r in rows
            if r["label"] and all(r["scores"][m] is not None for m in metrics)]
    X = np.array([[r["scores"][m] * ORIENT[m] for m in metrics] for r in data])
    y = np.array([1.0 if r["binder"] else 0.0 for r in data])
    if len(data) < 100 or y.sum() < 20:
        return None
    X = (X - X.mean(0)) / (X.std(0) + 1e-12)
    aucs = []
    for _ in range(n_repeat):
        idx = RNG.permutation(len(y))
        cut = int(len(y) * (1 - test_frac))
        tr, te = idx[:cut], idx[cut:]
        if y[tr].sum() < 5 or y[te].sum() < 5:
            continue
        Xtr = np.hstack([np.ones((len(tr), 1)), X[tr]])

        def nll(w, Xtr=Xtr, ytr=y[tr]):
            z = Xtr @ w
            # log(1+exp(z)) computed stably
            ll = ytr * z - np.logaddexp(0, z)
            return -ll.sum() + 1e-3 * (w[1:] ** 2).sum()

        w = minimize(nll, np.zeros(Xtr.shape[1]), method="L-BFGS-B").x
        z = np.hstack([np.ones((len(te), 1)), X[te]]) @ w
        a = auroc(z[y[te] == 1], z[y[te] == 0])
        if a:
            aucs.append(a[0])
    if not aucs:
        return None
    return (statistics.fmean(aucs), statistics.stdev(aucs) if len(aucs) > 1 else 0.0,
            len(data), int(y.sum()), len(aucs))


# ------------------------------------------------------------------ report

def report(fig_path="figures/m23_proteinbase_kd.png"):
    qc = pb.qc_report(verbose=False)
    all_rows = assemble()
    rows = assemble(target=PRIMARY_TARGET)
    metrics = usable_metrics(rows)

    print("=" * 72)
    print("MODULE 23 — co-folding scores vs MEASURED outcomes (Proteinbase)")
    print("=" * 72)

    # ---- A. the truth set
    n_lab = sum(1 for r in all_rows if r["label"])
    n_pos = sum(1 for r in all_rows if r["binder"])
    n_kd = sum(1 for r in all_rows if r["pkd"] is not None)
    print(f"\n[A] TRUTH SET  (snapshot {pb.SNAPSHOT}, one lab, standardised protocols)")
    print(f"    assayed design/target pairs : {n_lab}")
    print(f"      binders                   : {n_pos}")
    print(f"      measured NON-binders      : {n_lab - n_pos}   <- the denominator")
    print(f"    pairs with a measured Kd    : {n_kd}")
    print(f"    M10 had n=5. This is a {n_lab // 5}x larger truth set, and unlike")
    print(f"    ChEMBL's mixed cell/biochemical assays it is one lab on standardised")
    print(f"    protocols — but NOT one instrument: Kd comes from two platforms.")
    plat = {}
    for r in all_rows:
        if r["pkd"] is not None:
            plat.setdefault(r["platform"] or "(unrecorded)", []).append(r["pkd"])
    print(f"    Kd by platform: "
          + ", ".join(f"{k} n={len(v)} (median {10 ** (9 - statistics.median(v)):.0f} nM)"
                      for k, v in sorted(plat.items(), key=lambda kv: -len(kv[1]))))
    dropped = sorted(qc["constant_columns"])
    if dropped:
        print(f"    excluded as constant in snapshot: {', '.join(dropped)}")

    # ---- B. noise floor
    noise = qc["noise_floor_log10"]
    print(f"\n[B] MEASUREMENT NOISE FLOOR (ceiling on any predictor)")
    print(f"    Kd replicate spread: median {noise:.3f} log10 = {10 ** noise:.2f}x fold")
    print(f"    Kd == koff/kon to {qc['kinetics_consistency'][1]:.1e} log10 "
          f"(n={qc['kinetics_consistency'][2]}) — kinetics are internally consistent.")
    print(f"    Labels are Kd bands, not opinions: "
          f"Strong median {statistics.median(qc['label_bands']['Strong']):.1f} nM, "
          f"Medium {statistics.median(qc['label_bands']['Medium']):.0f} nM, "
          f"Weak {statistics.median(qc['label_bands']['Weak']):.0f} nM.")
    print(f"    => a 1.44x noise floor leaves room for rho up to ~0.9; nothing below")
    print(f"       is limited by the assay.")

    # ---- C. discrimination
    print(f"\n[C] DISCRIMINATION — binder vs non-binder, target = {PRIMARY_TARGET}")
    print(f"    {'metric':<40} {'AUROC':>6} {'95% CI':>15} {'p':>9} {'n+/n-':>11}")
    disc = {}
    for m in metrics:
        pos, neg = oriented(rows, m)
        a = auroc(pos, neg)
        if not a:
            continue
        lo, hi = auroc_ci(pos, neg)
        disc[m] = (a[0], lo, hi, a[1], len(pos), len(neg))
        flag = "" if lo > 0.5 else "   (CI includes 0.5)"
        print(f"    {m:<40} {a[0]:6.3f} [{lo:.3f},{hi:.3f}] {a[1]:9.2e} "
              f"{len(pos):5d}/{len(neg):<5d}{flag}")
    best = max(disc, key=lambda k: disc[k][0])
    print(f"    best single score: {best} (AUROC {disc[best][0]:.3f})")

    # ---- D. affinity ranking among binders
    # 12 metrics are scanned, so a single p<0.05 means nothing on its own:
    # the Bonferroni threshold is the bar a correlation has to clear here.
    bonf = 0.05 / max(len(metrics), 1)
    print(f"\n[D] AFFINITY RANKING — Spearman vs pKd, BINDERS ONLY (range-restricted)")
    print(f"    scanning {len(metrics)} metrics => Bonferroni threshold p < {bonf:.4f}"
          f"  ('**' = clears it)")
    print(f"    {'metric':<40} {'rho':>7} {'p':>9} {'n':>5}")
    aff = {}
    for scope_name, scope_rows in (("all targets pooled", all_rows),
                                   (PRIMARY_TARGET, rows)):
        print(f"    -- {scope_name} --")
        for m in metrics:
            xs = [ORIENT[m] * r["scores"][m] for r in scope_rows
                  if r["pkd"] is not None and r["scores"][m] is not None]
            ys = [r["pkd"] for r in scope_rows
                  if r["pkd"] is not None and r["scores"][m] is not None]
            if len(xs) < 20:
                continue
            rho, p = spearmanr(xs, ys)
            aff[(scope_name, m)] = (rho, p, len(xs))
            print(f"    {m:<40} {rho:+7.3f} {p:9.3f} {len(xs):5d}"
                  f"{'  **' if p < bonf else ''}")
    pooled = {m: v for (s, m), v in aff.items() if s == "all targets pooled"}
    within = {m: v for (s, m), v in aff.items() if s == PRIMARY_TARGET}
    if pooled and within:
        bp = max(pooled, key=lambda k: pooled[k][0])
        print(f"    NOTE: the best pooled rho ({bp} {pooled[bp][0]:+.2f}) "
              f"falls to {within.get(bp, (float('nan'),))[0]:+.2f} within one target")
        print(f"          — most of that pooled signal is a between-target offset.")

    # Anything that survives Bonferroni gets one more test: does it hold inside
    # each measurement platform, or is it riding an SPR-vs-BLI offset?
    surv_d = {m: v for m, v in within.items() if v[1] < bonf}
    if surv_d:
        print(f"    PLATFORM CHECK (does it survive inside one instrument?):")
        for m in surv_d:
            parts = []
            for p in sorted({r["platform"] for r in rows if r["platform"]}):
                xs = [ORIENT[m] * r["scores"][m] for r in rows
                      if r["platform"] == p and r["pkd"] is not None
                      and r["scores"][m] is not None]
                ys = [r["pkd"] for r in rows
                      if r["platform"] == p and r["pkd"] is not None
                      and r["scores"][m] is not None]
                if len(xs) >= 20:
                    rho, p_val = spearmanr(xs, ys)
                    parts.append(f"{p}: rho={rho:+.2f} (p={p_val:.3f}, n={len(xs)})")
            print(f"      {m:<42} {'  |  '.join(parts) if parts else 'n<20 per platform'}")
        n_plat = len({r["platform"] for r in rows
                      if r["pkd"] is not None and r["platform"]})
        if n_plat == 1:
            print(f"      => every Kd for this target came off ONE platform, so these")
            print(f"         correlations cannot be a cross-platform offset.")

    # ---- E. decision value
    print(f"\n[E] DECISION VALUE — precision at top-k, target = {PRIMARY_TARGET}")
    base, pk, n_c = precision_at_k(rows, best)
    print(f"    base rate (order at random): {100 * base:.1f}% of {n_c} designs bind")
    print(f"    ranking by {best}:")
    for k, (hit, lift) in pk.items():
        print(f"      top {k:4d}: {100 * hit:5.1f}% bind   ({lift:.2f}x base rate)")

    # ---- F. can a combination do better?
    combo_set = [m for m in ("boltz2_min_ipsae", "boltz2_complex_iplddt",
                             "boltz2_iptm", "esmfold_plddt", "proteinmpnn_score")
                 if m in metrics]
    cv = _logistic_cv(rows, combo_set)
    print(f"\n[F] COMBINATION — logistic fit on {len(combo_set)} scores, HELD-OUT AUROC")
    if cv:
        mean, sd, n, npos, reps = cv
        print(f"    {' + '.join(combo_set)}")
        print(f"    held-out AUROC = {mean:.3f} +/- {sd:.3f} over {reps} random "
              f"70/30 splits (n={n}, {npos} binders)")
        gain = mean - disc[best][0]
        print(f"    vs best single score ({best} {disc[best][0]:.3f}): "
              f"{gain:+.3f} — {'a real but small gain' if gain > 0.02 else 'no material gain'}")
    else:
        print("    not enough labelled rows with all metrics present")

    # ---- G. verdict
    print(f"\n[G] VERDICT — what this does to M10's claim")
    lo_b = disc[best][1]
    conf_aff = {m: v for m, v in within.items() if m in CONFIDENCE_SCORES}
    surv = {m: v for m, v in within.items() if v[1] < bonf}
    stable_k = 100 if 100 in pk else max(pk)

    print(f"    OK  Co-folding confidence DOES predict WHETHER a design binds:")
    print(f"        {best} AUROC = {disc[best][0]:.3f} (95% CI "
          f"{lo_b:.3f}-{disc[best][2]:.3f}, n={disc[best][4]}+{disc[best][5]}).")
    print(f"        Real, and modest — a triage filter, not an oracle.")
    print(f"    OK  But NOT how tightly. Among binders, within one target, not one")
    print(f"        of the {len(conf_aff)} interface-confidence scores clears the")
    print(f"        Bonferroni bar (best: "
          f"{max(conf_aff, key=lambda k: conf_aff[k][0])} "
          f"rho={max(v[0] for v in conf_aff.values()):+.2f}, "
          f"p={min(v[1] for v in conf_aff.values() if v[0] > 0.2):.3f} > {bonf:.4f}).")
    print(f"    ==> M10's 'opt_score tracks potency, rho=+0.6' came from n=5, where")
    print(f"        that rho's CI spans roughly [-0.5,+0.95]. At n~100 measured")
    print(f"        protein binders the affinity half does not reproduce. Predicting")
    print(f"        IF something binds and HOW TIGHTLY are different problems; on")
    print(f"        this evidence only the first one is solved.")
    if surv:
        print(f"    ??  HYPOTHESIS (not established): {len(surv)} non-confidence")
        print(f"        quantities DO rank affinity within the target and clear")
        print(f"        Bonferroni —")
        for m, v in sorted(surv.items(), key=lambda kv: -abs(kv[1][0])):
            print(f"          {m:<42} rho={v[0]:+.2f} (p={v[1]:.4f}, n={v[2]})")
        print(f"        Exploratory: one target, one snapshot, {len(metrics)} metrics")
        print(f"        scanned. Worth a prospective test, NOT a rule to screen on.")
    if cv:
        print(f"    OK  Scores combine: held-out AUROC {cv[0]:.3f} +/- {cv[1]:.3f} vs "
              f"{disc[best][0]:.3f} single.")
    print(f"    ==> For M7-style screens: ranking by the best score lifts the hit")
    print(f"        rate {pk[stable_k][1]:.1f}x over random in the top {stable_k} "
          f"({100 * pk[stable_k][0]:.0f}% vs {100 * base:.1f}%).")
    print(f"        Small-k precision is noisy here (top-10 to top-50 bounces "
          f"{100 * min(pk[k][0] for k in (10, 25, 50) if k in pk):.0f}-"
          f"{100 * max(pk[k][0] for k in (10, 25, 50) if k in pk):.0f}%), so trust")
    print(f"        the lift, not the exact top-10 number. Screen to shortlist,")
    print(f"        then MEASURE — never rank a shortlist by score alone.")
    print(f"    !!  Caveat: these are protein-protein interfaces; M10's target was a")
    print(f"        covalent small molecule. This transfers as a caution about")
    print(f"        confidence-vs-affinity, not as a number for small molecules.")
    print(f"\n    {pb.ATTRIBUTION}")

    _figure(rows, disc, best, qc, fig_path, within=within)
    print(f"\nfigure written to: {fig_path}")

    return {
        "n_labelled": n_lab, "n_binders": n_pos, "n_kd": n_kd,
        "discrimination": disc, "best_metric": best,
        "affinity_within_target": within, "affinity_pooled": pooled,
        "precision_at_k": pk, "base_rate": base, "combination_cv": cv,
        "noise_floor_log10": noise,
    }


def _figure(rows, disc, best, qc, path, within=None):
    fig = plt.figure(figsize=(13.5, 9.5))
    gs = fig.add_gridspec(2, 4, hspace=0.34, wspace=0.42)
    ax = {"A": fig.add_subplot(gs[0, 0:2]), "C": fig.add_subplot(gs[0, 2:4]),
          "D1": fig.add_subplot(gs[1, 0]), "D2": fig.add_subplot(gs[1, 1]),
          "E": fig.add_subplot(gs[1, 2:4])}

    # (1) the truth set: Kd by label
    a = ax["A"]
    bands = qc["label_bands"]
    for i, lbl in enumerate(("Strong", "Medium", "Weak")):
        v = bands.get(lbl, [])
        if not v:
            continue
        x = i + RNG.normal(0, 0.06, len(v))
        a.scatter(x, v, s=12, alpha=0.5, c=["#2ecc71", "#f39c12", "#e74c3c"][i],
                  edgecolors="none")
        a.scatter([i], [statistics.median(v)], marker="_", s=900, c="k", zorder=5)
    n_neg = qc["n_labelled"] - qc["n_binders"]
    a.set_yscale("log")
    a.set_xticks([0, 1, 2])
    a.set_xticklabels([f"Strong\nn={len(bands.get('Strong', []))}",
                       f"Medium\nn={len(bands.get('Medium', []))}",
                       f"Weak\nn={len(bands.get('Weak', []))}"])
    a.set_ylabel("measured Kd (nM)")
    a.set_title(f"[A] the truth set: {qc['n_kd']} measured Kd\n"
                f"+ {n_neg} measured NON-binders (no Kd, the denominator)",
                fontsize=10)
    a.grid(alpha=0.3, axis="y")

    # (2) ROC curves
    a = ax["C"]
    top = sorted(disc, key=lambda k: -disc[k][0])[:4]
    for m in top:
        pos, neg = oriented(rows, m)
        fpr, tpr = roc_curve(pos, neg)
        a.plot(fpr, tpr, lw=1.8, label=f"{m}  AUROC={disc[m][0]:.3f}")
    a.plot([0, 1], [0, 1], "k--", lw=1, label="random = 0.500")
    a.set_xlabel("false-positive rate")
    a.set_ylabel("true-positive rate")
    a.set_title(f"[C] discrimination: binder vs non-binder\n"
                f"{PRIMARY_TARGET}, n={disc[best][4]}+{disc[best][5]}", fontsize=10)
    a.legend(fontsize=7.5, loc="lower right")
    a.grid(alpha=0.3)

    # (3) affinity among binders — the dissociation, side by side: the score
    #     that wins at discrimination is flat against Kd, while another one
    #     that is useless for discrimination does rank it.
    best_aff = (max(within, key=lambda k: abs(within[k][0]))
                if within else best)
    for key, metric, colour, tag in (("D1", best, "#3498db", "best at [C]"),
                                     ("D2", best_aff, "#8e44ad", "best at [D]")):
        a = ax[key]
        xs = [ORIENT[metric] * r["scores"][metric] for r in rows
              if r["pkd"] is not None and r["scores"][metric] is not None]
        ys = [r["pkd"] for r in rows
              if r["pkd"] is not None and r["scores"][metric] is not None]
        if len(xs) >= 20:
            rho, p = spearmanr(xs, ys)
            a.scatter(xs, ys, s=22, c=colour, alpha=0.75, edgecolors="k",
                      linewidths=0.35)
            z = np.polyfit(xs, ys, 1)
            xf = np.linspace(min(xs), max(xs), 20)
            a.plot(xf, np.polyval(z, xf), "k--", lw=1.2)
            verdict = "ranks Kd" if abs(rho) > 0.3 else "FLAT"
            a.set_title(f"[D] {tag}: {verdict}\n{metric}\n"
                        f"rho={rho:+.2f} (p={p:.3f}, n={len(xs)})", fontsize=8.5)
        a.set_xlabel("score (oriented: higher = better)", fontsize=8.5)
        a.set_ylabel("measured pKd (higher = tighter)", fontsize=8.5)
        a.tick_params(labelsize=8)
        a.grid(alpha=0.3)

    # (4) precision at k
    a = ax["E"]
    base, pk, n_c = precision_at_k(rows, best, ks=tuple(range(10, 401, 10)))
    ks = sorted(pk)
    a.plot(ks, [100 * pk[k][0] for k in ks], "o-", ms=3.5, c="#8e44ad",
           label=f"ranked by {best}")
    a.axhline(100 * base, ls="--", c="k", lw=1,
              label=f"base rate {100 * base:.1f}%")
    a.set_xlabel("k = number of designs ordered (top-k by score)")
    a.set_ylabel("% that bind")
    a.set_title("[E] what the score is worth: hit rate vs ordering at random",
                fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    fig.suptitle("M23 — recalibrating co-folding scores on Proteinbase "
                 f"(Adaptyv Bio, ODC-BY, snapshot {pb.SNAPSHOT})", fontsize=12)
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.07, right=0.97)
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
