"""Does the in-silico (RibonanzaNet) OpenKnot score predict the experimental one?

Compares, per design: in-silico OpenKnot score (RibonanzaNet-predicted 2A3 -> OpenKnot
score) against the experimentally measured OpenKnot score from the published
OpenKnotBench SHAPE data.
"""
import pickle, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

SUCCESS = 90.0  # the paper's success cutoff

res = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else "/home/user/work/rnet_preds_r13.pkl", "rb"))
for r in res:  # GC content: the obvious confound, since GC-rich designs both look
    r["gc"] = (r["seq"].count("G") + r["seq"].count("C")) / len(r["seq"])  # unreactive AND fold better
df = pd.DataFrame([{k: v for k, v in r.items()
                    if k not in ("pred_react", "exp_react", "seq", "struct")} for r in res])

# ---------- 0. implementation check: does my scorer reproduce the published score? ----------
ok = df["oks_exp_published"].notna()
print("=" * 78)
print("0. SCORER REPRODUCTION CHECK (my OpenKnot implementation on the *experimental* data)")
print(f"   n = {ok.sum()}   Pearson r = {pearsonr(df.oks_exp_recomputed[ok], df.oks_exp_published[ok])[0]:.5f}"
      f"   MAE = {np.abs(df.oks_exp_recomputed[ok] - df.oks_exp_published[ok]).mean():.3f}")

# ---------- 1. profile-level accuracy ----------
pc = []
for r in res:
    e = r["exp_react"]; m = np.isfinite(e)
    if m.sum() > 20 and r["sn_filter"] == 1:
        pc.append(pearsonr(r["pred_react"][m], e[m])[0])
print("\n1. PER-DESIGN SHAPE PROFILE ACCURACY (predicted vs measured 2A3, SN_filter=1)")
print(f"   n = {len(pc)}   mean Pearson r = {np.mean(pc):.3f}   median = {np.median(pc):.3f}"
      f"   frac r>0.7 = {np.mean(np.array(pc) > 0.7):.2f}")

# ---------- 2. score-level agreement ----------
d = df[(df.sn_filter == 1) & df.oks_exp_published.notna()].copy()
d["exp"] = d.oks_exp_published
d["hit"] = d.exp > SUCCESS

print(f"\n2. IN-SILICO vs EXPERIMENTAL OpenKnot SCORE   (SN_filter=1, n={len(d)})")
for rnd, g in d.groupby("round"):
    rho = spearmanr(g.oks_insilico, g.exp).correlation
    per_t = [spearmanr(gg.oks_insilico, gg.exp).correlation
             for _, gg in g.groupby("puzzle") if len(gg) >= 20]
    print(f"   Round {rnd}: n={len(g):5d}  pooled Spearman={rho:.3f}   "
          f"within-target Spearman (median of {len(per_t)} targets) = {np.nanmedian(per_t):.3f}"
          f"  [range {np.nanmin(per_t):.2f}..{np.nanmax(per_t):.2f}]")

# ---------- 3. the decision that matters: does top-k in-silico enrich for real successes? ----------
print(f"\n3. SCREENING ENRICHMENT  (pick top-k designs per target by in-silico score;")
print(f"   how many are experimental successes, OpenKnot score > {SUCCESS:.0f}?)")
for rnd, g in d.groupby("round"):
    print(f"\n   --- Round {rnd} ---")
    base = g.hit.mean()
    print(f"   base rate (all designs)              : {base * 100:5.1f}%   n={len(g)}")
    for k in (1, 5, 10, 25):
        hits, tot, per_target_success = 0, 0, []
        for _, gg in g.groupby("puzzle"):
            if len(gg) < k:
                continue
            top = gg.nlargest(k, "oks_insilico")
            hits += top.hit.sum(); tot += len(top)
            per_target_success.append(bool(top.hit.any()))
        print(f"   top-{k:<3d} by in-silico score           : {hits / tot * 100:5.1f}%   "
              f"(enrichment x{hits / tot / base:.2f})   "
              f"targets with >=1 success in top-{k}: {np.mean(per_target_success) * 100:.0f}%")
    # random baseline for the "targets solved" column
    rng = np.random.default_rng(0)
    for k in (10,):
        rand = []
        for _ in range(200):
            s = [gg.sample(min(k, len(gg)), random_state=rng.integers(1e9)).hit.any()
                 for _, gg in g.groupby("puzzle")]
            rand.append(np.mean(s))
        print(f"   random-{k} baseline                   : targets with >=1 success: "
              f"{np.mean(rand) * 100:.0f}% +/- {np.std(rand) * 100:.0f}%")

# ---------- 4. compare predictors ----------
def auroc(score, label):
    score, label = np.asarray(score, float), np.asarray(label, bool)
    m = np.isfinite(score)
    score, label = score[m], label[m]
    if label.all() or not label.any():
        return np.nan
    order = np.argsort(score)
    ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score) + 1)
    n1, n0 = label.sum(), (~label).sum()
    return (ranks[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

print("\n4. PREDICTOR COMPARISON — AUROC for 'design will be an experimental success'")
print("   (pooled within round; RNet_F1 is the filter reported in the paper;")
print("    GC content is the trivial baseline any sequence model could be exploiting)")
print(f"   {'round':>6} {'n':>6} {'in-silico OKS':>14} {'RNet_F1':>9} {'RNet_F1_cp':>11} {'GC only':>9}")
for rnd, g in d.groupby("round"):
    print(f"   {rnd:>6} {len(g):>6} {auroc(g.oks_insilico, g.hit):>14.3f} "
          f"{auroc(g.rnet_f1, g.hit):>9.3f} {auroc(g.rnet_f1_cp, g.hit):>11.3f} "
          f"{auroc(g.gc, g.hit):>9.3f}")

print("\n4b. IS THE MODEL BEATING GC CONTENT?  AUROC of the in-silico score computed")
print("    inside narrow GC strata (GC held ~constant, so it cannot be the explanation)")
for rnd, g in d.groupby("round"):
    out = []
    for lo in np.arange(0.45, 0.80, 0.05):
        s = g[(g.gc >= lo) & (g.gc < lo + 0.05)]
        if len(s) >= 200 and 0.05 < s.hit.mean() < 0.95:
            out.append((f"{lo:.2f}-{lo + 0.05:.2f}", len(s), auroc(s.oks_insilico, s.hit)))
    print(f"   Round {rnd}: " + "   ".join(f"GC {a} (n={b}) AUROC={c:.3f}" for a, b, c in out))

# ---------- 5. per-method view ----------
print("\n5. BY DESIGN METHOD (Round 3, SN_filter=1)")
g3 = d[d["round"] == 3]
rows = []
for meth, gg in g3.groupby("method"):
    if len(gg) < 50:
        continue
    rows.append(dict(method=meth, n=len(gg),
                     exp_success_pct=100 * gg.hit.mean(),
                     mean_exp=gg.exp.mean(), mean_insilico=gg.oks_insilico.mean(),
                     spearman=spearmanr(gg.oks_insilico, gg.exp).correlation,
                     auroc=auroc(gg.oks_insilico, gg.hit)))
print(pd.DataFrame(rows).sort_values("exp_success_pct", ascending=False).to_string(index=False, float_format="%.3f"))

d.to_csv("/home/user/work/analysis_table.csv", index=False)
print("\nsaved per-design table -> /home/user/work/analysis_table.csv")
