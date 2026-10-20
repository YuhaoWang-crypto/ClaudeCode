"""
End-to-end feasibility run for the three DCSC brain-holes.

    python3 -m dcsc.run_demo            # full run (~3-5 min on 4 cores)
    python3 -m dcsc.run_demo --quick    # smaller sizes, smoke test

Writes dcsc/RESULTS.md. Every number is computed against synthetic ground
truth (see synth.py for what is and is not realistic).
"""
import argparse
import os
import time

import numpy as np

from . import codec as C
from . import delta as DL
from . import index as IX
from . import partial as PV
from . import synth as S

OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def table(rows, cols, fmt=None):
    fmt = fmt or {}
    say("| " + " | ".join(cols) + " |")
    say("|" + "---|" * len(cols))
    for r in rows:
        say("| " + " | ".join(fmt.get(c, "{}").format(r[c]) for c in cols) + " |")
    say()


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


# ---------------------------------------------------------------------------
def part1_codec(world, enc, rng, n_cells, n_samples):
    say("## 1. Factorized State Barcode (64 B)\n")
    A, M, E = S.random_factors(n_cells, rng)
    leaf = rng.integers(0, len(S.LEAVES), n_cells)
    X = world.sample_norm(leaf, A, M, E, rng.integers(0, world.n_batches, n_cells), rng)
    out = enc.encode(X)
    Z = np.hstack([A, M, E])
    lin_true = world.leaf_lineage[leaf]
    lin_hat = world.leaf_lineage[out["leaf"]]
    say(f"Identity tree (single cells, n={n_cells}): lineage acc = {np.mean(lin_true == lin_hat):.3f}, "
        f"leaf acc = {np.mean(leaf == out['leaf']):.3f}\n")

    Xs, ls, Zs = PV.sample_states(world, rng, n_samples, 100)
    zs = enc.encode(Xs, n_cells=100)["z"]
    overlap = [("TNF_NFkB", 0), ("Heat_shock", 0), ("EMT", 1), ("Prolif_G2M", 1), ("Oxphos", 2)]
    rows = []
    nv = enc.naive_scores(X, out["leaf"])
    nvs = enc.naive_scores(Xs, enc.classify(Xs))
    for name, est, truth in [("joint posterior (DCSC), single cell", out["z"], Z),
                             ("naive gene-set score, single cell", nv, Z),
                             ("joint posterior (DCSC), pseudobulk n=100", zs, Zs),
                             ("naive gene-set score, pseudobulk n=100", nvs, Zs)]:
        r = [corr(est[:, j], truth[:, j]) for j in range(S.K_LATENT)]
        leak = [abs(corr(est[:, S.PROGRAMS.index(p)], truth[:, S.P + S.NM + e])) for p, e in overlap]
        rows.append(dict(method=name, mod=np.mean(r[:S.P]), mem=np.mean(r[S.P:S.P + S.NM]),
                         env=np.mean(r[S.P + S.NM:]), leak=np.mean(leak),
                         null=np.sqrt(2 / np.pi) / np.sqrt(len(est)),
                         tern=np.mean(C.ternary(est[:, :S.P]) == C.ternary(truth[:, :S.P]))))
    say("Recovery of each block (Pearson r with truth, mean over coordinates) and leakage:\n")
    table(rows, ["method", "mod", "mem", "env", "leak", "null", "tern"],
          {k: "{:.3f}" for k in ["mod", "mem", "env", "leak", "null", "tern"]})
    say("`leak` = mean |r| between an estimated program and the TRUE artifact it shares 12 genes "
        "with by construction (TNF_NFkB/Heat_shock↔Dissociation, EMT/Prolif_G2M↔Culture_2D, "
        "Oxphos↔Mito); `null` = expected |r| from sampling noise alone. "
        "`tern` = ternary-state accuracy over 32 programs.\n")

    raw = world.G * 4
    say(f"Size: {C.CODE_BYTES} B/cell vs {raw:,} B for a dense float32 profile of {world.G} genes "
        f"(×{raw // C.CODE_BYTES}; ×{20000 * 4 // C.CODE_BYTES} for 20k genes). "
        f"100 M cells = {100e6 * C.CODE_BYTES / 1e9:.1f} GB → fits in RAM.\n")
    say("Example barcodes:\n")
    for i in range(3):
        say(f"- `{C.describe(out['code'][i])}`")
    say()

    # query modes
    codes = out["code"]
    true_t = C.ternary(A)
    stats = {"mod": [], "same_id": []}
    for qi in rng.choice(n_cells, 200, replace=False):
        for mode in stats:
            idx, _ = IX.search(codes, codes[qi], k=11, mode=mode)
            idx = idx[idx != qi][:10]
            stats[mode].append((np.abs(true_t[idx] - true_t[qi]).sum(1).mean(),
                                np.mean(world.leaf_lineage[leaf[idx]] != world.leaf_lineage[leaf[qi]]),
                                np.mean(leaf[idx] == leaf[qi])))
    rnd = np.abs(true_t[rng.integers(0, n_cells, 2000)] - true_t[rng.integers(0, n_cells, 2000)]).sum(1).mean()
    rows = [dict(mode=m, prog_dist=np.mean([s[0] for s in v]), other_lineage=np.mean([s[1] for s in v]),
                 same_leaf=np.mean([s[2] for s in v])) for m, v in stats.items()]
    rows.append(dict(mode="random pairs", prog_dist=rnd, other_lineage=2 / 3, same_leaf=1 / 9))
    say("Query modes (top-10 neighbours of 200 query cells): true program-state distance "
        "(lower = same functional state), share from another lineage, share of same leaf:\n")
    table(rows, ["mode", "prog_dist", "other_lineage", "same_leaf"],
          {k: "{:.2f}" for k in ["prog_dist", "other_lineage", "same_leaf"]})


# ---------------------------------------------------------------------------
def synthetic_codes(n, rng, noise=0.3):
    A, M, E = S.random_factors(n, rng)
    z = np.hstack([A, M, E]) + rng.normal(0, noise, (n, S.K_LATENT))
    leaf = rng.integers(0, len(S.LEAVES), n)
    return C.pack(leaf, rng.normal(0, 30, (n, 12)), z)


def part2_speed(rng, n_db, n_hnsw):
    say("## 2. Search speed on 64-B codes (single core)\n")
    db = synthetic_codes(n_db, rng)
    qs = synthetic_codes(64, rng)
    rows = []

    t = time.perf_counter()
    for q in qs[:8]:
        IX.mod_distance(db, q)
    dt = (time.perf_counter() - t) / 8
    rows.append(dict(engine="numpy bitwise_count", n=n_db, per_query_ms=dt * 1e3, cmp_per_s=n_db / dt))

    if IX.c_kernel():
        IX.mod_distance_c(db, qs[0])
        t = time.perf_counter()
        for q in qs[:16]:
            IX.mod_distance_c(db, q)
        dt = (time.perf_counter() - t) / 16
        assert np.array_equal(IX.mod_distance_c(db, qs[1]), IX.mod_distance(db, qs[1]))
        rows.append(dict(engine="C POPCNT kernel", n=n_db, per_query_ms=dt * 1e3, cmp_per_s=n_db / dt))
    c_rate = rows[-1]["cmp_per_s"]

    try:
        import faiss
        faiss.omp_set_num_threads(1)
        fl = IX.faiss_mod_index(db, "flat")
        t = time.perf_counter()
        fl.search(IX.mod_block(qs), 10)
        dt = (time.perf_counter() - t) / len(qs)
        rows.append(dict(engine="faiss IndexBinaryFlat", n=n_db, per_query_ms=dt * 1e3, cmp_per_s=n_db / dt))

        sub = db[:n_hnsw]
        t = time.perf_counter()
        hn = IX.faiss_mod_index(sub, "hnsw")
        build = time.perf_counter() - t
        hn.hnsw.efSearch = 128
        t = time.perf_counter()
        Dh, _ = hn.search(IX.mod_block(qs), 10)
        dt = (time.perf_counter() - t) / len(qs)
        Df, _ = IX.faiss_mod_index(sub, "flat").search(IX.mod_block(qs), 10)
        recall = np.mean(Dh <= Df[:, -1:])     # ties are common: recall by distance
        rows.append(dict(engine=f"faiss IndexBinaryHNSW (build {build:.0f}s, recall@10 {recall:.2f})",
                         n=n_hnsw, per_query_ms=dt * 1e3, cmp_per_s=float("nan")))
    except ImportError:
        say("(faiss not installed — skipped)\n")
    table(rows, ["engine", "n", "per_query_ms", "cmp_per_s"],
          {"n": "{:,}", "per_query_ms": "{:.3f}", "cmp_per_s": "{:.3g}"})
    say(f"Extrapolation: brute-force scan of 100 M codes ≈ {100e6 / c_rate:.2f} s/query/core with the "
        "C kernel, i.e. ~0.1 s on 8 cores; HNSW gives ms-level latency at some recall cost. "
        "Ternary Hamming has only 0–64 distinct values, so ties are frequent: re-rank the "
        "top few thousand with the continuous posterior or Mem/Env blocks.\n")


# ---------------------------------------------------------------------------
def part3_delta(world, enc, rng, n_q, n_cells):
    say("## 3. Δ-state signatures + residual pursuit (cocktail recipe)\n")
    L, Lsig = DL.build_library(world, enc, rng, leaf=0, n=300)
    eq = DL.equivalence(L)
    Ln = L / np.linalg.norm(L, axis=1, keepdims=True)
    G = Ln @ Ln.T
    np.fill_diagonal(G, 0)
    pairs = [(DL.NAMES[i], DL.NAMES[j]) for i, j in zip(*np.where(np.triu(eq, 1)))]
    say(f"Library: {len(DL.NAMES)} stimuli measured in {S.LEAVES[0][1]} (pseudobulk 300 vs 300 cells). "
        f"Mutual coherence max cos = {G.max():.2f}; near-duplicate pairs (cos>0.9): "
        f"{', '.join(f'{a}~{b}' for a, b in pairs) or 'none'}.\n")
    E = DL.effect_matrix()
    settings = [
        ("1 stimulus, same cell type", 1, 0, None, 0.0),
        ("2 stimuli, same cell type", 2, 0, None, 0.0),
        ("3 stimuli, same cell type", 3, 0, None, 0.0),
        ("2 stimuli, other lineage (CD4_T)", 2, 3, None, 0.0),
        ("2 stimuli + dissociation/2D artifact", 2, 0, np.array([1.2, 1.0, 0.0, 0.3]), 0.0),
        ("2 stimuli + synergy (non-additive)", 2, 0, None, 1.0),
    ]
    rows = []
    for name, k, leaf, env, syn in settings:
        res = dict(setting=name, exact=0, equiv=0, top_k_cos=0, dose_mae=[], unexpl=[], flagged=0)
        naive_ok = 0
        for _ in range(n_q):
            ids = rng.choice(len(DL.NAMES), k, replace=False)
            dose = rng.uniform(0.6, 1.4, k)
            a = dose @ E[ids]
            if syn:
                v = np.zeros(S.P)
                v[rng.choice(S.P, 3, replace=False)] = rng.choice([-1, 1], 3) * syn
                a = a + dose.prod() * v
            d, d_naive = DL.pseudobulk_delta(world, enc, leaf, a, rng, n_cells, env_shift=env)
            out = DL.residual_pursuit(d[C.SL_MOD], L, Lsig)
            got = out["stimuli"]
            truth = set(ids.tolist())
            res["exact"] += set(got) == truth
            ok = len(got) == k and all(any(eq[g, t] for g in got) for t in truth)
            res["equiv"] += ok
            res["top_k_cos"] += set(np.argsort(-(Ln @ d[C.SL_MOD]))[:k].tolist()) == truth
            if set(got) == truth:
                dmap = dict(zip(got, out["dose"]))
                res["dose_mae"].append(np.mean([abs(dmap[i] - x) for i, x in zip(ids, dose)]))
            res["unexpl"].append(out["unexplained"])
            res["flagged"] += out["unexplained"] > 0.35
            if env is not None:
                naive_ok += set(DL.residual_pursuit(d_naive[C.SL_MOD], L, Lsig)["stimuli"]) == truth
        for key in ["exact", "equiv", "top_k_cos", "flagged"]:
            res[key] /= n_q
        res["dose_mae"] = np.mean(res["dose_mae"]) if res["dose_mae"] else float("nan")
        res["unexpl"] = np.mean(res["unexpl"])
        if env is not None:
            res["setting"] += f" (naive-score Δ: {naive_ok / n_q:.2f})"
        rows.append(res)
    say(f"{n_q} random queries per row; target and baseline = {n_cells} cells each; doses U(0.6,1.4):\n")
    table(rows, ["setting", "exact", "equiv", "top_k_cos", "dose_mae", "unexpl", "flagged"],
          {k: "{:.2f}" for k in ["exact", "equiv", "top_k_cos", "dose_mae", "unexpl", "flagged"]})
    say("`exact` = recovered set equals truth; `equiv` = equal up to near-duplicate stimuli; "
        "`top_k_cos` = baseline that just takes the k best single matches (k given!); "
        "`unexpl` = |residual|/|ΔS|; `flagged` = share with unexplained > 0.35 "
        "(reported as 'missing stimulus or non-additive').\n")


# ---------------------------------------------------------------------------
def part4_partial(world, enc, rng, n_test, n_cells):
    say("## 4. Asymmetric partial-view query\n")
    targets = ["IFNa", "IFNg", "TNF_NFkB", "Hypoxia"]
    tgt = [S.PROGRAMS.index(p) for p in targets]
    leaf = 0
    # calibration plate: samples with known state, measured on every modality
    Xp, lp, Zp = PV.sample_states(world, rng, 150, n_cells, leaf=leaf)
    Rp = Xp - enc.centroid[lp]
    assay = PV.CytokineAssay(rng)
    zfull_p = enc.encode(Xp, leaf=lp, n_cells=n_cells)["z"]
    y_p = assay.measure(Zp, rng)
    assay.fit_on_paired(zfull_p, y_p)

    Xt, lt, Zt = PV.sample_states(world, rng, n_test, n_cells, leaf=leaf)
    Rt = Xt - enc.centroid[lt]
    snr = np.linalg.norm(enc.W / np.sqrt(enc.D + enc.D_sys)[:, None], axis=1)
    panels = {
        "full transcriptome": np.arange(world.G),
        "targeted 800": np.argsort(-snr)[:800],
        "qPCR 12 (top-loaded)": PV.top_loaded(enc, targets, 3),
        "qPCR 12 (A-optimal design)": PV.design_panel(enc, targets, 12, n_cells),
    }
    rows, views = [], {}
    for name, genes in panels.items():
        v = C.LinearGaussianView(enc.W, enc.D, genes, n_cells, name=name, D_floor=enc.D_sys)
        v.calibrate(Rp[:, genes], Zp, C.SL_MOD)
        r, _, _ = PV.eval_view(v, Rt[:, genes], Zt, tgt)
        rows.append(r)
        views[name] = (v, genes)
    cv = assay.view()
    cv.calibrate(y_p, Zp, C.SL_MOD)
    r, _, _ = PV.eval_view(cv, assay.measure(Zt, rng), Zt, [S.PROGRAMS.index(p) for p in assay.PROGS])
    r["panel"] = "4 secreted cytokines"
    rows.append(r)
    say(f"{n_test} test samples ({n_cells}-cell pseudobulk, CD14_Mono, random programs/memory/artifacts). "
        "Same loadings W for every gene panel; only rows are sliced. Error bars = linear-Gaussian "
        "posterior (prior sd 1, single-cell noise/n + between-batch floor) × conformal factor κ "
        "fitted on a separate 150-sample calibration plate:\n")
    table(rows, ["panel", "n_obs", "n_identifiable", "rmse", "rmse_targets", "kappa", "coverage90",
                 "frac_determined", "tern_acc_determined", "tern_acc_all", "r95_targets", "cond"],
          {"rmse": "{:.2f}", "rmse_targets": "{:.2f}", "kappa": "{:.2f}", "coverage90": "{:.2f}",
           "frac_determined": "{:.2f}", "tern_acc_determined": "{:.3f}", "tern_acc_all": "{:.3f}",
           "r95_targets": "{:.2f}", "cond": "{:.3g}"})
    say("`n_identifiable` = programs whose posterior sd < 0.5 × prior sd; `rmse_targets` = the 4 "
        "programs the panel was built for (cytokine row: IFNg/TNF/IL6/IL10); `kappa` = error-bar "
        "inflation needed on the calibration plate (1 = model noise already honest); `coverage90` = "
        "share of TEST errors inside the 90% interval (target ≈0.90); `frac_determined` = share of "
        "program bits whose 90% interval sits in one ternary bin (only those bits are searched); "
        "`r95_targets` = longest semi-axis of the 95% confidence ellipsoid over the 4 target programs "
        "(programs the panel does not see fall back to the prior: radius √χ²₃₂ ≈ 6.8); `cond` = "
        "condition number of the Fisher matrix (inf = fewer observations than latent coordinates — "
        "the prior regularises, the mask reports it).\n")

    # which library stimulus is this sample?
    E = DL.effect_matrix()
    eq = DL.equivalence(E, 0.9)
    Xl, ll, _ = PV.sample_states(world, rng, len(DL.NAMES), 300, leaf=leaf, a=E, mem_on=False)
    lib = enc.encode(Xl, leaf=ll, n_cells=300)["code"]
    reps = 3
    Xq, lq, Zq = PV.sample_states(world, rng, len(DL.NAMES) * reps, n_cells, leaf=leaf,
                                  a=np.repeat(E, reps, 0), mem_on=False)
    truth = np.repeat(np.arange(len(DL.NAMES)), reps)
    Rq = Xq - enc.centroid[lq]
    rows = []
    for name in ["full transcriptome", "qPCR 12 (A-optimal design)", "4 secreted cytokines"]:
        if name.startswith("4 secreted"):
            mean, sd = cv.posterior(assay.measure(Zq, rng))
            seen = [S.PROGRAMS.index(p) for p in assay.PROGS]
        else:
            v, genes = views[name]
            mean, sd = v.posterior(Rq[:, genes])
            seen = list(range(S.P)) if name.startswith("full") else tgt
        covered = (np.abs(E[:, seen]) > C.TAU).any(1)[truth]
        codes = C.pack(lq, np.zeros((len(lq), 12)), mean)
        masks = C.mask_bits(C.determined_mask(mean[:, C.SL_MOD], sd[:, C.SL_MOD]))
        outcome = []
        for i in range(len(codes)):
            d = IX.mod_distance(lib, codes[i], masks[i])
            best = np.where(d == d.min())[0]
            if len(best) > 2:
                outcome.append("ambiguous")
            elif any(eq[b, truth[i]] for b in best):
                outcome.append("correct")
            else:
                outcome.append("wrong")
        outcome = np.array(outcome)
        for grp, sel in [("covered", covered), ("not covered", ~covered)]:
            if sel.sum() == 0:
                continue
            o = outcome[sel]
            rows.append(dict(panel=name, stimuli=f"{grp} ({sel.sum() // reps})",
                             correct=np.mean(o == "correct"), ambiguous=np.mean(o == "ambiguous"),
                             wrong=np.mean(o == "wrong")))
    say(f"Which of the {len(DL.NAMES)} library stimuli produced this sample? Library = full-transcriptome "
        f"codes; query = panel posterior, undetermined bits masked; {reps} samples per stimulus. "
        "'covered' = the stimulus moves a program the panel can see. 'ambiguous' = >2 library entries "
        "tie at the best masked distance, so the system answers 'cannot tell' instead of guessing:\n")
    table(rows, ["panel", "stimuli", "correct", "ambiguous", "wrong"],
          {k: "{:.2f}" for k in ["correct", "ambiguous", "wrong"]})


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    q = args.quick
    rng = np.random.default_rng(7)
    t0 = time.time()
    world = S.World(seed=0)
    enc = C.DCSCEncoder().fit(world, rng, n_ref=200 if q else 400, n_cal=100 if q else 200)
    say("# DCSC feasibility run — results\n")
    say(f"Synthetic world: {world.G} genes, {len(S.LEAVES)} cell types in a 2-level tree, "
        f"{S.P} programs, {S.NM} memory and {S.NE} environment factors, {world.n_batches} batches, "
        "NB counts (~2500 UMI/cell). Encoder fitted from simulated reference + calibration "
        "contrasts only.\n")
    part1_codec(world, enc, rng, 1500 if q else 4000, 60 if q else 200)
    part2_speed(rng, 200_000 if q else 2_000_000, 50_000 if q else 300_000)
    part3_delta(world, enc, rng, 10 if q else 60, 150 if q else 300)
    part4_partial(world, enc, rng, 40 if q else 150, 200)
    say(f"_Run time {time.time() - t0:.0f} s{' (quick mode)' if q else ''}._")
    if not q:
        with open(os.path.join(os.path.dirname(__file__), "RESULTS.md"), "w") as f:
            f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
