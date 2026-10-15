"""Build and gate the R2 conformational-state generator from the CG trajectories.

Seeds 101, 102 build (state centres, pi, generator); seed 103 is the held-out validation seed for
the Chapman-Kolmogorov / direct-propagation test; after validation the generator is re-estimated on
all seeds (101-104) with the same states and lag, and both records are kept.

Gate operationalization (pre-declared thresholds from PROTOCOL.md; concrete rules below):
  timescales  : slowest two timescales resolvable at the lag (t_k > tau) change <= 20 % from tau to 2tau and 4tau
  CK          : for 2tau, 4tau, 8tau, max |P_pred - P_emp| over initial states with pi >= 0.05 is <= 0.10 AND
                compatible with seed-103 block-bootstrap uncertainty (deviation <= 3 SE, or <= 0.03 absolute)
  sampling    : states covering 95 % of sum(pi_i w_i): ESS >= 200, >= 50 committed entries (residence >= tau),
                committed entries present in both build seeds
Writes results/r2_model_dt{dt}.npz and results/r2_gates_dt{dt}.json.
"""
import json
import os
import sys

import numpy as np

import md_analysis as ma
import msm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
LAGS = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)   # frames (10 ps each); 500, 1000 added (deviation: tau_z ~ 20 ns)
STATE_COUNTS = (4, 8, 12)
BETAS = (2.0, 5.0, 10.0)


def replica_runs(name):
    r = ma.load(name)
    out = []
    for j in range(r["z"].shape[1]):
        out.append(dict(z=r["z"][:, j], rg=r["rg"][:, j], r_end=r["r_end"][:, j], nc=r["nc"][:, j]))
    return out, r["t_ps"][1] - r["t_ps"][0]


def weights_per_state(reps, dtrajs, k, beta):
    num = np.zeros(k)
    cnt = np.zeros(k)
    for r, s in zip(reps, dtrajs):
        num += np.bincount(s, weights=np.exp(-beta * (r["z"] - ma.Z_MIN)), minlength=k)
        cnt += np.bincount(s, minlength=k)
    return num / np.maximum(cnt, 1)


def evaluate(sm, build, val, frame_ps, rng):
    k = sm.k
    dB = [sm.assign(r) for r in build["reps"]]
    dV = [sm.assign(r) for r in val["reps"]]
    pi = msm.stationary_from(dB, k)
    w5 = weights_per_state(build["reps"], dB, k, 5.0)
    contrib = pi * w5 / np.sum(pi * w5)
    order = np.argsort(contrib)[::-1]
    cover = order[: int(np.searchsorted(np.cumsum(contrib[order]), 0.95) + 1)]
    seed_of = np.repeat(np.arange(len(build["n_per_seed"])), build["n_per_seed"])
    ess = np.zeros(k)
    for sd in range(len(build["n_per_seed"])):
        idx = np.where(seed_of == sd)[0]
        ess += msm.ess_and_entries([dB[i] for i in idx], k, 1)[0]
    results = []
    for lag in LAGS:
        gens = {}
        for m in (1, 2, 4):
            C = msm.count_matrix(dB, lag * m, k)
            g = msm.fit_generator(C, pi, lag * m * frame_ps)
            g.dt = frame_ps
            gens[m] = g
        ts = {m: gens[m].timescales() for m in gens}
        res_n = int(np.sum(ts[1][:2] > lag * frame_ps))
        ts_change = [float(np.max(np.abs(ts[m][:max(res_n, 1)] / ts[1][:max(res_n, 1)] - 1))) for m in (2, 4)]
        ts_ok = res_n == 0 or max(ts_change) <= 0.20
        ck = msm.ck_test(gens[1], dV, lag, mults=(2, 4, 8), rng=rng)
        ck_ok = all(c["max_dev"] <= 0.10 and (c["max_dev"] <= 3 * c["se_at_max"] or c["max_dev"] <= 0.03) for c in ck)
        # sampling gate: per-seed counting so that "from >= 2 seeds" can be checked
        ent = np.zeros(k, int)
        seeds_with = np.zeros(k, int)
        for sd in range(len(build["n_per_seed"])):
            idx = np.where(seed_of == sd)[0]
            n_ent = msm.committed_entries([dB[i] for i in idx], k, lag)
            ent += n_ent
            seeds_with += (n_ent > 0).astype(int)
        samp_ok = bool(all(ess[i] >= 200 and ent[i] >= 50 and seeds_with[i] >= 2 for i in cover))
        results.append(dict(lag_frames=lag, lag_ps=lag * frame_ps, timescales_ps={m: ts[m][:3].tolist() for m in ts},
                            resolvable=res_n, ts_change=ts_change, ts_ok=bool(ts_ok), ck=ck, ck_ok=bool(ck_ok),
                            cover_states=cover.tolist(), ess=ess.tolist(), entries=ent.tolist(),
                            seeds_with_entries=seeds_with.tolist(), sampling_ok=samp_ok,
                            passed=bool(ts_ok and ck_ok and samp_ok), loglik=gens[1].loglik,
                            converged=all(g.converged for g in gens.values())))
        if results[-1]["passed"]:
            break
    return pi, contrib, results


def equilibrium_blocks(names, dt_frame_ps):
    """Seed x half blocks: contact occupancy, <w>, tau_z, z-PMF agreement over the important region."""
    blocks = []
    zs = []
    for n in names:
        r = ma.load(n)
        h = r["z"].shape[0] // 2
        for part, sl in (("a", slice(0, h)), ("b", slice(h, None))):
            z = r["z"][sl]
            w = np.exp(-5 * (z - ma.Z_MIN))
            taus = [ma.tau_int(z[:, j], dt_frame_ps) for j in range(0, z.shape[1], 8)]
            blocks.append(dict(block=f"{n}{part}", p_contact=float((z < ma.Z_MIN + 0.2).mean()),
                               w_mean=float(w.mean()), tau_z_ps=float(np.mean(taus))))
            zs.append(z.ravel())
    edges = np.arange(0.2, max(z.max() for z in zs) + 0.1, 0.1)
    P = np.array([np.histogram(z, edges)[0] for z in zs], float)
    Ppool = P.sum(0) / P.sum()
    cen = 0.5 * (edges[1:] + edges[:-1])
    wz = np.exp(-5 * (cen - ma.Z_MIN))
    imp = np.zeros(len(cen), bool)
    o = np.argsort(Ppool)[::-1]
    imp[o[: np.searchsorted(np.cumsum(Ppool[o]), 0.99) + 1]] = True
    pw = Ppool * wz / np.sum(Ppool * wz)
    o = np.argsort(pw)[::-1]
    imp[o[: np.searchsorted(np.cumsum(pw[o]), 0.95) + 1]] = True
    Fpool = -np.log(Ppool[imp])
    dmax = []
    for b in range(len(zs)):
        Pb = P[b] / P[b].sum()
        if (Pb[imp] == 0).any():
            dmax.append(np.inf)
            continue
        Fb = -np.log(Pb[imp])
        off = np.sum(Ppool[imp] * (Fb - Fpool)) / np.sum(Ppool[imp])
        dmax.append(float(np.max(np.abs(Fb - Fpool - off))))
    for b, d in zip(blocks, dmax):
        b["pmf_max_dev_kT"] = d
    region = (float(cen[imp].min()), float(cen[imp].max()))
    return blocks, region, bool(max(dmax) <= 1.0)


def main(dt):
    rng = np.random.default_rng(11)
    names = [f"dt{dt}_s{s}" for s in (101, 102, 103, 104)]
    reps = {}
    frame_ps = None
    for n in names:
        reps[n], frame_ps = replica_runs(n)
    build = dict(reps=reps[names[0]] + reps[names[1]], n_per_seed=[len(reps[names[0]]), len(reps[names[1]])])
    val = dict(reps=reps[names[2]])
    report = dict(dt_ps=dt, frame_ps=frame_ps, build=names[:2], validation=names[2], extra=names[3], candidates=[])
    chosen = None
    for k in STATE_COUNTS:
        sm = msm.StateModel(build["reps"], k, seed=k, stride=20)
        print(f"k={k}: states built", flush=True)
        pi, contrib, res = evaluate(sm, build, val, frame_ps, rng)
        report["candidates"].append(dict(k=k, pi_build=pi.tolist(), current_share=contrib.tolist(),
                                         centres_std=sm.centres.tolist(), lags=res))
        print(f"k={k}: " + "; ".join(f"lag {r['lag_ps']:.0f}ps ts{'+' if r['ts_ok'] else '-'} ck{'+' if r['ck_ok'] else '-'} "
                                      f"samp{'+' if r['sampling_ok'] else '-'} maxCK {max(c['max_dev'] for c in r['ck']):.3f}"
                                      for r in res), flush=True)
        if res[-1]["passed"]:
            chosen = (k, sm, res[-1]["lag_frames"])
            break
    report["eligible"] = chosen is not None
    if chosen is None:
        # fallback (not ranked): 12 states, lag with the smallest worst-case CK error
        cand = report["candidates"][-1]
        best = min(cand["lags"], key=lambda r: max(c["max_dev"] for c in r["ck"]))
        sm = msm.StateModel(build["reps"], cand["k"], seed=cand["k"], stride=20)
        chosen = (cand["k"], sm, best["lag_frames"])
        report["fallback_note"] = "no candidate passed all gates; conditional (unvalidated) generator"
    k, sm, lag = chosen
    report["chosen"] = dict(k=k, lag_frames=lag, lag_ps=lag * frame_ps)
    # build-only generator (validation record) and pooled re-estimate
    allreps = sum((reps[n] for n in names), [])
    dA = [sm.assign(r) for r in allreps]
    dB = dA[: len(build["reps"])]
    out = {}
    for tag, d, rr in (("build", dB, build["reps"]), ("pooled", dA, allreps)):
        pi = msm.stationary_from(d, k)
        g = msm.fit_generator(msm.count_matrix(d, lag, k), pi, lag * frame_ps)
        out[tag] = dict(pi=pi, L_per_s=g.L * 1e12, timescales_ns=g.timescales() / 1e3,
                        w={b: weights_per_state(rr, d, k, b) for b in BETAS})
    blocks, region, pmf_ok = equilibrium_blocks(names, frame_ps)
    report["equilibrium_blocks"] = blocks
    report["pmf_region_nm"] = region
    report["pmf_ok"] = pmf_ok
    report["pi_by_seed"] = [msm.stationary_from(dA[i * 64:(i + 1) * 64], k).tolist() for i in range(4)]
    report["timescales_ns_build"] = out["build"]["timescales_ns"][:4].tolist()
    report["timescales_ns_pooled"] = out["pooled"]["timescales_ns"][:4].tolist()
    os.makedirs(RES, exist_ok=True)
    np.savez(os.path.join(RES, f"r2_model_dt{dt}.npz"), k=k, lag_ps=lag * frame_ps,
             **{f"{t}_pi": out[t]["pi"] for t in out}, **{f"{t}_L": out[t]["L_per_s"] for t in out},
             **{f"{t}_w{b:g}": out[t]["w"][b] for t in out for b in BETAS}, eligible=report["eligible"])
    json.dump(report, open(os.path.join(RES, f"r2_gates_dt{dt}.json"), "w"), indent=1, default=float)
    print("chosen", report["chosen"], "eligible", report["eligible"], "pmf_ok", pmf_ok)
    return report


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.5)
