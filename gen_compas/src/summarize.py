"""Collect every number that goes into REPORT.md, straight from results/."""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def brute_force():
    out = []
    for f in sorted(glob.glob(os.path.join(common.RESULTS, "bruteforce_r*.npz"))):
        d = np.load(f)
        phi = d["phi"]
        core = common.which_core(phi)
        cc = core[core >= 0]
        n_switch = int(np.sum(np.diff(cc) != 0)) if len(cc) > 1 else 0
        out.append(dict(file=os.path.basename(f), ns=len(phi) / 1000.0,
                        transitions=n_switch,
                        frac_A=float((cc == 0).mean()) if len(cc) else 0.0))
    return out


def reference():
    from analysis import metad_fes, state_free_energies, barrier_height, \
        free_energy_along_phi
    f, edges = metad_fes()
    if f is None:
        return None
    centers = 0.5 * (edges[1:] + edges[:-1])
    fphi = free_energy_along_phi(f, edges)
    core = common.which_core(centers)
    return dict(dG=float(state_free_energies(f, edges)),
                barrier=barrier_height(f, edges),
                fphi_A_min=float(np.nanmin(fphi[core == 0])),
                fphi_B_min=float(np.nanmin(fphi[core == 1])))


def convergence(f0, f1):
    """Change in the phi profile between two metadynamics snapshots."""
    from analysis import KT
    a, b = np.load(f0).T, np.load(f1).T
    def prof(F):
        p = np.exp(-F / KT).sum(1)
        g = -KT * np.log(p)
        return g - g.min()
    pa, pb = prof(a), prof(b)
    sel = pb < 12.0
    return float(np.abs(pa[sel] - pb[sel]).max())


def main():
    print("=" * 72)
    print("REFERENCE  (well-tempered metadynamics on phi,psi; same force field)")
    print("=" * 72)
    r = reference()
    if r:
        print(f"  dG(B - A)                      {r['dG']:+.2f} kcal/mol")
        print(f"  phi~0 saddle above basin A     {r['barrier']['from_A']:.2f} kcal/mol")
        print(f"  phi~0 saddle above basin B     {r['barrier']['from_B']:.2f} kcal/mol")
    snaps = sorted(glob.glob(os.path.join(common.RESULTS,
                                          "metad_snapshot_*.npy")))
    if len(snaps) >= 2:
        print(f"  drift of the phi profile between the last two snapshots: "
              f"{convergence(snaps[0], snaps[-1]):.2f} kcal/mol")

    print()
    print("=" * 72)
    print("BRUTE-FORCE BASELINE  (plain unbiased MD, same force field)")
    print("=" * 72)
    for b in brute_force():
        print(f"  {b['file']}: {b['ns']:.1f} ns, {b['transitions']} A<->B "
              f"transitions, {100 * b['frac_A']:.0f}% of time in A")

    print()
    for tag in ("main", "rep1"):
        pa = os.path.join(common.RESULTS, f"analysis_{tag}.json")
        pl = os.path.join(common.RESULTS, f"gencompas_{tag}_log.json")
        if not os.path.exists(pa):
            continue
        with open(pa) as fh:
            a = json.load(fh)
        print("=" * 72)
        print(f"GEN-COMPAS  [{tag}]")
        print("=" * 72)
        if os.path.exists(pl):
            with open(pl) as fh:
                L = json.load(fh)
            print(f"  total MD                       {L['total_ns']:.2f} ns")
            print(f"  wall clock (1 CPU core)        {L['wall_s'] / 60:.0f} min")
            print("  per iteration:")
            print("    it   ns    cum ns   targets  TMD  shot pts  "
                  "on-separatrix  <q_emp>")
            for r_ in L["log"]:
                print(f"    {r_['iteration']:2d}  {r_['ns_this_iteration']:5.2f}"
                      f"  {r_['ns_cumulative']:6.2f}   {r_['n_targets']:5d}"
                      f"  {r_['n_tmd']:4d}  {r_['n_shot_points']:7d}"
                      f"  {100 * r_['frac_separatrix']:11.0f}%"
                      f"  {r_['empirical_q_mean']:8.2f}")
        print(f"  frames collected               {a['n_frames']}")
        print(f"  dG(B - A)                      {a['dG_gencompas']:+.2f} kcal/mol"
              + (f"   [reference {a['dG_metad']:+.2f}]" if "dG_metad" in a else ""))
        bg = a["barrier_gencompas"]
        print(f"  phi~0 saddle above A           {bg['from_A']:.2f} kcal/mol"
              + (f"   [reference {a['barrier_metad']['from_A']:.2f}]"
                 if "barrier_metad" in a else ""))
        if "fel_rmse" in a:
            print(f"  landscape RMSE vs reference    {a['fel_rmse']:.2f} kcal/mol")
        if "mfpt_ab_ns" in a:
            print(f"  MSM first passage time A->B    {a['mfpt_ab_ns']:.0f} ns")
            print(f"  MSM first passage time B->A    {a['mfpt_ba_ns']:.0f} ns")
        print(f"  transition-state ensemble      {a['tse_n']} frames, "
              f"phi(10/50/90) = "
              f"{[round(v) for v in a['tse_phi']] if a['tse_phi'] else '-'}")
        print(f"  pathways found                 {len(a['paths'])}")
        for p in a["paths"]:
            mid = [x for x in p["points"] if abs(x["q"] - 0.5) < 0.1]
            if mid:
                print(f"    channel {p['channel']}: {p['n_frames']} frames, "
                      f"q=1/2 at phi={mid[0]['phi']:.0f}, psi={mid[0]['psi']:.0f}")

        pv = os.path.join(common.RESULTS, f"validation_{tag}.json")
        if os.path.exists(pv):
            with open(pv) as fh:
                v = json.load(fh)
            print(f"  --- committor validation ({v['ns']:.1f} ns of fresh MD, "
                  f"{v['n_shots']} shots per structure) ---")
            print(f"  mean |q_pred - q_measured|     {v['mae']:.3f} "
                  f"(binomial noise floor {v['expected_mae']:.3f})")
            print(f"  correlation                    {v['corr']:.3f}")
            if v.get("tse_mean") is not None:
                print(f"  measured committor of the predicted TSE  "
                      f"{v['tse_mean']:.3f} +/- {v['tse_std']:.3f}")
        print()


if __name__ == "__main__":
    main()
