"""
Module 5 - the result that traditional MD analysis does not surface.

Question: does a distant point mutation change the population of the CRYPTIC
ALLOSTERIC POCKET, 30 A away from where the mutation sits?

Design note.  The 32 apo runs are a balanced 4 variants x 8 crystal starts
grid, so the variant comparison is PAIRED on the starting structure.  That
matters: the opening is rare and bursty, so an unpaired test on 8 noisy runs
would mostly measure which crystal structure you happened to start from.  The
paired Wilcoxon signed-rank test on the 8 matched starts is the right one, and
its floor p-value at n=8 one-sided is 1/2^6 = 0.0156.

Everything in this module is a DIRECT OBSERVABLE - occupancies and committed
crossing counts - so none of it depends on the MSM converging (M3).
"""
import numpy as np
from scipy import stats

from .common import (APO_STARTS, HOLO_START, LIGANDS, LIGAND_NAMES, TRAJ_NS,
                     VARIANTS, apo_indices, banner, holo_indices, load)
from .m2_cryptic_pocket import (OPEN_THR, collective_variables,
                                committed_crossings)


def _by_start(cv, names, thr=OPEN_THR):
    idx = {n: k for k, n in enumerate(names)}
    return {v: np.array([(cv[idx[f"{s}_{v}"]] > thr).mean() for s in APO_STARTS])
            for v in VARIANTS}


def main():
    banner("M5  does a distal mutation open the cryptic pocket?  (direct observables)")
    _, names = load()
    _, cv1, _ = collective_variables()

    print(f"  open-state occupancy, mouth distance L220-N276 > {OPEN_THR} A\n")
    print(f"  {'variant':8s} {'mutations':22s} {'open %':>8s} {'crossings':>10s} "
          f"{'runs w/ event':>14s}")
    muts = {"TEM-1": "(wild type)", "TEM-2": "Q39K",
            "TEM-10": "R164S, E240K", "TEM-52": "E104K, M182T, G238S"}
    occ, crs = {}, {}
    for v in VARIANTS:
        ks = apo_indices(names, v)
        occ[v] = np.array([(cv1[k] > OPEN_THR).mean() for k in ks])
        crs[v] = np.array([committed_crossings(cv1[k])[0] for k in ks])
        print(f"  {v:8s} {muts[v]:22s} {100*occ[v].mean():8.4f} {crs[v].sum():10d} "
              f"{int((crs[v]>0).sum()):>10d}/{len(ks)}")

    ratio = occ["TEM-10"].mean() / occ["TEM-1"].mean()
    print(f"\n  TEM-10 opens the pocket {ratio:.0f}x more than wild-type TEM-1.")

    print("\n  PAIRED by crystal starting structure (the correct test):")
    M = _by_start(cv1, names)
    print(f"  {'start':6s}" + "".join(f"{v:>10s}" for v in VARIANTS))
    for i, s in enumerate(APO_STARTS):
        print(f"  {s:6s}" + "".join(f"{100*M[v][i]:10.3f}" for v in VARIANTS))
    print()
    for v in ["TEM-1", "TEM-2", "TEM-52"]:
        w = stats.wilcoxon(M["TEM-10"], M[v], alternative="greater")
        hi = int((M["TEM-10"] > M[v]).sum())
        print(f"    TEM-10 > {v:7s}  Wilcoxon p = {w.pvalue:.4f}   higher in {hi}/8 starts")
    print("    (0.0156 is the floor p-value for a one-sided n=8 signed-rank test)")

    print(f"""
  INTERPRETATION
    R164S sits in the omega loop, at the active site.  E240K sits on the
    beta-turn.  Neither touches the H11/H12 pocket - that pocket is ~30 A away
    on the far face of the protein.  Yet the TEM-10 pair shifts its open
    population by {ratio:.0f}x.  That is long-range allosteric coupling, read out as a
    change in the POPULATION of a rare state rather than as a change in any
    average structure.

    Testable prediction: the non-competitive H11/H12 inhibitors (the 1PZO /
    1PZP chemotype of Horn & Shoichet) should be MORE potent against TEM-10
    than against TEM-1, because their binding site is pre-populated.  A simple
    IC50 comparison on the two purified enzymes settles it.""")

    # ---- does substrate at the active site change the distal pocket? -------
    print("\n  SUBSTRATE EFFECT (apo vs holo) - and why it is NOT conclusive here")
    idx = {n: k for k, n in enumerate(names)}
    for v in VARIANTS:
        kh = holo_indices(names, v)
        holo = np.array([(cv1[k] > OPEN_THR).mean() for k in kh])
        matched = np.array([(cv1[idx[f'{s}_{v}']] > OPEN_THR).mean()
                            for s in sorted({HOLO_START[l] for l in LIGANDS})])
        allapo = occ[v]
        print(f"    {v:7s} apo(all 8) {100*allapo.mean():7.4f} %  "
              f"apo(matched starts) {100*matched.mean():7.4f} %  "
              f"holo {100*holo.mean():7.4f} %")
    a = np.array([(cv1[idx[f'{s}_TEM-10']] > OPEN_THR).mean() for s in ["1XPB", "1HTZ"]])
    h = np.array([(cv1[k] > OPEN_THR).mean() for k in holo_indices(names, "TEM-10")])
    u = stats.mannwhitneyu(a, h, alternative="greater")
    print(f"""
    Taken at face value TEM-10 looks like it loses the open pocket on substrate
    binding ({100*occ['TEM-10'].mean():.3f} % -> {100*h.mean():.3f} %), which would be textbook
    negative coupling between the active site and the distal pocket.  But the
    holo runs were only ever started from 1XPB and 1HTZ, and those two starts
    are among the LOW openers in the apo set.  Against the matched apo baseline
    the effect shrinks to {a.mean()/max(h.mean(),1e-9):.1f}x with n = 2 vs 4 and p = {u.pvalue:.2f}.

    So: suggestive, underpowered, and confounded by starting structure.  It
    needs holo runs launched from the same 8 starts as the apo set before it
    can be called an effect.  Reported here as a hypothesis, not a result.""")
    return {"ratio_tem10_tem1": float(ratio),
            "wilcoxon_p": float(stats.wilcoxon(M["TEM-10"], M["TEM-1"],
                                               alternative="greater").pvalue)}


if __name__ == "__main__":
    main()
