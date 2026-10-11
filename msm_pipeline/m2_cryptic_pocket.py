"""
Module 2 - define the cryptic allosteric pocket from CRYSTALLOGRAPHY, not from
the simulation, then ask how often the simulation goes there.

Why this ordering matters: if you read the reaction coordinate off the same
trajectories you then score, you have fitted the answer.  Here the opening
axis is fixed by two crystal structures that were solved independently of
this dataset --

  1BTL  apo TEM-1, cryptic pocket closed
  1PZO  TEM-1 with N,N-bis(4-chlorobenzyl)-tetrazol-5-amine bound in the
        H11/H12 pocket, a core-disrupting NON-COMPETITIVE inhibitor

-- and the MD is projected onto it.  The pocket is the one identified
experimentally by Horn & Shoichet (JMB 2004) and recovered computationally by
Bowman et al. (PNAS 2015, doi:10.1073/pnas.1417811112) and Porter et al.
(Biophys J 2019, doi:10.1016/j.bpj.2018.11.3144).

Outputs (cached): proj.npy, cv_220_276.npy, cv_223_280.npy
"""
import numpy as np

from .common import (CACHE, VARIANTS, XTAL_CLOSED, XTAL_OPEN, XTALS,
                     amb_to_idx, apo_indices, banner, cached, crystal_ca,
                     kabsch_rmsd, load)

# C-alpha pairs spanning the H10-H11 loop / H11 face against H12.
AXIS_PAIRS = [(220, 276), (223, 280), (219, 273), (223, 277), (216, 279), (218, 272)]
# residues lining the pocket, used for the RMSD validation
POCKET_RES = list(range(214, 232)) + list(range(269, 287))

OPEN_THR = 9.5      # angstrom on the L220-N276 mouth distance
CLOSED_THR = 7.0    # core-set definition for committed-crossing counting


def opening_axis():
    """Unit vector in AXIS_PAIRS distance space pointing closed -> inhibitor-bound."""
    c, o = crystal_ca(XTAL_CLOSED), crystal_ca(XTAL_OPEN)
    vc = np.array([np.linalg.norm(c(i) - c(j)) for i, j in AXIS_PAIRS])
    vo = np.array([np.linalg.norm(o(i) - o(j)) for i, j in AXIS_PAIRS])
    d = vo - vc
    return vc, d / np.linalg.norm(d), float(np.linalg.norm(d))


def _compute():
    X, _ = load()
    vc, axis, _ = opening_axis()
    ii = [amb_to_idx(a) for a, _ in AXIS_PAIRS]
    jj = [amb_to_idx(b) for _, b in AXIS_PAIRS]
    proj = np.empty((X.shape[0], X.shape[1]), np.float32)
    for k in range(X.shape[0]):
        x = np.asarray(X[k])
        proj[k] = (np.linalg.norm(x[:, ii] - x[:, jj], axis=-1) - vc) @ axis
    return proj


def collective_variables():
    """-> proj (48,T) along the crystal opening axis, plus the two raw mouth CVs."""
    X, _ = load()
    proj = cached("proj", _compute)
    cv1 = cached("cv_220_276", lambda: np.stack(
        [np.linalg.norm(np.asarray(X[k])[:, amb_to_idx(220)]
                        - np.asarray(X[k])[:, amb_to_idx(276)], axis=1)
         for k in range(X.shape[0])]))
    cv2 = cached("cv_223_280", lambda: np.stack(
        [np.linalg.norm(np.asarray(X[k])[:, amb_to_idx(223)]
                        - np.asarray(X[k])[:, amb_to_idx(280)], axis=1)
         for k in range(X.shape[0])]))
    return proj, cv1, cv2


def committed_crossings(x, lo=CLOSED_THR, hi=OPEN_THR):
    """Core-set (milestoning) count of committed closed->open transitions."""
    s = np.full(len(x), -1, np.int8)
    s[x < lo] = 0
    s[x > hi] = 1
    seq = s[s >= 0]
    if len(seq) < 2:
        return 0, 0
    d = np.diff(seq)
    return int((d == 1).sum()), int((d == -1).sum())


def main():
    banner("M2  cryptic allosteric pocket, defined by crystallography")
    X, names = load()
    vc, axis, L = opening_axis()
    print(f"  opening axis from {XTAL_CLOSED} (closed) -> {XTAL_OPEN} (inhibitor-bound)")
    for (i, j), a, b in zip(AXIS_PAIRS, vc, vc + axis * L):
        print(f"    {i:3d}-{j:3d}   {a:5.2f} -> {b:5.2f} A   ({b-a:+.2f})")
    print(f"  |axis| = {L:.2f} A   (1BTL = 0, 1PZO = {L:.2f})")

    print("\n  pocket RMSD between reference crystals "
          f"({len(POCKET_RES)} C-alpha, Kabsch):")
    ref = {p: np.array([crystal_ca(p)(a) for a in POCKET_RES]) for p in XTALS}
    print("          " + "".join(f"{p:>8s}" for p in XTALS))
    for a in XTALS:
        print(f"  {a:6s}  " + "".join(f"{kabsch_rmsd(ref[a], ref[b]):8.2f}" for b in XTALS))
    print("  -> the two inhibitor-bound structures sit ~2.2 A away from all three")
    print("     apo structures: that gap IS the cryptic opening.")

    proj, cv1, _ = collective_variables()
    apo = apo_indices(names)
    print(f"\n  apo frames with mouth distance L220-N276 > {OPEN_THR} A: "
          f"{100*(cv1[apo] > OPEN_THR).mean():.3f} %")
    print("  (a state this rare is invisible to an average structure - see M6)")

    print("\n  how far along the crystal axis does the apo ensemble travel?")
    print(f"  {'variant':8s} {'mean':>6s} {'p99':>6s} {'max':>7s}  frac > half-way")
    for v in VARIANTS:
        p = proj[apo_indices(names, v)]
        print(f"  {v:8s} {p.mean():6.2f} {np.percentile(p,99):6.2f} {p.max():7.2f}"
              f"  {100*(p > L/2).mean():8.3f} %")

    tot = sum(committed_crossings(cv1[k])[0] for k in range(len(names)))
    print(f"\n  committed closed(<{CLOSED_THR}A) -> open(>{OPEN_THR}A) crossings, "
          f"all 48 runs: {tot}")
    print("  -> enough events to compare variants (M5); NOT enough for converged")
    print("     absolute rates (M3/M4).")
    return {"axis_len": L, "crossings": tot}


if __name__ == "__main__":
    main()
