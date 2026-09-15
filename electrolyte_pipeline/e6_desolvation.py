"""
E6 — desolvation: one coordinate, three different "energy differences"   [digest p11]

The digest's table distinguishes
   (1) binding electronic energy      — E5 gives it (complex vs. fragments);
   (2) reaction electronic-energy barrier along a path in ONE model;
   (3) process free-energy barrier along a DEFINED coordinate in an environment.

Here the *same* coordinate — the Li…O(DME) distance — is followed in two
different models so that (2) and (3) can be put side by side without being
confused with each other:

   A. gas-phase relaxed scan, Li+(DME), B3LYP/def2-SVP (+ GFN2-xTB):
      Li–O1 constrained from 1.8 to 5.0 Å, everything else relaxed.
      This is an *electronic energy profile* of the 1:1 complex in vacuum.
      Result labels: reaction-electronic-energy along that coordinate.

   B. liquid-phase potential of mean force from the E2 Li–O(DME) RDF,
      W(r) = −kT ln g(r), classical force field, all other DME/TFSI present.
      This is a *free energy* along the same coordinate, but only where the
      unbiased trajectory sampled it; the large-r plateau is the bulk
      reference and the first minimum is the coordinated state.
      ⚠️ No umbrella sampling / metadynamics was run: barrier regions with
      g(r) ≈ 0 are not resolved (they show as gaps / large W).

Neither number is "the desolvation energy" of the electrolyte; the plot says
what each one is.
"""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG
from electrolyte_pipeline.e5_qc_clusters import (starting_structures, _mol, XC, OPT_BASIS,
                                                 HARTREE_KCAL, _xtb_energy, _write_xyz, OUT)

KT_KCAL = 0.0019872041 * 298.15


def gas_scan(distances=(1.8, 2.0, 2.2, 2.5, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0), force=False) -> dict:
    cache = os.path.join(OUT, "e6_gas_scan.json")
    if os.path.exists(cache) and not force:
        return json.load(open(cache))
    from pyscf import dft
    from pyscf.geomopt.geometric_solver import optimize
    sym, xyz = starting_structures()["DME_bidentate"]
    li = len(sym) - 1
    o_idx = [i for i, s in enumerate(sym) if s == "O"]
    # start from the optimised bidentate complex if E5 has produced it
    opt_xyz = os.path.join(OUT, "DME_bidentate_opt.xyz")
    if os.path.exists(opt_xyz):
        xyz = np.loadtxt(opt_xyz, skiprows=2, usecols=(1, 2, 3))
    o1 = min(o_idx, key=lambda i: np.linalg.norm(xyz[i] - xyz[li]))
    E, Ex, geoms, d2 = [], [], [], []
    cur = xyz.copy()
    for d in distances:
        # displace Li along the O1->Li direction to the target distance, then relax with the constraint
        v = cur[li] - cur[o1]; v /= np.linalg.norm(v)
        cur[li] = cur[o1] + d * v
        mol = _mol(sym, cur, +1, OPT_BASIS)
        mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = 3
        # geomeTRIC constraint file (1-based indices)
        cpath = os.path.join(OUT, "e6_constraint.txt")
        with open(cpath, "w") as fh:
            fh.write(f"$freeze\n$set\ndistance {o1+1} {li+1} {d:.3f}\n")
        mol_eq = optimize(mf, constraints=cpath, maxsteps=120, convergence_energy=1e-5,
                          convergence_grms=5e-4, convergence_gmax=1e-3, convergence_drms=2e-3, convergence_dmax=4e-3)
        cur = mol_eq.atom_coords(unit="Angstrom")
        mf2 = dft.RKS(mol_eq); mf2.xc = XC; mf2.grids.level = 3
        E.append(mf2.kernel())
        Ex.append(_xtb_energy(sym, cur, +1))
        other = [i for i in o_idx if i != o1][0]
        d2.append(float(np.linalg.norm(cur[li] - cur[other])))
        geoms.append(cur.copy())
        print(f"[E6 scan] Li–O1 = {d:.2f} Å   E = {(E[-1]-E[0])*HARTREE_KCAL:7.2f} kcal/mol rel.   Li–O2 = {d2[-1]:.2f} Å")
    E = np.array(E); rel = (E - E.min()) * HARTREE_KCAL
    relx = ((np.array(Ex) - min(Ex)) * HARTREE_KCAL).tolist() if None not in Ex else None
    res = {"coordinate": "Li…O1(DME) distance, all else relaxed", "level": f"{XC}/{OPT_BASIS}, gas phase",
           "d_A": list(distances), "E_rel_kcal": rel.tolist(), "xtb_E_rel_kcal": relx,
           "Li_O2_A": d2, "E_abs_Ha": E.tolist(),
           "note": "electronic energy along a constrained coordinate in vacuum: a (2)-type quantity, not a free energy"}
    for d, g in zip(distances, geoms):
        _write_xyz(os.path.join(OUT, f"e6_scan_{d:.1f}.xyz"), sym, g)
    json.dump(res, open(cache, "w"), indent=1)
    return res


def liquid_pmf(label="LiTFSI_r10", workdir=WORK) -> dict | None:
    e2 = os.path.join(workdir, label, "e2_rdf.json")
    if not os.path.exists(e2):
        return None
    p = json.load(open(e2))["pairs"]["O(DME)"]
    r, g = np.array(p["curve"]["r"]), np.array(p["curve"]["g"])
    with np.errstate(divide="ignore"):
        W = -KT_KCAL * np.log(g)
    ok = g > 0
    W[~ok] = np.nan
    # reference: bulk plateau (r > 7 Å) -> W = 0
    ref = np.nanmean(W[(r > 7.0)])
    W -= ref
    i_min = int(np.nanargmin(np.where(r < p["r_min"], W, np.nan)))
    # barrier between the first minimum and r_min region (maximum of W in [r_peak, r_min+1])
    sel = (r > p["r_peak"]) & (r < p["r_min"] + 1.0) & ok
    i_bar = int(np.where(sel)[0][np.nanargmax(W[sel])]) if sel.any() else i_min
    res = {"label": label, "coordinate": "Li…O(DME) distance in the liquid (classical FF, all species present)",
           "r_A": r.tolist(), "W_kcal": np.where(np.isnan(W), None, W).tolist(),
           "W_min_kcal": float(W[i_min]), "r_W_min_A": float(r[i_min]),
           "W_barrier_kcal": float(W[i_bar] - W[i_min]), "r_barrier_A": float(r[i_bar]),
           "unsampled_fraction": float(1 - ok[(r > 1.5) & (r < 7)].mean()),
           "note": "⚠️ W = -kT ln g(r) from unbiased MD: a free energy along r where sampled; no enhanced sampling; "
                   "the barrier value is a lower bound if g≈0 bins are missing"}
    return res


def report(force=False) -> dict:
    print("E6  desolvation along Li…O(DME): gas-phase electronic scan vs. liquid PMF (same coordinate, different quantities)")
    gas = gas_scan(force=force)
    liq = liquid_pmf()
    d, e = np.array(gas["d_A"]), np.array(gas["E_rel_kcal"])
    print(f"    gas scan: minimum at Li–O1 = {d[np.argmin(e)]:.1f} Å; E(5.0 Å) − E(min) = {e[-1]-e.min():.1f} kcal/mol "
          f"(electronic, vacuum, 1:1 complex, O2 still bound at {gas['Li_O2_A'][-1]:.2f} Å)")
    if liq:
        print(f"    liquid PMF ({liq['label']}): W_min = {liq['W_min_kcal']:.2f} kcal/mol at {liq['r_W_min_A']:.2f} Å; "
              f"barrier to leave the shell ≈ {liq['W_barrier_kcal']:.2f} kcal/mol at {liq['r_barrier_A']:.2f} Å "
              f"({liq['unsampled_fraction']*100:.0f}% of bins unsampled)")
    _plot(gas, liq)
    return {"gas": gas, "liquid": liq}


def _plot(gas, liq):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(gas["d_A"], gas["E_rel_kcal"], "o-", label=f"{gas['level']}")
    if gas.get("xtb_E_rel_kcal"):
        ax[0].plot(gas["d_A"], gas["xtb_E_rel_kcal"], "s--", ms=4, label="GFN2-xTB")
    ax[0].set_xlabel("Li…O1 (Å), constrained"); ax[0].set_ylabel("ΔE electronic (kcal/mol)")
    ax[0].set_title("(2) gas-phase Li+(DME) electronic-energy profile", fontsize=9); ax[0].legend(fontsize=8)
    if liq:
        r = np.array(liq["r_A"]); W = np.array([np.nan if w is None else w for w in liq["W_kcal"]])
        ax[1].plot(r, W, "k-")
        ax[1].axvline(liq["r_W_min_A"], ls=":", color="gray"); ax[1].axvline(liq["r_barrier_A"], ls=":", color="r")
        ax[1].set_xlim(1.5, 8); ax[1].set_ylim(-3, 4)
        ax[1].set_title(f"(3) liquid PMF W=-kT ln g, {liq['label']} (classical FF)", fontsize=9)
    ax[1].set_xlabel("Li…O(DME) (Å)"); ax[1].set_ylabel("W (kcal/mol), bulk = 0")
    fig.suptitle("E6  same coordinate, different quantities: not interchangeable 'desolvation energies'", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e6_desolvation.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(force="--force" in sys.argv)
