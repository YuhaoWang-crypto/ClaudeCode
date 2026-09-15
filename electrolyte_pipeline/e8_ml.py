"""
E8 — AI4Electrolytes: two different model classes, two different validations   [digest p12]

Part A — potential model (energies & forces).  GFN2-xTB is used as a stand-in
  for a "cheap potential not trained on this environment" and tested against
  B3LYP/def2-SVP on Li+ first-shell clusters cut from the E1 liquid (E2's
  shell extraction, several snapshots).  Reported: energy RMSE on *relative*
  energies, force RMSE and force-component correlation, and the same
  quantities on the 1:1 Li+–solvent complexes of E5.  The point of the digest:
  "训练数据若不含目标反应、组成或界面环境，不能默认模型可靠外推" — a potential that is
  fine for 1:1 gas-phase complexes may not be fine for the condensed-phase shell.

Part B — property model (X -> y).  A per-ion dataset is assembled from the E1
  trajectories: X = (local shell composition of one Li+ at time t, salt ratio),
  y = displacement of that Li+ over the next 100 ps.  The same regressor is
  evaluated with
     (i)  a random frame split       — frames of the same trajectory land in
                                       train AND test (the leak the digest warns about),
     (ii) leave-one-concentration-out — the honest split for "screen a new formulation",
  and against a trivial baseline (predict the training mean).
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, SERIES_BY_LABEL
from electrolyte_pipeline.e5_qc_clusters import _mol, XC, OPT_BASIS, HARTREE_KCAL, OUT

BOHR = 0.52917721092


# --------------------------------------------------------------------------
# Part A
# --------------------------------------------------------------------------
def _dft_energy_forces(sym, xyz, charge):
    from pyscf import dft
    mol = _mol(sym, xyz, charge, OPT_BASIS)
    mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = 3
    e = mf.kernel()
    g = mf.nuc_grad_method().kernel()                 # Ha/Bohr
    return e * HARTREE_KCAL, -g * HARTREE_KCAL / BOHR  # kcal/mol, kcal/mol/Å


def _xtb_energy_forces(sym, xyz, charge):
    from xtb.interface import Calculator, Param
    Z = {"H": 1, "Li": 3, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16}
    calc = Calculator(Param.GFN2xTB, np.array([Z[s] for s in sym]), np.array(xyz) / BOHR, charge=float(charge))
    calc.set_verbosity(0)
    r = calc.singlepoint()
    return r.get_energy() * HARTREE_KCAL, -r.get_gradient() * HARTREE_KCAL / BOHR


def extract_shells(label="LiTFSI_r10", workdir=WORK, n=12, r_cut=3.0, seed=0):
    """Cut n Li+ first shells (whole DME/TFSI molecules with an O within r_cut) from random frames."""
    from electrolyte_pipeline.traj import Traj
    from MDAnalysis.lib.distances import distance_array
    tr = Traj(label, workdir)
    rng = np.random.default_rng(seed)
    li = tr.sel(species="Li"); O = tr.sel(element="O")
    shells = []
    frames = rng.choice(np.arange(tr.n_frames), n, replace=False)
    for f in sorted(frames):
        ts = tr.u.trajectory[tr.start + f]
        c = rng.choice(li)
        d = distance_array(tr.u.atoms.positions[[c]], tr.u.atoms.positions[O], box=ts.dimensions)[0]
        mols = np.unique(tr.molid[O[d < r_cut]])
        idx = np.concatenate([[c], np.where(np.isin(tr.molid, mols))[0]])
        pos = tr.u.atoms.positions[idx] - tr.u.atoms.positions[c]
        box = tr.box(ts); pos -= box * np.round(pos / box)
        # make each molecule whole around its first atom
        for m in mols:
            sel = np.where(tr.molid[idx] == m)[0]
            dd = pos[sel] - pos[sel[0]]; dd -= box * np.round(dd / box); pos[sel] = pos[sel[0]] + dd
        sym = [s.capitalize() for s in tr.element[idx]]
        n_tfsi = int(np.sum(tr.species[idx] == "TFSI") // 15)
        shells.append({"sym": sym, "xyz": pos.tolist(), "charge": 1 - n_tfsi, "frame": int(f),
                       "n_DME": int(np.sum(tr.species[idx] == "DME") // 16), "n_TFSI": n_tfsi})
    return shells


def part_a(force=False) -> dict:
    cache = os.path.join(OUT, "e8_potential_check.json")
    if os.path.exists(cache) and not force:
        return json.load(open(cache))
    shells = extract_shells()
    rows = []
    for s in shells:
        t0 = time.time()
        e_d, f_d = _dft_energy_forces(s["sym"], np.array(s["xyz"]), s["charge"])
        e_x, f_x = _xtb_energy_forces(s["sym"], np.array(s["xyz"]), s["charge"])
        rows.append({"n_atoms": len(s["sym"]), "n_DME": s["n_DME"], "n_TFSI": s["n_TFSI"], "charge": s["charge"],
                     "E_dft": e_d, "E_xtb": e_x, "F_dft": f_d.tolist(), "F_xtb": f_x.tolist(),
                     "wall_s": round(time.time() - t0, 1)})
        print(f"[E8-A] shell {len(rows):2d}: {len(s['sym'])} atoms ({s['n_DME']} DME, {s['n_TFSI']} TFSI)  "
              f"F rmse {np.sqrt(np.mean((f_d-f_x)**2)):.2f} kcal/mol/Å   {rows[-1]['wall_s']} s")
    # relative energies are only comparable within the same composition; group by (n_DME, n_TFSI)
    groups = {}
    for r in rows:
        groups.setdefault((r["n_DME"], r["n_TFSI"]), []).append(r)
    e_pairs = []
    for g in groups.values():
        if len(g) > 1:
            ed = np.array([r["E_dft"] for r in g]); ex = np.array([r["E_xtb"] for r in g])
            e_pairs += list(zip(ed - ed.mean(), ex - ex.mean()))
    e_pairs = np.array(e_pairs) if e_pairs else np.zeros((0, 2))
    F_d = np.concatenate([np.array(r["F_dft"]).ravel() for r in rows]); F_x = np.concatenate([np.array(r["F_xtb"]).ravel() for r in rows])
    li_mask = np.concatenate([np.repeat(np.array([s == "Li" for s in sh["sym"]]), 3) for sh in shells])
    res = {"reference": f"{XC}/{OPT_BASIS} single points", "candidate": "GFN2-xTB (untrained on this system)",
           "n_shells": len(rows), "compositions": {f"{k[0]}DME_{k[1]}TFSI": len(v) for k, v in groups.items()},
           "force_rmse_kcal_mol_A": float(np.sqrt(np.mean((F_d - F_x) ** 2))),
           "force_rmse_on_Li_kcal_mol_A": float(np.sqrt(np.mean((F_d[li_mask] - F_x[li_mask]) ** 2))),
           "force_pearson": float(np.corrcoef(F_d, F_x)[0, 1]),
           "force_rms_reference": float(np.sqrt(np.mean(F_d ** 2))),
           "relative_energy_rmse_kcal": float(np.sqrt(np.mean((e_pairs[:, 0] - e_pairs[:, 1]) ** 2))) if len(e_pairs) else None,
           "n_relative_energy_pairs": int(len(e_pairs)),
           "rows": [{k: v for k, v in r.items() if not k.startswith("F_")} for r in rows],
           "scatter": {"F_dft": F_d[::7].round(3).tolist(), "F_xtb": F_x[::7].round(3).tolist(),
                       "E_rel": e_pairs.round(3).tolist()}}
    # 1:1 complexes from E5 for the "in-distribution" comparison
    e5 = os.path.join(OUT, "e5_pairs.json")
    if os.path.exists(e5):
        p = json.load(open(e5))["structures"]
        res["pair_complexes"] = {k: {"dft_vertical": v["dE_vertical_kcal"], "xtb_vertical": v["xtb_dE_vertical_kcal"]} for k, v in p.items()}
    json.dump(res, open(cache, "w"), indent=1)
    return res


# --------------------------------------------------------------------------
# Part B
# --------------------------------------------------------------------------
def build_dataset(workdir=WORK, horizon_ps=100.0, r_cut=3.0, stride=5):
    """Per-Li rows: features = shell composition at t, ratio; target = |Δr| over horizon."""
    from electrolyte_pipeline.traj import Traj
    from MDAnalysis.lib.distances import distance_array
    X, y, group = [], [], []
    for lab, comp in SERIES_BY_LABEL.items():
        if "Li" not in comp.counts or not os.path.exists(os.path.join(workdir, lab, "prod.dcd")):
            continue
        tr = Traj(lab, workdir)
        li = tr.sel(species="Li"); od = tr.sel(species="DME", element="O"); ot = tr.sel(species="TFSI", element="O")
        h = int(horizon_ps / tr.frame_ps)
        ratio = comp.counts["DME"] / comp.counts["Li"]
        # unwrapped Li positions
        P = []
        prev = None
        for ts in tr.frames():
            p = tr.u.atoms.positions[li].astype(float)
            if prev is not None:
                d = p - prev; p -= tr.box(ts) * np.round(d / tr.box(ts))
            prev = p.copy(); P.append(p)
        P = np.array(P)
        for f in range(0, tr.n_frames - h, stride):
            ts = tr.u.trajectory[tr.start + f]
            dd = distance_array(tr.u.atoms.positions[li], tr.u.atoms.positions[od], box=ts.dimensions)
            dt = distance_array(tr.u.atoms.positions[li], tr.u.atoms.positions[ot], box=ts.dimensions)
            n_od = (dd < r_cut).sum(1); n_ot = (dt < r_cut).sum(1)
            n_dme = np.array([len(np.unique(tr.molid[od][dd[i] < r_cut])) for i in range(len(li))])
            n_tfsi = np.array([len(np.unique(tr.molid[ot][dt[i] < r_cut])) for i in range(len(li))])
            disp = np.linalg.norm(P[f + h] - P[f], axis=1)
            for i in range(len(li)):
                X.append([n_od[i], n_ot[i], n_dme[i], n_tfsi[i], ratio]); y.append(disp[i]); group.append(lab)
    return np.array(X, float), np.array(y), np.array(group)


def part_b(workdir=WORK) -> dict | None:
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import KFold
    from sklearn.metrics import r2_score, mean_absolute_error
    X, y, g = build_dataset(workdir)
    if len(np.unique(g)) < 3:
        return None
    rng = np.random.default_rng(0)
    model = lambda: GradientBoostingRegressor(n_estimators=200, max_depth=3, learning_rate=0.05, random_state=0)
    # (i) random split
    r2_rand, mae_rand = [], []
    for tr_i, te_i in KFold(5, shuffle=True, random_state=0).split(X):
        m = model().fit(X[tr_i], y[tr_i]); p = m.predict(X[te_i])
        r2_rand.append(r2_score(y[te_i], p)); mae_rand.append(mean_absolute_error(y[te_i], p))
    # (ii) leave-one-concentration-out
    loco = {}
    for lab in np.unique(g):
        te = g == lab
        m = model().fit(X[~te], y[~te]); p = m.predict(X[te])
        base = np.full(te.sum(), y[~te].mean())
        loco[lab] = {"r2": float(r2_score(y[te], p)), "mae": float(mean_absolute_error(y[te], p)),
                     "mae_baseline_train_mean": float(mean_absolute_error(y[te], base)),
                     "y_mean_true": float(y[te].mean()), "y_mean_pred": float(p.mean())}
    # (ii') same but WITHOUT the concentration feature: what does local structure alone carry?
    loco_nofeat = {}
    for lab in np.unique(g):
        te = g == lab
        m = model().fit(X[~te][:, :4], y[~te]); p = m.predict(X[te][:, :4])
        loco_nofeat[lab] = float(r2_score(y[te], p))
    res = {"n_rows": int(len(y)), "groups": {l: int((g == l).sum()) for l in np.unique(g)},
           "features": ["n_O(DME)", "n_O(TFSI)", "n_DME_molecules", "n_TFSI_molecules", "DME:Li ratio"],
           "target": "|Δr(Li)| over 100 ps (Å)",
           "random_split": {"r2_mean": float(np.mean(r2_rand)), "mae_mean": float(np.mean(mae_rand))},
           "leave_one_concentration_out": loco, "loco_without_ratio_feature_r2": loco_nofeat,
           "loco_r2_mean": float(np.mean([v["r2"] for v in loco.values()])),
           "note": "random split mixes frames of the same trajectory into train and test; LOCO is the split that "
                   "matches 'predict a new formulation'. Negative R² = worse than predicting the mean."}
    json.dump(res, open(os.path.join(OUT, "e8_property_split.json"), "w"), indent=1)
    return res


def report(force=False) -> dict:
    print("E8  AI4Electrolytes — potential model vs property model, validated differently")
    a = part_a(force)
    print(f"  A. GFN2-xTB vs {a['reference']} on {a['n_shells']} Li+ shells cut from the liquid "
          f"({a['compositions']}):")
    print(f"     force RMSE {a['force_rmse_kcal_mol_A']:.2f} kcal/mol/Å (reference RMS force {a['force_rms_reference']:.2f}), "
          f"on Li {a['force_rmse_on_Li_kcal_mol_A']:.2f}; force Pearson r = {a['force_pearson']:.3f}; "
          f"relative-energy RMSE {a['relative_energy_rmse_kcal']} kcal/mol over {a['n_relative_energy_pairs']} same-composition pairs")
    if "pair_complexes" in a:
        for k, v in a["pair_complexes"].items():
            print(f"     1:1 {k:16s} ΔE_vertical DFT {v['dft_vertical']:7.2f}  xtb {v['xtb_vertical']:7.2f} kcal/mol")
    b = part_b()
    if b:
        print(f"  B. Li displacement over 100 ps from shell composition, {b['n_rows']} rows, {len(b['groups'])} concentrations:")
        print(f"     random 5-fold split : R² = {b['random_split']['r2_mean']:.3f}, MAE {b['random_split']['mae_mean']:.2f} Å")
        for lab, v in b["leave_one_concentration_out"].items():
            print(f"     hold out {lab:11s}: R² = {v['r2']:+.3f}, MAE {v['mae']:.2f} Å (baseline train-mean MAE {v['mae_baseline_train_mean']:.2f}); "
                  f"without ratio feature R² = {b['loco_without_ratio_feature_r2'][lab]:+.3f}")
        print(f"     {b['note']}")
    _plot(a, b)
    return {"A": a, "B": b}


def _plot(a, b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(a["scatter"]["F_dft"], a["scatter"]["F_xtb"], ".", ms=2, alpha=0.5)
    lim = np.max(np.abs(a["scatter"]["F_dft"])) * 1.05
    ax[0].plot([-lim, lim], [-lim, lim], "k-", lw=0.5); ax[0].set_xlabel("F DFT (kcal/mol/Å)"); ax[0].set_ylabel("F GFN2-xTB")
    ax[0].set_title(f"A  forces on Li+ shells: RMSE {a['force_rmse_kcal_mol_A']:.1f}, r={a['force_pearson']:.2f}", fontsize=9)
    if a["scatter"]["E_rel"]:
        e = np.array(a["scatter"]["E_rel"])
        ax[1].plot(e[:, 0], e[:, 1], "o"); l = np.abs(e).max() * 1.1
        ax[1].plot([-l, l], [-l, l], "k-", lw=0.5)
        ax[1].set_xlabel("ΔE DFT (kcal/mol, same composition)"); ax[1].set_ylabel("ΔE xtb")
        ax[1].set_title(f"A  relative energies: RMSE {a['relative_energy_rmse_kcal']:.1f} kcal/mol", fontsize=9)
    if b:
        labs = list(b["leave_one_concentration_out"])
        r2 = [b["leave_one_concentration_out"][l]["r2"] for l in labs]
        ax[2].bar(np.arange(len(labs)), r2, color="C1", label="leave-one-concentration-out")
        ax[2].axhline(b["random_split"]["r2_mean"], color="C0", ls="--", label=f"random split R²={b['random_split']['r2_mean']:.2f}")
        ax[2].axhline(0, color="k", lw=0.5)
        ax[2].set_xticks(np.arange(len(labs))); ax[2].set_xticklabels(labs, fontsize=7, rotation=20)
        ax[2].set_ylabel("R² (test)"); ax[2].legend(fontsize=7); ax[2].set_ylim(min(-1, min(r2) - 0.1), 1)
        ax[2].set_title("B  property model: the split decides the score", fontsize=9)
    fig.suptitle("E8  potential validation (forces/energies on the target environment) vs property-model validation (split design)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e8_ml.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(force="--force" in sys.argv)
