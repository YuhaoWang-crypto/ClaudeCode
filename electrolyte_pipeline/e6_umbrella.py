"""
E6b — umbrella sampling of the Li…O(DME) coordinate in the liquid   [digest p11, "自由能的来源要说清楚"]

The E6 liquid PMF, W = −kT ln g(r), is only defined where the unbiased run
sampled r.  Here the same coordinate is driven with harmonic windows so the
barrier region is sampled on purpose, and the PMF is recovered with MBAR.

Reported, as the digest asks: bias form and force constant, window centres,
samples per window, overlap between neighbouring windows, and the MBAR
statistical uncertainty.  The coordinate is the distance between ONE Li+ and
ONE O atom of ONE DME molecule (chosen from a coordinated pair in final.pdb),
so the PMF is that of *that* pair separating while everything else — other
DME, TFSI, the second O of the same DME — relaxes freely.  The g(r)-based PMF
averages over all Li–O(DME) pairs; the two are compared, not equated.

Free energies are relative to the first minimum; the large-r plateau is not
reached (windows stop at 7 Å) so the absolute solvation free energy is NOT
claimed.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, DynamicsSpec

KT_KCAL = 0.0019872041 * 298.15
K_KJ_MOL_NM2 = 4000.0                 # ~ 9.6 kcal/mol/Å² : σ_window ≈ 0.25 Å at 298 K
CENTRES_A = np.arange(1.8, 7.01, 0.25)


def _pick_pair(sysdir):
    """One Li and one DME O currently within 2.5 Å of it (a coordinated pair)."""
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array
    u = mda.Universe(os.path.join(sysdir, "final.pdb"))
    atoms = json.load(open(os.path.join(sysdir, "atoms.json")))
    sp = np.array([a["species"] for a in atoms]); el = np.array([a["element"] for a in atoms])
    li = np.where(sp == "Li")[0]; od = np.where((sp == "DME") & (el == "O"))[0]
    d = distance_array(u.atoms.positions[li], u.atoms.positions[od], box=u.dimensions)
    i, j = np.unravel_index(np.argmin(np.abs(d - 2.1)), d.shape)
    return int(li[i]), int(od[j])


def run(label: str = "LiTFSI_r10", workdir: str = WORK, ps_per_window: float = 400.0,
        equil_ps: float = 100.0, sample_ps: float = 0.2) -> dict:
    import openmm
    from openmm import app, unit
    from electrolyte_pipeline.e1_run import _platform
    dyn = DynamicsSpec()
    sysdir = os.path.join(workdir, label)
    out = os.path.join(sysdir, "umbrella"); os.makedirs(out, exist_ok=True)
    li, o = _pick_pair(sysdir)
    system = openmm.XmlSerializer.deserialize(open(os.path.join(sysdir, "system.xml")).read())
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i), openmm.MonteCarloBarostat):
            system.removeForce(i)
    bias = openmm.CustomBondForce("0.5*k*(r-r0)^2")
    bias.addGlobalParameter("k", K_KJ_MOL_NM2); bias.addGlobalParameter("r0", 0.21)
    bias.addBond(li, o, []); bias.setUsesPeriodicBoundaryConditions(True)
    system.addForce(bias)
    pdb = app.PDBFile(os.path.join(sysdir, "final.pdb"))
    integ = openmm.LangevinMiddleIntegrator(dyn.temperature_K * unit.kelvin, dyn.friction_per_ps / unit.picosecond,
                                            dyn.timestep_fs * unit.femtosecond)
    plat, props = _platform()
    sim = app.Simulation(pdb.topology, system, integ, plat, props)
    sim.context.setPositions(pdb.positions); sim.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
    sim.context.setVelocitiesToTemperature(dyn.temperature_K * unit.kelvin, 11)
    spp = int(1000 / dyn.timestep_fs)
    samples = {}
    t0 = time.time()
    # walk the windows outward from the coordinated state so each starts from the previous one
    for c in CENTRES_A:
        sim.context.setParameter("r0", c / 10.0)
        sim.step(int(equil_ps * spp))
        rs = []
        for _ in range(int((ps_per_window - equil_ps) / sample_ps)):
            sim.step(int(sample_ps * spp))
            st = sim.context.getState(getPositions=True, enforcePeriodicBox=False)
            p = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
            box = np.array([st.getPeriodicBoxVectors()[i][i].value_in_unit(unit.angstrom) for i in range(3)])
            d = p[o] - p[li]; d -= box * np.round(d / box)
            rs.append(float(np.linalg.norm(d)))
        samples[float(c)] = rs
        print(f"[E6b] window r0={c:.2f} Å: <r>={np.mean(rs):.2f} sd {np.std(rs):.2f}  ({time.time()-t0:.0f} s)")
    rec = {"label": label, "pair": {"Li_index": li, "O_index": o}, "bias": "0.5 k (r-r0)^2",
           "k_kJ_mol_nm2": K_KJ_MOL_NM2, "centres_A": CENTRES_A.tolist(), "ps_per_window": ps_per_window,
           "equil_ps": equil_ps, "sample_ps": sample_ps, "platform": plat.getName(),
           "samples": samples, "wall_s": round(time.time() - t0, 1)}
    json.dump(rec, open(os.path.join(out, "umbrella_samples.json"), "w"))
    return {k: v for k, v in rec.items() if k != "samples"}


def analyse(label: str = "LiTFSI_r10", workdir: str = WORK, nbins: int = 52) -> dict | None:
    from pymbar import MBAR, timeseries
    path = os.path.join(workdir, label, "umbrella", "umbrella_samples.json")
    if not os.path.exists(path):
        return None
    rec = json.load(open(path))
    centres = np.array(rec["centres_A"]); k = rec["k_kJ_mol_nm2"] / 100.0 / 4.184   # kcal/mol/Å²
    K = len(centres)
    # subsample each window by its statistical inefficiency
    r_k, N_k, g_k = [], [], []
    for c in centres:
        x = np.array(rec["samples"][str(float(c))]) if str(float(c)) in rec["samples"] else np.array(rec["samples"][f"{c}"])
        g = timeseries.statistical_inefficiency(x)
        idx = timeseries.subsample_correlated_data(x, g=g)
        r_k.append(x[idx]); N_k.append(len(idx)); g_k.append(float(g))
    N_k = np.array(N_k); r_all = np.concatenate(r_k)
    u_kn = np.array([0.5 * k * (r_all - c) ** 2 / KT_KCAL for c in centres])          # reduced bias energies
    edges = np.linspace(1.6, 7.2, nbins + 1); mid = 0.5 * (edges[1:] + edges[:-1])
    bin_n = np.digitize(r_all, edges) - 1

    def pmf_from(u, N, bins, rng=None):
        """MBAR unbiased weights -> histogram PMF (kcal/mol, min = 0)."""
        if rng is not None:                                   # bootstrap: resample within each window
            idx, start = [], 0
            for n in N:
                idx.append(start + rng.integers(0, n, n)); start += n
            idx = np.concatenate(idx); u = u[:, idx]; bins = bins[idx]
        m = MBAR(u, N)
        f_k = m.f_k
        log_w = -np.logaddexp.reduce(np.log(N)[:, None] + f_k[:, None] - u, axis=0)   # unbiased (u_0 = 0)
        w = np.exp(log_w - log_w.max())
        hist = np.zeros(len(mid))
        okb = (bins >= 0) & (bins < len(mid))
        np.add.at(hist, bins[okb], w[okb])
        with np.errstate(divide="ignore"):
            Wb = -KT_KCAL * np.log(hist / hist.sum())
        return Wb - np.nanmin(Wb[np.isfinite(Wb)])

    W = pmf_from(u_kn, N_k, bin_n)
    rng = np.random.default_rng(0)
    boots = np.array([pmf_from(u_kn, N_k, bin_n, rng) for _ in range(20)])
    dW = np.nanstd(np.where(np.isfinite(boots), boots, np.nan), axis=0)
    W = np.where(np.isfinite(W), W, np.nan)
    # Jacobian: the 1-D distance PMF as sampled is W(r); to compare with -kT ln g(r) add 2 kT ln r
    W_g = W + 2 * KT_KCAL * np.log(mid)   # convert W(r) -> radial-density PMF comparable to -kT ln g(r)
    W_g -= np.nanmin(W_g)
    # overlap: fraction of neighbouring-window sample ranges that overlap
    lo = np.array([np.percentile(x, 5) for x in r_k]); hi = np.array([np.percentile(x, 95) for x in r_k])
    overlap = float(np.mean(hi[:-1] > lo[1:]))
    i_min = int(np.nanargmin(np.where(mid < 3.0, W_g, np.nan)))
    sel = (mid > mid[i_min]) & (mid < 5.0)
    i_bar = int(np.where(sel)[0][np.nanargmax(W_g[sel])])
    res = {"label": label, "r_A": mid.tolist(), "W_r_kcal": W.tolist(), "dW_kcal": dW.tolist(),
           "W_radial_kcal": W_g.tolist(),
           "barrier_kcal": float(W_g[i_bar] - W_g[i_min]), "dbarrier_kcal": float(np.sqrt(dW[i_bar] ** 2 + dW[i_min] ** 2)),
           "r_min_A": float(mid[i_min]), "r_barrier_A": float(mid[i_bar]),
           "windows": K, "k_kcal_mol_A2": k, "samples_per_window_uncorrelated": N_k.tolist(),
           "statistical_inefficiency": g_k, "neighbour_overlap_fraction": overlap,
           "note": ("one specific Li–O(DME) pair; W(r) is the 1-D distance PMF, W_radial adds 2kT ln r for "
                    "comparison with -kT ln g(r); plateau not reached (windows end at 7 Å)")}
    json.dump(res, open(os.path.join(workdir, label, "umbrella", "e6b_pmf.json"), "w"))
    return res


def report(workdir: str = WORK):
    r = analyse(workdir=workdir)
    if r is None:
        print("E6b umbrella: no samples (run `modal run electrolyte_pipeline/modal_run.py --stage umbrella`)")
        return None
    print(f"E6b umbrella sampling, {r['windows']} windows, k = {r['k_kcal_mol_A2']:.1f} kcal/mol/Å², "
          f"neighbour overlap {r['neighbour_overlap_fraction']*100:.0f} %, uncorrelated samples/window "
          f"{min(r['samples_per_window_uncorrelated'])}–{max(r['samples_per_window_uncorrelated'])}")
    print(f"    minimum at {r['r_min_A']:.2f} Å; barrier to {r['r_barrier_A']:.2f} Å = "
          f"{r['barrier_kcal']:.2f} ± {r['dbarrier_kcal']:.2f} kcal/mol (MBAR)")
    _plot(r, workdir)
    return r


def _plot(r, workdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from electrolyte_pipeline.e6_desolvation import liquid_pmf
    fig, ax = plt.subplots(figsize=(6, 3.8))
    liq = liquid_pmf(r["label"], workdir)
    if liq:
        rr = np.array(liq["r_A"]); W = np.array([np.nan if w is None else w for w in liq["W_kcal"]])
        W -= np.nanmin(W[(rr > 1.8) & (rr < 3)])
        ax.plot(rr, W, "k-", lw=1, label="−kT ln g(r), all pairs, unbiased 10 ns")
    ax.errorbar(r["r_A"], r["W_radial_kcal"], yerr=r["dW_kcal"], fmt="o-", ms=3, color="C3",
                label=f"umbrella + MBAR, one pair ({r['windows']} windows)")
    ax.set_xlim(1.6, 7.2); ax.set_ylim(-0.5, 8); ax.set_xlabel("Li…O(DME) (Å)"); ax.set_ylabel("W (kcal/mol), minimum = 0")
    ax.set_title(f"E6b  barrier {r['barrier_kcal']:.1f} ± {r['dbarrier_kcal']:.1f} kcal/mol at {r['r_barrier_A']:.1f} Å", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e6b_umbrella.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    if "--run" in sys.argv:
        run(ps_per_window=4.0, equil_ps=2.0)
    report()
