r"""
E10 - Real signalling pathways in one pool: what pattern, under what conditions?
================================================================================

Twelve literature-parameterised modules (``plife.bio.pathways.ATLAS``) are put
through the four closed-form conditions and then dumped into a single shared
pool.

Panels
------
1. The valence gate.  Threshold site concentration against ``(n_A-1)(n_B-1)``,
   with every atlas module placed on it.  The gate is a wall, not a slope.
2. Threshold vs available material: how far each module sits from switching on,
   and how far a 2-D membrane boost moves it.
3. The stoichiometry window: percolation region in the ``(c_A, c_B)`` plane,
   showing both edges of the re-entrant boundary.
4. Miscibility: cross-Kd against the harmonic mean of the self-Kd's, with the
   coexistence calculation checking which side of the line each point lands on.
5. The pool at bulk cytosolic concentration: which modules are on, which off.
6. Mode structure of the pool -- condensation modes vs demixing modes, the same
   coherence order parameter used on the Particle Life dispersion relation.
7. Orthogonality: how tight cross-pathway binding has to get before the sorted
   pool collapses into one blob.
8. The reciprocity boundary: what an ATP cycle has to outrun to make any of this
   move rather than just sit there.
"""

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt

from plife.bio import pathways as PW, phases as P, sites as SI, thermo as T, units as U

BOOST = float(PW.membrane_boost())


# --------------------------------------------------------------------------- #
def main():
    out = {}
    fig, ax = plt.subplots(2, 4, figsize=(19.5, 9.2))

    # ---------------------------------------------------------------- 1 gate
    a = ax[0, 0]
    G = np.geomspace(0.05, 200, 400)
    for kd, col in zip((1.0, 10.0, 100.0), ACC):
        p = 1.0 / np.sqrt(G)
        S = np.where(p < 1, kd * p / np.maximum((1 - p) ** 2, 1e-30), np.nan)
        a.plot(G, S, color=col, lw=1.8, label=f"Kd = {kd:g} uM")
    a.axvline(1.0, color=ACC[6], lw=2)
    a.text(0.03, 0.52, "valence gate\n$(n_A\\!-\\!1)(n_B\\!-\\!1)=1$\n\nleft of it:\nno network\nat any Kd",
           transform=a.transAxes, color=ACC[6], fontsize=7.5, ha="left", va="center")
    a.axvspan(0.05, 1.0, color=ACC[6], alpha=0.10)

    rows = []
    for s in PW.ATLAS:
        r = PW.assess(s)
        r["need"] = PW.boost_to_percolate(s)
        rows.append(r)
        if np.isfinite(r["S_gel_uM"]) and r["G"] > 0:
            a.plot(r["G"], r["S_gel_uM"], "o", color=ACC[2], ms=6, mec=ACC[1], mew=0.8)
            a.annotate(r["key"], (r["G"], r["S_gel_uM"]), fontsize=7,
                       xytext=(4, 3), textcoords="offset points", color=ACC[1])
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("branching capacity  $(n_A-1)(n_B-1)$")
    a.set_ylabel("gel-point site concentration  $S^*$   uM")
    a.set_title("1  the valence gate\nleft of the wall: no network at any Kd", fontsize=9.5)
    a.legend(fontsize=7, loc="upper right")

    # ------------------------------------------------ 2 how far from switching
    a = ax[0, 1]
    live = [r for r in rows]
    order = np.argsort([r["branching"] for r in live])
    y = np.arange(len(order))
    lam0 = np.array([live[i]["branching"] for i in order])
    lamM = np.array([PW.assess(PW.get(live[i]["key"]), local_boost=BOOST)["branching"]
                     for i in order])
    a.barh(y, lam0, color=ACC[2], alpha=0.85, label="bulk cytosolic")
    a.barh(y, np.maximum(lamM - lam0, 0), left=lam0, color=ACC[0], alpha=0.55,
           label=f"+ {BOOST:.0f}x local (membrane)")
    a.axvline(1.0, color=ACC[6], lw=2)
    a.text(1.0, len(order) - 0.2, " percolation", color=ACC[6], fontsize=8, va="top")
    a.set_yticks(y, [live[i]["key"] for i in order], fontsize=8)
    a.set_xscale("log"); a.set_xlim(1e-3, 30); a.set_ylim(-0.9, len(order) - 0.2)
    a.set_xlabel("branching eigenvalue  $\\varrho(M)$")
    a.set_title("2  distance to the switch\nmost modules are OFF until localised",
                fontsize=9.5)
    a.legend(fontsize=7, loc="lower right")
    a.grid(axis="y", alpha=0)

    # ----------------------------------------------- 3 stoichiometry window
    a = ax[0, 2]
    nA, nB, kd = 4.0, 4.0, 10.0
    m = SI.SiteMixture.from_valence(["A", "B"], [5.0, 5.0], [nA, nB],
                                    np.array([[np.inf, kd], [kd, np.inf]]) * 1e-6)
    cs = np.geomspace(0.3, 3e4, 90)
    lam = np.array([[SI.percolates([x, yv], m)[1] for yv in cs] for x in cs])
    a.contourf(cs, cs, lam.T, levels=[1.0, 1e9], colors=[ACC[0]], alpha=0.35)
    a.contour(cs, cs, lam.T, levels=[1.0], colors=[ACC[0]], linewidths=1.8)
    Gw = PW.stoichiometry_window(nA, nB)
    for f, ls in ((Gw, "--"), (1.0 / Gw, "--")):
        a.plot(cs, cs * (nA / nB) / f, ls, color=ACC[1], lw=1.3)
    a.text(0.04, 0.93, f"dashed: $n_Ac_A/n_Bc_B = G^{{\\pm1}}$, $G={Gw:.0f}$",
           transform=a.transAxes, fontsize=7.5, color=ACC[1])
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(cs[0], cs[-1]); a.set_ylim(cs[0], cs[-1])
    a.set_xlabel("[A]  uM"); a.set_ylabel("[B]  uM")
    a.set_title(f"3  stoichiometry window\n$n_A=n_B={nA:.0f}$, Kd = {kd:g} uM",
                fontsize=9.5)

    # ------------------------------------------------------- 4 miscibility
    a = ax[0, 3]
    kd_self = 30.0
    cross = np.geomspace(1.0, 1000.0, 13)
    n_dense, asym = [], []
    cmax = U.density_to_conc(0.55 / ((4 / 3) * np.pi * 2.5 ** 3))
    grid = np.geomspace(0.3, cmax, 34)
    for kc in cross:
        mx = T.Mixture(names=["A", "B"], diameters=[5.0, 5.0], valence=[4.0, 4.0],
                       kd=np.array([[kd_self, kc], [kc, kd_self]]) * 1e-6)
        rho, shape = P.grid2d(grid, grid)
        res = T.coexistence(mx, rho)
        i = int(np.argmin(np.abs(grid - 0.35 * cmax)))
        sp = P._facet_split(res, i * shape[1] + i, rho)  # noqa: SLF001
        if sp is None:
            n_dense.append(1); asym.append(0.0); continue
        v = sp[0][sp[1] > 1e-3]
        v = v[np.argsort(v.sum(axis=1))]
        n_dense.append(len(v))
        d = v[-1]
        asym.append(float(abs(d[0] - d[1]) / max(d.max(), 1e-300)))
    a.plot(cross, asym, "o-", color=ACC[0], ms=5)
    a.axvline(kd_self, color=ACC[6], lw=2)
    a.text(kd_self * 1.1, 0.5, "$K_d^{AB}=\\mathrm{harm}(K_d^{self})$",
           color=ACC[6], fontsize=8, rotation=90, va="center")
    a.set_xscale("log")
    a.set_xlabel("cross-affinity  $K_d^{AB}$   uM")
    a.set_ylabel("composition asymmetry of the dense phase")
    a.set_title("4  one body or two\nleft: shared condensate.  right: immiscible",
                fontsize=9.5)
    out["miscibility"] = dict(cross_uM=cross.tolist(), n_dense=n_dense, asym=asym)

    # ------------------------------------------------------------- 5 the pool
    pool = PW.pool_analysis()
    poolM = PW.pool_analysis(local_boost=BOOST)
    a = ax[1, 0]
    keys = [r["key"] for r in rows]
    on0 = np.array([k in pool["live"] for k in keys])
    onM = np.array([k in poolM["live"] for k in keys])
    y = np.arange(len(keys))
    a.barh(y, np.where(on0, 1, 0), color=ACC[0], alpha=0.9, label="on in bulk cytosol")
    a.barh(y, np.where(onM & ~on0, 1, 0), left=np.where(on0, 1, 0), color=ACC[3],
           alpha=0.75, label="on only once localised")
    a.barh(y, np.where(~onM, 1, 0), color=ACC[6], alpha=0.35, label="never (gate shut)")
    a.set_yticks(y, keys, fontsize=8)
    a.set_xticks([]); a.set_xlim(0, 1.75)
    a.set_title(f"5  the shared pool\n{len(pool['live'])}/{len(keys)} on in bulk, "
                f"{len(poolM['live'])}/{len(keys)} once localised", fontsize=9.5)
    a.legend(fontsize=7, loc="center right")
    a.grid(alpha=0)

    # -------------------------------------------------------- 6 mode structure
    a = ax[1, 1]
    for tag, p, col, mk in (("bulk", pool, ACC[2], "o"), (f"{BOOST:.0f}x", poolM, ACC[0], "s")):
        if p["modes"]:
            e = [-m["eig"] for m in p["modes"]]
            c = [m["coherence"] for m in p["modes"]]
            a.scatter(e, c, s=46, marker=mk, color=col, alpha=0.85, label=tag)
    a.axhline(0.7, color=ACC[6], ls="--", lw=1.2)
    a.axhline(0.3, color=ACC[6], ls="--", lw=1.2)
    a.text(0.02, 0.86, "condensation (all species together)", transform=a.transAxes,
           fontsize=7.5, color=ACC[6])
    a.text(0.02, 0.06, "demixing (sorts into separate bodies)", transform=a.transAxes,
           fontsize=7.5, color=ACC[6])
    a.set_xscale("log"); a.set_ylim(-0.03, 1.03)
    a.set_xlabel("growth rate of the unstable mode  $-\\lambda$")
    a.set_ylabel("mode coherence  $|\\sum v_a| / \\sum |v_a|$")
    a.set_title("6  what the pool's unstable modes do\nsame order parameter as the "
                "Particle Life\ndispersion analysis", fontsize=9.5)
    a.legend(fontsize=7)

    # ------------------------------------------------------- 7 orthogonality
    a = ax[1, 2]
    kds = np.geomspace(0.02, 300, 16)
    pairs = {("psd", "spop"): None, ("spop", "nucleolus"): None,
             ("psd", "sg"): None, ("sg", "cgas"): None}
    nb = []
    for kd_c in kds:
        nb.append(PW.pool_analysis(cross={k: kd_c for k in pairs},
                                   modes=False)["n_bodies"])
    a.step(kds, nb, where="mid", color=ACC[0], lw=2)
    a.fill_between(kds, 0, nb, step="mid", color=ACC[0], alpha=0.2)
    harm = 2.0 / (1 / 1.0 + 1 / 10.0)
    a.axvline(harm, color=ACC[6], lw=2)
    a.text(harm * 1.15, 1.4, "$\\mathrm{harm}(K_d^{self})$", color=ACC[6], fontsize=8,
           rotation=90)
    a.set_xscale("log"); a.set_ylim(0, max(nb) + 0.6)
    a.set_xlabel("cross-pathway $K_d$   uM")
    a.set_ylabel("distinct bodies in the pool")
    a.set_title("7  orthogonality is the condition\ntighten cross-talk and the sorted\n"
                "pool collapses into one blob", fontsize=9.5)
    out["orthogonality"] = dict(cross_kd_uM=kds.tolist(), n_bodies=nb, harm_uM=harm)

    # --------------------------------------------------------- 8 reciprocity
    a = ax[1, 3]
    kd_range = np.geomspace(1e-3, 1e3, 200)          # uM
    k_on = 1e8 * 1e-6                                 # 1/(uM s), diffusion-limited-ish
    k_off = kd_range * k_on                           # 1/s
    a.loglog(kd_range, k_off, color=ACC[2], lw=2, label="$k_{off}$ of the bond")
    for kcat_E, col, lab in ((0.1, ACC[0], "slow cycle  $k_{cat}[E] = 0.1$/s"),
                             (10.0, ACC[3], "fast cycle  $k_{cat}[E] = 10$/s")):
        a.axhline(kcat_E, color=col, lw=1.6, ls="--", label=lab)
        kd_x = kcat_E / k_on
        a.plot([kd_x], [kcat_E], "o", color=col, ms=7)
    a.fill_between(kd_range, 1e-4, k_off, color=ACC[6], alpha=0.10)
    a.text(2e-3, 3e-4, "modification outruns the bond:\nnon-reciprocal, motion possible",
           fontsize=7.5, color=ACC[6])
    a.set_xlabel("bond $K_d$   uM"); a.set_ylabel("rate   s$^{-1}$")
    a.set_ylim(1e-4, 1e6)
    a.set_title("8  the reciprocity boundary\nno Kd makes a pattern move;\n"
                "only $k_{cat}[E] \\gtrsim k_{off}$ does", fontsize=9.5)
    a.legend(fontsize=7, loc="upper left")

    fig.suptitle("E10  Signalling-pathway modules in a shared pool: "
                 "one body per module, and what it takes to break that", fontsize=12.5)
    fig.tight_layout()
    save(fig, "e10_signaling_pool.png")

    # ------------------------------------------------------------------ dump
    out["boost"] = BOOST
    out["systems"] = [
        dict(key=r["key"], name=r["name"], family=r["family"],
             valence=list(np.asarray(r["valence"], float)),
             load_bearing=list(r["tightest_pair"]), kd_uM=r["kd_uM"],
             p_star=r["p_star"], S_gel_uM=r["S_gel_uM"], site_conc_uM=r["site_conc_uM"],
             margin=r["margin"], G=r["G"], site_ratio=r["site_ratio"],
             branching=r["branching"], percolates=r["percolates"],
             boost_needed=r["need"], membrane=r["membrane"], source=r["source"])
        for r in rows]
    out["pool_bulk"] = dict(live=pool["live"], n_bodies=pool["n_bodies"],
                            groups=pool["groups"], n_unstable=pool["n_unstable"],
                            modes=[{k: v for k, v in m.items() if k != "leaders"}
                                   for m in pool["modes"]])
    out["pool_local"] = dict(live=poolM["live"], n_bodies=poolM["n_bodies"],
                             groups=poolM["groups"], n_unstable=poolM["n_unstable"])
    dump(out, "e10_signaling_pool.json")

    # ------------------------------------------------------------- summary
    print("\n  module        load-bearing bond              Kd   S*    sites  branch  need")
    for r in rows:
        need = "never" if not np.isfinite(r["need"]) else f"{r['need']:.0f}x"
        S = "inf" if not np.isfinite(r["S_gel_uM"]) else f"{r['S_gel_uM']:.1f}"
        print(f"  {r['key']:<13} {r['tightest_pair'][0]}~{r['tightest_pair'][1]:<18}"
              f"{r['kd_uM']:7.2f}{S:>6}{r['site_conc_uM']:8.2f}{r['branching']:8.3f}{need:>7}")
    print(f"\n  pool, bulk cytosolic : {len(pool['live'])} modules on -> "
          f"{pool['n_bodies']} distinct bodies")
    print(f"  pool, localised {BOOST:.0f}x  : {len(poolM['live'])} modules on -> "
          f"{poolM['n_bodies']} distinct bodies")


if __name__ == "__main__":
    main()
