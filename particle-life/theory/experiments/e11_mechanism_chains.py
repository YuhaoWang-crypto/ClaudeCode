r"""
E11 - From pattern to function: what each morphology lets a cell DO.
=====================================================================

E10 said which pattern a protein group makes.  This one closes the chain

    valence + Kd + concentration -> model quantity -> morphology -> capability

with a computed number at every arrow.

Panels
------
1. Ultrasensitivity from percolation alone: the gel fraction is exactly zero
   below threshold and near-vertical above it, giving Hill coefficients of 4-6
   out of bonds that have **no cooperativity whatsoever**.
2. Concentration buffering: above threshold, added protein goes into the network,
   not into free monomer, so the free pool is clamped.
3. Exclusion-driven bistability: phosphorylation builds the network, the network
   excludes the phosphatase, and the loop closes into a hysteretic switch.
4. How much exclusion is needed: the bistable region in (phosphatase access,
   kinase:phosphatase) space.
5. Polymer length as valence: the minimum dsDNA length that switches cGAS on,
   and how it moves with cGAS abundance.
6. Stoichiometry, not affinity, limits LAT: more SOS1 builds the network, more
   Grb2 destroys it.
7. Co-concentration catalysis, with the volume and viscosity penalties that stop
   the naive K_E*K_S figure from being the answer.
8. Where the real atlas modules sit on the steepness axis.
"""

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt

from plife.bio import function as F, pathways as PW, sites as SI

BOOST = float(PW.membrane_boost())


def cross(kd_uM, nA, nB, d=5.0):
    return SI.SiteMixture.from_valence(["A", "B"], [d, d], [nA, nB],
                                       np.array([[np.inf, kd_uM],
                                                 [kd_uM, np.inf]]) * 1e-6)


def main():
    out = {}
    fig, ax = plt.subplots(2, 4, figsize=(19.5, 9.2))

    # ------------------------------------------------- 1 ultrasensitivity
    a = ax[0, 0]
    out["hill"] = {}
    for n, col in zip((3, 4, 6, 8), ACC):
        m = cross(10.0, n, n)
        sc = np.geomspace(0.02, 300, 320)
        dr = F.dose_response(m, [1.0, 1.0], sc)
        c = dr["conc"][:, 0]
        nH = F.hill_coefficient(c, dr["gel"][:, 0], at=0.5)
        cstar = float(PW.gel_site_conc(10.0, n, n)) / n
        a.plot(c / cstar, dr["gel"][:, 0], color=col, lw=1.9,
               label=f"n = {n},  $n_H$ = {nH:.1f}")
        out["hill"][n] = dict(nH=nH, c_star_uM=cstar)
    x = np.geomspace(0.05, 20, 200)
    a.plot(x, x / (1 + x), "--", color=ACC[6], lw=1.4, label="hyperbolic,  $n_H$ = 1")
    a.plot(x, x**2.8 / (1 + x**2.8), ":", color=ACC[6], lw=1.4,
           label="haemoglobin,  $n_H$ = 2.8")
    a.axvline(1.0, color=ACC[5], lw=1.0, alpha=0.6)
    a.set_xscale("log"); a.set_xlim(0.05, 20); a.set_ylim(-0.03, 1.03)
    a.set_xlabel("concentration / percolation threshold")
    a.set_ylabel("fraction in the network  $g$")
    a.set_title("1  ultrasensitivity with NO cooperativity\n"
                "one Kd, no allostery - the sharpness\nis the transition itself", fontsize=9.5)
    a.legend(fontsize=7, loc="upper left")

    # ------------------------------------------------------- 2 buffering
    a = ax[0, 1]
    m = cross(10.0, 4, 4)
    cstar = float(PW.gel_site_conc(10.0, 4, 4)) / 4
    sc = np.geomspace(2e-3, 80, 340)
    dr = F.dose_response(m, [1.0, 1.0], sc)
    tot, free = dr["conc"][:, 0], dr["free"][:, 0]
    a.loglog(tot, free, color=ACC[0], lw=2, label="free (all sites unbonded)")
    a.loglog(tot, tot, "--", color=ACC[6], lw=1.3, label="no buffering (free = total)")
    a.axvline(cstar, color=ACC[5], lw=1.4)
    a.text(0.62, 0.05, "percolation\nthreshold", transform=a.transAxes,
           color=ACC[5], fontsize=7.5, ha="right")
    sl = F.buffering_quality(tot, free)
    lo = sl[tot < 0.05 * cstar].mean()
    hi = sl[tot > 6 * cstar].mean()
    a.text(0.03, 0.95, f"$d\\ln c_{{free}}/d\\ln c_{{tot}}$\n"
                       f"  dilute:  {lo:+.2f}\n  above:   {hi:+.2f}",
           transform=a.transAxes, fontsize=8, va="top", color=ACC[3])
    a.set_xlabel("total protein  uM"); a.set_ylabel("free protein  uM")
    a.set_title("2  the free pool is clamped\nexpression noise stops propagating\n"
                "once the network forms", fontsize=9.5)
    a.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.30, 0.99))
    out["buffering"] = dict(slope_dilute=float(lo), slope_above=float(hi),
                            c_star_uM=cstar)

    # --------------------------------------------------- 3 bistability
    N_MAX, N_ADAPT, C, KD, ACCESS = 10.0, 4.0, 10.0, 10.0, 0.05

    def gel_of_p(p):
        n = max(N_MAX * p, 1e-6)
        sm = SI.SiteMixture.from_valence(["S", "A"], [5.0, 5.0], [n, N_ADAPT],
                                         np.array([[np.inf, KD], [KD, np.inf]]) * 1e-6)
        return F.gel_fraction([C, C], sm)[0]

    Rs = np.geomspace(0.02, 20, 220)
    bi = F.phospho_bistability(gel_of_p, Rs, access_in_gel=ACCESS)
    a = ax[0, 2]
    a2 = a.twinx()
    a2.plot(bi["p_grid"], bi["gel"], color=ACC[2], lw=1.4, alpha=0.75)
    a2.plot(bi["p_grid"], bi["access"], color=ACC[3], lw=1.4, alpha=0.75)
    a2.set_ylim(-0.03, 1.03); a2.set_ylabel("$g(p)$ and $\\alpha(p)$", fontsize=8)
    a2.tick_params(labelsize=7); a2.grid(alpha=0)
    a2.text(0.55, 0.90, "gel fraction $g$", color=ACC[2], fontsize=7.5,
            transform=a2.transAxes)
    a2.text(0.55, 0.80, "phosphatase access $\\alpha$", color=ACC[3], fontsize=7.5,
            transform=a2.transAxes)
    a.plot(bi["p_grid"], bi["rhs"], color=ACC[0], lw=2.2, zorder=5)
    st = [o for o in bi["states"] if o["bistable"]]
    if st:
        mid = st[len(st) // 2]
        Rmid = mid["R"]
        a.axhline(Rmid, color=ACC[6], lw=1.5, ls="--", zorder=4)
        for r, s in zip(mid["roots"], mid["stable"]):
            a.plot([r], [Rmid], "o" if s else "x", color=ACC[6], ms=9, mew=2.2, zorder=6)
        a.text(0.03, 0.10, f"$R = p\\,\\alpha(p)/(1-p)$  (green)\n"
               f"at $R$ = {Rmid:.2f}: three roots\n(o stable, x unstable)",
               transform=a.transAxes, color=ACC[6], fontsize=7.5, va="bottom")
    a.set_yscale("log"); a.set_ylim(1e-3, 30); a.set_xlim(0, 1)
    a.set_zorder(a2.get_zorder() + 1); a.patch.set_visible(False)
    a.set_xlabel("phosphorylation fraction  $p$")
    a.set_ylabel("kinase : phosphatase  $R$")
    a.set_title("3  bistability EMERGES\nphospho builds the net, the net\n"
                "excludes the phosphatase", fontsize=9.5)
    out["bistability"] = dict(n_bistable=len(st),
                              R_range=[float(min(o["R"] for o in st)),
                                       float(max(o["R"] for o in st))] if st else None,
                              p_off=float(min(o["off"] for o in st)) if st else None,
                              p_on=float(max(o["on"] for o in st)) if st else None,
                              access=ACCESS)

    # ------------------------------------ 4 how much exclusion is needed
    a = ax[0, 3]
    accs = np.geomspace(0.005, 1.0, 26)
    lo_hi = []
    for acc in accs:
        b = F.phospho_bistability(gel_of_p, Rs, access_in_gel=float(acc), n_p=161)
        s = [o["R"] for o in b["states"] if o["bistable"]]
        lo_hi.append((min(s), max(s)) if s else (np.nan, np.nan))
    lo_hi = np.array(lo_hi)
    ok = np.isfinite(lo_hi[:, 0])
    a.fill_between(accs[ok], lo_hi[ok, 0], lo_hi[ok, 1], color=ACC[0], alpha=0.35)
    a.plot(accs[ok], lo_hi[ok, 0], color=ACC[0], lw=1.8)
    a.plot(accs[ok], lo_hi[ok, 1], color=ACC[0], lw=1.8)
    for k, lab, col, ha in ((0.0025, "8 nm probe", ACC[1], 0.06),
                            (0.059, "3 nm probe", ACC[3], 0.62)):
        a.axvline(k, color=col, lw=1.4, ls="--")
        a.text(ha, 0.03, lab, color=col, fontsize=7.5, transform=a.transAxes,
               va="bottom")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("phosphatase access inside the network  $a$")
    a.set_ylabel("kinase : phosphatase  $R$")
    a.set_title("4  how much exclusion it takes\nshaded = bistable.  measured\n"
                "size-exclusion values marked", fontsize=9.5)
    out["bistable_region"] = dict(access=accs.tolist(), R_lo=lo_hi[:, 0].tolist(),
                                  R_hi=lo_hi[:, 1].tolist())

    # ------------------------------------------- 5 DNA length as valence
    a = ax[1, 0]

    def build_dna(L_bp):
        return SI.SiteMixture(
            names=["cGAS", "dsDNA"], diameters=[5.2, 2.0],
            sites=[SI.Site("cGAS", "DNA", 2.0),
                   SI.Site("dsDNA", "turn", max(L_bp / 20.0, 1.0))],
            kd=np.array([[np.inf, 0.5], [0.5, np.inf]]) * 1e-6)

    L = np.geomspace(20, 6000, 130)
    out["dna_length"] = {}
    for cg, col in zip((0.5, 0.2, 0.05, 0.02), ACC):
        r = F.min_polymer_length(build_dna, L, [cg, cg / 10])
        a.plot(L, r["branching"], color=col, lw=1.8, label=f"[cGAS] = {cg:g} uM")
        if np.isfinite(r["L_min"]):
            a.plot([r["L_min"]], [1.0], "o", color=col, ms=6)
        out["dna_length"][cg] = r["L_min"]
    a.axhline(1.0, color=ACC[6], lw=1.8)
    a.text(25, 1.12, "percolation", color=ACC[6], fontsize=8)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("dsDNA length  bp"); a.set_ylabel("branching eigenvalue")
    a.set_title("5  ligand length IS the valence\nself/non-self discrimination\n"
                "by DNA length, from arithmetic", fontsize=9.5)
    a.legend(fontsize=7, loc="lower right")

    # ---------------------------------------- 6 LAT is stoichiometry-limited
    a = ax[1, 1]
    sysL = PW.get("lat")
    lat = PW.to_mixture(sysL)
    c0 = sysL.conc("mid") * BOOST
    mult = np.geomspace(0.05, 40, 90)
    for idx, lab, col in ((2, "SOS1", ACC[0]), (1, "Grb2", ACC[1]), (0, "pLAT", ACC[2])):
        g = []
        for f in mult:
            c = c0.copy(); c[idx] *= f
            g.append(F.gel_fraction(c, lat)[0])
        a.plot(mult, g, color=col, lw=2, label=f"vary {lab}")
    a.axvline(1.0, color=ACC[6], lw=1.2, ls="--")
    a.set_xscale("log")
    a.set_xlabel("fold change of one component (others at membrane level)")
    a.set_ylabel("fraction of pLAT in the network")
    a.set_title("6  which knob actually works\nmore SOS1 builds it,\n"
                "more Grb2 DESTROYS it", fontsize=9.5)
    a.legend(fontsize=7.5, loc="upper right")
    out["lat_knobs"] = dict(
        mult=mult.tolist(),
        gel_even=float(F.gel_fraction(c0, lat)[0]),
        gel_sos6=float(F.gel_fraction(c0 * np.array([1, 1, 6.0]), lat)[0]),
        gel_grb6=float(F.gel_fraction(c0 * np.array([1, 6.0, 1]), lat)[0]))

    # ------------------------------------------------ 7 rate enhancement
    a = ax[1, 2]
    phi = np.geomspace(1e-4, 0.3, 200)
    for K, col in zip((10.0, 100.0, 1000.0), ACC):
        for eta, ls, alp in ((1.0, "-", 1.0), (100.0, "--", 0.75)):
            w, _ = F.rate_enhancement(K, K, phi, viscosity_penalty=eta)
            a.loglog(phi, w, ls, color=col, lw=1.7, alpha=alp,
                     label=f"$K_E=K_S=${K:g}" if eta == 1 else None)
    a.axhline(1.0, color=ACC[6], lw=1.2)
    a.text(0.97, 0.05, "solid: no viscosity penalty\ndashed: 100x slower interior",
           transform=a.transAxes, fontsize=7.5, ha="right", color=ACC[6])
    a.set_xlabel("condensate volume fraction  $\\phi$")
    a.set_ylabel("whole-cell flux, fold change")
    a.set_title("7  co-concentration catalysis\nthe $K_EK_S$ headline is local;\n"
                "volume and viscosity take it back", fontsize=9.5)
    a.legend(fontsize=7, loc="upper left")

    # ----------------------------------- 8 steepness of the real modules
    a = ax[1, 3]
    rows = []
    for s in PW.ATLAS:
        sm = PW.to_mixture(s)
        base = s.conc("mid")
        sc = np.geomspace(1e-3, 3000, 300)
        dr = F.dose_response(sm, base, sc)
        g = dr["gel"].max(axis=1)
        if g.max() < 0.05:
            continue
        nH = F.hill_coefficient(sc, g / max(g.max(), 1e-9), at=0.5)
        on = sc[np.argmax(g > 0)] if (g > 0).any() else np.nan
        rows.append((s.key, nH, on, g.max()))
    rows.sort(key=lambda r: -(r[1] if np.isfinite(r[1]) else -1))
    y = np.arange(len(rows))
    a.barh(y, [r[1] for r in rows], color=ACC[2], alpha=0.9)
    a.axvline(1.0, color=ACC[6], lw=1.6)
    a.axvline(2.8, color=ACC[5], lw=1.4, ls="--")
    a.text(1.05, -0.75, "hyperbolic", color=ACC[6], fontsize=7.5, va="bottom")
    a.text(2.9, -0.75, "haemoglobin", color=ACC[5], fontsize=7.5, va="bottom")
    a.set_ylim(-1.0, len(rows) - 0.4)
    a.set_yticks(y, [r[0] for r in rows], fontsize=8)
    a.set_xlabel("effective Hill coefficient of the switch")
    a.set_title("8  every real module is a sharp switch\n"
                "and none of them uses cooperative binding", fontsize=9.5)
    a.grid(axis="y", alpha=0)
    out["module_steepness"] = [dict(key=r[0], nH=r[1], switch_at=r[2], gel_max=r[3])
                               for r in rows]

    fig.suptitle("E11  Pattern -> capability: ultrasensitivity, buffering, bistability, "
                 "length sensing - all from connectivity", fontsize=12.5)
    fig.tight_layout()
    save(fig, "e11_mechanism_chains.png")
    dump(out, "e11_mechanism_chains.json")

    print("\n  capability                          computed value")
    print(f"  ultrasensitivity (n=4, Kd 10 uM)    n_H = {out['hill'][4]['nH']:.2f}"
          f"   (hyperbolic = 1, haemoglobin = 2.8)")
    print(f"  buffering above threshold           dln free/dln total = "
          f"{out['buffering']['slope_above']:+.2f}")
    if out["bistability"]["R_range"]:
        r = out["bistability"]
        print(f"  bistable window                     R = {r['R_range'][0]:.2f}-"
              f"{r['R_range'][1]:.2f}  ({r['R_range'][1] / r['R_range'][0]:.1f}x hysteresis), "
              f"p = {r['p_off']:.2f} vs {r['p_on']:.2f}")
    print(f"  cGAS length threshold               "
          + ",  ".join(f"{k:g} uM -> {v:.0f} bp" for k, v in out["dna_length"].items()))
    k = out["lat_knobs"]
    print(f"  LAT: 6x SOS1 vs 6x Grb2             gel {k['gel_even']:.3f} -> "
          f"{k['gel_sos6']:.3f} (SOS1) vs {k['gel_grb6']:.3f} (Grb2)")


if __name__ == "__main__":
    main()
