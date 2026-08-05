r"""
E09 - Protein mixtures: forward prediction and inverse inference.
=================================================================

Forward
-------
Given measured Kd's, valences, sizes and concentrations, what forms in the
cytoplasm?  Panels 1-4 map the (c_A, c_B) plane, the titration curve, the
partition coefficients, and the valence/affinity chart that says whether a pair
*can* phase separate at all.

Inverse
-------
Given observations, what can be inferred?  Panels 5-8 fit synthetic (noisy)
data back to the interaction matrix, compare the information content of binary
"does it condense" calls against measured partition coefficients, run the
sloppy-model eigen-analysis, and answer the practical question: **given [A],
what range of [B] is compatible with the observed structure?**

Everything uses the Wertheim TPT1 + BMCSL free energy of ``plife.bio.thermo``.
"""

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from plife.bio import inverse as I, phases as P, thermo as T, units as U

# --------------------------------------------------------------------------- #
#  a concrete two-protein system, specified the way a bench experiment would be
# --------------------------------------------------------------------------- #
TRUE = dict(kd_aa=300e-6, kd_ab=30e-6, kd_bb=np.inf)
MW = (60.0, 35.0)                    # kDa
VAL = (5.0, 3.0)                     # binding sites per molecule


def make_mix(kd_aa=TRUE["kd_aa"], kd_ab=TRUE["kd_ab"], kd_bb=TRUE["kd_bb"]):
    return T.Mixture(
        names=["A (scaffold)", "B (partner)"],
        diameters=[2 * float(U.radius_from_mw(MW[0])), 2 * float(U.radius_from_mw(MW[1]))],
        valence=list(VAL),
        kd=np.array([[kd_aa, kd_ab], [kd_ab, kd_bb]]),
    )


CMAP = ListedColormap(["#161b22", "#39d353", "#58a6ff", "#ff5ec4"])


def main():
    mix = make_mix()
    out = {}

    # ------------------------------------------------------------ forward
    cA = np.geomspace(0.03, 3000, 64)
    cB = np.geomspace(0.03, 3000, 64)
    pd = P.phase_diagram(mix, cA, cB)

    cB_tit = np.geomspace(0.01, 3000, 40)
    csat, cA_axis, _, _ = I.csat_curve(mix, cB_tit)

    print("system:", mix.names, " MW", MW, "kDa  valence", VAL)
    print(f"  diameters {np.round(mix.diameters,2)} nm   "
          f"Kd_AA={TRUE['kd_aa']*1e6:.0f} uM  Kd_AB={TRUE['kd_ab']*1e6:.0f} uM  Kd_BB=inf")
    print(f"\n  c_sat(A) with no partner : {csat[0]:.3f} uM")
    print(f"  minimum c_sat(A)          : {np.nanmin(csat):.3f} uM "
          f"at [B] = {cB_tit[int(np.nanargmin(csat))]:.2f} uM")
    print("  re-entrance: c_sat rises again at high [B] "
          f"({csat[-1]:.0f} uM at [B]={cB_tit[-1]:.0f} uM)")

    # Q3: given [A], which [B] gives a condensate?
    windows = {}
    for q in (0.3, 1.0, 3.0, 10.0, 100.0):
        w = P.partner_window(pd, q, want=(1, 2, 3))
        windows[q] = w
    print("\n  given [A], the [B] window that still gives a condensate:")
    for q, w in windows.items():
        txt = ", ".join(f"{lo:.2f}-{hi:.0f} uM" for lo, hi in w) if w else "none"
        print(f"    [A] = {q:6.1f} uM  ->  [B] in {txt}")
    out["windows"] = {str(k): v for k, v in windows.items()}

    # valence / affinity chart
    kds = np.geomspace(1e-7, 3e-3, 22)
    vals = np.array([1, 2, 3, 4, 5, 6, 8])
    chart = np.full((len(vals), len(kds)), np.inf)
    for i, n in enumerate(vals):
        for j, kd in enumerate(kds):
            m = T.Mixture(names=["A"], diameters=[mix.diameters[0]],
                          valence=[float(n)], kd=np.array([[kd]]))
            r = T.csat_scan(m, 0, lo_uM=1e-3, hi_uM=1e4, n=160)
            chart[i, j] = r["csat_uM"]
    out["valence_chart"] = dict(kd_M=kds.tolist(), valence=vals.tolist(),
                                csat_uM=np.where(np.isfinite(chart), chart, -1).tolist())
    print("\n  valence gate (homotypic only):  n=1 never condenses at any Kd -> "
          f"{np.isfinite(chart[0]).sum()}/{len(kds)} finite c_sat")
    print(f"                                  n=3 condenses for "
          f"{np.isfinite(chart[2]).sum()}/{len(kds)} of the Kd range")

    # ------------------------------------------------------------ inverse
    rng = np.random.default_rng(0)
    probe_B = np.array([0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])
    truth, _, _, _ = I.csat_curve(mix, probe_B)
    noisy = truth * 10 ** rng.normal(0, 0.12, truth.shape)
    obs_titration = [I.Observation("csat", value=float(v), partner_uM=float(b),
                                   sd_log10=0.12)
                     for b, v in zip(probe_B, noisy) if np.isfinite(v)]
    obs_single = [obs_titration[3]]                     # one point only

    grid_aa = np.geomspace(3e-5, 3e-3, 13)
    grid_ab = np.geomspace(3e-6, 3e-4, 13)
    grid_bb = np.array([np.inf])                        # B has no self-affinity (known)

    fits = {}
    for tag, obs in (("one point", obs_single), ("titration", obs_titration)):
        f = I.fit_from_titration(make_mix(), obs, grid_aa, grid_ab, grid_bb)
        fits[tag] = f
        m = I.marginals(f)
        ci_aa = I.credible_interval(grid_aa, m[0])
        ci_ab = I.credible_interval(grid_ab, m[1])
        print(f"\n  inverse [{tag}]  n_obs = {len(obs)}")
        print(f"    MAP   Kd_AA = {f['map'][0]*1e6:8.1f} uM   (true {TRUE['kd_aa']*1e6:.0f})")
        print(f"          Kd_AB = {f['map'][1]*1e6:8.1f} uM   (true {TRUE['kd_ab']*1e6:.0f})")
        print(f"    95%   Kd_AA in [{ci_aa[0]*1e6:.1f}, {ci_aa[1]*1e6:.1f}] uM "
              f"({np.log10(ci_aa[1]/ci_aa[0]):.2f} decades)")
        print(f"          Kd_AB in [{ci_ab[0]*1e6:.1f}, {ci_ab[1]*1e6:.1f}] uM "
              f"({np.log10(ci_ab[1]/ci_ab[0]):.2f} decades)")
        out[f"fit_{tag.replace(' ','_')}"] = dict(
            map=list(f["map"]), ci_aa=list(ci_aa), ci_ab=list(ci_ab), n_obs=len(obs))

    out["sloppiness"] = {}
    for tag in ("one point", "titration"):
        sl = I.sloppiness(fits[tag])
        print(f"\n  sloppy-model analysis [{tag}]:")
        for sd, vec in zip(sl["sd_decades"], sl["eigenvectors"].T):
            combo = " ".join(f"{c:+.2f}*{l.split()[-1]}" for c, l in zip(vec, sl["labels"])
                             if abs(c) > 0.15)
            print(f"    +-{sd:7.3f} decades along  {combo}")
        print(f"    condition number = {sl['condition']:.1f}")
        out["sloppiness"][tag] = dict(sd_decades=sl["sd_decades"].tolist(),
                                      eigenvectors=sl["eigenvectors"].tolist(),
                                      condition=sl["condition"],
                                      labels=sl["labels"])

    # ------------------------------------------------------------ mRNA band
    tpm_probe = np.array([100.0, 400.0, 1600.0])
    med, lo, hi, samples = U.protein_from_mrna(tpm_probe, n_samples=4000,
                                               rng=np.random.default_rng(1))
    cs_A = float(csat[0])
    print(f"\n  mRNA propagation (c_sat = {cs_A:.3f} uM):")
    mrows = []
    for t, m_, l_, h_, sm in zip(tpm_probe, med, lo, hi, samples):
        p_ps = float((sm > cs_A).mean())
        print(f"    {t:6.0f} TPM -> [A] median {m_:7.3f} uM  (95% band {l_:.3f}-{h_:.3f})"
              f"   P(condensate) = {p_ps:.2f}")
        mrows.append(dict(tpm=float(t), median_uM=float(m_), lo_uM=float(l_),
                          hi_uM=float(h_), p_condensate=p_ps))
    print("    -> the answer from transcript data is a probability, never a yes/no")
    out["mrna"] = dict(csat_uM=cs_A, rows=mrows)

    # ------------------------------------------------------------ figures
    fig, ax = plt.subplots(2, 4, figsize=(17.0, 7.6))

    a = ax[0, 0]
    im = a.pcolormesh(cB, cA, pd["label"], cmap=CMAP,
                      norm=BoundaryNorm(np.arange(-0.5, 4), CMAP.N), shading="auto")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("[B]  uM"); a.set_ylabel("[A]  uM")
    a.set_title("phase diagram", fontsize=10)
    cb = fig.colorbar(im, ax=a, ticks=range(4), fraction=0.046)
    cb.ax.set_yticklabels([P.LABELS[i] for i in range(4)], fontsize=6.5)

    a = ax[0, 1]
    a.loglog(cB_tit, csat, "o-", color=ACC[0], ms=4)
    a.axhline(csat[0], color="#8b949e", ls=":", lw=1)
    jm = int(np.nanargmin(csat))
    a.plot([cB_tit[jm]], [csat[jm]], "*", color=ACC[1], ms=13)
    a.set_xlabel("[B]  uM"); a.set_ylabel("c$_{sat}$(A)  uM")
    a.set_title("titration: re-entrant boundary", fontsize=10)
    a.annotate("excess B\ndissolves it", xy=(cB_tit[-4], csat[-4]),
               xytext=(0.35, 0.75), textcoords="axes fraction", fontsize=8,
               color="#e6edf3", arrowprops=dict(arrowstyle="->", color="#8b949e"))

    a = ax[0, 2]
    with np.errstate(invalid="ignore"):
        im = a.pcolormesh(cB, cA, np.log10(np.where(pd["K_B"] > 0, pd["K_B"], np.nan)),
                          cmap="magma", shading="auto")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("[B]  uM"); a.set_ylabel("[A]  uM")
    a.set_title(r"log$_{10}$ partition coefficient of B", fontsize=10)
    fig.colorbar(im, ax=a, fraction=0.046)

    a = ax[0, 3]
    finite = np.isfinite(chart)
    plotc = np.where(finite, np.log10(np.where(finite, chart, 1)), np.nan)
    cm = plt.get_cmap("viridis").copy(); cm.set_bad("#161b22")
    im = a.pcolormesh(kds * 1e6, vals, np.ma.masked_invalid(plotc), cmap=cm,
                      shading="auto")
    a.text(0.04, 0.9, "n = 1: never,\nat any affinity", transform=a.transAxes,
           fontsize=7.5, color="#8b949e")
    a.set_xscale("log")
    a.set_xlabel("homotypic K$_d$  uM"); a.set_ylabel("valence n")
    a.set_title(r"log$_{10}$ c$_{sat}$;  blank = never condenses", fontsize=10)
    fig.colorbar(im, ax=a, fraction=0.046)

    for col, tag in enumerate(("one point", "titration")):
        a = ax[1, col]
        f = fits[tag]
        post2 = f["post"][:, :, 0]
        im = a.pcolormesh(grid_ab * 1e6, grid_aa * 1e6, post2, cmap="magma", shading="auto")
        a.plot([TRUE["kd_ab"] * 1e6], [TRUE["kd_aa"] * 1e6], "*", color=ACC[0], ms=15)
        a.plot([f["map"][1] * 1e6], [f["map"][0] * 1e6], "o", mfc="none",
               mec=ACC[1], ms=11, mew=2)
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel(r"K$_d^{AB}$  uM"); a.set_ylabel(r"K$_d^{AA}$  uM")
        a.set_title(f"posterior, {tag} ({len(obs_single if tag=='one point' else obs_titration)} obs)"
                    "\n★ truth   ○ MAP", fontsize=9.5)
        fig.colorbar(im, ax=a, fraction=0.046)

    a = ax[1, 2]
    for tag, c in (("one point", ACC[2]), ("titration", ACC[0])):
        m = I.marginals(fits[tag])
        a.semilogx(grid_aa * 1e6, m[0] / m[0].max(), "-o", color=c, ms=3,
                   label=f"{tag}: Kd_AA")
        a.semilogx(grid_ab * 1e6, m[1] / m[1].max(), "--s", color=c, ms=3, alpha=0.6,
                   label=f"{tag}: Kd_AB")
    a.axvline(TRUE["kd_aa"] * 1e6, color="#8b949e", ls=":", lw=1)
    a.axvline(TRUE["kd_ab"] * 1e6, color="#8b949e", ls=":", lw=1)
    a.set_xlabel("K$_d$  uM"); a.set_ylabel("marginal posterior")
    a.set_title("what the data pins down", fontsize=10)
    a.legend(fontsize=6.5)

    a = ax[1, 3]
    qs = sorted(windows)
    los = [windows[q][0][0] if windows[q] else np.nan for q in qs]
    his = [windows[q][0][1] if windows[q] else np.nan for q in qs]
    a.fill_between(qs, los, his, color=ACC[0], alpha=0.3, step="mid")
    a.plot(qs, los, "o-", color=ACC[0], ms=5, label="lower [B]")
    a.plot(qs, his, "s-", color=ACC[1], ms=5, label="upper [B]")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("measured [A]  uM"); a.set_ylabel("compatible [B]  uM")
    a.set_title("given [A], the [B] window\nthat still gives a condensate", fontsize=9.5)
    a.legend(fontsize=8)

    fig.suptitle("E09  Protein mixtures: Kd + concentration -> structure, and back again",
                 fontsize=12)
    fig.tight_layout()
    save(fig, "e09_protein_phases.png")
    dump(out, "e09_protein.json")


if __name__ == "__main__":
    main()
