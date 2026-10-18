"""
LF7 — Isothermal microcalorimetry: can total heat read out the "entropy
production of a cell state"?

Published inputs (PMID)
  resting human lymphocytes  2.32 ± 0.10 pW/cell     Valdemarsson & Monti 1994 (8156102)
  normal vs CLL lymphocytes  2.6 vs 1.8 pW/cell       Brandt 1979 (432550)
  NHL cells, progressive vs responding  5.5 vs 3.1 pW/cell   Monti 1981 (7346996)
  E. coli, single cell       1.3–4.5 pW               Vehusheia 2023 (37180454)
  TAM III detection limit    0.02 µW                  Braissant 2010 Sensors (22163413)
  oxycaloric equivalent      −470 kJ/mol O2           Gnaiger & Kemp 1990 (2184896)

Questions
 (1) heat budget: how large is a cell's total dissipation in k_B units, and
     how big a share could a state-defining gene circuit (the flux "curl" of
     LF2) possibly be?   [order-of-magnitude, assumptions stated]
 (2) detectability: cells needed to see the published state differences.
 (3) can single-cell fluctuations (critical slowing, "thermal fluctuation
     frequency") survive population averaging inside a calorimeter?
"""
import numpy as np
from . import _data

KB = 1.380649e-23
T = 310.15
KT = KB * T
NA = 6.02214076e23

LIT = {
    "resting lymphocyte": 2.32e-12,
    "normal lymphocyte (Brandt)": 2.6e-12, "CLL lymphocyte": 1.8e-12,
    "NHL responding": 3.1e-12, "NHL progressive": 5.5e-12,
}
TAM_LIMIT = 0.02e-6       # W
OXYCAL = 470e3            # J per mol O2
P_O2 = 5.0                # ATP per O2 (≈ P/O 2.5 per O atom), textbook value


def heat_budget():
    q = LIT["resting lymphocyte"]
    kT_per_s = q / KT
    o2 = q / OXYCAL * NA                 # O2 molecules per second
    atp = o2 * P_O2
    # assumed state-defining circuit: 10 TFs × 1e4 copies, 400 aa, 4 ATP/aa,
    # 1 h protein lifetime  +  1e4 kinase/phosphatase cycles per second (1 ATP each)
    tf_atp = 10 * 1e4 * 400 * 4 / 3600
    sig_atp = 1e4
    atp_kT = 20.0                        # ~ΔG of ATP hydrolysis in vivo ≈ 20 k_B T
    circuit_kT = (tf_atp + sig_atp) * atp_kT
    return dict(q=q, kT_per_s=kT_per_s, o2=o2, atp=atp, circuit_kT=circuit_kT,
                frac=circuit_kT / kT_per_s, tf_atp=tf_atp, sig_atp=sig_atp)


def detectability(snr=3.0):
    out = {}
    out["n_at_limit"] = TAM_LIMIT / LIT["resting lymphocyte"]
    for a, b in (("normal lymphocyte (Brandt)", "CLL lymphocyte"),
                 ("NHL responding", "NHL progressive")):
        dq = abs(LIT[a] - LIT[b])
        out[(a, b)] = snr * TAM_LIMIT / dq
    return out


def fluctuation_averaging(N=1e5, cv_cell=0.3, tau_cell=600.0, tau_instr=300.0,
                          dt=10.0, t_end=4 * 3600, seed=0):
    """N independent cells, each heat output an OU process (mean q, CV cv_cell,
    correlation time tau_cell). The calorimeter sees the SUM, low-pass filtered
    with time constant tau_instr (assumed; isothermal instruments integrate over
    minutes). Returns relative fluctuation of the measured signal and the
    synchronised alternative where all cells share one OU process."""
    rng = np.random.default_rng(seed)
    q = LIT["resting lymphocyte"]
    n = int(t_end / dt)
    a = np.exp(-dt / tau_cell)
    # sum of N independent OUs is an OU with sd sqrt(N)·σ
    x_ind = np.zeros(n); x_syn = np.zeros(n)
    s = cv_cell * q
    for k in range(1, n):
        x_ind[k] = a * x_ind[k - 1] + np.sqrt(1 - a * a) * s * np.sqrt(N) * rng.normal()
        x_syn[k] = a * x_syn[k - 1] + np.sqrt(1 - a * a) * s * N * rng.normal()
    b = np.exp(-dt / tau_instr)

    def lowpass(x):
        y = np.zeros_like(x)
        for k in range(1, len(x)):
            y[k] = b * y[k - 1] + (1 - b) * x[k]
        return y
    mean = N * q
    y_ind = lowpass(x_ind)[n // 4:]; y_syn = lowpass(x_syn)[n // 4:]
    return dict(N=N, mean_W=mean, rel_ind=y_ind.std() / mean, rel_syn=y_syn.std() / mean,
                abs_ind=y_ind.std(), abs_syn=y_syn.std(), t=np.arange(n)[n // 4:] * dt,
                y_ind=y_ind, y_syn=y_syn)


def run(make_fig=True):
    print("=" * 72); print("LF7  Microcalorimetry: heat as a state readout"); print("=" * 72)
    B = heat_budget()
    print(f"[1] resting lymphocyte: {B['q']*1e12:.2f} pW = {B['kT_per_s']:.1e} k_B T/s "
          f"≈ {B['o2']:.1e} O2/s ≈ {B['atp']:.1e} ATP/s")
    print(f"    assumed state-defining circuit (10 TFs×1e4 copies + 1e4 phospho-cycles/s): "
          f"{B['circuit_kT']:.1e} k_B T/s = {B['frac']:.2%} of total heat")
    print("    -> the calorimeter reports TOTAL metabolism; the LF2-type circuit flux is "
          "a sub-percent term inside it")
    D = detectability()
    print(f"\n[2] TAM III limit {TAM_LIMIT*1e9:.0f} nW = {D['n_at_limit']:.0f} resting lymphocytes")
    for k, v in D.items():
        if isinstance(k, tuple):
            print(f"    {k[0]} vs {k[1]}: ≈ {v:.1e} cells to resolve the difference at SNR 3")
    F = fluctuation_averaging()
    print(f"\n[3] N = {F['N']:.0e} cells, single-cell CV 30%, τ_cell 10 min, instrument τ 5 min:")
    print(f"    independent cells   -> signal fluctuation {F['rel_ind']:.2%} of mean "
          f"({F['abs_ind']*1e9:.2f} nW; below the {TAM_LIMIT*1e9:.0f} nW limit)")
    print(f"    synchronised cells  -> signal fluctuation {F['rel_syn']:.1%} of mean "
          f"({F['abs_syn']*1e9:.0f} nW; visible)")
    print("    -> only collective (population-wide) transitions are visible; single-cell "
          "critical slowing is averaged away by 1/√N")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
        names = list(LIT); vals = [LIT[k] * 1e12 for k in names]
        ax[0].barh(names, vals, color=["gray", "tab:blue", "tab:orange", "tab:green", "tab:red"])
        ax[0].set_xlabel("heat output (pW / cell)"); ax[0].set_title("(a) published per-cell heat")
        parts = {"total heat\n(measured)": B["kT_per_s"], "TF circuit\n(assumed)": B["tf_atp"] * 20,
                 "phospho-cycles\n(assumed)": B["sig_atp"] * 20}
        ax[1].bar(list(parts), list(parts.values()), color=["k", "tab:purple", "tab:pink"])
        ax[1].set_yscale("log"); ax[1].set_ylabel("dissipation (k_B T / s per cell)")
        ax[1].set_title("(b) heat budget of one lymphocyte")
        th = (F["t"] - F["t"][0]) / 3600
        ax[2].plot(th, (F["mean_W"] + F["y_syn"]) * 1e6, label="synchronised")
        ax[2].plot(th, (F["mean_W"] + F["y_ind"]) * 1e6, label="independent")
        ax[2].axhspan(F["mean_W"] * 1e6 - TAM_LIMIT * 1e6, F["mean_W"] * 1e6 + TAM_LIMIT * 1e6,
                      color="gray", alpha=0.2, label="± detection limit")
        ax[2].set_xlabel("time (h)"); ax[2].set_ylabel("measured heat (µW)")
        ax[2].set_title(f"(c) 1e5 cells in a calorimeter"); ax[2].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(_data.fig("lf7_calorimetry.png"), dpi=130)
        plt.close(fig)
    return dict(budget=B, detect=D, fluct=F)


if __name__ == "__main__":
    run()
