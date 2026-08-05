r"""
E08 - Pattern cards: one parameter set per morphology, with its model.
======================================================================

For each morphology the platform can produce, this experiment fixes an explicit,
minimal parameter set (the "card"), runs it, and reports side by side:

  * the **card**       - A, alpha, beta, friction, density: exactly what to type in
  * the **model**      - which term of the theory is responsible
  * the **prediction** - k*, l* = 2 pi / k*, growth rate, regime, mode coherence,
                         and for the moving ones the cluster's internal force and
                         torque turned into a predicted drift speed and spin
  * the **measurement** - motility, spin, D2, psi6, cluster statistics, S(k) peak
  * a **sensitivity sweep** on the single knob that destroys the pattern

The three diagnostics that do most of the classifying:

    nu_f        0 -> frozen (theorem), > 0 -> can move
    coherence   1 -> all species pile together (droplets)
                0 -> species separate in antiphase (core-shell, membranes)
    D2          ~1 -> filaments, ~2 -> blobs/sheets
"""

import multiprocessing as mp

import numpy as np
from common import ACC, dump, save, species_colors
import matplotlib.pyplot as plt

from plife import kernels as K, matrices as M, model, observables as OB
from plife import stability as ST, twobody as TB

MAIN = dict(N=3000, L=1200.0, T_EQ=35.0, T_MEAS=8.0)
SWEEP = dict(N=2000, L=1000.0, T_EQ=22.0, T_MEAS=5.0)


# --------------------------------------------------------------------------- #
#  helpers for building the cards
# --------------------------------------------------------------------------- #
def cyclic(S, self_a, fwd, bwd, other=0.0):
    """A_ii = self_a, A_{i,i+1} = fwd, A_{i,i-1} = bwd, everything else = other."""
    A = np.full((S, S), float(other))
    for i in range(S):
        A[i, i] = self_a
        A[i, (i + 1) % S] = fwd
        A[i, (i - 1) % S] = bwd
    return A


def wavefield(S, self_a, amp):
    A = np.full((S, S), 0.0)
    for i in range(S):
        for j in range(S):
            A[i, j] = self_a if i == j else amp * np.sin(2 * np.pi * (j - i) / S)
    return A


def P(A, alpha, beta, friction=0.30, L=None, repel=1.0, force_factor=1.0):
    """Assemble Params with scalar or matrix radii."""
    S = len(A)
    A = np.asarray(A, float)
    al = np.full((S, S), alpha, float) if np.isscalar(alpha) else np.asarray(alpha, float)
    be = np.full((S, S), beta, float) if np.isscalar(beta) else np.asarray(beta, float)
    LL = MAIN["L"] if L is None else L
    return model.Params(A=A, alpha=al, beta=be, L=(LL, LL), friction=friction,
                        repel=repel, force_factor=force_factor)


# --------------------------------------------------------------------------- #
#  the cards
# --------------------------------------------------------------------------- #
def card_gas(x=-0.6, L=None):
    return P(np.full((3, 3), x), 15.0, 45.0, friction=0.30, L=L)


def card_blobs(alpha=15.0, L=None):
    A = np.full((3, 3), -0.30)
    np.fill_diagonal(A, 0.90)
    return P(A, alpha, 3.0 * alpha, friction=0.30, L=L)


def card_crystal(repel=4.0, L=None):
    # one effective species (all pairs identical), narrow attraction band and a
    # stiff core -> every bond locks at r = alpha -> hexagonal packing
    A = np.full((3, 3), 0.90)
    return P(A, 16.0, 26.0, friction=0.30, repel=repel, L=L)


def card_cells(a11=-0.50, L=None):
    # species 0 = core (self-attracting), species 1 = membrane (self-repelling,
    # bound to the core).  A symmetric -> reciprocal -> the cells stay put.
    A = np.array([[0.90, 0.45], [0.45, a11]])
    al = np.array([[12.0, 22.0], [22.0, 16.0]])
    be = np.array([[42.0, 60.0], [60.0, 48.0]])
    return P(A, al, be, friction=0.30, L=L)


def card_chains(fwd=0.35, L=None):
    return P(cyclic(5, 1.0, fwd, fwd, other=-0.9), 15.0, 45.0, friction=0.30, L=L)


def card_snakes(theta=1.0, L=None):
    """Chain backbone + a one-directional link along the type cycle."""
    base = cyclic(5, 1.0, 0.35, 0.35, other=-0.9)
    anti = cyclic(5, 0.0, 0.45, -0.45, other=0.0)
    return P(np.clip(base + theta * anti, -1, 1), 15.0, 45.0, friction=0.25, L=L)


def card_chasers(friction=0.10, L=None):
    # A_S = diag(0.8) exactly; A_N = +-0.9 off-diagonal: pure chase, no cross bonding
    A = np.array([[0.80, 0.90], [-0.90, 0.80]])
    return P(A, 14.0, 42.0, friction=friction, L=L)


def card_rotors(bwd=-0.20, L=None):
    # all three species co-locate (self 0.9, forward 0.6) so the cluster stays
    # compact, while forward != backward gives it a net internal torque
    return P(cyclic(3, 0.90, 0.60, bwd), 14.0, 42.0, friction=0.20, L=L)


def card_waves(friction=0.05, L=None):
    return P(wavefield(5, 0.60, 0.80), 15.0, 45.0, friction=friction, L=L)


def card_dimers(L=None, a10=26.0):
    # strong self-repulsion keeps like particles apart, so what is left is a gas
    # of isolated bound pairs, each obeying the exact two-body result
    # V = mu_disc f_A(r*).  The knob that matters here is the DENSITY: the
    # two-body formula is a dilute-limit statement.
    A = np.array([[-0.90, 0.90], [-0.40, -0.90]])
    al = np.array([[15.0, 15.0], [a10, 15.0]])
    be = np.full((2, 2), 45.0)
    return P(A, al, be, friction=0.10, L=L if L is not None else 1500.0)


CARDS = [
    dict(key="gas", zh="均匀气体", en="Gas",
         build=card_gas, sweep=("A_ij  (all pairs)", [-0.6, -0.2, 0.1, 0.4]),
         mech="所有 k 上 Re λ < 0：均匀态线性稳定，没有不稳定波段。"),
    dict(key="blobs", zh="静态团块 / 液滴", en="Static droplets",
         build=card_blobs, sweep=(r"min radius $\alpha$", [10.0, 15.0, 22.0, 30.0]),
         mech="对角自吸引 + 非对角互斥，A 对称 ⇒ ν_f = 0 ⇒ 定理保证冻结。"
              "最不稳定模是 coherence ≈ 0 的**分相模**（各物种分别成滴，而不是混在一起），"
              "尺度 ℓ* ∝ α。"),
    dict(key="crystal", zh="晶体团簇", en="Crystalline clusters",
         build=card_crystal, sweep=("core repulsion R", [0.5, 1.5, 3.0, 6.0]),
         mech="窄吸引带 (β/α = 1.67) + 硬核 ⇒ 键长严格锁在 α，六角密堆，|ψ₆| 高。"),
    dict(key="cells", zh="核-壳细胞", en="Core-shell cells",
         build=card_cells, sweep=(r"membrane self-coupling $A_{11}$",
                                  [-0.9, -0.5, -0.1, 0.3]),
         mech="不稳定本征向量分量反号（coherence ≈ 0）⇒ 分相模：核与膜互相排开又互相束缚。"
              "α₀₁ > α₀₀ 把膜推到外圈。"),
    dict(key="chains", zh="链 / 网络", en="Chains & networks",
         build=card_chains, sweep=("cycle coupling w = b", [0.0, 0.2, 0.35, 0.8]),
         mech="类型空间的环图 ⇒ 实空间的一维流形；A 对称 ⇒ 静止。粗细 = α。"),
    dict(key="snakes", zh="爬行的蛇", en="Crawling snakes",
         build=card_snakes, sweep=(r"non-reciprocity $\theta$", [0.0, 0.35, 0.7, 1.0]),
         mech="同样的环图，但正向/反向耦合不等 ⇒ 每条键都是一次追逐 ⇒ "
              "F_int 沿链方向不为零，整条蛇自己往前爬。"),
    dict(key="chasers", zh="追逐者 / 迁移细胞", en="Chasers", clust_cut=2.8,
         build=card_chasers, sweep=(r"friction $\lambda$", [0.30, 0.15, 0.10, 0.04]),
         mech="A_S = diag(0.8) 负责成团，A_N = ±0.9 负责追逐。交叉的对称部分为 0，"
              "所以两个物种**不共处一团**——绿团追着品红团跑，驱动力是**团簇之间**的，"
              "不是团簇内部的（因此用较大的聚类半径 2.8α 才能把追逐对当成一个单元）。"
              "摩擦越低，行波不稳定性越强（判据 σᵢ² > γ²σᵣ/(60κk²)）。"),
    dict(key="rotors", zh="自旋细胞 / 旋涡", en="Rotating cells",
         build=card_rotors, sweep=("backward coupling b", [0.60, 0.20, -0.20, -0.60]),
         mech="三个物种互相吸引所以待在同一团里，但正向 0.6 ≠ 反向 b ⇒ "
              "团簇净内力矩 τ_int ≠ 0 ⇒ 自旋。b = 0.6（= 正向）时矩阵变成互易，"
              "τ_int 精确归零，转动消失。"),
    dict(key="waves", zh="行波 / 传送带", en="Travelling waves",
         build=card_waves, sweep=(r"friction $\lambda$", [0.30, 0.15, 0.08, 0.03]),
         mech="纯反对称 sin 耦合 ⇒ A_S 只剩对角 ⇒ 几乎全部动力学都是非变分的；"
              "低摩擦下均匀态直接失稳成行波。"),
    dict(key="dimers", zh="自驱动二聚体气体", en="Self-propelled dimers", clust_cut=1.45,
         N=250, L=1800.0, dimer_probe=True,
         # value = box size; N chosen so the sweep spans rho*beta^2 = 2.25 -> 0.06
         sweep_N={900.0: 900, 1200.0: 600, 1500.0: 400, 1800.0: 250, 2200.0: 150},
         build=card_dimers,
         sweep=(r"box size L  (density $\rho\beta^2$: 2.3 $\to$ 0.06)",
                [900.0, 1200.0, 1500.0, 1800.0, 2200.0]),
         mech="自斥阻止聚团，剩下孤立的束缚对。半径矩阵不对称 ⇒ f_A(r*) ≠ 0 ⇒ "
              "每个二聚体以 V = μ_disc f_A(r*) 匀速漂移。"),
]


# --------------------------------------------------------------------------- #
def analyse(p, N, L, T_EQ, T_MEAS, seed=0, want_snapshot=False, nk=350,
            clust_cut=1.35, t_early=6.0, probe=False):
    sim = model.ParticleLife(p, N=N, seed=seed)
    # l* from linear theory describes the FIRST structures, so the fair
    # comparison is an early-time S(k) peak, before nonlinear coarsening.
    sim.run(t_early)
    k_early = OB.peak_wavenumber(sim.pos, np.array([L, L]), ngrid=256,
                                 kmax=3 * np.pi / float(np.median(p.alpha)))
    sim.run(T_EQ - t_early)
    o1 = OB.observables(sim)
    sim.run(T_MEAS)
    o2 = OB.observables(sim)
    o = {k: (0.5 * (o1[k] + o2[k]) if isinstance(o1[k], float) else o1[k]) for k in o1}

    rho = np.full(p.S, N / p.S / L**2)
    lin = ST.scan(p, rho, nk=nk)
    mc = ST.mode_character(p, rho, k=lin.k_star)
    cl = TB.cluster_prediction(sim, cutoff_scale=clust_cut)

    out = dict(
        nu_A=M.nonreciprocity(p.A),
        nu_f=M.nonreciprocity_forces(p.A, p.alpha, p.beta, p.repel),
        k_star=lin.k_star, length_pred=lin.length, growth=lin.growth,
        osc_growth=lin.osc_growth, osc_speed=lin.osc_speed, regime=lin.regime(),
        coherence=mc["coherence"],
        motility=o["motility"], speed=o["speed"], spin=o["spin"],
        abs_spin=o["abs_spin"], dim2=o["dim2"], psi6=o["psi6"],
        n_clusters=o["n_clusters"], mean_cluster=o["mean_cluster"],
        length_meas=o["length"], label=OB.classify(o),
        length_early=float(2 * np.pi / k_early) if k_early == k_early and k_early > 0 else np.nan,
        psi6_global=o.get("psi6_global", np.nan),
        cluster=cl, gamma=p.gamma,
    )
    if probe:
        out["dimer"] = dimer_probe(sim, clust_cut)
    if want_snapshot:
        out["pos"] = sim.pos.tolist()
        out["types"] = sim.types.tolist()
    return out


def dimer_probe(sim, cut):
    r"""Isolated bound pairs: measured bond length and drift vs the exact two-body result."""
    p = sim.p
    box = np.asarray(p.L, float)
    lab, _ = OB.clusters(sim.pos, box, cut * float(np.median(p.alpha)))
    cnt = np.bincount(lab)
    speeds, rs = [], []
    for c in np.where(cnt == 2)[0]:
        idx = np.where(lab == c)[0]
        if sim.types[idx[0]] == sim.types[idx[1]]:
            continue                      # like pairs carry no non-reciprocal drive
        speeds.append(float(np.hypot(*sim.vel[idx].mean(0))))
        d = sim.pos[idx[1]] - sim.pos[idx[0]]
        d -= box * np.round(d / box)
        rs.append(float(np.hypot(*d)))
    th = TB.dimer(p, 0, 1)
    return dict(n_pairs=len(speeds),
                r_theory=float(th["r_star"]) if th["r_star"] else np.nan,
                r_meas=float(np.mean(rs)) if rs else np.nan,
                V_theory=float(abs(th["V"])),
                V_meas=float(np.mean(speeds)) if speeds else np.nan,
                rho_beta2=float(sim.N / (box[0] * box[1]) * float(p.beta.max()) ** 2))


def job(args):
    kind, key, value = args
    card = next(c for c in CARDS if c["key"] == key)
    cut = card.get("clust_cut", 1.35)
    if kind == "main":
        L = card.get("L", MAIN["L"])
        N = card.get("N", MAIN["N"])
        p = card["build"](L=L) if card.get("N") else card["build"]()
        r = analyse(p, N, L, MAIN["T_EQ"], MAIN["T_MEAS"],
                    want_snapshot=True, clust_cut=cut, probe=card.get("dimer_probe"))
        r.update(kind="main", key=key, value=None, box=L,
                 A=p.A.tolist(), alpha=p.alpha.tolist(), beta=p.beta.tolist(),
                 friction=p.friction, repel=p.repel, S=p.S, N=N, rho=N / L ** 2)
        return r
    if "sweep_N" in card:                     # value is the box size
        L = float(value)
        N = card["sweep_N"][value]
        p = card["build"](L=L)
    else:
        L, N = SWEEP["L"], SWEEP["N"]
        p = card["build"](value, L=L)
    r = analyse(p, N, L, SWEEP["T_EQ"] + 14.0, SWEEP["T_MEAS"], nk=250,
                clust_cut=cut, t_early=5.0, probe=card.get("dimer_probe"))
    r.update(kind="sweep", key=key, value=float(value), box=L, N=N)
    return r


# --------------------------------------------------------------------------- #
def draw_card(card, main, sweeps, path):
    fig = plt.figure(figsize=(12.4, 3.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.05, 1, 1, 1.05], wspace=0.32)
    S = main["S"]
    cols = species_colors(S)

    # --- snapshot
    ax = fig.add_subplot(gs[0, 0])
    pos, ty = np.array(main["pos"]), np.array(main["types"])
    for s in range(S):
        m = ty == s
        ax.scatter(pos[m, 0], pos[m, 1], s=1.1, c=cols[s], linewidths=0, alpha=0.9)
    LL = main.get("box", MAIN["L"])
    ax.set_xlim(0, LL); ax.set_ylim(0, LL)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    ax.set_title(card["en"], fontsize=10)

    # --- interaction matrix
    ax = fig.add_subplot(gs[0, 1])
    A = np.array(main["A"])
    im = ax.imshow(A, cmap="PiYG", vmin=-1, vmax=1)
    for i in range(S):
        for j in range(S):
            ax.text(j, i, f"{A[i, j]:.2f}", ha="center", va="center", fontsize=6.6,
                    color="#0d1117" if abs(A[i, j]) > 0.45 else "#e6edf3")
    ax.set_xticks(range(S)); ax.set_yticks(range(S))
    ax.set_xticklabels(range(S), fontsize=7); ax.set_yticklabels(range(S), fontsize=7)
    ax.set_xlabel("exerted by", fontsize=8); ax.set_ylabel("felt by", fontsize=8)
    ax.grid(False)
    ax.set_title(r"$A$    " + rf"$\nu_f$={main['nu_f']:.2f}", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046)

    # --- dispersion
    ax = fig.add_subplot(gs[0, 2])
    p = card["build"](L=main.get("box", MAIN["L"])) if card.get("N") else card["build"]()
    rho = np.full(S, main["N"] / S / main.get("box", MAIN["L"]) ** 2)
    r = ST.scan(p, rho, nk=400)
    re = r.lam.real.reshape(len(r.k), -1).max(axis=1)
    im_ = np.abs(r.lam.imag.reshape(len(r.k), -1)).max(axis=1)
    ax.plot(r.k, re, color=ACC[2], lw=1.6)
    osc = (im_ > 1e-9) & (re > 0)
    if osc.any():
        ax.plot(r.k[osc], re[osc], ".", color=ACC[1], ms=3.2, label="oscillatory")
        ax.legend(fontsize=7)
    ax.axhline(0, color="#8b949e", lw=0.6)
    if r.growth > 0:
        ax.axvline(r.k_star, color=ACC[0], ls="--", lw=1)
    ax.set_xlabel("k", fontsize=8); ax.set_ylabel(r"max Re $\lambda$", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(rf"$\ell^*$={main['length_pred']:.0f} (meas {main['length_early']:.0f})"
                 rf", coh={main['coherence']:.2f}"
                 f"\n{main['regime']}", fontsize=9)

    # --- sensitivity sweep
    ax = fig.add_subplot(gs[0, 3])
    lab, vals = card["sweep"]
    xs = [s["value"] for s in sweeps]
    ax.plot(xs, [s["motility"] for s in sweeps], "o-", color=ACC[1], ms=5,
            label="motility")
    ax2 = ax.twinx()
    ax2.plot(xs, [s["dim2"] for s in sweeps], "s--", color=ACC[2], ms=4, label=r"$D_2$")
    ax2.plot(xs, [s["nu_f"] for s in sweeps], "^:", color=ACC[3], ms=4, label=r"$\nu_f$")
    ax2.grid(False); ax2.tick_params(labelsize=7)
    ax.set_xlabel(lab, fontsize=8); ax.tick_params(labelsize=7)
    ax.set_ylabel("motility", color=ACC[1], fontsize=8)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=6.5, loc="best")
    ax.set_title("sensitivity", fontsize=9)

    fig.tight_layout()
    save(fig, path)


if __name__ == "__main__":
    jobs = [("main", c["key"], None) for c in CARDS]
    for c in CARDS:
        for v in c["sweep"][1]:
            jobs.append(("sweep", c["key"], v))
    with mp.Pool(min(4, mp.cpu_count())) as pool:
        res = pool.map(job, jobs)

    mains = {r["key"]: r for r in res if r["kind"] == "main"}
    sweeps = {}
    for r in res:
        if r["kind"] == "sweep":
            sweeps.setdefault(r["key"], []).append(r)
    for k in sweeps:
        sweeps[k].sort(key=lambda x: x["value"])

    hdr = (f"{'pattern':22s} {'nu_f':>5s} {'coh':>5s} {'l*':>6s} {'l_early':>8s} "
           f"{'l_late':>7s} {'regime':>14s} {'motil':>8s} {'|spin|':>7s} "
           f"{'D2':>5s} {'psi6':>5s}  label")
    print(hdr); print("-" * len(hdr))
    for c in CARDS:
        m = mains[c["key"]]
        print(f"{c['zh']:22s} {m['nu_f']:5.2f} {m['coherence']:5.2f} "
              f"{m['length_pred']:6.0f} {m['length_early']:8.0f} "
              f"{m['length_meas']:7.0f} {m['regime']:>14s} "
              f"{m['motility']:8.4f} {m['abs_spin']:7.3f} {m['dim2']:5.2f} "
              f"{m['psi6']:5.2f}  {m['label']}")

    print(f"\n{'pattern':22s} {'N_C':>5s} {'|F_int|':>9s} {'V_pred':>9s} {'V_meas':>9s} "
          f"{'tau_int':>9s} {'w_pred':>9s} {'w_meas':>9s}")
    for c in CARDS:
        cl = mains[c["key"]]["cluster"]
        print(f"{c['zh']:22s} {cl['size']:5d} {cl['F_int']:9.4f} {cl['V_pred']:9.4f} "
              f"{cl['V_meas']:9.4f} {cl['tau_int']:9.4f} {cl['omega_pred']:+9.4f} "
              f"{cl['omega_meas']:+9.4f}")

    for c in CARDS:
        if "dimer" in mains[c["key"]]:
            print(f"\nisolated-dimer probe ({c['zh']}):")
            print(f"{'rho*beta^2':>11s} {'#pairs':>7s} {'r_th':>7s} {'r_meas':>7s} "
                  f"{'V_th':>7s} {'V_meas':>7s} {'ratio':>7s}")
            for row in [mains[c["key"]]["dimer"]] + [s_["dimer"] for s_ in sweeps[c["key"]]]:
                ratio = row["V_meas"] / row["V_theory"] if row["V_theory"] else np.nan
                print(f"{row['rho_beta2']:11.3f} {row['n_pairs']:7d} {row['r_theory']:7.3f} "
                      f"{row['r_meas']:7.3f} {row['V_theory']:7.3f} {row['V_meas']:7.3f} "
                      f"{ratio:7.3f}")

    for c in CARDS:
        draw_card(c, mains[c["key"]], sweeps[c["key"]], f"e08_card_{c['key']}.png")

    dump(dict(
        cards=[dict(key=c["key"], zh=c["zh"], en=c["en"], mech=c["mech"],
                    sweep_label=c["sweep"][0],
                    main={k: v for k, v in mains[c["key"]].items()
                          if k not in ("pos", "types")},
                    sweep=[{k: v for k, v in s.items()} for s in sweeps[c["key"]]])
               for c in CARDS],
        settings=dict(main=MAIN, sweep=SWEEP),
    ), "e08_cards.json")
