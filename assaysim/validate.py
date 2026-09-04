"""
验证层 —— 每一条都是可证伪的检查，失败会让退出码非零。

分三类:
  [解析]  与闭式解 / 特征值 / 守恒律比对，容差到数值精度
  [交叉]  与独立第三方实现比对 (Biopython)
  [真值]  与 ChEMBL 真实实验数据比对 (留出集)

运行:  python3 -m assaysim.validate
"""

from __future__ import annotations

import math
import sys

import numpy as np

from . import neutralization as neu
from . import viral_dynamics as vd
from .nn_thermo import duplex_thermo, mismatch_penalty_dG

R_CAL = 1.98720425864083
T37 = 310.15

_FAILURES: list[str] = []
_CHECKS = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global _CHECKS
    _CHECKS += 1
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        _FAILURES.append(name)


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


# --- M1 ------------------------------------------------------------------

def validate_nn_thermo() -> None:
    print("\nM1 最近邻热力学")

    # [交叉] 与 Biopython 的 Allawi & SantaLucia 1997 表比对
    try:
        import random

        from Bio.SeqUtils import MeltingTemp as mt

        random.seed(0)
        seqs = [
            "".join(random.choice("ACGT") for _ in range(n))
            for n in (18, 20, 22, 25)
            for _ in range(50)
        ]
        diffs = [
            duplex_thermo(s, ct_nM=250.0, salt_correct=False).tm
            - mt.Tm_NN(s, nn_table=mt.DNA_NN3, dnac1=125.0, dnac2=125.0,
                       saltcorr=0, Na=1000)
            for s in seqs
        ]
        md = max(abs(d) for d in diffs)
        check("[交叉] Tm 与 Biopython DNA_NN3 (Allawi97) 一致",
              md < 0.01, f"max|Δ| = {md:.5f} °C, n={len(seqs)}")
    except ImportError:
        check("[交叉] Biopython 可用", False, "未安装 biopython")

    # [解析] 热力学自洽 ΔG = ΔH - TΔS
    t = duplex_thermo("GTCATGACGTCATGAC")
    dg = t.dH - T37 * t.dS / 1000.0
    check("[解析] ΔG37 = ΔH − TΔS 自洽", close(dg, t.dG37, 1e-10),
          f"Δ = {abs(dg - t.dG37):.2e} kcal/mol")

    # [解析] 3' 端错配的惩罚必须重于 5' 端
    L = 20
    p3 = mismatch_penalty_dG([L - 1], L)
    p5 = mismatch_penalty_dG([0], L)
    check("[解析] 错配惩罚 3' 端 > 5' 端", p3 > p5 * 3,
          f"3'={p3:.2f} vs 5'={p5:.2f} kcal/mol")

    # [解析] GC 含量高 -> Tm 高 (同长度)
    gc_rich = duplex_thermo("GCGCGCGCGCGCGCGCGCGC").tm
    at_rich = duplex_thermo("ATATATATATATATATATAT").tm
    check("[解析] GC 富集序列 Tm 更高", gc_rich > at_rich + 20,
          f"GC20={gc_rich:.1f} °C, AT20={at_rich:.1f} °C")


# --- M3 ------------------------------------------------------------------

def validate_neutralization() -> None:
    print("\nM3 多击中占据模型")

    # [解析] k=1 闭式解 vs 数值根求解
    worst = 0.0
    for n in (1, 2, 5, 14, 25, 100, 340):
        num = neu.nt50_from_kd(1.0, neu.Virion("t", n, 1))
        ana = neu.nt50_closed_form_single_hit(1.0, n)
        worst = max(worst, abs(num - ana) / ana)
    check("[解析] k=1 数值解 == 闭式解 Kd·(2^(1/n)−1)", worst < 1e-12,
          f"max 相对误差 = {worst:.2e}")

    # [解析] n=1, k=1 -> NT50 == Kd
    v = neu.Virion("single", 1, 1)
    check("[解析] n=1,k=1 时 NT50 == Kd",
          close(neu.nt50_from_kd(7.3, v), 7.3, 1e-12),
          f"NT50 = {neu.nt50_from_kd(7.3, v):.12f} (Kd=7.3)")

    # [解析] 报告里 n=100 的 1/144 断言
    amp = neu.amplification_factor(neu.Virion("n100", 100, 1))
    check("[解析] n=100,k=1 时 NT50/Kd = 2^0.01−1 ≈ 1/144",
          close(amp, 2 ** 0.01 - 1, 1e-14) and close(1 / amp, 143.77, 0.01),
          f"NT50/Kd = {amp:.6f}  即 Kd/{1 / amp:.2f}")

    # [解析] 固定 k=1 时，n 越大 NT50 越低 (多价放大) —— 严格单调
    amps = [neu.amplification_factor(neu.Virion("v", n, 1)) for n in range(1, 200)]
    check("[解析] k=1 时 NT50/Kd 随 n 严格单调下降",
          all(amps[i] > amps[i + 1] for i in range(len(amps) - 1)),
          f"n=1: {amps[0]:.3f} -> n=199: {amps[-1]:.5f}")

    # [解析] 占据模型的边界行为
    ok = (close(neu.residual_infectivity(0.0, 50, 3), 1.0, 1e-12)
          and close(neu.residual_infectivity(1.0, 50, 3), 0.0, 1e-12))
    check("[解析] P_inf(θ=0)=1, P_inf(θ=1)=0", ok)

    # [解析] k 越大越难中和 -> NT50 越高，严格单调
    n = 50
    nts = [neu.nt50_from_kd(1.0, neu.Virion("v", n, k)) for k in range(1, n + 1)]
    check("[解析] 固定 n 时 NT50 随击中阈值 k 严格单调上升",
          all(nts[i] < nts[i + 1] for i in range(len(nts) - 1)),
          f"k=1: {nts[0]:.4f} -> k={n}: {nts[-1]:.2f} (×Kd)")

    # [解析] 从模拟曲线反解 k —— 参数可辨识性的直接检验
    true_k, kd, nsp = 7, 3.0, 50
    concs = np.logspace(-3, 3, 40)
    resid = neu.neutralization_curve(concs, kd, neu.Virion("v", nsp, true_k))
    khat = neu.fit_k_from_curve(concs, resid, kd, nsp)
    check("[解析] 无噪声曲线可唯一反解击中阈值 k", khat == true_k,
          f"真值 k={true_k}, 拟合 k={khat}")


# --- M4 ------------------------------------------------------------------

def validate_viral_dynamics() -> None:
    print("\nM4 靶细胞受限病毒动力学")

    par = vd.REFERENCE_PARAMS["influenza_MDCK"]

    # [解析] 质量守恒 T+E+I+D 恒等于 T0
    tr = vd.simulate(par, moi=0.01, t_end_h=120)
    total = tr.T + tr.E + tr.I + tr.D
    err = float(np.max(np.abs(total - par.T0)) / par.T0)
    check("[解析] 细胞数守恒 T+E+I+D = T0", err < 1e-6,
          f"max 相对偏差 = {err:.2e}")

    # [解析] 指数期增长率 == 线性化系统最大特征值
    # [解析] 特征方程根 == 直接对雅可比做特征分解 (两条独立路径)
    lam = par.growth_rate()
    bT0 = par.beta * par.T0
    J = np.array([[-par.k_eclipse, 0.0, bT0],
                  [par.k_eclipse, -par.delta, 0.0],
                  [0.0, par.p, -(par.c + bT0)]])
    lam_eig = float(np.max(np.linalg.eigvals(J).real))
    check("[解析] 特征方程求根 == 雅可比特征分解",
          abs(lam - lam_eig) / lam_eig < 1e-9,
          f"brentq {lam:.9f} vs eig {lam_eig:.9f} /h")

    # [解析] 模拟的渐近指数斜率 == λ
    # 接种后有瞬态 (初期只有 V、尚无 E/I)，必须在瞬态之后取窗口
    tr2 = vd.simulate(par, moi=1e-9, t_end_h=40, n_points=4000)
    m = (tr2.V > 0) & (tr2.t > 20) & (tr2.t < 30)
    slope = float(np.polyfit(tr2.t[m], np.log(tr2.V[m]), 1)[0])
    check("[解析] 模拟渐近增长率 == 特征方程根 (λ+k)(λ+δ)(λ+c+βT0)=kpβT0",
          abs(slope - lam) / lam < 0.02,
          f"模拟 {slope:.5f} /h  vs 解析 {lam:.5f} /h  (R0={par.R0:.1f})")

    # [解析] R0 阈值：R0<1 感染熄灭，R0>1 扩增
    import dataclasses
    sub = dataclasses.replace(par, beta=par.beta * 0.5 / par.R0)   # R0 = 0.5
    sup = dataclasses.replace(par, beta=par.beta * 2.0 / par.R0)   # R0 = 2.0
    t_sub = vd.simulate(sub, moi=1e-4, t_end_h=200)
    t_sup = vd.simulate(sup, moi=1e-4, t_end_h=200)
    check("[解析] R0<1 感染熄灭 / R0>1 扩增",
          t_sub.V[-1] < t_sub.V[0] and t_sup.peak_titer() > 10 * t_sup.V[0],
          f"R0=0.5: V末/V初={t_sub.V[-1] / t_sub.V[0]:.2e}; "
          f"R0=2.0: 峰值/初值={t_sup.peak_titer() / t_sup.V[0]:.1f}")

    # [解析] 快清除极限下最终规模 -> Kermack-McKendrick 超越方程
    #        ln(T0/T∞) = R0 (1 − T∞/T0)
    fast = dataclasses.replace(par, c=par.c * 3000.0, p=par.p * 3000.0)  # 保持 R0
    trf = vd.simulate(fast, moi=1e-6, t_end_h=4000, n_points=4000)
    s_inf = trf.T[-1] / fast.T0
    lhs = -math.log(s_inf)
    rhs = fast.R0 * (1.0 - s_inf)
    check("[解析] 快清除极限的最终规模满足 KM 超越方程",
          abs(lhs - rhs) / rhs < 0.02,
          f"−ln(S∞)={lhs:.4f} vs R0(1−S∞)={rhs:.4f}, S∞={s_inf:.4f}, R0={fast.R0:.2f}")

    # [解析] 剂量-反应：孔水平表观 EC50 与分子层 EC50 不相等
    drug = vd.Drug("示例复制抑制剂", moa="replication", ec50_mol=100.0, hill=1.0)
    res = vd.apparent_ec50(par, drug, moi=0.01, readout_h=72)
    check("[解析] 剂量-反应曲线可拟合出表观 EC50",
          res["ec50_apparent"] == res["ec50_apparent"],
          f"表观 EC50 = {res['ec50_apparent']:.2f} nM  vs 分子层 {res['ec50_molecular']:.0f} nM  "
          f"(位移 {res['shift_log10']:+.2f} log)")

    # [解析] 表观 EC50 随 MOI 与读板时间系统性变化 (模型的可检验预言)
    shifts = {}
    for moi in (0.001, 0.01, 0.1):
        r = vd.apparent_ec50(par, drug, moi=moi, readout_h=72)
        shifts[("moi", moi)] = r["ec50_apparent"]
    spread_moi = max(shifts.values()) / min(shifts.values())
    check("[解析] 表观 EC50 随 MOI 变化 (同一分子、同一药)",
          spread_moi > 1.2,
          "MOI 0.001/0.01/0.1 -> EC50 " +
          " / ".join(f"{v:.1f}" for v in shifts.values()) + f" nM (跨度 {spread_moi:.2f}×)")

    times = {}
    for th in (24.0, 48.0, 72.0):
        r = vd.apparent_ec50(par, drug, moi=0.01, readout_h=th)
        times[th] = r["ec50_apparent"]
    spread_t = max(times.values()) / min(times.values())
    check("[解析] 表观 EC50 随读板时间变化",
          spread_t > 1.2,
          "24/48/72 h -> EC50 " +
          " / ".join(f"{v:.1f}" for v in times.values()) + f" nM (跨度 {spread_t:.2f}×)")

    # [解析] 中和抗体经占据模型进入 ODE，孔水平 NT50 应远低于 Kd
    ab = vd.Antibody("示例中和抗体", kd_nM=10.0,
                     virion=neu.KNOWN_VIRIONS["influenza_A"])
    r = vd.apparent_ec50(par, ab, moi=0.01, readout_h=72)
    check("[解析] 抗体孔水平 NT50 << Kd (多价放大穿透到细胞层)",
          r["ec50_apparent"] < ab.kd_nM,
          f"孔水平 NT50 = {r['ec50_apparent']:.4f} nM  vs Kd = {ab.kd_nM} nM  "
          f"({ab.kd_nM / r['ec50_apparent']:.0f}× 更低)")


# --- M6 ------------------------------------------------------------------

def validate_bridge() -> None:
    print("\nM6 非细胞 -> 细胞 桥接 (ChEMBL 真实数据)")
    try:
        from .bridge import evaluate, fit_linear, run_analysis
        from .chembl_data import build_pairs
    except ImportError as e:  # pragma: no cover
        check("[真值] 依赖可用", False, str(e))
        return

    try:
        hiv = build_pairs("HIV1_RT", verbose=False)
        flu = build_pairs("FLU_NA", verbose=False)
    except Exception as e:
        check("[真值] ChEMBL 数据可用", False, f"{type(e).__name__}: {e}")
        return

    check("[真值] HIV-1 RT 配对数据量充足", len(hiv) > 500, f"N = {len(hiv)}")
    check("[真值] 流感 NA 配对数据量充足", len(flu) > 50, f"N = {len(flu)}")

    a = run_analysis(hiv, seed=0)
    res = {r["kind"]: r for r in a["results"]}
    check("[真值] 线性桥接优于 identity (留出骨架划分)",
          res["linear"]["rmse"] < res["identity"]["rmse"],
          f"linear RMSE {res['linear']['rmse']:.3f} < identity {res['identity']['rmse']:.3f}")

    floor = a["noise_cell"]["approx_sd_log10"]
    check("[真值] 预测误差已接近细胞法自身噪声下限",
          res["linear"]["rmse"] < 2.0 * floor,
          f"RMSE {res['linear']['rmse']:.3f} vs 噪声下限 ≈ {floor:.2f} log "
          f"(比值 {res['linear']['rmse'] / floor:.2f})")

    # [真值] 跨体系迁移必须变差 —— "参数不可迁移" 的直接检验
    cross = evaluate(fit_linear(hiv), flu)
    native = evaluate(fit_linear(flu), flu)
    check("[真值] HIV 标定的桥接迁移到流感后显著变差",
          cross["rmse"] > native["rmse"] * 1.2 and abs(cross["bias"]) > 0.5,
          f"迁移 RMSE {cross['rmse']:.3f} (bias {cross['bias']:+.2f}) "
          f"vs 自标定 {native['rmse']:.3f} (bias {native['bias']:+.2f})")


def validate_agid() -> None:
    print("\nM5 琼脂免疫扩散 (AGID) 反应-扩散 PDE")
    from . import agid

    # [文献] Stokes-Einstein 对三个实测蛋白的复现
    worst = 0.0
    detail = []
    for name, p in agid.REFERENCE_PROTEINS.items():
        D = agid.stokes_einstein_D(p["r_h_nm"], T_K=293.15)
        rel = abs(D - p["D_measured_cm2_s"]) / p["D_measured_cm2_s"]
        worst = max(worst, rel)
        detail.append(f"{name} {D:.2e} vs 实测 {p['D_measured_cm2_s']:.2e}")
    check("[文献] Stokes-Einstein 复现 IgG/BSA/溶菌酶的实测扩散系数",
          worst < 0.08, f"最大相对偏差 {worst:.1%} · " + " · ".join(detail))

    # [文献] MW -> R_h 的经验式对三个实测半径的复现
    worst_r = 0.0
    dr = []
    for name, p in agid.REFERENCE_PROTEINS.items():
        rh = agid.hydrodynamic_radius_from_MW(p["MW_kDa"] * 1000.0, p["shape_factor"])
        rel = abs(rh - p["r_h_nm"]) / p["r_h_nm"]
        worst_r = max(worst_r, rel)
        dr.append(f"{name} {rh:.2f} vs 实测 {p['r_h_nm']}nm")
    check("[文献] MW -> R_h 经验式复现实测流体力学半径",
          worst_r < 0.10, f"最大相对偏差 {worst_r:.1%} · " + " · ".join(dr))

    # [解析] PDE 求解器 vs erfc 解析解 (无反应)
    D = 6.0e-7
    t = 3600.0
    x, num = agid.diffuse_numeric(L_cm=0.5, n=401, t_s=t, D=D, C0=1.0)
    ana = agid.diffuse_semi_infinite_analytic(x, t, D, 1.0)
    err = float(np.max(np.abs(num - ana)))
    check("[解析] 扩散求解器 == erfc 解析解", err < 2e-3,
          f"最大绝对偏差 {err:.2e}（浓度归一到 1）")

    # [解析] 扩散前沿 ∝ √t
    def front(tt):
        xx, c = agid.diffuse_numeric(L_cm=1.0, n=401, t_s=tt, D=D, C0=1.0)
        return float(np.interp(-0.5, -c, xx))  # c 单调下降，找 c=0.5 的位置

    f1, f4 = front(1800.0), front(7200.0)  # 时间 ×4 -> 前沿应 ×2
    check("[解析] 扩散前沿位置 ∝ √t", abs(f4 / f1 - 2.0) < 0.03,
          f"t×4 后前沿比 {f4 / f1:.4f}（理论 2）")

    # [解析] 等价带窗函数在 r=1 取最大，两侧衰减
    # 等价点在 ν·A = B，即 ν=2、B=1 时 A=0.5
    w_eq = float(agid.equivalence_window(np.array([0.5]), np.array([1.0]), 2.0, 0.9)[0])
    w_ag = float(agid.equivalence_window(np.array([5.0]), np.array([1.0]), 2.0, 0.9)[0])
    w_ab = float(agid.equivalence_window(np.array([0.05]), np.array([1.0]), 2.0, 0.9)[0])
    check("[解析] 等价带窗函数在化学计量等价点最大、抗原/抗体过量两侧衰减",
          w_eq > 0.99 and w_ag < 0.1 and w_ab < 0.1,
          f"νA/B=1 时 {w_eq:.3f} · 抗原过量(νA/B=10) {w_ag:.4f} · 抗体过量(νA/B=0.1) {w_ab:.4f}")

    # [解析] 沉淀线：确实形成一条局域的带，而不是铺满整个凝胶
    s = agid.AgidSetup()
    r = agid.simulate_agid(s, t_end_s=12 * 3600.0, n_out=12)
    sd, pos = r.band_sharpness(), r.band_position()
    check("[解析] 形成局域沉淀带（宽度远小于孔间距）",
          sd == sd and sd < 0.25 * s.L_cm and 0 < pos < s.L_cm,
          f"带心 {pos:.3f} cm · 带宽 SD {sd:.3f} cm · 孔间距 {s.L_cm} cm")

    # [解析] 沉淀带随时间向抗体侧推进（抗原扩散更快）
    p_early, p_late = r.band_position(3), r.band_position(-1)
    check("[解析] 快扩散组分把沉淀带推向慢扩散一侧",
          p_late > p_early,
          f"{r.t[3] / 3600:.1f}h {p_early:.3f} cm -> {r.t[-1] / 3600:.1f}h {p_late:.3f} cm "
          f"(D_ag={s.D_ag:.1e} > D_ab={s.D_ab:.1e})")

    # [解析] 双向扩散：沉淀线偏向**较稀**的一侧（Ouchterlony 教科书行为）
    #        注意这里不该出现前带/后带 —— 两孔是恒浓度储库，凝胶内浓度比
    #        从 ∞ 连续扫到 0，等价带总能找到位置，线只会移动不会消失。
    tc = agid.titration_curve(agid.AgidSetup(B_well=1.0),
                              [0.1, 0.3, 1.0, 3.0, 10.0], t_end_s=6 * 3600.0)
    pos = [d["band_x"] for d in tc]
    check("[解析] 抗原越浓，沉淀线越靠近抗体孔（偏向较稀一侧）",
          all(pos[i] < pos[i + 1] for i in range(len(pos) - 1)),
          "A/B=0.1→10 时带心 " + " → ".join(f"{p:.3f}" for p in pos) + " cm")

    # [解析] 试管沉淀（封闭体系）才有 Heidelberger-Kendall 钟形曲线
    amounts = [0.005, 0.02, 0.1, 0.5, 2.0, 10.0, 50.0]
    hk = agid.heidelberger_kendall_curve(amounts, B_total=1.0, nu=2.0)
    ppt = [d["precipitate"] for d in hk]
    imax = int(np.argmax(ppt))
    check("[解析] 试管沉淀呈 Heidelberger-Kendall 钟形（前带/后带）",
          0 < imax < len(ppt) - 1 and ppt[imax] > 5 * ppt[0] and ppt[imax] > 5 * ppt[-1],
          f"峰在 A/B={hk[imax]['ratio']:g}；沉淀量 " +
          " ".join(f"{v:.3g}" for v in ppt))

    # [解析] 钟形曲线的峰必须落在化学计量等价点 A/B = 1/ν 附近
    check("[解析] 沉淀峰位于化学计量等价点 A/B = 1/ν",
          abs(math.log10(hk[imax]["ratio"] * 2.0)) < 0.4,
          f"峰 A/B={hk[imax]['ratio']:g}，理论 1/ν={1 / 2.0:g}"
          f"（等价比 r={hk[imax]['equivalence_r']:.2f}，理论 1）")

    # [解析] 抗原过量端：抗体被消耗进可溶复合物而非沉淀
    tail = hk[-1]
    check("[解析] 抗原过量时抗体主要进入可溶复合物（后带的机制）",
          tail["soluble"] > 5 * tail["precipitate"],
          f"A/B={tail['ratio']:g}: 可溶 {tail['soluble']:.3g} vs 沉淀 {tail['precipitate']:.3g}")


def validate_agglutination() -> None:
    print("\nM7 平板凝集 Smoluchowski 聚集")
    from . import agglutination as ag

    K, n0, t_end = 1e-3, 1.0e3, 6.0
    r = ag.simulate_aggregation(K, n0, t_end, k_max=90, n_out=30)

    # [解析] 截断误差先确认足够小，后面的比对才有意义
    loss = float(r.mass_loss[-1])
    check("[解析] 截断质量流失可忽略（k_max=90）", loss < 1e-6,
          f"末时刻流失 {loss:.2e}")

    # [解析] 总粒子数 vs 闭式解 N(t)=n0/(1+Kn0t/2)
    ana_N = np.array([ag.total_particles_analytic(t, K, n0) for t in r.t])
    relN = float(np.max(np.abs(r.N_total - ana_N) / ana_N))
    check("[解析] 总粒子数 == 闭式解 n₀/(1+Kn₀t/2)", relN < 1e-6,
          f"最大相对偏差 {relN:.2e}（{len(r.t)} 个时间点）")

    # [解析] 逐个簇尺寸的分布 vs 闭式解
    worst, wk = 0.0, 0
    for ti in (5, 15, 29):
        ana = ag.smoluchowski_analytic(r.k[:30], r.t[ti], K, n0)
        rel = np.abs(r.n[ti, :30] - ana) / np.maximum(ana, 1e-12)
        if rel.max() > worst:
            worst, wk = float(rel.max()), int(np.argmax(rel)) + 1
    check("[解析] 簇尺寸分布 n_k(t) == 闭式解", worst < 1e-5,
          f"最大相对偏差 {worst:.2e}（出现在 k={wk}）")

    # [解析] 质量守恒 Σk·n_k = n₀
    relM = float(np.max(np.abs(r.mass - n0) / n0))
    check("[解析] 质量守恒 Σ k·n_k = n₀", relM < 1e-6,
          f"最大相对偏差 {relM:.2e}")

    # [解析] 桥联效率在 θ=0.5 取极大
    ths = np.linspace(0.01, 0.99, 99)
    effs = [ag.bridging_efficiency(t) for t in ths]
    imax = int(np.argmax(effs))
    check("[解析] 桥联效率 4θ(1−θ) 在 θ=0.5 取最大",
          abs(ths[imax] - 0.5) < 0.02 and abs(effs[imax] - 1.0) < 1e-3,
          f"峰在 θ={ths[imax]:.2f}，值 {effs[imax]:.4f}")

    # [解析] 最优抗体浓度 = Kd
    kd = 7.5
    check("[解析] 使桥联最大的抗体浓度 == Kd",
          abs(ag.epitope_occupancy(ag.optimal_ab_concentration(kd), kd) - 0.5) < 1e-12,
          f"[Ab]=Kd={kd} 时 θ={ag.epitope_occupancy(kd, kd):.6f}")

    # [解析] 抗体滴定呈钩状：两端弱、中间强，峰在 Kd 附近
    concs = np.logspace(-3, 3, 25) * kd
    pc = ag.prozone_curve(concs, kd=kd, K0=2e-3, n0=1e3, t_read_s=6.0,
                          k_max=40, min_size=4)
    vis = [d["visible"] for d in pc]
    ip = int(np.argmax(vis))
    check("[解析] 抗体滴定出现前带（钩状）：高抗体端凝集塌陷",
          0 < ip < len(vis) - 1 and vis[-1] < 0.3 * vis[ip] and vis[0] < 0.3 * vis[ip],
          f"峰在 [Ab]/Kd={pc[ip]['ab'] / kd:.3g}（θ={pc[ip]['theta']:.2f}），"
          f"可见分数 低端 {vis[0]:.3f} / 峰 {vis[ip]:.3f} / 高端 {vis[-1]:.3f}")

    # [解析] 单调性：聚集过程中总粒子数只减不增
    check("[解析] 总粒子数随时间单调下降",
          bool(np.all(np.diff(r.N_total) <= 1e-9)),
          f"N: {r.N_total[0]:.1f} -> {r.N_total[-1]:.1f}")


def main() -> int:
    print("=" * 72)
    print("assaysim 验证套件")
    print("=" * 72)
    validate_nn_thermo()
    validate_neutralization()
    validate_viral_dynamics()
    validate_agid()
    validate_agglutination()
    validate_bridge()
    print("\n" + "=" * 72)
    if _FAILURES:
        print(f"{len(_FAILURES)}/{_CHECKS} 项失败: " + ", ".join(_FAILURES))
        return 1
    print(f"全部 {_CHECKS} 项通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
