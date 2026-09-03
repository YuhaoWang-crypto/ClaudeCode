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


def main() -> int:
    print("=" * 72)
    print("assaysim 验证套件")
    print("=" * 72)
    validate_nn_thermo()
    validate_neutralization()
    validate_viral_dynamics()
    validate_bridge()
    print("\n" + "=" * 72)
    if _FAILURES:
        print(f"{len(_FAILURES)}/{_CHECKS} 项失败: " + ", ".join(_FAILURES))
        return 1
    print(f"全部 {_CHECKS} 项通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
