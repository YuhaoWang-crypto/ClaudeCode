"""
跨实现交叉验证 — assaysim (Python) vs av_demos_bundle 的 VC-NEU / VC-AGD (浏览器 JS)

为什么值得做
------------
两套代码是**独立写的**：不同语言、不同 lgamma 实现（我用 Lanczos g=7，
他们用 Numerical Recipes 系数）、不同求根方式、不同单位约定。
如果两边在同一批输入上逐位吻合，那么剩下的共同错误只可能是**模型本身**的错，
而不是某一次实现手滑。这比任何一方的自我验证都强。

用法
----
    python3 -m assaysim.crossvalidate_av /path/to/av_demos_bundle/

需要 playwright + 一个 chromium。找不到浏览器时会跳过并说明原因。
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

from . import agid
from .neutralization import Virion, nt50_from_kd

# (Kd nM, n, k) —— 含 ASFV p72 的真实化学计量
NEUT_GRID = [
    (5.0, 24, 1), (1.0, 100, 1), (10.0, 340, 1), (2.5, 50, 7),
    (0.5, 200, 20), (7.3, 14, 1), (100.0, 25, 3), (0.01, 600, 60),
    (50.0, 277, 1), (3.0, 2760, 1), (3.0, 8280, 1),
]
MW_GRID = [14300, 66500, 73179, 150000, 219537, 900000]


def find_chromium() -> str | None:
    for pat in ("/opt/pw-browsers/chromium-*/chrome-linux/chrome",
                "/opt/pw-browsers/chromium/chrome-linux/chrome"):
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    return None


def main(bundle: str) -> int:
    root = Path(bundle).resolve()
    neut, agd = root / "av_neut.html", root / "av_agid.html"
    for f in (neut, agd):
        if not f.exists():
            print(f"找不到 {f}")
            return 2

    exe = find_chromium()
    if exe is None:
        print("未找到 chromium，跳过跨实现验证（assaysim 自身的 44 项验证不受影响）")
        return 0

    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("未安装 playwright，跳过跨实现验证")
        return 0

    worst_neut = worst_rh = worst_d = 0.0
    with sync_playwright() as pw:
        b = pw.chromium.launch(executable_path=exe, args=["--no-sandbox"])
        pg = b.new_page()

        # --- VC-NEU vs M3 占据模型 ---
        pg.goto(neut.as_uri())
        pg.wait_for_timeout(2000)
        print("VC-NEU (他们, JS) vs assaysim M3 (我, Python) — 多击中占据模型")
        print(f"  {'Kd(nM)':>8} {'n':>6} {'k':>5} | {'他们':>15} {'我':>15} {'相对偏差':>10}")
        for kd, n, k in NEUT_GRID:
            t = pg.evaluate("a=>nt50(a[0],a[1],a[2])", [kd, n, k])
            m = nt50_from_kd(kd, Virion("v", n, k))
            rel = abs(t - m) / m
            worst_neut = max(worst_neut, rel)
            tag = "  ← ASFV p72 真实化学计量" if n in (2760, 8280) else ""
            print(f"  {kd:8.2f} {n:6d} {k:5d} | {t:15.8g} {m:15.8g} {rel:10.2e}{tag}")
        print(f"  最大相对偏差 = {worst_neut:.3e}\n")

        # --- VC-AGD vs M5 的扩散层 ---
        # 单位对齐：他们 mm²/day, 298.15 K, η=1e-3 Pa·s, 凝胶因子 0.75, shape 1.20
        pg.goto(agd.as_uri())
        pg.wait_for_timeout(2000)
        print("VC-AGD (他们) vs assaysim M5 (我) — Stokes-Einstein 扩散层（单位与常数已对齐）")
        print(f"  {'MW(Da)':>9} | {'他们 R_h':>10} {'我 R_h':>10} {'偏差':>9} | "
              f"{'他们 D':>10} {'我 D':>10} {'偏差':>9}")
        for mw in MW_GRID:
            tr = pg.evaluate("m=>Rh(m)", mw)
            mr = agid.hydrodynamic_radius_from_MW(mw, 1.20)
            td = pg.evaluate("m=>Dcoef(m)", mw)
            md = agid.D_from_MW(mw, 1.20, T_K=298.15, eta_Pa_s=1e-3,
                                gel_factor=0.75) * 1e2 * 86400  # cm²/s -> mm²/day
            er, ed = abs(tr - mr) / mr, abs(td - md) / md
            worst_rh, worst_d = max(worst_rh, er), max(worst_d, ed)
            tag = "  ← p72 单体" if mw == 73179 else ("  ← p72 三聚体" if mw == 219537 else "")
            print(f"  {mw:9,d} | {tr:10.4f} {mr:10.4f} {er:9.2e} | "
                  f"{td:10.4f} {md:10.4f} {ed:9.2e}{tag}")
        print(f"  最大偏差: R_h {worst_rh:.2e}   D {worst_d:.2e}")

        # --- VC-VK vs M4 靶细胞受限 ODE (Baccam 2006 参数) ---
        vk = next(iter(sorted(Path(bundle).rglob("*vc_vk.html"))), None)
        worst_vk = None
        if vk is not None:
            import numpy as np

            from . import viral_dynamics as vd

            pg.goto(vk.as_uri())
            pg.wait_for_timeout(2000)
            bac = vd.REFERENCE_PARAMS["baccam2006_H1N1"]
            V0 = 7.5e-2
            worst_vk = 0.0
            print("\nVC-VK (他们) vs assaysim M4 (我) — Baccam 2006 靶细胞受限 ODE")
            print("  ⚠️ 只有 absorption=False 才能吃这套文献参数；含吸附项时 R0 会掉到 0.009")
            print(f"  {'场景':<24} {'他们 峰值':>12} {'我 峰值':>12} {'偏差':>9}")
            for label, C, ec, ts in [("无药", None, None, None),
                                     ("100nM / EC50 20nM", 100.0, 20.0, 1.0),
                                     ("10nM / EC50 100nM", 10.0, 100.0, 1.0),
                                     ("1µM / EC50 50nM, t=2d", 1000.0, 50.0, 2.0)]:
                if C is None:
                    js = ("()=>{const s=VCVK.simulate(null,null,8.0,0.0005);"
                          "return VCVK.metrics(s.t,s.V,1.0).peak_titer;}")
                    tr = vd.simulate(bac, moi=V0 / bac.T0, t_end_h=8.0, n_points=16001)
                else:
                    js = (f"()=>{{const e=VCVK.constantEpsFn({C},{ec},{ts});"
                          f"const s=VCVK.simulate(null,e,8.0,0.0005);"
                          f"return VCVK.metrics(s.t,s.V,1.0).peak_titer;}}")
                    d = vd.Drug("d", moa="replication", ec50_mol=ec, hill=1.0)
                    tr = vd.simulate(bac, moi=V0 / bac.T0, t_end_h=8.0, n_points=16001,
                                     agent=d, dose=C, dose_start_h=ts)
                t = pg.evaluate(js)
                m = tr.peak_titer()
                rel = abs(t - m) / m
                worst_vk = max(worst_vk, rel)
                print(f"  {label:<24} {t:12.4e} {m:12.4e} {rel:9.2e}")
            print(f"  最大峰值相对偏差 = {worst_vk:.2e}")
            print("  （延迟给药一档偏差较大是因为峰值恰落在给药阶跃处，"
                  "定步长 RK4 与变步长 LSODA 对间断的处理不同）")

        b.close()

    ok = (worst_neut < 1e-9 and worst_rh < 1e-6 and worst_d < 1e-6
          and (worst_vk is None or worst_vk < 2e-3))
    print("\n" + ("跨实现一致 ✅" if ok else "存在超出容差的差异 ⚠️  —— 需逐项排查约定差异"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
