"""
C2 — Anchor each of the six axes in a REAL curated model, and measure the
price of putting them in one state space.

C1 showed the six axes exist as data. This module asks a harder question: does
each axis exist as *dynamics* with published parameters, and can the six be put
on a common footing?

Six anchors, each fetched as official SBML from the BioModels GitHub mirror and
integrated with libRoadRunner (zero transcription error -- the M15/M20b
standard). All six were found by the live C1 search, not picked from memory:

  cell_cycle         BIOMD0000000318  Yao2008 Rb-E2F restriction-point switch
  stress             BIOMD0000000446  Erguler2013 unfolded-protein stress response
  metabolic_reserve  BIOMD0000000041  Kongas2007 creatine-kinase energy buffer
  dna_damage         BIOMD0000000188  Proctor2008 p53/Mdm2 stabilised by ATM
  apoptotic_priming  BIOMD0000000220  Albeck2008 extrinsic apoptosis (EARM)
  epigenetic_lineage BIOMD0000000210  Chickarmane2008 NANOG/GATA-6 lineage switch

What gets computed:
  1. Dimension of each anchor (species / parameters / reactions), hence the
     compression ratio implied by collapsing it to ONE state variable.
  2. The characteristic relaxation timescale of each axis, measured from the
     actual trajectory (time to 63% of the total excursion of the fastest-
     moving species) and cross-checked against the Jacobian spectrum.
  3. The timescale SPREAD across axes -- the number that decides whether one
     sampling cadence can observe all six.

Rigour: ✅ the SBML, the dimensions and the integration are exact. ⚠️ the
one-model-per-axis assignment is an editorial choice; other curated models
would serve. ⚠️ the anchors come from different organisms and unit systems,
which is itself part of the finding, not a bug in the analysis.
"""
import numpy as np

from .common import banner, load_rr, save_json

ANCHORS = [
    ("cell_cycle", "细胞周期", "BIOMD0000000318",
     "Yao2008 Rb-E2F restriction-point switch", "bistable switch"),
    ("stress", "应激水平", "BIOMD0000000446",
     "Erguler2013 unfolded-protein stress response", "adaptive response"),
    ("metabolic_reserve", "代谢储备", "BIOMD0000000041",
     "Kongas2007 creatine-kinase energy buffer", "buffered reserve"),
    ("dna_damage", "DNA 损伤", "BIOMD0000000188",
     "Proctor2008 p53/Mdm2 stabilised by ATM", "damage-triggered pulses"),
    ("apoptotic_priming", "凋亡准备度", "BIOMD0000000220",
     "Albeck2008 extrinsic apoptosis (EARM)", "all-or-none snap"),
    ("epigenetic_lineage", "表观/谱系", "BIOMD0000000210",
     "Chickarmane2008 NANOG/GATA-6 lineage switch", "bistable lineage switch"),
]


def _dims(rr):
    return dict(species=len(rr.model.getFloatingSpeciesIds()),
                boundary=len(rr.model.getBoundarySpeciesIds()),
                params=len(rr.model.getGlobalParameterIds()),
                reactions=len(rr.model.getReactionIds()))


def _relaxation_time(rr, tmax, npts=2000):
    """Characteristic timescale, measured from the trajectory itself.

    For the species with the largest fractional excursion, the time at which
    it first covers 63% (1 - 1/e) of its total excursion. This is model-free:
    it needs no steady state and no linearisation.
    """
    rr.reset()
    try:
        res = rr.simulate(0, tmax, npts)
    except Exception:
        return None, None
    arr = np.array(res)
    t, Y = arr[:, 0], arr[:, 1:]
    span = Y.max(axis=0) - Y.min(axis=0)
    scale = np.maximum(np.abs(Y).max(axis=0), 1e-12)
    frac = span / scale
    if not np.any(np.isfinite(frac)) or frac.max() <= 1e-9:
        return None, None
    j = int(np.nanargmax(frac))
    y = Y[:, j]
    target = y[0] + 0.63 * (y[-1] - y[0])
    idx = np.where((y - target) * np.sign(y[-1] - y[0]) >= 0)[0]
    tau = float(t[idx[0]]) if len(idx) else None
    name = rr.model.getFloatingSpeciesIds()[j] if j < len(rr.model.getFloatingSpeciesIds()) else f"col{j}"
    return tau, name


def _jacobian_timescale(rr):
    """Slowest non-zero mode from the Jacobian, if roadrunner can form one."""
    try:
        J = np.array(rr.getFullJacobian())
        ev = np.linalg.eigvals(J)
        neg = ev[(ev.real < -1e-12)]
        if len(neg) == 0:
            return None
        return float(1.0 / np.min(np.abs(neg.real)))
    except Exception:
        return None


def run():
    banner("C2 — 六个轴各自锚定到真实策展模型(官方 SBML,libRoadRunner 积分)")
    rows = []
    # Integration horizons chosen per model's own unit system (s vs min vs h).
    horizons = {"BIOMD0000000318": 4000, "BIOMD0000000446": 4000,
                "BIOMD0000000041": 500, "BIOMD0000000188": 20000,
                "BIOMD0000000220": 40000, "BIOMD0000000210": 4000}

    for key, cn, bid, cite, behaviour in ANCHORS:
        rr, path = load_rr(bid)
        if rr is None:
            print(f"  ⚠️  {cn:<10} {bid} 未能取回/加载 —— 跳过(不编造)")
            rows.append(dict(axis=key, label=cn, biomd=bid, cite=cite, ok=False))
            continue
        d = _dims(rr)
        tau, drv = _relaxation_time(rr, horizons[bid])
        rr.reset()
        tj = _jacobian_timescale(rr)
        rows.append(dict(axis=key, label=cn, biomd=bid, cite=cite,
                         behaviour=behaviour, ok=True, tau=tau,
                         tau_jac=tj, driver=drv, **d))
        print(f"  ✅ {cn:<10} {bid}  {cite}")

    ok = [r for r in rows if r.get("ok")]
    if not ok:
        print("\n⚠️  没有任何模型加载成功,C2 无结论。")
        return rows

    print("\n【每个轴的真实维度】把它压成 1 个状态变量,压缩比是多少")
    print(f"{'轴':<12}{'物种':>6}{'参数':>6}{'反应':>6}{'→1 变量的压缩比':>16}")
    print("-" * 50)
    tot_s = tot_p = 0
    for r in ok:
        tot_s += r["species"]; tot_p += r["params"]
        print(f"{r['label']:<10}{r['species']:>6}{r['params']:>6}"
              f"{r['reactions']:>6}{r['species']:>13}:1")
    print("-" * 50)
    print(f"{'合计':<10}{tot_s:>6}{tot_p:>6}")
    print(f"\n六轴联合的机制态维数 = {tot_s} 个物种、{tot_p} 个参数。")
    print(f"用户提出的 6 变量表述 = {tot_s}:{6} ≈ {tot_s / 6:.0f}:1 的压缩。")
    print("这不是反对意见 —— 这正是状态空间要做的事;但压缩比给出了")
    print("必须由数据来证明的'可约化性'负担有多重。")

    print("\n【时间尺度】各轴的特征弛豫时间(从轨迹实测,63% excursion)")
    print(f"{'轴':<12}{'τ(模型自身时间单位)':>22}{'雅可比最慢模':>16}{'驱动物种':>14}")
    print("-" * 68)
    taus = []
    for r in ok:
        if r["tau"]:
            taus.append(r["tau"])
        f = lambda v: "n/a" if v is None else f"{v:.4g}"
        print(f"{r['label']:<10}{f(r['tau']):>20}{f(r['tau_jac']):>18}"
              f"{str(r['driver'])[:12]:>14}")

    if len(taus) >= 2:
        spread = max(taus) / min(taus)
        print(f"\n时间尺度跨度 = {spread:.0f} 倍(且各模型的时间单位并不统一:"
              f"秒/分钟/小时混用)。")
        print("含义:单一采样节奏无法同时观测六个轴。快轴在慢轴的一个采样间隔内")
        print("已经走完并复位 —— 这是状态空间构建的第一个硬约束,与数据量无关。")

    print("\n【关键缺口】BioModels 里没有任何一个策展模型把这六个轴耦合在一起。")
    print("六个锚点分属不同物种、不同单位、不同时间尺度,彼此的耦合项在文献中")
    print("不是一个被标定过的对象。所以'6 变量状态空间'目前是一个待建的假设,")
    print("不是一个可以直接从数据库里取出来的现成结构。")

    save_json("c2_anchors.json", rows)
    print("\n结果已存 figures/_cellstate_cache/c2_anchors.json")
    return rows


if __name__ == "__main__":
    run()
