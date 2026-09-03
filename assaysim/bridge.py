"""
M6 — "非细胞 -> 细胞" 桥接模型与它的诚实误差界

要回答的问题
------------
    只测非细胞数据 (纯酶 IC50 / SPR Kd)，能把原本必须上细胞才能得到的
    EC50 预测到什么精度？

做法
----
1. 从 ChEMBL 取同一分子的 (生化 pIC50, 细胞 pEC50) 真实配对
2. 用**骨架划分** (Murcko scaffold) 切训练/测试集 —— 随机划分会把同一
   化学系列的类似物拆到两边，导致乐观偏差
3. 比三个预测器：
      identity  : p_cell = p_biochem          ("生化就当细胞用")
      offset    : p_cell = p_biochem - median(Δ)
      linear    : p_cell = a·p_biochem + b
4. 把误差和**测量本身的噪声下限**对比 —— 任何预测器都不可能优于目标量
   自身的重复性

噪声下限
--------
同一分子在不同论文里被反复测量。这些重复值的离散度就是该终点的
真实重现性；预测误差低于它是不可能的，接近它就说明模型已经用尽了
非细胞数据里的信息。
"""

from __future__ import annotations

import math
import random
import statistics as st
from dataclasses import dataclass

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:  # pragma: no cover
    HAVE_RDKIT = False


# --- 划分 ----------------------------------------------------------------

def murcko_scaffold(smiles: str | None) -> str | None:
    if not HAVE_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return None


def scaffold_split(pairs: list[dict], test_frac: float = 0.3,
                   seed: int = 0) -> tuple[list[dict], list[dict], str]:
    """按 Murcko 骨架分组划分：同一骨架的分子整组进训练或测试。

    这样测试集里不会出现训练集分子的近似物，是化学信息学的标准做法。
    RDKit 不可用时回退到随机划分，并在返回值里标明。
    """
    groups: dict[str, list[dict]] = {}
    n_no_scaffold = 0
    for p in pairs:
        s = murcko_scaffold(p.get("smiles"))
        if s is None:
            n_no_scaffold += 1
            s = f"__none_{n_no_scaffold}"
        groups.setdefault(s, []).append(p)

    mode = "scaffold" if HAVE_RDKIT else "random(no-rdkit)"
    keys = sorted(groups)
    rng = random.Random(seed)
    rng.shuffle(keys)

    target = int(round(test_frac * len(pairs)))
    test: list[dict] = []
    for k in keys:
        if len(test) >= target:
            break
        test.extend(groups[k])
    test_ids = {id(x) for x in test}
    train = [p for p in pairs if id(p) not in test_ids]
    return train, test, f"{mode}, {len(groups)} 组"


# --- 预测器 --------------------------------------------------------------

@dataclass(frozen=True)
class Bridge:
    """p_cell_pred = slope * p_biochem + intercept"""

    slope: float
    intercept: float
    n_train: int
    kind: str

    def predict(self, p_biochem):
        x = np.asarray(p_biochem, dtype=float)
        return self.slope * x + self.intercept


def fit_identity(train: list[dict]) -> Bridge:
    return Bridge(1.0, 0.0, len(train), "identity")


def fit_offset(train: list[dict]) -> Bridge:
    med = st.median([p["delta"] for p in train])
    return Bridge(1.0, -med, len(train), "offset")


def fit_linear(train: list[dict]) -> Bridge:
    x = np.array([p["p_biochem"] for p in train])
    y = np.array([p["p_cell"] for p in train])
    slope, intercept = np.polyfit(x, y, 1)
    return Bridge(float(slope), float(intercept), len(train), "linear")


# --- 评估 ----------------------------------------------------------------

def evaluate(bridge: Bridge, test: list[dict]) -> dict:
    x = np.array([p["p_biochem"] for p in test])
    y = np.array([p["p_cell"] for p in test])
    yhat = bridge.predict(x)
    err = yhat - y

    ss_res = float(np.sum(err ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))

    return {
        "kind": bridge.kind,
        "n_test": len(test),
        "rmse": float(np.sqrt(np.mean(err ** 2))),
        "mae": float(np.mean(np.abs(err))),
        "bias": float(np.mean(err)),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "pearson": float(np.corrcoef(x, y)[0, 1]),
        "within_0.5log": float(np.mean(np.abs(err) <= 0.5)),
        "within_1.0log": float(np.mean(np.abs(err) <= 1.0)),
        "slope": bridge.slope,
        "intercept": bridge.intercept,
    }


def noise_floor(pairs: list[dict], side: str = "cell", min_n: int = 3) -> dict:
    """用文献重复测量估计该终点自身的重现性 (log10 单位)。

    ⚠️ 这是**跨论文**重现性，包含实验室间差异、细胞系差异、方案差异，
       因此它是"同一分子在不同报道间的离散度"，不是单一实验室的精度。
       但它正是任何跨论文训练的预测器所面对的噪声上限。
    """
    key = "spread_cell" if side == "cell" else "spread_biochem"
    nkey = "n_cell" if side == "cell" else "n_biochem"

    spreads = [p[key] for p in pairs if p[nkey] >= min_n]
    if not spreads:
        return {"n": 0}

    # 极差 -> SD 的换算依赖 n；用无量纲的稳健做法：报告极差的中位数，
    # 并给出把极差近似为 ~2.5 SD (n≈5 时的期望) 的粗略 SD 估计
    med_range = st.median(spreads)
    return {
        "n": len(spreads),
        "median_range_log10": med_range,
        "approx_sd_log10": med_range / 2.5,
        "note": "极差 -> SD 用 n≈5 的期望因子 2.5 换算，仅为量级估计",
    }


def run_analysis(pairs: list[dict], *, test_frac: float = 0.3,
                 seed: int = 0) -> dict:
    """完整跑一遍：划分 -> 拟合三个预测器 -> 留出评估 -> 噪声下限。"""
    train, test, split_desc = scaffold_split(pairs, test_frac, seed)

    results = []
    for fitter in (fit_identity, fit_offset, fit_linear):
        b = fitter(train)
        results.append(evaluate(b, test))

    return {
        "n_pairs": len(pairs),
        "n_train": len(train),
        "n_test": len(test),
        "split": split_desc,
        "results": results,
        "noise_cell": noise_floor(pairs, "cell"),
        "noise_biochem": noise_floor(pairs, "biochem"),
    }


def format_analysis(name: str, a: dict) -> str:
    lines = [
        f"=== {name} ===",
        f"配对 N={a['n_pairs']}  训练 {a['n_train']} / 测试 {a['n_test']}  划分: {a['split']}",
        "",
        f"  {'预测器':<10} {'RMSE':>7} {'MAE':>7} {'bias':>7} {'R²':>7} {'r':>6} {'±0.5log':>8} {'±1log':>7}",
    ]
    for r in a["results"]:
        lines.append(
            f"  {r['kind']:<10} {r['rmse']:7.3f} {r['mae']:7.3f} {r['bias']:+7.3f} "
            f"{r['r2']:7.3f} {r['pearson']:6.3f} {r['within_0.5log']:8.1%} {r['within_1.0log']:7.1%}"
        )
    nc, nb = a["noise_cell"], a["noise_biochem"]
    lines.append("")
    if nc.get("n"):
        lines.append(
            f"  噪声下限 (细胞法, n>=3 的 {nc['n']} 个分子): "
            f"文献极差中位数 {nc['median_range_log10']:.2f} log  -> SD ≈ {nc['approx_sd_log10']:.2f} log"
        )
    if nb.get("n"):
        lines.append(
            f"  噪声下限 (生化法, n>=3 的 {nb['n']} 个分子): "
            f"文献极差中位数 {nb['median_range_log10']:.2f} log  -> SD ≈ {nb['approx_sd_log10']:.2f} log"
        )
    return "\n".join(lines)
