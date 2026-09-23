"""指标实现的自检。

写这个文件的直接原因：BEDROC 的归一化因子最初写反了（`rie / fac` 而非 `rie * fac`），
完美排序下返回 17.0 而不是 1.0。这类错误不会抛异常，只会让报告里的数字悄悄失真，
所以用"已知输入 -> 已知输出"的方式钉住。

    python -m mers_pipeline.test_metrics
"""
from __future__ import annotations

import numpy as np

from .a5_docking import bedroc, enrichment_factor


def test_bedroc_bounds_and_extremes() -> None:
    n, n_act = 79, 19
    y = np.zeros(n, dtype=int)
    y[:n_act] = 1
    perfect = np.concatenate([np.linspace(10, 9, n_act), np.linspace(1, 0, n - n_act)])

    assert abs(bedroc(y, perfect, 20.0) - 1.0) < 1e-6, "完美排序的 BEDROC 应为 1"
    assert abs(bedroc(y, -perfect, 20.0) - 0.0) < 1e-6, "最差排序的 BEDROC 应为 0"

    # 随机排序的 BEDROC 期望值 ≈ 活性物比例 Ra
    rng = np.random.default_rng(0)
    vals = np.array([bedroc(y, rng.normal(size=n), 20.0) for _ in range(3000)])
    assert np.all((vals >= 0) & (vals <= 1)), "BEDROC 必须落在 [0, 1]"
    assert abs(vals.mean() - n_act / n) < 0.02, "随机排序的 BEDROC 均值应接近 Ra"
    print("✓ BEDROC: 完美=1, 最差=0, 随机均值 %.3f ≈ Ra %.3f，全部落在 [0,1]"
          % (vals.mean(), n_act / n))


def test_enrichment_factor() -> None:
    n, n_act = 79, 19
    y = np.zeros(n, dtype=int)
    y[:n_act] = 1
    perfect = np.concatenate([np.linspace(10, 9, n_act), np.linspace(1, 0, n - n_act)])

    # 前 10% (k=8) 全是活性物 -> EF = 1 / Ra
    assert abs(enrichment_factor(y, perfect, 0.10) - n / n_act) < 1e-6
    # n=79 时 EF1% 只覆盖 1 个化合物，应拒绝给出数字而不是给个退化值
    assert np.isnan(enrichment_factor(y, perfect, 0.01)), "k 过小时 EF 应返回 NaN"
    print("✓ EF: 完美排序 EF10%% = %.2f = 1/Ra；EF1%% 在 n=79 时正确返回 NaN"
          % enrichment_factor(y, perfect, 0.10))


def test_enrichment_factor_random_is_one() -> None:
    n, n_act = 400, 100
    y = np.zeros(n, dtype=int)
    y[:n_act] = 1
    rng = np.random.default_rng(1)
    vals = [enrichment_factor(y, rng.normal(size=n), 0.10) for _ in range(2000)]
    assert abs(np.mean(vals) - 1.0) < 0.05, "随机排序的 EF 期望值应为 1"
    print("✓ EF: 随机排序 EF10%% 均值 %.3f ≈ 1" % np.mean(vals))


if __name__ == "__main__":
    test_bedroc_bounds_and_extremes()
    test_enrichment_factor()
    test_enrichment_factor_random_is_one()
    print("\n全部指标自检通过。")
