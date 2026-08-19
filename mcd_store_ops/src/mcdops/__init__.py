"""mcd_store_ops —— 一家快餐门店的日常运营报表 / 指标矩阵 / 看板 package。

典型用法：

    from datetime import date
    from mcdops.metrics import MetricEngine, load_metrics, load_profile
    from mcdops.simulate import simulate_period
    from mcdops import matrix, report

    profile = load_profile()
    data = simulate_period(profile, date(2026, 8, 18), days=14)
    engine = MetricEngine(data, profile, load_metrics())

    print(report.build_console_brief(engine, date(2026, 8, 18)))
"""

__version__ = "0.1.0"
