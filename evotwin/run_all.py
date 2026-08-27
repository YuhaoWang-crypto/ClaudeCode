"""EvoTwin 位置/结构/遗传图谱扩展 —— 四个模块一次跑完。

    python3 -m evotwin.run_all
"""
from . import e1_target_supply, e2_liquidity, e3_recombination, e4_avoidance


def main():
    for mod in (e1_target_supply, e2_liquidity, e3_recombination, e4_avoidance):
        mod.run()
        print()


if __name__ == "__main__":
    main()
