"""共享配置：路径、ChEMBL 靶标、PDB 受体定义。

所有模块从这里取路径，保证结果目录结构一致。
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PDB_DIR = DATA / "pdb"
CHEMBL_DIR = DATA / "chembl"
SCIPLEX_DIR = DATA / "sciplex"
RESULTS = ROOT / "results"
RES_A = RESULTS / "trackA"
RES_B = RESULTS / "trackB"
FIGS = ROOT / "figures" / "mers"

for _d in (DATA, PDB_DIR, CHEMBL_DIR, SCIPLEX_DIR, RESULTS, RES_A, RES_B, FIGS):
    _d.mkdir(parents=True, exist_ok=True)

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

# ---------------------------------------------------------------- 靶标定义
TARGETS = {
    # MERS-CoV replicase pp1ab —— Mpro(3CLpro) 与 PLpro 的活性都挂在这个靶标下，
    # 需要按 assay_description 再拆分。
    "MERS_PP1AB": "CHEMBL4295557",
    # MERS-CoV 细胞水平（Vero / Huh-7 等）抗病毒数据
    "MERS_CELL": "CHEMBL4296578",
    # 跨冠状病毒迁移训练集
    "SARS2_MPRO": "CHEMBL4523582",
    "SARS_PP1AB_1": "CHEMBL3927",
    "SARS_PP1AB_2": "CHEMBL5118",
}

# assay_description 关键词 -> 酶。顺序重要：PLpro 关键词先匹配。
PLPRO_KEYS = ("papain-like", "papain like", "plpro", "pl protease", "pl-pro")
MPRO_KEYS = ("3cl", "3c-like", "3c like", "mpro", "main protease", "ктс")

# ---------------------------------------------------------------- 受体定义
# box_pad: 以共晶配体重原子外接盒各方向外扩的 Å 数
RECEPTORS = {
    "MPRO_5WKK": dict(
        pdb="5WKK", enzyme="Mpro", resolution=1.55,
        ligand_code=None,   # 运行时从结构里自动挑选最大的非聚合物配体
        chain=None, box_pad=5.0,
        note="GC813 共晶，首选 Mpro 受体",
    ),
    "MPRO_5WKJ": dict(
        pdb="5WKJ", enzyme="Mpro", resolution=2.05,
        ligand_code=None, chain=None, box_pad=5.0,
        note="GC376 共晶，备选 Mpro 受体",
    ),
    "PLPRO_4RF1": dict(
        pdb="4RF1", enzyme="PLpro", resolution=2.15,
        ligand_code=None, chain=None, box_pad=5.0,
        note="带小分子配体，运行时核实配体身份与位点",
    ),
    "PLPRO_4RNA": dict(
        pdb="4RNA", enzyme="PLpro", resolution=1.79,
        ligand_code=None, chain=None, box_pad=6.0,
        note="apo 回退受体；若 4RF1 配体不适合 redock 则启用",
    ),
}

# ---------------------------------------------------------------- 运行参数
N_CPU = int(os.environ.get("MERS_NCPU", "4"))
VINA_EXHAUSTIVENESS = int(os.environ.get("MERS_EXH", "8"))
VINA_N_POSES = 5
RANDOM_SEED = 42

# drug-like 过滤（遵循 binding-affinity 建模惯例）
MAX_MW = 650.0
MIN_MW = 150.0


def log(msg: str) -> None:
    print(f"[mers] {msg}", flush=True)


def setup_cjk_fonts() -> str | None:
    """让 matplotlib 能渲染中文标签；返回选中的字体名。"""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager, rcParams

    candidates = [
        "Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Serif CJK SC",
        "WenQuanYi Zen Hei", "Source Han Sans SC", "Droid Sans Fallback",
    ]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            rcParams["axes.unicode_minus"] = False
            return name
    log("警告：未找到中文字体，图中中文可能显示为方框")
    return None
