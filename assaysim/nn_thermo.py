"""
M1 — 最近邻 (nearest-neighbour) 双链热力学

实现 Allawi & SantaLucia (1997) / SantaLucia (1998) 的"统一"NN 参数集，
计算 DNA 双链的 ΔH°, ΔS°, ΔG°37 与解链温度 Tm。

参数来源
--------
Allawi HT, SantaLucia J (1997) Biochemistry 36:10581  — 统一 NN 参数
SantaLucia J (1998) PNAS 95:1460                      — 统一参数汇总
SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415  — 综述

盐校正
------
Owczarzy R et al. (2004) Biochemistry 43:3537 — 1/Tm 对 ln[Na+] 的经验式

接口
----
    duplex_thermo(seq)               -> Thermo(dH, dS, dG37, ...)
    tm(seq, ct_nM, na_mM, mg_mM)     -> float, °C
    mismatch_penalty_dG(...)         -> float, kcal/mol  (给 M2 用)

单位约定 (全局一致，验证里会检查)
    dH   kcal/mol
    dS   cal/(mol·K)
    dG   kcal/mol
    Tm   °C
"""

from __future__ import annotations

import math
from dataclasses import dataclass

R_CAL = 1.98720425864083  # cal/(mol·K)
T37 = 310.15  # K

# --- Allawi & SantaLucia 1997 统一最近邻参数 (1 M NaCl) -------------------
# key = 5'->3' 二核苷酸 (顶链)；ΔH kcal/mol, ΔS cal/(mol·K)
NN_UNIFIED: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# 起始项 (SantaLucia 1998)
INIT_GC = (0.1, -2.8)  # 末端为 G 或 C
INIT_AT = (2.3, 4.1)  # 末端为 A 或 T
SYMMETRY_DS = -1.4  # 自互补分子的对称性校正 (cal/mol/K)

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def is_self_complementary(seq: str) -> bool:
    return seq.upper() == revcomp(seq.upper())


@dataclass(frozen=True)
class Thermo:
    """一个双链的热力学量。dH kcal/mol, dS cal/(mol·K), dG kcal/mol, tm °C."""

    seq: str
    dH: float
    dS: float
    dG37: float
    tm: float
    gc_frac: float
    length: int

    def __str__(self) -> str:
        return (
            f"{self.seq}  len={self.length}  GC={self.gc_frac:.1%}  "
            f"dH={self.dH:.2f} kcal/mol  dS={self.dS:.2f} cal/mol/K  "
            f"dG37={self.dG37:.2f} kcal/mol  Tm={self.tm:.2f} C"
        )


def _nn_sum(seq: str) -> tuple[float, float]:
    """对完全互补双链累加 NN 项 + 起始项 + 对称校正。返回 (dH, dS)。"""
    s = seq.upper()
    if len(s) < 2:
        raise ValueError("序列至少需要 2 nt")
    bad = set(s) - set("ACGT")
    if bad:
        raise ValueError(f"含非 ACGT 字符: {sorted(bad)}")

    dH = dS = 0.0
    for i in range(len(s) - 1):
        h, sdt = NN_UNIFIED[s[i : i + 2]]
        dH += h
        dS += sdt

    # 两端各加一个起始项
    for end in (s[0], s[-1]):
        h, sdt = INIT_GC if end in "GC" else INIT_AT
        dH += h
        dS += sdt

    if is_self_complementary(s):
        dS += SYMMETRY_DS

    return dH, dS


def salt_correction_owczarzy(tm_1M_K: float, gc_frac: float, n_bp: int,
                             na_mM: float) -> float:
    """Owczarzy 2004 单价盐校正。输入 1 M Na+ 下的 Tm(K)，返回校正后的 Tm(K)。

    1/Tm(Na) = 1/Tm(1M) + (4.29*fGC - 3.95)e-5 * ln(Na) + 9.40e-6 * ln(Na)^2
    """
    if na_mM <= 0:
        raise ValueError("Na+ 浓度必须为正")
    ln_na = math.log(na_mM / 1000.0)
    inv = (
        1.0 / tm_1M_K
        + (4.29 * gc_frac - 3.95) * 1e-5 * ln_na
        + 9.40e-6 * ln_na * ln_na
    )
    return 1.0 / inv


def duplex_thermo(seq: str, ct_nM: float = 250.0, na_mM: float = 50.0,
                  salt_correct: bool = True) -> Thermo:
    """完全互补双链的热力学。

    ct_nM: 总链浓度 (nM)。非自互补时有效浓度取 Ct/4 (两条链等量的惯例)。
    na_mM: 单价阳离子浓度 (mM)。salt_correct=False 时返回 1 M Na+ 的 Tm。
    """
    s = seq.upper()
    dH, dS = _nn_sum(s)

    ct = ct_nM * 1e-9
    x = 1.0 if is_self_complementary(s) else 4.0

    # Tm = dH / (dS + R ln(Ct/x))；dH 由 kcal 转 cal
    tm_K = (dH * 1000.0) / (dS + R_CAL * math.log(ct / x))

    gc = (s.count("G") + s.count("C")) / len(s)
    if salt_correct:
        tm_K = salt_correction_owczarzy(tm_K, gc, len(s), na_mM)

    dG37 = dH - T37 * dS / 1000.0

    return Thermo(seq=s, dH=dH, dS=dS, dG37=dG37, tm=tm_K - 273.15,
                  gc_frac=gc, length=len(s))


def tm(seq: str, ct_nM: float = 250.0, na_mM: float = 50.0) -> float:
    return duplex_thermo(seq, ct_nM, na_mM).tm


# --- 错配惩罚 (给 M2 in-silico PCR 用) -----------------------------------
#
# 完整的内部错配 NN 参数需要 Allawi & SantaLucia 1997-1998 四篇论文的全部表格。
# 这里用一个位置加权的经验惩罚：错配的自由能代价本身在 0.5-4 kcal/mol 量级，
# 但对 PCR 而言真正决定成败的是错配离 3' 端的距离 —— 聚合酶延伸需要 3' 端配对。
#
# 该惩罚函数是 ⚠️ 经验模型，不是实测热力学。它的用途是给引物-模板对做
# **排序**，validate.py 只验证它的单调性与边界行为，不声称它给出真实 ΔΔG。

MISMATCH_BASE_DG = 3.0  # kcal/mol, 内部错配的典型代价量级


def mismatch_penalty_dG(mismatch_positions: list[int], primer_len: int,
                        three_prime_weight: float = 6.0,
                        decay_nt: float = 4.0) -> float:
    """错配的经验自由能惩罚 (kcal/mol, 正值 = 不利)。

    mismatch_positions: 0-based，从引物 5' 端计数。
    引物 3' 末端 (position = primer_len-1) 的错配惩罚最重，向 5' 端指数衰减。
    """
    total = 0.0
    for pos in mismatch_positions:
        dist_from_3p = (primer_len - 1) - pos
        weight = 1.0 + (three_prime_weight - 1.0) * math.exp(-dist_from_3p / decay_nt)
        total += MISMATCH_BASE_DG * weight
    return total
