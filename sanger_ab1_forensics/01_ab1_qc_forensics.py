#!/usr/bin/env python3
"""ABI .ab1 质控 + 真实性取证。

对一个 Sanger .ab1 文件给出两类结论：
  (1) 常规 QC —— 读长、QV 分布、信噪比、简并碱基；
  (2) 真实性取证 —— 分析道 (DATA9-12) 是否为真实电泳信号。

真实的毛细管电泳数据必然带有：峰高随迁移时间指数衰减、峰形不对称拖尾、
峰间距随读长漂移、基线有噪音、读段两端 QV 低。合成/重绘的色谱图没有这些性质。

用法: python3 01_ab1_qc_forensics.py <file.ab1>
"""
import sys
import numpy as np
from Bio import SeqIO


def load(path):
    ab = SeqIO.read(path, "abi").annotations["abif_raw"]
    dec = lambda v: v.decode("latin-1") if isinstance(v, bytes) else v
    fwo = dec(ab["FWO_1"])                       # 染料-碱基顺序, 通常 'GATC'
    return ab, dec, fwo


def qc(ab, dec, fwo):
    seq = dec(ab["PBAS1"]).upper()
    qv = np.array([ord(c) for c in dec(ab["PCON1"])])
    loc = np.array(ab["PLOC1"])
    n = len(seq)
    print("=" * 72)
    print("【1】常规质控")
    print("=" * 72)
    print(f"  仪器            {dec(ab['HCFG3'])}   聚合物 {dec(ab['GTyp1']).strip()}   "
          f"化学 {dec(ab['DySN1'])}")
    print(f"  运行日期        {dec(ab['RUND1'])} {dec(ab['RUNT1'])}   泳道 {ab['LANE1']}")
    print(f"  Basecaller      {dec(ab['SVER2'])}")
    print(f"  读长            {n} bp")
    print(f"  QV              mean={qv.mean():.1f}  median={np.median(qv):.0f}  "
          f"min={qv.min()}  max={qv.max()}  唯一取值={len(np.unique(qv))}")
    for t in (20, 30, 40):
        print(f"    QV>={t}        {(qv >= t).sum():4d} ({(qv >= t).sum() / n * 100:5.1f}%)")
    amb = [(i + 1, c) for i, c in enumerate(seq) if c not in "ACGT"]
    print(f"  简并/N 碱基     {len(amb)}  {amb[:20]}")

    sn, nois = ab["S/N%1"], ab["NOIS1"]
    print("  各通道信号/噪音:")
    for i, b in enumerate(fwo):
        print(f"    {b}  signal={sn[i]:6.0f}  noise={nois[i]:5.2f}  S/N={sn[i] / nois[i]:6.0f}")

    sp = np.diff(loc)
    print(f"  峰间距          实测 mean={sp.mean():.2f} sd={sp.std():.3f} "
          f"[{sp.min()}, {sp.max()}]   头文件 SPAC={ab['SPAC1']:.2f}")
    return seq, qv, loc, sp


def forensics(ab, dec, fwo, seq, qv, loc, sp):
    """五项独立检验。每项返回 (是否异常, 说明)。"""
    ana = {fwo[i]: np.array(ab[f"DATA{9 + i}"], float) for i in range(4)}
    n = len(seq)
    tests = []

    # T1 峰间距恒定 —— 真实数据中峰间距随读长单调漂移
    tests.append((sp.std() < 0.5,
                  f"峰间距标准差 = {sp.std():.3f} scan (真实数据应 >1，随读长漂移)"))

    # T2 分析道长度 == 碱基数 x 固定间距 —— 合成器的痕迹
    L = len(ana[fwo[0]])
    tests.append((L == n * int(round(sp.mean())) and sp.std() < 0.5,
                  f"分析道长度 {L} == 碱基数 {n} × 固定间距 {int(round(sp.mean()))}"))

    # T3 峰高恒定 —— 真实数据峰高跨读长衰减数倍
    h = np.array([ana[b][p] for p, b in zip(loc, seq)])
    cv = h.std() / h.mean()
    ratio = h[:50].mean() / max(h[-50:].mean(), 1e-9)
    tests.append((cv < 0.05,
                  f"主峰高变异系数 CV = {cv:.4f}（真实 0.3~0.8）；"
                  f"首50峰/末50峰 = {ratio:.2f}（真实应显著>1）；唯一峰高值 {len(np.unique(h))} 个"))

    # T4 峰形对称 —— 真实毛细管峰有拖尾, 峰两侧同距点强度不等
    sk = []
    for p, b in zip(loc[10:-10], seq[10:-10]):
        y = ana[b]
        for k in (1, 2, 3, 4):                 # 只比同一通道、峰心等距的两点
            lo, hi = y[p - k], y[p + k]
            if lo + hi > 0:
                sk.append(abs(lo - hi) / (lo + hi))
    sk = float(np.mean(sk))
    tests.append((sk < 0.01,
                  f"峰心等距两点强度差 = {sk:.6f}（真实拖尾峰应 >0.03，完美高斯=0）"))

    # T5 基线 —— 真实数据经基线扣除后必有负值与小幅涨落, 不会硬性归零
    allv = np.concatenate([ana[b] for b in fwo])
    zero, neg = (allv == 0).mean(), (allv < 0).mean()
    tests.append((neg == 0.0 and zero > 0.3,
                  f"分析道取值：恰好为 0 的点占 {zero * 100:.1f}%，负值点占 {neg * 100:.2f}%"
                  f"（真实 KB 分析道基线扣除后必然出现负值）"))

    # T6 QV 恒定 —— 真实读段起始/末端 QV 必然下降
    tests.append((len(np.unique(qv)) == 1,
                  f"QV 唯一取值 {len(np.unique(qv))} 个（真实读段两端 QV 必然衰减）"))

    # T7 内部数组长度自洽性
    lens = {t: len(dec(ab[t]) if isinstance(ab[t], (bytes, str)) else ab[t])
            for t in ("PBAS1", "PCON1", "PLOC1", "P1AM1", "P2AM1", "P2BA1", "P1WD1")
            if t in ab}
    tests.append((len(set(lens.values())) > 1,
                  f"每碱基数组长度不自洽: {lens}"))

    print()
    print("=" * 72)
    print("【2】真实性取证")
    print("=" * 72)
    bad = 0
    for i, (flag, msg) in enumerate(tests, 1):
        mark = "✗ 异常" if flag else "✓ 正常"
        bad += flag
        print(f"  T{i} [{mark}] {msg}")
    print()
    print(f"  异常项 {bad}/{len(tests)}")
    print("  结论: " + ("色谱图为合成/重绘，非真实电泳信号 —— 该文件不能作为测序凭证。"
                        if bad >= 4 else "未检出合成痕迹。"))
    return bad


if __name__ == "__main__":
    path = sys.argv[1]
    ab, dec, fwo = load(path)
    print(f"文件: {path}\n")
    seq, qv, loc, sp = qc(ab, dec, fwo)
    forensics(ab, dec, fwo, seq, qv, loc, sp)
