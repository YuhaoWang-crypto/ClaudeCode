#!/usr/bin/env python3
"""从 .ab1 的原始道 (DATA1-4) 独立反演序列。

当分析道 (DATA9-12) 与判读序列 PBAS1 可疑时，原始道是唯一未被改写的证据。
原始道未做光谱解混和迁移率校正，因此需要：
  1. 滑动分位数扣基线；
  2. 由数据自估 4x4 染料串扰矩阵（对高强度点按主导通道聚类取归一化均值）；
  3. 求逆解混，逐通道寻峰，按 scan 位置合并成序列。

得到的序列约有 10% 错误（无 .mob 迁移率校正），足以做同源比对判定，
不能用作最终序列。

用法: python3 02_raw_trace_basecall.py <file.ab1> [out.fasta]
"""
import sys
import numpy as np
from Bio import SeqIO

BLOB_END = 3100          # 起始游离染料/引物峰区，跳过


def rolling_baseline(y, win=301, q=15, step=25):
    n = len(y)
    b = np.empty(n)
    for i in range(0, n, step):
        a, z = max(0, i - win // 2), min(n, i + win // 2)
        b[i:i + step] = np.percentile(y[a:z], q)
    return b


def estimate_color_matrix(Xb):
    """列 = 染料，行 = 检测通道。对高强度 scan 按主导通道聚类取归一化均值。"""
    tot = Xb.sum(0)
    sel = tot > np.percentile(tot, 80)
    Xn = Xb[:, sel] / (Xb[:, sel].sum(0) + 1e-9)
    dom, pur = Xn.argmax(0), Xn.max(0)
    M = np.zeros((4, 4))
    for k in range(4):
        m = (dom == k) & (pur > np.percentile(pur[dom == k], 75))
        M[:, k] = Xn[:, m].mean(1)
    return M / M.sum(0)


def channel_peaks(y, thr):
    out = []
    for i in range(2, len(y) - 2):
        if y[i] > thr and y[i] >= y[i - 1] and y[i] >= y[i - 2] \
                and y[i] > y[i + 1] and y[i] >= y[i + 2]:
            out.append(i)
    return out


def basecall_raw(path):
    ab = SeqIO.read(path, "abi").annotations["abif_raw"]
    fwo = ab["FWO_1"].decode() if isinstance(ab["FWO_1"], bytes) else ab["FWO_1"]
    X = np.vstack([np.array(ab[f"DATA{1 + i}"], float) for i in range(4)])
    Xb = np.clip(np.vstack([X[i] - rolling_baseline(X[i]) for i in range(4)]), 0, None)

    M = estimate_color_matrix(Xb)
    print("自估染料串扰矩阵 (列=染料, 行=通道 %s):" % fwo)
    print(np.round(M, 3))

    S = np.clip(np.linalg.solve(M, Xb), 0, None)
    S = np.vstack([np.convolve(S[i], np.ones(3) / 3, "same") for i in range(4)])

    peaks = []
    for k, b in enumerate(fwo):
        thr = np.percentile(S[k], 70) * 0.6 + 50
        peaks += [(p, b, S[k][p]) for p in channel_peaks(S[k], thr)]
    peaks.sort()

    merged = []
    for p, b, h in peaks:                       # 同一碱基位置只保留最强通道
        if merged and p - merged[-1][0] < 5:
            if h > merged[-1][2]:
                merged[-1] = (p, b, h)
        else:
            merged.append((p, b, h))
    merged = [m for m in merged if m[0] > BLOB_END]

    seq = "".join(b for _, b, _ in merged)
    pos = [p for p, _, _ in merged]
    print(f"\n反演峰数 = {len(seq)}   跨度 scan {pos[0]}-{pos[-1]}")
    print("峰间距(每60峰):",
          [round(float(np.mean(np.diff(pos[i:i + 60]))), 1) for i in range(0, len(pos) - 60, 60)])
    return seq


if __name__ == "__main__":
    s = basecall_raw(sys.argv[1])
    print("\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)))
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as fh:
            fh.write(">raw_trace_reconstructed\n"
                     + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n")
        print(f"\n已写入 {sys.argv[2]}")
