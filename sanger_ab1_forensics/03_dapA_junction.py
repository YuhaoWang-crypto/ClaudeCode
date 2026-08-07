#!/usr/bin/env python3
"""把送检读段定位到 E. coli 染色体，解析 dapA 缺失断点。

野生型 E. coli 该位点基因顺序为  purC - bamC(nlpB/dapX) - dapA - gcvR。
若菌株敲除了 dapA，跨断点的读段会出现：左半段比对到 gcvR 侧、右半段比对到
bamC 侧，而两者在野生型基因组上相隔 ~880 bp。本脚本把两臂分别做局部比对，
求出各自在读段与基因组上的精确边界，从而给出缺失区间与连接处微同源。

需要联网 (NCBI E-utilities)。
用法: python3 03_dapA_junction.py <read.fasta>
"""
import sys
import urllib.parse
import urllib.request

from Bio import Align
from Bio.Seq import Seq

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REF = "NC_000913.3"          # E. coli K-12 MG1655, 野生型参照
WIN_START, WIN_STOP = 2597000, 2601000
# MG1655 该区段注释 (均在负链)
FEATURES = {"bamC": (2597831, 2598865), "dapA": (2598882, 2599760),
            "gcvR": (2599906, 2600478), "purC": (2596905, 2597618)}


def efetch_fasta(acc, start, stop):
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "fasta",
                                "retmode": "text", "seq_start": start, "seq_stop": stop})
    txt = urllib.request.urlopen(f"{EUTILS}?{q}", timeout=180).read().decode()
    return "".join(txt.split("\n")[1:]).upper()


def read_fasta(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def aligner():
    a = Align.PairwiseAligner()
    a.mode, a.match_score, a.mismatch_score = "local", 2, -3
    a.open_gap_score, a.extend_gap_score = -10, -4
    return a


def map_arm(al, wt, off, sub, sub_off, lo, hi, label):
    aln = al.align(wt[lo - off:hi - off], sub)[0]
    t, q = aln.aligned
    rs, re_ = q[0][0] + sub_off + 1, q[-1][1] + sub_off
    ws, we = t[0][0] + lo, t[-1][1] + lo - 1
    print(f"  {label}: 读段 {rs}-{re_}  <->  {REF} {ws:,}-{we:,}   score={aln.score:.0f}")
    return rs, re_, ws, we


def main(path):
    read = read_fasta(path)
    wt = efetch_fasta(REF, WIN_START, WIN_STOP)
    al = aligner()

    # 读段方向未知，取与参考正链同向的一条
    fwd, rev = read, str(Seq(read).reverse_complement())
    r = fwd if al.score(wt, fwd) >= al.score(wt, rev) else rev
    print(f"读段长 {len(read)} bp；与 {REF} 正链同向的是 "
          f"{'原始方向' if r is fwd else '反向互补'}\n")

    print("两臂精确边界:")
    _, r_end, _, w_end = map_arm(al, wt, WIN_START, r[:420], 0, 2598300, 2599100, "右臂 (bamC 侧)")
    r_sta, _, w_sta, _ = map_arm(al, wt, WIN_START, r[300:], 300, 2599600, 2600400, "左臂 (gcvR 侧)")

    ov = r_end - r_sta + 1
    print()
    if ov > 0:
        print(f"  两臂在读段上重叠 {ov} bp -> 连接处微同源: {r[r_sta - 1:r_end]}")
    else:
        print(f"  两臂之间插入 {-ov} bp: {r[r_end:r_sta - 1]}")

    dele = w_sta - 1 - w_end
    print(f"\n【缺失区间】{REF} {w_end + 1:,} .. {w_sta - 1:,}  = {dele} bp")
    for g, (s, e) in sorted(FEATURES.items(), key=lambda x: x[1][0]):
        inside = max(0, min(e, w_sta - 1) - max(s, w_end + 1) + 1)
        print(f"    {g:5s} {s:,}-{e:,} ({e - s + 1} bp)  落入缺失区 {inside} bp "
              f"({inside / (e - s + 1) * 100:.1f}%)")

    print(f"\n断点邻近序列 (小写=bamC 侧 / | / 大写=gcvR 侧):")
    print(f"  ...{r[max(0, r_end - 45):r_end].lower()} | {r[r_end:r_end + 45]}...")


if __name__ == "__main__":
    main(sys.argv[1])
