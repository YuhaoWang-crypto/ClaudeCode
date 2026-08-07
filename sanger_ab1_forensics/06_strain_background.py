#!/usr/bin/env python3
"""判定读段的菌株背景。

跨缺失断点的读段不能整条比对到任何野生型基因组（断点两侧在野生型上相隔 ~880 bp），
所以"全长最佳命中"会被唯一含同一断点的基因组吃掉，从而误判背景。
正确做法是把两臂拆开，各自对候选菌株做限定物种的 BLAST，比较一致性。

用法: python3 06_strain_background.py <read.fasta> [断点在读段上的位置, 默认 365]
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

from Bio.Seq import Seq

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
STRAINS = [
    ("BL21", "Escherichia coli BL21[Organism]"),
    ("Nissle 1917", "Escherichia coli Nissle 1917[Organism]"),
    ("K-12 MG1655", "Escherichia coli str. K-12 substr. MG1655[Organism]"),
]


def _http(fn, *a, tries=4):
    for i in range(tries):
        try:
            return fn(*a)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def put(seq, entrez):
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "nt", "QUERY": seq,
         "MEGABLAST": "on", "HITLIST_SIZE": "5", "ENTREZ_QUERY": entrez,
         "FORMAT_TYPE": "Text"}
    r = _http(lambda: urllib.request.urlopen(
        URL, urllib.parse.urlencode(p).encode(), timeout=120).read().decode("utf-8", "replace"))
    m = re.search(r"RID = (\S+)", r)
    return m.group(1) if m else None


def fetch(rid):
    p = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text"}
    for _ in range(30):
        t = _http(lambda: urllib.request.urlopen(
            URL + "?" + urllib.parse.urlencode(p), timeout=180).read().decode("utf-8", "replace"))
        if "Status=WAITING" in t:
            time.sleep(12)
            continue
        return t
    return ""


def main(path, cut=365):
    read = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    r = str(Seq(read).reverse_complement())
    arms = {f"bamC 臂 ({cut} bp)": r[:cut], f"gcvR 臂 ({len(r) - cut + 1} bp)": r[cut - 1:]}

    jobs = {}
    for an, aseq in arms.items():
        for sn, eq in STRAINS:
            jobs[(an, sn)] = put(aseq, eq)
            time.sleep(2)

    print(f"{'':22s}" + "".join(f"{sn:>20s}" for sn, _ in STRAINS))
    for an, aseq in arms.items():
        row = f"{an:22s}"
        for sn, _ in STRAINS:
            t = fetch(jobs[(an, sn)])
            m = re.search(r"Identities = (\d+)/(\d+)\s*\((\d+)%\)", t)
            row += f"{(m.group(1) + '/' + m.group(2) + ' ' + m.group(3) + '%') if m else '无命中':>20s}"
        print(row)
    print("\n注：该位点 BL21 与 K-12 序列相同，故本检验只能排除某个背景，不能区分 BL21 与 K-12。")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 365)
