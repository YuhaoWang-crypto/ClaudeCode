#!/usr/bin/env python3
"""
分析 Boltz-2 预测的 KRAS–配体复合物：配体落在哪个口袋？与第 12 位残基接触吗？

用法: python3 analyze_pocket.py structures/g12d_mrtx1133.cif
无第三方依赖，手工解析 mmCIF 的 atom_site 环。
距离单位 Å。
"""
import sys, math

# KRAS 功能区（1-based 残基编号，人源 KRAS4B 1-169）
REGIONS = {
    "P-loop(10-17)":      set(range(10, 18)),   # 核苷酸磷酸结合环（内源 GTP/GDP 位点一部分）
    "Switch-I(30-38)":    set(range(30, 39)),
    "Switch-II(60-76)":   set(range(60, 77)),   # S-IIP 隐匿口袋主体
    "a3-helix(87-104)":   set(range(87, 105)),  # 含 His95 次级沟槽
    "NKxD(116-119)":      set(range(116, 120)), # 鸟嘌呤识别基序（核苷酸口袋）
}
MUT_POS = 12  # 突变热点

def parse_cif(path):
    atoms = []
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("_atom_site."):
            # 收集列名
            cols = []
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip().split(".", 1)[1])
                i += 1
            idx = {c: k for k, c in enumerate(cols)}
            # 读取数据行直到下一个循环/块
            while i < len(lines):
                row = lines[i]
                s = row.strip()
                if s == "" or s.startswith("#") or s.startswith("loop_") or s.startswith("_") or s.startswith("data_"):
                    break
                parts = s.split()
                if len(parts) < len(cols):
                    i += 1
                    continue
                def g(name, default="?"):
                    return parts[idx[name]] if name in idx else default
                try:
                    atoms.append({
                        "group": g("group_PDB"),
                        "atom":  g("label_atom_id").strip('"'),
                        "comp":  g("label_comp_id"),
                        "chain": g("label_asym_id"),
                        "seq":   g("label_seq_id"),
                        "x": float(g("Cartn_x")), "y": float(g("Cartn_y")), "z": float(g("Cartn_z")),
                    })
                except (ValueError, KeyError):
                    pass
                i += 1
            continue
        i += 1
    return atoms

def dist(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def main(path):
    atoms = parse_cif(path)
    # 配体 = HETATM（非聚合物），或链不是蛋白主链；Boltz 把蛋白放 ATOM、配体放 HETATM
    lig = [a for a in atoms if a["group"] == "HETATM"]
    prot = [a for a in atoms if a["group"] == "ATOM"]
    if not lig:
        # 退化处理：按链 L 识别
        lig = [a for a in atoms if a["chain"] == "L"]
        prot = [a for a in atoms if a["chain"] != "L"]
    # 剔除配体里的金属/水等（保留有坐标的全部重原子即可，这里不特判）
    print(f"\n=== {path} ===")
    print(f"蛋白原子 {len(prot)}，配体原子 {len(lig)}")

    # 每个蛋白残基到配体的最近距离
    res_min = {}   # (seq)->(min_dist, comp)
    for a in prot:
        try:
            seq = int(a["seq"])
        except ValueError:
            continue
        dmin = min(dist(a, l) for l in lig)
        if seq not in res_min or dmin < res_min[seq][0]:
            res_min[seq] = (dmin, a["comp"])

    CUT = 4.5
    contacts = sorted([(s, d, c) for s, (d, c) in res_min.items() if d <= CUT], key=lambda t: t[1])
    print(f"\n口袋接触残基（最近重原子距离 ≤ {CUT} Å）：")
    region_hits = {k: [] for k in REGIONS}
    for s, d, c in contacts:
        tags = [name for name, members in REGIONS.items() if s in members]
        tagstr = (" ["+", ".join(tags)+"]") if tags else ""
        star = " ★突变位点" if s == MUT_POS else ""
        print(f"  {c}{s:<4} {d:4.2f} Å{tagstr}{star}")
        for t in tags:
            region_hits[t].append(s)

    print("\n命中功能区统计：")
    for name, members in REGIONS.items():
        hit = sorted(set(region_hits[name]))
        print(f"  {name:18s}: {len(hit)} 残基接触 {hit if hit else ''}")

    # 第 12 位残基的最近接触（突变位点）
    r12 = [a for a in prot if a["seq"] == str(MUT_POS)]
    if r12 and lig:
        comp12 = r12[0]["comp"]
        d12 = min(dist(a, l) for a in r12 for l in lig)
        print(f"\n第12位残基 = {comp12}12，到配体最近重原子 = {d12:.2f} Å")
        # 若为 ASP，专门看羧基 OD1/OD2（盐桥供体）
        if comp12 == "ASP":
            od = [a for a in r12 if a["atom"] in ("OD1", "OD2")]
            if od:
                dod = min(dist(a, l) for a in od for l in lig)
                # 找最近的配体氮原子（二胺质子化端）
                lig_n = [l for l in lig if l["atom"].startswith("N")]
                dN = min((dist(a, l) for a in od for l in lig_n), default=float("nan"))
                print(f"  Asp12 羧基(OD1/OD2) ↔ 配体最近 = {dod:.2f} Å；↔ 配体最近 N = {dN:.2f} Å （盐桥判据 <4 Å）")
        if comp12 == "CYS":
            sg = [a for a in r12 if a["atom"] == "SG"]
            if sg:
                dsg = min(dist(a, l) for a in sg for l in lig)
                print(f"  Cys12 硫(SG) ↔ 配体最近重原子 = {dsg:.2f} Å （共价弹头靶点；Boltz 未显式建键）")

if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
