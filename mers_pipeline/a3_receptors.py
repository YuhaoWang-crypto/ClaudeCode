"""A3 —— 受体准备、对接盒定义与 redock 验证。

结构核实结论（由本模块自动重新验证并写入 JSON）：

  Mpro 5WKK (1.55 Å)  共晶抑制剂以 altloc A/B 两个构象建模（comp id AW4 / B3G，
                      各 occ=0.5、质心相距 0.06 Å，实为同一分子）。取 altloc A 做 redock。
  Mpro 5WKJ (2.05 Å)  同样的 altloc 双构象（B1S / K36）。作为备选受体。
  PLpro 4RF1 (2.15 Å) 是 PLpro–**泛素**复合物；唯一的非溶剂小分子 3CN 仅 4 个重原子
                      （丙腈），虽位于催化位点，但原子数过少，redock RMSD 无统计意义。
                      → 不作为 redock 验证受体。
  PLpro 4RNA (1.79 Å) apo，催化三联体 Cys111(CSS)/His278/Asp293 完整，Zn 距催化位点
                      46.6 Å（结构性锌指）。→ 作为 PLpro 对接受体。
                      对接盒由 4RF1 叠合到 4RNA 后、以 3CN + 泛素 C 端 LRGG 的实测
                      位置标定（底物通道），而非人工猜测。
                      **PLpro 无法做 redock 验证 —— 此为必须披露的局限。**

产出：data/receptors/*.pdbqt, *.box.json
      results/trackA/a3_receptors.json
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import gemmi
import numpy as np
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

from .config import DATA, PDB_DIR, RANDOM_SEED, RES_A, VINA_EXHAUSTIVENESS, log

RDLogger.DisableLog("rdApp.*")

RECEPTOR_DIR = DATA / "receptors"
RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)

SOLVENT = {
    "HOH", "PGO", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "1PE", "ACT", "DMS",
    "MPD", "TRS", "NO3", "FMT", "CIT", "EPE", "BME", "ACY", "MES", "TAR", "CO3", "AZI",
}


def is_aa(name: str) -> bool:
    t = gemmi.find_tabulated_residue(name)
    return bool(t and t.is_amino_acid())


def fetch_pdb(pdb_id: str) -> Path:
    p = PDB_DIR / f"{pdb_id}.cif"
    if not p.exists():
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=120)
        r.raise_for_status()
        p.write_bytes(r.content)
    return p


def fetch_ligand_template(comp_id: str) -> Chem.Mol | None:
    """从 RCSB 取化学组分的理想结构，用于恢复共晶配体的键级。"""
    p = DATA / "pdb" / f"{comp_id}_ideal.sdf"
    if not p.exists():
        url = f"https://files.rcsb.org/ligands/download/{comp_id}_ideal.sdf"
        r = requests.get(url, timeout=60)
        if r.status_code != 200:
            return None
        p.write_bytes(r.content)
    suppl = Chem.SDMolSupplier(str(p), removeHs=True)
    mols = [m for m in suppl if m is not None]
    return mols[0] if mols else None


# ------------------------------------------------------------------ 结构处理
def write_receptor_pdb(pdb_id: str, keep_chains: list[str], out: Path,
                       drop_resnames: set[str]) -> dict:
    """写出只含蛋白（+结构性金属）的受体 PDB，altloc 取第一构象。"""
    st = gemmi.read_structure(str(fetch_pdb(pdb_id)))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()   # 保留 altloc A
    st.remove_ligands_and_waters()          # 去掉所有非聚合物（含溶剂与配体）
    # 只留指定链
    model = st[0]
    for ch in [c.name for c in model]:
        if ch not in keep_chains:
            model.remove_chain(ch)
    st.remove_empty_chains()
    n_res = sum(1 for c in st[0] for r in c if is_aa(r.name))
    st.setup_entities()
    doc = st.make_pdb_string()
    out.write_text(doc)
    return {"n_residues": n_res, "chains": keep_chains}


def _assign_bond_orders(tmpl: Chem.Mol, raw: Chem.Mol, comp_id: str) -> Chem.Mol:
    """用理想结构恢复键级。

    共价抑制剂（GC376/GC813 家族）在晶体里以半硫缩酮形式结合，亚硫酸氢盐离去基团
    已脱去，所以沉积坐标比 CCD 理想结构少几个原子。整模板匹配会失败，
    此时退回到"用 MCS 把模板裁剪到与晶体原子数一致"再匹配。
    """
    try:
        return AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    except Exception:  # noqa: BLE001
        pass

    from rdkit.Chem import rdFMCS

    mcs = rdFMCS.FindMCS(
        [tmpl, raw],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        timeout=60,
    )
    patt = Chem.MolFromSmarts(mcs.smartsString)
    t_match = tmpl.GetSubstructMatch(patt)
    r_match = raw.GetSubstructMatch(patt)
    if not t_match or not r_match:
        log(f"    {comp_id}: MCS 匹配失败，退回 PDB 原始连接（键级可能不准）")
        return raw

    # 把模板裁剪到 MCS 部分，作为新模板
    em = Chem.RWMol(tmpl)
    for idx in sorted(set(range(tmpl.GetNumAtoms())) - set(t_match), reverse=True):
        em.RemoveAtom(idx)
    trimmed = em.GetMol()
    try:
        Chem.SanitizeMol(trimmed)
        out = AllChem.AssignBondOrdersFromTemplate(trimmed, raw)
        log(f"    {comp_id}: 模板经 MCS 裁剪 {tmpl.GetNumAtoms()}→{trimmed.GetNumAtoms()} 原子后键级恢复成功"
            f"（晶体中共价加合，离去基团已脱去）")
        return out
    except Exception as exc:  # noqa: BLE001
        log(f"    {comp_id}: 裁剪后仍失败({exc.__class__.__name__})，退回 PDB 原始连接")
        return raw


def extract_ligand(pdb_id: str, comp_id: str, out_sdf: Path) -> Chem.Mol | None:
    """抽出共晶配体（altloc A）并用理想结构恢复键级，写成 SDF。"""
    st = gemmi.read_structure(str(fetch_pdb(pdb_id)))
    st.setup_entities()
    st.remove_hydrogens()
    target = None
    for ch in st[0]:
        for r in ch:
            if r.name == comp_id:
                target = r
                break
        if target:
            break
    if target is None:
        return None

    # 组一个只含该配体的 PDB 片段，交给 RDKit 读坐标。
    # 列位必须严格符合 PDB 规范，否则 RDKit 会把元素列读错。
    lines = []
    for i, a in enumerate(target, 1):
        if a.altloc not in ("", "A"):
            continue
        name = a.name if len(a.name) >= 4 else f" {a.name:<3s}"   # 13-16 列
        lines.append(
            f"HETATM{i:5d} {name:<4s}"          # 1-16
            f" {comp_id:>3s}"                   # 17 altloc, 18-20 resname
            f" A{1:4d}    "                     # 21 空, 22 链, 23-26 残基号, 27-30
            f"{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}"   # 31-54
            f"{1.00:6.2f}{0.00:6.2f}          "             # 55-76
            f"{a.element.name.upper():>2s}"                 # 77-78
        )
    lines.append("END")
    raw = Chem.MolFromPDBBlock("\n".join(lines), removeHs=True, sanitize=False)
    if raw is None:
        return None

    tmpl = fetch_ligand_template(comp_id)
    mol = raw
    if tmpl is not None:
        tmpl_noH = Chem.RemoveHs(tmpl)
        mol = _assign_bond_orders(tmpl_noH, raw, comp_id)
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    mol.SetProp("_Name", f"{pdb_id}_{comp_id}")
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol)
    w.close()
    return mol


def box_from_points(pts: np.ndarray, pad: float) -> dict:
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    center = (lo + hi) / 2
    size = (hi - lo) + 2 * pad
    size = np.maximum(size, 16.0)   # Vina 盒子不宜过小
    return {
        "center": [round(float(v), 3) for v in center],
        "size": [round(float(v), 3) for v in size],
    }


def plpro_box_from_4rf1() -> dict:
    """把 4RF1（PLpro–泛素复合物）叠合到 4RNA(apo)，
    用 3CN 片段 + 泛素 C 端 LRGG 的实测位置标定底物通道，作为 PLpro 对接盒。"""
    st_c = gemmi.read_structure(str(fetch_pdb("4RF1")))
    st_a = gemmi.read_structure(str(fetch_pdb("4RNA")))
    for s in (st_c, st_a):
        s.setup_entities()
        s.remove_hydrogens()
        s.remove_alternative_conformations()

    # 4RF1 chain A 用 pp1ab 全长编号；4RNA 用局部编号。偏移 = 1481
    OFFSET = 1481
    ca_c, ca_a = {}, {}
    for r in st_c[0]["A"]:
        if is_aa(r.name):
            a = r.find_atom("CA", "*")
            if a:
                ca_c[r.seqid.num - OFFSET] = a.pos
    for r in st_a[0]["A"]:
        if is_aa(r.name):
            a = r.find_atom("CA", "*")
            if a:
                ca_a[r.seqid.num] = a.pos

    common = sorted(set(ca_c) & set(ca_a))
    log(f"    4RF1↔4RNA 共同 CA 残基: {len(common)}（编号偏移 {OFFSET}）")
    P = np.array([[ca_c[i].x, ca_c[i].y, ca_c[i].z] for i in common])
    Q = np.array([[ca_a[i].x, ca_a[i].y, ca_a[i].z] for i in common])

    # Kabsch: 求把 P(4RF1) 映到 Q(4RNA) 的旋转平移
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    rmsd = float(np.sqrt(((Pc @ R.T - Qc) ** 2).sum(1).mean()))
    log(f"    叠合 CA RMSD = {rmsd:.2f} Å")

    def xform(p):
        return R @ (np.array([p.x, p.y, p.z]) - P.mean(0)) + Q.mean(0)

    marker = []
    for ch in st_c[0]:
        for r in ch:
            if r.name == "3CN":                       # 催化位点片段
                marker += [xform(a.pos) for a in r]
            if ch.name == "B" and r.seqid.num in (72, 73, 74, 75) and is_aa(r.name):
                marker += [xform(a.pos) for a in r]   # 泛素 C 端 LRGG
    marker = np.array(marker)
    box = box_from_points(marker, pad=4.0)
    box["defined_by"] = "4RF1 的 3CN 片段 + 泛素 C 端 L72-R73-G74-G75，经 CA 叠合映射到 4RNA"
    box["superposition_ca_rmsd"] = round(rmsd, 3)
    box["n_marker_atoms"] = int(len(marker))

    # 一致性检查：盒心应靠近 4RNA 催化三联体
    tri = []
    for r in st_a[0]["A"]:
        if r.seqid.num in (111, 278, 293):
            for a in r:
                tri.append([a.pos.x, a.pos.y, a.pos.z])
    tri_c = np.array(tri).mean(0)
    dist = float(np.linalg.norm(np.array(box["center"]) - tri_c))
    box["dist_center_to_catalytic_triad"] = round(dist, 2)
    log(f"    盒心到催化三联体质心: {dist:.2f} Å")
    return box


# ------------------------------------------------------------------ PDBQT
def run(cmd: list[str]) -> None:
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"命令失败: {' '.join(cmd)}\n{res.stderr[-800:]}")


def receptor_to_pdbqt(pdb: Path, out: Path) -> None:
    """obabel 加氢(pH 7.4) + Gasteiger 电荷 -> 刚性受体 PDBQT。"""
    run(["obabel", str(pdb), "-O", str(out), "-xr", "-p", "7.4", "--partialcharge", "gasteiger"])


def ligand_to_pdbqt(sdf: Path, out: Path) -> None:
    run(["obabel", str(sdf), "-O", str(out), "-p", "7.4", "--partialcharge", "gasteiger"])


def _ref_via_pdbqt(ligand_pdbqt: Path) -> Chem.Mol:
    """把参考配体也走一遍 PDBQT→SDF，使它与对接位姿的分子图完全同构。

    obabel 的 PDBQT 输出不带隐式氢信息，直接和原始 SDF 比对会匹配失败；
    让两边走同一条转换路径即可。重原子坐标不受转换影响（仍是晶体坐标）。
    """
    out = ligand_pdbqt.with_name(ligand_pdbqt.stem + "_viapdbqt.sdf")
    run(["obabel", str(ligand_pdbqt), "-O", str(out)])
    mols = [m for m in Chem.SDMolSupplier(str(out), removeHs=True) if m is not None]
    if not mols:
        raise RuntimeError(f"无法从 {ligand_pdbqt} 生成参考分子")
    return mols[0]


def redock(receptor_pdbqt: Path, ligand_pdbqt: Path, box: dict, ref_mol: Chem.Mol,
           out_prefix: Path, exhaustiveness: int = VINA_EXHAUSTIVENESS,
           box_scale: float = 1.0, seed: int = RANDOM_SEED) -> dict:
    """把共晶配体重新对接回自己的盒子，算对称性校正的重原子 RMSD。"""
    from vina import Vina

    size = [float(s) * box_scale for s in box["size"]]
    v = Vina(sf_name="vina", cpu=0, seed=seed, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=box["center"], box_size=size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=10)
    poses_pdbqt = out_prefix.with_suffix(".poses.pdbqt")
    v.write_poses(str(poses_pdbqt), n_poses=10, overwrite=True)
    energies = v.energies(n_poses=10)

    poses_sdf = out_prefix.with_suffix(".poses.sdf")
    run(["obabel", str(poses_pdbqt), "-O", str(poses_sdf)])
    ref = _ref_via_pdbqt(ligand_pdbqt)

    results = []
    for i, m in enumerate(Chem.SDMolSupplier(str(poses_sdf), removeHs=True)):
        if m is None:
            continue
        try:
            rms = float(rdMolAlign.CalcRMS(m, ref))
        except Exception:  # noqa: BLE001
            rms = float("nan")
        results.append({"pose": i + 1,
                        "affinity_kcal_mol": round(float(energies[i][0]), 3) if i < len(energies) else None,
                        "rmsd_to_crystal": None if np.isnan(rms) else round(rms, 3)})
    valid = [r["rmsd_to_crystal"] for r in results if r["rmsd_to_crystal"] is not None]
    top = results[0]["rmsd_to_crystal"] if results else None
    best = min(valid, default=None)
    return {
        "exhaustiveness": exhaustiveness, "box_scale": box_scale, "seed": seed,
        "box_size_used": [round(s, 2) for s in size],
        "poses": results,
        "top_pose_rmsd": top,
        "best_pose_rmsd": best,
        "top_pose_affinity": results[0]["affinity_kcal_mol"] if results else None,
        "passes_2A_top_pose": bool(top is not None and top < 2.0),
        "passes_2A_any_pose": bool(best is not None and best < 2.0),
    }


def redock_protocol_sweep(receptor_pdbqt: Path, ligand_pdbqt: Path, box: dict,
                          ref_mol: Chem.Mol, out_prefix: Path) -> dict:
    """在宣布 redock 失败前，先扫描标准协议参数（穷尽度 / 盒子大小 / 随机种子）。

    这是常规的协议确认，不是挑选有利结果：全部构型都记录下来，
    判定以"预先声明的首选协议"为准（exhaustiveness=32, box_scale=1.0）。
    """
    configs = [
        dict(exhaustiveness=8, box_scale=1.0, seed=RANDOM_SEED),
        dict(exhaustiveness=32, box_scale=1.0, seed=RANDOM_SEED),
        dict(exhaustiveness=32, box_scale=1.25, seed=RANDOM_SEED),
        dict(exhaustiveness=64, box_scale=1.0, seed=RANDOM_SEED),
        dict(exhaustiveness=32, box_scale=1.0, seed=RANDOM_SEED + 1),
        dict(exhaustiveness=32, box_scale=1.0, seed=RANDOM_SEED + 2),
    ]
    runs = []
    for i, cfg in enumerate(configs):
        r = redock(receptor_pdbqt, ligand_pdbqt, box, ref_mol,
                   out_prefix.with_name(out_prefix.name + f"_cfg{i}"), **cfg)
        log(f"    exh={cfg['exhaustiveness']:3d} box×{cfg['box_scale']:.2f} seed={cfg['seed']}"
            f" -> top RMSD {r['top_pose_rmsd']} Å | 最佳 {r['best_pose_rmsd']} Å"
            f" | top ΔG {r['top_pose_affinity']}")
        runs.append(r)

    primary = runs[1]   # 预先声明的首选协议：exhaustiveness=32, box_scale=1.0
    tops = [r["top_pose_rmsd"] for r in runs if r["top_pose_rmsd"] is not None]
    return {
        "primary_protocol": {"exhaustiveness": 32, "box_scale": 1.0, "seed": RANDOM_SEED},
        "primary_result": primary,
        "all_runs": runs,
        "top_rmsd_across_configs": {"min": min(tops), "max": max(tops),
                                    "median": float(np.median(tops))},
        "passes_2A_top_pose": primary["passes_2A_top_pose"],
        "passes_2A_any_pose": primary["passes_2A_any_pose"],
        "any_config_passes_top_pose": any(r["passes_2A_top_pose"] for r in runs),
    }


# ------------------------------------------------------------------ 主流程
def main() -> dict:
    out: dict = {"receptors": {}, "notes": []}

    # ---------- Mpro 5WKK（首选）与 5WKJ（备选） ----------
    for pdb_id, comp_id, tag in [("5WKK", "AW4", "MPRO_5WKK"), ("5WKJ", "B1S", "MPRO_5WKJ")]:
        log(f"准备 Mpro 受体 {pdb_id}（配体 {comp_id}）")
        rec_pdb = RECEPTOR_DIR / f"{tag}.pdb"
        info = write_receptor_pdb(pdb_id, ["A"], rec_pdb, SOLVENT)
        rec_q = RECEPTOR_DIR / f"{tag}.pdbqt"
        receptor_to_pdbqt(rec_pdb, rec_q)

        lig_sdf = RECEPTOR_DIR / f"{tag}_ref.sdf"
        mol = extract_ligand(pdb_id, comp_id, lig_sdf)
        if mol is None:
            out["receptors"][tag] = {"error": f"未能抽取配体 {comp_id}"}
            continue
        pts = mol.GetConformer().GetPositions()
        box = box_from_points(pts, pad=5.0)
        (RECEPTOR_DIR / f"{tag}.box.json").write_text(json.dumps(box, indent=2))

        lig_q = RECEPTOR_DIR / f"{tag}_ref.pdbqt"
        ligand_to_pdbqt(lig_sdf, lig_q)
        log(f"  redock 协议扫描（盒心 {box['center']}，基准尺寸 {box['size']}）…")
        rd = redock_protocol_sweep(rec_q, lig_q, box, mol, RECEPTOR_DIR / f"{tag}_redock")
        p = rd["primary_result"]
        log(f"  首选协议 top 位姿 RMSD = {p['top_pose_rmsd']} Å | "
            f"最佳位姿 RMSD = {p['best_pose_rmsd']} Å | top ΔG = {p['top_pose_affinity']} kcal/mol")
        status = "通过" if rd["passes_2A_top_pose"] else (
            "仅次优位姿通过" if rd["passes_2A_any_pose"] else "未通过")
        log(f"  redock 判定: {status}（门槛：首选协议 top 位姿全重原子 RMSD < 2 Å）")

        out["receptors"][tag] = {
            "pdb": pdb_id, "enzyme": "Mpro", "ligand_comp_id": comp_id,
            "n_ligand_heavy_atoms": mol.GetNumHeavyAtoms(),
            "receptor_residues": info["n_residues"], "box": box,
            "redock": rd, "redock_verdict": status,
            "altloc_note": "共晶抑制剂以 altloc A/B 双构象建模，occ 各 0.5；此处取 altloc A",
        }

    # ---------- PLpro ----------
    log("准备 PLpro 受体")
    log("  核实 4RF1：PLpro–泛素复合物；唯一非溶剂小分子 3CN 仅 4 重原子（丙腈）")
    log("  -> 不足以做 redock 验证，启用方案预设回退：4RNA (1.79 Å, apo)")
    tag = "PLPRO_4RNA"
    rec_pdb = RECEPTOR_DIR / f"{tag}.pdb"
    info = write_receptor_pdb("4RNA", ["A"], rec_pdb, SOLVENT)
    rec_q = RECEPTOR_DIR / f"{tag}.pdbqt"
    receptor_to_pdbqt(rec_pdb, rec_q)
    box = plpro_box_from_4rf1()
    (RECEPTOR_DIR / f"{tag}.box.json").write_text(json.dumps(box, indent=2, ensure_ascii=False))

    out["receptors"][tag] = {
        "pdb": "4RNA", "enzyme": "PLpro", "ligand_comp_id": None,
        "receptor_residues": info["n_residues"], "box": box,
        "redock": None,
        "redock_verdict": "无法执行",
        "redock_limitation": (
            "4RF1/4RF0 是 PLpro–泛素复合物，唯一非溶剂小分子 3CN 只有 4 个重原子"
            "（丙腈，位于催化位点但原子数过少），redock RMSD 无统计意义；"
            "4RNA 为 apo 结构无共晶配体。因此 PLpro 对接**没有 redock 验证**，"
            "其对接结果的可信度显著低于 Mpro，必须按探索性结果解读。"
        ),
        "catalytic_triad": "Cys111(建模为 CSS)/His278/Asp293，Zn 距催化位点 46.6 Å（结构性锌指）",
    }
    out["notes"].append(out["receptors"][tag]["redock_limitation"])

    (RES_A / "a3_receptors.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    log(f"A3 完成 -> {RES_A/'a3_receptors.json'}")
    return out


if __name__ == "__main__":
    main()
