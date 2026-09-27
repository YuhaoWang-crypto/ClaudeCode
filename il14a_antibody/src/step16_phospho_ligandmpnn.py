"""Step 16 - The phospho arm: put a real phosphate on S515 and redesign the
paratope around it with LigandMPNN.

WHY LIGANDMPNN AND NOT PROTEINMPNN
----------------------------------
ProteinMPNN sees N, CA, C and O and nothing else, so a phosphate is invisible
to it - asking it to "design for the phosphorylated peptide" would return the
same sequences as for the unphosphorylated one and the difference would be
cosmetic. LigandMPNN was built for exactly this: prody selects everything that
is "not protein and not water" as context atoms, so a phosphate group written
as a HETATM residue is seen, and the designed residues are conditioned on it.

WHERE THE PHOSPHATE GOES, AND WHY IT IS NOT ARBITRARY
-----------------------------------------------------
Step 14 chose TXLNA 514-522 for the SM3 scaffold on backbone grounds alone -
it is the window whose prolines align with both of MUC1's. S515 lands at
window position 2 by consequence, not by design, and that position carries 5
heavy-atom contacts with the paratope in the template. So the one annotated
phosphosite in the C-terminal region falls on a contacted position of the
best-threading window. That is worth testing rather than asserting.

THE CONTROL THAT DECIDES WHETHER ANY OF THIS IS REAL
-----------------------------------------------------
Every design is run twice: once with --ligand_mpnn_use_atom_context 1 and once
with 0, same seed, same positions. If the phosphate is genuinely shaping the
paratope, the two sequence sets must differ, and they must differ in the
direction chemistry predicts - more basic residues facing the phosphate. If
they come out the same, "phospho-aware design" was decoration and this script
says so.
"""
import json
import math
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from step13_template_mining import SOLVENT, get_structure
from step14_graft_and_design import THREE2ONE, ONE2THREE, chain_residues, read_atoms_named
from step15_mpnn_rounds import paratope, expand, liabilities

LMPNN = Path("/tmp/claude-0/tools/LigandMPNN")
WORK = Path("/tmp/claude-0/lmpnn_work")

TEMPLATE = {"pdb": "1SM3", "ab": ["H", "L"], "pep": "P"}

# WHICH WINDOW THE PHOSPHO ARM USES IS NOT THE ONE THE PAN ARM USES.
# Step 14 picked TXLNA 514-522 for the unmodified peptide because it aligns
# both of MUC1's prolines. But that puts S515 at a position where all 36
# phosphate rotamers are clash-free - which sounds encouraging and means the
# opposite: the phosphate points into solvent and the paratope never touches
# it. A phospho-specific antibody needs the phosphate BURIED, so the phospho
# arm is chosen on burial and pocket chemistry instead. Both windows below are
# backbone-admissible by step 14's proline test.
WINDOWS = [
    {"id": "pS515_512_520", "window": "APSSPRVTE", "start": 512, "idx": 4,
     "why": "the template residue at this position is Asp, so the pocket "
            "already binds an anion, and it carries H52 Arg"},
    {"id": "pS515_513_521", "window": "PSSPRVTEA", "start": 513, "idx": 3,
     "why": "deepest burial of the phosphate (173 antibody heavy atoms "
            "within 6 A), in the CDR-L3 aromatic box"},
]

# Ideal geometry. Ser chi1 rotamers are the three staggered wells; chi2 (the
# CB-OG-P torsion) is essentially free in solution, so it is scanned.
D_CB_OG, A_CA_CB_OG = 1.417, 110.8
D_OG_P, A_CB_OG_P = 1.60, 120.5
D_P_O, A_OG_P_O = 1.50, 110.0
CHI1_WELLS = (-177.0, -65.0, 62.0)

HARD_CLASH = 3.0     # against main chain - no mutation can relieve it
SOFT_CLASH = 3.2     # against a side chain - a mutation can


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def vdot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vnorm(a):
    n = math.sqrt(vdot(a, a))
    return (a[0] / n, a[1] / n, a[2] / n)


def place(a, b, c, bond, angle, dihedral):
    """NeRF: position D given A-B-C and internal coordinates (degrees)."""
    ang, dih = math.radians(angle), math.radians(dihedral)
    bc = vnorm(vsub(c, b))
    n = vnorm(vcross(vsub(b, a), bc))
    m = vcross(n, bc)
    d2 = (-bond * math.cos(ang),
          bond * math.cos(dih) * math.sin(ang),
          bond * math.sin(dih) * math.sin(ang))
    return tuple(c[i] + d2[0] * bc[i] + d2[1] * m[i] + d2[2] * n[i]
                 for i in range(3))


def build_phosphate(res, chi1, chi2):
    """OG + PO3 on a residue that already has N, CA, CB."""
    N, CA, CB = res["at"]["N"], res["at"]["CA"], res["at"]["CB"]
    OG = place(N, CA, CB, D_CB_OG, A_CA_CB_OG, chi1)
    P = place(CA, CB, OG, D_OG_P, A_CB_OG_P, chi2)
    oxy = [place(CB, OG, P, D_P_O, A_OG_P_O, d) for d in (60.0, 180.0, 300.0)]
    return {"OG": OG}, {"P": P, "O1P": oxy[0], "O2P": oxy[1], "O3P": oxy[2]}


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def clash_profile(group, ab_atoms):
    hard = soft = 0
    near = Counter()
    for _, xyz in group.items():
        for a in ab_atoms:
            d = dist(xyz, a["xyz"])
            if d <= 4.5:
                near[(a["chain"], a["resseq"], a["comp"])] += 1
            if a["name"] in ("N", "CA", "C", "O"):
                if d < HARD_CLASH:
                    hard += 1
            elif d < SOFT_CLASH:
                soft += 1
    return hard, soft, near


def write_pdb(path, prot_chains, het):
    n = 0
    with open(path, "w") as fh:
        for cid, res in prot_chains:
            for i, r in enumerate(res, 1):
                for nm, (x, y, z) in r["at"].items():
                    n += 1
                    atom = nm if len(nm) == 4 else f" {nm:<3s}"
                    fh.write(f"ATOM  {n:5d} {atom} {r['comp']:>3s} {cid}"
                             f"{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                             f"  1.00  0.00          {nm[0]:>2s}\n")
            fh.write("TER\n")
        for (cid, rnum, rname, atoms) in het:
            for nm, (x, y, z) in atoms.items():
                n += 1
                atom = nm if len(nm) == 4 else f" {nm:<3s}"
                fh.write(f"HETATM{n:5d} {atom} {rname:>3s} {cid}"
                         f"{rnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                         f"  1.00  0.00          {nm[0]:>2s}\n")
        fh.write("END\n")
    return n


def run(cmd, cwd=None):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True,
                       text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:])
        print(r.stderr[-2500:])
        raise SystemExit("LigandMPNN failed")
    return r.stdout


def read_designs(folder):
    fa = list((Path(folder) / "seqs").glob("*.fa"))
    out = []
    for f in fa:
        hdr = None
        for line in f.read_text().splitlines():
            if line.startswith(">"):
                hdr = line
            elif hdr:
                kv = dict(re.findall(r"([\w_]+)=([-\d.]+)", hdr))
                # The first record in a LigandMPNN fasta is the INPUT
                # sequence, not a design. Only designs carry id=. Reading it
                # as a design puts SM3's own CDRs at the top of the list.
                if "id" not in kv:
                    hdr = None
                    continue
                if "T" in hdr and "overall_confidence" in kv or "id" in kv:
                    pass
                out.append({
                    "header": hdr, "seqs": line.strip().split(":"),
                    "score": float(kv.get("overall_confidence", "nan")),
                    # reported separately by LigandMPNN and the more relevant
                    # number here: confidence restricted to residues near the
                    # context atoms, i.e. the ones facing the phosphate
                    "ligand_confidence": float(
                        kv.get("ligand_confidence", "nan")),
                    "seq_recovery": float(kv.get("seq_rec", "nan"))})
                hdr = None
    return out


def burial(group, ab_atoms, cutoff=6.0):
    return sum(1 for _, x in group.items() for a in ab_atoms
               if dist(x, a["xyz"]) <= cutoff)


def scan_rotamers(site, ab_atoms):
    rot = []
    for chi1 in CHI1_WELLS:
        for chi2 in range(-180, 180, 30):
            og, po3 = build_phosphate(site, chi1, float(chi2))
            grp = dict(og)
            grp.update(po3)
            hard, soft, near = clash_profile(grp, ab_atoms)
            rot.append({"chi1": chi1, "chi2": float(chi2), "hard": hard,
                        "soft": soft, "burial": burial(po3, ab_atoms),
                        "og": og, "po3": po3,
                        "near": sorted({f"{c}{r}{n}" for (c, r, n) in near})})
    # a main-chain clash is disqualifying; among the rest, prefer the rotamer
    # that buries the phosphate most, and break ties on fewest soft clashes.
    rot.sort(key=lambda r: (r["hard"], -r["burial"], r["soft"]))
    return rot


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    atoms = read_atoms_named(get_structure(TEMPLATE["pdb"]))
    pep0 = [r for r in chain_residues(atoms, TEMPLATE["pep"])
            if r["comp"] not in SOLVENT and {"N", "CA", "C"} <= set(r["at"])]
    ab = {c: [r for r in chain_residues(atoms, c)
              if r["comp"] not in SOLVENT and r["comp"] in THREE2ONE
              and {"N", "CA", "C"} <= set(r["at"])] for c in TEMPLATE["ab"]}
    ab_atoms = [a for a in atoms if a["chain"] in TEMPLATE["ab"]
                and a["comp"] not in SOLVENT and a["elem"] != "H"]
    sel = {c: expand(paratope(ab[c], pep0), len(ab[c])) for c in TEMPLATE["ab"]}
    redes = " ".join(f"{c}{p}" for c in TEMPLATE["ab"] for p in sel[c])

    out_all = []
    for W in WINDOWS:
        assert len(W["window"]) == len(pep0)
        site = pep0[W["idx"] - 1]
        print("=" * 100)
        print(f"{W['id']}   TXLNA {W['start']}-{W['start']+len(pep0)-1} "
              f"{W['window']}   pS515 at window position {W['idx']} "
              f"(template {site['comp']})")
        print(f"  chosen because: {W['why']}")
        print("=" * 100)
        if "CB" not in site["at"]:
            print("  template residue has no CB - cannot build a side chain")
            continue
        rot = scan_rotamers(site, ab_atoms)
        feas = [r for r in rot if r["hard"] == 0]
        print(f"  {len(feas)}/{len(rot)} rotamers free of main-chain clash "
              f"(a main-chain clash cannot be relieved by any mutation)")
        best = feas[0]
        print(f"  best rotamer chi1={best['chi1']:.0f} chi2={best['chi2']:.0f}"
              f"  burial={best['burial']} atoms within 6 A"
              f"  soft clashes={best['soft']} (relievable by design)")
        print(f"  pocket: {', '.join(best['near'])}")
        basic = [x for x in best["near"] if x.endswith(("ARG", "LYS", "HIS"))]
        print(f"  basic residues already lining it: "
              f"{', '.join(basic) if basic else 'NONE'}")

        pep_g = [{"comp": ONE2THREE[aa], "resseq": r["resseq"],
                  "at": {k: v for k, v in r["at"].items()
                         if k in ("N", "CA", "C", "O", "CB")}}
                 for r, aa in zip(pep0, W["window"])]
        chains = [(c, ab[c]) for c in TEMPLATE["ab"]] + [(TEMPLATE["pep"], pep_g)]
        pdb_p = WORK / f"{W['id']}_withPO4.pdb"
        pdb_0 = WORK / f"{W['id']}_noPO4.pdb"
        write_pdb(pdb_p, chains, [("P", 900, "PO4", dict(best["po3"]))])
        write_pdb(pdb_0, chains, [])

        runs = {}
        for tag, pdb, ctx in (("pS_context_ON", pdb_p, 1),
                              ("pS_context_OFF", pdb_p, 0),
                              ("no_phosphate", pdb_0, 1)):
            o = WORK / f"{W['id']}_{tag}"
            shutil.rmtree(o, ignore_errors=True)
            run([sys.executable, LMPNN / "run.py",
                 "--model_type", "ligand_mpnn",
                 "--checkpoint_ligand_mpnn",
                 LMPNN / "model_params" / "ligandmpnn_v_32_010_25.pt",
                 "--pdb_path", pdb, "--out_folder", o,
                 "--redesigned_residues", redes,
                 "--ligand_mpnn_use_atom_context", str(ctx),
                 # an unpaired cysteine in a CDR is a manufacturability
                 # problem, and the smoke test produced one (CDR-H3
                 # TFGGTMCYF), so the alphabet excludes it up front rather
                 # than filtering it out afterwards
                 "--omit_AA", "C",
                 "--batch_size", "10", "--number_of_batches", "10",
                 "--temperature", "0.2", "--seed", "37"], cwd=LMPNN)
            runs[tag] = read_designs(o)
            print(f"  {tag:16} {len(runs[tag])} designs")

        def profile(recs):
            cnt = {(c, p): Counter() for c in TEMPLATE["ab"] for p in sel[c]}
            for r in recs:
                for ci, c in enumerate(TEMPLATE["ab"]):
                    s_ = r["seqs"][ci]
                    for p in sel[c]:
                        cnt[(c, p)][s_[p - 1]] += 1
            return cnt

        on, off = profile(runs["pS_context_ON"]), profile(runs["pS_context_OFF"])
        diffs, tot_on, tot_off = [], 0, 0
        for k in on:
            n_on = sum(on[k].values()) or 1
            n_off = sum(off[k].values()) or 1
            b_on = sum(on[k][x] for x in "RKH") / n_on
            b_off = sum(off[k][x] for x in "RKH") / n_off
            tot_on += sum(on[k][x] for x in "RKH")
            tot_off += sum(off[k][x] for x in "RKH")
            a = on[k].most_common(1)[0][0]
            b = off[k].most_common(1)[0][0]
            if a != b or abs(b_on - b_off) > 0.15:
                diffs.append({"position": f"{k[0]}{k[1]}",
                              "with_phosphate": a, "without": b,
                              "basic_frac_with": round(b_on, 3),
                              "basic_frac_without": round(b_off, 3),
                              "in_pocket": any(f"{k[0]}{k[1]}" == x[:len(f'{k[0]}{k[1]}')]
                                               for x in best["near"])})
        print(f"\n  ATOM-CONTEXT CONTROL (same structure, same seed, "
              f"phosphate visible vs hidden)")
        print(f"    positions whose consensus or basic content moved: "
              f"{len(diffs)}/{len(on)}")
        for d in diffs:
            print(f"      {d['position']:6} {d['without']} -> "
                  f"{d['with_phosphate']}   basic "
                  f"{d['basic_frac_without']:.2f} -> {d['basic_frac_with']:.2f}")
        print(f"    basic residues sampled across the paratope: "
              f"{tot_on} with context vs {tot_off} without "
              f"({tot_on - tot_off:+d})")
        verdict = ("the phosphate is shaping the design"
                   if diffs and tot_on > tot_off else
                   "phosphate changes the consensus but not the basic content"
                   if diffs else
                   "IDENTICAL - the phosphate is NOT influencing the design; "
                   "calling this phospho-aware would be decoration")
        print(f"    -> {verdict}")

        out_all.append({
            "window": W, "rotamers_total": len(rot),
            "rotamers_no_mainchain_clash": len(feas),
            "chosen_rotamer": {k: v for k, v in best.items()
                               if k not in ("og", "po3")},
            "pocket_basic": basic,
            "redesigned_positions": redes.split(),
            "n_designs": {k: len(v) for k, v in runs.items()},
            "context_control": {"positions_changed": diffs,
                                "basic_with_context": tot_on,
                                "basic_without_context": tot_off,
                                "verdict": verdict},
            "designs": {k: [{"header": r["header"], "seqs": r["seqs"]}
                            for r in v[:20]] for k, v in runs.items()},
        })
        print()

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step16_phospho.json").write_text(json.dumps(out_all, indent=2))
    print("wrote results/step16_phospho.json")


if __name__ == "__main__":
    main()
