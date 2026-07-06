# Allosteric-Pocket / Metal-Switch Demo — Target Structures & Candidate-Finding Pipeline

计算 demo 的目标：对一组酶体系，用 **active/inactive（或 apo/ligand-bound）结构对做结构叠合**，
定位活性中心结合状态改变后在 **远端（distal）** 出现的构象响应（Cα 位移、hinge/界面运动），
把 5–20 个远端候选区域排出来，再围绕最高优先级的 hinge/界面 loop 设计
**His-pair / His-triplet 小文库**，筛选 Zn²⁺ / Ni²⁺ / Co²⁺ 对 kcat、Km、kcat/Km 的影响。

This repository packages the **structures/sequences** for the seven candidate systems and a
reproducible **alignment → distal-hotspot → His-site** pipeline.

## Layout

```
data/
  targets.json            # manifest: 7 systems, PDB pairs, roles, references, UniProt IDs
  pdb/                    # 13 downloaded RCSB coordinate files (*.pdb)
  fasta/all_targets.fasta # canonical sequences for every chain (for construct design)
scripts/
  download_structures.py  # fetch all PDBs + FASTA from RCSB (idempotent, retrying)
  find_allosteric_candidates.py  # robust superposition -> distal candidate regions + His-pairs
results/
  <SYS>_candidates.csv    # per-residue Ca displacement & distance-to-active-site
  <SYS>_regions.json      # ranked distal regions + His-pair/triplet geometric suggestions
  summary.json
docs/
  TARGETS.md              # per-system rationale, PDB pairs, wet-lab notes, assays
  EXPERIMENTAL_ROUTE.md   # the alignment -> library -> metal-screen protocol
```

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/download_structures.py            # -> data/pdb/, data/fasta/
python3 scripts/find_allosteric_candidates.py     # -> results/  (GCK, PTP1B, AdK)
```

## Target systems (recommended order)

计算 demo 顺序：**GCK → PTP1B → AdK → PFK → ATCase**
湿实验工程顺序：**AdK 或 PTP1B 优先；PFK 次之；cpTEM-1 作为文献复现型项目**

| System | PDB pair (unbound → bound) | Oligomer | Role |
|--------|---------------------------|----------|------|
| **GCK** human glucokinase | 1V4T → 1V4S | monomer | large domain motion; small-molecule activator site (clean positive control) |
| **PTP1B** | 1SUG → 1T49 (+1PTY) | monomer | known distal allosteric inhibitor site ~20 Å; WPD-loop coupling |
| **AdK** *E. coli* adenylate kinase | 4AKE → 1AKE | monomer | large open/closed; clear hinge → **de novo Zn-switch** target |
| **PFK** *E. coli* phosphofructokinase | 2PFK → 1PFK | tetramer | natural effector (ADP/Mg²⁺); oligomeric-allostery benchmark |
| **TEM-1** β-lactamase | 1BTL | monomer | scaffold for cpTEM-1 metal-switch literature replication |
| **ATCase** *E. coli* | 6AT1 (T) → 1D09 (R, PALA) | 12-mer | classic strong T/R allostery benchmark |
| **GP** human liver glycogen phosphorylase | 3CEH | dimer | AMP-site / allosteric-inhibitor drug-target benchmark |

See `docs/TARGETS.md` for per-system detail and `docs/EXPERIMENTAL_ROUTE.md` for the full protocol.

## Method summary (`find_allosteric_candidates.py`)

1. Match chain-A Cα atoms between the two states by residue number.
2. **Robust core superposition** — iterative outlier rejection locks the reference frame onto
   the rigid domain, so hinge/mobile regions stand out instead of being averaged away.
3. Per-residue Cα displacement + distance to the catalytic site (ligand centroid or catalytic residue).
4. Residues that are **high-displacement AND distal** → clustered into contiguous regions, ranked.
5. Geometric scan flags residue pairs/triples with Cα–Cα spacing (4.5–12 Å) compatible with an
   engineered His-pair / His-triplet metal site.

### Current results (monomeric positive controls)

| System | core-RMSD | overall-RMSD | max Cα disp | top distal region(s) |
|--------|-----------|--------------|-------------|----------------------|
| GCK   | 1.05 Å | 10.89 Å | 39.5 Å | 85–147 (small domain); 65–77 & 446–461 bracket the **allosteric activator site** |
| PTP1B | 0.18 Å | 0.97 Å | 6.1 Å | 185–186 (WPD edge), 239–241, 279–282 (subtle — read loop dynamics, not just max disp) |
| AdK   | 2.02 Å | 8.01 Å | 24.1 Å | **117–163 (LID)**, **31–70 (NMP)** — hinge/counterweight Zn-switch targets |

The GCK run recovers the known small-molecule activator pocket among the top-ranked distal regions —
validating the core hypothesis that active-site state change produces a detectable distal response.
