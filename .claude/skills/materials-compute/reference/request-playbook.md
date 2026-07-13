# Request playbook — client chat → demo + HPC scaffold

Worked mappings from real client requests to the chosen local proxy and the
production scaffold. Use as pattern-matching precedent for new asks.

---

## 1. Quantum-dot polarization / fine-structure splitting
**Ask:** "can you do this theoretical calculation (Fig 2/3)?" — exciton FSS,
polarization angle, uniaxial stress (Gong et al., PRL 106, 227401).
**Local proxy:** `qd_fss_model.py` — the analytical two-level bright-exciton
Hamiltonian, pure numpy. Reproduces the paper's limits exactly
(β=0 ⇒ FSS_min=2|κ|; α=0 ⇒ FSS_min=2|δ|).
**Why local works:** the model IS analytical; the paper only uses atomistic
pseudopotentials to *confirm* it. Full atomistic version → HPC.

## 2. Cr-doped single-atom biomimetic catalyst (peroxidase mimic)
**Ask:** screen metal centers by structural stability / electronic structure /
charge (initial screen), then barriers/selectivity (fine screen). Real system
Cr-Co₃O₄; reference used ~10 metal centers on an M-MnOC scaffold.
**Local proxy:** `metal_center_screen.py` — [M–OH] active-site model, PySCF UKS
DFT, scans spin states, reports gap / metal charge / spin population across
Cr/Mn/Fe/Co/Ni/Cu. (Cr came out with strong Lewis acidity + redox-active oxo +
decent gap — supports the client's Cr hypothesis, qualitatively.)
**HPC production:** periodic Cr-Co₃O₄ slab in VASP/QE, geometry opt, real
coordination, then NEB reaction barriers for the peroxidase cycle.

## 3. Scintillator / self-trapped-exciton (STE) paper — reviewer rebuttal
**Ask (reviewers):** STE mechanism under-supported; theory too shallow; only
ground-state; no phonon or electron-phonon calculations for the "strong e-ph
coupling" claim.
**Local proxy:** `demo3_method_checklist.md` — maps each reviewer point to
concrete calc + method + software + cost.
**HPC production:** `hpc_templates/phonon_qe/` — DFPT phonons (`ph.in`) and the
high-ROI **Huang-Rhys factor + configuration-coordinate diagram**
(`huang_rhys.md`) as the first, cheapest, most persuasive补算; EPW / GW-BSE if
reviewers push further. Priority order: phonon DOS → S factor → ΔSCF STE
geometry → (optional) GW-BSE.

## 4. Molten-salt machine-learning potential
**Ask (郭硕):** 50-100 atom molten salt, AIMD at 6 temperatures, DFT, want a
DeePMD training set and ultimately density/viscosity.
**Local proxy:** `moltensalt_mlp_pipeline.py` (LJ liquid MD → g(r), diffusion
Arrhenius, and **DeePMD/npy packaging** — the client's deliverable format,
validated) + `moltensalt_builder.py` (real FLiNaK / LiCl-KCl / NaCl initial
structures at target density → POSCAR + xyz).
**HPC production:** `hpc_templates/moltensalt_cp2k|vasp/` (AIMD) →
`hpc_templates/deepmd/` (train.json + `workflow.md`: dpdata → dp train → LAMMPS
NPT/Green-Kubo for density & viscosity). Consider DP-GEN active learning.
**Honesty flag:** LJ forces are NOT DFT — the proxy proves plumbing + format;
real physics needs the AIMD.

## 5. Protein Cd²⁺ binding sites → mutants for MST
**Ask (bizlikery):** find which sites bind divalent Cd²⁺, design point mutants
for MST; had used Discovery Studio, asks about AutoDock.
**Local proxy:** `cd_binding_site_predictor.py` — HSAB-weighted donor atoms
(Cys-S > His-N > carboxylate-O) + maximal-clique clustering; ranks sites +
mutation targets. Validated: top site on carbonic anhydrase (1CA2) = the real
His94/His96/His119 metal site.
**HPC production:** `hpc_templates/cd_docking/` — metal-aware docking (AutoDock4
Cd type / GOLD, not vanilla Vina) then QM/MM coordination-sphere refinement
(`qmmm_cluster.md`). Feed the client their protein's PDB to the predictor first.

## 6. Chiral SERS — L vs D amino acid on gold
**Ask (旺仔秋秋唐):** L- vs D-arginine adsorption on Au nanoparticles; DFT
geometry, adsorption energy, vibrational/Raman spectra; SERS difference.
**Local proxy:** `sers_chirality_dft.py` — L-alanine proxy, PySCF DFT; optimizes
L, mirror-constructs D, computes frequencies; proves E(L)=E(D) to ~1e-12 Ha and
spectra identical to ~0.002 cm⁻¹ ⇒ the SERS difference cannot be a free-molecule
effect, it must come from adsorption geometry on the surface.
**HPC production:** `hpc_templates/sers_au_slab/` — Au(111) slab, D3 dispersion,
relax both enantiomers' best poses, finite-difference frequencies + Raman; the
L/D difference in adsorbed spectra is the SERS observable.

---

## Environment reality check
This container has Python 3.11 + pip only. Installable & used: numpy, scipy,
matplotlib, pyscf (molecular DFT), ase, dpdata, biopython, rdkit, pyberny.
NOT present: VASP, Quantum ESPRESSO, CP2K, LAMMPS, EPW, phonopy, ORCA, AutoDock,
and any HPC scheduler. The attached MCP servers are bio/pharma (ChEMBL, PubMed,
bioRxiv, Boltz, ClinicalTrials) — none are materials engines.
