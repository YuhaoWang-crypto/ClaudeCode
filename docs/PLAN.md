# Technical plan — AI-guided lipid & targeting-ligand library screening

This document turns the request into a concrete, staged engineering plan and
records the key scientific reframing that makes the project tractable.

---

## 0. The core reframing: two problems, not one

The original idea describes a single pipeline:

> lipid in the LNP → recognized by a receptor (specificity) → maybe needs to be
> extracted from the particle → target endocytosis → after internalization,
> pH change → release.

This mixes up **two mechanistically separate design problems**. Separating them
is the single most important decision in the project, because they need
different data, different models, and different validation.

### Problem A — ionizable-lipid design (what the two papers do)

The ionizable lipid is one of 4 components of an LNP (ionizable lipid + helper
phospholipid + cholesterol + PEG-lipid). It is **not a receptor ligand** and is
**not extracted** from the particle to be "recognized". Its jobs are:

1. **Encapsulate** mRNA/siRNA (positive charge at low pH during formulation).
2. **Endosomal escape**: after the whole particle is endocytosed, the endosome
   acidifies (pH ~6.5→5); the lipid re-protonates, becomes fusogenic, disrupts
   the endosomal membrane, and releases cargo to the cytosol. This is the
   "pH change → release" step — it is a **membrane-biophysics** event, not a
   receptor event.
3. **Passive organ tropism** via the **protein corona**: the lipid's surface
   chemistry / 3D shape determines which serum proteins adsorb (e.g. ApoE →
   liver LDLR; IgM → spleen), and *that* determines which organ takes up the
   particle. This is the mechanism Su 2026 demonstrates.

→ Modeling target: **structure → property** (potency, organ). No docking. Graph
neural nets (LiON/Chemprop), descriptor QSAR, or MD-shape features (Su).

### Problem B — active-targeting ligand design (your CD-marker idea)

Here you attach a *separate* recognition molecule (small molecule, peptide,
antibody fragment, aptamer) to the **particle surface**, almost always via the
**PEG-lipid** (e.g. DSPE-PEG-maleimide → conjugate ligand). That ligand binds a
cell-surface receptor ectodomain (e.g. a CD marker, TfR, ASGPR, FRα, PSMA) and
triggers **receptor-mediated endocytosis** of the whole particle.

→ Modeling target: **protein–ligand binding**. This *is* a docking / affinity /
contrastive-retrieval problem. This is where DrugCLIP, Boltz, ChEMBL, and the
PPI databases belong.

**Both problems can feed one LNP**: a Problem-A optimized ionizable lipid for
potency + a Problem-B ligand on the PEG for active targeting. But they are built
and validated on separate tracks.

The receptor table in your first screenshot (ASGPR, TfR/CD71, megalin, CD206,
DEC-205, Stabilin, FRα, PSMA, CD19/CD20/BCMA …) is a **Problem-B target list** —
each row is "receptor → what ligand class already hits it → what payload". It is
the shortlist of receptors to design surface ligands against.

---

## 1. What the two attached papers actually give us (Problem A)

### Witten et al. 2024, *Nat Biotech* — "LiON"
- **Model**: directed message-passing NN (D-MPNN) via **Chemprop 1.7.0**.
- **Input**: ionizable-lipid structure (as a molecular graph) **+ metadata**
  (formulation molar ratios, cargo type, target cell/organ).
- **Output**: a single scalar delivery-potency prediction.
- **Data**: >9,000 LNP activity measurements across 20 screens.
- **Reuse**: code is **MIT-licensed** at `github.com/jswitten/LNP_ML`, ships
  `data/all_data.csv` + `main_script.py` (split / train / predict / screen).
  No pretrained checkpoint, but training is a single command and reproducible.
  → **This is our Phase-1 backbone. We can reuse the data and retrain directly.**

### Su et al. 2026, *Nat Biomed Eng* — spatial-conformation ML
- **Library**: 1,408 lipids = 14 amine heads × {ester, amide} linkers × 16 tails.
- **Pipeline**: MD (>2,000 conformations) → 2D density maps → 28 features
  (22 spatial + 6 chemical) → SISSO symbolic regression + ML → predict potency;
  proteomics ties shape → corona → organ.
- **Reuse**: the *concept* (shape features, cone-shape ⇒ good endosomal escape)
  is reusable; full reproduction needs MD (GPU). Extractable now: the library
  enumeration (heads/linkers/tails) and the SISSO feature idea. Treat MD-shape
  features as an **optional Phase-3 enrichment**, not a Phase-1 dependency.

**Bottom line on "extract all their data":** LNP_ML's `all_data.csv` is directly
harvestable and is the highest-value dataset in the project. Su's raw MD data may
not be fully released; we reproduce its *method* only if shape features pay off.

---

## 2. Where every resource you mentioned fits

| Resource | Problem | Role |
|---|---|---|
| LNP_ML (Witten) | A | Training data + baseline D-MPNN model |
| Su 2026 spatial ML | A | Optional MD-shape feature enrichment |
| SMILES embeddings | A & B | Featurize any molecule library (this repo, Phase 1) |
| DrugCLIP / drug-the-whole-genome (#1 PPI, #4 ligand) | **B** | Contrastive pocket↔molecule retrieval: given a receptor ectodomain, retrieve candidate binders from a library |
| ChEMBL / BindingDB (MCP available here) | B | Known actives per target → training/eval pairs, decoys |
| Boltz (MCP available here) | B | Structure + **binding-affinity** co-folding for validation of ligand↔receptor and peptide↔receptor |
| humanPPI (prodata.swmed.edu) | B | Protein–protein interactions; useful for peptide/PPI-derived targeting, secondary |
| Delivery-kinetics ODE/stochastic models (Mihaila 2017/2019, Müller 2024) | A | Mechanistic layer to turn "uptake + escape" into a delivery/expression number — a **complementary simulator**, not an ML model |
| Modal / your GPU | infra | Run Chemprop training, DrugCLIP, Boltz-batch, DiffDock, MD |

---

## 3. Staged roadmap

Each phase is a shippable increment. Start with Phase 1; A and B can then proceed
in parallel.

### Phase 0 — scaffold (this turn)
- Repo structure, featurization module, reference-lipid dataset builder,
  resource inventory, fetch helpers. ✅ (in progress)

### Phase 1 — SMILES embedding + reproduce LiON (Problem A)  ← recommended start
1. Featurize any SMILES library: ECFP/Morgan, MACCS, RDKit 2D descriptors
   (done here); add learned embeddings (Chemprop D-MPNN, MolFormer) later.
2. Clone LNP_ML, load `all_data.csv`, reproduce the train/val/test split, retrain
   the D-MPNN, confirm we match reported correlations.
3. Wrap prediction so we can score **our own** enumerated lipid library.
   Deliverable: `predict_potency(smiles, formulation, target) → score`.

### Phase 2 — build & screen an ionizable-lipid library (Problem A)
1. Enumerate a combinatorial library (heads × linkers × tails), Su-style, as
   SMILES via reaction templates (RDKit `ReactionFromSmarts`).
2. Score with the Phase-1 model; rank; filter by synthesizability + pKa proxy.
3. Optional: add MD-shape features (Su) for the top-k on GPU/Modal.

### Phase 3 — active-targeting ligand library (Problem B)
1. Pick receptors from the screenshot table (start with 1–2, e.g. ASGPR for
   liver, TfR for BBB, or a CD marker of interest). Get the **ectodomain**
   structure (PDB / AlphaFold).
2. Assemble known binders + decoys from ChEMBL/BindingDB per receptor.
3. Retrieve candidate binders from a molecule/peptide library with **DrugCLIP**
   (embedding-based screen) and/or dock with DiffDock/Vina.
4. **Validate** top hits with **Boltz** (structure + affinity co-folding).
   Deliverable: ranked ligand shortlist per receptor, each with a predicted
   binding pose + affinity, ready for PEG-lipid conjugation.

### Phase 4 — delivery-efficiency estimate (mechanistic, both)  ✅ DONE
- `lipidlib/kinetics.py` + `analysis/delivery_kinetics.py`: 4-compartment
  uptake→escape→translation ODE cascade (Müller 2024 / Mihaila 2017-19), with a
  closed-form AUC (numeric match <1%). Reproduces the ~1.8% escape bottleneck and
  ~13 h protein peak. Links to Track A by mapping LiON potency → escape rate →
  expression dynamics (lipid sets amplitude, cargo stability sets timing). See
  [`docs/DELIVERY_KINETICS.md`](DELIVERY_KINETICS.md).

### Modal — verified working (2026-07)

`modal_app/lion.py` runs end-to-end on a Modal A10G GPU. Verified:
`modal run modal_app/lion.py::demo_screen` builds the image, trains a 3-epoch
toy ensemble, and produces a real LiON screen (`pred_file.csv`, 186 rows, with
`cv_0..cv_4_pred_delivery` + `avg_pred_delivery`), persisted to the `lion-models`
volume.

Build gotchas solved (all encoded in the image):
- Chemprop 1.7.0 hard-pins Python `>=3.7,<3.9`; Modal's builder needs `>=3.10`.
  → function runs in 3.10, Chemprop lives in a `uv`-built 3.8 venv, invoked via
  subprocess.
- `uv` venvs omit setuptools → hyperopt's `pkg_resources` fails → install
  `setuptools wheel` into the 3.8 venv.
- RDKit `Draw` needs `libxrender1 libxext6 libsm6 libglib2.0-0` (apt).
- LiON writes screen output to `<split>_preds/<library>/pred_file.csv`.
- This environment's proxy needs `modal[api-proxy-support]` / `python-socks`.

### Two training tiers (both verified on Modal)

| mode | command | folds × epochs | GPU | time | use |
|---|---|---|---|---|---|
| **lite** | `modal run modal_app/lion.py::train_lite` | 1 × 15 | T4 | **~6 min** | fast, cheap; one usable model |
| **full** | `modal run modal_app/lion.py::train --epochs 30` | 5 × 30 | A10G | ~40–60 min | paper-faithful ensemble + error bars |

The lite tier trains on the **full ~13k-point dataset** (single fold, fewer epochs,
cheaper GPU) — ~10× cheaper than the ensemble. Implemented by making LiON's
hardcoded `cv_num = 5` read `LION_CV_NUM` (a one-line sed patch in the image), so
both `train`/`train_lite` and `screen(folds=…)` honour the fold count.

Verified lite run: 1-fold, 15-epoch train on `all_random_split_for_paper` → test
RMSE 0.90 (z-scored delivery) in 6 min on T4; screening the 5 clinical ionizable
lipids ranked MC3 highest in the HeLa/in-vitro context. Screen a lite model with
`screen(..., folds=1)`.

Next Modal step (optional, more GPU-hours): the full 5-fold ensemble + `analyze`
for held-out Pearson/Spearman/Kendall, matching the paper.

### Phase 5 — deploy at scale
- Modal apps for: Chemprop training/inference, DrugCLIP screen, Boltz batch,
  DiffDock, and (if needed) MD. This session is CPU-only + ephemeral, so heavy
  compute is authored here and executed on Modal or your GPU.

---

## 4. Compute reality

- **This session**: CPU-only, ephemeral container. Good for featurization,
  library enumeration, data wrangling, small sklearn baselines, and authoring
  Modal/GPU jobs. **Not** for MD, large NN training, or batch docking.
- **Boltz** and **ChEMBL** are available *right now* as hosted MCP tools — we can
  validate individual designs without any GPU.
- **Your GPU / Modal**: for Chemprop training, DrugCLIP inference over large
  libraries, DiffDock, and MD. All the above tools are open-source and
  Modal-deployable; Phase 5 provides the wrappers.

---

## 5. Open decisions (need your steer)

1. **Priority**: start Problem A (ionizable-lipid potency/organ — the two papers)
   or Problem B (targeting-ligand for a specific CD/receptor)? *Recommendation:
   Phase 1 (A) first — it's the most reproducible and both papers hand us the
   data; B can start in parallel once a target receptor is chosen.*
2. **Compute target for heavy jobs**: your own GPU, or Modal? (Changes how Phase
   5 wrappers are written.)
3. **For Problem B**: which receptor(s) first, from the screenshot table?
