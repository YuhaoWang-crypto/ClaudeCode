---
name: enzyme-activity-engineering
description: "Use this skill to computationally engineer an enzyme (or any protein) for improved activity/stability/binding and to validate the candidates with physics. It runs a staged pipeline on Modal GPU: (1) zero-shot single-mutation scoring across ESM-1v, ESM2, ProteinMPNN, and ESM-IF with consensus ranking; (2) structure-based annotation (catalytic core, RNA/DNA/interface distances, burial) from a PDB/cryo-EM structure; (3) epistasis-aware combination-mutation scoring; (4) fold self-consistency via Boltz-2; (5) MM-GBSA binding ΔΔG (OpenMM+AmberTools); (6) QM active-site cluster and QM/MM reaction-barrier attempts (xtb+Psi4); (7) relative binding FEP/TI (perses); plus structural figures and an interactive 3Dmol viewer. Triggers: 'engineer this enzyme', 'improve activity', 'predict mutations', 'ESM/ProteinMPNN/ESM-IF variant effects', 'rank mutations', 'MM-GBSA / FEP / QM-MM binding or barrier', 'is this mutation stabilizing/tighter-binding'. Also use when asked to reproduce or adapt this exact protein-engineering workflow. Requires MODAL_TOKEN_ID/SECRET (GPU) and HF_TOKEN; reads structures from RCSB/NCBI/UniProt."
license: Proprietary.
---

# Enzyme activity engineering (screen → rank → physics-validate)

A staged, honesty-first pipeline. Cheap ML screening first, expensive physics only where it
answers the specific mechanism — and **every physics result passes a sanity gate or is reported
as a negative**. Do not manufacture numbers: if a calculation fails a gate, say so.

## Prerequisites (verify first)
- **Modal GPU**: `MODAL_TOKEN_ID` + `MODAL_TOKEN_SECRET` set. `pip install 'modal[api-proxy-support]' python-socks` (the sandbox proxy needs `python-socks` or every `modal run` errors — see references/LESSONS.md #1).
- **HF_TOKEN** for model weights. Local: `pip install biopython numpy matplotlib`.
- Smoke-test GPU once (see `scripts/` header comments) before long runs.

## Adapting to a new target
Every script has a **`CONFIGURE`** block at the top. For a new protein set: the `SEQ` (1-letter),
the `PDB` id + which chains are protein / RNA / DNA / other, the `CATALYTIC` residue set (never
mutate these), and the mutation shortlist. Get the sequence+structure from RCSB (`files.rcsb.org`),
NCBI E-utilities, or UniProt; **confirm the sequence with the user** if it came from a file you
could not read (numbering must match their construct).

## Phase 0 — identify & fetch (local)
Pull the structure (`curl files.rcsb.org/download/<PDB>.pdb`) and the sequence. Classify chains
(protein vs RNA `A/U/G/C` vs DNA `DA/DT/DG/DC`), find the catalytic residues by 3D-clustering the
conserved acidic/nucleophilic residues around any catalytic metal, and extract a clean single
protein chain for the structure-based scorers.

## Phase 1 — single-mutation scoring (Modal GPU, ~30–60 min)
Run all four, each writes a per-position score matrix (LLR = logP(mut)−logP(wt); + favors mutant):
- `scripts/modal_esm.py` — ESM-1v + ESM2 masked-marginal (sequence PLMs).
- `scripts/modal_mpnn.py` — ProteinMPNN `unconditional_probs` on chain A (structure).
- `scripts/modal_esmif.py` — ESM-IF1 inverse folding on chain A (structure). Pins torch 2.2.1+cu121 for PyG wheels.
Then `scripts/annotate_structure.py` (distances/burial/tags) and `scripts/rank.py` (percentile
consensus across ≥3 methods, catalytic residues excluded). Pick top candidates for **mechanistic
diversity** (fold/yield, RNA affinity, DNA affinity, turnover), not just top score.

## Phase 2 — combinations (Modal GPU, minutes)
`scripts/modal_combo.py` — epistasis-aware ESM: for each combo, score the full mutant sequence
in the *mutant context* and subtract the additive baseline. `scripts/combo_rank.py` ranks and
flags synergy/interference. Check pairwise 3D distances first: sites >12 Å apart are ~additive.

## Phase 3 — fold self-consistency (Boltz-2 via Boltz_API MCP, ~1 min/seq, ~$0.05)
Fold WT and each designed sequence; superpose to the reference chain by **domain** (a global RMSD
can be large just from apo-vs-complex hinge motion — judge per-domain). Conclusion you want:
catalytic domain reproduced (~2 Å) and mutant ≈ WT with unchanged pLDDT → mutations don't break
the fold. This is the RFdiffusion-style design self-consistency check applied to point mutants.

## Phase 4 — binding ΔΔG (choose by rigor)
- **MM-GBSA** (`scripts/modal_md.py`, OpenMM+AmberTools, ~1 hr): explicit-solvent MD of
  protein+nucleic vs mutant, single-trajectory OBC2. **Direction only** — it overestimates
  charge-adding mutations. Split protein/ligand by **atom index**, not residue name (LESSONS #2).
- **FEP/TI** (`scripts/modal_fep.py`, perses): rigorous ΔΔG_bind = ΔG_mut(complex)−ΔG_mut(apo).
  **Blocked without an OpenEye license** in a bare sandbox (LESSONS #5). Charge-neutral mutations
  (e.g. Phe→Tyr) are the clean case; charge-changing (Ala→Lys) needs
  `transform_waters_into_ions_for_charge_changes=True`. Open-source alternative: GROMACS + pmx.

## Phase 5 — catalysis / kcat (hardest; expect to gate-fail on a bare cluster)
- `scripts/modal_qm.py` — QM active-site cluster (GFN2-xTB opt + Psi4 DFT single point). POC only.
- `scripts/build_reactant.py` + `scripts/modal_barrier2.py` — build a pre-cleavage Michaelis
  complex (re-ligate scissile bond, add 2nd metal at the canonical site, in-line nucleophile) and
  run a gated reaction-coordinate scan. **A gas-phase cluster cannot restrain a two-metal site —
  the metals collapse** (LESSONS #4). A trustworthy barrier needs real QM/MM (protein-embedded,
  MD-equilibrated, DFT TS with frequency verification). Only worth it for genuine active-site /
  second-shell mutations (check distance to the scissile atom first, not just to the catalytic residues).

## Phase 6 — figures & interactive viewer
`scripts/make_figures.py` / `make_methods_fig.py` / `make_combo_fig.py` / `make_validation_fig.py`
(matplotlib, real coordinates). Interactive: inline `3Dmol-min.js` (cdnjs) + the PDB into a
self-contained HTML artifact (CSP blocks external scripts, so inline everything). See
references/METHODS.md for the viewer assembly.

## Trustworthiness tiers (state these honestly in any report)
1. **High** — PLM/inverse-folding consensus (screening); Boltz fold self-consistency.
2. **Medium/directional** — MM-GBSA ΔΔG_bind (sign, not magnitude).
3. **Needs dedicated/licensed env** — FEP/TI (OpenEye or GROMACS+pmx); QM/MM kcat barrier.
Final activity confirmation is always the **wet-lab assay**; computation prioritizes and de-risks.

Read `references/LESSONS.md` (the environment/pipeline pitfalls that cost real debugging time) and
`references/METHODS.md` (per-method rationale, gates, and the reaction-modeling details) before running.
