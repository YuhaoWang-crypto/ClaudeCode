# Methods — what each predicts, key parameters, and how to read it

## Scoring convention
All single-mutation scores are **LLR = log P(mut) − log P(wt)**; positive favors the mutant.
Normalize each method to a within-method percentile, then average across methods for a consensus.

| Method | Kind | Predicts | Notes |
|---|---|---|---|
| ESM-1v (650M) | sequence PLM, masked-marginal | evolutionary fitness/tolerance | `esm1v_t33_650M_UR90S_1` |
| ESM2 (650M) | sequence PLM, masked-marginal | independent sequence likelihood | `esm2_t33_650M_UR50D` |
| ProteinMPNN | inverse folding (structure) | foldability/packing | `unconditional_probs`, vanilla `v_48_020.pt` |
| ESM-IF1 (142M) | inverse folding (GVP) | structure-conditioned residue pref | teacher-forced logits; torch 2.2.1+cu121 for PyG |

Two sequence + two structure methods → agreement across the pair-of-pairs is the strong signal.
**Structure methods see protein atoms only** — nucleic-acid effects come from the annotation
(distances measured on the experimental complex), not the models.

## What "improved activity" means (map method → mechanism)
- **Fold/yield** (buried/surface, distal): PLM+inverse-folding favorable → more folded/soluble
  enzyme → more active molecules. Validate with Boltz self-consistency + (rigorously) FEP ΔΔG_fold.
- **Guide/substrate affinity** (RNA- or DNA-contact residue): favorable + adds an H-bond/charge to
  the nucleic backbone → tighter binding (lower Km). Validate with MM-GBSA direction, then FEP.
- **Turnover/kcat** (active-site or true second-shell): needs QM/MM barrier — hard; only attempt
  when the residue is close to the *scissile atom*.

## Consensus ranking (`rank.py`)
Exclude catalytic residues. Require ≥3 methods scoring a position. Rank by mean percentile; also
report per-method sign agreement and the worst-method percentile (an agreement floor). Deduplicate
to best-substitution-per-position for a diverse shortlist. Then pick for **mechanistic diversity**.

## Combination FEP-free epistasis (`modal_combo.py`)
For a mutation set S, build the full mutant sequence; for each mutated position i, mask it *in the
mutant context* and read logP(mut_i)−logP(wt_i). Sum = epistasis-aware joint score. Epistasis =
(that sum) − (additive single-mutant sum). Positive = synergy. Sites >12 Å apart are ~additive.

## Boltz-2 fold self-consistency (Boltz_API MCP)
`boltz_start_structure_and_binding` with the protein sequence (num_samples=1, ~$0.05). Download the
CIF, Kabsch-superpose Cα onto the reference chain **per domain**. Report: WT-vs-experiment,
mutant-vs-experiment, mutant-vs-WT (all per catalytic domain), and pLDDT/pTM. A large *global*
RMSD with a small *catalytic-domain* RMSD = inter-domain hinge (apo vs complex), not a fold failure.

## MM-GBSA (`modal_md.py`)
ff19SB / OL21-DNA / OL3-RNA / TIP3P, PME, 4 fs HMR, ~1.5 ns, OBC2 single-trajectory. ΔΔG_bind =
ΔG(mutant) − ΔG(WT). See LESSONS #2 for the atom-index split and the charge-mutation caveat.

## QM cluster + QM/MM barrier (`modal_qm.py`, `build_reactant.py`, `modal_barrier2.py`)
Cluster = catalytic metal(s) + first-shell carboxylates + waters + scissile phosphate (+ the
mutated residue, so its effect is "seen"). GFN2-xTB geometry, Psi4 B3LYP/def2-SVP single points.
Barrier = relaxed scan of the scissile P–O distance (cleavage/re-ligation through the shared
pentacovalent TS). Gate the profile (LESSONS #4). A converged barrier realistically needs true
QM/MM (protein-embedded) — the cluster is a POC / negative-result demonstrator.

## FEP/TI (`modal_fep.py`)
perses `PointMutationExecutor` builds the hybrid topology; `HybridRepexSampler` (12 λ states) +
`MultiStateSamplerAnalyzer` (MBAR). ΔΔG_bind = ΔG_mut(complex) − ΔG_mut(apo), run as separate legs
and subtracted; the complex leg uses the nucleic partner the residue actually contacts (RNA for an
RNA-contact residue, DNA for a DNA-contact one). Blocked without OpenEye (LESSONS #5) — else use
GROMACS+pmx. Always run a 5-iteration smoke test first; check MBAR overlap + fwd/bwd convergence.

## Interactive 3Dmol viewer (self-contained artifact)
1. Trim the PDB to one protomer + the contacted nucleic acids (keeps the artifact ~1 MB).
2. `curl` `3Dmol-min.js` from cdnjs; inline it (`<script>…</script>`) and inline the PDB in a
   `<script type="text/plain">`. Embed candidate data as a JS array.
3. Style: protein cartoon; RNA/DNA cartoon+stick (semantic colors); catalytic residues red sticks;
   candidate sites lime spheres+labels; click a card → `zoomTo({chain,resi})`.
4. Validate headless (Playwright, `/opt/pw-browsers/chromium`), capture console errors, then publish.
Reference implementation of the assembly is in the project's `analysis/` history (build step that
replaces `__PDBDATA__` / `__THREEDMOL__` / `__CANDIDATES__` placeholders in an HTML template).
