# Tier 2: running the shortlist through real models

The candidates in this directory were selected by an **uncalibrated structural
triage heuristic**. Nothing here is an affinity prediction yet. These are the
inputs for the tools that can produce one.

Run at least two orthogonal methods and treat disagreement as the signal it is.

## 1. Co-folding (needs a GPU)

AlphaFold3 is the better choice for TCR-pMHC and should be run first if you
have access to the weights:

```bash
python run_alphafold.py --json_path=alphafold3/D_Q30E.json \
  --model_dir=$AF3_WEIGHTS --output_dir=af3_out
```

Boltz inputs are provided as an open-weights fallback:

```bash
boltz predict boltz/ --use_msa_server --out_dir boltz_out
```

**Read this before trusting either.** Two independent 2025 benchmarks on
TCR-pMHC complexes, both including structures deposited after the models'
training cutoffs, found AlphaFold3 clearly ahead, with roughly 90% of models
at acceptable quality or better. The same benchmark reported that Boltz-2 and
tFold-TCR produced **zero** medium- or high-quality models across 20
complexes: they recover the fold (TCR interface RMSD under 2 A) but place the
receptor in the wrong docking geometry, which is precisely the thing interface
engineering depends on. Treat Boltz output here as a structural sanity check,
not as a ranking signal.

Second caveat: **co-folding confidence is not affinity.** ipTM measures
modelling self-consistency, not Kd. A Boltz-2 affinity-head adaptation scored
about r = 0.15 against TCR3d affinity data, and a 2026 analysis argues that the
apparent predictive power of AlphaFold-derived scores tracks how the TCR is
positioned over the pMHC rather than the chemistry of the contacts. Use
co-folding to discard mutants that break the interface, not to rank survivors.

## 2. Rosetta flex ddG (CPU, best-validated physics option)

Resfiles are in `rosetta/`. The protocol capture is not vendored here; get it
from the Kortemme lab repository (`flex_ddG_tutorial`) and point it at the
structure plus each resfile. Expect roughly 35 backrub models per mutation and
a few CPU-hours per mutation for converged numbers.

flex ddG was benchmarked on 1,240 mutations of general protein-protein
interfaces, with its largest gains on small-to-large substitutions and on
antibody-antigen interfaces. **No TCR-specific flex ddG benchmark appears to
exist**, so its accuracy on these interfaces is an extrapolation. TCR-pMHC
interfaces are unusually flat and low-affinity, which is the regime where ddG
methods tend to do worst.

## 3. MM/PBSA (the best TCR-specific energetics number available)

If you can run MD, this is the method with real TCR validation behind it.
Benchmarked against measured mutant affinities in the ATLAS database it reached
Pearson r = 0.9 across 7 mutants with solved structures and r = 0.74 across 35
modelled mutants, with balanced accuracy 0.79 for calling ddG above
0.5 kcal/mol. Those sample sizes are small, and the protocol omitted explicit
interface waters and entropy corrections.

## 4. FoldX BuildModel (fast, orthogonal, less accurate)

```bash
foldx --command=BuildModel --pdb=1ao7.pdb \
  --mutant-file=foldx/individual_list.txt --numberOfRuns=5
```

Repair the structure first (`--command=RepairPDB`). Average over runs. FoldX is
quicker than flex ddG and noticeably less accurate on charged and buried-polar
substitutions.

## 5. Inverse folding for loop redesign

For redesigning whole CDR loops rather than single positions:

```bash
python protein_mpnn_run.py --pdb_path 1ao7.pdb \
  --pdb_path_chains "D E" \
  --chain_id_jsonl proteinmpnn_chains.jsonl \
  --fixed_positions_jsonl proteinmpnn_fixed_positions.jsonl \
  --num_seq_per_target 64 --sampling_temp 0.2 --out_folder mpnn_out
```

On TCR-pMHC interfaces specifically, ProteinMPNN and ESM-IF1 recover 44% and
50% of native CDR3-interface residues respectively, and removing the pMHC from
the input costs about 20% of that, so always keep the pMHC chains as context.
That benchmark had **no wet-lab validation**, and it found that randomly chosen
dissimilar sequences sometimes still scored high model confidence. Confidence
alone is not a filter. ProteinMPNN also optimises sequence for a fixed
backbone, so it will propose sequences that need a different loop conformation
than the one it was given; re-fold anything it suggests.

## 6. The step that is not optional

Whatever survives, screen it for cross-reactivity properly. The structural
screen in this package is a cheap first pass on one crystal structure. The
established cascade is:

1. Determine the recognition motif of the *actual engineered receptor*
   experimentally: a glycine/alanine or combinatorial X-scan, or randomised-CDR3
   yeast selection. Residue-importance rankings from such selections correlate
   with a Rosetta computational alanine scan, so run that too.
2. Search the motif broadly against the human proteome. CrossDome's
   TCR-centred mode, which weights positions by TCR contact hotspots in the
   same spirit as the screen in this package, raised top-50 enrichment of
   validated cross-reactive peptides from 63% to 82% across 16 known cases, and
   moved the titin peptide from rank 27 to rank 6 of 36,000 in the MAGE-A3
   case. Extend the search to common SNP variants and layer tissue expression.
3. Run an unbiased experimental screen. A yeast-display pHLA library of around
   10^8 peptides is the only method that has prospectively caught a peptide
   that killed patients.
4. Confirm function in physiologically realistic cells from the at-risk
   tissue. In the MAGE-A3 case titin was undetectable in cardiac-derived cell
   lines and absent from mouse models; it showed up only in beating
   iPSC-derived cardiomyocytes.

**Do not substitute a sequence-based TCR-epitope predictor for any of this.**
On the task that off-target screening actually requires, scoring the effect of
single amino-acid changes in a peptide, these models perform at chance: across
eight predictors tested on 21 TCRs against 172 single-residue mutants of one
epitope, none achieved a significant correlation with measured binding.
