# GPU dispatch recipes (Modal) — Stages 4, 6, 7

All GPU work runs on the user's Modal via `host.compute.create("byoc:modal")` in
the **repl** tool (or `compute_provider` for image builds). Weights are cached on
persistent volumes across sessions — check before re-uploading.

## Antibody campaign (Stage 6) — RFantibody
- Build CUDA-devel image ONLY in a compute_provider cell (cold build > 30-min clamp):
  `nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04`, torch 2.3.1 cu118, dgl 2.4.0+cu118,
  RFantibody + USalign. Weights from `files.ipd.uw.edu/pub/RFantibody/`.
- CLI is CLICK flags, NOT Hydra: `rfdiffusion -t <antigen> -f <framework> -o <prefix>
  -n <N> -l "H1:7,H3:5-13" -h "T218,T222,T226" -w <weights>`; `proteinmpnn -i -o -n -l -w`;
  `rf2 -i -o -w`. Frameworks: `h-NbBCII10.pdb` (VHH, cheaper) or `hu-4D5-8_Fv.pdb` (VH/VL).
- Target PDB in HLT format: antigen chain relabeled to chain T, hotspots as T/K/E residues.
- Cap 25 backbones x 4 seqs = 100 designs. Parse CDR3 with regex `Y[YF]C(.*?)WG.G`.
- pLDDT often recoverable only from the log TAIL (n<100) — report n_recoverable honestly.

## Boltz-2 co-fold (Stages 6,7) — mini-binders + small molecules
- Env image `proteomics_boltz_gpu`, A100-80GB, volume `/root/.boltz:claude-science-boltz-cache`.
- YAML per complex: protein sequence + ligand SMILES (or peptide), affinity head on.
  Read `affinity_pred_value` = log10(IC50 uM), `affinity_probability_binary` for ranking,
  `iptm` for interface confidence (>0.5 pass). 128-atom affinity cap; `--no_kernels` fallback.
- Batch 4 co-folds per job. For covalent warheads, disable the reactive-group SMARTS filter.

## Evo2 genomic scoring (Stage 4)
- Env `genomics_evo2_gpu` (evo2_7b ~22GB VRAM, H100/A100). Weights cached at
  `/cache/hf/hub/models--arcinstitute--evo2_7b` on the proto-cache volume — mount
  `/cache`, set `HF_HOME=/cache/hf HF_HUB_OFFLINE=1`. Score = mean per-token
  log-likelihood; compare synthetic vs dinucleotide-shuffled control + natural positive control.

## Synthetic cassette (Stage 4) — Proto
- `synbio-cassette-designer` skill in the `proto` conda env. Design ARE/state-sensor
  cassette, QC (splice sites, GC, restriction sites, homopolymers), then validate
  naturalness with Evo2.

## General
- Submit via `c.submit_job(command=, intent=, inputs=[...], timeout=)`; results arrive
  as a `compute_done` notification with `featured_files` ready for save_artifacts.
- Egress may be unrestricted (legacy) — offline weights preferred. Always verify weights
  are warm-cached before dispatch to avoid a 15GB re-download.
