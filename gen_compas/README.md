# Gen-COMPAS reproduction

An independent, from-scratch implementation of **Gen-COMPAS** (generative
committor-guided path sampling), the framework introduced in

> Chenyu Tang, Mayank Prakash Pandey, Cheng Giuseppe Chen, Alberto Megías,
> François Dehez, Christophe Chipot,
> *Breaking the Timescale Barrier: Generative Discovery of Conformational
> Free-Energy Landscapes and Transition Pathways*, arXiv:2510.24979v1 (2025).

The paper ships no code, no data and no supplementary information, so nothing
here is a port: the pipeline was rebuilt from the description in the main
text and applied to the system the authors themselves use as their
proof-of-concept, **NANMA** (N-acetyl-N'-methylalaninamide, alanine
dipeptide) in vacuum.

## What is and is not reproduced

| Paper system | Status here | Why |
|---|---|---|
| NANMA / trialanine in vacuum (proof of concept) | **reproduced** | tractable on CPU |
| Trp-cage folding, explicit solvent | not attempted | needs a GPU and the 208 µs DESRES reference trajectory |
| Ribose-binding protein, binding-upon-folding | not attempted | ~10⁵ atoms, GPU-days |
| Mitochondrial ADP/ATP carrier in a membrane | not attempted | ~10⁵ atoms, GPU-days |

This machine has four CPU cores and no GPU. The three protein systems are out
of reach by orders of magnitude, so the claim tested here is the *method*, on
the *system the paper itself uses to establish it*.

## The pipeline

Following Fig. 1 of the paper:

0. **Seed.** 1 ns of unbiased MD in each of the two metastable states.
1. **Generate.** A denoising diffusion model (DDPM, cosine schedule) is
   trained on every conformation sampled so far, in a
   translation/rotation-invariant heavy-atom Cartesian representation.
   Intermediates are produced by deterministic DDIM inversion of one
   end-state structure and one product structure followed by spherical
   interpolation of the two latents and DDIM sampling back.
2. **Filter.** A committor network `q(x)`, trained on the accumulated
   commitment outcomes in the same feature space, scores every generated
   intermediate; those near the `q = 1/2` separatrix are kept as targets.
3. **Project.** Targeted MD from a state-A structure and from a state-B
   structure onto each target, with a best-fit heavy-atom restraint whose
   reference is re-superposed every window, then released for a short
   unbiased relaxation. This is what turns a generated geometry into a
   physically realisable one.
4. **Shoot.** Fixed-length unbiased trajectories from each projected
   structure, with fresh Maxwell-Boltzmann velocities.
5. **Repeat.** Everything goes back into steps 1-2.

No collective variable is used anywhere in the loop. `phi` and `psi` appear
only to define the reactant and product macrostates and, afterwards, to plot
the results.

## Files

| file | role |
|---|---|
| `src/build_system.py` | builds ACE-ALA-NME from a Z-matrix and minimises it (AMBER14) |
| `src/common.py` | system, CVs, state definitions, invariant featurisation |
| `src/models.py` | DDPM with DDIM inversion/sampling; committor network |
| `src/mdops.py` | MD engine, targeted MD, committor shooting |
| `src/gen_compas.py` | the iterative loop |
| `src/msm.py` | Markov-state model (reversible MLE) over the collected sampling |
| `src/analysis.py` | free-energy landscape, committor map, TSE, pathways |
| `src/validate_committor.py` | committor test by fresh shooting |
| `src/reference_metad.py` | reference free energy by well-tempered metadynamics |
| `src/reference_bruteforce.py` | brute-force MD cost baseline |
| `src/compare_generators.py` | diffusion model vs. plain linear interpolation |
| `src/summarize.py` | collects every reported number from `results/` |
| `src/figures.py` | figures |

## Reproducing

```bash
pip install numpy scipy matplotlib openmm torch
python3 src/build_system.py
python3 src/reference_metad.py 0 2 25 &     # reference, 2 walkers
python3 src/reference_metad.py 1 2 25 &
python3 src/reference_bruteforce.py 0 120 & # cost baseline
python3 src/gen_compas.py --iterations 5 --seed-ns 1.0 --tag main
python3 src/analysis.py --tag main
python3 src/validate_committor.py --tag main
python3 src/compare_generators.py main
python3 src/figures.py main
python3 src/summarize.py
```

To run the loop with the generative model removed (straight-line
interpolation instead), add `--generator linear`.

Results, numbers and caveats are in `REPORT.md`.
