# Reproducing Gen-COMPAS on NANMA

## The question

Can the framework of arXiv:2510.24979v1 be reproduced?

The paper describes Gen-COMPAS in prose only. There is no code, no data, no
supplementary information (it is referenced repeatedly but not attached to the
preprint), and no software-availability statement. So "reproduce" cannot mean
"run the authors' pipeline". It can only mean: **rebuild the method from the
description and check that it does what the paper claims, on a system small
enough to run here.**

## Scope

The paper demonstrates Gen-COMPAS on five systems of increasing size. Four of
them need GPU-days to GPU-weeks, and the Trp-cage comparison additionally
needs the 208 µs DESRES reference trajectory, which is available only by
request from D. E. Shaw Research. This machine has four CPU cores and no GPU.

What is reproduced here is the system the paper itself uses to establish the
method: **NANMA (N-acetyl-N'-methylalaninamide, alanine dipeptide) in
vacuum** — AMBER14/ff14SB, 300 K, Langevin dynamics, 2 fs with constrained
X-H bonds. The rare event is the rotation of the backbone dihedral phi
between the negative-phi basin (C7eq / alpha_R) and the positive-phi basin
(C7ax, continuing into an extended region past phi = 180).

Everything in the loop is CV-free, as the paper requires: both networks act
on rotation- and translation-invariant heavy-atom Cartesian coordinates. The
dihedrals phi and psi are used only to say which macrostate a structure is in
and, after the fact, to draw the figures.

## How the method was rebuilt

The paper's Fig. 1 gives the loop; the details below are choices I had to make
because the text does not specify them.

**Generating intermediates.** The paper says a diffusion model trained on the
accumulated sampling "produces intermediate conformations connecting the
states", and cites both DDPM and DDIM. DDIM is only needed if you want a
deterministic, invertible map, so the natural reading is: invert an A
structure and a B structure to their latents, interpolate, and decode. That is
what is implemented — spherical interpolation, which is the right geometry for
Gaussian latents. Measured round-trip fidelity of the inversion is 0.5 Å RMSD,
well inside thermal spread.

**Targeted MD.** Schlitter's TMD restrains the RMSD to a reference. OpenMM
cannot update the reference of an `RMSDForce` nested in a `CustomCVForce`, so
the restraint is written as a harmonic pull toward reference positions that
are re-superposed onto the instantaneous structure at every steering window.
A distance-matrix restraint was tried first and **fails**: the heavy-atom
distance matrix cannot tell a structure from its mirror image, and flipping
phi is close to a reflection of this skeleton, so the restraint was satisfied
without any transition.

**Committor.** Learned by weighted logistic regression on commitment outcomes.
Each frame is labelled with the state *its own future* reaches first, not with
the parent trajectory's first commitment — a run that touches A at 5 ps and B
at 18 ps has frames in between whose future is B.

**Free energy.** The paper reweights the cumulative sampling. Since every
trajectory Gen-COMPAS collects is unbiased, the cleanest estimator is a
Markov-state model: k-means in the same CV-free feature space, a reversible
maximum-likelihood transition matrix, and the resulting stationary
distribution as frame weights. Shooting runs are therefore **fixed-length**
rather than stopped on commitment; stopping on a state would censor the data
and bias the model.

## What had to be fixed to make it work

None of these are in the paper. Each was found by a result not making sense,
and each is recorded here because anyone rebuilding the method will hit it.

1. **A distance-matrix targeted-MD restraint fails outright.** The heavy-atom
   distance matrix is invariant to reflection, and flipping phi is close to a
   reflection of this skeleton, so the restraint was satisfied without any
   transition: the molecule reached the target's distance matrix as its own
   mirror image. The restraint has to be chirality-sensitive.

2. **A single idealised starting geometry is not good enough.** Building
   NANMA at (phi, psi) = (-82, 73) from a Z-matrix and minimising lands in a
   side minimum 11 kcal/mol above the bottom of the basin. Releasing that much
   energy into 22 atoms is comparable to the barrier itself: three of three
   test trajectories started this way left state A within 75-500 ps, while
   trajectories started from a properly equilibrated structure stayed for tens
   of nanoseconds. Seeding now scans the whole Ramachandran plane, minimises
   every point, and takes the deepest minimum inside the target macrostate.

3. **State definitions have to follow the actual landscape.** Defining state A
   as "all negative phi" puts the boundary on the wrong side of a second,
   13 kcal/mol saddle near phi = -127, and state B genuinely wraps across
   phi = 180, so its membership test must be periodic.

4. **Shooting from the end of the targeted-MD path stalls the loop.** With the
   shooting point taken at the end of the steering, every point committed
   decisively to one basin, no sampling reached the separatrix, and the
   committor never acquired the data it needed to improve. Two things were
   wrong. The barrier crossing occupies one or two steering windows out of
   sixty, and a structure on the barrier falls into a basin within about a
   picosecond, so the 2 ps relaxation after releasing the restraint guaranteed
   the point was already committed. Keeping the structures sampled along the
   steering and picking the one whose predicted committor is closest to 1/2
   fixes it: bracketing a crossing along a one-dimensional path is something
   the committor can do reliably even while it is still badly calibrated in
   the full conformational space. That is the bootstrap the loop needs, and
   it moved the mean measured committor at the shooting points from 0.0-0.25
   to 0.48-0.58.

5. **Commitment cores and thermodynamic basins are not the same sets.** A
   boundary placed where it belongs for integrating a free energy sits about
   25 degrees from the saddle, close enough that ordinary libration carries a
   structure still genuinely on the barrier across it inside one saved frame.
   The commitment test needs strict cores, here at least 45 degrees from the
   saddle.

6. **The transition path time sets the shot length.** The phi crossing is
   quasi-ballistic and commitment happens within a few hundred femtoseconds,
   so the initial 20 ps shots spent about 95% of their cost watching the
   molecule sit in a basin while resolving the crossing with roughly one
   frame. Shots of 6 ps saved every 0.1 ps resolve it and cost three times
   less.

## Reference calculations

Two independent references, same force field, same thermostat:

- **Well-tempered metadynamics** biasing (phi, psi) directly — the kind of
  predefined-CV method Gen-COMPAS is meant to avoid, used here only to have a
  converged answer to check against. Two walkers, bias factor 8.
- **Brute-force unbiased MD** — the cost baseline.

## Results

Two independent replicas, five iterations each. All numbers are regenerated
by `src/summarize.py` into `results/summary.txt`; figures are in `figures/`.

### Reference

| quantity | value |
|---|---|
| free-energy difference between the basins | −2.13 kcal/mol |
| barrier over the phi ~ 0 saddle, from state A | 6.16 kcal/mol |
| drift of the phi profile over the last 4 ns of metadynamics | 0.15 kcal/mol |
| transitions in 90 ns of brute-force unbiased MD | 0 |

### Gen-COMPAS

| quantity | replica 1 | replica 2 | reference |
|---|---|---|---|
| unbiased MD used | 8.62 ns | 8.36 ns | — |
| landscape RMSE below 8 kcal/mol | 1.24 | 1.20 | — |
| free-energy difference between basins | +0.06 | +0.17 | −2.13 |
| barrier from state A | 4.89 | 5.65 | 6.16 |
| committor correlation, fresh shooting | 0.836 | 0.825 | — |
| committor error / binomial noise floor | 0.184 / 0.063 | 0.151 / 0.047 | — |

The free-energy landscape is recovered: both basins and the extended region
appear in the right places, and the barrier is 0.5 to 1.3 kcal/mol too low.
The free-energy *difference* between the basins is off by about 2 kcal/mol in
both replicas, consistently in the same direction.

**The committor is the part that holds up best.** Fitted on commitment
outcomes alone, its predictions track the frequencies actually observed, bin
by bin:

| phi window | frames | observed | predicted |
|---|---|---|---|
| [−45, −20) | 949 | 0.00 | 0.00 |
| [−20, −10) | 26 | 0.08 | 0.09 |
| [−10, 0) | 27 | 0.33 | 0.34 |
| [0, 10) | 18 | 0.83 | 0.82 |
| [10, 20) | 24 | 1.00 | 0.99 |
| [20, 55) | 413 | 1.00 | 1.00 |

Tested against 15 fresh trajectories per structure, launched with new
velocities from structures across the whole committor range, the correlation
between predicted and measured committor is 0.83. Structures the method picks
out near phi = −4 measure 0.27 to 0.67 — genuinely on the separatrix.

**The loop does concentrate sampling where it claims to.** The fraction of
shooting points with a measured committor between 0.2 and 0.8 rises from
10-21% at the first iteration to 29-32% by the fifth, in both replicas.

**The generative step earns its place, but not in the way I expected.**
Replacing the diffusion model with straight-line interpolation between the
same end-state pairs, everything else identical:

| | physically realisable | realisable *and* near the saddle |
|---|---|---|
| diffusion model | 96.4% | 3.3% |
| linear interpolation | 40.0% | 12.1% |

The diffusion model wins decisively on realism, and inside the loop 92-99% of
its proposals survive the geometry check against 51% for the interpolator.
But it places *fewer* proposals near the barrier, because it generates on the
manifold it was trained on and that manifold is concentrated in the basins.
Naive interpolation lands closer to the saddle precisely because it ignores
the manifold — it just moves atoms along chords, which is also why most of
its output is not a molecule.

So in this implementation the generative model supplies plausible targets,
and the barrier-region structures come from the targeted-MD steering that
chases those targets, not from the generator itself. That is worth saying
plainly: the paper presents the generative step as what "produces physically
realistic intermediates", and it does, but on this system the intermediates
that matter are found by the steering.

## What this does and does not establish

Reproduced, on NANMA:

- the framework runs as described, with no collective variable in the loop;
- the free-energy landscape is recovered to about 1.2 kcal/mol from 8.5 ns of
  unbiased sampling, against a brute-force run that saw no transition at all
  in 90 ns;
- the committor learned in conformational space is quantitatively calibrated
  and identifies real transition-state structures;
- the iteration measurably concentrates sampling on the separatrix;
- the generative model contributes something a trivial interpolator does not.

Not established:

- **Kinetics.** The Markov model's slowest implied timescale grows
  monotonically with lag and never plateaus, and the two replicas disagree by
  a factor of 50 on the first passage time (77 ns against 3681 ns). The
  passage time is not determined by this data and should not be quoted as a
  result. The paper's kinetic claims are untested here.
- **Quantitative thermodynamics.** A systematic ~2 kcal/mol error in the
  basin free-energy difference, in the same direction in both replicas,
  points at a bias in the Markov-model weighting rather than noise.
- **Transition-state ensemble precision.** The committor ranks structures
  well but does not pin the q = 1/2 level set tightly: structures predicted
  at 1/2 measure 0.60 and 0.67 on average, with a spread that covers most of
  the range.
- **Everything about the protein systems.** Trp-cage, the ribose-binding
  protein and the ADP/ATP carrier carry the paper's actual claims to novelty,
  and none of them is tested here.

One result differs from the paper by way of the system rather than the
method: the paper finds two competing channels for Trp-cage and for the
ribose-binding protein, whereas NANMA in this force field has one. The
reference free energy puts the alternative route, round the far side through
phi = −127, about 13 kcal/mol above the basin against 6 for the phi ~ 0
saddle, so it carries no measurable flux and a two-channel analysis fits
noise.
