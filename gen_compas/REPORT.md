# Reproducing Gen-COMPAS on NANMA

*Status: results section filled in from `results/` after the production run.*

## The question

Can the framework of arXiv:2510.24979v1 be reproduced? The paper describes
Gen-COMPAS in prose only: there is no code, no data, no supplementary
information, and no software availability statement. So "reproduce" here can
only mean: **rebuild the method from the description and check that it does
what the paper claims, on a system small enough to run on this machine.**

## Scope

The paper demonstrates Gen-COMPAS on five systems. Four of them
(trialanine aside) need GPU-days to GPU-weeks and, for Trp-cage, the 208 µs
DESRES reference trajectory. This machine has four CPU cores and no GPU.

What is tested here is the system the paper itself uses to establish the
method: **NANMA (alanine dipeptide) in vacuum**, AMBER14/ff14SB, 300 K,
Langevin dynamics, 2 fs with constrained X-H bonds.

The transition studied is the rotation of the backbone dihedral phi between
the negative-phi basin (C7eq / alpha_R) and the positive-phi basin (C7ax and
the extended region beyond it).

## Reference

Two independent references were computed with the same force field:

| quantity | value |
|---|---|
| free-energy difference between the basins | filled from results |
| barrier over the phi ~ 0 saddle, from A | filled from results |
| transitions seen in brute-force unbiased MD | filled from results |

## Result

filled from results

## Caveats

filled from results
