"""
recorder_pipeline — Disease Recorder biomarker framework, made computable.

Companion to `grn_pipeline`. Where grn_pipeline asks "which network node is
irreducible / near a tipping point", this package asks the IVD-facing question:

    given that a disease process happened, WHERE was it written, HOW LONG does
    the writing survive, and WHAT should it be divided by?

Four modules, each producing computed numbers (not asserted ones):

  r1_kernel        persistence kernels; effective memory of a recorder is set
                   by the carrier's AGE distribution, not its half-life.
                   Reproduces HbA1c's clinical 8-12 week window from the RBC
                   lifespan alone.
  r2_pairing       when does a denominator help? Exact admissibility criterion
                   rho > kappa/2, and the residual bound Var -> (1-rho^2).
  r3_individuality index of individuality + reference change value; which
                   existing analytes are structurally wasted by a population
                   cutoff.
  r4_catalog       machine-readable catalog of recorder classes and candidate
                   assays, scored by an explicit (heuristic) design rubric.

Rigour labels used throughout, following REPORT.md convention:
  [OK]   rigorous  — derived or computed, reproducible from the code
  [LIT]  literature — input value taken from a cited publication
  [HYP]  hypothesis — a proposal this framework generates, not yet evidence
"""

__all__ = ["r1_kernel", "r2_pairing", "r3_individuality", "r4_catalog"]
