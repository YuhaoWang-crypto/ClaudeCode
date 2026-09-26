# Boltz-2 configs (not submitted)

Each file is a request body for the named Boltz MCP tool. Always run the
matching `boltz_estimate_*` tool with the same body first and check the cost.

Order of execution:
  1. stage1_screen_*        - does anything we already have touch the epitope?
  2. stage2_design_*        - constrained de novo nanobody generation
  3. stage3_cdr3_redesign_* - parent-anchored CDR3 redesign (preferred path)

The critical difference from the earlier campaign: every target block sets
`epitope_residues` AND `non_binding_residues`. The earlier runs applied no
epitope or contact constraints and their replicate contact sets did not
reproduce (Jaccard 0.012). Treat any run whose replicate contact sets disagree
as a failed run, not as a low-affinity result.

Acceptance gate before believing any ranking:
  * two independent samples of the SAME input must agree on the contact set
    (target-residue Jaccard >= 0.6), and
  * the contacted residues must lie within the intended epitope.
Only then are the relative scores worth reading.
