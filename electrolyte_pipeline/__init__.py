"""
electrolyte_pipeline — a runnable, honesty-labelled reproduction of the
methods surveyed in the "电解液计算专题" digest of
Yao et al., Chem. Rev. 2022, 122, 10970 (DOI 10.1021/acs.chemrev.1c00904).

Module map (page of the digest in brackets):

  e0_systems      compositions, force-field spec, experimental reference table     [p4]
  e1_build        packmol box -> OpenFF Sage parametrisation -> OpenMM System      [p4]
  e1_run          minimise / NPT equilibrate / production with convergence checks  [p4]
  e2_rdf_cn       RDF, r_min, density-weighted N(r), CN distribution P(n)          [p5]
  e3_clusters     SSIP / CIP / AGG with explicit criteria, cluster-size network    [p6]
  e4_properties   density, MSD/D, conductivity, Green–Kubo viscosity, dielectric   [p8]
  e5_qc_clusters  Li+–EC / Li+–DME binding energies, multi-start, BSSE             [p2]
  e6_desolvation  electronic-energy scan vs. classical PMF on one coordinate       [p11]
  e7_interface    graphite slab + electrolyte, z-resolved number densities         [p9]
  e8_ml           property regression (random vs grouped split), GFN2 vs DFT       [p12]
  e9_reactive     why reactive MD is NOT run here (scope note, no numbers)         [p10]

Every number written by these modules carries one of three labels:
  ✅ converged   — statistics checked (block errors / integral plateau)
  ⚠️ demo        — pipeline runs, statistics NOT sufficient for a quantitative claim
  ❌ not run     — out of scope for the force field / compute used
"""
