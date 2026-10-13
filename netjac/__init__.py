"""
netjac — a network-Jacobian diagnostics platform.

Turns a curated kinetic model (any SBML, e.g. a BioModels entry) into the
three objects that network-dynamics diagnostics need, and keeps them apart:

  * the STATE Jacobian      J = d(dx/dt)/dx      -> eigenvalues, stability class,
                                                      relaxation time, fragility
  * the ELASTICITY matrix   eps = dv/dx          -> the "edges" (enzyme <- metabolite,
                                                      including allosteric regulation)
  * the CONTROL matrix      C^J = dlnFlux/dlnE   -> the user's "dFlux_i/dEnzyme_j"
                                                      (metabolic control analysis)

They are linked exactly by  C^S = -B (B^T N eps B)^{-1} B^T N,  C^J = I + eps C^S.

Modules
  core          KineticModel: SBML -> steady state, J, eps, C^J, theorem checks
  spectrum      eigen-classification, relaxation time, eigenvalue sensitivities
                (edge / enzyme fragility) with finite-difference verification
  fluctuation   linear-noise covariance (Lyapunov eq.) forward + inverse
                (recover the Jacobian from metabolomics-style covariance)
  cohort        virtual cohorts (inter-individual enzyme variation + assay noise)
  relax         perturbation-relaxation: pulse response, decay-rate fitting
  discriminate  AUC with bootstrap CIs, cross-validated classifiers
  data.cgmacros real CGM data (PhysioNet CGMacros) via HTTP range reads
  statmech      micro -> macro: random-matrix nulls, collective modes, shrinkage
                couplings, all-pairs rewiring, fluctuation-response, MRA
  demos/        D1 edge biomarkers, D2 perturbation-relaxation, D3 fragility
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures", "netjac")
CACHE = os.path.join(ROOT, "figures", "_sbml_cache")
DATA_CACHE = os.path.join(ROOT, "figures", "_data_cache")
