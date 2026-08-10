"""Multi-centre reproduction and honest re-appraisal of MRI-based prediction of
CD8+ TIL spatial distribution in breast cancer.

The package is organised around one claim and one doubt:

* the claim -- DCE-MRI radiomics can recover the CD8+ T-cell spatial
  immunophenotype (immune-desert / immune-excluded / inflamed);
* the doubt -- the published AUCs of 0.97-0.99 come from a single centre,
  n=182, with feature selection performed before the train/validation split
  and no external validation.

``mri_cd8_til.experiments.optimism`` measures how much of that AUC the protocol
manufactures on its own; ``mri_cd8_til.cohorts`` records the public data that
can replace it with a multi-centre estimate.
"""

__version__ = "0.1.0"
