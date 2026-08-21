"""Personalized neoantigen selection and mRNA-construct assembly.

Public-workflow reimplementation of the individualized neoantigen-therapy
pipeline (the class of product exemplified by intismeran autogene /
mRNA-4157/V940): tumor+normal DNA and tumor RNA -> somatic variants ->
HLA type -> presentation prediction -> ranked, constrained selection of
<=34 neoantigens -> one concatemeric mRNA construct.

The proprietary parts of any commercial product (its exact scoring model,
training data, thresholds and construct rules) are NOT reproduced here.
This package supplies an explicit, editable, literature-grounded stand-in
for those, so the selection layer can be inspected and benchmarked.
"""

__version__ = "1.0.0"
