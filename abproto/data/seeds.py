"""A small seed list of known antibody-antigen complexes for the Stage 1 demo.

Mix of soluble/ectodomain and GPCR targets, including the paper's GPCR refs.
The classifier validates each; non-complexes are skipped and logged, so a wrong
guess here is harmless. For a full dataset use rcsb.search_antibody_antigen().
"""

SEED_COMPLEXES = [
    # soluble / ectodomain
    "1N8Z",  # trastuzumab Fab / HER2
    "1VFB",  # D1.3 Fv / lysozyme
    "1NCA",  # NC41 / influenza neuraminidase
    "3HFM",  # HyHEL-10 / lysozyme
    "1MLC",  # D44.1 / lysozyme
    "1FDL",  # HyHEL-5 / lysozyme
    "2DQJ",  # HyHEL-10 variant / lysozyme
    "1YQV",  # HyHEL-5 / lysozyme
    "3GBN",  # CR6261 / influenza HA
    "4FQI",  # CR9114 / influenza HA
    "1BJ1",  # bevacizumab Fab / VEGF
    "1CZ8",  # VEGF / Fab
    # GPCR (hard targets; paper class)
    "8U4S",  # REGN7663 Fab / CXCR4  (confirmed)
    "8TLM",  # Fab / class A GPCR  (receptor identity: verify — RCSB title is generic)
]
