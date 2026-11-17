"""Stage 1 data pipeline: build a leakage-controlled antibody-antigen dataset.

Primary source is the RCSB PDB (SAbDab is derived from it). The new SAbDab site
is a single-page app with no stable TSV endpoint, so we go to RCSB directly via
its GraphQL batch API (verified reachable) and keep a loader for a user-supplied
SAbDab summary dump.

Flow:
    seed list / RCSB search
      -> rcsb.fetch_entries        (dates, resolution, sequences)
      -> classify.to_complexes     (heavy / light / antigen; GPCR flag)
      -> dedup.cluster_sequences   (identity clustering; mmseqs2 in production)
      -> splits.make_splits        (time split + homology split, cluster-disjoint)
      -> dataset manifest + stats
"""

from .rcsb import fetch_entries
from .classify import to_complexes, AbAgComplex
from .splits import make_splits

__all__ = ["fetch_entries", "to_complexes", "AbAgComplex", "make_splits"]
