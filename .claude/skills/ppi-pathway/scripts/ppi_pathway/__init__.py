"""
ppi_pathway — download STRING & BioGRID protein-interaction data and run
omics / pathway analysis on top of it.

Two entry points:

  * string_api  — ONLINE. Talks to the STRING REST API. Only needs `requests`.
                  Best for functional/pathway enrichment of a gene list and for
                  pulling a small interaction network around a few genes.

  * network     — OFFLINE. Loads the bulk STRING / BioGRID flat files (fetched
                  once by `download.py` and cached on disk) into a networkx
                  graph, so you can work genome-wide without hitting the API.

See SKILL.md for the full workflow.
"""

__all__ = ["download", "string_api", "network", "subnetwork", "enrich"]
__version__ = "0.1.0"

# Human by default; every function takes `taxon=` so other species work too.
HUMAN_TAXON = 9606
STRING_VERSION = "12.0"
