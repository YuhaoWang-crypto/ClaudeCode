"""Multi-model epitope / immunogenicity automation for DTU Health Tech tools.

This package wraps the *stand-alone command-line* versions of the DTU Health
Tech (former CBS) predictors behind a single orchestrator so that one protein
can be pushed through many models in parallel and the results merged into one
normalized table.

Why CLI and not a REST API: services.healthtech.dtu.dk does NOT expose a
public REST/JSON API. The historical SOAP/WSDL endpoints were retired. The
officially supported way to automate is to download the free-for-academia
stand-alone packages and run them locally. This package is the glue around
those binaries.
"""

__version__ = "0.1.0"
