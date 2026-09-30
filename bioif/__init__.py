"""
bioif -- a typed interchange layer for composing heterogeneous biological
models into longer causal chains.

The thesis: cross-model composition fails at the *semantics* of the seam, not
at the data format. So the package supplies a six-field contract (entity,
quantity, context, estimate, domain, provenance), adapters that must declare
their domain and their assumptions, a registry that routes chains
automatically, and an executor that propagates distributions, permits
refusal, and attributes the endpoint's variance back to individual links.

See INTEROP.md for the bottleneck analysis and the ranked plan; run
`python3 -m bioif.demo` for the worked example and `python3 -m bioif.selftest`
for the tests.
"""
from .core import (Claim, Context, Entity, Estimate, Provenance, Quantity,
                   Verdict, check_identity,
                   IN_DOMAIN, EXTRAPOLATE, REFUSE)
from .adapter import Adapter, Coercion, BRIDGE, COERCION, EMPIRICAL
from .registry import Registry
from .chain import ChainResult, Experiment, rank_experiments, run_chain

__all__ = [
    "Claim", "Context", "Entity", "Estimate", "Provenance", "Quantity",
    "Verdict", "check_identity", "IN_DOMAIN", "EXTRAPOLATE", "REFUSE",
    "Adapter", "Coercion", "BRIDGE", "COERCION", "EMPIRICAL",
    "Registry", "ChainResult", "Experiment", "rank_experiments", "run_chain",
]
