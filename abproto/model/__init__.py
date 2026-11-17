"""Model core: the mixed continuous-discrete generative math.

Two pieces, coupled by a shared time t and (in production) a shared network:

  interpolants.py : continuous Stochastic-Interpolant flow matching over
                    all-atom coordinates (backbone + side chains).
  ctmc.py         : discrete Continuous-Time-Markov-Chain sequence generation
                    with EXACT motif hard constraints baked into the rate matrix.

Everything here is a dependency-free numpy reference implementation: the paths,
targets, losses and samplers are mathematically faithful and runnable/testable
today. The learned velocity field and denoiser are passed in as callables, so a
torch network drops straight in (see the ``*Field`` / ``Denoiser`` type aliases).
"""

from .interpolants import (
    LinearInterpolant,
    flow_matching_loss,
    sample_ode,
    coupled_loss,
    structure_first_time,
)
from .ctmc import ConstrainedCTMC, MotifConstraint, DEFAULT_LIABILITIES

__all__ = [
    "LinearInterpolant",
    "flow_matching_loss",
    "sample_ode",
    "coupled_loss",
    "structure_first_time",
    "ConstrainedCTMC",
    "MotifConstraint",
    "DEFAULT_LIABILITIES",
]
