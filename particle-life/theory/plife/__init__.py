"""
plife - the physics of Particle Life.

A faithful re-implementation of the sandbox-science.com/particle-life engine plus
the continuum theory that explains which structures it produces and why.

Modules
-------
``kernels``      the exact force law, its potential and Hankel transform
``matrices``     the platform's interaction-matrix generators + reciprocity algebra
``model``        the N-body integrator (bit-for-bit the platform's update order)
``twobody``      exactly solvable dimer -> origin of self-propulsion
``stability``    continuum linear stability -> which structures, what size, moving or not
``observables``  order parameters + morphology classifier
"""

from . import kernels, matrices, model, observables, stability, twobody  # noqa: F401
from .model import Params, ParticleLife, default_params  # noqa: F401

__all__ = [
    "kernels", "matrices", "model", "observables", "stability", "twobody",
    "Params", "ParticleLife", "default_params",
]
__version__ = "1.0.0"
