"""
genotox — virtual in-vitro genotoxicity assays, built as four swappable layers.

    chemistry            biology              measurement          decision
    ---------            -------              -----------          --------
    DamageSource  ->  DamageFlux  ->  SignalCore  ->  Readout  ->  call_result
    (damage.py)       (damage.py)     (core.py)      (readouts.py) (doseresponse)

Step 1 (implemented): the umu test — SOS/umuDC-lacZ core read out as
beta-galactosidase ONPG absorbance, with the growth-factor validity gate.

The layer boundaries are the deliverable.  The damage object is a per-channel
vector rather than a scalar precisely so that later cores (p53/GADD45a
reporter, comet, micronucleus) can read different channels of the same
upstream prediction; see :mod:`genotox.damage` for why that matters.
"""
from .assay import VirtualAssay
from .core import Exposure, SOSCore
from .damage import (CHANNELS, Compound, DamageFlux, DamageSource,
                     DEMO_COMPOUNDS, QSARSource, TabulatedSource)
from .doseresponse import call_result, dose_series, ec_ir, log_doses
from .readouts import REGISTRY as READOUTS

__all__ = [
    "VirtualAssay", "Exposure", "SOSCore",
    "CHANNELS", "Compound", "DamageFlux", "DamageSource", "DEMO_COMPOUNDS",
    "QSARSource", "TabulatedSource",
    "call_result", "dose_series", "ec_ir", "log_doses", "READOUTS",
]
