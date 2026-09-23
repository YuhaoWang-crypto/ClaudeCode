"""
genotox — virtual in-vitro genotoxicity assays, built as four swappable layers.

    chemistry            biology              measurement          decision
    ---------            -------              -----------          --------
    DamageSource  ->  DamageFlux  ->  SignalCore  ->  Readout  ->  call_result
    (damage.py)       (damage.py)     (core.py)      (readouts.py) (doseresponse)

Step 1: the umu test — SOS/umuDC-lacZ core read out as beta-galactosidase
ONPG absorbance, with the growth-factor validity gate.
Step 2: a mammalian GADD45a-GFP reporter line — p53/Mdm2 delayed-feedback
core, GFP readout, relative-cell-density gate.  Reuses step 1's upstream,
readout registry and decision layer unchanged; only the core is new.
Step 3: comet AND micronucleus — one cytogenetic core read by two
instruments at two exposure times, sharing the p53 sub-model with step 2.
Centromere status is tracked so the endpoint can name a mechanism
(clastogenic vs aneugenic), which no transcriptional reporter can.

The layer boundaries are the deliverable.  The damage object is a per-channel
vector rather than a scalar precisely so that later cores (p53/GADD45a
reporter, comet, micronucleus) can read different channels of the same
upstream prediction; see :mod:`genotox.damage` for why that matters.
"""
from .assay import VirtualAssay
from .core import Exposure, SOSCore
from .cytogenetic import CytogeneticCore
from .p53 import P53Core, P53Mdm2
from .damage import (CHANNELS, Compound, DamageFlux, DamageSource,
                     DEMO_COMPOUNDS, QSARSource, TabulatedSource)
from .doseresponse import (COMET, GADD45A_GFP, MN_CBMN, UMU, Protocol,
                           call_result, dose_series, ec_ir, log_doses)
from .readouts import REGISTRY as READOUTS

__all__ = [
    "VirtualAssay", "Exposure", "SOSCore", "P53Core", "P53Mdm2",
    "CytogeneticCore",
    "Protocol", "UMU", "GADD45A_GFP", "COMET", "MN_CBMN",
    "CHANNELS", "Compound", "DamageFlux", "DamageSource", "DEMO_COMPOUNDS",
    "QSARSource", "TabulatedSource",
    "call_result", "dose_series", "ec_ir", "log_doses", "READOUTS",
]
