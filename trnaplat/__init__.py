"""Platform layer: the three demos, built on `fftplsr` and `pylrs`.

| Module | Demo | What it decides |
|---|---|---|
| `demo1_multisite` | 1 | Does activity engineering pay off, and is the ceiling the synthetase or the tRNA pool? |
| `demo2_scaffold`  | 2 | Is picking a better starting scaffold worth more than mutating one? |
| `demo3_evidence`  | 3 | For a new ncAA, how much should the prediction be trusted? |

Shared convention, inherited from `fftplsr` and `pylrs`: ✅ marks a number this
code computed (stating the data and the split), ⚠️ marks an extrapolation, and
❌ marks a published claim that did not reproduce.
"""
