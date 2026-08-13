
## Context zero-shot (the challenge's setting)

*2053 knockdowns, mean over 4 held-out cell lines*

### Headline metrics

| model | discrimination | DE overlap@100 | MAE | VCC score | balanced |
|---|---|---|---|---|---|
| control (Δ=0) | 0.5010 | 0.0047 | 0.0570 | 0.009 | -0.020 |
| global mean *(baseline)* | 0.5011 | 0.0716 | 0.0576 | 0.000 | +0.000 |
| naive transfer | 0.6495 | 0.2489 | 0.0613 | 0.163 | +0.139 |
| **ContextTransfer** | 0.6235 | 0.2432 | 0.0551 | 0.158 | +0.158 |

### Supplementary metrics

| model | DE direction | DE LFC ρ | Pearson (effect) | MAE (effect) | VCC score | balanced |
|---|---|---|---|---|---|---|
| control (Δ=0) | 0.0000 | 0.0000 | 0.0000 | 0.0570 | 0.009 | -0.020 |
| global mean *(baseline)* | 0.6783 | 0.2683 | 0.2171 | 0.0576 | 0.000 | +0.000 |
| naive transfer | 0.8191 | 0.4664 | 0.2781 | 0.0613 | 0.163 | +0.139 |
| **ContextTransfer** | 0.8410 | 0.4921 | 0.3278 | 0.0551 | 0.158 | +0.158 |

### Discrimination, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.501 | 0.501 | 0.501 | 0.501 | 0.501 |
| global mean *(baseline)* | 0.502 | 0.501 | 0.501 | 0.501 | 0.501 |
| naive transfer | 0.726 | 0.568 | 0.665 | 0.639 | 0.649 |
| **ContextTransfer** | 0.676 | 0.567 | 0.616 | 0.636 | 0.624 |

### DE overlap@100, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.005 | 0.004 | 0.003 | 0.007 | 0.005 |
| global mean *(baseline)* | 0.050 | 0.058 | 0.129 | 0.050 | 0.072 |
| naive transfer | 0.261 | 0.255 | 0.270 | 0.209 | 0.249 |
| **ContextTransfer** | 0.281 | 0.207 | 0.279 | 0.206 | 0.243 |

### MAE, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.042 | 0.070 | 0.060 | 0.056 | 0.057 |
| global mean *(baseline)* | 0.046 | 0.068 | 0.059 | 0.058 | 0.058 |
| naive transfer | 0.053 | 0.069 | 0.061 | 0.062 | 0.061 |
| **ContextTransfer** | 0.043 | 0.065 | 0.056 | 0.057 | 0.055 |

### Hyperparameters chosen per fold

```
K562     temp=0.1  smooth=0.15  n_neighbors=5  shrink=0.0  gamma=0.25  mod_clip=5.0  beta=0.8  use_global=0.0  rank=80  rank_mix=0.5  unseen_k=25  renorm=0.0
RPE1     temp=0.1  smooth=0.15  n_neighbors=5  shrink=0.5  gamma=0.5  mod_clip=3.0  beta=1.25  use_global=0.0  rank=80  rank_mix=0.5  unseen_k=25  renorm=0.0
HepG2    temp=0.1  smooth=0.15  n_neighbors=5  shrink=0.0  gamma=0.5  mod_clip=5.0  beta=0.8  use_global=0.0  rank=80  rank_mix=0.5  unseen_k=25  renorm=0.0
Jurkat   temp=0.1  smooth=0.0  n_neighbors=5  shrink=0.0  gamma=0.5  mod_clip=3.0  beta=1.0  use_global=0.0  rank=80  rank_mix=0.75  unseen_k=25  renorm=0.0
```