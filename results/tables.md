
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

## Double-blind (knockdown unseen everywhere too)

*2053 knockdowns, mean over 4 held-out cell lines*

### Headline metrics

| model | discrimination | DE overlap@100 | MAE | VCC score | balanced |
|---|---|---|---|---|---|
| control (Δ=0) | 0.5010 | 0.0047 | 0.0570 | 0.005 | -0.020 |
| global mean *(baseline)* | 0.5011 | 0.0640 | 0.0572 | 0.000 | +0.000 |
| naive transfer | 0.5011 | 0.0703 | 0.0572 | 0.003 | +0.002 |
| **ContextTransfer** | 0.5019 | 0.1729 | 0.0565 | 0.045 | +0.044 |

### Supplementary metrics

| model | DE direction | DE LFC ρ | Pearson (effect) | MAE (effect) | VCC score | balanced |
|---|---|---|---|---|---|---|
| control (Δ=0) | 0.0000 | 0.0000 | 0.0000 | 0.0570 | 0.005 | -0.020 |
| global mean *(baseline)* | 0.6578 | 0.2392 | 0.1901 | 0.0572 | 0.000 | +0.000 |
| naive transfer | 0.6746 | 0.2566 | 0.2001 | 0.0572 | 0.003 | +0.002 |
| **ContextTransfer** | 0.7078 | 0.2988 | 0.2156 | 0.0565 | 0.045 | +0.044 |

### Discrimination, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.501 | 0.501 | 0.501 | 0.501 | 0.501 |
| global mean *(baseline)* | 0.501 | 0.501 | 0.501 | 0.501 | 0.501 |
| naive transfer | 0.501 | 0.501 | 0.501 | 0.501 | 0.501 |
| **ContextTransfer** | 0.502 | 0.501 | 0.503 | 0.501 | 0.502 |

### DE overlap@100, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.005 | 0.004 | 0.003 | 0.007 | 0.005 |
| global mean *(baseline)* | 0.052 | 0.050 | 0.110 | 0.044 | 0.064 |
| naive transfer | 0.053 | 0.050 | 0.133 | 0.046 | 0.070 |
| **ContextTransfer** | 0.179 | 0.153 | 0.237 | 0.122 | 0.173 |

### MAE, per cell line

| model | K562 | RPE1 | HepG2 | Jurkat | mean |
|---|---|---|---|---|---|
| control (Δ=0) | 0.042 | 0.070 | 0.060 | 0.056 | 0.057 |
| global mean *(baseline)* | 0.044 | 0.068 | 0.059 | 0.058 | 0.057 |
| naive transfer | 0.044 | 0.068 | 0.059 | 0.058 | 0.057 |
| **ContextTransfer** | 0.042 | 0.068 | 0.058 | 0.057 | 0.056 |

### Hyperparameters chosen per fold

```
K562     temp=0.6  smooth=0.0  n_neighbors=15  shrink=2.0  gamma=0.25  mod_clip=5.0  beta=0.6  use_global=0.35  rank=80  rank_mix=0.0  unseen_k=15  renorm=0.0
RPE1     temp=inf  smooth=0.3  n_neighbors=30  shrink=1.0  gamma=0.75  mod_clip=3.0  beta=0.6  use_global=0.35  rank=20  rank_mix=0.25  unseen_k=15  renorm=0.0
HepG2    temp=0.25  smooth=0.3  n_neighbors=30  shrink=0.0  gamma=0.5  mod_clip=5.0  beta=0.6  use_global=0.15  rank=80  rank_mix=0.0  unseen_k=25  renorm=0.0
Jurkat   temp=inf  smooth=0.3  n_neighbors=30  shrink=1.0  gamma=0.5  mod_clip=2.0  beta=0.6  use_global=0.35  rank=80  rank_mix=0.0  unseen_k=50  renorm=0.0
```

## Context ablation

*VCC score against the number of source cell lines, hyperparameters held fixed*

| held-out line | 1 source line | 2 source lines | 3 source lines |
|---|---|---|---|
| K562 | 0.179 | 0.215 | 0.235 |
| RPE1 | 0.084 | 0.099 | 0.108 |
| HepG2 | 0.119 | 0.146 | 0.165 |
| Jurkat | 0.124 | 0.142 | 0.153 |
| **mean** | **0.126** | **0.150** | **0.165** |

## By measured perturbation strength

```
### performance by measured perturbation strength

stratum                 model                                  n      PDS  ovl@100      dir  pearson    score
-------------------------------------------------------------------------------------------------------------
silent (<10 DE genes)   global mean [challenge baseline]    1263    0.503    0.001    0.633    0.087    0.000
silent (<10 DE genes)   naive transfer                      1263    0.557    0.286    0.861    0.113    0.131
silent (<10 DE genes)   ContextTransfer (ours)              1263    0.543    0.279    0.866    0.142    0.138

weak (10-100)           global mean [challenge baseline]     803    0.504    0.047    0.677    0.226    0.000
weak (10-100)           naive transfer                       803    0.632    0.182    0.801    0.287    0.133
weak (10-100)           ContextTransfer (ours)               803    0.602    0.173    0.821    0.334    0.114

moderate (100-500)      global mean [challenge baseline]     623    0.505    0.127    0.725    0.326    0.000
moderate (100-500)      naive transfer                       623    0.699    0.219    0.805    0.394    0.166
moderate (100-500)      ContextTransfer (ours)               623    0.661    0.241    0.836    0.450    0.160

strong (>500)           global mean [challenge baseline]     511    0.505    0.217    0.730    0.399    0.000
strong (>500)           naive transfer                       511    0.761    0.297    0.822    0.524    0.218
strong (>500)           ContextTransfer (ours)               511    0.716    0.338    0.857    0.582    0.226

### source-line weights chosen from control profiles alone

  target K562     RPE1=0.24  HepG2=0.31  Jurkat=0.45
  target RPE1     K562=0.29  HepG2=0.41  Jurkat=0.30
  target HepG2    K562=0.35  RPE1=0.39  Jurkat=0.26
  target Jurkat   K562=0.48  RPE1=0.27  HepG2=0.25
```