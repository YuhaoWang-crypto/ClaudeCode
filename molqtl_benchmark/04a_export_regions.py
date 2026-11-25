"""Export real genotype blocks (cis windows of 150 random genes) for the fine-mapping / colocalisation simulations."""
import numpy as np, pandas as pd
from common import *
G, bim, E, genes, sets = load(); rng = np.random.default_rng(1)
eur, yri = (sets.set != "YRI").values, (sets.set == "YRI").values
f_e, f_y = G[eur].mean(0) / 2, G[yri].mean(0) / 2
ok = (np.minimum(f_e, 1 - f_e) >= 0.05) & (np.minimum(f_y, 1 - f_y) >= 0.05)     # common in both ancestries
os.makedirs(f"{WORK}/regions", exist_ok=True); n = 0
for g in rng.permutation(genes.index.values):
    ix = np.where((bim.chr.values == str(genes.chr[g])) & (np.abs(bim.pos.values - genes.tss[g]) <= 200_000) & ok)[0]
    if len(ix) < 200: continue
    if len(ix) > 1000: ix = ix[:1000]
    pd.DataFrame(G[:, ix].astype(np.int8), index=sets.set.values).to_csv(f"{WORK}/regions/r{n:03d}.csv.gz", header=False)
    n += 1
    if n == 150: break
print(n, "regions")
