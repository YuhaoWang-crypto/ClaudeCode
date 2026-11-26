"""Predict mean ribosome load (MRL) with the pretrained Optimus 5-Prime CNN (Sample et al. 2019),
on its own held-out random-library test set and on real human 5'UTRs it never trained on."""
import os, numpy as np, pandas as pd
from tensorflow import keras
N = os.environ["NB_DATA"]; OUT = "results"
nuc = {'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1], 'n': [0, 0, 0, 0]}
def onehot(seqs, L=50):
    V = np.zeros((len(seqs), L, 4), np.float32)
    for i, s in enumerate(seqs):
        s = s.lower()[:L].rjust(L, 'n') if len(s) < L else s.lower()[:L]
        for j, x in enumerate(s): V[i, j] = nuc.get(x, [0, 0, 0, 0])
    return V
m = keras.models.load_model(f"{N}/repo/modeling/saved_models/main_MRL_model.hdf5", compile=False)

# random library: top 280k by reads, test = top 20k (Sample's split); model trained on the remaining 260k
rnd = pd.read_csv(f"{N}/GSM3130435_egfp_unmod_1.csv.gz")
rnd = rnd.sort_values("total_reads", ascending=False).reset_index(drop=True).iloc[:280000]
rtest = rnd.iloc[:20000].copy()
rtest["pred_cnn"] = m.predict(onehot(rtest.utr.values), verbose=0).ravel()
rtest[["utr", "rl", "total_reads", "pred_cnn"]].to_csv(f"{OUT}/pred_random_test.csv", index=False)

# human + designed library (never seen by the model)
des = pd.read_csv(f"{N}/GSM3130443_designed_library.csv.gz")
des = des[des.utr.str.len() == 50].copy()
des["pred_cnn"] = m.predict(onehot(des.utr.values), verbose=0).ravel()
des[["utr", "rl", "r_total", "library", "pred_cnn"]].to_csv(f"{OUT}/pred_designed.csv", index=False)
print("random test", len(rtest), "| designed", len(des), "| library counts:", des.library.value_counts().to_dict())
