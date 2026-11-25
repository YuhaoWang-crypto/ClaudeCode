"""Pre-train scVI -> scANVI on the reference, then map each unseen query with scArches surgery.

Query labels are never shown to the model. Outputs per-cell predictions + uncertainty scores.
"""
import os, sys, time, numpy as np, pandas as pd, anndata as ad, torch, scvi
from sklearn.neighbors import KNeighborsClassifier
from common import *

torch.set_num_threads(4); scvi.settings.seed = 0
EP_SCVI, EP_SCANVI, EP_QUERY = int(os.environ.get("EP_SCVI", 40)), 15, 30
mdir = f"{WORK}/models"; os.makedirs(mdir, exist_ok=True)
kw = dict(accelerator="cpu", enable_progress_bar=False, check_val_every_n_epoch=None, train_size=1.0)

ref = ad.read_h5ad(f"{WORK}/ref.h5ad")
if not os.path.exists(f"{mdir}/scanvi/model.pt"):
    t = time.time()
    scvi.model.SCVI.setup_anndata(ref, batch_key="batch", labels_key="fine")
    if os.path.exists(f"{mdir}/scvi/model.pt"):          # resume after an interrupted run
        m = scvi.model.SCVI.load(f"{mdir}/scvi", ref)
    else:
        m = scvi.model.SCVI(ref, n_layers=2, n_latent=30, gene_likelihood="nb", encode_covariates=True,
                            use_layer_norm="both", use_batch_norm="none", dropout_rate=0.2)
        m.train(max_epochs=EP_SCVI, batch_size=256, **kw)
        m.save(f"{mdir}/scvi", overwrite=True)
    print("scVI done", time.time() - t, flush=True)
    lm = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown", labels_key="fine")
    lm.train(max_epochs=EP_SCANVI, n_samples_per_label=100, batch_size=256, **kw)
    lm.save(f"{mdir}/scanvi", overwrite=True)
    print("scANVI done", time.time() - t, flush=True)
lm = scvi.model.SCANVI.load(f"{mdir}/scanvi", ref)
z_ref = lm.get_latent_representation(ref)
ref_pred = lm.predict(ref)
print("reference training accuracy (fine):", (ref_pred == ref.obs.fine.values).mean(), flush=True)
knn = KNeighborsClassifier(n_neighbors=15, weights="distance").fit(z_ref, ref.obs.fine.values)
# reference-internal kNN distance scale (leave-self-out) for calibrating the novelty score
d_ref = knn.kneighbors(z_ref, n_neighbors=16)[0][:, 1:].mean(1)

for name in ["q_triana_wta", "q_ts_bm"]:
    t = time.time()
    q = ad.read_h5ad(f"{WORK}/{name}.h5ad")
    truth = q.obs[["fine", "coarse", "batch"]].copy()
    q.obs["fine"] = "Unknown"                       # hide labels from the model
    scvi.model.SCANVI.prepare_query_anndata(q, lm)
    qm = scvi.model.SCANVI.load_query_data(q, lm, freeze_dropout=True)
    qm.train(max_epochs=EP_QUERY, plan_kwargs=dict(weight_decay=0.0), batch_size=256, **kw)
    z = qm.get_latent_representation(q)
    soft = qm.predict(q, soft=True)
    out = truth.copy()
    out["scanvi_fine"] = soft.idxmax(1).values
    out["scanvi_maxp"] = soft.max(1).values
    coarse_p = soft.T.groupby(soft.columns.map(FINE2COARSE)).sum().T
    out["scanvi_coarse"] = coarse_p.idxmax(1).values
    out["scanvi_coarse_maxp"] = coarse_p.max(1).values
    out["knn_fine"] = knn.predict(z)
    out["knn_dist"] = knn.kneighbors(z)[0].mean(1) / np.median(d_ref)
    out.to_csv(f"{OUT}/pred_scanvi_{name}.csv")
    np.save(f"{WORK}/z_{name}.npy", z)
    print(name, "mapped", time.time() - t, "coarse acc (incl. held-out):",
          (out.scanvi_coarse == out.coarse).mean(), flush=True)
np.save(f"{WORK}/z_ref.npy", z_ref)
