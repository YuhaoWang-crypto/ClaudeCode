"""
GPU half of Module 23 — Geneformer in-silico deletion on Modal.

Run (from the repo root, needs a Modal token):
    modal run -m grn_pipeline.m23_geneformer_modal               # full run
    modal run -m grn_pipeline.m23_geneformer_modal --n-cells 24  # smoke test

It writes grn_pipeline/data/m23_geneformer_raw.json.gz, which
grn_pipeline/m23_geneformer.py turns into statistics + a figure on CPU.

What runs on the GPU (all inside one Modal function):
  1. download Geneformer-V2-104M + its gc104M dictionaries from Hugging Face,
     pinned to a fixed repo revision;
  2. download the 10x Genomics HEK293T scRNA-seq dataset (2,885 cells);
  3. tokenise each cell exactly like geneformer/tokenizer.py (counts / n_counts
     * 1e4 / gene median, rank the non-zero genes, truncate to 4094,
     wrap in <cls> ... <eos>);
  4. for every sampled cell, delete MAPK-graph genes singly and in pairs and
     read out the <cls> cell embedding of the penultimate layer (Geneformer's
     default emb_layer=-1); each real deletion is mirrored by deletions of
     RANK-MATCHED null genes from the same cell (a random non-MAPK gene within
     +/-3 rank positions), the confound-control that Geneformer deletion
     magnitudes need because they scale with rank position;
  5. a no-op re-run of the unperturbed cell measures the numerical noise
     floor of the forward pass.
"""
import gzip
import json
import os
import time

import modal

APP_NAME = "grn-m23-geneformer"
HF_REPO = "ctheodoris/Geneformer"
HF_REVISION = "1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5"    # 2026-09-12 main
MODEL_SUBDIR = "Geneformer-V2-104M"
DICTS = {"token": "geneformer/token_dictionary_gc104M.pkl",
         "median": "geneformer/gene_median_dictionary_gc104M.pkl",
         "mapping": "geneformer/ensembl_mapping_dict_gc104M.pkl"}
DATA_URL = ("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/"
            "293t_filtered_gene_bc_matrices.tar.gz")
INPUT_SIZE = 4096            # V2: 4094 genes + <cls> + <eos>
MAX_SHIFT = 3                # rank window for the matched null genes

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                   "m23_geneformer_raw.json.gz")

app = modal.App(APP_NAME)
cache = modal.Volume.from_name("grn-m23-cache", create_if_missing=True)
CACHE = "/cache"
HF_CACHE = f"{CACHE}/hf"
DATA_TAR = f"{CACHE}/293t_filtered_gene_bc_matrices.tar.gz"
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("torch==2.5.1", "transformers==4.46.3",
                      "safetensors==0.4.5", "huggingface_hub==0.26.2",
                      "numpy==1.26.4", "scipy==1.14.1"))


# ---------------------------------------------------------------- container
def _retry(fn, tries=5):
    for i in range(tries):
        try:
            return fn()
        except Exception as e:                       # flaky CDN / HF hub
            if i == tries - 1:
                raise RuntimeError(f"download failed: {e!r}") from None
            time.sleep(2 ** (i + 1))


def _hf_paths(local_only):
    from huggingface_hub import hf_hub_download, snapshot_download
    kw = dict(revision=HF_REVISION, cache_dir=HF_CACHE,
              local_files_only=local_only)
    dicts = {k: hf_hub_download(HF_REPO, v, **kw) for k, v in DICTS.items()}
    mdir = snapshot_download(HF_REPO, allow_patterns=[
        f"{MODEL_SUBDIR}/*.json", f"{MODEL_SUBDIR}/model.safetensors"], **kw)
    return dicts, os.path.join(mdir, MODEL_SUBDIR)


@app.function(image=image, volumes={CACHE: cache}, timeout=1800)
def prefetch():
    """Download model, dictionaries and dataset once into the shared volume,
    so the GPU shards never hit Hugging Face / the 10x CDN concurrently."""
    import urllib.request
    _retry(lambda: _hf_paths(local_only=False))
    if not os.path.exists(DATA_TAR):
        def get():
            req = urllib.request.Request(DATA_URL,
                                         headers={"User-Agent": "curl/8.5.0"})
            raw = urllib.request.urlopen(req, timeout=300).read()
            with open(DATA_TAR + ".part", "wb") as fh:
                fh.write(raw)
            os.replace(DATA_TAR + ".part", DATA_TAR)
        _retry(get)
    cache.commit()


def _load_dataset(workdir):
    import tarfile
    import numpy as np
    import scipy.io
    with tarfile.open(DATA_TAR) as tf:
        tf.extractall(workdir, filter="data")
    d = os.path.join(workdir, "filtered_matrices_mex", "hg19")
    X = scipy.io.mmread(os.path.join(d, "matrix.mtx")).T.tocsr()  # cells x genes
    genes = [l.split("\t") for l in
             open(os.path.join(d, "genes.tsv")).read().splitlines()]
    ens = np.array([g[0] for g in genes]); sym = np.array([g[1] for g in genes])
    return X.astype(np.float64), ens, sym


def _tokenise(X, ens, tok, med, mapping):
    """Return per-cell token arrays (with <cls>/<eos>) and the collapsed-gene
    normalised matrix used to build them."""
    import numpy as np
    import scipy.sparse as sp
    n_counts = np.asarray(X.sum(1)).ravel()
    coll = np.array([mapping.get(e) for e in ens], dtype=object)
    keep = np.array([c is not None and c in tok and c in med for c in coll])
    ids = sorted(set(coll[keep]))
    col = {c: i for i, c in enumerate(ids)}
    rows = np.where(keep)[0]
    agg = sp.csr_matrix((np.ones(len(rows)), (rows, [col[c] for c in coll[rows]])),
                        shape=(len(ens), len(ids)))          # sums duplicates
    Xc = (X @ agg).tocsr()
    norm = np.array([med[c] for c in ids])
    tokens = np.array([tok[c] for c in ids])
    Xn = sp.diags(1e4 / n_counts) @ Xc @ sp.diags(1.0 / norm)
    Xn = Xn.tocsr()
    cells = []
    for i in range(Xn.shape[0]):
        row = Xn.getrow(i)
        order = np.argsort(-row.data, kind="stable")
        t = tokens[row.indices[order]][:INPUT_SIZE - 2]
        cells.append(np.concatenate([[tok["<cls>"]], t, [tok["<eos>"]]]))
    return cells, Xc, ids, n_counts


@app.function(image=image, gpu="H100", volumes={CACHE: cache},
              timeout=4 * 3600)
def run(genes, n_cells=None, n_null=2, seed=0, batch=32, precision="fp32",
        shard=0, n_shards=1):
    import pickle
    import tempfile
    import numpy as np
    import torch
    from scipy.sparse import csr_matrix as sp_csr
    from scipy.stats import spearmanr
    from transformers import BertModel

    t0 = time.time()
    rng = np.random.default_rng([seed, shard])
    dict_paths, model_dir = _hf_paths(local_only=True)
    dicts = {k: pickle.load(open(v, "rb")) for k, v in dict_paths.items()}
    tok, med, mapping = dicts["token"], dicts["median"], dicts["mapping"]
    model = BertModel.from_pretrained(model_dir,
                                      add_pooling_layer=False,
                                      attn_implementation="sdpa")
    n_layers = model.config.num_hidden_layers
    # Geneformer emb_layer=-1 = hidden_states[n_layers-1] = output of the
    # penultimate block: drop the last block so last_hidden_state is exactly it.
    model.encoder.layer = model.encoder.layer[:n_layers - 1]
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    dtype = {"fp32": torch.float32, "bf16": torch.bfloat16}[precision]
    model = model.to("cuda", dtype=dtype).eval()

    with tempfile.TemporaryDirectory() as wd:
        X, ens, sym = _load_dataset(wd)
    cells, Xc, ids, n_counts = _tokenise(X, ens, tok, med, mapping)

    # MAPK genes -> tokens (symbol -> Ensembl in this dataset -> collapsed id)
    sym2tok = {}
    for g in genes:
        hits = [mapping.get(e) for e in ens[sym == g]]
        hits = [h for h in hits if h in tok]
        if hits:
            sym2tok[g] = tok[hits[0]]
    tok2sym = {t: g for g, t in sym2tok.items()}
    mapk_tokens = set(sym2tok.values())

    # gene-level descriptors on ALL cells (for the baselines)
    col_of = {tok[c]: i for i, c in enumerate(ids)}
    lognorm = sp_csr(Xc.multiply(1e4 / n_counts[:, None])).log1p().tocsc()
    expr = {g: np.asarray(lognorm[:, col_of[t]].todense()).ravel()
            for g, t in sym2tok.items()}
    detect = {g: float((v > 0).mean()) for g, v in expr.items()}
    mean_expr = {g: float(np.asarray(Xc[:, col_of[t]].todense()).mean())
                 for g, t in sym2tok.items()}
    gl = sorted(sym2tok)
    coexpr = {}
    for i, a in enumerate(gl):
        for b in gl[i + 1:]:
            coexpr[f"{a}|{b}"] = float(spearmanr(expr[a], expr[b]).statistic)

    @torch.no_grad()
    def embed(seqs):
        """<cls> embedding of every sequence + full hidden states of seqs[0]."""
        cls_out, full0 = [], None
        for s in range(0, len(seqs), batch):
            chunk = seqs[s:s + batch]
            L = max(len(x) for x in chunk)
            ids_t = torch.zeros((len(chunk), L), dtype=torch.long)
            att = torch.zeros((len(chunk), L), dtype=torch.long)
            for j, x in enumerate(chunk):
                ids_t[j, :len(x)] = torch.as_tensor(x); att[j, :len(x)] = 1
            h = model(input_ids=ids_t.cuda(), attention_mask=att.cuda()).last_hidden_state
            if s == 0:
                full0 = h[0, :len(chunk[0])].float().cpu().numpy()
            cls_out.append(h[:, 0].float().cpu().numpy())
        return np.concatenate(cls_out), full0

    order = np.random.default_rng(seed).permutation(len(cells))[:n_cells]
    order = order[shard::n_shards]
    t_loop = time.time()
    single = {g: {"real": [], "null": []} for g in gl}
    pair = {}
    gemb = {}
    noise = []
    n_pass = 0

    def unit(v):
        return v / np.linalg.norm(v)

    for ci, c in enumerate(order):
        seq = cells[c]
        body = seq[1:-1]
        pos = {tok2sym[t]: p for p, t in enumerate(body) if t in mapk_tokens}
        present = sorted(pos)
        # rank-matched null(s) for each present MAPK gene
        nulls = {}
        for g in present:
            cand = [p for p in range(pos[g] - MAX_SHIFT, pos[g] + MAX_SHIFT + 1)
                    if 0 <= p < len(body) and p != pos[g]
                    and body[p] not in mapk_tokens]
            nulls[g] = list(rng.choice(cand, size=min(n_null, len(cand)),
                                       replace=False)) if cand else []
        pairs = [(a, b) for i, a in enumerate(present) for b in present[i + 1:]]
        jobs = [("orig",), ("noise",)]
        jobs += [("s", g) for g in present]
        jobs += [("sn", g, k) for g in present for k in range(len(nulls[g]))]
        jobs += [("d", a, b) for a, b in pairs]
        for a, b in pairs:
            for k in range(min(len(nulls[a]), len(nulls[b]))):
                if nulls[a][k] != nulls[b][k]:
                    jobs.append(("dn", a, b, k))

        def deleted(job):
            kind = job[0]
            if kind in ("orig", "noise"):
                idx = []
            elif kind == "s":
                idx = [pos[job[1]]]
            elif kind == "sn":
                idx = [nulls[job[1]][job[2]]]
            elif kind == "d":
                idx = [pos[job[1]], pos[job[2]]]
            else:
                idx = [nulls[job[1]][job[3]], nulls[job[2]][job[3]]]
            keep = np.ones(len(body), bool); keep[idx] = False
            return np.concatenate([[seq[0]], body[keep], [seq[-1]]])

        # 'noise' sits in a different batch slot / padding than 'orig'
        seqs = [deleted(j) for j in jobs]
        perm = [0] + list(rng.permutation(np.arange(1, len(seqs))))
        cls_perm, h0 = embed([seqs[p] for p in perm])
        n_pass += len(seqs)
        cls = {jobs[p]: unit(cls_perm[slot]) for slot, p in enumerate(perm)}
        u0 = cls[("orig",)]
        noise.append(float(1 - u0 @ cls[("noise",)]))

        # contextual gene embeddings from the unperturbed pass
        ge = {g: unit(h0[pos[g] + 1]) for g in present}
        for a, b in pairs:
            gemb.setdefault(f"{a}|{b}", []).append(float(ge[a] @ ge[b]))

        disp = {j: v - u0 for j, v in cls.items()}
        for g in present:
            single[g]["real"].append(float(1 - u0 @ cls[("s", g)]))
            nv = [float(1 - u0 @ cls[("sn", g, k)]) for k in range(len(nulls[g]))]
            single[g]["null"].append(float(np.mean(nv)) if nv else None)

        def metrics(da, db, dab):
            na, nb = np.linalg.norm(da), np.linalg.norm(db)
            align = float(da @ db / (na * nb))
            s = da + db
            epi = float((dab - s) @ (s / np.linalg.norm(s)) / (na + nb))
            return align, epi

        for a, b in pairs:
            key = f"{a}|{b}"
            rec = pair.setdefault(key, {"align_real": [], "align_null": [],
                                        "epi_real": [], "epi_null": []})
            al, ep = metrics(disp[("s", a)], disp[("s", b)], disp[("d", a, b)])
            nal, nep = [], []
            for k in range(min(len(nulls[a]), len(nulls[b]))):
                if ("dn", a, b, k) in disp:
                    x, y = metrics(disp[("sn", a, k)], disp[("sn", b, k)],
                                   disp[("dn", a, b, k)])
                    nal.append(x); nep.append(y)
            if nal:
                rec["align_real"].append(al); rec["epi_real"].append(ep)
                rec["align_null"].append(float(np.mean(nal)))
                rec["epi_null"].append(float(np.mean(nep)))
        if ci % 50 == 0:
            print(f"cell {ci}/{len(order)}  passes={n_pass}  "
                  f"{n_pass / (time.time() - t_loop):.1f} seq/s", flush=True)

    return {"meta": {"model": f"{HF_REPO}/{MODEL_SUBDIR}",
                     "revision": HF_REVISION, "dataset": DATA_URL,
                     "n_cells_total": len(cells), "n_cells": int(len(order)),
                     "shard": shard, "n_shards": n_shards,
                     "n_null": n_null, "seed": seed, "precision": precision,
                     "layer": f"hidden_states[{n_layers - 1}] (penultimate), <cls>",
                     "rank_window": MAX_SHIFT, "forward_passes": n_pass,
                     "seconds": round(time.time() - t0, 1),
                     "seq_per_s": round(n_pass / (time.time() - t_loop), 2),
                     "median_genes_per_cell": float(np.median([len(x) - 2 for x in cells]))},
            "genes": gl, "detect": detect, "mean_expr": mean_expr,
            "single": single, "pair": pair, "gene_emb_cos": gemb,
            "coexpr": coexpr, "noise": noise}


# -------------------------------------------------------------------- local
def _merge(parts):
    """Concatenate the per-shard cell-level records."""
    res = parts[0]
    for p in parts[1:]:
        for g, v in p["single"].items():
            for k in ("real", "null"):
                res["single"][g][k] += v[k]
        for key, v in p["pair"].items():
            rec = res["pair"].setdefault(key, {k: [] for k in v})
            for k in v:
                rec[k] += v[k]
        for key, v in p["gene_emb_cos"].items():
            res["gene_emb_cos"].setdefault(key, []).extend(v)
        res["noise"] += p["noise"]
    m = res["meta"]
    m["n_cells"] = sum(p["meta"]["n_cells"] for p in parts)
    m["forward_passes"] = sum(p["meta"]["forward_passes"] for p in parts)
    m["seconds"] = max(p["meta"]["seconds"] for p in parts)
    m["seq_per_s"] = round(sum(p["meta"]["seq_per_s"] for p in parts), 2)
    m.pop("shard")
    return res


@app.local_entrypoint()
def main(n_cells: int = 0, n_null: int = 2, seed: int = 0,
         precision: str = "fp32", n_shards: int = 6, out: str = OUT):
    from grn_pipeline.m11_fibration import (build_expanded_mapk,
                                            fibration_partition)
    G = build_expanded_mapk()
    colour = fibration_partition(G, role_seed=True)
    genes = sorted(G.nodes())
    kw = dict(n_cells=n_cells or None, n_null=n_null, seed=seed,
              precision=precision, n_shards=n_shards)
    prefetch.remote()
    calls = [run.spawn(genes, shard=i, **kw) for i in range(n_shards)]
    parts = [c.get() for c in calls]
    res = _merge(sorted(parts, key=lambda p: p["meta"]["shard"]))
    res["fiber"] = {g: str(colour[g]) for g in genes}
    res["role"] = {g: G.nodes[g]["role"] for g in genes}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with gzip.open(out, "wt") as fh:
        json.dump(res, fh)
    m = res["meta"]
    print(f"wrote {out}: {m['n_cells']} cells, {m['forward_passes']} passes, "
          f"{m['seconds']} s wall on {n_shards} GPUs")
