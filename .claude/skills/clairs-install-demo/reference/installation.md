# Installing ClairS

Four routes exist upstream. Pick by what the machine actually has.

## 1. Docker (upstream default)

```bash
docker run -it -v "$IN":"$IN" -v "$OUT":"$OUT" hkubal/clairs:latest \
  /opt/bin/run_clairs -T "$IN/tumor.bam" -N "$IN/normal.bam" -R "$IN/ref.fa" \
  -t 8 -p ont_r10_dorado_sup_5khz_ssrs -o "$OUT"
```

Absolute paths are required for both mounts. Everything (models, Clair3,
compiled helpers) is baked into the image.

## 2. Singularity

```bash
singularity pull docker://hkubal/clairs:latest
singularity exec -B /path:/path clairs_latest.sif \
  /opt/bin/run_clairs -T tumor.bam -N normal.bam -R ref.fa -t 8 -p ilmn_ssrs -o out
```

## 3. Conda / micromamba from source — `scripts/install_clairs.sh`

Use when there is no Docker daemon (most agent sandboxes, many CI runners:
`docker` the client exists, `/var/run/docker.sock` does not). This is the route
validated end-to-end in this repo.

```bash
bash scripts/install_clairs.sh /opt/clairs
source /opt/clairs/env.sh
```

What it does, and why each step matters:

| Step | Detail |
|---|---|
| micromamba | single static binary from `micro.mamba.pm`; no conda install needed |
| env | `clair3=1.2.0 pytorch torchinfo tqdm scipy scikit-learn libboost-headers` from `conda-forge` + `bioconda`, Python 3.10 |
| clone | `git clone --depth 1 https://github.com/HKU-BAL/ClairS.git` |
| compile | `realigner` and `debruijn_graph` in `src/realign/` |
| models | `clairs_models.tar.gz` (832 MB) → `$CONDA_PREFIX/bin/clairs_models/` |

### The clair3 version trap

This is the one that costs an hour. `run_clairs` shells out to
`$CONDA_PREFIX/bin/run_clair3.sh` for germline calling on both BAMs, and points
it at a Clair3 model **inside ClairS's own bundle**:

```
$CONDA_PREFIX/bin/clairs_models/clair3_models/ont_r104_e81_sup_g5015/
    pileup.index  pileup.data-00000-of-00002  ...      # TensorFlow SavedModel
```

- Clair3 **1.x** is TensorFlow-based and reads exactly these. Last 1.x: `1.2.0`.
- Clair3 **2.x** is PyTorch-based and demands `pileup.pt` / `full_alignment.pt`.

So an unpinned `clair3` (today: 2.0.2) installs cleanly, starts running, gets
through candidate extraction, and then dies inside the Clair3 sub-run:

```
FileNotFoundError: [Errno 2] No such file or directory:
  '.../bin/clairs_models/clair3_models/ont_r104_e81_sup_g5015/pileup.pt'
```

Pin `clair3=1.2.0`. Two further notes:

- `clair3=1.0.11` (the version ClairS's changelog names) **cannot be solved
  today** on linux-64: it pins `libcurl 7.88.1` → `openssl 1.1.1w`, which no
  current `pysam` build accepts. 1.2.0 is on openssl 3 and solves fine.
- The env ends up with TensorFlow 2.15 (for Clair3) *and* PyTorch 2.1 (for
  ClairS). That is expected, not a mis-solve.

### Boost

`debruijn_graph.cpp` includes `boost/graph/adjacency_list.hpp`. Upstream's
Dockerfile apt-installs `libboost-graph-dev`; from conda:

```bash
micromamba install -c conda-forge libboost-headers
g++ -std=c++11 -O3 -shared -fPIC -o debruijn_graph \
    -I"$CONDA_PREFIX/include" debruijn_graph.cpp
```

Header-only is enough — nothing links against a Boost library.

### http vs https on bio8.cs.hku.hk

Every upstream `wget` uses `http://www.bio8.cs.hku.hk/...`. The same paths are
served over `https://`, and only the https form gets through egress policies
that restrict outbound traffic to 443. If a model or demo download returns a
bare `403 Forbidden` with no body, that is the proxy, not the file server.

### Verdict (optional, off by default)

`--enable_verdict` classifies calls as germline / somatic / subclonal from a CNV
profile and purity estimate. It needs two extra things:

```bash
ENABLE_VERDICT=1 bash scripts/install_clairs.sh /opt/clairs
```

- `reference_files.tar.gz` (230 MB) → `$CONDA_PREFIX/bin/cnv_data/`
- `src/verdict/allele_counter/setup.sh` — builds alleleCount, which compiles its
  own htslib/libdeflate and needs perl plus zlib/bz2/lzma/curl dev headers
- `scikit-learn` (already in the env spec) — imported by `src/verdict/correct_logr.py`

Recommended only when tumour purity < 0.8.

## 4. Docker build from the repo Dockerfile

`docker build -f ./Dockerfile -t hkubal/clairs:latest .` — same content as
route 1, useful only when patching ClairS itself. Note the Dockerfile is
`FROM ubuntu:16.04` and creates the conda env with an unpinned `clair3`, so a
rebuild today hits the version trap above; add the pin before building.

## Verifying an install

```bash
run_clairs --version
ls "$CONDA_PREFIX/bin/clairs_models" | head          # platform dirs
ls "$CONDA_PREFIX/bin/clairs_models/clair3_models"   # TF Clair3 models
ls "$CLAIRS_HOME/src/realign/realigner" "$CLAIRS_HOME/src/realign/debruijn_graph"
```

Then run `scripts/run_demo.sh ont` — hitting 28 TP / 0 FP / 1 FN on the ONT
quick demo is the real proof the install matches upstream.

## Resource notes

- Env ≈ 4.5 GB, models ≈ 3 GB unpacked (832 MB download), demo data ≈ 100 MB
- Sub-100 GB RAM even for high-coverage whole genomes (fixed in v0.1.0)
- CPU-only is fine; the models are small and run on CPU by default
