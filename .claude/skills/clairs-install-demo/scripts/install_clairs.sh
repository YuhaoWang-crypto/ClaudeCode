#!/usr/bin/env bash
# Install ClairS (HKU-BAL) from source into a self-contained micromamba env.
#
# Why not the Docker image? `hkubal/clairs:latest` is the upstream-recommended
# route, but it needs a Docker daemon. In sandboxes/CI without one, this script
# reproduces the Dockerfile's environment with micromamba instead.
#
# Usage:  install_clairs.sh [PREFIX]        (default PREFIX=$HOME/clairs)
# Result: $PREFIX/env.sh  -- source it to get `run_clairs` on PATH.
set -euo pipefail

PREFIX="${1:-$HOME/clairs}"
ENV_NAME=clairs
CLAIR3_PIN="${CLAIR3_PIN:-1.2.0}"   # MUST be 1.x -- see reference/installation.md
# Upstream docs use http:// for this host; https:// is used here because some
# egress policies allow only 443. Both serve the same files.
BIO8=https://www.bio8.cs.hku.hk

export MAMBA_ROOT_PREFIX="$PREFIX/mamba"
MM="$PREFIX/bin/micromamba"
ENV_DIR="$MAMBA_ROOT_PREFIX/envs/$ENV_NAME"

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }

mkdir -p "$PREFIX/bin" "$PREFIX/models"

# ---------------------------------------------------------------- micromamba
if [[ ! -x "$MM" ]]; then
  say "Fetching micromamba"
  curl -fsSL --retry 4 --retry-delay 2 \
    https://micro.mamba.pm/api/micromamba/linux-64/latest -o "$PREFIX/mm.tar.bz2"
  tar -xjf "$PREFIX/mm.tar.bz2" -C "$PREFIX" bin/micromamba
  rm -f "$PREFIX/mm.tar.bz2"
fi
"$MM" --version

# ------------------------------------------------------------- ClairS source
# Cloned before the env solve so a bad clone fails in seconds, not 10 minutes.
if [[ ! -d "$PREFIX/ClairS" ]]; then
  say "Cloning ClairS"
  git clone --depth 1 https://github.com/HKU-BAL/ClairS.git "$PREFIX/ClairS"
fi

# ----------------------------------------------------------------- conda env
if [[ ! -d "$ENV_DIR" ]]; then
  say "Creating conda env (clair3=$CLAIR3_PIN + pytorch); ~10 min, ~6 GB"
  # clair3 pinned to 1.x on purpose: 1.x is TensorFlow-based and reads the
  # TF-format Clair3 models shipped inside clairs_models.tar.gz. clair3 >=2.0
  # is PyTorch-based and demands *.pt checkpoints that the bundle does not
  # contain -- it fails at runtime with
  #   FileNotFoundError: .../clair3_models/ont_r104_e81_sup_g5015/pileup.pt
  # ClairS itself is PyTorch; both frameworks coexist in this one env.
  "$MM" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    "clair3=$CLAIR3_PIN" pytorch torchinfo tqdm scipy scikit-learn libboost-headers
else
  say "Conda env already present at $ENV_DIR"
fi

say "Compiling C++ helpers (realigner, debruijn_graph)"
# debruijn_graph needs Boost graph headers; upstream's Dockerfile gets them from
# apt (libboost-graph-dev). Here they come from the conda env, hence the -I.
pushd "$PREFIX/ClairS/src/realign" >/dev/null
g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
g++ -std=c++11 -O3 -shared -fPIC -o debruijn_graph -I"$ENV_DIR/include" debruijn_graph.cpp
popd >/dev/null

# ------------------------------------------------------------------- models
# run_clairs resolves models as $CONDA_PREFIX/bin/clairs_models/<platform>/ so
# the tarball has to land exactly there.
if [[ ! -d "$ENV_DIR/bin/clairs_models" ]]; then
  say "Downloading pre-trained models (~832 MB)"
  curl -fsSL --retry 4 --retry-delay 2 \
    -o "$PREFIX/models/clairs_models.tar.gz" "$BIO8/clairs/models/clairs_models.tar.gz"
  mkdir -p "$ENV_DIR/bin/clairs_models"
  tar -xzf "$PREFIX/models/clairs_models.tar.gz" -C "$ENV_DIR/bin/clairs_models"
  rm -f "$PREFIX/models/clairs_models.tar.gz"
fi

# Verdict (--enable_verdict) additionally needs a 230 MB CNV reference bundle and
# a compiled allele counter (builds its own htslib). Off by default because the
# core caller does not need it; opt in with ENABLE_VERDICT=1.
if [[ "${ENABLE_VERDICT:-0}" == "1" && ! -d "$ENV_DIR/bin/cnv_data" ]]; then
  say "Downloading Verdict CNV reference files (~230 MB)"
  mkdir -p "$ENV_DIR/bin/cnv_data"
  curl -fsSL --retry 4 --retry-delay 2 -o "$PREFIX/models/reference_files.tar.gz" \
    "$BIO8/clairs/data/reference_files.tar.gz" &&
    tar -xzf "$PREFIX/models/reference_files.tar.gz" -C "$ENV_DIR/bin/cnv_data" &&
    rm -f "$PREFIX/models/reference_files.tar.gz" ||
    echo "[WARN] Verdict CNV data unavailable; --enable_verdict will not work"
  ( cd "$PREFIX/ClairS/src/verdict/allele_counter" && chmod +x setup.sh &&
    bash setup.sh "$PWD" ) || echo "[WARN] allele_counter build failed; Verdict disabled"
fi

# ---------------------------------------------------------------- activation
cat > "$PREFIX/env.sh" <<EOF
# source this file to use ClairS
export MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX"
export CLAIRS_HOME="$PREFIX/ClairS"
eval "\$("$MM" shell hook --shell bash)"
micromamba activate $ENV_NAME
export PATH="\$CLAIRS_HOME:\$PATH"
EOF

say "Verifying"
"$MM" run -n "$ENV_NAME" "$PREFIX/ClairS/run_clairs" --version || true
cat <<EOF

Installed. To use:

    source $PREFIX/env.sh
    run_clairs --tumor_bam_fn T.bam --normal_bam_fn N.bam --ref_fn ref.fa \\
               --threads 4 --platform ont_r10_dorado_sup_5khz_ssrs --output_dir out

Next: scripts/fetch_demo_data.sh ont <dir>   # download the HCC1395 demo data
EOF
