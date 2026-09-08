#!/usr/bin/env bash
# Fetch the RibonanzaNet / Struct2SeQ source and the model weights.
#
# The official checkpoints live on Kaggle (account required). This pulls the
# HuggingFace mirrors instead and converts the secondary-structure one into the
# original layout; see scripts/convert_rnet_weights.py, which verifies the
# conversion against the base model before trusting it.
#
# Set HF_TOKEN if the HuggingFace endpoints rate-limit anonymous requests.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THIRD_PARTY="$ROOT/third_party"
WEIGHTS="$ROOT/weights"

# Pinned commits: the model code must match the checkpoints.
STRUCT2SEQ_COMMIT="381b7f6aa33e0ba64b436294de4cccbc721ffdcd"
RIBONANZANET_COMMIT="ace422b629f98154f78bb409b7b354a523aa08b8"

clone_pinned() {
  local url="$1" dest="$2" commit="$3"
  if [[ -d "$dest/.git" ]]; then
    echo "have    $(basename "$dest")"
    return
  fi
  echo "clone   $(basename "$dest")"
  GIT_LFS_SKIP_SMUDGE=1 git clone -q "$url" "$dest"
  git -C "$dest" checkout -q "$commit"
}

mkdir -p "$THIRD_PARTY" "$WEIGHTS"
clone_pinned https://github.com/Shujun-He/Struct2SeQ "$THIRD_PARTY/Struct2SeQ" "$STRUCT2SEQ_COMMIT"
clone_pinned https://github.com/Shujun-He/RibonanzaNet "$THIRD_PARTY/RibonanzaNet" "$RIBONANZANET_COMMIT"

fetch_hf() {
  local repo="$1" file="$2" out="$3"
  if [[ -s "$out" ]]; then
    echo "have    $(basename "$out")"
    return
  fi
  echo "fetch   $(basename "$out")"
  local auth=()
  [[ -n "${HF_TOKEN:-}" ]] && auth=(-H "Authorization: Bearer $HF_TOKEN")
  curl -fsSL "${auth[@]}" -o "$out" "https://huggingface.co/$repo/resolve/main/$file"
}

# Primary source: the official checkpoints, re-hosted alongside gRNAde.
fetch_hf chaitjo/gRNAde ribonanzanet/ribonanzanet.pt "$WEIGHTS/RibonanzaNet.pt"
fetch_hf chaitjo/gRNAde ribonanzanet_sec_struct/ribonanzanet_ss.pt "$WEIGHTS/RibonanzaNet-SS.pt"

# Independent copies, used by convert_rnet_weights.py to cross-check the above.
fetch_hf roos23/RibonanzaNet RibonanzaNet.pt "$WEIGHTS/roos23_RibonanzaNet.pt"
fetch_hf multimolecule/ribonanzanet pytorch_model.bin "$WEIGHTS/multimolecule_ribonanzanet.bin"
fetch_hf multimolecule/ribonanzanet-ss pytorch_model.bin "$WEIGHTS/multimolecule_ribonanzanet_ss.bin"

python3 "$ROOT/scripts/convert_rnet_weights.py"
echo "weights in $WEIGHTS"
