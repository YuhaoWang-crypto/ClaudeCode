#!/usr/bin/env bash
# Demo: the ProteinMPNN kit from anthropics/uplifting-biomolecular-modeling, on a CPU-only Linux box.
# Clones the kit tree, builds a Python 3.11 venv (pinned lock, CPU torch 2.5.1 instead of the cu124 wheel),
# fetches the pinned ProteinMPNN weights, then runs: check, stock design (off), exact (refused without GPU),
# exact --allow-partial, and the kit's unit tests.
set -uo pipefail

KITS="${KITS:-$HOME/uplifting-biomolecular-modeling}"
VENV="${VENV:-$HOME/venv-mpnn}"
export MPNN_DIR="${MPNN_DIR:-$HOME/ProteinMPNN}"
export PYTHONHASHSEED=0

[ -d "$KITS" ] || GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 https://github.com/anthropics/uplifting-biomolecular-modeling "$KITS"
cd "$KITS/proteinmpnn"

if [ ! -x "$VENV/bin/python" ]; then
  uv venv --seed --managed-python --python 3.11.5 "$VENV"
  # every pin of the lock except the CUDA wheels; torch 2.5.1 comes from the CPU index instead
  grep -v -E '^(torch @|nvidia-|triton)' environment/requirements-mpnn.lock > "$VENV/lock-cpu.txt"
  "$VENV/bin/python" -m pip install -q --no-deps -r "$VENV/lock-cpu.txt"
  "$VENV/bin/python" -m pip install -q --no-deps torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu
fi
. "$VENV/bin/activate"

bash run.sh install --weights "$MPNN_DIR"

X="--input $MPNN_DIR/inputs/PDB_monomers/pdbs"
S="--seed 37 --num_seq_per_target 8 --batch_size 8"

echo "== 1. check (dry run): expect exit 3, no CUDA device"
bash run.sh check --config h100 --mode exact; echo "exit=$?"

echo "== 2. stock ProteinMPNN (--mode off)"
bash run.sh design --config h100 --mode off $X --out out/off $S; echo "exit=$?"

echo "== 3. exact mode: refused on CPU (exit 3), never falls back silently"
bash run.sh design --config h100 --mode exact $X --out out/exact $S; echo "exit=$?"

echo "== 4. exact --allow-partial: the CPU-capable levers only"
bash run.sh design --config h100 --mode exact --allow-partial $X --out out/partial $S; echo "exit=$?"
for f in 5L33 6MRR; do
  cmp -s out/off/seqs/$f.fa out/partial/seqs/$f.fa && echo "$f: identical to stock" || echo "$f: differs from stock"
done

echo "== 5. kit unit tests"
python -m pip install -q pytest
python -m pytest -q opt/proteinmpnn_opt/tests
