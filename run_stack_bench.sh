#!/bin/bash
# Generate Stack predictions for each leave-one-cell-line-out fold.
# ~69 s per split measured, 101 splits per fold => ~2 h per fold.
set -u
cd /home/user/ClaudeCode
export OMP_NUM_THREADS=4 PYTHONPATH=/home/user/ClaudeCode
PY=/home/user/venv/bin/python
BENCH=/home/user/stack_bench

for T in k562 rpe1 hepg2 jurkat; do
  OUT=$BENCH/pred_$T
  if [ -s "$OUT/non-targeting.h5ad" ]; then echo "== have $T"; continue; fi
  echo "########## FOLD $T ##########"
  [ -s "$BENCH/context_$T.h5ad" ] || $PY -m virtualcell.prep_stack_inputs \
      --target $T --n-perts 100 --n-control 500
  PERTS=$(tr '\n' ' ' < "$BENCH/perts_$T.txt")
  cd $BENCH
  $PY -m virtualcell.patch_stack -- stack-generation \
    --checkpoint /home/user/stack_data/bc_large_aligned.ckpt \
    --base-adata context_$T.h5ad --test-adata test_$T.h5ad \
    --genelist /home/user/stack_data/basecount_genes.pkl \
    --gene-name-col gene_name --split-column gene \
    --split-values $PERTS non-targeting \
    --batch-size 16 --num-steps 5 --output-dir pred_$T \
    --prompt-ratio 0.25 --context-ratio 0.4 > gen_$T.log 2>&1
  echo "  fold $T exit=$? files=$(ls pred_$T 2>/dev/null | wc -l)"
  cd /home/user/ClaudeCode
done
echo "########## STACK GENERATION DONE ##########"
