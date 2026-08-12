#!/bin/bash
set -u
cd /home/user/ClaudeCode
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/user/venv/bin/python
echo "########## CONTEXT REGIME ##########"
$PY -m virtualcell.run --regime context --tune-eval 300 --out results/context.json
echo "########## DOUBLE-BLIND REGIME ##########"
$PY -m virtualcell.run --regime double --tune-eval 300 --out results/double.json
echo "########## DONE ##########"
