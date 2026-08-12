#!/bin/bash
set -u
cd /home/user/ClaudeCode
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/user/venv/bin/python
until [ -f results/double.json ]; do sleep 20; done
echo "########## CONTEXT ABLATION ##########"
$PY -m virtualcell.run --ablation --n-eval 600
echo "########## STRATIFIED ANALYSIS ##########"
$PY -m virtualcell.analysis
echo "########## FIGURES ##########"
$PY -m virtualcell.figures
echo "########## ALL DONE ##########"
