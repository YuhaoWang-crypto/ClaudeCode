#!/bin/bash
# Everything that runs once both benchmark regimes have landed.
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

echo "########## DATA PACKAGE ##########"
$PY -m virtualcell.export
mkdir -p virtualcell/datapackage/results
cp -f results/*.json results/*.txt virtualcell/datapackage/results/ 2>/dev/null || true
rm -f virtualcell/datapackage/results/*.partial.json

echo "########## RESULT TABLES ##########"
$PY -m virtualcell.report

echo "########## ALL DONE ##########"
