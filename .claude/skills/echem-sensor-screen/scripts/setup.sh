#!/usr/bin/env bash
# Install the open-source compute stack for fairchem UMA + PySCF + Bader.
set -e
pip install --quiet torch --index-url https://download.pytorch.org/whl/cpu
# SETUPTOOLS_USE_DISTUTILS=stdlib works around a Debian setuptools/antlr4 build quirk
SETUPTOOLS_USE_DISTUTILS=stdlib pip install --quiet fairchem-core ase pyscf rdkit matplotlib
python -c "import torch, ase, fairchem.core, pyscf, rdkit; print('stack OK')"
echo "Set HF_TOKEN (access to facebook/UMA) before running. For Bader:"
echo "  curl -O http://theory.cm.utexas.edu/henkelman/code/bader/download/bader_lnx_64.tar.gz && tar xzf bader_lnx_64.tar.gz"
