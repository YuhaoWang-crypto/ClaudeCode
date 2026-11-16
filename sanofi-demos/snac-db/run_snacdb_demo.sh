#!/usr/bin/env bash
# SNAC-DB demo: install env, run unit tests, curate Ab/Nb complexes from real PDBs,
# then structurally cluster their antigens with Foldseek.
#
# Requires: micromamba (or conda/mamba). Adjust MAMBA below if needed.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HERE/work}
MAMBA=${MAMBA:-micromamba}
mkdir -p "$WORK" && cd "$WORK"

# 1. Clone
[ -d SNAC-DB ] || git clone --depth 1 https://github.com/Sanofi-Public/SNAC-DB

# 2. Environment: Python 3.10 + ANARCI, HMMER, Foldseek (conda-only tools)
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$WORK/mamba}
"$MAMBA" create -y -n snacdb -c conda-forge -c bioconda \
    python=3.10 anarci hmmer foldseek pip
RUN="$MAMBA" ; ENV_BIN="$MAMBA_ROOT_PREFIX/envs/snacdb/bin"
export PATH="$ENV_BIN:$PATH"

# 3. Install package + patch ANARCI to tolerate unknown 'X' residues
cd SNAC-DB
pip install -r requirements.txt .
snacdb-patch-anarci

# 4. Unit tests
python -m pytest -q unit_tests/unit_test_pipeline.py

# 5. Download a handful of real structures and curate them
PDBS="$WORK/pdbs"; mkdir -p "$PDBS"; cd "$PDBS"
for id in 1N8Z 7BWJ 6WAQ 5JDS 1MLC 7KGK 6M0J; do
    [ -f "$id.pdb1" ] || { curl -sSL -o "$id.pdb1.gz" "https://files.rcsb.org/download/$id.pdb1.gz" && gunzip -f "$id.pdb1.gz"; }
done
cd "$WORK/SNAC-DB"
bash data_curation.sh "$PDBS" True   # True = keep intermediate outputs

# 6. Structurally cluster antigens of the curated complexes with Foldseek
bash finding_hits.sh "$PDBS"_curated antigen

echo "Done. Curated complexes: $PDBS"_curated
echo "Curation summary:        $WORK/pdbs_curation_summary.csv"
echo "Antigen clusters:        $WORK/pdbs_curated_antigen_all_cluster.tsv"
