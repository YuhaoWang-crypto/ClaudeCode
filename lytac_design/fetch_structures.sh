#!/usr/bin/env bash
# Fetch the experimental structures the stage-0 LTR epitope spec is built from.
#
#   3F6K  sortilin Vps10p + neurotensin (2.0 A)        -- primary sortilin template
#   4PO7  sortilin + neurotensin, excess ligand         -- resolves a SECOND NT site
#   6X3L  sortilin + compound UMJ                       -- progranulin-site ligand
#   5JQ1  ASGPR CRD + galactosamine mimic ZPF + Ca      -- marks the glycan site
#   6UM2  IGF2R full ectodomain + IGF2 (bovine, cryoEM) -- only D6+D11 in one frame
#   1GP0  IGF2R domain 11, human (1.4 A)                -- correct-species D11
#
# Usage: ./fetch_structures.sh [out-dir]
set -euo pipefail

out_dir="${1:-structures_raw}"
mkdir -p "$out_dir"

for pdb_id in 3F6K 4PO7 6X3L 5JQ1 6UM2 1GP0; do
    dest="$out_dir/$pdb_id.cif"
    if [[ -s "$dest" ]]; then
        echo "$pdb_id  already present"
        continue
    fi
    echo -n "$pdb_id  "
    curl -sSf -o "$dest" "https://files.rcsb.org/download/$pdb_id.cif" \
        -w "%{http_code} %{size_download} bytes\n"
done

echo
echo "Next: python3 ltr_epitope_spec.py --struct-dir $out_dir --out-dir specs"
