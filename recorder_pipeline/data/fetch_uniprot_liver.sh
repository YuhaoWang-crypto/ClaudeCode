#!/bin/sh
# Regenerate the R8 inputs.
# 1) UniProt sequences for the site audit
for acc in P02461 P02768 P02647 P14136; do
  curl -sSL "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "uniprot/${acc}.fasta"
done
# 2) fibrotic-liver mesenchyme co-expression (CELLxGENE Census 2025-11-08)
#    datasets: 02792605-... and 1873a18a-... (PSC + PBC, 13 donors, no normal arm)
#    see r8_audit.part_d; the aggregated outputs are liver_mes_expr.tsv and
#    liver_coexpr.tsv, produced by filtering tissue_general=='liver',
#    is_primary_data==True, cell_type in {hepatic stellate cell, fibroblast}.
