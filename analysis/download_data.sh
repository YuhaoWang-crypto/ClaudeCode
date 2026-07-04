#!/usr/bin/env bash
# Download the Sade-Feldman et al. (Cell 2018) melanoma scRNA-seq data
# from GEO series GSE120575. Files are large and are NOT committed
# (see ../.gitignore); run this once before analysis/run_scrna_analysis.py.
set -euo pipefail

DEST="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$DEST"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl"

echo "Downloading per-cell clinical metadata (small)..."
curl -fsSL -o "$DEST/GSE120575_patient_ID_single_cells.txt.gz" \
  "$BASE/GSE120575_patient_ID_single_cells.txt.gz"

echo "Downloading single-cell TPM matrix (~120 MB gz)..."
curl -fsSL -o "$DEST/GSE120575_TPM.txt.gz" \
  "$BASE/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"

# --- External validation cohorts ---
echo "Downloading Riaz anti-PD-1 bulk RNA-seq (GSE91061)..."
RIAZ="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061"
curl -fsSL -o "$DEST/GSE91061_fpkm.csv.gz" \
  "$RIAZ/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
curl -fsSL -o "$DEST/GSE91061_series_matrix.txt.gz" \
  "$RIAZ/matrix/GSE91061_series_matrix.txt.gz"

echo "Downloading Yost anti-PD-1 scRNA + TCR (GSE123813)..."
YOST="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl"
for f in bcc_tcell_metadata.txt.gz bcc_tcr.txt.gz bcc_scRNA_counts.txt.gz; do
  curl -fsSL -o "$DEST/GSE123813_$f" "$YOST/GSE123813_$f"
done

echo "Extracting Yost signature-gene rows + per-cell library sizes..."
zcat "$DEST/GSE123813_bcc_scRNA_counts.txt.gz" | awk -F'\t' '
  BEGIN{ split("TCF7 IL7R SELL CCR7 CD28 LEF1 CD27 GZMK PDCD1 HAVCR2 LAG3 TIGIT CTLA4 ENTPD1 CD38 TOX CD8A CD8B MKI67 GZMB GNLY", g, " ");
         for(i in g) want[g[i]]=1; }
  NR==1{ print > "'"$DEST"'/yost_sig_genes.tsv"; next }
  { for(j=2;j<=NF;j++) col[j]+=$j; ncol=NF; if($1 in want) print >> "'"$DEST"'/yost_sig_genes.tsv"; }
  END{ for(j=2;j<=ncol;j++) printf "%d\n", col[j] > "'"$DEST"'/yost_libsize.tsv" }'

echo "Done. Files in $DEST:"
ls -lh "$DEST"
