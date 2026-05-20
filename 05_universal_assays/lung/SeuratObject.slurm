#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --job-name=LungUniAssay
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load r/4.4.2

# 1. Paths
UP="/QRISdata/Q8448/Mouse_disease_data/universal_peaks_v5/consensus_regions_v5.bed.gz"
OBJ="/QRISdata/Q8448/Mouse_disease_data/Lung/lung_integrated_clean_annotated.rds"
OUT="/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_v5.rds"
FRAG_TPL="/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/{id}/outs/fragments.tsv.gz"
SCRIPT="/home/s4869245/scripts/integration_analysis/Lung_seurat_obj/SeuratObject.R"
SAMPLE_KEY="folder_id"

if [[ ! -f "$OBJ" ]]; then
  echo "[ERROR] Object not found: $OBJ"
  exit 1
fi

mkdir -p "$(dirname "${OUT}")"

echo "Running script..."
echo "Sample Key Column: ${SAMPLE_KEY}"

# < [CRITICAL FIX] Completed the Rscript command with all necessary flags
Rscript "${SCRIPT}" \
  --obj "${OBJ}" \
  --objtype rds \
  --up "${UP}" \
  --out "${OUT}" \
  --assay_in peaks \
  --assay_out peaks_universal \
  --mode per_sample \
  --sample_key "${SAMPLE_KEY}" \
  --frag_tpl "${FRAG_TPL}"

echo "Rscript finished. Running quick check..."

# < Quick check to ensure the matrix is healthy and NOT empty
Rscript -e '
x <- readRDS("'"${OUT}"'");
cat("DefaultAssay:", DefaultAssay(x), "\n");
mat <- GetAssayData(x, assay="peaks_universal", layer="counts");
cat("Universal Matrix Dimensions:", nrow(mat), "rows x", ncol(mat), "cols\n");
cat("Total Reads in Matrix:", sum(mat, na.rm=TRUE), "\n");
if(sum(mat, na.rm=TRUE) == 0) {
  stop("CRITICAL ERROR: Matrix is completely empty (0 reads)!")
} else {
  cat("SUCCESS: Data successfully extracted!\n")
}
'
