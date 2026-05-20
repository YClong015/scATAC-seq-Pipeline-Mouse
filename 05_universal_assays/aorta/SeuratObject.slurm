#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --job-name=AortaUniAssay
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load r/4.4.2

# Universal Peaks BED file
UP="/QRISdata/Q8448/Mouse_disease_data/universal_peaks_v5/consensus_regions_v5.bed.gz"

# Input Seurat Object
OBJ="/QRISdata/Q8448/Mouse_disease_data/Aorta/aorta_filtered_obj.rds"

# Output Seurat Object
OUT="/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_universal_v5.rds"

# Fragments template
FRAG_TPL="/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac/{id}/outs/fragments.tsv.gz"

# Path to the NEW R script
SCRIPT="/home/s4869245/scripts/integration_analysis/Aortic_seurat_obj/SeuratObject.R"

echo "Starting Rscript..."

# < [CRITICAL FIX] Removed the trailing backslash, and made sure --frag_tpl is passed!
Rscript "${SCRIPT}" \
  --obj "${OBJ}" \
  --objtype rds \
  --up "${UP}" \
  --out "${OUT}" \
  --assay_in peaks \
  --assay_out peaks_universal \
  --mode per_sample \
  --sample_key sample_id \
  --frag_tpl "${FRAG_TPL}"

echo "Rscript finished. Running quick check..."

# Quick check to ensure the object is healthy and NOT empty
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
