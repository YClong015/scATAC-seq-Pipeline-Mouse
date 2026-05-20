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

# ============================================================
# Stage 05 — Lung (step 1 of 2):
# Re-quantify the cleaned + annotated Lung Seurat object
# against the OLD universal peak set.
# Output -> lung_universal.rds (will be pruned to
# lung_universal_new_pruned.rds by Diagnose_lung.R in step 2).
# ============================================================

set -euo pipefail
module load r/4.4.2

# OLD universal peak set (the one DAR was actually run against).
UP="/QRISdata/Q8448/Mouse_disease_data/Universal_peaks/Kidney_Aorta_Tcells_Lung/consensus_peak_calling/consensus_regions.bed.gz"

# Input: cleaned + annotated Lung Seurat object built by
# 01_preprocessing/lung/atac_Lung.Rmd.
OBJ="/QRISdata/Q8448/Mouse_disease_data/Lung/lung_integrated_clean_annotated.rds"

# Output BEFORE pruning. Diagnose_lung.R in this same dir then prunes
# the resulting object to match Kidney's feature count (667,459 peaks)
# and writes lung_universal_new_pruned.rds (which is the canonical DAR input).
OUT="/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal.rds"

FRAG_TPL="/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/{id}/outs/fragments.tsv.gz"
SAMPLE_KEY="folder_id"

if [[ ! -f "$OBJ" ]]; then
  echo "[ERROR] Object not found: $OBJ"; exit 1
fi
mkdir -p "$(dirname "${OUT}")"

Rscript SeuratObject.R \
  --obj "${OBJ}" \
  --objtype rds \
  --up "${UP}" \
  --out "${OUT}" \
  --assay_in peaks \
  --assay_out peaks_universal \
  --mode per_sample \
  --sample_key "${SAMPLE_KEY}" \
  --frag_tpl "${FRAG_TPL}"

echo "Step 1 done. Output: ${OUT}"
echo "Next: Rscript Diagnose_lung.R   # prunes to lung_universal_new_pruned.rds"
