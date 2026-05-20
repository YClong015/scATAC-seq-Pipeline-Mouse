#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --job-name=TcellsUni
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
module load r/4.4.2

# 1. Universal Peak file (v5)
UP="/QRISdata/Q8448/Mouse_disease_data/Universal_peaks/Kidney_Aorta_Tcells_Lung/consensus_peak_calling/consensus_regions.bed.gz"

# 2. Input
OBJ="/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_processed.rds"

# 3. Output
OUT="/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
mkdir -p "$(dirname "${OUT}")"

# 4. Fragment file (the only one)
FRAG_FILE="/QRISdata/Q8448/Mouse_disease_data/Tcells/atac_fragments.tsv.gz"

# 5. Script path
SCRIPT="./SeuratObject.R"

echo "Processing Tcells using deMultliplex2_final_mapped..."


Rscript "${SCRIPT}" \
  --obj "${OBJ}" \
  --objtype rds \
  --up "${UP}" \
  --out "${OUT}" \
  --assay_in ATAC \
  --assay_out peaks_universal \
  --mode per_sample \
  --sample_key deMultliplex2_final_mapped \
  --frag_tpl "${FRAG_FILE}"

# Validate the output
Rscript -e '
x <- readRDS("'"${OUT}"'");
cat("----------------------------------\n");
cat("Success! Tcells universal assay added.\n");
cat("Features: ", nrow(x[["peaks_universal"]]), "\n");
cat("Assay Key used: deMultliplex2_final_mapped\n");
cat("----------------------------------\n");
'
