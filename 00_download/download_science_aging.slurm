#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=download_science_aging
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Mouse aging chromatin atlas, Lu et al. 2026 Science (GSE283134).
# Kidney/Lung peak-count h5ad files are downloaded manually from GEO or epiage.net.

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/QRISdata/Q8448/Mouse_disease_data}"
SCIENCE_DIR="${DATA_ROOT}/DAR/DAR_science_comparison"
mkdir -p "${SCIENCE_DIR}/kidney_processed" "${SCIENCE_DIR}/lung_processed"

echo "Download from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283134) or https://epiage.net/ and place at:"
echo "  ${SCIENCE_DIR}/kidney_processed/GSM8774007_Kidney_peak_count.h5ad"
echo "  ${SCIENCE_DIR}/lung_processed/GSM8774006_Lung_peak_count.h5ad"
