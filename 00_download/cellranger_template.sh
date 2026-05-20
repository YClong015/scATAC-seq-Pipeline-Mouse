#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --job-name=cellranger
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd

# ============================================================
# Template — run Cell Ranger ATAC on one sample.
# Submit one job per sample. Set SAMPLE, FASTQ_DIR, OUT_DIR.
# ============================================================

set -euo pipefail

# === EDIT THESE ===
SAMPLE="SRR27367347_Kidney_atac"
FASTQ_DIR="/QRISdata/Q8448/Mouse_disease_data/Kidney/SRR27367347_Fastq/SRR27367347"
OUT_DIR="/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/${SAMPLE}"
REF="/scratch/user/${USER}/cellranger_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"

# === Run ===
module load cellranger-atac/2.1.0  # adjust to your HPC's module
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

cellranger-atac count \
  --id="${SAMPLE}" \
  --reference="${REF}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE%%_*}" \
  --localcores=16 \
  --localmem=120

echo "Done: ${SAMPLE}"
echo "Outputs at ${OUT_DIR}/${SAMPLE}/outs/"
echo "  fragments.tsv.gz, singlecell.csv, filtered_peak_bc_matrix.h5, peaks.bed"
