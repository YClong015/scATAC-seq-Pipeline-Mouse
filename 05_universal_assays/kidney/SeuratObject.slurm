#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --job-name=KidneyUni
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
module load r/4.4.2

UP="/QRISdata/Q8448/Mouse_disease_data/Universal_peaks/Kidney_Aorta_Tcells_Lung/consensus_peak_calling/consensus_regions.bed.gz"

OBJ="/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/kidney_merged_annotated.rds"

OUT="/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds"

FRAG_TPL="/QRISdata/Q8448/Mouse_disease_data/Kidney/\
cellranger_unpacked_data/{id}_Kidney_atac/{id}_Kidney_atac/outs/fragments.tsv.gz"

Rscript SeuratObject.R \
  --obj "$OBJ" \
  --objtype rds \
  --up "$UP" \
  --out "$OUT" \
  --assay_in peaks \
  --assay_out peaks_universal \
  --mode per_sample \
  --sample_key dataset \
  --frag_tpl "$FRAG_TPL"

