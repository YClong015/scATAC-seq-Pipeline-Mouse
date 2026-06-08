#!/bin/bash
#SBATCH --job-name=Kidney_ATAC_Merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --output=log_%x_%j.out
#SBATCH --error=log_%x_%j.err
#SBATCH --account=a_imb_ccbcd

echo "Job started at: $(date)"
module purge
module load r/4.4.2

Rscript Kidney_scATAC_Combine.R

echo "Job finished at: $(date)"
