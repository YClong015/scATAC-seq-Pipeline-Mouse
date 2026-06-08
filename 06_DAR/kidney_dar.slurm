#!/bin/bash
#SBATCH --job-name=Kidney_DAR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=Kidney_DAR_%j.out
#SBATCH --error=Kidney_DAR_%j.err

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"

module load r/4.4.2

Rscript Kidney_pseudo-bulk_DAR.R

echo "Job finished at: $(date)"
