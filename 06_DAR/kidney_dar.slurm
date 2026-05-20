#!/bin/bash
#SBATCH --job-name=Kidney_DAR          # Job name
#SBATCH --nodes=1                      # Request 1 node
#SBATCH --ntasks=1                     # 1 task
#SBATCH --cpus-per-task=16             # Request 16 CPU cores (must match workers = 16 in the R code)
#SBATCH --mem=256G                     # Request 150GB RAM (to avoid OOM when running 16 parallel workers)
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=Kidney_DAR_%j.out     # Stdout log file (%j will be replaced by the job ID)
#SBATCH --error=Kidney_DAR_%j.err      # Stderr log file

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"

module load r/4.4.2

Rscript kidney_dar.R

echo "Job finished at: $(date)"
