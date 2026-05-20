#!/bin/bash
#SBATCH --job-name=Kidney_ATAC_Merge      # Job name
#SBATCH --nodes=1                         # Number of nodes (R usually requires 1)
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=8                 # CPU cores per task (Recommended 4-8 for Seurat)
#SBATCH --mem=256G                        # Memory allocation (256GB - Essential for large merges)
#SBATCH --time=24:00:00                   # Time limit (hh:mm:ss)
#SBATCH --partition=general               # Partition name (e.g., general, highmem - Check your HPC)
#SBATCH --output=log_%x_%j.out            # Standard output log (%x=job_name, %j=job_id)
#SBATCH --error=log_%x_%j.err             # Standard error log
#SBATCH --account=a_imb_ccbcd

# ==========================================
# 1. Load Environment
# ==========================================
echo "Job started at: $(date)"

module purge                              # Clear existing modules to avoid conflicts

# Load R module
# NOTE: Check your HPC specific module name using 'module avail R'
# Common examples: R/4.3.0, r/4.2.2, foss/2022b
module load r

# Optional: If you use Conda/Mamba, uncomment the line below:
# source activate your_env_name

# ==========================================
# 2. Run R Script
# ==========================================
echo "Running Seurat/Signac Analysis..."

# IMPORTANT: Replace 'Your_Analysis_Script.R' with your actual R filename
Rscript Kidney_scATAC\(Combine\).R 

echo "Job finished at: $(date)"
