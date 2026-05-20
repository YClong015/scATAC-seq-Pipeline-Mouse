#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --job-name=AtlasMerge
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load r/4.4.2

echo "Starting Integration Pipeline..."
Rscript /home/s4869245/scripts/integration_analysis/Merge_and_integrate/Merge_and_integrate.R
