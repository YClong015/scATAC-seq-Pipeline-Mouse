#!/bin/bash
#SBATCH --job-name=homer_tcells
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=homer_tcells_%j.out
#SBATCH --error=homer_tcells_%j.err

./run_homer_all.sh
