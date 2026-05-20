#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=build_univ_scpm
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

echo "$(date)  Starting build_universal_scpm.py (v5)"

PYTHON="/home/s4869245/.conda/envs/scenicplus/bin/python"

# ── Environment variables ──────────────────────────────────────────────────
export BLACKLIST="/scratch/user/s4869245/pycisTopic/blacklist/mm10-blacklist.v2.bed"
export OUTDIR="/QRISdata/Q8448/Mouse_disease_data/universal_peaks_v5"
export CHROMSIZES_FILE="/scratch/user/s4869245/pycisTopic/chromsizes/mm10.chrom.sizes"
export PEAK_HALF_WIDTH="250"
export SCPM_CUTOFF="1.0"

# ── Run ───────────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"

${PYTHON} /home/s4869245/scripts/universal_peaks/build_universal_scpm.py

echo "$(date)  Done."
