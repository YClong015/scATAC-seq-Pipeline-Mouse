#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=PeakMerge_Universal
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load anaconda3/2023.09-0
source "$EBROOTANACONDA3/etc/profile.d/conda.sh"
conda activate scenicplus

BASE="/QRISdata/Q8448/Mouse_disease_data"
UROOT="${BASE}/Universal_peaks/Kidney_Aorta_Tcells_Lung"
OUTDIR="${UROOT}/consensus_peak_calling"
mkdir -p "${OUTDIR}"

export PEAK_HALF_WIDTH="250"
export BLACKLIST="/scratch/user/s4869245/pycisTopic/blacklist/mm10-blacklist.v2.bed"
export CHROMSIZES_FILE="/scratch/user/s4869245/mouse_mm10_ref/chrom.sizes"

PEAKPATH="${UROOT}/tissue_consensus_paths.txt"
: > "${PEAKPATH}"

printf "%s: %s\n" \
  "Kidney" \
  "${BASE}/Kidney/peak_merging/consensus_peak_calling/consensus_regions.bed.gz" \
  >> "${PEAKPATH}"

printf "%s: %s\n" \
  "Aorta" \
  "${BASE}/Aorta/peak_merging/consensus_peak_calling/consensus_regions.bed.gz" \
  >> "${PEAKPATH}"

printf "%s: %s\n" \
  "Tcells" \
  "${BASE}/Tcells/peak_merging/consensus_peak_calling/consensus_regions.bed.gz" \
  >> "${PEAKPATH}"

printf "%s: %s\n" \
  "Lung" \
  "${BASE}/Lung/peak_merging/consensus_peak_calling/consensus_regions.bed.gz" \
  >> "${PEAKPATH}"

echo "[INFO] PEAKPATH:"
cat "${PEAKPATH}"

export PEAKPATH
export OUTDIR

PY="$HOME/.conda/envs/scenicplus/bin/python"
SCRIPT="${SLURM_SUBMIT_DIR}/merge_universal_mm10.py"
"$PY" "$SCRIPT"

echo -n "[INFO] Lines: "
zcat "${OUTDIR}/consensus_regions.bed.gz" | wc -l

conda deactivate
echo "Universal peaks done."

