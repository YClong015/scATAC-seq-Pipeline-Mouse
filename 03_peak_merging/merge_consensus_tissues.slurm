#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=PeakMerge_Tissues
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load anaconda3/2023.09-0
source "$EBROOTANACONDA3/etc/profile.d/conda.sh"
conda activate scenicplus

BASE="/QRISdata/Q8448/Mouse_disease_data"
SCRIPT="${SLURM_SUBMIT_DIR}/consensus_mm10.py"

export PEAK_HALF_WIDTH="250"
export BLACKLIST="/scratch/user/s4869245/pycisTopic/blacklist/mm10-blacklist.v2.bed"
export CHROMSIZES_FILE="/scratch/user/s4869245/mouse_mm10_ref/chrom.sizes"

PY="$HOME/.conda/envs/scenicplus/bin/python"

TISSUES=("Kidney" "Aorta" "Tcells" "Lung")

for TIS in "${TISSUES[@]}"; do
  echo "Tissue: ${TIS}"

  IN_ROOT="${BASE}/${TIS}/peaks"
  OUT_ROOT="${BASE}/${TIS}/peak_merging"
  OUT_DIR="${OUT_ROOT}/consensus_peak_calling"
  mkdir -p "${OUT_DIR}"

  PEAKPATH="${OUT_ROOT}/narrow_peak_paths.txt"
  : > "${PEAKPATH}"

  find "${IN_ROOT}" -type f \
    -path "*/sample_peak/*_peaks.narrowPeak" | sort \
    | while read -r p; do
        b="$(basename "${p}")"
        k="${b%_peaks.narrowPeak}"
        printf "%s: %s\n" "${k}" "${p}" >> "${PEAKPATH}"
      done

  echo "[INFO] Peak list: ${PEAKPATH}"
  echo -n "[INFO] #narrowPeak: "
  wc -l "${PEAKPATH}" | awk '{print $1}'

  export PEAKPATH
  export OUTDIR="${OUT_DIR}"

  "$PY" "$SCRIPT"

  echo -n "[INFO] Lines: "
  (zcat "${OUT_DIR}/consensus_regions.bed.gz" 2>/dev/null || \
   cat "${OUT_DIR}/consensus_regions.bed.gz") | wc -l
done

conda deactivate
echo "All tissue consensus peak merging done."

