#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=MACS_Lung
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=1-15

set -euo pipefail

module load anaconda3/2023.09-0
source "$EBROOTANACONDA3/etc/profile.d/conda.sh"
conda activate scenicplus

BED_DIR="/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_fragments_file"
LIST="${BED_DIR}/DataList.txt"

DATA="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")"
DATA_BASE="${DATA%.bed}"

INPUT_BED="${BED_DIR}/${DATA}"
if [[ ! -f "${INPUT_BED}" ]]; then
  echo "ERROR: missing bed file: ${INPUT_BED}"
  exit 1
fi

OUT_ROOT="/QRISdata/Q8448/Mouse_disease_data/Lung/peaks"
FINAL_OUT_DIR="${OUT_ROOT}/${DATA_BASE}/sample_peak"
mkdir -p "${FINAL_OUT_DIR}"

BEDPATH_FILE="${TMPDIR}/${DATA_BASE}_bedPaths.txt"
echo "${DATA_BASE}: ${INPUT_BED}" > "${BEDPATH_FILE}"

export OUTDIR="${FINAL_OUT_DIR}"
export BEDPATH="${BEDPATH_FILE}"
export TMPDIR="${TMPDIR}"

PY="$HOME/.conda/envs/scenicplus/bin/python"
"$PY" /home/s4869245/scripts/Lung/peak_calling_for_lung/peak_calling.py

conda deactivate


