#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=MACS
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail

module load anaconda3/2023.09-0
source "$EBROOTANACONDA3/etc/profile.d/conda.sh"
conda activate scenicplus

BED_DIR="/QRISdata/Q8448/Mouse_disease_data/Kidney"
BED_DIR="${BED_DIR}/fragment_files_split_by_celltype"
LIST="${BED_DIR}/DataList.txt"

DATA="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")"
DATA_BASE="${DATA%.bed}"

INPUT_BED="${BED_DIR}/${DATA}"
if [[ ! -f "${INPUT_BED}" ]]; then
  echo "ERROR: missing bed file: ${INPUT_BED}"
  exit 1
fi

OUT_DIR="/QRISdata/Q8448/Mouse_disease_data/Kidney"
OUT_DIR="${OUT_DIR}/peaks/${DATA_BASE}/sample_peak"
mkdir -p "${OUT_DIR}"

BEDPATH_FILE="${TMPDIR}/${DATA_BASE}_bedPaths.txt"
echo "${DATA_BASE}: ${INPUT_BED}" > "${BEDPATH_FILE}"

export OUTDIR="${OUT_DIR}"
export BEDPATH="${BEDPATH_FILE}"
export TMPDIR="${TMPDIR}"

PY="$HOME/.conda/envs/scenicplus/bin/python"
"$PY" /home/s4869245/scripts/peak_calling_for_kidney/peak_calling_kidney.py

conda deactivate

