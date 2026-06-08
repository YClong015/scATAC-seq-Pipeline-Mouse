#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=MACS_Tcells
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=1-9

set -euo pipefail

module load macs2/2.2.9.1
export MACS2=$(which macs2)

PY="/home/s4869245/.conda/envs/scenicplus/bin/python"

BED_DIR="/QRISdata/Q8448/Mouse_disease_data/Tcells/fragment_files_split_by_celltype_final"
LIST="${BED_DIR}/DataList.txt"

DATA="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")"
DATA_BASE="${DATA%.bed}"

INPUT_BED="${BED_DIR}/${DATA}"
if [[ ! -f "${INPUT_BED}" ]]; then
  echo "ERROR: missing fragment file: ${INPUT_BED}"
  exit 1
fi

OUT_ROOT="/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks"
FINAL_OUT_DIR="${OUT_ROOT}/${DATA_BASE}/sample_peak"
mkdir -p "${FINAL_OUT_DIR}"

BEDPATH_FILE="${TMPDIR}/${DATA_BASE}_bedPaths.txt"
echo "${DATA_BASE}: ${INPUT_BED}" > "${BEDPATH_FILE}"

export OUTDIR="${FINAL_OUT_DIR}"
export BEDPATH="${BEDPATH_FILE}"
export TMPDIR="${TMPDIR}"

"${PY}" ${SLURM_SUBMIT_DIR}/peak_calling_Tcell.py
