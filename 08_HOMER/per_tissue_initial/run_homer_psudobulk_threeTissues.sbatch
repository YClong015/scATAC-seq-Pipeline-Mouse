#!/bin/bash
#SBATCH --job-name=homer_pseudobulk
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=homer_pseudobulk_%j.out
#SBATCH --error=homer_pseudobulk_%j.err

set -euo pipefail

HOMER_HOME="/scratch/user/s4869245/homer"
export PATH="${PATH}:${HOMER_HOME}/bin"

BASE="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2"
GENOME="mm10"
CPU=10
MIN_PEAKS=20

TISSUES=(
  "DAR_pseudobulk_Aorta_DESeq2"
  "DAR_pseudobulk_Lung_DESeq2"
  "DAR_pseudobulk_Tcells_DESeq2"
)

echo "========================================"
echo "HOMER batch run (pseudo-bulk DAR BED)"
echo "BASE      : ${BASE}"
echo "GENOME    : ${GENOME}"
echo "CPU       : ${CPU}"
echo "MIN_PEAKS : ${MIN_PEAKS}"
echo "HOMER_HOME: ${HOMER_HOME}"
echo "========================================"

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then
  echo "ERROR: findMotifsGenome.pl not found in PATH"
  echo "PATH=${PATH}"
  exit 1
fi

for T in "${TISSUES[@]}"; do
  TDIR="${BASE}/${T}"
  BED_DIR="${TDIR}/DAR_BED"
  OUT_DIR="${TDIR}/HOMER"

  echo ""
  echo "----------------------------------------"
  echo "Tissue folder: ${TDIR}"
  echo "BED dir      : ${BED_DIR}"
  echo "HOMER out    : ${OUT_DIR}"
  echo "----------------------------------------"

  if [ ! -d "${BED_DIR}" ]; then
    echo "WARNING: BED dir not found, skip: ${BED_DIR}"
    continue
  fi

  mkdir -p "${OUT_DIR}"

  shopt -s nullglob
  beds=( "${BED_DIR}"/*.bed )

  if [ ${#beds[@]} -eq 0 ]; then
    echo "WARNING: No BED files found in ${BED_DIR}"
    continue
  fi

  for bed in "${beds[@]}"; do
    base_name="$(basename "${bed}" .bed)"
    out_sub="${OUT_DIR}/${base_name}"

    n_peaks="$(awk 'NF>0 {n++} END {print n+0}' "${bed}")"
    echo "BED: ${base_name}  peaks=${n_peaks}"

    if [ "${n_peaks}" -lt "${MIN_PEAKS}" ]; then
      echo "  -> Skip (<${MIN_PEAKS} peaks)"
      continue
    fi

    mkdir -p "${out_sub}"

    findMotifsGenome.pl \
      "${bed}" \
      "${GENOME}" \
      "${out_sub}" \
      -size given \
      -mask \
      -nomotif \
      -p "${CPU}" \
      > "${out_sub}/homer.stdout.log" \
      2> "${out_sub}/homer.stderr.log" || true

    if [ -f "${out_sub}/knownResults.txt" ]; then
      echo "  -> OK: knownResults.txt generated"
    else
      echo "  -> WARNING: knownResults.txt missing (check logs)"
    fi
  done
done

echo ""
echo "========================================"
echo "HOMER batch finished."
echo "========================================"
