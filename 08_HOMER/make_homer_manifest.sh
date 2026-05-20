#!/bin/bash
set -euo pipefail

BED_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_Motif_Kidney/DAR_BED"
OUT_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_Motif_Kidney/HOMER"
MANIFEST="${OUT_DIR}/homer_bed_manifest.txt"
MIN_PEAKS=20

mkdir -p "${OUT_DIR}"
: > "${MANIFEST}"

echo "Scanning BED files in: ${BED_DIR}"
echo "Minimum peaks required: ${MIN_PEAKS}"
echo

count_total=0
count_keep=0
count_skip=0

shopt -s nullglob
for bed in "${BED_DIR}"/*.bed; do
  ((count_total+=1))
  n_peaks=$(wc -l < "${bed}" | tr -d ' ')
  base=$(basename "${bed}" .bed)

  if [[ "${n_peaks}" -lt "${MIN_PEAKS}" ]]; then
    echo "SKIP  ${base}  (${n_peaks} peaks)"
    ((count_skip+=1))
    continue
  fi

  echo "${bed}" >> "${MANIFEST}"
  echo "KEEP  ${base}  (${n_peaks} peaks)"
  ((count_keep+=1))
done

echo
echo "Done."
echo "Total BEDs : ${count_total}"
echo "Kept       : ${count_keep}"
echo "Skipped    : ${count_skip}"
echo "Manifest   : ${MANIFEST}"
