#!/usr/bin/env bash
set -euo pipefail

export HOMER_HOME="/scratch/user/s4869245/homer"
export PATH="${HOMER_HOME}/bin:${PATH}"
hash -r

GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-10}"
MIN_PEAKS=20

BASE_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2"
BED_DIR="${BASE_DIR}/DAR_BED_NS"
OUT_DIR="${BASE_DIR}/HOMER_NS_bg"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

run_homer_if_needed () {
  local fg="$1"
  local bg="$2"
  local out_sub="$3"
  local log_file="$4"

  if [ ! -s "${fg}" ] || [ ! -s "${bg}" ]; then
    echo "Skip (missing/empty file): ${out_sub}"
    return
  fi

  local n_fg
  local n_bg
  n_fg="$(wc -l < "${fg}" | tr -d ' ')"
  n_bg="$(wc -l < "${bg}" | tr -d ' ')"

  if [ "${n_fg}" -lt "${MIN_PEAKS}" ] || [ "${n_bg}" -lt "${MIN_PEAKS}" ]; then
    echo "Skip (<${MIN_PEAKS} peaks): ${out_sub} | fg=${n_fg}, bg=${n_bg}"
    return
  fi

  if [ -s "${out_sub}/knownResults.txt" ]; then
    echo "Done already: $(basename "${out_sub}")"
    return
  fi

  rm -rf "${out_sub}"

  echo "RUN: $(basename "${out_sub}")"
  echo "  FG: ${fg}"
  echo "  BG: ${bg}"

  findMotifsGenome.pl "${fg}" "${GENOME}" "${out_sub}" \
    -bg "${bg}" \
    -size given \
    -p "${CPU}" \
    -nomotif \
    2>&1 | tee "${log_file}"
}

for ns in "${BED_DIR}"/*Case_vs_Control*__NS.bed; do
  [ -f "${ns}" ] || continue

  base="$(basename "${ns}")"
  prefix="${base%__NS.bed}"

  opening="${BED_DIR}/${prefix}__opening.bed"
  closing="${BED_DIR}/${prefix}__closing.bed"

  echo "======================================================"
  echo "Processing: ${prefix}"

  run_homer_if_needed \
    "${opening}" \
    "${ns}" \
    "${OUT_DIR}/${prefix}__opening_vs_NS" \
    "${LOG_DIR}/${prefix}__opening_vs_NS.log"

  run_homer_if_needed \
    "${closing}" \
    "${ns}" \
    "${OUT_DIR}/${prefix}__closing_vs_NS" \
    "${LOG_DIR}/${prefix}__closing_vs_NS.log"

  run_homer_if_needed \
    "${ns}" \
    "${opening}" \
    "${OUT_DIR}/${prefix}__NS_vs_opening" \
    "${LOG_DIR}/${prefix}__NS_vs_opening.log"

  run_homer_if_needed \
    "${ns}" \
    "${closing}" \
    "${OUT_DIR}/${prefix}__NS_vs_closing" \
    "${LOG_DIR}/${prefix}__NS_vs_closing.log"
done

echo "Lung HOMER NS resume finished."
echo "Output dir: ${OUT_DIR}"
