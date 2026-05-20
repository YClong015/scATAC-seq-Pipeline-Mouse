#!/usr/bin/env bash
set -euo pipefail

export HOMER_HOME="/scratch/user/s4869245/homer"
export PATH="${HOMER_HOME}/bin:${PATH}"
hash -r

GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-10}"
MIN_PEAKS=20

BASE_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2"
BED_DIR="${BASE_DIR}/DAR_BED_stable"
OUT_DIR="${BASE_DIR}/HOMER_stable_bg"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

for stable in "${BED_DIR}"/*__stable.bed; do
  [ -f "${stable}" ] || continue

  base="$(basename "${stable}")"
  prefix="${base%__stable.bed}"

  opening="${BED_DIR}/${prefix}__opening.bed"
  closing="${BED_DIR}/${prefix}__closing.bed"

  n_stable="$(wc -l < "${stable}" | tr -d ' ')"
  n_opening=0
  n_closing=0

  [ -f "${opening}" ] && n_opening="$(wc -l < "${opening}" | tr -d ' ')"
  [ -f "${closing}" ] && n_closing="$(wc -l < "${closing}" | tr -d ' ')"

  echo "======================================================"
  echo "Comparison: ${prefix}"
  echo "opening=${n_opening} closing=${n_closing} stable=${n_stable}"

  # -------------------------
  # 1. opening vs stable
  # -------------------------
  if [ -f "${opening}" ] && [ "${n_opening}" -ge "${MIN_PEAKS}" ] && \
     [ "${n_stable}" -ge "${MIN_PEAKS}" ]; then

    out_sub="${OUT_DIR}/${prefix}__opening_vs_stable"
    log_file="${LOG_DIR}/${prefix}__opening_vs_stable.log"

    if [ ! -s "${out_sub}/knownResults.txt" ]; then
      rm -rf "${out_sub}"
      echo "RUN opening_vs_stable: ${prefix}"
      findMotifsGenome.pl "${opening}" "${GENOME}" "${out_sub}" \
        -bg "${stable}" \
        -size given \
        -p "${CPU}" \
        -nomotif \
        2>&1 | tee "${log_file}"
    else
      echo "Done already: ${prefix}__opening_vs_stable"
    fi
  else
    echo "Skip opening_vs_stable: ${prefix}"
  fi

  # -------------------------
  # 2. closing vs stable
  # -------------------------
  if [ -f "${closing}" ] && [ "${n_closing}" -ge "${MIN_PEAKS}" ] && \
     [ "${n_stable}" -ge "${MIN_PEAKS}" ]; then

    out_sub="${OUT_DIR}/${prefix}__closing_vs_stable"
    log_file="${LOG_DIR}/${prefix}__closing_vs_stable.log"

    if [ ! -s "${out_sub}/knownResults.txt" ]; then
      rm -rf "${out_sub}"
      echo "RUN closing_vs_stable: ${prefix}"
      findMotifsGenome.pl "${closing}" "${GENOME}" "${out_sub}" \
        -bg "${stable}" \
        -size given \
        -p "${CPU}" \
        -nomotif \
        2>&1 | tee "${log_file}"
    else
      echo "Done already: ${prefix}__closing_vs_stable"
    fi
  else
    echo "Skip closing_vs_stable: ${prefix}"
  fi

  # -------------------------
  # 3. stable vs opening
  # -------------------------
  if [ -f "${opening}" ] && [ "${n_opening}" -ge "${MIN_PEAKS}" ] && \
     [ "${n_stable}" -ge "${MIN_PEAKS}" ]; then

    out_sub="${OUT_DIR}/${prefix}__stable_vs_opening"
    log_file="${LOG_DIR}/${prefix}__stable_vs_opening.log"

    if [ ! -s "${out_sub}/knownResults.txt" ]; then
      rm -rf "${out_sub}"
      echo "RUN stable_vs_opening: ${prefix}"
      findMotifsGenome.pl "${stable}" "${GENOME}" "${out_sub}" \
        -bg "${opening}" \
        -size given \
        -p "${CPU}" \
        -nomotif \
        2>&1 | tee "${log_file}"
    else
      echo "Done already: ${prefix}__stable_vs_opening"
    fi
  else
    echo "Skip stable_vs_opening: ${prefix}"
  fi

  # -------------------------
  # 4. stable vs closing
  # -------------------------
  if [ -f "${closing}" ] && [ "${n_closing}" -ge "${MIN_PEAKS}" ] && \
     [ "${n_stable}" -ge "${MIN_PEAKS}" ]; then

    out_sub="${OUT_DIR}/${prefix}__stable_vs_closing"
    log_file="${LOG_DIR}/${prefix}__stable_vs_closing.log"

    if [ ! -s "${out_sub}/knownResults.txt" ]; then
      rm -rf "${out_sub}"
      echo "RUN stable_vs_closing: ${prefix}"
      findMotifsGenome.pl "${stable}" "${GENOME}" "${out_sub}" \
        -bg "${closing}" \
        -size given \
        -p "${CPU}" \
        -nomotif \
        2>&1 | tee "${log_file}"
    else
      echo "Done already: ${prefix}__stable_vs_closing"
    fi
  else
    echo "Skip stable_vs_closing: ${prefix}"
  fi
done

echo "All HOMER jobs finished."
echo "Output dir: ${OUT_DIR}"
