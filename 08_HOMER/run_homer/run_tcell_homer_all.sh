#!/usr/bin/env bash
set -euo pipefail

export HOMER_HOME="/scratch/user/s4869245/homer"
export PATH="${HOMER_HOME}/bin:${PATH}"
hash -r

GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-10}"
MIN_PEAKS=20

BASE_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Tcells_DESeq2"
DAR_BED_DIR="${BASE_DIR}/DAR_BED"
DAR_TAB_DIR="${BASE_DIR}/DAR_tables"

OUT_BASE="${BASE_DIR}/HOMER"
BG_DIR="${OUT_BASE}/background_bed"
LOG_DIR="${OUT_BASE}/logs"

mkdir -p "${OUT_BASE}" "${BG_DIR}" "${LOG_DIR}"

beds=(
  "${DAR_BED_DIR}"/*__opening.bed
  "${DAR_BED_DIR}"/*__closing.bed
)

for bed in "${beds[@]}"; do
  [ -f "${bed}" ] || continue

  bn="$(basename "${bed}")"
  stub="${bn%.bed}"

  if [[ "${stub}" == *"__opening" ]]; then
    direction="opening"
    prefix="${stub%__opening}"
  elif [[ "${stub}" == *"__closing" ]]; then
    direction="closing"
    prefix="${stub%__closing}"
  else
    echo "Skip: ${bn}"
    continue
  fi

  n_peaks="$(wc -l < "${bed}" | tr -d ' ')"
  if [ "${n_peaks}" -lt "${MIN_PEAKS}" ]; then
    echo "Skip (<${MIN_PEAKS} peaks): ${bn} (${n_peaks})"
    continue
  fi

  all_match=( "${DAR_TAB_DIR}/${prefix}"__*_DESeq2_all.tsv )
  if [ "${#all_match[@]}" -eq 0 ]; then
    echo "WARN: missing all.tsv for ${prefix}"
    continue
  fi

  all_tsv="${all_match[0]}"
  if [ ! -s "${all_tsv}" ]; then
    echo "WARN: empty all.tsv for ${prefix}"
    continue
  fi

  bg_bed="${BG_DIR}/${prefix}__background_tested.bed"
  if [ ! -s "${bg_bed}" ]; then
    awk -F'\t' '
BEGIN {OFS="\t"}
NR==1 {
  peak_col=0
  for (i=1; i<=NF; i++) {
    if ($i=="peak") {
      peak_col=i
      break
    }
  }
  if (peak_col==0) {
    print "ERROR: peak column not found in " FILENAME > "/dev/stderr"
    exit 1
  }
  next
}
{
  split($peak_col, a, "-")
  if (length(a) >= 3) {
    print a[1], a[2], a[3]
  }
}
' "${all_tsv}" | sort -k1,1 -k2,2n > "${bg_bed}"
  fi

  out_dir="${OUT_BASE}/${prefix}__${direction}"
  log_file="${LOG_DIR}/${prefix}__${direction}.log"

  if [ -d "${out_dir}" ] && [ -s "${out_dir}/knownResults.txt" ]; then
    echo "Done already: ${prefix} ${direction}"
    continue
  fi

  echo "RUN: ${prefix} ${direction} | peaks=${n_peaks}"
  echo "  bed: ${bed}"
  echo "  bg : ${bg_bed}"
  echo "  out: ${out_dir}"

  rm -rf "${out_dir}"

  findMotifsGenome.pl "${bed}" "${GENOME}" "${out_dir}" \
    -bg "${bg_bed}" \
    -size given \
    -p "${CPU}" \
    -nomotif \
    2>&1 | tee "${log_file}"
done

echo "All HOMER jobs finished."
echo "Output dir: ${OUT_BASE}"
