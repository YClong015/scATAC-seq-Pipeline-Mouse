#!/usr/bin/env bash
set -euo pipefail

export HOMER_HOME="/scratch/user/s4869245/homer"
export PATH="${HOMER_HOME}/bin:${PATH}"
hash -r

GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-10}"
MIN_PEAKS=20

BASES=(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2"
)

run_one_mode () {
  local base_dir="$1"
  local mode="$2"

  local bed_dir=""
  local out_dir=""
  local bg_name=""

  if [ "${mode}" = "stable" ]; then
    bed_dir="${base_dir}/DAR_BED_stable"
    out_dir="${base_dir}/HOMER_stable_bg"
    bg_name="stable"
  elif [ "${mode}" = "NS" ]; then
    bed_dir="${base_dir}/DAR_BED_NS"
    out_dir="${base_dir}/HOMER_NS_bg"
    bg_name="NS"
  else
    echo "Unknown mode: ${mode}"
    exit 1
  fi

  local log_dir="${out_dir}/logs"
  mkdir -p "${out_dir}" "${log_dir}"

  for bg in "${bed_dir}"/*__${bg_name}.bed; do
    [ -f "${bg}" ] || continue

    local base
    local prefix
    local opening
    local closing
    local n_bg
    local n_opening
    local n_closing

    base="$(basename "${bg}")"
    prefix="${base%__${bg_name}.bed}"

    opening="${bed_dir}/${prefix}__opening.bed"
    closing="${bed_dir}/${prefix}__closing.bed"

    n_bg="$(wc -l < "${bg}" | tr -d ' ')"
    n_opening=0
    n_closing=0

    [ -f "${opening}" ] && n_opening="$(wc -l < "${opening}" | tr -d ' ')"
    [ -f "${closing}" ] && n_closing="$(wc -l < "${closing}" | tr -d ' ')"

    echo "======================================================"
    echo "BASE: ${base_dir}"
    echo "MODE: ${mode}"
    echo "COMP: ${prefix}"
    echo "opening=${n_opening} closing=${n_closing} ${bg_name}=${n_bg}"

    # opening vs background
    if [ -f "${opening}" ] && [ "${n_opening}" -ge "${MIN_PEAKS}" ] && \
       [ "${n_bg}" -ge "${MIN_PEAKS}" ]; then
      out_sub="${out_dir}/${prefix}__opening_vs_${bg_name}"
      log_file="${log_dir}/${prefix}__opening_vs_${bg_name}.log"
      if [ ! -s "${out_sub}/knownResults.txt" ]; then
        rm -rf "${out_sub}"
        findMotifsGenome.pl "${opening}" "${GENOME}" "${out_sub}" \
          -bg "${bg}" -size given -p "${CPU}" -nomotif \
          2>&1 | tee "${log_file}"
      fi
    fi

    # closing vs background
    if [ -f "${closing}" ] && [ "${n_closing}" -ge "${MIN_PEAKS}" ] && \
       [ "${n_bg}" -ge "${MIN_PEAKS}" ]; then
      out_sub="${out_dir}/${prefix}__closing_vs_${bg_name}"
      log_file="${log_dir}/${prefix}__closing_vs_${bg_name}.log"
      if [ ! -s "${out_sub}/knownResults.txt" ]; then
        rm -rf "${out_sub}"
        findMotifsGenome.pl "${closing}" "${GENOME}" "${out_sub}" \
          -bg "${bg}" -size given -p "${CPU}" -nomotif \
          2>&1 | tee "${log_file}"
      fi
    fi

    # background vs opening
    if [ -f "${opening}" ] && [ "${n_opening}" -ge "${MIN_PEAKS}" ] && \
       [ "${n_bg}" -ge "${MIN_PEAKS}" ]; then
      out_sub="${out_dir}/${prefix}__${bg_name}_vs_opening"
      log_file="${log_dir}/${prefix}__${bg_name}_vs_opening.log"
      if [ ! -s "${out_sub}/knownResults.txt" ]; then
        rm -rf "${out_sub}"
        findMotifsGenome.pl "${bg}" "${GENOME}" "${out_sub}" \
          -bg "${opening}" -size given -p "${CPU}" -nomotif \
          2>&1 | tee "${log_file}"
      fi
    fi

    # background vs closing
    if [ -f "${closing}" ] && [ "${n_closing}" -ge "${MIN_PEAKS}" ] && \
       [ "${n_bg}" -ge "${MIN_PEAKS}" ]; then
      out_sub="${out_dir}/${prefix}__${bg_name}_vs_closing"
      log_file="${log_dir}/${prefix}__${bg_name}_vs_closing.log"
      if [ ! -s "${out_sub}/knownResults.txt" ]; then
        rm -rf "${out_sub}"
        findMotifsGenome.pl "${bg}" "${GENOME}" "${out_sub}" \
          -bg "${closing}" -size given -p "${CPU}" -nomotif \
          2>&1 | tee "${log_file}"
      fi
    fi
  done
}

for base in "${BASES[@]}"; do
  run_one_mode "${base}" "stable"
  run_one_mode "${base}" "NS"
done

echo "All HOMER jobs finished."
