#!/bin/bash
#SBATCH --job-name=homer_kidney
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --output=homer_kidney_%j.out
#SBATCH --error=homer_kidney_%j.err

set -euo pipefail

# -----------------------------
# 0) User settings
# -----------------------------
HOMER_HOME="/scratch/user/s4869245/homer"
GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-10}"
MIN_PEAKS=20

BASE_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2"
DAR_BED_DIR="${BASE_DIR}/DAR_BED"
HOMER_OUT_DIR="${BASE_DIR}/HOMER"

# -----------------------------
# 1) Setup
# -----------------------------
export PATH="${PATH}:${HOMER_HOME}/bin"
mkdir -p "${HOMER_OUT_DIR}"

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then
  echo "ERROR: findMotifsGenome.pl not found in PATH"
  echo "PATH=${PATH}"
  exit 1
fi

if [ ! -d "${DAR_BED_DIR}" ]; then
  echo "ERROR: BED directory not found: ${DAR_BED_DIR}"
  exit 1
fi

# Deterministic BED order
mapfile -t bed_files < <(find "${DAR_BED_DIR}" -maxdepth 1 -type f -name "*.bed" | sort)

if [ "${#bed_files[@]}" -eq 0 ]; then
  echo "ERROR: No BED files found in ${DAR_BED_DIR}"
  exit 1
fi

echo "========================================"
echo "HOMER Kidney batch (no array)"
echo "BED count : ${#bed_files[@]}"
echo "GENOME    : ${GENOME}"
echo "CPU       : ${CPU}"
echo "MIN_PEAKS : ${MIN_PEAKS}"
echo "========================================"

ok_n=0
skip_n=0
fail_n=0

# -----------------------------
# 2) Loop all BED files
# -----------------------------
for bed in "${bed_files[@]}"; do
  base_name="$(basename "${bed}" .bed)"
  out_dir="${HOMER_OUT_DIR}/${base_name}"

  echo ""
  echo "----------------------------------------"
  echo "Processing: ${base_name}"
  echo "BED: ${bed}"
  echo "OUT: ${out_dir}"
  echo "----------------------------------------"

  # Count valid BED lines
  n_peaks="$(awk '!/^#/ && !/^track/ && !/^browser/ && NF>=3 {n++} END {print n+0}' "${bed}")"
  echo "Peaks: ${n_peaks}"

  if [ "${n_peaks}" -lt "${MIN_PEAKS}" ]; then
    echo "Skip: < ${MIN_PEAKS} peaks"
    skip_n=$((skip_n + 1))
    continue
  fi

  mkdir -p "${out_dir}"

  if findMotifsGenome.pl \
    "${bed}" \
    "${GENOME}" \
    "${out_dir}" \
    -size given \
    -mask \
    -nomotif \
    -p "${CPU}" \
    > "${out_dir}/homer.stdout.log" \
    2> "${out_dir}/homer.stderr.log"; then

    # Known motif result check (different HOMER versions may differ)
    if [ -s "${out_dir}/knownResults/knownResults.txt" ] || \
       [ -s "${out_dir}/knownResults.txt" ] || \
       [ -s "${out_dir}/knownResults.html" ] || \
       [ -s "${out_dir}/knownResults/knownResults.html" ]; then
      echo "Done: ${base_name}"
      ok_n=$((ok_n + 1))
    else
      echo "FAILED (no known motif result file): ${base_name}"
      fail_n=$((fail_n + 1))
    fi
  else
    echo "FAILED (command error): ${base_name}"
    fail_n=$((fail_n + 1))
  fi
done

echo ""
echo "========================================"
echo "Summary:"
echo "  Success : ${ok_n}"
echo "  Skipped : ${skip_n}"
echo "  Failed  : ${fail_n}"
echo "========================================"

if [ "${fail_n}" -gt 0 ]; then
  exit 1
fi

