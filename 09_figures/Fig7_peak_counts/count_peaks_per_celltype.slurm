#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=count_peaks
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Count cell-type-specific peaks per tissue (PRE-consensus)
#
# Counts lines in every {CellType}_peaks.narrowPeak file produced by MACS2 in
# the per-tissue peak-calling pipeline. These are the inputs to consensus_mm10.py
# (i.e., the peaks BEFORE per-tissue consensus merging).
#
# Outputs:
#   - peaks_per_celltype.csv          : tissue, cell_type, n_peaks (one row per cell type)
#   - peaks_per_celltype_summary.csv  : per-tissue totals + post-consensus count
#
# Usage:
#   sbatch count_peaks_per_celltype.slurm

set -euo pipefail

BASE="/QRISdata/Q8448/Mouse_disease_data"
# New output folder for the pooled-T-cell version; leaves the original
# peaks_per_celltype/ (old per-subtype version) completely untouched.
OUT_DIR="${BASE}/QC_figures/peaks_per_celltype_pooled"
mkdir -p "${OUT_DIR}"

RAW_CSV="${OUT_DIR}/peaks_per_celltype.csv"
SUM_CSV="${OUT_DIR}/peaks_per_celltype_summary.csv"

TISSUES=("Kidney" "Lung" "Aorta" "Tcells")

echo "tissue,cell_type,n_peaks" > "${RAW_CSV}"
echo "tissue,n_cell_types,total_peaks_precons,median_peaks,min_peaks,max_peaks,n_peaks_consensus" > "${SUM_CSV}"

for TIS in "${TISSUES[@]}"; do
  # T cells: subtype annotation lacks an independent reference, so peaks are
  # called from ALL T cells pooled into one group (peaks_pooled/). A single
  # pooled peak set is reported rather than per-subtype counts.
  # Other tissues use their per-cell-type peaks (peaks/).
  if [[ "${TIS}" == "Tcells" ]]; then
    PEAK_ROOT="${BASE}/${TIS}/peaks_pooled"
  else
    PEAK_ROOT="${BASE}/${TIS}/peaks"
  fi
  if [[ ! -d "${PEAK_ROOT}" ]]; then
    echo "  WARNING: directory missing for ${TIS}: ${PEAK_ROOT}" >&2
    continue
  fi

  echo "Processing ${TIS} ..."

  ## per cell type peak counts
  find "${PEAK_ROOT}" -type f -path "*/sample_peak/*_peaks.narrowPeak" | sort \
    | while read -r f; do
        ct=$(basename "${f}" "_peaks.narrowPeak")
        n=$(wc -l < "${f}")
        echo "${TIS},${ct},${n}" >> "${RAW_CSV}"
      done

  ## per tissue summary
  TISSUE_CONSENSUS="${BASE}/${TIS}/peak_merging/consensus_peak_calling/consensus_regions.bed.gz"
  if [[ -f "${TISSUE_CONSENSUS}" ]]; then
    n_cons=$(zcat "${TISSUE_CONSENSUS}" | wc -l)
  else
    n_cons="NA"
  fi

  # Aggregate stats from the rows just appended
  awk -F',' -v tis="${TIS}" -v ncons="${n_cons}" '
    BEGIN { n=0; tot=0; mn=0; mx=0 }
    $1==tis {
      n++; v=$3+0; vals[n]=v; tot+=v
      if (n==1) { mn=v; mx=v } else { if (v<mn) mn=v; if (v>mx) mx=v }
    }
    END {
      # median
      asort(vals)
      if (n==0) { print tis",0,0,0,0,0,"ncons; exit }
      if (n%2==1) med=vals[(n+1)/2]
      else        med=(vals[n/2]+vals[n/2+1])/2
      printf "%s,%d,%d,%g,%d,%d,%s\n", tis, n, tot, med, mn, mx, ncons
    }
  ' "${RAW_CSV}" >> "${SUM_CSV}"

done

echo ""
echo "=== Output ==="
echo "Raw counts:  ${RAW_CSV}"
echo "Summary:     ${SUM_CSV}"
echo ""
echo "Preview:"
cat "${SUM_CSV}" | column -t -s ','
