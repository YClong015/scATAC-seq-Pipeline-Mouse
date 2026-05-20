#!/bin/bash
# ============================================================
# Lung COPD scATAC-seq — Zhang Q. et al. 2025 PLOS ONE
# Source: CNGBdb project CNP0004399
#   https://db.cngb.org/search/project/CNP0004399/
# 6 BGI samples, naming:
#   Control_F2 -> CL100168054_L01
#   Control_M1 -> CL100167942_L01
#   Case_F1    -> CL100168078_L02
#   Case_F3    -> CL100168054_L02
#   Case_M2    -> CL100167942_L02
#   Case_M3    -> CL100168078_L01
# ============================================================

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/QRISdata/Q8448/Mouse_disease_data}"
OUTDIR="${DATA_ROOT}/Lung/raw_fastq"
mkdir -p "${OUTDIR}"

echo "=================================================="
echo "Lung COPD scATAC-seq download (CNGBdb CNP0004399)"
echo "Output: ${OUTDIR}"
echo "=================================================="
echo ""
echo "CNGB does not provide a stable wget URL pattern for these"
echo "raw FASTQ; you must download manually using the CNGBdb"
echo "data delivery tool (cngbdb-cli) or the project web interface:"
echo ""
echo "    https://db.cngb.org/search/project/CNP0004399/"
echo ""
echo "After download, the directory structure should be:"
echo "  ${DATA_ROOT}/Lung/Lung_cellatac/{CL100168054_L01,CL100167942_L01,...}/outs/"
echo "    fragments.tsv.gz"
echo "    singlecell.csv"
echo "    filtered_peak_bc_matrix.h5"
echo "    peaks.bed"
echo ""
echo "These are the standard Cell Ranger ATAC outs/. The CNGB project"
echo "may provide already-processed Cell Ranger outputs alongside FASTQ."
echo ""
echo "Sample → BGI ID mapping (used throughout the pipeline):"
cat <<'EOF'
  Control_F2 -> CL100168054_L01
  Control_M1 -> CL100167942_L01
  Case_F1    -> CL100168078_L02
  Case_F3    -> CL100168054_L02
  Case_M2    -> CL100167942_L02
  Case_M3    -> CL100168078_L01
EOF
