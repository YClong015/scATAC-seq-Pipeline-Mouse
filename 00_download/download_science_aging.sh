#!/bin/bash
# ============================================================
# Mouse aging chromatin atlas — Lu et al. 2026 Science
# DOI: 10.1126/science.adw6273
# Portal: https://epiage.net/
#
# We use the two GSM-level peak-count h5ad files matching our
# disease tissues:
#   - GSM8774007_Kidney_peak_count.h5ad
#   - GSM8774006_Lung_peak_count.h5ad
# ============================================================

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/QRISdata/Q8448/Mouse_disease_data}"
SCIENCE_DIR="${DATA_ROOT}/DAR/DAR_science_comparison"
mkdir -p "${SCIENCE_DIR}/kidney_processed" "${SCIENCE_DIR}/lung_processed"

echo "=================================================="
echo "Lu et al. 2026 Science aging chromatin atlas"
echo "Output: ${SCIENCE_DIR}"
echo "=================================================="
echo ""
echo "Download the peak-count h5ad files from GEO:"
echo "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283134"
echo ""
echo "Or from the epiage.net portal:"
echo "  https://epiage.net/  →  Downloads  →  ATAC peak counts"
echo ""
echo "Place files at:"
echo "  ${SCIENCE_DIR}/kidney_processed/GSM8774007_Kidney_peak_count.h5ad"
echo "  ${SCIENCE_DIR}/lung_processed/GSM8774006_Lung_peak_count.h5ad"
echo ""
echo "After download, inspect the file structure with:"
echo "  bash ../09_aging_comparison/01_explore_h5ad.sh"
echo ""
echo "Each h5ad provides cells × peaks count matrix + per-cell"
echo "Main_cell_type, Age (Young/Adult/Aged), Sample, Gender."
