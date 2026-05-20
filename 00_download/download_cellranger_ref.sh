#!/bin/bash
# ============================================================
# Download the 10x Genomics Cell Ranger ATAC mm10 reference.
# Required for processing raw FASTQ into fragments.tsv.gz.
# ============================================================

set -euo pipefail

REF_ROOT="${REF_ROOT:-/scratch/user/$USER/cellranger_ref}"
mkdir -p "${REF_ROOT}"
cd "${REF_ROOT}"

REF_NAME="refdata-cellranger-arc-mm10-2020-A-2.0.0"

if [ ! -d "${REF_NAME}" ]; then
  echo "Downloading Cell Ranger ARC mm10 reference (~9 GB)..."
  wget "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"
  tar -xzf "refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"
  rm "refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"
else
  echo "Already present: ${REF_ROOT}/${REF_NAME}"
fi

echo ""
echo "=== Done ==="
echo "Reference path: ${REF_ROOT}/${REF_NAME}"
echo ""
echo "Pass this to cellranger-atac count via --reference."
