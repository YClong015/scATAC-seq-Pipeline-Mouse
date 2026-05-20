#!/bin/bash
# ============================================================
# Download mm10 reference files used throughout the pipeline:
#   - mm10.chrom.sizes (UCSC)
#   - ENCODE mm10 blacklist v2
# ============================================================

set -euo pipefail

REF_ROOT="${REF_ROOT:-/scratch/user/$USER/mm10_ref}"
mkdir -p "${REF_ROOT}"

# ── mm10 chromosome sizes ──────────────────────────────────────
if [ ! -f "${REF_ROOT}/mm10.chrom.sizes" ]; then
  echo "Downloading mm10.chrom.sizes..."
  wget -O "${REF_ROOT}/mm10.chrom.sizes" \
    "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
else
  echo "Already present: ${REF_ROOT}/mm10.chrom.sizes"
fi

# ── ENCODE mm10 blacklist v2 ───────────────────────────────────
if [ ! -f "${REF_ROOT}/mm10-blacklist.v2.bed" ]; then
  echo "Downloading mm10 blacklist v2..."
  wget -O "${REF_ROOT}/mm10-blacklist.v2.bed.gz" \
    "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz"
  gunzip "${REF_ROOT}/mm10-blacklist.v2.bed.gz"
else
  echo "Already present: ${REF_ROOT}/mm10-blacklist.v2.bed"
fi

echo ""
echo "=== Done ==="
echo "REF_ROOT=${REF_ROOT}"
ls -la "${REF_ROOT}"
echo ""
echo "Set this environment variable so pipeline scripts can find these files:"
echo "  export CHROMSIZES_FILE=${REF_ROOT}/mm10.chrom.sizes"
echo "  export BLACKLIST=${REF_ROOT}/mm10-blacklist.v2.bed"
