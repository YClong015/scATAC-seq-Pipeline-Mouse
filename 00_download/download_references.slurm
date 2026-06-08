#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=download_references
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Download mm10 chrom.sizes (UCSC) and ENCODE blacklist v2

set -euo pipefail

REF_ROOT="${REF_ROOT:-/scratch/user/$USER/mm10_ref}"
mkdir -p "${REF_ROOT}"

if [ ! -f "${REF_ROOT}/mm10.chrom.sizes" ]; then
  echo "Downloading mm10.chrom.sizes..."
  wget -O "${REF_ROOT}/mm10.chrom.sizes" \
    "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
else
  echo "Already present: ${REF_ROOT}/mm10.chrom.sizes"
fi

if [ ! -f "${REF_ROOT}/mm10-blacklist.v2.bed" ]; then
  echo "Downloading mm10 blacklist v2..."
  wget -O "${REF_ROOT}/mm10-blacklist.v2.bed.gz" \
    "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz"
  gunzip "${REF_ROOT}/mm10-blacklist.v2.bed.gz"
else
  echo "Already present: ${REF_ROOT}/mm10-blacklist.v2.bed"
fi

echo "Done. CHROMSIZES_FILE=${REF_ROOT}/mm10.chrom.sizes  BLACKLIST=${REF_ROOT}/mm10-blacklist.v2.bed"
ls -la "${REF_ROOT}"
