#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=download_cellranger_ref
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Download the 10x Cell Ranger ARC mm10 reference

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

echo "Reference path: ${REF_ROOT}/${REF_NAME}"
