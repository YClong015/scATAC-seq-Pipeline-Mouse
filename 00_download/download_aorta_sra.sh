#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD_aorta
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --array=0-1
#SBATCH --output=DWLD_aorta_%A_%a.out
#SBATCH --error=DWLD_aorta_%A_%a.err

# ============================================================
# Aorta AAD scATAC-seq — Zhang C. et al. 2023 ATVB
# 2 samples: Control (SRR21686724), Challenge (SRR21686722)
# ============================================================

set -euo pipefail

SRAIDS=(SRR21686724 SRR21686722)
SRAID=${SRAIDS[$SLURM_ARRAY_TASK_ID]}

DATA_ROOT="${DATA_ROOT:-/QRISdata/Q8448/Mouse_disease_data}"
DATADIR="${DATA_ROOT}/Aorta"
mkdir -p "${DATADIR}"

module load sra-toolkit

cd $TMPDIR
mkdir -p ${SRAID}_Fastq
cd ${SRAID}_Fastq

OUTDIR=${TMPDIR}/${SRAID}_Fastq/${SRAID}
fasterq-dump -e 8 -O ${OUTDIR} --progress --split-files --include-technical ${SRAID}

cd ${SRAID}
gzip *

cd ../../
tar -zcvf ${SRAID}_Fastq.tar.gz ${SRAID}_Fastq
cp ${SRAID}_Fastq.tar.gz ${DATADIR}/${SRAID}_Fastq.tar.gz

echo "Done: ${SRAID}"
