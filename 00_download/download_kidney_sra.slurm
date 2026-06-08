#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD_kidney
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --array=0-8
#SBATCH --output=DWLD_kidney_%A_%a.out
#SBATCH --error=DWLD_kidney_%A_%a.err

# Kidney IRI scATAC-seq, Muto et al. 2024 (GSE197391), 9 samples (Sham/Day14/Day42)

set -euo pipefail

SRAIDS=(SRR27367347 SRR27367332 SRR27367344 \
        SRR27367351 SRR27367349 SRR27367346 \
        SRR27367331 SRR27367330 SRR27367340)
SRAID=${SRAIDS[$SLURM_ARRAY_TASK_ID]}

DATA_ROOT="${DATA_ROOT:-/QRISdata/Q8448/Mouse_disease_data}"
DATADIR="${DATA_ROOT}/Kidney"
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
