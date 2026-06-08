#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH --array=0
#SBATCH --output=DWLD_%j.out
#SBATCH --error=DWLD_%j.err

URLS=(
"ftp://ftp2.cngb.org/pub/CNSA/data6/CNP0004399/CNS0809369/CNX0739875/CNR0841131/CL100168078_L02_read_2.fq.gz"
)
URL=${URLS[$SLURM_ARRAY_TASK_ID]}
FILE=$(basename "$URL")
BASE=${FILE%.fq.gz}

DATADIR=/QRISdata/Q8448/Mouse_disease_data/Lung

cd $TMPDIR
mkdir ${BASE}_Fastq
cd ${BASE}_Fastq
wget -O "$FILE" "$URL"


cd ../
tar -zcvf ${BASE}_Fastq.tar.gz ${BASE}_Fastq

cp ${BASE}_Fastq.tar.gz ${DATADIR}/${BASE}_Fastq.tar.gz
