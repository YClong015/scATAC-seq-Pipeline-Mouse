#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=dnbc4tools_Lung
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH --output=MGI_%j.out
#SBATCH --error=MGI_%j.err
#SBATCH --array=1-1

DATA=Lung
SAMPLE=CL100168054_L01
SOURCE_DIR="/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_raw_data/CNX0739876/CNR0841132"
FQ_R1="CL100168054_L01_read_1.fq.gz"
FQ_R2="CL100168054_L01_read_2.fq.gz"

DNBC4TOOLS_EXE="/scratch/user/s4869245/dnbc4tools/dnbc4tools3.0beta_v5/dnbc4tools"
GENOME_REF_DIR="/scratch/user/s4869245/mouse_mm10_ref"
OUTPUT_DIR="/QRISdata/Q8448/Mouse_disease_cellatac/${DATA}_dnbc4tools"

mkdir -p $OUTPUT_DIR

cd $TMPDIR

# --- Data Preparation ---
cp ${SOURCE_DIR}/${FQ_R1} ${TMPDIR}/${SAMPLE}_R1.fastq.gz
cp ${SOURCE_DIR}/${FQ_R2} ${TMPDIR}/${SAMPLE}_R2.fastq.gz

# --- Core Analysis ---
echo "Starting dnbc4tools atac run for sample ${SAMPLE} at $(date)"
$DNBC4TOOLS_EXE atac run \
        --name ${SAMPLE} \
        --fastq1 ${TMPDIR}/${SAMPLE}_R1.fastq.gz \
        --fastq2 ${TMPDIR}/${SAMPLE}_R2.fastq.gz \
        --genomeDir ${GENOME_REF_DIR} \
        --threads 32 \
	--chemistry scATAC

# --- Packaging ---
echo "Packaging results for sample ${SAMPLE} at $(date)"
tar -zcvf $TMPDIR/${SAMPLE}.tar.gz ${TMPDIR}/${SAMPLE}
cp ${TMPDIR}/${SAMPLE}.tar.gz ${OUTPUT_DIR}/${DATA}_${SAMPLE}.tar.gz
 
