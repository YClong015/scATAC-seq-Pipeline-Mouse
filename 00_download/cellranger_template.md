# Cell Ranger ATAC (Kidney + Aorta only — NOT Lung)

Lung uses MGI/dnbc4tools instead — see `01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_*.slurm`.

Lung does NOT need Cell Ranger.

## Template

```bash
#!/bin/bash --login
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd

set -euo pipefail

# EDIT per sample
SAMPLE="SRR27367347_Kidney_atac"
FASTQ_DIR="${DATA_ROOT}/Kidney/SRR27367347_Fastq/SRR27367347"
OUT_DIR="${DATA_ROOT}/Kidney/cellranger_unpacked_data/${SAMPLE}"
REF="/scratch/user/${USER}/cellranger_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"

module load cellranger-atac/2.1.0
mkdir -p "${OUT_DIR}" && cd "${OUT_DIR}"

cellranger-atac count \
  --id="${SAMPLE}" \
  --reference="${REF}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE%%_*}" \
  --localcores=16 \
  --localmem=120
```

Outputs under `${OUT_DIR}/${SAMPLE}/outs/`:
- `fragments.tsv.gz`
- `singlecell.csv`
- `filtered_peak_bc_matrix.h5`
- `peaks.bed`

Submit one SLURM job per sample. Kidney has 9 samples (SRR27367330-32, 40, 44, 46, 47, 49, 51); Aorta has 2 (SRR21686722, SRR21686724).

## mm10 reference (one-time)

```bash
mkdir -p /scratch/user/$USER/cellranger_ref && cd $_
wget https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
tar -xzf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
```
