# Environment setup

## R (4.4.2)

```bash
Rscript environment/R_packages.R
```

Installs Seurat 5, Signac 1.10+, DESeq2, ComplexHeatmap, EnsDb.Mmusculus.v79, BSgenome, scDblFinder, harmony, and plotting helpers.

## Python (conda)

```bash
conda env create -f environment/python_env.yml
conda activate scenicplus
```

Provides `pycisTopic 1.0.3` (peak merging), `macs2`, `scanpy`, `pydeseq2` (aging-DAR re-calling), and `pyranges`.

## HPC modules

See `hpc_modules.txt`. The pipeline assumes a SLURM cluster with these modules available:
- `r/4.4.2`, `anaconda3/2023.09-0`, `sra-toolkit`, `bedtools/2.30.0`, `cellranger-atac/2.0.0`

## HOMER

Install HOMER once per user account:
```bash
mkdir -p /scratch/user/${USER}/homer && cd /scratch/user/${USER}/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
perl configureHomer.pl -install mm10
export HOMER_HOME=/scratch/user/${USER}/homer
export PATH=${HOMER_HOME}/bin:${PATH}
```
Most HOMER-running SLURM scripts in `07_HOMER/` set `HOMER_HOME` at the top - adjust there if your install path differs.

## Required environment variables

The scripts read these for portability:

| Variable | Default in scripts | Meaning |
|---|---|---|
| `DATA_ROOT` | `/QRISdata/Q8448/Mouse_disease_data` | Root of all raw + processed data |
| `REF_ROOT` | `/scratch/user/$USER/mm10_ref` | mm10 chrom.sizes + ENCODE blacklist |
| `HOMER_HOME` | `/scratch/user/$USER/homer` | HOMER install root |
| `CHROMSIZES_FILE` | `${REF_ROOT}/mm10.chrom.sizes` | UCSC mm10 chrom sizes |
| `BLACKLIST` | `${REF_ROOT}/mm10-blacklist.v2.bed` | ENCODE blacklist v2 |

Most scripts have these hard-coded to the original development paths. If you run from a different cluster account, edit the path constants at the top of each script (search for `/QRISdata/` and `/scratch/user/`).
