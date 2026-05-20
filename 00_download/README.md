# Stage 00 — Data Download

All raw data sources used in this project. Run these scripts on the HPC node that has internet access.

| Tissue / Dataset | Source | Accession | Sequencing platform | Download script |
|---|---|---|---|---|
| Kidney IRI scATAC-seq | NCBI SRA | GSE197391 (9 SRR runs) | 10x Chromium | `download_kidney_sra.sh` |
| Lung COPD scATAC-seq | CNGBdb | CNP0004399 (6 samples) | **MGI / BGI DNBSEQ** | `lung_cngb/download_full_directory.sh` (recursive) or `lung_cngb/download_per_sample.sh` (one URL at a time) |
| Aorta AAD scATAC-seq | NCBI SRA | SRR21686722, SRR21686724 | 10x Chromium | `download_aorta_sra.sh` |
| T-cell exhaustion multiome | Nefzger Lab (private) | embargoed | 10x Chromium Multiome | See `tcells_HANDOVER.md` |
| Mouse aging chromatin atlas | GEO / epiage.net | GSM8774006 (Lung), GSM8774007 (Kidney) | EasySci-ATAC (Lu et al. 2026) | `download_science_aging.sh` |
| mm10 reference + blacklist | UCSC / ENCODE | mm10.chrom.sizes, ENCODE blacklist v2 | n/a | `download_references.sh` |
| dnbc4tools mm10 reference | local build from GENCODE vM25 | n/a | for MGI Lung data | `mkref_MGI.sh` |

---

## ⚠ Important: Lung uses MGI sequencing, not 10x Chromium

The lung COPD dataset (Zhang Q. et al., 2025 PLOS ONE; original paper PDF at `/paper/Lung_mice_paper.pdf`) was generated on the **MGI/BGI DNBSEQ platform** using **dnbc4tools** (BGI's analog of Cell Ranger) — NOT 10x Genomics.

This means:
- The 6 Lung samples are named with BGI run IDs: `CL100168054_L01`, `CL100167942_L01`, `CL100168054_L02`, `CL100168078_L02`, `CL100167942_L02`, `CL100168078_L01`
- FASTQ → fragments/peak-matrix is done by `dnbc4tools atac run` (see `01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_*.slurm`), NOT `cellranger-atac count`
- The mm10 reference also needs to be built via `dnbc4tools atac mkref` (see `mkref_MGI.sh`)

All Kidney, Aorta, and T-cell datasets use 10x Chromium and standard Cell Ranger ATAC (run separately — see Step 4 below).

---

## Quick start

```bash
# Set the data root (canonical: /QRISdata/Q8448/Mouse_disease_data on UQ HPC)
export DATA_ROOT=/QRISdata/Q8448/Mouse_disease_data
export REF_ROOT=/scratch/user/$USER/mm10_ref

# 1. References
bash   00_download/download_references.sh
sbatch 00_download/mkref_MGI.sh                  # dnbc4tools mm10 reference (~30 min)

# 2. Raw data per tissue
sbatch 00_download/download_kidney_sra.sh        # 9 SRR (SLURM array)
sbatch 00_download/download_aorta_sra.sh         # 2 SRR (SLURM array)
sbatch 00_download/lung_cngb/download_full_directory.sh   # CNGB recursive wget (6 samples × R1/R2)
bash   00_download/download_science_aging.sh     # Lu 2026 h5ad (manual — see script)

# T cells: request from Patrick / Nefzger lab — see tcells_HANDOVER.md
```

### Step 4 — FASTQ → fragments

After raw FASTQ is on disk, run the per-platform alignment:

| Tissue | Aligner | Reference | Script |
|---|---|---|---|
| Kidney, Aorta | `cellranger-atac count` v2.1.0 | 10x Genomics mm10 reference (refdata-cellranger-arc-mm10-2020-A-2.0.0) | run manually per sample — see `cellranger_template.md` |
| Lung | `dnbc4tools atac run` v3.0 | locally built via `mkref_MGI.sh` | `01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_*.slurm` |
| T cells | already pre-processed | n/a | obj handover via Patrick Lab |
