# Stage 00 — Data Download

All raw data sources used in this project. Run these scripts on the HPC node that has internet access.

| Tissue / Dataset | Source | Accession | Download script |
|---|---|---|---|
| Kidney IRI scATAC-seq | NCBI SRA | GSE197391 (9 SRR runs) | `download_kidney_sra.sh` |
| Lung COPD scATAC-seq | CNGB | CNP0004399 (6 samples) | `download_lung_cngb.sh` |
| Aorta AAD scATAC-seq | NCBI SRA | SRR21686722, SRR21686724 | `download_aorta_sra.sh` |
| T-cell exhaustion multiome | Nefzger/Patrick lab (private) | See `tcells_HANDOVER.md` | not applicable |
| Mouse aging chromatin atlas | Lu et al., 2026 *Science* | GSM8774006 (Lung), GSM8774007 (Kidney) | `download_science_aging.sh` |
| mm10 reference + blacklist | UCSC / ENCODE | mm10.chrom.sizes, ENCODE blacklist v2 | `download_references.sh` |
| Cell Ranger ATAC mm10 reference | 10x Genomics | refdata-cellranger-arc-mm10-2020-A-2.0.0 | `download_cellranger_ref.sh` |

After SRA download you still need to run **Cell Ranger ATAC** on the FASTQ files to produce `fragments.tsv.gz` + `filtered_peak_bc_matrix.h5`. See `cellranger_template.sh`.

---

## Quick start

```bash
# Set the data root (canonical: /QRISdata/Q8448/Mouse_disease_data on UQ HPC)
export DATA_ROOT=/QRISdata/Q8448/Mouse_disease_data

# Download SRA tissues (parallel SLURM array jobs)
sbatch 00_download/download_kidney_sra.sh
sbatch 00_download/download_aorta_sra.sh
bash   00_download/download_lung_cngb.sh    # CNGB — single download, large
bash   00_download/download_science_aging.sh
bash   00_download/download_references.sh
bash   00_download/download_cellranger_ref.sh

# Once FASTQ is on disk, run cellranger-atac per sample. Template:
bash 00_download/cellranger_template.sh
```

For the T-cell dataset, request `Tcells_Seurat_filtered.RData` from the Patrick lab — see `tcells_HANDOVER.md`. This is the entry point for `01_preprocessing/tcells/`.
