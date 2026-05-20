# Chromatin Accessibility and Transcription Factor Regulation in Age-Related Diseases

Cross-tissue single-cell ATAC-seq analysis of four mouse models of age-related disease — kidney ischemia-reperfusion injury (IRI), lung COPD, aortic aneurysm/dissection (AAD), and chronic T-cell stimulation — integrated against the organism-wide mouse aging chromatin atlas (Lu et al., 2026, *Science*).

This repository accompanies the BIOX7026 Master's thesis of **Yanchen Zheng** (UQ IMB, Christian Nefzger Lab supervised by Ralph Patrick, 2026).

---

## Central finding

| Layer | Programme | Driver TFs |
|---|---|---|
| Universal disease opening | Inflammatory enhancer opening | **AP-1** (Fos, Jun, Atf3, BATF) — shared across all four tissues |
| Tissue-specific disease closing | Loss of cell identity | HNF1b/HNF4a/PPARα (kidney PT); Smad2/HIF-1α (lung AT2); KLF4/WT1 (aorta SMC); TCF7/LEF1 (T cells) |
| Aging concordance | Disease-aging chromatin overlap | Kidney PT shows strongest concordance (Pearson r = 0.45–0.62) with Lu et al. 2026 aging atlas |

Disease chromatin remodelling is consistent with the SIPHON model (Patrick et al., 2024, *Cell Metabolism*): age-related disease accelerates the AP-1-driven aging chromatin trajectory in the same cell type.

---

## Pipeline overview

```
                                ┌─────────────────────────┐
                                │ 00_download             │  Raw FASTQ + Lu 2026 h5ad + mm10 ref
                                └────────────┬────────────┘
                                             │
                Kidney + Aorta ──► Cell Ranger ATAC                  Lung ──► dnbc4tools (MGI/BGI)
                T cells (private) ──► Ralph Patrick's pre-processed Tcells_Seurat_filtered.RData
                                             │
                                             ▼
                                ┌─────────────────────────┐
                                │ 01_preprocessing        │  QC + Harmony + cell-type
                                │  (one R script/tissue)  │  annotation + SplitFragments
                                └────────────┬────────────┘
                                             ▼
                                ┌─────────────────────────┐
                                │ 02_peak_calling         │  MACS2 (pycisTopic) per cell type
                                └────────────┬────────────┘  → narrowPeak
                                             ▼
                                ┌─────────────────────────┐
                                │ 03_peak_merging         │  Per-tissue consensus peaks
                                └────────────┬────────────┘
                                             ▼
                                ┌─────────────────────────┐
                                │ 04_universal_peaks      │  Universal peak set
                                └────────────┬────────────┘  (667,473 peaks; merge_universal_mm10.py)
                                             ▼
                                ┌─────────────────────────┐
                                │ 05_universal_assays     │  Per-tissue Seurat re-quantified
                                └────────────┬────────────┘  against the universal peak set
                                             ▼
                                ┌─────────────────────────┐
                                │ 06_DAR                  │  Pseudo-bulk DESeq2 (DATesting.R)
                                └────────────┬────────────┘  per cell type per tissue
                                             ▼
                                ┌─────────────────────────┐
                                │ 07_HOMER                │  Motif enrichment (NS + Stable bg)
                                └────────────┬────────────┘  4 backgrounds × 2 directions
                                             ▼
                                ┌─────────────────────────┐
                                │ 08_aging_comparison     │  bedtools intersect vs Lu 2026
                                │ (Aim 3)                 │  Fisher OR + Pearson r
                                └────────────┬────────────┘
                                             ▼
                                ┌─────────────────────────┐
                                │ 09_figures              │  All 21 thesis figures
                                └─────────────────────────┘
```

---

## ⚠ Tissue-specific platform note

**Lung uses MGI/BGI DNBSEQ sequencing — NOT 10x Chromium.** The 6 lung samples are aligned with `dnbc4tools atac run` (BGI's analog of Cell Ranger), not `cellranger-atac count`. See `00_download/lung_cngb/README.md` and `01_preprocessing/lung/dnbc4tools_per_sample/`.

Kidney + Aorta use 10x Chromium → Cell Ranger ATAC. T cells use 10x Multiome (Esmaeili PR2 PDF for wet-lab details), entry point is Ralph Patrick's pre-processed Seurat object.

---

## Quickstart

```bash
# 1. Clone repo
git clone https://github.com/YClong015/scatac-aging-disease.git
cd scatac-aging-disease

# 2. Install environments
Rscript environment/R_packages.R
conda env create -f environment/python_env.yml

# 3. Set data paths
export DATA_ROOT=/QRISdata/Q8448/Mouse_disease_data
export REF_ROOT=/scratch/user/$USER/mm10_ref
export HOMER_HOME=/scratch/user/$USER/homer

# 4. Download raw data (see 00_download/README.md)
bash   00_download/download_references.sh
sbatch 00_download/mkref_MGI.sh                          # dnbc4tools reference for Lung
sbatch 00_download/download_kidney_sra.sh
sbatch 00_download/download_aorta_sra.sh
sbatch 00_download/lung_cngb/download_full_directory.sh  # MGI / CNGB
bash   00_download/download_science_aging.sh             # manual
# T cells: request from r.patrick@uq.edu.au (see 00_download/tcells_HANDOVER.md)

# 5. Per-platform alignment
# Kidney + Aorta — Cell Ranger ATAC (template at 00_download/cellranger_template.md)
# Lung — dnbc4tools (run scripts at 01_preprocessing/lung/dnbc4tools_per_sample/)

# 6. Run pipeline stages 01-09 in order (each subdir's README has commands)
```

---

## Data sources

| Tissue | Source | Accession | Platform | Reference |
|---|---|---|---|---|
| Kidney IRI | NCBI SRA | GSE197391 | 10x Chromium | Muto et al., 2024 *Sci. Adv.* (doi:10.1126/sciadv.adk8845) |
| Lung COPD | CNGBdb | CNP0004399 | **MGI/BGI DNBSEQ** | Zhang Q. et al., 2025 *PLOS ONE* (doi:10.1371/journal.pone.0322538) |
| Aorta AAD | NCBI SRA | SRR21686722, SRR21686724 | 10x Chromium | Zhang C. et al., 2023 *ATVB* (doi:10.1161/ATVBAHA.122.318135) |
| T cells | Nefzger Lab (in-house) | private | 10x Multiome | Esmaeili et al., in prep — see `00_download/tcells_HANDOVER.md` |
| Mouse aging atlas | epiage.net / GEO | GSM8774006, GSM8774007 | EasySci-ATAC | Lu et al., 2026 *Science* (doi:10.1126/science.adw6273) |

---

## Repo layout

```
00_download/             — Raw data download (SRA, CNGB MGI, Science h5ad, mm10 + dnbc4tools refs)
01_preprocessing/        — Per-tissue alignment outs/ → Seurat (QC, Harmony, annotation)
   ├── kidney/Kidney_scATAC_Combine.R
   ├── lung/atac_Lung.Rmd  +  dnbc4tools_per_sample/Run_MGI_75..80.slurm
   ├── aorta/Aortic_scATAC.R
   └── tcells/Tcell_scATAC.R + Tcell_reannotate.R   (entry: Ralph's pre-processed RData)
02_peak_calling/         — pycisTopic + MACS2 per cell type
03_peak_merging/         — pycisTopic consensus per tissue
04_universal_peaks/      — Cross-tissue universal peak set (merge_universal_mm10.py)
05_universal_assays/     — Re-quantify each tissue's Seurat obj against universal peaks
                            (Lung gets an extra prune step via Diagnose_lung.R)
06_DAR/                  — Pseudo-bulk DESeq2 — one script per tissue, all using their own obj
07_HOMER/                — Motif enrichment (NS + Stable backgrounds)
08_aging_comparison/     — Aim 3 — bedtools intersect against Lu 2026 + Fisher + Pearson r
09_figures/              — All 21 thesis figures + supplementary
environment/             — R + Python + HPC module specs
```

Note: the 4-tissue Harmony integration step (used only to make the integrated UMAP figure for Aim 1) is not part of the canonical DAR pipeline. The script that produced the integrated UMAP lives under `09_figures/Fig7_UMAP/Fig_integration_UMAP.R`.

---

## Reproducibility notes

- All randomised steps use `set.seed(2024)` (R) or `random_state=2024` (Python).
- Software versions: R 4.4.2, Python 3.10.12, Seurat 5.0.1, Signac 1.10.0, DESeq2 1.40.2, pycisTopic 1.0.3, HOMER 4.11.1, MACS2 2.2.7.1, bedtools 2.30.0, Cell Ranger ATAC 2.1.0, dnbc4tools 3.0 (Lung MGI alignment).
- All file paths in scripts were developed on UQ Bunya HPC. Edit the constants at the top of each script (search `/QRISdata/` and `/scratch/user/`) to match your cluster.

## License

Code: MIT. Documentation: CC-BY-4.0. Data: subject to the licences of the underlying primary datasets.

## Contact

Yanchen Zheng (`yanchenzheng34@gmail.com`) · Patrick Lab, UQ IMB
