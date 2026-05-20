# Chromatin Accessibility and Transcription Factor Regulation in Age-Related Diseases

Cross-tissue single-cell ATAC-seq analysis of four mouse models of age-related disease — kidney ischemia-reperfusion injury (IRI), lung COPD, aortic aneurysm/dissection (AAD), and chronic T-cell stimulation — integrated against the organism-wide mouse aging chromatin atlas (Lu et al., 2026, *Science*).

This repository accompanies the BIOX7026 Master's thesis of **Yanchen Zheng** (UQ IMB, Christian Lab supervised by Ralph Patrick, 2026).

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
                                         ▼
   Kidney  Lung  Aorta             ┌─────────────────────────┐
   (Cell Ranger ATAC outs/)        │ 01_preprocessing        │  QC + Harmony + cell-type
                  +                │  (one R script/tissue)  │  annotation + SplitFragments
   T cells: Tcells_Seurat_         └────────────┬────────────┘
   filtered.RData (Ralph)                       │
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
                            │ 04_universal_peaks      │  Build universal peak set
                            └────────────┬────────────┘  (SCPM-filtered, v5: 667,473)
                                         ▼
                            ┌─────────────────────────┐
                            │ 05_universal_assays     │  Per-tissue Seurat re-quantified
                            └─────┬───────────────┬───┘  against the universal peak set
                                  │               │
                                  ▼               ▼
                            ┌──────────┐    ┌────────────────────┐
                            │ 06_integ │    │ 07_DAR             │  Pseudo-bulk DESeq2
                            │ Harmony  │    │ per cell type      │  (DATesting.R)
                            └────┬─────┘    └──────────┬─────────┘
                                 ▼                     ▼
                            ┌─────────────────────────┐ Opening / closing DARs
                            │ 08_HOMER                │  4 backgrounds × 2 directions
                            └────────────┬────────────┘  per cell type per tissue
                                         ▼
                            ┌─────────────────────────┐
                            │ 09_aging_comparison     │  bedtools intersect vs Lu 2026
                            │ (Aim 3)                 │  Fisher OR + Pearson r
                            └────────────┬────────────┘
                                         ▼
                            ┌─────────────────────────┐
                            │ 10_figures              │  All 21 thesis figures
                            └─────────────────────────┘
```

---

## Quickstart

```bash
# 1. Clone repo
git clone https://github.com/<user>/scatac-aging-disease.git
cd scatac-aging-disease

# 2. Install environments
Rscript environment/R_packages.R
conda env create -f environment/python_env.yml

# 3. Set data paths (or just edit the constants at the top of each script)
export DATA_ROOT=/QRISdata/Q8448/Mouse_disease_data
export REF_ROOT=/scratch/user/$USER/mm10_ref
export HOMER_HOME=/scratch/user/$USER/homer

# 4. Download raw data (see 00_download/README.md)
bash   00_download/download_references.sh
bash   00_download/download_cellranger_ref.sh
sbatch 00_download/download_kidney_sra.sh
sbatch 00_download/download_aorta_sra.sh
bash   00_download/download_lung_cngb.sh         # manual download from CNGBdb
bash   00_download/download_science_aging.sh
# T cells: request from r.patrick@uq.edu.au (see 00_download/tcells_HANDOVER.md)

# 5. Cell Ranger ATAC per sample
sbatch 00_download/cellranger_template.sh         # edit SAMPLE / FASTQ_DIR per sample

# 6. Run pipeline stages in order (see each subdir's README.md)
#    Each stage is fully scripted via SLURM .slurm / .sbatch files.
```

---

## Data sources

| Tissue | Source | Accession | Reference |
|---|---|---|---|
| Kidney IRI | NCBI SRA | GSE197391 | Muto et al., 2024 *Sci. Adv.* (doi:10.1126/sciadv.adk8845) |
| Lung COPD | CNGB | CNP0004399 | Zhang Q. et al., 2025 *PLOS ONE* (doi:10.1371/journal.pone.0322538) |
| Aorta AAD | NCBI SRA | SRR21686722, SRR21686724 | Zhang C. et al., 2023 *ATVB* (doi:10.1161/ATVBAHA.122.318135) |
| T cells | Patrick Lab (in-house) | private | Esmaeili et al., in prep |
| Mouse aging atlas | epiage.net / GEO | GSM8774006, GSM8774007 | Lu et al., 2026 *Science* (doi:10.1126/science.adw6273) |

---

## Repo layout

```
00_download/             — Data download scripts (SRA, CNGB, references)
01_preprocessing/        — Per-tissue cellranger-out → Seurat (QC, Harmony, annotation)
   ├── kidney/Kidney_scATAC_Combine.R
   ├── aorta/Aortic_scATAC.R
   └── tcells/Tcell_scATAC.R + Tcell_reannotate.R
       (no lung script: Lung pre-processed obj loaded as-is — see 05_universal_assays/lung/)
02_peak_calling/         — pycisTopic + MACS2 per cell type (per-tissue, plus 9 SRR per-sample for kidney)
03_peak_merging/         — pycisTopic consensus per tissue
04_universal_peaks/      — Two-stage merge + SCPM>1 filter → universal v5
05_universal_assays/     — Re-quantify each tissue's Seurat obj against universal v5 peaks
06_integration/          — 4-tissue Harmony integration + re-annotation + UMAP
07_DAR/                  — Pseudo-bulk DESeq2 DAR calling via DATesting.R
08_HOMER/                — Motif enrichment (8 comparisons × cell type)
   ├── prepare_bed/      — Generate opening/closing/stable/NS BED files
   ├── run_homer/        — Canonical 4-tissue + per-tissue HOMER runs
   └── per_tissue_initial/ — Initial whole-DAR HOMER runs per tissue
09_aging_comparison/     — Aim 3 — bedtools intersect against Lu 2026 atlas
   ├── 02_aging_DAR/     — Re-call aging DARs (Aged-vs-Young + 21mo-vs-5mo variants)
   ├── 03_overlap/       — Fisher OR + Pearson r + scatter plots
   ├── 04_region_classification.sh — disease-specific / shared / aging-specific
   ├── methodology_tests/  — sex covariate / CPM-filter sensitivity tests
   └── replots/          — Final replots
10_figures/              — All 21 thesis figures + supplementary figures
environment/             — R + Python + HPC module specs
thesis_reference/        — Thesis markdown + figure structure outline
```

---

## Citation

If you use this pipeline or its results, please cite the thesis (forthcoming) and the underlying datasets:

> Zheng Y. (2026). *Chromatin Accessibility and Transcription Factor Regulation in Age-Related Diseases*. MBioinformatics thesis, University of Queensland.

Please also cite Muto et al. 2024, Zhang Q. 2025, Zhang C. 2023, Lu et al. 2026, and Patrick et al. 2024 — see `thesis_reference/Thesis_Yanchen_Zheng_FULL.md` §References for full bibliographic detail.

---

## Reproducibility notes

- All randomised steps use `set.seed(2024)` (R) or `random_state=2024` (Python).
- Software versions: R 4.4.2, Python 3.10.12, Seurat 5.0.1, Signac 1.10.0, DESeq2 1.40.2, pycisTopic 1.0.3, HOMER 4.11.1, MACS2 2.2.7.1, bedtools 2.30.0, Cell Ranger ATAC 2.1.0.
- All file paths in scripts were developed on UQ Bunya HPC. Edit the constants at the top of each script (search `/QRISdata/` and `/scratch/user/`) to match your cluster.

## License

Code: MIT. Documentation: CC-BY-4.0. Data: subject to the licences of the underlying primary datasets.

## Contact

Yanchen Zheng (`yanchenzheng34@gmail.com`) · Patrick Lab, UQ IMB
