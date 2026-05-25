# Stage 01 — Per-tissue preprocessing

Per-tissue FASTQ / pre-processed obj → annotated, QC-passed Seurat object → per-cell-type fragment BED files (via `SplitFragments()`).

## Tissue-specific alignment

| Tissue | Aligner | Why |
|---|---|---|
| Kidney | Cell Ranger ATAC | 10x Chromium data |
| Aorta | Cell Ranger ATAC | 10x Chromium data |
| **Lung** | **`dnbc4tools atac run` (BGI)** | **MGI/DNBSEQ sequencing platform** (Zhang Q. et al. 2025 PLOS ONE — see `paper/Lung_mice_paper.pdf`) |
| T cells | n/a (pre-processed) | Starts from Ralph Patrick's `Tcells_Seurat_filtered.RData` |

## Inputs

| Tissue | Required input(s) |
|---|---|
| Kidney | `${DATA_ROOT}/Kidney/cellranger_unpacked_data/{SRR}_Kidney_atac/{SRR}_Kidney_atac/outs/` — 9 samples (3 Sham, 3 Day14, 3 Day42) |
| Lung | `${DATA_ROOT}/Lung/Lung_cellatac/{CL_id}/outs/` (built by dnbc4tools — 6 samples) — see `dnbc4tools_per_sample/Run_MGI_*.slurm` for FASTQ→outs |
| Aorta | `${DATA_ROOT}/Aorta/Aorta_cellranger_atac/{SRR}_Aorta_atac/outs/` — 2 samples (Control, Challenge) |
| T cells | `${DATA_ROOT}/Tcells/Tcells_Seurat_filtered.RData` + `atac_fragments.tsv.gz` (from Patrick Lab — see `00_download/tcells_HANDOVER.md`) |

## Run order

```bash
# Lung — special: first run dnbc4tools on each FASTQ pair to produce outs/
# (six SLURM jobs, one per sample — submit all in parallel)
for i in 75 76 77 78 79 80; do
  sbatch 01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_${i}.slurm
done

# Then build the merged Lung Seurat object + cell-type annotation + SplitFragments
# (run interactively via RStudio "Knit" or Rscript -e 'rmarkdown::render(...)')
Rscript -e 'rmarkdown::render("01_preprocessing/lung/atac_Lung.Rmd")'

# Kidney (Combine pipeline: per-sample LSI + Harmony integration → 11 cell types)
sbatch 01_preprocessing/kidney/run_kidney_combine.sh

# Aorta (rLSI integration anchors → 6 cell types)
Rscript 01_preprocessing/aorta/Aortic_scATAC.R

# T cells — two-step
Rscript 01_preprocessing/tcells/Tcell_scATAC.R       # initial QC, Harmony, clustering, annotation
Rscript 01_preprocessing/tcells/Tcell_reannotate.R   # corrected annotations + remove B/LowQ
```

## Outputs

| Tissue | Output | Used by |
|---|---|---|
| Kidney | `kidney_merged_annotated.rds` + `fragment_files_split_by_celltype/*.bed` | 02_peak_calling/kidney, 05_universal_assays/kidney |
| Lung | `lung_integrated_clean_annotated.rds` + `Lung_fragments_file/*.bed` | 02_peak_calling/lung, 05_universal_assays/lung |
| Aorta | `aortic_integrated_res0.6_up2k_seed*.rds` + `fragment_files/*.bed` | 02_peak_calling/aorta, 05_universal_assays/aorta |
| T cells | `tcells_processed.rds` + `fragment_files_split_by_celltype/*.bed` | 02_peak_calling/tcells, 05_universal_assays/tcells |

The `SplitFragments()` step (inside each tissue's R script) writes one per-cell-type `.bed.gz` per tissue, which is what MACS2 reads in stage 02.

## Cell-type vocabulary

| Tissue | Cell types after annotation |
|---|---|
| Kidney | PCT, PST, Injured_PT, TAL, DCT_CNT, DTL_ATL, PC_URO, IC, EC, FIB, PODO_PEC, LEUK |
| Lung | AT2, B, Ciliated, EC-vasc, Eosinophils, Fib, Mac-alv, Mac-inter, Mesothelial, Mo-Ly6c+, NK, Pen, SMCs, T (Low Quality cluster removed) |
| Aorta | Endothelial, Fibroblast, Mac, Pericyte, SMC, T-cell |
| T cells | Naive_T, Naive_CD8_T, Effector_CD8_T, Cytotoxic_CD8_T, CD8_Eff, Memory_CD8_T, Treg, Tfh_like_T, NK (B_cell + Low_quality removed) |
