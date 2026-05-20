# Stage 01 — Per-tissue preprocessing

Per-tissue Cell Ranger outputs → annotated, QC-passed Seurat object → per-cell-type fragment BED files (via `SplitFragments()`).

## Inputs

| Tissue | Required input |
|---|---|
| Kidney | `${DATA_ROOT}/Kidney/cellranger_unpacked_data/{SRR}_Kidney_atac/{SRR}_Kidney_atac/outs/` — 9 samples (3 Sham, 3 Day14, 3 Day42) |
| Lung | `${DATA_ROOT}/Lung/Lung_cellatac/{CL...}/outs/` — 6 samples |
| Aorta | `${DATA_ROOT}/Aorta/Aorta_cellranger_atac/{SRR}_Aorta_atac/outs/` — 2 samples (Control, Challenge) |
| T cells | `${DATA_ROOT}/Tcells/Tcells_Seurat_filtered.RData` + `atac_fragments.tsv.gz` (from Patrick Lab — see `00_download/tcells_HANDOVER.md`) |

> **Lung note:** there is no preprocessing script for lung in this repo. The lung Seurat object was constructed externally (by the Zhang et al. 2025 authors / Patrick Lab); we load `lung_universal_new_pruned.rds` directly at stage 05.

## Run order

```bash
# Kidney (uses Combine pipeline: per-sample LSI + Harmony integration → 11 cell types)
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
| Kidney | `kidney_merged_annotated.rds` + `fragment_files_split_by_celltype/*.bed` | 02_peak_calling/kidney |
| Aorta | `aortic_integrated_res0.6_up2k_seed*.rds` + `fragment_files/*.bed` | 02_peak_calling/aorta |
| T cells | `tcells_processed.rds` + `fragment_files_split_by_celltype/*.bed` | 02_peak_calling/tcells |
| Lung | (n/a — pre-built `lung_universal_new_pruned.rds`) | 05_universal_assays/lung |

The `SplitFragments()` step writes one per-cell-type `.bed.gz` per tissue, which is what MACS2 reads in stage 02.

## Cell-type vocabulary

| Tissue | Cell types after annotation |
|---|---|
| Kidney | PT, Injured_PT, TAL, DCT_CNT, DTL_ATL, PC_URO, IC, EC, Pen, PODO_PEC, LEUK |
| Lung | AT2, B, Ciliated, EC-vasc, Eosinophils, Fib, Mac, Mac-alv, Mac-inter, Mesothelial, Mo-Ly6c+, NK, Pen, SMCs, T |
| Aorta | Endothelial, Fibroblast, Mac, Pericyte, SMC, T-cell |
| T cells | Naive_T, Naive_CD8_T, Effector_CD8_T, Cytotoxic_CD8_T, CD8_Eff, Memory_CD8_T, Treg, Tfh_like_T, NK (B_cell + Low_quality removed) |
