# Stage 01 - Per-tissue preprocessing

Per-tissue FASTQ / pre-processed obj -> annotated, QC-passed Seurat object -> per-cell-type fragment BED files (via `SplitFragments()`).

## Tissue-specific alignment

| Tissue | Aligner | Why |
|---|---|---|
| Kidney | Cell Ranger ATAC | 10x Chromium data |
| Aorta | Cell Ranger ATAC | 10x Chromium data |
| **Lung** | **`dnbc4tools atac run` (BGI)** | **MGI/DNBSEQ sequencing platform** (Zhang Q. et al. 2025 PLOS ONE - see `paper/Lung_mice_paper.pdf`) |
| T cells | n/a (pre-processed) | Starts from Ralph Patrick's `Tcells_Seurat_filtered.RData` |

## Inputs

| Tissue | Required input(s) |
|---|---|
| Kidney | `${DATA_ROOT}/Kidney/cellranger_unpacked_data/{SRR}_Kidney_atac/{SRR}_Kidney_atac/outs/` - 9 samples (3 Sham, 3 Day14, 3 Day42) |
| Lung | `${DATA_ROOT}/Lung/Lung_cellatac/{CL_id}/outs/` (built by dnbc4tools - 6 samples) - see `dnbc4tools_per_sample/Run_MGI_*.slurm` for FASTQ->outs |
| Aorta | `${DATA_ROOT}/Aorta/Aorta_cellranger_atac/{SRR}_Aorta_atac/outs/` - 2 samples (Control, Challenge) |
| T cells | `${DATA_ROOT}/Tcells/Tcells_Seurat_filtered.RData` + `atac_fragments.tsv.gz` (from Patrick Lab - see `00_download/tcells_HANDOVER.md`) |

## Run order

```bash
# Lung - special: first run dnbc4tools on each FASTQ pair to produce outs/
# (six SLURM jobs, one per sample - submit all in parallel)
for i in 75 76 77 78 79 80; do
  sbatch 01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_${i}.slurm
done

# Then build the merged Lung Seurat object + cell-type annotation + SplitFragments
# (run interactively via RStudio "Knit" or Rscript -e 'rmarkdown::render(...)')
Rscript 01_preprocessing/lung/atac_Lung.R

# Kidney (Combine pipeline: per-sample LSI + Harmony integration -> 11 cell types)
sbatch 01_preprocessing/kidney/run_kidney_combine.slurm

# Aorta (rLSI integration anchors -> 6 cell types)
Rscript 01_preprocessing/aorta/Aortic_scATAC.R

# T cells - initial annotation, then contamination cleanup + final annotation
Rscript 01_preprocessing/tcells/Tcell_scATAC.R               # initial QC, Harmony, clustering, annotation
Rscript 01_preprocessing/tcells/Tcell_clean_recluster.R      # remove B/myeloid/endo contamination, re-embed
Rscript 01_preprocessing/tcells/Tcell_final_annotate.R       # final Masopust cell-type assignment + UMAP
Rscript 01_preprocessing/tcells/Tcell_split_fragments_final.R # split fragments by final cell type (for peak calling)
```

## Outputs

| Tissue | Output | Used by |
|---|---|---|
| Kidney | `kidney_merged_annotated.rds` + `fragment_files_split_by_celltype/*.bed` | 02_peak_calling/kidney, 05_universal_assays/kidney |
| Lung | `lung_integrated_clean_annotated.rds` + `Lung_fragments_file/*.bed` | 02_peak_calling/lung, 05_universal_assays/lung |
| Aorta | `aortic_integrated_res0.6_up2k_seed*.rds` + `fragment_files/*.bed` | 02_peak_calling/aorta, 05_universal_assays/aorta |
| T cells | `tcells_processed.rds` (DAR input) -> 05_universal_assays/tcells; `tcells_final_annotated.rds` + `fragment_files_split_by_celltype_final/*.bed` -> 02_peak_calling/tcells | as noted |

The `SplitFragments()` step (inside each tissue's R script) writes one per-cell-type `.bed.gz` per tissue, which is what MACS2 reads in stage 02.

## Cell-type vocabulary

| Tissue | Cell types after annotation |
|---|---|
| Kidney | PCT, PST, Injured_PT, TAL, DCT_CNT, DTL_ATL, PC_URO, IC, EC, FIB, PODO_PEC, LEUK |
| Lung | AT2, B, Ciliated, EC-vasc, Eosinophils, Fib, Mac-alv, Mac-inter, Mesothelial, Mo-Ly6c+, NK, Pen, SMCs, T (Low Quality cluster removed) |
| Aorta | Endothelial, Fibroblast, Mac, Pericyte, SMC, T-cell |
| T cells | Naive_central_T, Activated_Cd69_T, Effector_SLEC_CD8, Exhausted_Tex_CD8, Treg, TRM_CD8_CD103, TRM_CD8_CD49a, Tpex_Tfh_like |

### T-cell re-annotation (QC update)

The initial T-cell annotation (`Tcell_scATAC.R`) was revised after marker-level
QC. Because the T-cell dataset is **10x Multiome (RNA + ATAC)**, the real GEX
(`SCT` assay; ~471 genes/cell) was used - rather than ATAC gene-activity - to
validate cluster identity against the consensus T-cell nomenclature
(Masopust et al., *Nat Rev Immunol* 2026). `FindAllMarkers` showed that several
clusters in the original annotation were contaminating non-T cells carried over
from the CD8 enrichment: two B-cell clusters (Pax5/Ebf1/Cd79a/b), one plasma-cell
cluster (Sdc1/Jchain/Igkc; originally mis-labelled "NK"), and myeloid/endothelial
residuals (Adgre1/Pparg; Vwf). These were removed and the ~6,800 bona fide T cells
re-clustered and re-annotated:

- `Tcell_clean_recluster.R` - remove contamination, re-embed, FindAllMarkers + Masopust-marker z-score
- `Tcell_final_annotate.R` - final cell-type assignment + annotated UMAP

Marker basis (Masopust 2026 backbone + additional primary sources where needed):
SLEC/effector (Klrg1/Cx3cr1/Tbx21/Gzmb; Zeb2 - Omilusik et al. 2015 *JEM*),
exhausted/Tex (Tox/TIM3/PD-1/CD39; Maf, Nr4a2 - Giordano et al. 2015 *EMBO J*),
Treg (Foxp3/CD25/CTLA4/Helios/GITR), TRM (CD103/CD49a/CXCR6).
