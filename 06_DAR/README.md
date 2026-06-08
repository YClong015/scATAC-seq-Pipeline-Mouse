# Stage 06 - Pseudo-bulk DAR calling (Aim 2 input)

Per-cell-type, per-contrast differential accessibility via pseudo-bulk DESeq2, using the `apply_DESeq2_test_seurat()` wrapper in `DATesting.R` (Patrick Lab standard, mirrors Squair et al. 2021 recommendations).

## Per-tissue scripts (consistent: each reads its OWN per-tissue universal-peak obj)

| Tissue | Script | Input object (from stage 05) | Output dir |
|---|---|---|---|
| Kidney | `Kidney_pseudo-bulk_DAR.R` | `kidney_merged_universal.rds` | `DAR_pseudobulk_Kidney_DESeq2/` |
| Lung | `Lung_pseudo-bulk_DAR.R` | `lung_universal_new_pruned.rds` | `DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/` |
| Aorta | `Aorta_pseudo-bulk_DAR.R` | `Aorta_integrated_universal.rds` | `DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2/` |
| T cells | `Tcells_pseudo-bulk_DAR.R` | `tcells_universal.rds` | `DAR_pseudobulk_Tcells_DESeq2/` |

All four tissues now read their own per-tissue object (no combined integrated obj, no v5 figure-only obj) - fully consistent.

## Run

```bash
# Kidney
sbatch 06_DAR/kidney_dar.slurm        # invokes Kidney_pseudo-bulk_DAR.R

# Lung
Rscript 06_DAR/Lung_pseudo-bulk_DAR.R

# Aorta
Rscript 06_DAR/Aorta_pseudo-bulk_DAR.R

# T cells (mounts fragments + DAR - 9 pairwise contrasts among 5 conditions)
Rscript 06_DAR/Tcells_pseudo-bulk_DAR.R
```

## Tissue-specific contrasts

| Tissue | `group_col` | Contrasts |
|---|---|---|
| Kidney | `condition` | Day14 vs Sham, Day42 vs Sham, Day42 vs Day14 |
| Lung | `Group` | Case vs Control |
| Aorta | `Group` | Challenge vs Control |
| T cells | `deMultliplex2_final_mapped` | 9 pairwise among {Aged, Juvenile, Young_acute, Young_chronic, Young_control} |

## DESeq2 parameters (`apply_DESeq2_test_seurat()`)

| Parameter | Default | Meaning |
|---|---|---|
| `exp.thresh` | 0.05 | Pre-filter: peak must be accessible in >=5% of cells in foreground OR background |
| `num.splits` | 10 | Each population randomly split into N pseudo-bulk pools |
| `padj.cutoff` | 0.05 | BH FDR threshold for "significant" |
| `min.cells.per.group` | 80 | Excluded cell types below threshold |

A peak is called a DAR if `padj < 0.05`. Opening = log2FC > 0; closing = log2FC < 0.

## Outputs per cell type x contrast

```
{output_dir}/DAR_tables/{cell_type}__{contrast}__005_DESeq2_all.tsv         # all tested peaks
{output_dir}/DAR_tables/{cell_type}__{contrast}__005_DESeq2_padj005.tsv     # padj<0.05
{output_dir}/DAR_tables/{cell_type}__{contrast}__005_opening.tsv            # padj<0.05 & lfc>0
{output_dir}/DAR_tables/{cell_type}__{contrast}__005_closing.tsv            # padj<0.05 & lfc<0
{output_dir}/DAR_BED/{cell_type}__{contrast}__opening.bed                   # for HOMER
{output_dir}/DAR_BED/{cell_type}__{contrast}__closing.bed                   # for HOMER
{output_dir}/QC/*.csv                                                       # per-CT cell counts
{output_dir}/Figures/*pvalue_hist.png                                       # p-value diagnostic
```

These feed into stage 08 (HOMER motif enrichment).
