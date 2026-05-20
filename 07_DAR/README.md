# Stage 07 — Pseudo-bulk DAR calling (Aim 2 input)

Per-cell-type, per-contrast differential accessibility via pseudo-bulk DESeq2, using the `apply_DESeq2_test_seurat()` wrapper in `DATesting.R` (Patrick Lab standard, mirrors Squair et al. 2021 recommendations).

## Per-tissue scripts and the objects they read

| Tissue | Script | Input object | Output dir |
|---|---|---|---|
| **Kidney** | `Kidney_specific_DAR_pseudo-bulk.R` | `kidney_universal_v5.rds` | `DAR_pseudobulk_Kidney_v5_DESeq2/` |
| **Lung** | `Lung_specific_DAR_pseudo-bulk.R` | `lung_universal_new_pruned.rds` ⚠ (older object) | `DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/` |
| **Aorta** | `Aorta_Lung_Tcells_pseudo-bulk_DAR.R` (the per-tissue loop subsets the integrated object) | `All_Tissues_Integrated_Annotated_Clean_for_DAR.rds` | `DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2/` |
| **T cells** | `Tcells_pseudo-bulk_DAR.R` | `tcells_universal.rds` (note: pre-v5 — see ⚠ below) | `DAR_pseudobulk_Tcells_DESeq2/` |

> ⚠ **Lung & T-cell DAR were performed against pre-v5 objects.** This was an intentional decision — the v5 re-quantification happened after the DAR calls were finalised for figure generation. If you re-run from scratch, either (a) re-call DAR on the v5 objects, or (b) keep this same configuration to reproduce the thesis figures exactly. The downstream HOMER / aging-comparison pipelines read from these specific output dirs.

## Run

```bash
# Kidney (recommended SLURM)
sbatch 07_DAR/kidney_dar.slurm        # invokes Kidney_specific_DAR_pseudo-bulk.R

# Lung (single R job, ~6 h)
Rscript 07_DAR/Lung_specific_DAR_pseudo-bulk.R

# Aorta (via the Aorta_Lung_Tcells per-tissue loop — restrict to Aorta if you only need Aorta)
Rscript 07_DAR/Aorta_Lung_Tcells_pseudo-bulk_DAR.R

# T cells (mounts fragments + DAR — 9 pairwise contrasts among 5 conditions)
Rscript 07_DAR/Tcells_pseudo-bulk_DAR.R
```

## DESeq2 parameters (`apply_DESeq2_test_seurat()`)

| Parameter | Default | Meaning |
|---|---|---|
| `exp.thresh` | 0.05 | Pre-filter: peak must be accessible in ≥5% of cells in foreground OR background |
| `num.splits` | 10 (5 for small Lung) | Each population randomly split into N pseudo-bulk pools |
| `padj.cutoff` | 0.05 | BH FDR threshold for "significant" |
| `min.cells.per.group` | 80 (Kidney/Aorta/Tcells); 50 (Lung) | Excluded cell types below threshold |

A peak is called a DAR if `padj < 0.05` AND `|log2FC| > 0.5`. Opening = log2FC > 0.5; closing = log2FC < -0.5.

## Outputs per cell type × contrast

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
