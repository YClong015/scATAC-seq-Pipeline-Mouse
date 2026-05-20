# Stage 09 — Aging chromatin comparison (Aim 3)

Re-call aging DARs from the Lu et al. (2026) mouse aging chromatin atlas using **the same pseudo-bulk DESeq2 pipeline** as for disease, intersect against the disease DARs (stage 07 output), and quantify overlap via Fisher's exact test + Pearson r.

## Aging-DAR call variants

| Variant | Contrast | Output dir | When to use |
|---|---|---|---|
| Standard | `Aged (21 mo) vs Young (1 mo)` | `aging_DARs/` | Default — broadest DAR coverage |
| 5mo-vs-21mo | `Aged (21 mo) vs Adult (5 mo)` | `aging_DARs_5movs21mo/` | "Pure aging" contrast (excludes maturation component) |

The 5mo-vs-21mo variant also adds a `sex` covariate when both sexes are balanced (≥2 per group), matching Lu et al. 2026 sex-stratified approach.

## Run

```bash
# 0. Inspect Lu 2026 h5ad structure (cell types, ages, samples)
sbatch 09_aging_comparison/01_explore_h5ad.sh

# 1. Re-call aging DARs (one job per tissue × variant)
sbatch 09_aging_comparison/02_aging_DAR/02_pseudobulk_DESeq2.sh        # original combined (kidney+lung)
# OR finer-grained variants:
sbatch 09_aging_comparison/02_aging_DAR/02a_pseudobulk_Kidney.sh
sbatch 09_aging_comparison/02_aging_DAR/02a_pseudobulk_Kidney_5movs21mo.sh
sbatch 09_aging_comparison/02_aging_DAR/02b_pseudobulk_Lung.sh
sbatch 09_aging_comparison/02_aging_DAR/02b_pseudobulk_Lung_5movs21mo.sh

# 2. Disease-vs-aging overlap (one variant per script)
sbatch 09_aging_comparison/03_overlap/03_overlap_and_plots.sh                       # combined kidney+lung
sbatch 09_aging_comparison/03_overlap/03a_overlap_Kidney_oldobj.sh                  # Kidney isolated
sbatch 09_aging_comparison/03_overlap/03a_overlap_Kidney_oldobj_5movs21mo.sh        # Kidney pure-aging
sbatch 09_aging_comparison/03_overlap/03b_overlap_and_plots_Lung.sh                 # Lung isolated
sbatch 09_aging_comparison/03_overlap/03b_overlap_and_plots_Lung_5movs21mo.sh       # Lung pure-aging

# 3. Region classification (disease-specific / shared / aging-specific)
sbatch 09_aging_comparison/04_region_classification.sh

# 4. (Optional) Final replots — used for Fig 13 and Aim 3 summary
Rscript 09_aging_comparison/replots/Fig13_replot.R
Rscript 09_aging_comparison/replots/replot_Kidney.R   # Aim 3 Kidney-only summary with OR≥1 filter
```

## Methodology sensitivity tests (`methodology_tests/`)

Two short scripts test alternative pseudo-bulk DAR strategies on the same kidney/lung cell types:
- `test_02_sex_only.sh` — V0 (baseline) vs V1 (add sex covariate)
- `test_02_variants.sh` — V0–V3 grid (sex covariate × CPM-top-25% pre-filter)

Both produce a CSV comparing DAR counts across the 2–4 variants. Used during method development; results documented in `${DATA_ROOT}/DAR/DAR_science_comparison/test_02_*.csv`.

## Statistical tests

**Fisher's exact test** (one-sided, "greater"):
- 2×2 table: disease-DAR × aging-DAR membership
- `N` = total peaks tested by DESeq2 in the aging atlas for that cell type
- Aging relaxed to `padj < 0.2` to maintain test power
- Output: OR (odds ratio) + p-value + significance star

**Pearson r**:
- Computed on the disease-significant peaks that overlap aging peaks
- Restricted to aging-significant (`padj < 0.05`) subset for robust correlation
- Concordance metric: % of disease-opening peaks where aging log2FC > 0; same for closing

## Cell-type mappings (disease → Lu 2026 `Main_cell_type`)

| Our label | Lu 2026 Main_cell_type |
|---|---|
| Kidney PT | Proximal tubule cells, Proximal tubule cells_S3T2 |
| Kidney TAL | Thick ascending limb of LOH cells |
| Kidney DCT | Distal convoluted tubule cells |
| Kidney PC | Principal cells |
| Kidney Macrophages | Myeloid cells_Macrophages |
| Lung AT2 | Type II alveolar epithelial cells |
| Lung Mac-alv | Myeloid cells_Alveolar macrophages |
| Lung Fib | Fibroblasts |
| (others) | See `02_pseudobulk_DESeq2.sh` for full map |

## Outputs

Per variant, each script produces:

```
${DATA_ROOT}/DAR/DAR_science_comparison/{variant}/
  results/overlap_stats.csv             # Fisher OR + p per cell type
  results/scatter_r_stats.csv           # Pearson r per cell type
  results/figures/DAR_overlap_enrichment.{pdf,png}
  results/figures/DAR_overlap_pvalue.{pdf,png}
  results/figures/DAR_scatter_combined.{pdf,png}
  results/figures/Aim3_DAR_overlap_summary.{pdf,png}
  classified_regions/{tissue}_{ct}_{disease|aging}_specific_{open|close}.bed
  classified_regions/{tissue}_{ct}_shared_{open|close}.bed
```

The `_shared_open.bed` / `_shared_close.bed` files are the candidate **aging-disease regulatory hubs** released as the supplementary peak track in the thesis Appendix B.
