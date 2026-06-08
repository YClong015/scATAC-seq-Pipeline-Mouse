# Stage 08 - Aging chromatin comparison (Aim 3)

Re-call aging DARs from the Lu et al. (2026) mouse aging chromatin atlas using **the same pseudo-bulk DESeq2 pipeline** as for disease, intersect against the disease DARs (stage 07 output), and quantify overlap via Fisher's exact test + Pearson r.

## Aging-DAR call

Aging DARs are re-called with the **Aged (21 mo) vs Adult (5 mo)** "pure aging" contrast (output `aging_DARs_5movs21mo/`), which excludes the early-life maturation component. A `sex` covariate is added when both sexes are balanced (>=2 per group), matching Lu et al. 2026's sex-stratified approach.

## Run

```bash
# 0. Inspect Lu 2026 h5ad structure (cell types, ages, samples)
sbatch 08_aging_comparison/01_explore_h5ad.slurm

# 1. Re-call aging DARs (Aged 21mo vs Adult 5mo), one job per tissue
sbatch 08_aging_comparison/02_aging_DAR/02a_pseudobulk_Kidney_5movs21mo.slurm
sbatch 08_aging_comparison/02_aging_DAR/02b_pseudobulk_Lung_5movs21mo.slurm

# 2. Disease-vs-aging overlap, one job per tissue
sbatch 08_aging_comparison/03_overlap/03a_overlap_Kidney_5movs21mo.slurm   # -> results_Kidney_5movs21mo/
sbatch 08_aging_comparison/03_overlap/03b_overlap_and_plots_Lung_5movs21mo.slurm  # -> results_Lung_5movs21mo/

# 3. Combined thesis figures live in 09_figures/ (they read the per-tissue 5mo results above):
#      09_figures/Fig13_aging_DAR_counts/Fig13_replot.R            # Fig 13 - ageing DAR counts
#      09_figures/Fig14_Fisher_OR/Fig14_Fisher_OR.R               # Fig 14 - Fisher OR (A Kidney, B Lung)
#      09_figures/Fig15_disease_aging_scatter/Fig15_disease_aging_scatter.R  # Fig 15 - log2FC scatter
```

## Methodology sensitivity tests (`methodology_tests/`)

Two short scripts test alternative pseudo-bulk DAR strategies on the same kidney/lung cell types:
- `test_02_sex_only.slurm` - V0 (baseline) vs V1 (add sex covariate)
- `test_02_variants.slurm` - V0-V3 grid (sex covariate x CPM-top-25% pre-filter)

Both produce a CSV comparing DAR counts across the 2-4 variants. Used during method development; results documented in `${DATA_ROOT}/DAR/DAR_science_comparison/test_02_*.csv`.

## Statistical tests

**Fisher's exact test** (one-sided, "greater"):
- 2x2 table: disease-DAR x aging-DAR membership
- `N` = total peaks tested by DESeq2 in the aging atlas for that cell type
- Both disease and aging significance use `padj < 0.05`
- Output: OR (odds ratio) + p-value + significance star

**Pearson r**:
- Computed on the disease-significant peaks that overlap aging peaks
- Restricted to aging-significant (`padj < 0.05`) subset for robust correlation
- Concordance metric: % of disease-opening peaks where aging log2FC > 0; same for closing

## Cell-type mappings (disease -> Lu 2026 `Main_cell_type`)

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
| (others) | See `02a_pseudobulk_Kidney_5movs21mo.slurm` for full map |

## Outputs

```
# Per-tissue overlap (03a / 03b)
${DATA_ROOT}/DAR/DAR_science_comparison/results_{Kidney,Lung}_5movs21mo/
  overlap_stats.csv      # Fisher OR + p per cell type   (read by Fig 14)
  scatter_r_stats.csv    # Pearson r per cell type
  paired_lfc/            # disease + aging log2FC per overlapping peak  (read by Fig 15)
  figures/               # per-tissue diagnostic plots

# Combined thesis figures (09_figures/Fig13-15, read the per-tissue results above)
${DATA_ROOT}/DAR/DAR_science_comparison/Fig14_output/Fig14_Fisher_OR_{Kidney,Lung}.{pdf,png}
${DATA_ROOT}/DAR/DAR_science_comparison/Fig15_output/Fig15_{solo_Kidney_TAL,...,combined}.{pdf,png}
```
