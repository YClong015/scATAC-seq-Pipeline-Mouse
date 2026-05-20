# Stage 10 — Figures

All 21 main thesis figures plus supplementary figures. Each subdirectory corresponds to one figure or set of related figures.

| Figure | Subdir | Script |
|---|---|---|
| Fig 3 (mice + cells per group) | `Fig3_4_QC/` | `QC_plots.sh`, `QC_plots_debug.R` (debug version has the cell-composition-per-sample panels) |
| Fig 4 (QC violin metrics) | `Fig3_4_QC/` | same |
| Fig 5 (universal peak bar chart) | `Fig5_6_universal_peaks/` | `Fig4_universal_peaks_v2.R` (reads `upset_input.csv` from stage 04) |
| Fig 6 (UpSet plot) | `Fig5_6_universal_peaks/` | same |
| Fig 7 (4-tissue UMAPs) | `Fig7_UMAP/` | `Fig5_UMAP_annotation.R`, `Fig5_Tcells_annotation.R`; integrated UMAP via `Fig_integration_UMAP.R` |
| Fig 8 (CoveragePlot at marker genes) | `Fig8_CoveragePlot/` | `link_fragments_and_save.R` (one-time setup) + `Fig_CoveragePlot.R` |
| Fig 9 (DAR counts diverging tile) | `Fig9_DAR_counts/` | `Fig9_DAR_counts.R` |
| Fig 10 (integrated HOMER heatmap) | `Fig10_HOMER_heatmap/` | `4tissues_integrate_heatmap.R` + `assemble_Fig7.R` (PNG composition) |
| Fig 11 (cell identity TF baseline) | `Fig11_identity_TF/` | `cell_identity_heatmap.sh` (CM-paper-style) OR `cell_identity_barplot.R` (the canonical v6 filter); `assemble_Fig6.R` to compose 6 panels |
| Fig 12 (shared opening TF motifs) | `Fig12_17_TF_motif_dotplots/` | `Fig10_11_16_dotplots.R` → `Fig10_Opening_shared_dotplot.png` |
| Fig 13 (Kidney TF motifs) | `per_tissue_HOMER_plots/` | `homer_kidney_plot.R` OR `Kidney_amalgamation.R` (true amalgamation across cell types) |
| Fig 14 (Lung TF motifs) | `per_tissue_HOMER_plots/` | `homer_lung_plot.R` |
| Fig 15 (Aorta TF motifs) | `per_tissue_HOMER_plots/` | `homer_aorta_plot.R` |
| Fig 16 (T cell TF motifs) | `per_tissue_HOMER_plots/` | `homer_tcells_plot.R` |
| Fig 17 (tissue-specific opening motifs) | `Fig12_17_TF_motif_dotplots/` | `Fig10_11_16_dotplots.R` → `Fig16_Opening_specific_dotplot.png` |
| Fig 18 (aging DAR counts) | (in `09_aging_comparison/replots/Fig13_replot.R`) | renamed Fig 18 in thesis |
| Fig 19 (Aim 3 OR + r summary) | `09_aging_comparison/03_overlap/03_overlap_and_plots.sh` → `Aim3_DAR_overlap_summary.{pdf,png}` | |
| Fig 20 (disease vs aging scatter) | `09_aging_comparison/03_overlap/03_overlap_and_plots.sh` → `DAR_scatter_combined.{pdf,png}` | |
| Fig 21 (region classification) | `09_aging_comparison/04_region_classification.sh` → `region_stacked_counts.{pdf,png}` etc. | |
| Supp Fig 7 (AllTissues comprehensive dotplot) | `Fig12_17_TF_motif_dotplots/` | `Fig10_11_16_dotplots.R` |
| Peak counts per cell type | `peak_counts/` | `count_peaks_per_celltype.sh` + `plot_peaks_per_celltype{,_per_tissue}.R` |

## Run order (after stages 07-09 complete)

```bash
# Fig 3 / 4 — QC composition + per-cell metrics
bash    10_figures/Fig3_4_QC/QC_plots.sh

# Fig 5 / 6 — universal peak set composition + UpSet
Rscript 10_figures/Fig5_6_universal_peaks/Fig4_universal_peaks_v2.R

# Fig 7 — per-tissue UMAPs and integrated UMAP
Rscript 10_figures/Fig7_UMAP/Fig5_UMAP_annotation.R
Rscript 10_figures/Fig7_UMAP/Fig5_Tcells_annotation.R
Rscript 10_figures/Fig7_UMAP/Fig_integration_UMAP.R

# Fig 8 — Coverage plots at marker genes
Rscript 10_figures/Fig8_CoveragePlot/link_fragments_and_save.R    # one-time prep
Rscript 10_figures/Fig8_CoveragePlot/Fig_CoveragePlot.R

# Fig 9 — DAR counts diverging tile
Rscript 10_figures/Fig9_DAR_counts/Fig9_DAR_counts.R

# Fig 10 — integrated 4-tissue HOMER heatmap
Rscript 10_figures/Fig10_HOMER_heatmap/4tissues_integrate_heatmap.R
Rscript 10_figures/Fig10_HOMER_heatmap/assemble_Fig7.R            # compose NS + Stable panels

# Fig 11 — cell-type identity TF baseline (CM-paper style or v6-filtered barplot)
sbatch  10_figures/Fig11_identity_TF/cell_identity_heatmap.sh
Rscript 10_figures/Fig11_identity_TF/cell_identity_barplot.R
Rscript 10_figures/Fig11_identity_TF/assemble_Fig6.R              # compose 6 cell-identity panels

# Fig 12 + 17 + Supp Fig 7 — shared / specific / all-tissues TF motif dotplots
Rscript 10_figures/Fig12_17_TF_motif_dotplots/Fig10_11_16_dotplots.R
sbatch  10_figures/Fig12_17_TF_motif_dotplots/run_TF_motif_dotplot.sh

# Fig 13-16 — per-tissue HOMER bubble plots
Rscript 10_figures/per_tissue_HOMER_plots/homer_kidney_plot.R
Rscript 10_figures/per_tissue_HOMER_plots/homer_lung_plot.R
Rscript 10_figures/per_tissue_HOMER_plots/homer_aorta_plot.R
Rscript 10_figures/per_tissue_HOMER_plots/homer_tcells_plot.R
Rscript 10_figures/per_tissue_HOMER_plots/Kidney_amalgamation.R   # alternative Kidney panel

# Fig 18-21 — see 09_aging_comparison/ (figures generated as part of the overlap pipeline)
```

All figures are written under `${DATA_ROOT}/QC_figures/` (Fig 3-8), `${DATA_ROOT}/DAR/Fig9_DAR_counts/` (Fig 9), `${DATA_ROOT}/DAR/Integrated_HOMER_Heatmaps/` (Fig 10), `${DATA_ROOT}/DAR/DAR_pseudobulk_DESeq2/.../HOMER_Plots/` (Fig 13-16), `${DATA_ROOT}/DAR/Combined_HOMER_Heatmap_Focused/TF_motif_plots/` (Fig 12, 17, Supp Fig 7), and `${DATA_ROOT}/DAR/DAR_science_comparison/results*/figures/` (Fig 19-21).
