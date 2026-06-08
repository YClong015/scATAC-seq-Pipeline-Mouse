# Stage 09 - Figures

Each subdirectory corresponds to one thesis figure. Figure numbers follow
ThesisDraft_V4 (results figures start at Fig 4).

| Figure | Subdir | Script |
|---|---|---|
| Fig 4 (QC violin metrics) | `Fig4_QC/` | `QC_plots.R` (cache-or-load: reads cached metadata CSVs if present, else loads RDS) |
| Fig 5 (4-tissue UMAP annotation) | `Fig5_UMAP/` | `Fig5_UMAP_annotation.R`, `Fig5_Tcells_annotation.R`; integrated UMAP via `Fig5_integration_UMAP.R` |
| Fig 6 (CoveragePlot at marker genes) | `Fig6_CoveragePlot/` | `Fig6_link_fragments.R` (one-time setup) + `Fig6_CoveragePlot.R` |
| Fig 7 (peak counts pre-consensus) | `Fig7_peak_counts/` | `count_peaks_per_celltype.slurm` + `plot_peaks_per_celltype{,_per_tissue}.R` |
| Fig 8 (universal peak set composition) | `Fig8_universal_peaks/` | `Fig8_universal_peaks.R` (bar chart + supplementary UpSet) |
| Fig 9 (DAR burden barplot) | `Fig9_DAR_burden/` | `Fig9_DAR_burden.R` |
| Fig 10 (integrated HOMER heatmap) | `Fig10_HOMER_heatmap/` | `Fig10_integrate_heatmap.R` + `assemble_Fig10.R` (PNG composition) |
| Fig 11 (4-tissue opening TF motifs, top 10) | `Fig11_opening_dotplot/` | `Fig11_opening_dotplot_4tissue.R` -> combined + 12A-D per-tissue solo panels |
| Fig 12 (cell-type identity TF, closing DARs) | `Fig12_identity_TF/` | `Fig12_cell_identity_barplot.R` (v7 filter, Kidney 11-type) |
| Fig 13 (ageing DARs) | `Fig13_aging_DAR_counts/` | `Fig13_replot.R` |
| Fig 14 (Fisher OR, disease vs ageing) | `Fig14_Fisher_OR/` | `Fig14_Fisher_OR.R` -> `Fig14_Fisher_OR_{Kidney,Lung}.{pdf,png}` |
| Fig 15 (disease-ageing concordance scatter) | `Fig15_disease_aging_scatter/` | `Fig15_disease_aging_scatter.R` -> `Fig15_combined.{pdf,png}` |

## Run order (after stages 07-09 complete)

```bash
# Fig 4 - QC violin metrics
Rscript 09_figures/Fig4_QC/QC_plots.R

# Fig 5 - per-tissue UMAPs and integrated UMAP
Rscript 09_figures/Fig5_UMAP/Fig5_UMAP_annotation.R
Rscript 09_figures/Fig5_UMAP/Fig5_Tcells_annotation.R
Rscript 09_figures/Fig5_UMAP/Fig5_integration_UMAP.R

# Fig 6 - Coverage plots at marker genes
Rscript 09_figures/Fig6_CoveragePlot/Fig6_link_fragments.R    # one-time prep
Rscript 09_figures/Fig6_CoveragePlot/Fig6_CoveragePlot.R

# Fig 7 - peak counts per cell type (pre-consensus)
sbatch 09_figures/Fig7_peak_counts/count_peaks_per_celltype.slurm
Rscript 09_figures/Fig7_peak_counts/plot_peaks_per_celltype.R
Rscript 09_figures/Fig7_peak_counts/plot_peaks_per_celltype_per_tissue.R

# Fig 8 - universal peak set composition
Rscript 09_figures/Fig8_universal_peaks/Fig8_universal_peaks.R

# Fig 9 - DAR burden barplot
Rscript 09_figures/Fig9_DAR_burden/Fig9_DAR_burden.R

# Fig 10 - integrated 4-tissue HOMER heatmap
Rscript 09_figures/Fig10_HOMER_heatmap/Fig10_integrate_heatmap.R
Rscript 09_figures/Fig10_HOMER_heatmap/assemble_Fig10.R            # compose NS + Stable panels

# Fig 11 + 12A-D - 4-tissue opening TF motif dotplot (combined + per-tissue solo panels)
Rscript 09_figures/Fig11_opening_dotplot/Fig11_opening_dotplot_4tissue.R

# Fig 12 - cell-type identity TF baseline (v7-filtered barplot, Kidney 11-type)
Rscript 09_figures/Fig12_identity_TF/Fig12_cell_identity_barplot.R

# Fig 13-16 - see 08_aging_comparison/ (figures generated as part of the overlap pipeline)
```

Figures are written under `${DATA_ROOT}/QC_figures/` (Fig 4-9), `${DATA_ROOT}/DAR/Fig9_DAR_burden/` (Fig 9), `${DATA_ROOT}/DAR/Integrated_HOMER_Heatmaps/` (Fig 10), `${DATA_ROOT}/DAR/Fig11_output/` (Fig 11 + 12A-D), `${DATA_ROOT}/DAR/DAR_closing_vs_opening/figures_filtered_v7_kidney11type/` (Fig 12), and `${DATA_ROOT}/DAR/DAR_science_comparison/` (Fig 13-16: per-tissue `results_*_5movs21mo/` + combined `Fig14_output/`, `Fig15_output/`).
