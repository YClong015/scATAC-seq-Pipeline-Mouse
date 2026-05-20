# Stage 06 — 4-tissue Harmony integration (Aim 1)

Merge the four per-tissue `*_universal_v5.rds` Seurat objects into one cross-tissue atlas, then run Harmony correction by `Tissue` to produce the integrated UMAP.

## Run order

```bash
# (a) Diagnostic — ensure all four objects share the same feature space.
#     Required because Lung pre-pruned object may have +14 features over Kidney.
#     This script prunes Lung to match Kidney's 667,459 features.
Rscript 06_integration/Diagnose_lung.R

# (b) Build merged + Harmony-integrated object.
sbatch 06_integration/Merge_and_integrate.slurm

# (c) Re-annotate clusters using the diagnostic plots.
#     Two-phase script: set PHASE <- 1, run, review plots, fill cluster_map, then PHASE <- 2.
Rscript 06_integration/Integration_reannotate.R     # phase 1: diagnostic plots
# (Manually fill cluster_map in the script using cluster_annotation_guide.csv)
Rscript 06_integration/Integration_reannotate.R     # phase 2: apply annotation

# (d) (Optional) Mount fragments + compute GeneActivity for downstream marker discovery.
Rscript 06_integration/UMAP_annotation_with_fragment_linking.R
```

## Outputs

| File | Contents |
|---|---|
| `Integrated/All_Tissues_Integrated.rds` | Merged + Harmony-corrected; `umap.lsi` and `umap.harmony` reductions present |
| `Integrated/All_Tissues_Integrated_Annotated.rds` | With `cell_type_annotated` metadata column (output of phase 2) |
| `Integrated/All_Tissues_Integrated_Annotated_Clean_for_DAR.rds` | Suspicious / OR-high clusters removed (used by `07_DAR/Aorta_Lung_Tcells_pseudo-bulk_DAR.R`) |
| `QC_figures/Fig_Integration_UMAP.{pdf,png}` | Thesis Fig 2 / supplementary |
| `QC_figures/annotation_diagnostics/cluster_annotation_guide.csv` | Used as input to manually fill cluster_map |

## Key parameters

- Universal peak assay: `peaks_universal`
- TF-IDF + SVD: dims 2:30 (dim 1 excluded — sequencing depth)
- Harmony: `group.by.vars = "Tissue"`, `dims.use = 2:30`
- UMAP after Harmony: dims 1:29
- Clustering: Leiden algorithm 3, resolution 0.5
- Quality filter: cells with `colSums(peaks_universal counts) <= 100` removed
