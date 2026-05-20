# Stage 05 — Re-quantify per-tissue Seurat objects on the universal peak set

Each tissue's annotated Seurat object (from stage 01) is re-quantified against the universal peak set (stage 04 output) using `FeatureMatrix()`. The output is a per-tissue Seurat object with a new `peaks_universal` assay — this is the canonical input to stage 07 (DAR).

## Why?

Different tissues use different peak conventions in their source Seurat objects (Kidney: `chrA-B-C`; Lung: `chrA:B-C`; Aorta: rLSI peaks; T cells: ATAC-original peaks). Re-quantifying everyone against the same universal peak set is required for cross-tissue comparison.

## Run

```bash
# Each tissue is independent and can run in parallel.
sbatch 05_universal_assays/kidney/SeuratObject.slurm
sbatch 05_universal_assays/aorta/SeuratObject.slurm
sbatch 05_universal_assays/tcells/SeuratObject.slurm

# Lung is two-step (the second step is needed because Lung emerged with +14 features
# over the Kidney/Aorta/Tcells feature count after re-quantification):
sbatch 05_universal_assays/lung/SeuratObject.slurm           # → lung_universal.rds
Rscript 05_universal_assays/lung/Diagnose_lung.R             # → lung_universal_new_pruned.rds
```

Each SLURM script invokes its tissue's `SeuratObject.R` with arguments:
- `--obj` — input Seurat .rds (from stage 01)
- `--up` — universal peak BED.gz (stage 04 output)
- `--out` — output .rds
- `--frag_tpl` — fragment-file path template, with `{id}` placeholder filled per-sample
- `--sample_key` — metadata column that maps cells to fragment files

## Tissue-specific barcode reconciliation

Each tissue's `SeuratObject.R` contains a tissue-specific `clean_barcodes()` / `strip_by_id()` function to convert Seurat cell names back into the raw barcode format the fragment file expects:

| Tissue | Barcode quirk | Reconciliation |
|---|---|---|
| Kidney | `SRR..._Kidney_atac_AAACG...-1` | `strip_by_id()` strips SRR prefix |
| Lung | `Control_F2_CELL771_N1` | regex `^.*?(CELL.*)` extracts the BGI raw barcode |
| Aorta | `AAACG...-1_1` or `_2` (Seurat merge suffix) | strip trailing `_N` |
| T cells | mixed: `Tcells_AAACG...-1` or raw | preferentially strip `Tcells_` then trailing `_N` |

## Outputs (canonical DAR substrates)

| Tissue | Output | Notes |
|---|---|---|
| Kidney | `${DATA_ROOT}/Kidney/kidney_merged_universal.rds` | Read by `07_DAR/Kidney_pseudo-bulk_DAR.R` |
| Lung | `${DATA_ROOT}/Lung/lung_universal_new_pruned.rds` | Pruned to match Kidney's 667,459 feature count via `Diagnose_lung.R`. Read by `07_DAR/Lung_specific_DAR_pseudo-bulk.R` |
| Aorta | `${DATA_ROOT}/Aorta/Aorta_integrated_universal.rds` | Read by `07_DAR/Aorta_pseudo-bulk_DAR.R` |
| T cells | `${DATA_ROOT}/Tcells/tcells_universal.rds` | Read by `07_DAR/Tcells_pseudo-bulk_DAR.R` |

All four objects share the same `peaks_universal` assay with the same peak coordinates (after Lung pruning), enabling cross-tissue comparison.

## Why Lung needs `Diagnose_lung.R` (the extra pruning step)

When `SeuratObject.R` re-quantifies Lung against the universal peak set, it occasionally produces +14 extra peak features beyond the Kidney/Aorta/Tcells count (667,473 vs 667,459). The reason is an artefact of the chromosome-set difference between the tissues' source fragment files. `Diagnose_lung.R` prunes Lung's feature space to exactly match Kidney's reference so that all four objects can be merged or compared peak-by-peak in downstream stages.

The output `lung_universal_new_pruned.rds` is what `07_DAR/Lung_specific_DAR_pseudo-bulk.R` reads.
