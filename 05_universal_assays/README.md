# Stage 05 — Re-quantify per-tissue Seurat objects on the universal peak set

Each tissue's original Seurat object is re-quantified against the universal peak set (`consensus_regions_v5.bed.gz`) using `FeatureMatrix()`. The output is a per-tissue `*_universal_v5.rds` with a new `peaks_universal` assay.

## Why?

Different tissues use different peak conventions in their source Seurat objects (Kidney: chrA-B-C; Lung: chrA:B-C; Aorta: rLSI peaks; T cells: ATAC-original peaks). Re-quantifying everyone against the same 667,473-peak universal set is required for cross-tissue comparison.

## Run

```bash
# Each tissue gets its own SLURM job — they're independent and can run in parallel.
sbatch 05_universal_assays/kidney/SeuratObject.slurm
sbatch 05_universal_assays/lung/SeuratObject.slurm
sbatch 05_universal_assays/aorta/SeuratObject.slurm
sbatch 05_universal_assays/tcells/SeuratObject.slurm
```

Each SLURM script invokes its tissue's `SeuratObject.R` with arguments:
- `--obj` — input Seurat .rds
- `--up` — universal peak BED.gz (`consensus_regions_v5.bed.gz`)
- `--out` — output .rds (`*_universal_v5.rds`)
- `--frag_tpl` — fragment-file path template, with `{id}` placeholder filled per-sample
- `--sample_key` — metadata column that maps cells to fragment files (varies per tissue)

## Tissue-specific barcode reconciliation

Each tissue's `SeuratObject.R` contains a tissue-specific `clean_barcodes()` / `strip_by_id()` function to convert Seurat cell names (`SampleID_BARCODE-1` or `SampleID_BARCODE-1_2` after merge) back into the raw barcode format the fragment file expects:

| Tissue | Barcode quirk | Reconciliation |
|---|---|---|
| Kidney | `SRR..._Kidney_atac_AAACG...-1` | `strip_by_id()` strips SRR prefix |
| Lung | `Control_F2_CELL771_N1` | regex `^.*?(CELL.*)` extracts the BGI raw barcode |
| Aorta | `AAACG...-1_1` or `_2` (Seurat merge suffix) | strip trailing `_N` |
| T cells | mixed: `Tcells_AAACG...-1` or raw | preferentially strip `Tcells_` then trailing `_N` |

## Outputs

```
${DATA_ROOT}/Kidney/kidney_universal_v5.rds
${DATA_ROOT}/Lung/lung_universal_v5.rds         ← used by 06_integration
${DATA_ROOT}/Lung/lung_universal_new_pruned.rds ← used by 07_DAR (canonical Lung DAR substrate)
${DATA_ROOT}/Aorta/Aorta_universal_v5.rds
${DATA_ROOT}/Tcells/tcells_universal_v5.rds
```

## ⚠ Lung object version note

There are two lung universal objects in active use:
1. `lung_universal_v5.rds` — built fresh by `05_universal_assays/lung/SeuratObject.R` against `consensus_regions_v5.bed.gz`. Used by **06_integration** (Merge_and_integrate.R).
2. `lung_universal_new_pruned.rds` — an earlier object built against an older Kidney+Aorta+Tcells+Lung universal set, then pruned to 667,459 peaks for feature-count compatibility with Kidney. **This is what 07_DAR/Lung_specific_DAR_pseudo-bulk.R reads** — so the canonical DAR-output-producing pipeline uses the older pruned object.

This dual-object setup is documented for clarity; if you're re-running from scratch, regenerate the pruned object via `06_integration/Diagnose_lung.R` after building both `kidney_universal_v5.rds` and `lung_universal_v5.rds`.
