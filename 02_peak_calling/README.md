# Stage 02 — Per-cell-type peak calling (MACS2 via pycisTopic)

After `SplitFragments()` produces one BED per cell type, MACS2 is run on each BED to call narrow peaks. We use the pycisTopic `peak_calling()` wrapper for consistent BEDPE-friendly defaults.

## MACS2 parameters

| Parameter | Value | Rationale |
|---|---|---|
| `shift` | -73 | Centre Tn5 fragments on accessibility site (Buenrostro 2015) |
| `ext_size` | 146 | Nucleosome footprint |
| `keep_dup` | all | scATAC pre-dedup'd by Cell Ranger |
| `q_value` | 0.05 | Default narrowPeak FDR |
| `genome_size` | 'mm' or '2.7e9' | mouse |
| `input_format` | BEDPE (kidney) or BED (lung/aorta/tcells) | depends on SplitFragments format |

## Run order

```bash
# Kidney — 12 cell types, BEDPE input
sbatch 02_peak_calling/kidney/peak_kidney.slurm   # SLURM array 1-12

# OR per-sample 9-SRR variants (one peak set per SRR, used for sample-level QC):
for srr in SRR27367330 SRR27367331 SRR27367332 SRR27367340 SRR27367344 SRR27367346 SRR27367347 SRR27367349 SRR27367351; do
  python 02_peak_calling/kidney/per_sample/$srr/peak_calling.py
done

# Lung — 15 cell types
sbatch 02_peak_calling/lung/peak_calling_for_lung.slurm   # SLURM array 1-15

# Aorta — 6 cell types
sbatch 02_peak_calling/aorta/peak_calling_for_aortic.slurm   # SLURM array 1-6

# T cells — 9 cell types
sbatch 02_peak_calling/tcells/peak_Tcell.slurm   # SLURM array 1-9
```

Each script reads a `DataList.txt` listing the per-cell-type BED filenames and dispatches one MACS2 job per file via the SLURM `--array` directive.

## Outputs

```
${DATA_ROOT}/{tissue}/peaks/{cell_type}/sample_peak/{cell_type}_peaks.narrowPeak
```

These narrowPeak files are the inputs to stage 03 (per-tissue consensus merging).
