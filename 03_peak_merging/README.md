# Stage 03 — Per-tissue consensus peak merging

Merge per-cell-type MACS2 narrowPeak calls into one consensus peak set per tissue using pycisTopic's `get_consensus_peaks()`. Each peak is extended ±250 bp around its summit before iterative greedy merging (Bravo González-Blas et al., 2023).

## Run

```bash
# Build per-tissue consensus for all 4 tissues in one SLURM job
sbatch 03_peak_merging/merge_consensus_tissues.slurm
```

This runs `consensus_mm10.py` once per tissue, reading the narrowPeak manifest at `${DATA_ROOT}/{tissue}/peak_merging/narrow_peak_paths.txt` and writing:

```
${DATA_ROOT}/{tissue}/peak_merging/consensus_peak_calling/consensus_regions.bed.gz
```

## Outputs (approximate peak counts)

| Tissue | Peaks (consensus) |
|---|---|
| Kidney | 581,481 |
| Lung | 245,604 |
| Aorta | 230,473 |
| T cells | 145,370 |

These per-tissue consensus BEDs feed into stage 04 to build the cross-tissue universal peak set.

## QC: count peaks per cell type

```bash
bash   09_figures/peak_counts/count_peaks_per_celltype.sh
Rscript 09_figures/peak_counts/plot_peaks_per_celltype.R
```
Produces `peaks_per_celltype.csv` and a faceted bar chart for Figure 5 / Supplementary.
