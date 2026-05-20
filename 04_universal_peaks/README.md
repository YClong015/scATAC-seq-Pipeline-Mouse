# Stage 04 — Universal peak set (v5)

Combines the four per-tissue consensus peak sets into a single 667,473-peak universal set spanning all four tissues. This is the substrate for all cross-tissue analysis.

## Three-stage construction (per `build_universal_scpm.py`)

1. **Stage 1.** Per-tissue consensus from raw narrowPeak calls (mirrors stage 03 but uses MACS2 narrowPeak directly, with `-log10(qval)` as score).
2. **Stage 2.** Cross-tissue union — narrow each 500 bp Stage-1 peak to its centre ±125 bp before cluster-merging via pyranges, max-score per cluster.
3. **Stage 3.** SCPM filter — peaks with signal-per-cell-per-million ≤ cutoff (default 1.0) are dropped.

## Run

```bash
sbatch 04_universal_peaks/run_build_universal_scpm.sh
```

Required env vars (set in `run_build_universal_scpm.sh`):
- `CHROMSIZES_FILE` — `mm10.chrom.sizes`
- `BLACKLIST` — ENCODE mm10 blacklist v2
- `OUTDIR` — output dir (`${DATA_ROOT}/universal_peaks_v5`)
- `PEAK_HALF_WIDTH=250`, `SCPM_CUTOFF=1.0`

## Outputs

| File | Contents |
|---|---|
| `consensus_regions_v5.bed.gz` | Final 667,473 universal peaks (BED + Score + SCPM) — used by 05_universal_assays |
| `universal_peaks_prescpm.bed.gz` | Pre-SCPM-filter version, useful for sensitivity analysis |
| `upset_input.csv` | Per-peak tissue membership matrix → used by `10_figures/Fig5_6_universal_peaks/Fig4_universal_peaks_v2.R` |
| `Kidney_consensus.bed.gz`, `Lung_consensus.bed.gz`, `Aorta_consensus.bed.gz`, `Tcells_consensus.bed.gz` | Cached per-tissue Stage-1 consensus (reused on subsequent runs) |

## Peak counts (final SCPM > 1)

- Universal: 667,473
- Kidney-only: 292,438
- Pan-tissue (all 4): 1,302
