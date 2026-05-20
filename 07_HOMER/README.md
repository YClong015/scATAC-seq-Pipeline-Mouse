# Stage 07 — HOMER motif enrichment

Run HOMER `findMotifsGenome.pl` on each cell type × direction × background combination to identify enriched TF motifs.

## Two backgrounds per direction

For each (cell type, contrast, direction), HOMER is run with two backgrounds:

- **NS background** (non-significant peaks): `padj >= 0.05` — everything DESeq2 tested that wasn't a DAR. Controls for accessibility-level bias.
- **Stable background** (truly invariant peaks): `padj > 0.9` AND `abs(log2FC) < 0.05`. Controls for both accessibility AND cell-type-specific peak architecture.

Plus the *reciprocal* tests — background vs opening/closing — giving **8 HOMER runs per cell type per contrast**.

## Pipeline

```bash
# 1. Generate the per-cell-type, per-direction BED files (opening / closing / stable / NS)
Rscript 07_HOMER/prepare_bed/mk_NS_Stable_files_4tissues.R

# 2. Run the canonical 4-tissue HOMER pipeline (stable + NS backgrounds, 8 comparisons per CT)
sbatch 07_HOMER/run_homer/run_4tissues_NS-stable_HOMER.sh

# 3. (Optional) Kidney-specific rerun with focused contrasts:
sbatch 07_HOMER/run_homer/kidney_rerun/run_kidney_homer_stable_bg.sh
sbatch 07_HOMER/run_homer/kidney_rerun/run_kidney_homer_NS_Day42_only.sh

# 4. (Optional) Closing-vs-opening identity-TF baseline (used for Fig 11)
sbatch 07_HOMER/run_homer/run_HOMER_closing_vs_opening.sh

# 5. (Optional) Resume incomplete Lung NS run (if it timed out)
bash 07_HOMER/run_homer/resume_lung_homer_NS.sh
```

The `per_tissue_initial/` directory contains the original whole-DAR HOMER runs (no background separation) — useful for the per-tissue dotplot figures (Fig 13-16) but **not** for the integrated heatmap (Fig 10), which requires the NS/Stable BG outputs.

## Outputs

For each cell type × contrast × direction × background:

```
${DATA_ROOT}/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_{tissue}_DESeq2/
  HOMER_stable_bg/{cell_type}__{contrast}__005__{direction}_vs_stable/
    knownResults.txt
    knownResults.html
    homer.log
  HOMER_NS_bg/    {cell_type}__{contrast}__005__{direction}_vs_NS/
    (same structure)
```

The `knownResults.txt` files feed every plotting script in `09_figures/Fig10_HOMER_heatmap/`, `09_figures/Fig11_identity_TF/`, `09_figures/Fig12_17_TF_motif_dotplots/`, and `09_figures/per_tissue_HOMER_plots/`.

## HOMER tuning

| Setting | Value | Why |
|---|---|---|
| Genome | `mm10` | Mouse |
| `-size given` | yes | Use BED coords as-is (peaks already fixed-width 501 bp) |
| `-nomotif` | yes | Skip de-novo motif discovery (slow; only known-motif scan needed for thesis) |
| `-mask` | yes (per_tissue_initial only) | Repeat-mask before scanning |
| Min peaks | 20 | Below this HOMER's hypergeometric test becomes unstable |
| Min DAR per direction | 20 | Excluded cell-type×contrast combos below threshold |

## Manual HOMER install required

```bash
mkdir -p /scratch/user/${USER}/homer && cd /scratch/user/${USER}/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
perl configureHomer.pl -install mm10
export HOMER_HOME=/scratch/user/${USER}/homer
```

Edit `HOMER_HOME` at the top of each `07_HOMER/run_homer/*.sh` to point at your install.
