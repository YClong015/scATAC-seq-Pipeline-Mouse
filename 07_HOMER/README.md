# Stage 07 - HOMER motif enrichment

One self-contained SLURM array job per tissue. Each task generates the opening / closing / stable / NS BED files for one unit, then runs all HOMER `findMotifsGenome.pl` comparisons.

## Backgrounds and comparisons

For each unit (cell type x contrast), HOMER is run against two backgrounds plus the closing-vs-opening pair:

- **NS background** (non-significant peaks, `padj >= 0.05`): everything DESeq2 tested that was not a DAR. Controls for accessibility-level bias.
- **Stable background** (truly invariant peaks, `padj > 0.9` and `abs(log2FC) < 0.05`): controls for accessibility and cell-type-specific peak architecture.
- **Closing vs opening**: closing DARs against opening DARs, to recover identity TFs (Fig 12).

Including the reciprocal tests, this is **10 HOMER runs per unit** (4 stable + 4 NS + 2 closing-vs-opening).

## Pipeline

Each tissue is one self-contained array job that builds its own BEDs and then runs HOMER, so there is no separate BED-prep step.

```bash
sbatch 07_HOMER/run_homer/Kidney_HOMER_full.slurm   # v5 11-type, array per cell type (Day42 vs Sham)
sbatch 07_HOMER/run_homer/Lung_HOMER_full.slurm     # array per cell type (Case vs Control)
sbatch 07_HOMER/run_homer/Aorta_HOMER_full.slurm    # array per cell type (Challenge vs Control)
sbatch 07_HOMER/run_homer/Tcells_HOMER_full.slurm   # pooled T cells, array per contrast (vs Young_control)
```

The cell-type list and contrast are set at the top of each script. T cells are a single pooled population, so that job arrays over disease contrasts instead of cell types.

## Outputs

For each unit (cell type x contrast):

```
${DAR_DIR}/
  DAR_BED_HOMER/{ct}__{contrast}__{opening,closing,stable,NS}.bed
  HOMER_stable_bg/{ct}__{contrast}__005__{direction}_vs_stable/knownResults.txt
  HOMER_NS_bg/    {ct}__{contrast}__005__{direction}_vs_NS/knownResults.txt
  HOMER_closing_vs_opening/{ct}__{contrast}__{closing_vs_opening,opening_vs_closing}/knownResults.txt
```

`${DAR_DIR}` is `.../DAR_pseudobulk_Kidney_v5_DESeq2` for Kidney and `.../DAR_pseudobulk_DESeq2/DAR_pseudobulk_{tissue}_DESeq2` for the others.

The `knownResults.txt` files feed the plotting scripts in `09_figures/Fig10_HOMER_heatmap/` (NS + stable), `09_figures/Fig11_opening_dotplot/` (NS opening), and `09_figures/Fig12_identity_TF/` (closing vs opening).

## HOMER tuning

| Setting | Value | Why |
|---|---|---|
| Genome | `mm10` | Mouse |
| `-size given` | yes | Use BED coordinates as-is (peaks already fixed-width 501 bp) |
| `-nomotif` | yes | Skip de-novo discovery; known-motif scan only |
| `-mask` | yes | Repeat-mask before scanning |
| `MIN_PEAKS` | 200 | Below this HOMER's hypergeometric test becomes unstable |

## Manual HOMER install required

```bash
mkdir -p /scratch/user/${USER}/homer && cd /scratch/user/${USER}/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
perl configureHomer.pl -install mm10
export HOMER_HOME=/scratch/user/${USER}/homer
```

Edit `HOMER_HOME` at the top of each `07_HOMER/run_homer/*_HOMER_full.slurm` to point at your install.
