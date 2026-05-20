#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --job-name=pseudobulk_deseq2
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'

import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings
warnings.filterwarnings('ignore')

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"
OUT_DIR     = f"{SCIENCE_DIR}/aging_DARs"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters ───────────────────────────────────────────────
MIN_CELLS_PER_SAMPLE = 10   # min cells per cell type per sample
MIN_SAMPLES          = 4    # min pseudobulk samples (total Young+Aged)
# No peak filtering: Science paper h5ad is already processed/filtered data;
# all peaks in the matrix are passed directly to DESeq2.

# ── Cell type mapping: our label → Science paper Main_cell_type ──
# Confirmed from explore output; add/remove as needed
kidney_map = {
    "PT":          ["Proximal tubule cells",
                    "Proximal tubule cells_S3T2"],   # S3 segment subtype, merge with PT
    "TAL":         ["Thick ascending limb of LOH cells"],
    "DCT":         ["Distal convoluted tubule cells"],
    "PC":          ["Principal cells"],
    "Macrophages": ["Myeloid cells_Macrophages"],
}
lung_map = {
    "AT2":         ["Type II alveolar epithelial cells"],
    "B":           ["Lymphoid cells_B cells"],
    "Ciliated":    ["Ciliated cells"],
    "EC-vasc":     ["Vascular endothelial cells",
                    "Vascular endothelial cells_Aerocytes"],
    "Eosinophils": ["Myeloid cells_Eosinophils"],
    "Fib":         ["Fibroblasts"],
    "Mac-alv":     ["Myeloid cells_Alveolar macrophages"],
    "Mac-inter":   ["Myeloid cells_Interstitial macrophages"],
    "NK":          ["Lymphoid cells_NK cells"],
    "Pen":         ["Neuroendocrine cells"],  # only 166 cells total, likely skipped
    "SMCs":        ["Smooth muscle cells"],
    "T":           ["Lymphoid cells_T cells"],
}

# ── Helper: pseudobulk one cell type ─────────────────────────
def pseudobulk(adata, ct_label, science_labels):
    """
    Subset to science_labels, filter Young+Aged, aggregate per Sample.
    Returns (counts_df [samples x peaks], meta_df) or None if skipped.
    """
    mask = adata.obs['Main_cell_type'].isin(science_labels)
    ct   = adata[mask].copy()
    # Keep only Young and Aged
    ct   = ct[ct.obs['Age'].isin(['Young', 'Aged'])].copy()

    print(f"    Cells (Young+Aged): {ct.shape[0]}")
    if ct.shape[0] < 50:
        print("    SKIP: too few cells"); return None, None

    pb_counts = {}
    pb_meta   = {}
    for samp in ct.obs['Sample'].unique():
        idx = ct.obs['Sample'] == samp
        n   = idx.sum()
        if n < MIN_CELLS_PER_SAMPLE:
            continue
        counts = np.asarray(ct.X[idx.values].sum(axis=0)).flatten()
        age    = ct.obs.loc[idx, 'Age'].iloc[0]
        pb_counts[samp] = counts
        pb_meta[samp]   = {'Age': age, 'n_cells': n}

    if len(pb_counts) < MIN_SAMPLES:
        print(f"    SKIP: only {len(pb_counts)} pseudo-bulk samples"); return None, None

    counts_df = pd.DataFrame(pb_counts, index=ct.var_names).T.astype(int)
    meta_df   = pd.DataFrame(pb_meta).T

    n_young = (meta_df['Age'] == 'Young').sum()
    n_aged  = (meta_df['Age'] == 'Aged').sum()
    if n_young == 0 or n_aged == 0:
        print("    SKIP: missing one age group"); return None, None

    print(f"    Pseudo-bulk: Young={n_young}, Aged={n_aged}, Peaks={counts_df.shape[1]}")

    # Remove peaks that are zero across ALL samples — DESeq2 cannot test these.
    keep = (counts_df > 0).any(axis=0)
    counts_df = counts_df.loc[:, keep]
    print(f"    Peaks (after removing all-zero): {counts_df.shape[1]}")

    if counts_df.shape[1] < 100:
        print("    SKIP: too few peaks after filtering"); return None, None

    return counts_df, meta_df

# ── Helper: run DESeq2 ────────────────────────────────────────
def run_deseq2(counts_df, meta_df, tissue, ct_label, out_dir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds  import DeseqStats

    meta_df = meta_df.copy()
    meta_df['condition'] = meta_df['Age'].map({'Young': 'Young', 'Aged': 'Aged'})

    try:
        dds = DeseqDataSet(
            counts         = counts_df,
            metadata       = meta_df[['condition']],
            design_factors = "condition",
            ref_level      = ["condition", "Young"],
            refit_cooks    = True,
            n_cpus         = 8
        )
        dds.deseq2()

        stat = DeseqStats(dds, contrast=["condition", "Aged", "Young"], n_cpus=8)
        stat.summary()
        res = stat.results_df.copy()
        res['peak']      = res.index          # chr1:start-end
        res['cell_type'] = ct_label
        res['tissue']    = tissue

        out_file = os.path.join(out_dir, f"{tissue}_{ct_label}_Aged_vs_Young_DAR.tsv")
        res.to_csv(out_file, sep='\t', index=False)

        n_open  = ((res['padj'] < 0.05) & (res['log2FoldChange'] > 0)).sum()
        n_close = ((res['padj'] < 0.05) & (res['log2FoldChange'] < 0)).sum()
        print(f"    Saved: {out_file}")
        print(f"    Significant (padj<0.05): opening={n_open}  closing={n_close}")
        return res

    except Exception as e:
        print(f"    ERROR in DESeq2: {e}")
        return None

# ════════════════════════════════════════════════════════════
# KIDNEY
# ════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("Loading Kidney h5ad...")
adata_k = sc.read_h5ad(
    f"{SCIENCE_DIR}/kidney_processed/GSM8774007_Kidney_peak_count.h5ad"
)
print(f"Shape: {adata_k.shape}")

# Print all available cell types for reference
print("\nAll Kidney Main_cell_type labels:")
for ct in sorted(adata_k.obs['Main_cell_type'].unique()):
    n = (adata_k.obs['Main_cell_type'] == ct).sum()
    print(f"  {n:>7,}  {ct}")

print("\nKidney Age breakdown:")
print(adata_k.obs['Age'].value_counts().to_string())

for our_ct, science_cts in kidney_map.items():
    print(f"\n--- Kidney: {our_ct} ({science_cts}) ---")
    counts_df, meta_df = pseudobulk(adata_k, our_ct, science_cts)
    if counts_df is not None:
        run_deseq2(counts_df, meta_df, "Kidney", our_ct, OUT_DIR)

del adata_k
print("\nKidney done, memory freed.")

# ════════════════════════════════════════════════════════════
# LUNG
# ════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("Loading Lung h5ad...")
adata_l = sc.read_h5ad(
    f"{SCIENCE_DIR}/lung_processed/GSM8774006_Lung_peak_count.h5ad"
)
print(f"Shape: {adata_l.shape}")

print("\nAll Lung Main_cell_type labels:")
for ct in sorted(adata_l.obs['Main_cell_type'].unique()):
    n = (adata_l.obs['Main_cell_type'] == ct).sum()
    print(f"  {n:>7,}  {ct}")

print("\nLung Age breakdown:")
print(adata_l.obs['Age'].value_counts().to_string())

for our_ct, science_cts in lung_map.items():
    print(f"\n--- Lung: {our_ct} ({science_cts}) ---")
    counts_df, meta_df = pseudobulk(adata_l, our_ct, science_cts)
    if counts_df is not None:
        run_deseq2(counts_df, meta_df, "Lung", our_ct, OUT_DIR)

del adata_l
print("\nLung done.")

print("\n" + "="*60)
print("All pseudobulk + DESeq2 done.")
print(f"Output: {OUT_DIR}")
PYEOF

# ════════════════════════════════════════════════════════════
# Fig 13: Aging DAR counts per cell type (Kidney / Lung)
# ════════════════════════════════════════════════════════════
module load r/4.4.2

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

AGING_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/aging_DARs"
OUT_DIR   <- AGING_DIR

# ── Read all TSV files and count DARs (padj < 0.05) ─────────
files <- list.files(AGING_DIR, pattern = "_Aged_vs_Young_DAR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No aging DAR TSV files found in: ", AGING_DIR)

counts <- lapply(files, function(f) {
  bname  <- basename(f)
  tissue <- sub("_.*", "", bname)
  ct     <- sub(paste0("^", tissue, "_"), "", sub("_Aged_vs_Young_DAR\\.tsv$", "", bname))

  tab <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(NULL)

  sig <- tab[!is.na(padj) & padj < 0.05]
  data.frame(
    tissue    = tissue,
    cell_type = ct,
    n_opening = sum(sig$log2FoldChange > 0, na.rm = TRUE),
    n_closing = sum(sig$log2FoldChange < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

write.csv(counts, file.path(OUT_DIR, "aging_DAR_counts_summary.csv"), row.names = FALSE)

# ── Long format for plotting ─────────────────────────────────
plot_df <- counts |>
  tidyr::pivot_longer(c(n_opening, n_closing),
                      names_to = "direction", values_to = "n_DAR") |>
  mutate(direction = ifelse(direction == "n_opening", "Opening", "Closing"))

dir_colors <- c(Opening = "#B2182B", Closing = "#2166AC")

make_panel <- function(df, tis, panel_label) {
  df <- df |> filter(tissue == tis) |>
    mutate(cell_type = factor(cell_type, levels = sort(unique(cell_type))))
  if (nrow(df) == 0)
    return(ggplot() + labs(title = panel_label) + theme_void())

  ggplot(df, aes(x = cell_type, y = n_DAR, fill = direction)) +
    geom_col(position = position_dodge(0.75), width = 0.65) +
    geom_text(aes(label = ifelse(n_DAR > 0, scales::comma(n_DAR), "")),
              position = position_dodge(0.75),
              vjust = -0.4, size = 3, color = "grey20") +
    scale_fill_manual(values = dir_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = panel_label,
         x = "Cell type", y = "Number of aging DARs (padj < 0.05)",
         fill = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 10),
      strip.text       = element_text(face = "bold"),
      legend.position  = "top"
    )
}

pA <- make_panel(plot_df, "Kidney", "A  Kidney")
pB <- make_panel(plot_df, "Lung",   "B  Lung")

fig13 <- (pA / pB) +
  plot_annotation(
    title    = "Fig. 13  Aging DARs from Science paper (Aged vs Young, padj < 0.05)",
    subtitle = "Pseudo-bulk DESeq2 | No peak pre-filtering",
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  )

ggsave(file.path(OUT_DIR, "Fig13_aging_DAR_counts.pdf"),
       fig13, width = 10, height = 9)
ggsave(file.path(OUT_DIR, "Fig13_aging_DAR_counts.png"),
       fig13, width = 10, height = 9, dpi = 300)

message("Saved: Fig13_aging_DAR_counts (.pdf + .png)")
message("Summary table: aging_DAR_counts_summary.csv")
REOF
