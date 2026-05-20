#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=pseudobulk_Lung_5movs21mo
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# ════════════════════════════════════════════════════════════════════
# Aging DAR re-call: Aged (21mo) vs Adult (5mo) — PURE aging contrast
# (Lung companion to 02a_pseudobulk_Kidney_5movs21mo.sh)
#
# ⚠️ AGE_LABEL ⚠️ — if Lu et al's middle-age label is not 'Adult'
# (could be '5mo', 'Middle', 'Adult_5mo' etc.), change the
# YOUNG_LABEL value at the top of the python block.
# ════════════════════════════════════════════════════════════════════

set -euo pipefail

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings, re
warnings.filterwarnings('ignore')

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"
OUT_DIR     = f"{SCIENCE_DIR}/aging_DARs_5movs21mo"           # ← new dir
os.makedirs(OUT_DIR, exist_ok=True)

# ⚠️ AGE_LABEL ⚠️
YOUNG_LABEL = "Adult"     # was "Young" (1mo) in the original run
AGED_LABEL  = "Aged"

MIN_CELLS_PER_SAMPLE = 10
MIN_SAMPLES          = 4
MIN_PEAK_FRAC        = 0.05
MIN_MEAN_COUNT       = 1

lung_map = {
    "AT2":         ["Type II alveolar epithelial cells"],
    "B":           ["Lymphoid cells_B cells"],
    "Ciliated":    ["Ciliated cells"],
    "EC-vasc":     ["Vascular endothelial cells", "Vascular endothelial cells_Aerocytes"],
    "Eosinophils": ["Myeloid cells_Eosinophils"],
    "Fib":         ["Fibroblasts"],
    "Mac-alv":     ["Myeloid cells_Alveolar macrophages"],
    "Mac-inter":   ["Myeloid cells_Interstitial macrophages"],
    "NK":          ["Lymphoid cells_NK cells"],
    "Pen":         ["Neuroendocrine cells"],
    "SMCs":        ["Smooth muscle cells"],
    "T":           ["Lymphoid cells_T cells"],
}

def pseudobulk(adata, ct_label, science_labels):
    mask = adata.obs['Main_cell_type'].isin(science_labels)
    ct   = adata[mask].copy()
    ct   = ct[ct.obs['Age'].isin([YOUNG_LABEL, AGED_LABEL])].copy()

    print(f"    Cells ({YOUNG_LABEL}+{AGED_LABEL}): {ct.shape[0]}")
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
        # Extract sex from Sample name, e.g. 'Adult_Female_3' -> 'Female'
        sex_match = re.search(r'(Female|Male)', str(samp))
        sex = sex_match.group(1) if sex_match else "Unknown"
        pb_counts[samp] = counts
        pb_meta[samp]   = {'Age': age, 'sex': sex, 'n_cells': n}

    if len(pb_counts) < MIN_SAMPLES:
        print(f"    SKIP: only {len(pb_counts)} pseudo-bulk samples"); return None, None

    counts_df = pd.DataFrame(pb_counts, index=ct.var_names).T.astype(int)
    meta_df   = pd.DataFrame(pb_meta).T

    n_young = (meta_df['Age'] == YOUNG_LABEL).sum()
    n_aged  = (meta_df['Age'] == AGED_LABEL).sum()
    if n_young == 0 or n_aged == 0:
        print("    SKIP: missing one age group"); return None, None

    print(f"    Pseudo-bulk: {YOUNG_LABEL}={n_young}, {AGED_LABEL}={n_aged}, Peaks={counts_df.shape[1]}")

    total_counts = counts_df.sum(axis=1)
    cpm_df       = counts_df.div(total_counts, axis=0) * 1e6
    max_cpm      = cpm_df.max(axis=0)
    cpm_thresh   = max_cpm.quantile(0.75)
    keep_cpm     = max_cpm >= cpm_thresh

    min_s      = max(2, int(MIN_PEAK_FRAC * len(pb_counts)))
    keep_frac  = (counts_df > 0).sum(axis=0) >= min_s

    keep_mean  = counts_df.mean(axis=0) >= MIN_MEAN_COUNT

    keep       = keep_cpm & keep_frac & keep_mean
    counts_df  = counts_df.loc[:, keep]
    print(f"    Peaks after CPM-top25% + frac>={min_s} + mean>={MIN_MEAN_COUNT}: {counts_df.shape[1]}")

    if counts_df.shape[1] < 100:
        print("    SKIP: too few peaks"); return None, None

    return counts_df, meta_df

def run_deseq2(counts_df, meta_df, tissue, ct_label, out_dir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds  import DeseqStats

    meta_df = meta_df.copy()
    meta_df['condition'] = meta_df['Age'].map({YOUNG_LABEL: YOUNG_LABEL, AGED_LABEL: AGED_LABEL})

    # Check if sex is balanced enough to be a covariate (need ≥2 each)
    n_female = (meta_df['sex'] == 'Female').sum()
    n_male   = (meta_df['sex'] == 'Male').sum()
    use_sex  = (n_female >= 2 and n_male >= 2)
    print(f"    Sex: Female={n_female}, Male={n_male}, use_sex_covariate={use_sex}")

    try:
        if use_sex:
            # design = ~ sex + condition  (matches Lu 2026 sex-stratified approach)
            dds = DeseqDataSet(
                counts         = counts_df,
                metadata       = meta_df[['sex', 'condition']],
                design_factors = ['sex', 'condition'],
                ref_level      = [['condition', YOUNG_LABEL], ['sex', 'Female']],
                refit_cooks    = True,
                n_cpus         = 8
            )
        else:
            # Fallback: no sex covariate
            dds = DeseqDataSet(
                counts         = counts_df,
                metadata       = meta_df[['condition']],
                design_factors = "condition",
                ref_level      = ["condition", YOUNG_LABEL],
                refit_cooks    = True,
                n_cpus         = 8
            )
        dds.deseq2()

        stat = DeseqStats(dds, contrast=["condition", AGED_LABEL, YOUNG_LABEL], n_cpus=8)
        stat.summary()
        res = stat.results_df.copy()
        res['peak']      = res.index
        res['cell_type'] = ct_label
        res['tissue']    = tissue

        out_file = os.path.join(out_dir,
            f"{tissue}_{ct_label}_{AGED_LABEL}_vs_{YOUNG_LABEL}_DAR.tsv")
        res.to_csv(out_file, sep='\t', index=False)

        n_open  = ((res['padj'] < 0.05) & (res['log2FoldChange'] > 0)).sum()
        n_close = ((res['padj'] < 0.05) & (res['log2FoldChange'] < 0)).sum()
        print(f"    Saved: {out_file}")
        print(f"    Significant (padj<0.05): opening={n_open}  closing={n_close}")
        return res

    except Exception as e:
        print(f"    ERROR in DESeq2: {e}"); return None

# ── Lung ─────────────────────────────────────────────────────
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

print("\nLung Age breakdown (verify '{YOUNG_LABEL}' and '{AGED_LABEL}' exist):"
      .format(YOUNG_LABEL=YOUNG_LABEL, AGED_LABEL=AGED_LABEL))
print(adata_l.obs['Age'].value_counts().to_string())

if YOUNG_LABEL not in adata_l.obs['Age'].unique():
    raise ValueError(
        f"'{YOUNG_LABEL}' not found in Age column. "
        f"Available: {sorted(adata_l.obs['Age'].unique())}. "
        f"Update YOUNG_LABEL at the top of the script."
    )

for our_ct, science_cts in lung_map.items():
    print(f"\n--- Lung: {our_ct} ---")
    counts_df, meta_df = pseudobulk(adata_l, our_ct, science_cts)
    if counts_df is not None:
        run_deseq2(counts_df, meta_df, "Lung", our_ct, OUT_DIR)

print(f"\nLung done. Outputs: {OUT_DIR}")
PYEOF

# ════════════════════════════════════════════════════════════
# Fig 13 (5mo vs 21mo version): Aging DAR counts per cell type
# ════════════════════════════════════════════════════════════
module load r/4.4.2

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

AGING_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/aging_DARs_5movs21mo"

# Match new filename pattern: e.g. Kidney_PT_Aged_vs_Adult_DAR.tsv
files <- list.files(AGING_DIR, pattern = "_Aged_vs_Adult_DAR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No aging DAR TSV files found in: ", AGING_DIR)

counts <- lapply(files, function(f) {
  bname  <- basename(f)
  tissue <- sub("_.*", "", bname)
  ct     <- sub(paste0("^", tissue, "_"), "", sub("_Aged_vs_Adult_DAR\\.tsv$", "", bname))
  tab    <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  sig <- tab[!is.na(padj) & padj < 0.05]
  data.frame(tissue = tissue, cell_type = ct,
             n_opening = sum(sig$log2FoldChange > 0, na.rm = TRUE),
             n_closing = sum(sig$log2FoldChange < 0, na.rm = TRUE),
             stringsAsFactors = FALSE)
}) |> bind_rows()

write.csv(counts, file.path(AGING_DIR, "aging_DAR_counts_summary_5movs21mo.csv"),
          row.names = FALSE)

plot_df <- counts |>
  tidyr::pivot_longer(c(n_opening, n_closing),
                      names_to = "direction", values_to = "n_DAR") |>
  mutate(direction = ifelse(direction == "n_opening", "Opening", "Closing"))

dir_colors <- c(Opening = "#B2182B", Closing = "#2166AC")

make_panel <- function(df, tis, panel_label) {
  df <- df |> filter(tissue == tis)
  if (nrow(df) == 0) return(ggplot() + labs(title = panel_label) + theme_void())

  totals <- df |>
    group_by(cell_type) |>
    summarise(total = sum(n_DAR), .groups = "drop") |>
    filter(total > 0) |>
    arrange(desc(total))
  if (nrow(totals) == 0) return(ggplot() + labs(title = panel_label) + theme_void())

  df <- df |>
    filter(cell_type %in% totals$cell_type) |>
    mutate(cell_type = factor(cell_type, levels = totals$cell_type))
  if (nrow(df) == 0) return(ggplot() + labs(title = panel_label) + theme_void())
  ggplot(df, aes(x = cell_type, y = n_DAR, fill = direction)) +
    geom_col(position = position_dodge(0.75), width = 0.65) +
    geom_text(aes(label = ifelse(n_DAR > 0, scales::comma(n_DAR), "")),
              position = position_dodge(0.75),
              vjust = -0.4, size = 3, color = "grey20") +
    scale_fill_manual(values = dir_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = panel_label,
         x = "Cell type", y = "Number of aging DARs (padj < 0.05)", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          axis.text.x     = element_text(angle = 40, hjust = 1, size = 10),
          legend.position = "top")
}

fig13 <- (make_panel(plot_df, "Kidney", "A  Kidney") /
           make_panel(plot_df, "Lung",   "B  Lung")) +
  plot_annotation(
    title    = "Fig. 13  Aging DARs (Aged 21mo vs Adult 5mo, padj < 0.05) — PURE aging contrast",
    subtitle = "Pseudo-bulk DESeq2 | non-zero in >=5% of samples",
    theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                  plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
  )

ggsave(file.path(AGING_DIR, "Fig13_aging_DAR_counts_5movs21mo.pdf"),
       fig13, width = 10, height = 9)
ggsave(file.path(AGING_DIR, "Fig13_aging_DAR_counts_5movs21mo.png"),
       fig13, width = 10, height = 9, dpi = 300)
message("Saved: Fig13_aging_DAR_counts_5movs21mo (.pdf + .png)")
REOF
