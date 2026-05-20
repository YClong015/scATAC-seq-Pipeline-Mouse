#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --job-name=QC_plots
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
module load r/4.4.2

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

OUT_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
META_DIR <- file.path(OUT_DIR, "metadata_cache")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load each tissue object and extract metadata ──────────────
paths <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
)

meta_list <- list()

for (tissue in names(paths)) {
  message("Loading: ", tissue)
  obj <- readRDS(paths[[tissue]])

  # Detect the correct assay name
  assay_name <- if ("peaks_universal" %in% names(obj@assays)) {
    "peaks_universal"
  } else if ("ATAC" %in% names(obj@assays)) {
    "ATAC"           # T cells use ATAC assay
  } else if ("peaks" %in% names(obj@assays)) {
    "peaks"
  } else {
    names(obj@assays)[1]
  }
  DefaultAssay(obj) <- assay_name

  md <- obj@meta.data
  md$Tissue <- tissue

  # Standardise column names across tissues
  # nCount: Kidney/Lung/Aorta → peak_region_fragments or nCount_peaks
  #         Tcells             → nCount_ATAC
  if ("peak_region_fragments" %in% colnames(md)) {
    md$nCount <- md$peak_region_fragments
  } else if ("nCount_ATAC" %in% colnames(md)) {
    md$nCount <- md$nCount_ATAC
  } else {
    cnt_col <- grep("^nCount", colnames(md), value = TRUE)[1]
    if (!is.na(cnt_col)) md$nCount <- md[[cnt_col]]
  }
  # TSS enrichment
  if (!"TSS.enrichment" %in% colnames(md)) {
    tss_col <- grep("TSS", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(tss_col)) md$TSS.enrichment <- md[[tss_col]]
  }
  # Nucleosome signal
  if (!"nucleosome_signal" %in% colnames(md)) {
    ns_col <- grep("nucleosome", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(ns_col)) md$nucleosome_signal <- md[[ns_col]]
  }
  # pct_reads_in_peaks (FRiP)
  if (!"pct_reads_in_peaks" %in% colnames(md)) {
    frip_col <- grep("pct_reads|FRiP|frip", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(frip_col)) md$pct_reads_in_peaks <- md[[frip_col]]
  }

  message(sprintf("  %s: %d cells | columns: %s",
                  tissue, nrow(md),
                  paste(intersect(c("nCount","TSS.enrichment",
                                    "nucleosome_signal","pct_reads_in_peaks"),
                                  colnames(md)), collapse=", ")))

  # Save metadata cache so QC_plots_debug.R can skip RDS loading next time
  write.csv(md, file.path(META_DIR, paste0(tissue, "_metadata.csv")))
  message("  Cache saved: ", tissue, "_metadata.csv")

  meta_list[[tissue]] <- md
  rm(obj); gc()
}

# ── Combine and factor tissue order ──────────────────────────
all_meta <- bind_rows(meta_list)
all_meta$Tissue <- factor(all_meta$Tissue,
                          levels = c("Kidney", "Lung", "Aorta", "Tcells"))

tissue_colors <- c(
  Kidney = "#2166AC",
  Lung   = "#B2182B",
  Aorta  = "#4DAC26",
  Tcells = "#D6604D"
)

# ── Helper: violin + boxplot ──────────────────────────────────
violin_box <- function(data, yvar, ylabel, log_scale = FALSE,
                       ylim = NULL, title = NULL) {
  p <- ggplot(data, aes(x = Tissue, y = .data[[yvar]], fill = Tissue)) +
    geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white",
                 outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = tissue_colors) +
    labs(x = NULL, y = ylabel, title = title) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 10),
      axis.text.y      = element_text(size = 9),
      axis.title.y     = element_text(size = 10),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    )
  if (log_scale)  p <- p + scale_y_log10()
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ── Panel A: Number of cells per sample ──────────────────────
# Count cells per tissue (bar chart)
cell_counts <- all_meta %>%
  group_by(Tissue) %>%
  summarise(n_cells = n(), .groups = "drop")

pA <- ggplot(cell_counts, aes(x = Tissue, y = n_cells, fill = Tissue)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n_cells)),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = tissue_colors) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of cells",
       title = "A  Number of high-quality cells") +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0, size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ── Panel B: Fragments per cell (nCount) ─────────────────────
pB <- NULL
if ("nCount" %in% colnames(all_meta)) {
  pB <- violin_box(
    all_meta %>% filter(!is.na(nCount) & nCount > 0),
    yvar      = "nCount",
    ylabel    = "Fragments in peaks (log10)",
    log_scale = TRUE,
    title     = "B  Fragments per cell"
  )
}

# ── Panel C: TSS enrichment ───────────────────────────────────
pC <- NULL
if ("TSS.enrichment" %in% colnames(all_meta)) {
  pC <- violin_box(
    all_meta %>% filter(!is.na(TSS.enrichment) & TSS.enrichment > 0),
    yvar   = "TSS.enrichment",
    ylabel = "TSS enrichment score",
    title  = "C  TSS enrichment score"
  )
}

# ── Panel D: FRiP (% reads in peaks) ─────────────────────────
pD <- NULL
if ("pct_reads_in_peaks" %in% colnames(all_meta)) {
  pD <- violin_box(
    all_meta %>% filter(!is.na(pct_reads_in_peaks) & pct_reads_in_peaks > 0),
    yvar   = "pct_reads_in_peaks",
    ylabel = "% reads in peaks (FRiP)",
    title  = "D  Fraction of reads in peaks"
  )
}

# ── Panel A standalone ─────────────────────────────────────────
ggsave(file.path(OUT_DIR, "QC_cell_counts.pdf"),
       pA, width = 4, height = 5)
ggsave(file.path(OUT_DIR, "QC_cell_counts.png"),
       pA, width = 4, height = 5, dpi = 300)
message("Saved: QC_cell_counts (.pdf + .png)")

# ── Panels B + C combined ──────────────────────────────────────
violin_panels <- Filter(Negate(is.null), list(pB, pC))
n_violin <- length(violin_panels)

if (n_violin > 0) {
  combined_violin <- wrap_plots(violin_panels, nrow = 1) +
    plot_annotation(
      title    = "scATAC-seq quality control metrics across tissues",
      subtitle = "Violin plots show per-cell distributions; box shows median and IQR",
      theme    = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
      )
    )
  ggsave(file.path(OUT_DIR, "QC_metrics.pdf"),
         combined_violin, width = n_violin * 4, height = 5)
  ggsave(file.path(OUT_DIR, "QC_metrics.png"),
         combined_violin, width = n_violin * 4, height = 5, dpi = 300)
  message("Saved: QC_metrics (.pdf + .png)")
}

# ── Also print summary table ──────────────────────────────────
summary_tbl <- all_meta %>%
  group_by(Tissue) %>%
  summarise(
    n_cells          = n(),
    median_nCount    = if ("nCount" %in% names(.)) median(nCount, na.rm=TRUE) else NA,
    median_TSS       = if ("TSS.enrichment" %in% names(.)) median(TSS.enrichment, na.rm=TRUE) else NA,
    median_nuc_sig   = if ("nucleosome_signal" %in% names(.)) median(nucleosome_signal, na.rm=TRUE) else NA,
    .groups = "drop"
  )
print(summary_tbl)
write.csv(summary_tbl, file.path(OUT_DIR, "QC_summary_table.csv"), row.names = FALSE)
message("Saved: QC_summary_table.csv")
REOF

echo "Done."
