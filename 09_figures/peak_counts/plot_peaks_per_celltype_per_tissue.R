#!/usr/bin/env Rscript
# ==============================================================================
# Per-tissue publication figures: cell type-specific peaks (one PDF per tissue)
#
# Companion to plot_peaks_per_celltype.R (which produces the combined
# facet_wrap version). This script generates standalone, thesis-figure-ready
# plots for each tissue, with in-plot statistics and a tissue-specific palette.
#
# Outputs (in IN_DIR/per_tissue/):
#   Fig_peaks_Kidney.pdf / .png
#   Fig_peaks_Lung.pdf   / .png
#   Fig_peaks_Aorta.pdf  / .png
#   Fig_peaks_Tcells.pdf / .png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
IN_DIR  <- ifelse(
  dir.exists("/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype"),
  "/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype",
  "."
)
OUT_DIR <- file.path(IN_DIR, "per_tissue")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

raw_csv <- file.path(IN_DIR, "peaks_per_celltype.csv")
sum_csv <- file.path(IN_DIR, "peaks_per_celltype_summary.csv")
stopifnot(file.exists(raw_csv), file.exists(sum_csv))

df_raw <- read.csv(raw_csv, stringsAsFactors = FALSE)
sm_raw <- read.csv(sum_csv, stringsAsFactors = FALSE)

# ---------------------------------------------------------------
# Filter stale / contaminating peak files
# (kept in sync with plot_peaks_per_celltype.R)
# ---------------------------------------------------------------
tcells_drop <- c("Aged", "Juvenile", "Young_acute", "Young_chronic", "Young_control")
lung_drop   <- c("Endothelial_cells")

drop_idx <- (df_raw$tissue == "Tcells" & df_raw$cell_type %in% tcells_drop) |
            (df_raw$tissue == "Lung"   & df_raw$cell_type %in% lung_drop)
if (any(drop_idx)) {
  message(sprintf("Filtered %d stale rows", sum(drop_idx)))
}
df <- df_raw[!drop_idx, , drop = FALSE]

# ---------------------------------------------------------------
# Tissue palette & labels
# ---------------------------------------------------------------
tissue_levels  <- c("Kidney", "Lung", "Aorta", "Tcells")
tissue_labels  <- c(Kidney = "Kidney", Lung = "Lung",
                    Aorta  = "Aorta",  Tcells = "T cells")
tissue_palette <- c(
  Kidney = "#51247A",   # UQ purple
  Lung   = "#0072B2",   # blue
  Aorta  = "#D55E00",   # vermillion
  Tcells = "#009E73"    # bluish green
)
darken <- function(hex, factor = 0.55) {
  v <- col2rgb(hex) / 255 * factor
  rgb(v[1], v[2], v[3])
}

# ---------------------------------------------------------------
# Recompute summary from filtered df
# ---------------------------------------------------------------
sm <- df %>%
  group_by(tissue) %>%
  summarise(
    n_cell_types        = dplyr::n(),
    total_peaks_precons = sum(n_peaks),
    median_peaks        = median(n_peaks),
    min_peaks           = min(n_peaks),
    max_peaks           = max(n_peaks),
    .groups = "drop"
  ) %>%
  left_join(sm_raw %>% select(tissue, n_peaks_consensus), by = "tissue")

# ---------------------------------------------------------------
# Plot function
# ---------------------------------------------------------------
make_tissue_plot <- function(tis) {
  d <- df %>%
    filter(tissue == tis) %>%
    arrange(desc(n_peaks)) %>%
    mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

  s <- sm %>% filter(tissue == tis)
  fill_col   <- unname(tissue_palette[tis])
  stroke_col <- darken(fill_col, 0.55)

  n_cons_num <- suppressWarnings(as.numeric(s$n_peaks_consensus))
  reduction  <- ifelse(is.na(n_cons_num) | n_cons_num == 0,
                       NA_real_, s$total_peaks_precons / n_cons_num)

  stats_lines <- c(
    sprintf("Cell types:  %d", s$n_cell_types),
    sprintf("Total peaks:  %s",
            formatC(s$total_peaks_precons, format = "d", big.mark = ",")),
    sprintf("Median:  %s",
            formatC(s$median_peaks, format = "d", big.mark = ",")),
    sprintf("Range:  %s — %s",
            formatC(s$min_peaks, format = "d", big.mark = ","),
            formatC(s$max_peaks, format = "d", big.mark = ",")),
    if (!is.na(reduction))
      sprintf("Consensus:  %s  (%.1f× reduction)",
              formatC(n_cons_num, format = "d", big.mark = ","), reduction)
    else
      sprintf("Consensus:  —")
  )
  stats_text <- paste(stats_lines, collapse = "\n")

  y_max <- max(d$n_peaks) / 1000
  y_breaks <- pretty(c(0, y_max), n = 6)

  p <- ggplot(d, aes(x = cell_type, y = n_peaks / 1000)) +
    geom_col(fill = fill_col, colour = stroke_col,
             width = 0.7, linewidth = 0.35) +
    geom_hline(yintercept = s$median_peaks / 1000,
               linetype = "dashed", colour = "grey55", linewidth = 0.45) +
    annotate("label",
             x = Inf, y = s$median_peaks / 1000,
             label = sprintf("median = %s",
                             formatC(s$median_peaks, format = "d", big.mark = ",")),
             hjust = 1.05, vjust = -0.25,
             label.size = 0,
             label.padding = unit(0.18, "lines"),
             fill = alpha("white", 0.88),
             size = 2.95,
             colour = "grey35", fontface = "italic") +
    geom_text(aes(label = formatC(n_peaks, format = "d", big.mark = ",")),
              vjust = -0.55, size = 2.95, colour = "grey15") +
    annotate("label",
             x = Inf, y = Inf,
             label = stats_text,
             hjust = 1.04, vjust = 1.05,
             label.size = 0.3,
             label.padding = unit(0.55, "lines"),
             label.r = unit(0.12, "lines"),
             fill = alpha("white", 0.92),
             colour = "grey25",
             size = 3.15, lineheight = 1.3, family = "sans") +
    scale_y_continuous(
      breaks = y_breaks,
      limits = c(0, max(y_breaks) * 1.22),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = expression("Number of peaks (" * 10^3 * ")"),
      title    = paste0(tissue_labels[tis], " — cell type-specific peaks"),
      subtitle = "MACS2 narrowPeak (q < 0.05), prior to per-tissue consensus merging"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", size = 14, hjust = 0,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(size = 10, colour = "grey35", hjust = 0,
                                       margin = margin(b = 10)),
      axis.text.x       = element_text(angle = 35, hjust = 1, size = 10,
                                       colour = "grey15"),
      axis.text.y       = element_text(size = 10, colour = "grey15"),
      axis.title.y      = element_text(size = 11.5, margin = margin(r = 7)),
      axis.line         = element_line(colour = "grey20", linewidth = 0.45),
      axis.ticks        = element_line(colour = "grey20", linewidth = 0.45),
      axis.ticks.length = unit(3, "pt"),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.margin       = margin(14, 16, 10, 14)
    )

  p
}

# ---------------------------------------------------------------
# Generate & save
# ---------------------------------------------------------------
for (tis in tissue_levels) {
  if (!(tis %in% df$tissue)) next

  p <- make_tissue_plot(tis)
  n_ct <- sum(df$tissue == tis)
  # Width scales with # cell types: ~0.45 in per bar, min 5.5"
  w <- max(5.5, 3.2 + 0.45 * n_ct)
  h <- 5.4

  pdf_path <- file.path(OUT_DIR, sprintf("Fig_peaks_%s.pdf", tis))
  png_path <- file.path(OUT_DIR, sprintf("Fig_peaks_%s.png", tis))

  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf)
  ggsave(png_path, p, width = w, height = h, dpi = 600)
  message(sprintf("Saved: %s  (%.1f x %.1f in., %d cell types)",
                  basename(pdf_path), w, h, n_ct))
}

cat("\nOutputs ->", OUT_DIR, "\n")
