#!/usr/bin/env Rscript
# Fig 7 - per-tissue cell-type peak counts, one panel per tissue + a 2x2 combined.
# Companion to plot_peaks_per_celltype.R (the faceted version).

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)   # for the 2x2 combined panel

## Paths
# Pooled-T-cell version: read + write under peaks_per_celltype_pooled/,
# leaving the original peaks_per_celltype/ (per-subtype version) untouched.
IN_DIR  <- ifelse(
  dir.exists("/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype_pooled"),
  "/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype_pooled",
  "."
)
OUT_DIR <- file.path(IN_DIR, "per_tissue")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

raw.csv <- file.path(IN_DIR, "peaks_per_celltype.csv")
sum.csv <- file.path(IN_DIR, "peaks_per_celltype_summary.csv")
stopifnot(file.exists(raw.csv), file.exists(sum.csv))

df.raw <- read.csv(raw.csv, stringsAsFactors = FALSE)
sm.raw <- read.csv(sum.csv, stringsAsFactors = FALSE)

## Filter stale / contaminating peak files
# (kept in sync with plot_peaks_per_celltype.R)
tcells.drop <- c("Aged", "Juvenile", "Young_acute", "Young_chronic", "Young_control")
lung.drop   <- c("Endothelial_cells")

drop.idx <- (df.raw$tissue == "Tcells" & df.raw$cell_type %in% tcells.drop) |
            (df.raw$tissue == "Lung"   & df.raw$cell_type %in% lung.drop)
if (any(drop.idx)) {
  message(sprintf("Filtered %d stale rows", sum(drop.idx)))
}
df <- df.raw[!drop.idx, , drop = FALSE]

## tissue palette & labels
tissue.levels  <- c("Kidney", "Aorta", "Lung", "Tcells")
tissue.labels  <- c(Kidney = "Kidney", Aorta = "Aorta",
                    Lung   = "Lung",   Tcells = "T cells")
tissue.palette <- c(
  Kidney = "#51247A",   # UQ purple
  Aorta  = "#D55E00",   # vermillion
  Lung   = "#0072B2",   # blue
  Tcells = "#009E73"    # bluish green
)
darken <- function(hex, factor = 0.55) {
  v <- col2rgb(hex) / 255 * factor
  rgb(v[1], v[2], v[3])
}

## Recompute summary from filtered df
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
  left_join(sm.raw %>% select(tissue, n_peaks_consensus), by = "tissue")

## Plot function
make.tissue.plot <- function(tis, panel_letter = NULL) {
  d <- df %>%
    filter(tissue == tis) %>%
    arrange(desc(n_peaks)) %>%
    # Tidy the pooled T-cell bar label
    mutate(cell_type = ifelse(cell_type == "Tcell_pooled",
                              "Pooled T cells", cell_type)) %>%
    mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

  s <- sm %>% filter(tissue == tis)
  fill.col   <- unname(tissue.palette[tis])
  stroke.col <- darken(fill.col, 0.55)

  n.cons.num <- suppressWarnings(as.numeric(s$n_peaks_consensus))
  reduction  <- ifelse(is.na(n.cons.num) | n.cons.num == 0,
                       NA_real_, s$total_peaks_precons / n.cons.num)

  # Single-group tissues (the pooled T-cell set) have no per-cell-type
  # spread, so median / range are not applicable.
  is.pooled <- s$n_cell_types == 1

  stats.lines <- c(
    if (is.pooled) "Pooled peak set (1 group)"
    else sprintf("Cell types:  %d", s$n_cell_types),
    sprintf("Total peaks:  %s",
            formatC(s$total_peaks_precons, format = "d", big.mark = ",")),
    if (is.pooled) "Median:  n/a"
    else sprintf("Median:  %s",
                 formatC(s$median_peaks, format = "d", big.mark = ",")),
    if (is.pooled) "Range:  n/a"
    else sprintf("Range:  %s - %s",
                 formatC(s$min_peaks, format = "d", big.mark = ","),
                 formatC(s$max_peaks, format = "d", big.mark = ",")),
    if (is.pooled)
      "Consensus:  n/a (single pooled set)"
    else if (!is.na(reduction))
      sprintf("Consensus:  %s  (%.1fx reduction)",
              formatC(n.cons.num, format = "d", big.mark = ","), reduction)
    else
      sprintf("Consensus:  -")
  )
  stats.text <- paste(stats.lines, collapse = "\n")

  y.max <- max(d$n_peaks) / 1000
  y.breaks <- pretty(c(0, y.max), n = 6)

  p <- ggplot(d, aes(x = cell_type, y = n_peaks / 1000)) +
    geom_col(fill = fill.col, colour = stroke.col,
             width = 0.7, linewidth = 0.4) +
    geom_hline(yintercept = s$median_peaks / 1000,
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    annotate("label",
             x = Inf, y = s$median_peaks / 1000,
             label = sprintf("median = %s",
                             formatC(s$median_peaks, format = "d", big.mark = ",")),
             hjust = 1.05, vjust = -0.3,
             label.size = 0,
             label.padding = unit(0.25, "lines"),
             fill = alpha("white", 0.92),
             size = 6.0,
             colour = "black", fontface = "italic",
             family = "sans") +
    annotate("label",
             x = Inf, y = Inf,
             label = stats.text,
             hjust = 1.02, vjust = 1.05,
             label.size = 0.5,
             label.padding = unit(0.7, "lines"),
             label.r = unit(0.18, "lines"),
             fill = alpha("white", 0.94),
             colour = "black",
             size = 6.8, lineheight = 1.3, family = "sans") +
    scale_y_continuous(
      breaks = y.breaks,
      limits = c(0, max(y.breaks) * 1.60),   # extra headroom so the inset clears the bars
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = expression("Number of peaks (" * 10^3 * ")"),
      # panel_letter (2x2 mode) prepended inline with the title
      title = if (is.null(panel_letter))
                paste0(tissue.labels[tis], ", cell-type-specific peaks")
              else
                paste0(panel_letter, "   ",
                       tissue.labels[tis], ", cell-type-specific peaks")
    ) +
    theme_classic(base_size = 20, base_family = "sans") +
    theme(
      plot.title        = element_text(face = "bold", size = 26, hjust = 0,
                                       margin = margin(b = 12),
                                       family = "sans", colour = "black"),
      axis.text.x       = element_text(angle = 35, hjust = 1, size = 20,
                                       colour = "black", family = "sans"),
      axis.text.y       = element_text(size = 18, colour = "black",
                                       family = "sans"),
      axis.title.y      = element_text(size = 22, margin = margin(r = 12),
                                       family = "sans", colour = "black"),
      axis.line         = element_line(colour = "black", linewidth = 0.8),
      axis.ticks        = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = unit(5, "pt"),   # outward-pointing ticks (Nature style)
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.margin       = margin(16, 20, 14, 18)
    )

  p
}

## Generate & save individual plots; collect for the 2x2 panel
plot.list <- list()

for (tis in tissue.levels) {
  if (!(tis %in% df$tissue)) next

  p <- make.tissue.plot(tis)
  plot.list[[tis]] <- p

  n.ct <- sum(df$tissue == tis)
  # Canvas scales with # cell types; floor + width-per-CT bumped to match
  # the next round of +4pt font growth without re-introducing collisions.
  w <- max(11, 5 + 0.8 * n.ct)
  h <- 8.5

  pdf.path <- file.path(OUT_DIR, sprintf("Fig7_peaks_%s.pdf", tis))
  png.path <- file.path(OUT_DIR, sprintf("Fig7_peaks_%s.png", tis))

  ggsave(pdf.path, p, width = w, height = h, device = cairo_pdf)
  ggsave(png.path, p, width = w, height = h, dpi = 600)
  message(sprintf("Saved: %s  (%.1f x %.1f in., %d cell types)",
                  basename(pdf.path), w, h, n.ct))
}

## combined 2x2 panel (patchwork), A/B/C/D tags inline with titles
if (length(plot.list) >= 1) {
  tagged.list <- list()
  ordered.names <- tissue.levels[tissue.levels %in% names(plot.list)]
  for (i in seq_along(ordered.names)) {
    tagged.list[[ordered.names[i]]] <-
      make.tissue.plot(ordered.names[i], panel_letter = LETTERS[i])
  }

  combined <- wrap_plots(tagged.list, ncol = 2)

  # Larger canvas: 26 x 18 in. (~66 x 46 cm); landscape-A3-friendly,
  # designed so each of the 4 panels retains thesis-readable fonts
  # after the latest font bump (axis 18-20pt, title 26pt).
  cw <- 26
  ch <- 18

  combo.pdf <- file.path(OUT_DIR, "Fig7_peaks_ALL_2x2.pdf")
  combo.png <- file.path(OUT_DIR, "Fig7_peaks_ALL_2x2.png")
  ggsave(combo.pdf, combined, width = cw, height = ch, device = cairo_pdf)
  ggsave(combo.png, combined, width = cw, height = ch, dpi = 400)
  message(sprintf("Saved 2x2 panel: %s  (%.0f x %.0f in., %d tissues)",
                  basename(combo.pdf), cw, ch, length(ordered)))

  # Display in RStudio Plots panel when sourced interactively.
  # No effect when run via Rscript (the device is null then).
  if (interactive()) print(combined)
}

cat("\nOutputs ->", OUT_DIR, "\n")
