#!/usr/bin/env Rscript
# Figure 14 - Enrichment of disease DARs within ageing DARs (Kidney + Lung)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Paths: two per-tissue Fisher overlap tables (from 03a/03b), merged in-script.
IN_KIDNEY <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Kidney_5movs21mo/overlap_stats.csv"
IN_LUNG   <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo/overlap_stats.csv"
OUT_DIR   <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/Fig14_output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(IN_KIDNEY), file.exists(IN_LUNG))


KEEP <- list(
  Kidney = c("DCT", "PC", "PT", "TAL"),
  Lung   = c("AT2", "B", "EC-vasc", "Mac-alv", "T")
)

# -- Load both CSVs and merge ------------------------------
read_overlap <- function(path, tissue_label) {
  tab <- read.csv(path, stringsAsFactors = FALSE)
  # Single-tissue CSVs may omit the 'tissue' column; backfill it.
  if (!"tissue" %in% names(tab)) tab$tissue <- tissue_label
  tab
}

fisher <- dplyr::bind_rows(
  read_overlap(IN_KIDNEY, "Kidney"),
  read_overlap(IN_LUNG,   "Lung")
) %>%
  filter(mapply(function(tis, ct) ct %in% KEEP[[tis]], tissue, cell_type))

message("Loaded ", nrow(fisher), " (tissue x cell-type) rows after whitelist.")

# p-value -> labeled string (specific numeric value above each bar)
fmt_pval <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ "",
    p == 0      ~ "P<1e-300",
    p < 0.001   ~ sprintf("P=%.1e", p),
    TRUE        ~ sprintf("P=%.3f", p)
  )
}

# Long format, capped OR at 12 (preserves visual hierarchy when one
# cell type has an extreme outlier OR; documented in y-axis label)
OR_CAP <- 12
fisher_long <- fisher %>%
  select(tissue, cell_type, fold_open, fold_close, pval_open, pval_close) %>%
  pivot_longer(cols = c(fold_open, fold_close),
               names_to = "direction", values_to = "OR") %>%
  mutate(
    direction  = ifelse(direction == "fold_open", "Opening", "Closing"),
    direction  = factor(direction, levels = c("Opening", "Closing")),
    pval       = ifelse(direction == "Opening", pval_open, pval_close),
    pval_label = fmt_pval(pval),
    OR_plot    = pmin(OR, OR_CAP, na.rm = TRUE)
  )

# -- Palette (Nature) ---------------------------------------
dir_colors <- c("Opening" = "#D62728", "Closing" = "#1F77B4")

# -- Per-panel plot helper ----------------------------------
make_panel <- function(tis, panel_letter = NULL) {
  d <- fisher_long %>%
    filter(tissue == tis) %>%
    mutate(cell_type = factor(cell_type, levels = KEEP[[tis]]))

  # Standalone: no panel letter. Combined layout: letter inlined with the title.
  title_str <- if (is.null(panel_letter)) tis
               else paste0(panel_letter, "   ", tis)

  if (nrow(d) == 0) {
    return(ggplot() + labs(title = title_str) + theme_void())
  }

  ggplot(d, aes(x = cell_type, y = OR_plot, fill = direction)) +
    geom_col(position = position_dodge(0.78), width = 0.7,
             colour = "white", linewidth = 0.35) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey25", linewidth = 1.0) +
    annotate("text", x = Inf, y = 1,
             label = "OR = 1",
             hjust = 1.05, vjust = -0.5,
             size = 5.5, fontface = "italic",
             colour = "grey25", family = "sans") +
    # Exact P-values above each bar (no significance stars / "ns")
    geom_text(aes(label = pval_label, group = direction),
              position = position_dodge(0.78),
              vjust = -0.4, size = 5.0,
              family = "sans", colour = "black") +
    scale_fill_manual(values = dir_colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    # Prevent the panel border from clipping the OR = 1 dashed line
    # and the "OR = 1" annotation that sits flush to the right edge.
    coord_cartesian(clip = "off") +
    labs(
      title = title_str,
      x = NULL,
      y = "Fisher exact OR (capped at 12)"
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
      # Full rectangular panel border; axis.line removed to avoid the
      # doubled left/bottom edge that would otherwise overlap it.
      panel.border      = element_rect(colour = "black", fill = NA,
                                       linewidth = 0.8),
      axis.line         = element_blank(),
      axis.ticks        = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = unit(5, "pt"),     # outward (Nature style)
      # Horizontal-only major grid (x-axis is categorical cell types).
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position   = "top",
      legend.justification = "center",
      legend.text       = element_text(size = 18, family = "sans",
                                       colour = "black"),
      legend.key.size   = unit(0.7, "cm"),
      legend.spacing.x  = unit(0.5, "cm"),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.margin       = margin(16, 20, 14, 18)
    )
}

# -- Build standalone single-tissue plots (no A/B prefix) --
pK_solo <- make_panel("Kidney")
pL_solo <- make_panel("Lung")

# -- Build tagged versions for the 2-panel combined figure --
pA <- make_panel("Kidney", "A")
pB <- make_panel("Lung",   "B")

combined <- (pA / pB) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.justification = "center")

# -- Save individual single-tissue PDFs --------------------
W_SOLO <- 12
H_SOLO <- 8

ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_Kidney.pdf"),
       pK_solo, width = W_SOLO, height = H_SOLO, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_Kidney.png"),
       pK_solo, width = W_SOLO, height = H_SOLO, dpi = 400)

ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_Lung.pdf"),
       pL_solo, width = W_SOLO, height = H_SOLO, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_Lung.png"),
       pL_solo, width = W_SOLO, height = H_SOLO, dpi = 400)

# -- Save combined 2-panel figure --------------------------

W_IN <- 13
H_IN <- 16

ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_KidneyLung.pdf"),
       combined, width = W_IN, height = H_IN, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig14_Fisher_OR_KidneyLung.png"),
       combined, width = W_IN, height = H_IN, dpi = 400)

message("Saved 3 figures (Kidney solo, Lung solo, Combined)  ->  ", OUT_DIR)

# Display in RStudio Plots panel when sourced interactively;
# silently skipped under Rscript (no display device).
if (interactive()) print(combined)
