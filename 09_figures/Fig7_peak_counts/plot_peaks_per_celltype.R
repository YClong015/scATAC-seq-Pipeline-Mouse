#!/usr/bin/env Rscript
# Fig 7 - cell-type peak counts per tissue (pre-consensus), faceted bar chart.
# Reads the CSVs written by count_peaks_per_celltype.slurm.

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)

## Paths (edit IN_DIR if running locally with downloaded CSVs)
IN_DIR  <- ifelse(
  dir.exists("/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype_pooled"),
  "/QRISdata/Q8448/Mouse_disease_data/QC_figures/peaks_per_celltype_pooled",
  "."
)
OUT_DIR <- IN_DIR

raw.csv <- file.path(IN_DIR, "peaks_per_celltype.csv")
sum.csv <- file.path(IN_DIR, "peaks_per_celltype_summary.csv")

stopifnot(file.exists(raw.csv), file.exists(sum.csv))

df.raw  <- read.csv(raw.csv,  stringsAsFactors = FALSE)
sm.raw  <- read.csv(sum.csv,  stringsAsFactors = FALSE)

## drop stale peak files (Tcells condition-pseudobulks; Lung pre-merge EC label)
tcells.drop <- c("Aged", "Juvenile", "Young_acute", "Young_chronic", "Young_control")
lung.drop   <- c("Endothelial_cells")

drop.idx <- (df.raw$tissue == "Tcells" & df.raw$cell_type %in% tcells.drop) |
            (df.raw$tissue == "Lung"   & df.raw$cell_type %in% lung.drop)
if (any(drop.idx)) {
  dropped <- paste(df.raw$tissue[drop.idx], df.raw$cell_type[drop.idx],
                   sep = ":", collapse = ", ")
  message(sprintf("Filtered %d stale rows: %s", sum(drop.idx), dropped))
}
df <- df.raw[!drop.idx, , drop = FALSE]

## Tidy
tissue.levels <- c("Kidney", "Lung", "Aorta", "Tcells")
tissue.labels <- c(Kidney = "Kidney", Lung = "Lung",
                   Aorta  = "Aorta",  Tcells = "T cells")

df <- df %>%
  filter(tissue %in% tissue.levels) %>%
  mutate(tissue = factor(tissue, levels = tissue.levels)) %>%
  arrange(tissue, desc(n_peaks)) %>%
  mutate(
    ct_key   = paste(tissue, cell_type, sep = "__"),
    ct_key   = factor(ct_key, levels = unique(ct_key))
  )

## Recompute summary from filtered df (do not trust sm.raw for Tcells)
sm <- df %>%
  group_by(tissue) %>%
  summarise(
    n_cell_types         = n(),
    total_peaks_precons  = sum(n_peaks),
    median_peaks         = median(n_peaks),
    min_peaks            = min(n_peaks),
    max_peaks            = max(n_peaks),
    .groups = "drop"
  ) %>%
  left_join(
    sm.raw %>% select(tissue, n_peaks_consensus),
    by = "tissue"
  )

ct.label.map <- setNames(df$cell_type, df$ct_key)

# Tissue palette (UQ purple + complementary, colour-blind friendly)
tissue.palette <- c(
  Kidney = "#51247A",   # UQ purple
  Lung   = "#0072B2",   # blue
  Aorta  = "#D55E00",   # vermillion
  Tcells = "#009E73"    # bluish green
)

## Main figure: bar chart, faceted by tissue
y.max <- max(df$n_peaks) / 1000
y.breaks <- pretty(c(0, y.max), n = 6)

p.main <- ggplot(df, aes(x = ct_key, y = n_peaks / 1000, fill = tissue)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.25) +
  geom_text(
    aes(label = formatC(n_peaks, format = "d", big.mark = ",")),
    vjust = -0.45, size = 2.7, colour = "grey15"
  ) +
  facet_wrap(~ tissue, scales = "free_x", nrow = 1,
             labeller = as_labeller(tissue.labels)) +
  scale_x_discrete(labels = ct.label.map) +
  scale_y_continuous(
    breaks = y.breaks,
    limits = c(0, max(y.breaks) * 1.12),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = tissue.palette, guide = "none") +
  labs(
    x = NULL,
    y = expression("Number of peaks (" * 10^3 * ")"),
    title    = "Cell type-specific peaks per tissue",
    subtitle = "MACS2 narrowPeak counts per cell type, prior to per-tissue consensus merging"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle   = element_text(size = 10, colour = "grey35", hjust = 0,
                                   margin = margin(b = 8)),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 11, margin = margin(4, 0, 4, 0)),
    panel.spacing.x  = unit(0.9, "lines"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9, colour = "grey15"),
    axis.text.y      = element_text(size = 9, colour = "grey15"),
    axis.title.y     = element_text(size = 11, margin = margin(r = 6)),
    axis.line        = element_line(colour = "grey25", linewidth = 0.4),
    axis.ticks       = element_line(colour = "grey25", linewidth = 0.4),
    plot.margin      = margin(10, 12, 8, 10)
  )

ggsave(file.path(OUT_DIR, "Fig7_peaks_per_celltype.pdf"),
       p.main, width = 11, height = 5.2, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7_peaks_per_celltype.png"),
       p.main, width = 11, height = 5.2, dpi = 600)
message("Saved: Fig7_peaks_per_celltype")

## Summary table (human-readable + as figure panel)
sm <- sm %>%
  mutate(tissue = factor(tissue, levels = tissue.levels)) %>%
  arrange(tissue)

fmt.int <- function(x) ifelse(is.na(x) | x == "NA",
                              "-",
                              formatC(as.numeric(x), format = "d", big.mark = ","))

sm.disp <- data.frame(
  Tissue              = tissue.labels[as.character(sm$tissue)],
  `Cell types`        = sm$n_cell_types,
  `Total peaks (pre-consensus)` = vapply(sm$total_peaks_precons, fmt.int, character(1)),
  `Median per cell type`        = vapply(sm$median_peaks,        fmt.int, character(1)),
  `Min`                         = vapply(sm$min_peaks,           fmt.int, character(1)),
  `Max`                         = vapply(sm$max_peaks,           fmt.int, character(1)),
  `Consensus peaks`             = vapply(sm$n_peaks_consensus,   fmt.int, character(1)),
  check.names = FALSE
)

write.csv(sm.disp,
          file.path(OUT_DIR, "peaks_per_celltype_summary_formatted.csv"),
          row.names = FALSE)
message("Saved: peaks_per_celltype_summary_formatted.csv")

# tableGrob for embedding in the figure
tt <- ttheme_minimal(
  core    = list(fg_params = list(cex = 0.78, hjust = 0.5, x = 0.5),
                 bg_params = list(fill = c("white", "grey97"))),
  colhead = list(fg_params = list(cex = 0.82, fontface = "bold", col = "white"),
                 bg_params = list(fill = "#51247A"))
)

tbl.grob <- tableGrob(sm.disp, rows = NULL, theme = tt)

# Combined figure: bar chart on top, summary table below
p.combined <- p.main /
  wrap_elements(full = tbl.grob) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, "Fig7_peaks_per_celltype_with_summary.pdf"),
       p.combined, width = 11, height = 7.6, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7_peaks_per_celltype_with_summary.png"),
       p.combined, width = 11, height = 7.6, dpi = 600)
message("Saved: Fig7_peaks_per_celltype_with_summary")

## Console summary
cat("\n=== Summary ===\n")
print(sm.disp, row.names = FALSE)
cat("\nOutputs ->", OUT_DIR, "\n")
