#!/usr/bin/env Rscript
# Figure 15 - Disease vs ageing log2FC concordance (Kidney tubular + Lung B)

library(ggplot2)
library(dplyr)
library(patchwork)

# -- Paths --------------------------------------------------
KIDNEY_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Kidney_5movs21mo/paired_lfc"
LUNG_DIR   <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo/paired_lfc"
OUT_DIR    <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/Fig15_output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(KIDNEY_DIR), dir.exists(LUNG_DIR))

# -- 4-panel definition (tissue, cell type, panel letter, display label) --
PANELS <- list(
  list(tissue = "Kidney", ct = "TAL", letter = "A", label = "Kidney_TAL"),
  list(tissue = "Kidney", ct = "PT",  letter = "B", label = "Kidney_PT"),
  list(tissue = "Kidney", ct = "DCT", letter = "C", label = "Kidney_DCT"),
  list(tissue = "Lung",   ct = "B",   letter = "D", label = "Lung_B")
)

# -- Colour & alpha palette --------------------------------
sig_levels <- c("Not sig in aging",
                "Aging trending (padj<0.1)",
                "Aging sig (padj<0.05)")
sig_colors <- c("Not sig in aging"          = "grey75",
                "Aging trending (padj<0.1)" = "#F4A582",
                "Aging sig (padj<0.05)"     = "#B2182B")
sig_alpha  <- c("Not sig in aging"          = 0.30,
                "Aging trending (padj<0.1)" = 0.75,
                "Aging sig (padj<0.05)"     = 0.95)

# -- Read + classify one paired BED -------------------------
read_paired <- function(f) {
  dat <- tryCatch(
    read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
               col.names = c("d_chr","d_start","d_end","d_lfc","d_padj",
                             "a_chr","a_start","a_end","a_lfc","a_padj")),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) < 10) return(NULL)

  dat %>%
    mutate(across(c(d_lfc, d_padj, a_lfc, a_padj), as.numeric)) %>%
    filter(!is.na(d_lfc), !is.na(a_lfc)) %>%
    mutate(
      aging_sig = case_when(
        !is.na(a_padj) & a_padj < 0.05 ~ "Aging sig (padj<0.05)",
        !is.na(a_padj) & a_padj < 0.1  ~ "Aging trending (padj<0.1)",
        TRUE                            ~ "Not sig in aging"
      ),
      aging_sig = factor(aging_sig, levels = sig_levels)
    ) %>%
    arrange(aging_sig)
}

# -- Per-panel plot helper ----------------------------------
make_scatter <- function(dat, panel_label, panel_letter = NULL) {
  aging_s <- dat %>% filter(aging_sig == "Aging sig (padj<0.05)")
  r_sig   <- if (nrow(aging_s) >= 5)
               cor(aging_s$d_lfc, aging_s$a_lfc, method = "pearson")
             else NA_real_

  r_str <- if (is.na(r_sig)) "Pearson r = NA"
           else sprintf("Pearson r = %.3f", r_sig)
  label_r <- sprintf("%s\n(n = %d, ageing padj < 0.05)",
                     r_str, nrow(aging_s))

  title_str <- if (is.null(panel_letter)) panel_label
               else paste0(panel_letter, "   ", panel_label)

  ggplot(dat, aes(x = d_lfc, y = a_lfc,
                  colour = aging_sig, alpha = aging_sig)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey25", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey25", linewidth = 0.8) +
    geom_point(size = 2.0, stroke = 0) +
    { if (nrow(aging_s) >= 5)
        geom_smooth(data = aging_s, aes(x = d_lfc, y = a_lfc),
                    method = "lm", se = TRUE,
                    colour = "#B2182B", fill = "#FDDBC7",
                    linewidth = 1.4, inherit.aes = FALSE) } +
    annotate("label", x = -Inf, y = Inf, label = label_r,
             hjust = -0.05, vjust = 1.05,
             size = 6.0, label.size = 0.5,
             label.padding = unit(0.6, "lines"),
             label.r = unit(0.15, "lines"),
             fill = alpha("white", 0.94),
             colour = "black", family = "sans") +
    scale_colour_manual(values = sig_colors, drop = FALSE, name = NULL) +
    scale_alpha_manual(values = sig_alpha, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(
      title = title_str,
      x = expression("Disease log"[2]*"FC (padj < 0.05)"),
      y = expression("Ageing log"[2]*"FC")
    ) +
    theme_classic(base_size = 20, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "bold", size = 26, hjust = 0,
                                        margin = margin(b = 12),
                                        family = "sans", colour = "black"),
      axis.text          = element_text(size = 18, colour = "black",
                                        family = "sans"),
      axis.title         = element_text(size = 22, family = "sans",
                                        colour = "black"),
      axis.title.y       = element_text(margin = margin(r = 12)),
      axis.title.x       = element_text(margin = margin(t = 10)),
      panel.border       = element_rect(colour = "black", fill = NA,
                                        linewidth = 0.8),
      axis.line          = element_blank(),
      axis.ticks         = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length  = unit(5, "pt"),
      panel.grid.major   = element_line(colour = "grey90", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.justification = "center",
      legend.text        = element_text(size = 20, family = "sans",
                                        colour = "black"),
      legend.key.size    = unit(0.9, "cm"),
      legend.spacing.x   = unit(0.6, "cm"),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(16, 20, 14, 18)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1)))
}

# -- Load each of the 4 panels ------------------------------
dat_by_panel <- list()
for (p in PANELS) {
  dir <- if (p$tissue == "Kidney") KIDNEY_DIR else LUNG_DIR
  pattern <- sprintf("^%s_%s_paired_sig\\.bed$", p$tissue, p$ct)
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    message("[SKIP] no paired BED for ", p$tissue, "_", p$ct, " in ", dir)
    next
  }
  d <- read_paired(files[1]); if (is.null(d)) next
  dat_by_panel[[p$label]] <- d
  message(sprintf("Loaded %s: %d paired peaks (%d ageing-sig)",
                  p$label, nrow(d),
                  sum(d$aging_sig == "Aging sig (padj<0.05)")))
}
if (length(dat_by_panel) == 0) stop("No panel data loaded.")

# -- Build solo + tagged versions ---------------------------
solo_list  <- list()
combo_list <- list()
for (p in PANELS) {
  if (is.null(dat_by_panel[[p$label]])) next
  solo_list[[p$label]]  <- make_scatter(dat_by_panel[[p$label]], p$label)
  combo_list[[p$label]] <- make_scatter(dat_by_panel[[p$label]], p$label,
                                        panel_letter = p$letter)
}

# -- Save individual single-CT PDFs (bigger than before) ---
W_SOLO <- 11
H_SOLO <- 10
for (lbl in names(solo_list)) {
  out_base <- sprintf("Fig15_solo_%s", lbl)
  ggsave(file.path(OUT_DIR, paste0(out_base, ".pdf")),
         solo_list[[lbl]], width = W_SOLO, height = H_SOLO, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(out_base, ".png")),
         solo_list[[lbl]], width = W_SOLO, height = H_SOLO, dpi = 400)
}

# -- Save combined 2x2 panel --------------------------------
if (length(combo_list) >= 1) {
  combined <- wrap_plots(combo_list, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.justification = "center")

  W_COMBO <- 20    # 2 columns x ~10 in
  H_COMBO <- 20    # 2 rows x ~10 in

  ggsave(file.path(OUT_DIR, "Fig15_combined.pdf"),
         combined, width = W_COMBO, height = H_COMBO, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, "Fig15_combined.png"),
         combined, width = W_COMBO, height = H_COMBO, dpi = 400)
}

message("Saved Fig 15 outputs (", length(solo_list),
        " solo + 1 combined)  ->  ", OUT_DIR)

if (interactive()) print(combined)
