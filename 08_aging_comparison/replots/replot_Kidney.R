#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Kidney"
FIG_DIR <- file.path(OUT_DIR, "figures")

fisher <- read.csv(file.path(OUT_DIR, "overlap_stats.csv"), stringsAsFactors = FALSE)
r_df   <- read.csv(file.path(OUT_DIR, "scatter_r_stats.csv"), stringsAsFactors = FALSE)

sig_label <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                   p < 0.05  ~ "*",   TRUE       ~ "ns")
}

ct_order <- paste0("Kidney\n", c("DCT", "PT", "TAL"))
dir_colors <- c("Opening" = "#B2182B", "Closing" = "#2166AC")

fmt_pval <- function(p) {
  dplyr::case_when(
    is.na(p) | p >= 1   ~ "",
    p == 0              ~ "p<1e-300",
    p < 0.001           ~ paste0("p=", formatC(p, format = "e", digits = 1)),
    TRUE                ~ paste0("p=", round(p, 3))
  )
}

# ── Panel A: Fisher OR (OR >= 1 only) ─────────────────────
fisher_long <- fisher %>%
  mutate(label     = paste0(tissue, "\n", cell_type),
         sig_open  = sig_label(pval_open),
         sig_close = sig_label(pval_close)) %>%
  select(label, tissue, fold_open, fold_close, sig_open, sig_close,
         pval_open, pval_close) %>%
  pivot_longer(cols = c(fold_open, fold_close),
               names_to = "direction", values_to = "OR") %>%
  mutate(
    direction  = ifelse(direction == "fold_open", "Opening", "Closing"),
    sig        = ifelse(direction == "Opening", sig_open, sig_close),
    pval       = ifelse(direction == "Opening", pval_open, pval_close),
    pval_label = fmt_pval(pval),
    label      = factor(label, levels = ct_order)
  ) %>%
  filter(is.na(OR) | OR >= 1) %>%          # ← remove OR < 1
  mutate(OR_plot = pmin(OR, 12, na.rm = TRUE))

pA <- ggplot(fisher_long, aes(x = label, y = OR_plot, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = sig, group = direction),
            position = position_dodge(0.8), vjust = -0.3, size = 4, fontface = "bold") +
  geom_text(aes(label = pval_label, group = direction),
            position = position_dodge(0.8), vjust = -1.6, size = 2.8, color = "grey30") +
  scale_fill_manual(values = dir_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "A  Enrichment of CKD DARs in aging DARs (Kidney)",
       subtitle = "OR < 1 (depletion) not shown",
       x = NULL, y = "Enrichment fold (OR)", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title  = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = "top")

# ── Panel B: Pearson r ─────────────────────────────────────
r_long <- r_df %>%
  mutate(label = factor(paste0(tissue, "\n", cell_type), levels = ct_order))

pB <- ggplot(r_long, aes(x = label, y = r_aging_sig)) +
  geom_col(fill = "#4D4D4D", width = 0.55) +
  geom_text(aes(label = sprintf("n=%d", n_aging_sig)),
            vjust = -0.4, size = 3, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  scale_y_continuous(limits = c(-0.1, 0.65),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "B  Correlation of log2FC (disease vs aging, Kidney)",
       x = NULL, y = "Pearson r (aging-significant peaks)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title  = element_text(face = "bold", size = 11))

combined <- pA / pB +
  plot_annotation(
    title    = "Kidney: CKD DARs overlap with aging DARs",
    subtitle = "Fisher's exact test (one-sided, enrichment only) | Pearson r on aging-significant peaks",
    theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                  plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"))
  )

ggsave(file.path(FIG_DIR, "Aim3_Kidney_OR1filter.pdf"),
       combined, width = 8, height = 9)
ggsave(file.path(FIG_DIR, "Aim3_Kidney_OR1filter.png"),
       combined, width = 8, height = 9, dpi = 300)
message("Saved: Aim3_Kidney_OR1filter (.pdf + .png)")
