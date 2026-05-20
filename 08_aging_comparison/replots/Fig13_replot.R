#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

AGING_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/aging_DARs"

files <- list.files(AGING_DIR, pattern = "_Aged_vs_Young_DAR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No aging DAR TSV files found in: ", AGING_DIR)

counts <- lapply(files, function(f) {
  bname  <- basename(f)
  tissue <- sub("_.*", "", bname)
  ct     <- sub(paste0("^", tissue, "_"), "", sub("_Aged_vs_Young_DAR\\.tsv$", "", bname))
  tab    <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  sig <- tab[!is.na(padj) & padj < 0.05]
  data.frame(tissue = tissue, cell_type = ct,
             n_opening = sum(sig$log2FoldChange > 0, na.rm = TRUE),
             n_closing = sum(sig$log2FoldChange < 0, na.rm = TRUE),
             stringsAsFactors = FALSE)
}) |> bind_rows()

write.csv(counts, file.path(AGING_DIR, "aging_DAR_counts_summary.csv"), row.names = FALSE)

plot_df <- counts |>
  tidyr::pivot_longer(c(n_opening, n_closing),
                      names_to = "direction", values_to = "n_DAR") |>
  mutate(direction = ifelse(direction == "n_opening", "Opening", "Closing"))

dir_colors <- c(Opening = "#B2182B", Closing = "#2166AC")

make_panel <- function(df, tis, panel_label) {
  df <- df |> filter(tissue == tis)
  if (nrow(df) == 0) return(ggplot() + labs(title = panel_label) + theme_void())

  # Remove cell types with 0 DARs in both directions, sort by total descending
  totals <- df |>
    group_by(cell_type) |>
    summarise(total = sum(n_DAR), .groups = "drop") |>
    filter(total >= 50) |>
    arrange(total)
  if (nrow(totals) == 0) return(ggplot() + labs(title = panel_label) + theme_void())

  df <- df |>
    filter(cell_type %in% totals$cell_type) |>
    mutate(cell_type = factor(cell_type, levels = totals$cell_type))

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
    title    = "Fig. 13  Aging DARs from Science paper (Aged vs Young, padj < 0.05)",
    subtitle = "Pseudo-bulk DESeq2 | non-zero in >=5% of samples",
    theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                  plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
  )

ggsave(file.path(AGING_DIR, "Fig13_aging_DAR_counts.pdf"), fig13, width = 10, height = 9)
ggsave(file.path(AGING_DIR, "Fig13_aging_DAR_counts.png"), fig13, width = 10, height = 9, dpi = 300)
message("Saved: Fig13_aging_DAR_counts (.pdf + .png)  →  ", AGING_DIR)
