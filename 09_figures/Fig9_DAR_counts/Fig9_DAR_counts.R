#!/usr/bin/env Rscript
# ============================================================
# Fig 9. DAR quantification across all four disease tissues
# All four panels use tile heatmaps (signed n_DARs)
# Layout: 2×2 grid (Kidney / Lung / Aorta / Tcells)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------
# 1) Paths
# ---------------------------------------------------------------
base_dar <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2"

dar_dirs <- list(
  Kidney = file.path(base_dar, "DAR_pseudobulk_Kidney_DESeq2", "DAR_tables"),
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Lung_DESeq2/DAR_tables",
  Aorta  = file.path(base_dar, "DAR_pseudobulk_Aorta_DESeq2",  "DAR_tables"),
  Tcells = file.path(base_dar, "DAR_pseudobulk_Tcells_DESeq2", "DAR_tables")
)

# Preferred cell-type and contrast orders per tissue
tissue_cfg <- list(
  Kidney = list(
    panel_label    = "A  Kidney",
    contrast_order = "Day42_vs_Sham",
    celltype_order = c("DCT", "Endothelial", "IC", "Macrophages",
                       "PC", "Podocytes", "PT", "TAL")
  ),
  Lung = list(
    panel_label    = "B  Lung",
    contrast_order = "Case_vs_Control",
    celltype_order = c("AT2", "B", "Ciliated", "EC-vasc", "Eosinophils",
                       "Fib", "Mac-alv", "Mac-inter", "Mo-Ly6c+",
                       "NK", "Pen", "SMCs", "T")
  ),
  Aorta = list(
    panel_label    = "C  Aorta",
    contrast_order = "Challenge_vs_Control",
    celltype_order = c("Macrophages", "Pericytes", "SMC")
  ),
  Tcells = list(
    panel_label    = "D  T cells",
    contrast_order = "Young_chronic_vs_Young_control",
    celltype_order = "Tcell"
  )
)

out_dir <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Fig9_DAR_counts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# 2) Helper: parse cell_type + contrast from filename
#    Supports both naming conventions:
#      Kidney:          CellType__Contrast__005_opening.tsv
#      Aorta/Lung/Tcell: CellType__Contrast__exp005_opening.tsv
# ---------------------------------------------------------------
parse_stem <- function(fname) {
  stem  <- tools::file_path_sans_ext(fname)
  parts <- strsplit(stem, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  list(cell_type = parts[1], contrast = parts[2], suffix = tolower(parts[3]))
}

# ---------------------------------------------------------------
# 3) Count DARs — strategy:
#    a) If *_opening.tsv / *_closing.tsv exist → count rows directly
#       (matches original pipeline: padj < 0.05, log2FC > 0 / < 0, no lfc_cut)
#    b) Fallback to *_DESeq2_padj005.tsv → split by log2FC sign
#       (also no lfc_cut, to match the original definition)
# ---------------------------------------------------------------
count_dars <- function(dar_dir, tissue_label) {
  all_files <- list.files(dar_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(all_files) == 0) {
    warning("No TSV files found in: ", dar_dir)
    return(NULL)
  }

  # --- path (a): pre-split opening/closing files ---
  open_files  <- all_files[grepl("_opening\\.tsv$",  basename(all_files))]
  close_files <- all_files[grepl("_closing\\.tsv$",  basename(all_files))]

  read_presplit <- function(fps, direction) {
    lapply(fps, function(fp) {
      meta <- parse_stem(basename(fp))
      if (is.null(meta)) return(NULL)
      if (!grepl(direction, meta$suffix)) return(NULL)
      n <- tryCatch(
        nrow(data.table::fread(fp, select = 1L, showProgress = FALSE)),
        error = function(e) NA_integer_
      )
      if (is.na(n) || n == 0) return(NULL)
      data.frame(cell_type = meta$cell_type, contrast = meta$contrast,
                 direction = if (direction == "opening") "Opening" else "Closing",
                 n_DARs = n, stringsAsFactors = FALSE)
    })
  }

  if (length(open_files) > 0 || length(close_files) > 0) {
    rows <- c(read_presplit(open_files, "opening"),
              read_presplit(close_files, "closing"))
    df <- bind_rows(Filter(Negate(is.null), rows))
    if (!is.null(df) && nrow(df) > 0)
      return(df %>% mutate(tissue = tissue_label))
  }

  # --- path (b): fallback — read DESeq2_padj005.tsv, split by log2FC sign ---
  padj_files <- all_files[grepl("DESeq2_padj005\\.tsv$", basename(all_files))]
  if (length(padj_files) == 0) {
    warning("No opening/closing or padj005 files in: ", dar_dir)
    return(NULL)
  }

  rows <- lapply(padj_files, function(fp) {
    meta <- parse_stem(basename(fp))
    if (is.null(meta)) return(NULL)

    tab <- tryCatch(
      data.table::fread(fp, showProgress = FALSE),
      error = function(e) NULL
    )
    if (is.null(tab) || nrow(tab) == 0) return(NULL)

    lfc_col <- intersect(c("log2FoldChange", "log2FC", "avg_log2FC"),
                         colnames(tab))[1]
    if (is.na(lfc_col)) return(NULL)

    lfc     <- tab[[lfc_col]]
    n_open  <- sum(lfc > 0, na.rm = TRUE)   # no lfc_cut — matches original
    n_close <- sum(lfc < 0, na.rm = TRUE)

    rows_out <- list()
    if (n_open  > 0) rows_out[[1]] <- data.frame(
      cell_type = meta$cell_type, contrast = meta$contrast,
      direction = "Opening", n_DARs = n_open, stringsAsFactors = FALSE)
    if (n_close > 0) rows_out[[2]] <- data.frame(
      cell_type = meta$cell_type, contrast = meta$contrast,
      direction = "Closing", n_DARs = n_close, stringsAsFactors = FALSE)
    bind_rows(rows_out)
  })

  df <- bind_rows(Filter(Negate(is.null), rows))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% mutate(tissue = tissue_label)
}

# ---------------------------------------------------------------
# 4) Tile heatmap — faceted Closing | Opening
#    Color intensity = n DARs; numbers inside tiles in black (always readable)
# ---------------------------------------------------------------
make_tile <- function(df, cfg) {

  df <- df %>%
    filter(contrast  %in% cfg$contrast_order,
           cell_type %in% cfg$celltype_order)

  if (nrow(df) == 0) return(ggplot() + labs(title = cfg$panel_label) + theme_void())

  ct_order <- intersect(cfg$celltype_order, unique(df$cell_type))

  df <- df %>%
    mutate(
      cell_type = factor(cell_type, levels = ct_order),
      direction = factor(direction, levels = c("Closing", "Opening")),
      signed_n  = ifelse(direction == "Opening", n_DARs, -n_DARs),
      label     = scales::comma(n_DARs, accuracy = 1)
    )

  abs_max <- max(abs(df$signed_n))

  ggplot(df, aes(x = contrast, y = cell_type, fill = signed_n)) +
    geom_tile(color = "white", linewidth = 0.6) +
    facet_wrap(~ direction, nrow = 1) +
    geom_text(aes(label = label),
              color = "black", size = 4.2, fontface = "bold") +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6191B",
      midpoint = 0,
      limits   = c(-abs_max, abs_max),
      labels   = function(x) scales::comma(abs(x)),
      name     = "n DAR\n(closing<0\nopening>0)"
    ) +
    scale_y_discrete(limits = rev(ct_order)) +
    labs(title = cfg$panel_label, x = NULL, y = "Cell type") +
    theme_bw(base_size = 13) +
    theme(
      plot.title        = element_text(face = "bold", size = 14, hjust = 0),
      axis.text.x       = element_text(size = 11, angle = 35, hjust = 1,
                                       colour = "grey20"),
      axis.text.y       = element_text(size = 12, colour = "grey10"),
      axis.title.y      = element_text(size = 12),
      strip.text        = element_text(size = 13, face = "bold"),
      strip.background  = element_rect(fill = "grey92", colour = "grey60"),
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 10),
      legend.key.height = unit(1.2, "cm"),
      panel.grid        = element_blank()
    )
}

# ---------------------------------------------------------------
# 5) Load data and build panels
# ---------------------------------------------------------------
panels <- lapply(names(dar_dirs), function(tis) {
  df <- count_dars(dar_dirs[[tis]], tis)
  if (is.null(df)) {
    message("WARNING: No data for ", tis)
    return(ggplot() +
             labs(title = paste0(tissue_cfg[[tis]]$panel_label,
                                 " — data not found")) +
             theme_void())
  }
  # Save per-tissue count table
  write.csv(df,
            file.path(out_dir, paste0(tis, "_DAR_counts_summary.csv")),
            row.names = FALSE)
  make_tile(df, tissue_cfg[[tis]])
})
names(panels) <- names(dar_dirs)

# ---------------------------------------------------------------
# 6) Assemble 2×2
# ---------------------------------------------------------------
fig9 <- (panels$Kidney | panels$Lung) /
        (panels$Aorta  | panels$Tcells) +
  plot_layout(heights = c(3, 1)) +   # Kidney/Lung taller; Aorta/Tcell compact
  plot_annotation(
    title = "Pseudo-bulk DAR counts across four disease tissues (DESeq2)",
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                            hjust = 0.5))
  )

# ---------------------------------------------------------------
# 7) Save
# ---------------------------------------------------------------
ggsave(file.path(out_dir, "Fig9_DAR_counts.pdf"),
       fig9, width = 20, height = 14)

ggsave(file.path(out_dir, "Fig9_DAR_counts.png"),
       fig9, width = 20, height = 14, dpi = 300)

message("Saved: Fig9_DAR_counts.pdf / .png  →  ", out_dir)
