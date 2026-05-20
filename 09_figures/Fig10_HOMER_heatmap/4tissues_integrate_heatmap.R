#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# Integrated 4-tissue (and 3-tissue, no-Lung) HOMER heatmaps
# Style: matches Cell Metabolism paper panel E/H
#   - column split titles at BOTTOM
#   - left-side AP-1 / CTCF row group labels
#   - fixed motif row order, no row clustering
# ============================================================

# -----------------------------
# 1) Paths
# -----------------------------
motif_file <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Heatmap_ralph/motif_names.txt"
out_dir    <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Integrated_HOMER_Heatmaps"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(motif_file)) stop("motif_names.txt not found: ", motif_file)

# -----------------------------
# 2) Tissue config  (full 4-tissue)
# -----------------------------
cfg_all <- list(
  list(
    tissue          = "Lung",
    contrast        = "Case_vs_Control",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/Lung_DAR_tmp",
    tissue_color    = "#1B9E77",
    preferred_order = c("AT2","B","Ciliated","EC-vasc","Eosinophils",
                        "Fib","Mac-alv","Mac-inter","Mo-Ly6c+","NK","Pen","SMCs","T")
  ),
  list(
    tissue          = "Aorta",
    contrast        = "Challenge_vs_Control",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2",
    tissue_color    = "#D95F02",
    preferred_order = c("Macrophages","Pericytes","SMC","SMCs")
  ),
  list(
    tissue          = "Kidney",
    contrast        = "Day42_vs_Sham",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2",
    tissue_color    = "#7570B3",
    preferred_order = c("DCT","Macrophages","PC","PT","TAL")
  ),
  list(
    tissue          = "Tcell",
    contrast        = "Young_chronic_vs_Young_control",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2",
    tissue_color    = "#E7298A",
    preferred_order = c("Tcell")
  )
)

# 3-tissue version (Lung removed)
cfg_no_lung <- cfg_all[sapply(cfg_all, function(x) x$tissue != "Lung")]

# -----------------------------
# 3) Settings
# -----------------------------
score_cap <- 50

# -----------------------------
# 4) Helper functions
# -----------------------------
short_motif <- function(x) {
  sapply(strsplit(as.character(x), "/"), `[`, 1)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", gsub("%", "", as.character(x)))))
}

normalize_cell <- function(x) sub("_+$", "", x)

parse_dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(NULL)
  data.frame(
    cell_type  = normalize_cell(parts[1]),
    contrast   = parts[2],
    thresh     = parts[3],
    comparison = parts[4],
    stringsAsFactors = FALSE
  )
}

detect_available_cells <- function(base_dir, mode, contrast_keep) {
  homer_dir <- file.path(base_dir,
                         ifelse(mode == "stable", "HOMER_stable_bg", "HOMER_NS_bg"))
  dirs <- list.dirs(homer_dir, recursive = FALSE, full.names = FALSE)
  out <- character(0)
  for (d in dirs) {
    if (d == "logs") next
    meta <- parse_dirname(d)
    if (is.null(meta) || meta$contrast != contrast_keep) next
    out <- c(out, meta$cell_type)
  }
  unique(out)
}

get_cell_order <- function(x, mode) {
  avail <- detect_available_cells(x$base_dir, mode, x$contrast)
  intersect(x$preferred_order, avail)   # strictly respect preferred_order whitelist
}

read_one_homer <- function(fp, motif_order, mode, bg_name,
                           contrast_keep, keep_cells) {
  meta <- parse_dirname(basename(dirname(fp)))
  if (is.null(meta))                          return(NULL)
  if (meta$contrast != contrast_keep)         return(NULL)
  if (!(meta$cell_type %in% keep_cells))      return(NULL)
  
  valid_comps <- if (mode == "stable")
    c("opening_vs_stable","closing_vs_stable","stable_vs_opening","stable_vs_closing")
  else
    c("opening_vs_NS","closing_vs_NS","NS_vs_opening","NS_vs_closing")
  
  if (!(meta$comparison %in% valid_comps)) return(NULL)
  
  tab <- tryCatch(
    read.delim(fp, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || nrow(tab) < 1)           return(NULL)
  if (!("Motif Name" %in% colnames(tab)))       return(NULL)
  if (!("P-value"    %in% colnames(tab)))       return(NULL)
  
  motif_vec <- as.character(tab[["Motif Name"]])
  p_vec     <- to_num(tab[["P-value"]])
  ok        <- !is.na(motif_vec) & !is.na(p_vec)
  motif_vec <- motif_vec[ok]; p_vec <- p_vec[ok]
  if (length(motif_vec) == 0) return(NULL)
  
  best_p <- tapply(p_vec, motif_vec, min, na.rm = TRUE)
  score  <- rep(0, length(motif_order))
  names(score) <- motif_order
  hit  <- intersect(motif_order, names(best_p))
  vals <- -log10(pmax(best_p[hit], 1e-300))
  
  # article-style sign: opening-axis = positive (red), closing-axis = negative (blue)
  if (meta$comparison %in% c(paste0("opening_vs_", bg_name),
                             paste0(bg_name, "_vs_opening"))) {
    vals <-  vals
  } else {
    vals <- -vals
  }
  score[hit] <- vals
  
  data.frame(
    motif_full = motif_order, score = score,
    cell_type = meta$cell_type, comparison = meta$comparison,
    stringsAsFactors = FALSE
  )
}

build_matrix_and_meta <- function(mode, motif_order, cfg) {
  bg_name    <- ifelse(mode == "stable", "stable", "NS")
  bg_pretty  <- ifelse(mode == "stable", "Stable", "NS")
  comp_order <- c(
    paste0("opening_vs_", bg_name), paste0("closing_vs_", bg_name),
    paste0(bg_name, "_vs_opening"), paste0(bg_name, "_vs_closing")
  )
  comp_pretty <- setNames(
    c(paste0("Opening vs\ncell type specific ", bg_pretty),
      paste0("Closing vs\ncell type specific ", bg_pretty),
      paste0(bg_pretty, " cell type\nspecific vs Opening"),
      paste0(bg_pretty, " cell type\nspecific vs Closing")),
    comp_order
  )
  
  score_store <- list(); col_meta <- list()
  
  for (x in cfg) {
    keep_cells <- get_cell_order(x, mode)
    homer_dir  <- file.path(x$base_dir,
                            ifelse(mode == "stable", "HOMER_stable_bg", "HOMER_NS_bg"))
    known_files <- list.files(homer_dir, pattern = "^knownResults\\.txt$",
                              recursive = TRUE, full.names = TRUE)
    if (length(known_files) == 0) next
    
    for (fp in known_files) {
      tmp <- read_one_homer(fp, motif_order, mode, bg_name, x$contrast, keep_cells)
      if (is.null(tmp)) next
      meta   <- parse_dirname(basename(dirname(fp)))
      col_id <- paste(x$tissue, meta$cell_type, meta$comparison, sep = "|")
      score_store[[col_id]] <- tmp$score
      col_meta[[col_id]] <- data.frame(
        col_id      = col_id,
        tissue      = x$tissue,
        contrast    = x$contrast,
        cell_type   = meta$cell_type,
        comparison  = meta$comparison,
        comp_pretty = comp_pretty[meta$comparison],
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(score_store) == 0)
    stop("No matching HOMER results found for mode = ", mode)
  
  # fixed column order: comparison block → tissue → cell type
  ordered_meta <- list()
  for (comp_use in comp_order) {
    for (x in cfg) {
      for (cell_use in get_cell_order(x, mode)) {
        col_id <- paste(x$tissue, cell_use, comp_use, sep = "|")
        if (!is.null(col_meta[[col_id]])) ordered_meta[[col_id]] <- col_meta[[col_id]]
      }
    }
  }
  
  meta_df <- do.call(rbind, ordered_meta)
  mat <- matrix(0, nrow = length(motif_order), ncol = nrow(meta_df),
                dimnames = list(motif_order, meta_df$col_id))
  for (cid in meta_df$col_id) mat[, cid] <- score_store[[cid]]
  
  mat[mat >  score_cap] <-  score_cap
  mat[mat < -score_cap] <- -score_cap
  mat[!is.finite(mat)]  <- 0
  
  list(mat = mat, meta = meta_df)
}

# ---------------------------------------------------------------
# Row mark annotation (AP-1 / CTCF) — does NOT reorder rows
# Uses anno_mark() to draw connecting lines to labels
# ---------------------------------------------------------------
make_row_mark_anno <- function(row_labels_vec) {
  ap1_pat  <- "AP-1|BATF|^JUN|JUNB|JUND|^FOS$|FOSL|IRF4|IRF8|NFAT|ETS|bZIP"
  ctcf_pat <- "CTCF|BORIS|YY1|CCCTC"
  
  ap1_idx  <- grep(ap1_pat,  row_labels_vec, ignore.case = TRUE)
  ctcf_idx <- grep(ctcf_pat, row_labels_vec, ignore.case = TRUE)
  
  # pick representative index for the label position (middle of each group)
  mark_idx    <- c(ap1_idx[ceiling(length(ap1_idx)/2)],
                   ctcf_idx[ceiling(length(ctcf_idx)/2)])
  mark_labels <- c("AP-1", "CTCF")
  
  # remove any NA (group not found)
  ok          <- !is.na(mark_idx)
  list(at = mark_idx[ok], labels = mark_labels[ok])
}

# ---------------------------------------------------------------
# Draw & save one heatmap  (paper style)
# ---------------------------------------------------------------
draw_one_heatmap <- function(obj, mode, cfg, out_dir, suffix = "") {
  mat     <- obj$mat
  meta_df <- obj$meta
  
  tissue_colors <- setNames(
    sapply(cfg, `[[`, "tissue_color"),
    sapply(cfg, `[[`, "tissue")
  )
  
  col_fun_signed <- colorRamp2(
    c(-score_cap, -10, 0, 10, score_cap),
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
  )
  
  # ---- top annotation: tissue color bar ----
  ha_top <- HeatmapAnnotation(
    Tissue = meta_df$tissue,
    col    = list(Tissue = tissue_colors),
    annotation_name_gp  = gpar(fontsize = 8, fontface = "bold"),
    annotation_name_side = "right",
    simple_anno_size     = unit(3, "mm"),
    show_legend          = TRUE,
    show_annotation_name = TRUE
  )
  
  # ---- column split: comparison blocks, labels at BOTTOM ----
  split_levels <- unique(meta_df$comp_pretty)
  split_factor <- factor(meta_df$comp_pretty, levels = split_levels)
  
  # ---- row labels (fixed order from motif_names.txt, never reordered) ----
  row_labels_short <- short_motif(rownames(mat))
  
  # ---- column labels: Tissue_CellType, but skip prefix if tissue == cell_type ----
  col_labels <- ifelse(
    meta_df$tissue == meta_df$cell_type,
    meta_df$tissue,
    paste0(meta_df$tissue, "_", meta_df$cell_type)
  )
  
  ht <- Heatmap(
    mat,
    name              = "signed\n-log10(P)",
    col               = col_fun_signed,
    top_annotation    = ha_top,
    
    # rows — strict fixed order from motif_names.txt
    cluster_rows      = FALSE,
    show_row_dend     = FALSE,
    show_row_names    = FALSE,
    
    # columns
    cluster_columns   = FALSE,
    show_column_dend  = FALSE,
    column_split      = split_factor,
    column_gap        = unit(1.5, "mm"),
    column_title_side = "top",
    column_title_gp   = gpar(fontsize = 8, fontface = "bold"),
    show_column_names = FALSE,
    
    heatmap_legend_param = list(
      title     = "signed -log10(P)",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      direction = "vertical"
    ),
    use_raster = TRUE
  )
  
  tag <- paste0(ifelse(suffix == "", "", paste0(gsub(" ", "_", trimws(suffix)), "_")), mode)
  
  # save outputs
  write.csv(mat,
            file.path(out_dir, paste0("Integrated_", tag, "_matrix.csv")))
  write.csv(meta_df,
            file.path(out_dir, paste0("Integrated_", tag, "_column_metadata.csv")),
            row.names = FALSE)
  
  pdf(file.path(out_dir, paste0("Integrated_", tag, "_heatmap.pdf")),
      width = 20, height = 16)
  draw(ht, merge_legend = TRUE,
       padding = unit(c(5, 120, 5, 2), "mm"))
  dev.off()

  png(file.path(out_dir, paste0("Integrated_", tag, "_heatmap.png")),
      width = 6500, height = 3500, res = 250)
  draw(ht, merge_legend = TRUE,
       padding = unit(c(5, 150, 5, 2), "mm"))
  dev.off()
  
  message("Saved: ", tag)
}

# -----------------------------
# 5) Run
# -----------------------------
motif_order <- readLines(motif_file)
motif_order <- motif_order[nzchar(motif_order)]

# ---- 4-tissue (all) ----
message("=== 4-tissue Stable ===")
stable_obj <- build_matrix_and_meta("stable", motif_order, cfg_all)
draw_one_heatmap(stable_obj, "Stable", cfg_all, out_dir, suffix = "")

message("=== 4-tissue NS ===")
ns_obj <- build_matrix_and_meta("NS", motif_order, cfg_all)
draw_one_heatmap(ns_obj, "NS", cfg_all, out_dir, suffix = "")

# ---- 3-tissue (no Lung) ----
message("=== 3-tissue (no Lung) Stable ===")
stable_no_lung <- build_matrix_and_meta("stable", motif_order, cfg_no_lung)
draw_one_heatmap(stable_no_lung, "Stable", cfg_no_lung, out_dir, suffix = " no_Lung")

message("=== 3-tissue (no Lung) NS ===")
ns_no_lung <- build_matrix_and_meta("NS", motif_order, cfg_no_lung)
draw_one_heatmap(ns_no_lung, "NS", cfg_no_lung, out_dir, suffix = " no_Lung")

message("Done. Output dir: ", out_dir)

