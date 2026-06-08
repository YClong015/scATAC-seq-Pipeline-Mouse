#!/usr/bin/env Rscript

library(ComplexHeatmap)
library(circlize)
library(grid)

# Integrated 4-tissue HOMER heatmaps, Cell Metabolism panel E/H style:
# column split titles at bottom, left-side AP-1 / CTCF row group labels,
# fixed motif row order, no row clustering.

## Paths
motif.file <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Heatmap_ralph/motif_names.txt"
out.dir    <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Integrated_HOMER_Heatmaps"

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(motif.file)) stop("motif_names.txt not found: ", motif.file)

## Tissue config (full 4-tissue)
cfg.all <- list(
  list(
    tissue          = "Lung",
    contrast        = "Case_vs_Control",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2",
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
    # v5 11-type, Day42_vs_Sham only (IC + DTL_ATL excluded)
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Kidney_v5_DESeq2",
    tissue_color    = "#7570B3",
    # 8 CTs (EC dropped to match Opening dotplot ordering)
    preferred_order = c("PCT","PST","Injured_PT","TAL","DCT_CNT","PC_URO","LEUK","FIB")
  ),
  list(
    tissue          = "Tcell",
    contrast        = "Young_chronic_vs_Young_control",
    base_dir        = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2",
    tissue_color    = "#E7298A",
    preferred_order = c("Tcell")
  )
)

## Settings
score.cap <- 50

## Helper functions
short.motif <- function(x) {
  sapply(strsplit(as.character(x), "/"), `[`, 1)
}

to.num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", gsub("%", "", as.character(x)))))
}

normalize.cell <- function(x) sub("_+$", "", x)

parse.dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(NULL)
  data.frame(
    cell_type  = normalize.cell(parts[1]),
    contrast   = parts[2],
    thresh     = parts[3],
    comparison = parts[4],
    stringsAsFactors = FALSE
  )
}

detect.available.cells <- function(base_dir, mode, contrast_keep) {
  homer.dir <- file.path(base_dir,
                         ifelse(mode == "stable", "HOMER_stable_bg", "HOMER_NS_bg"))
  dirs <- list.dirs(homer.dir, recursive = FALSE, full.names = FALSE)
  out <- character(0)
  for (d in dirs) {
    if (d == "logs") next
    meta <- parse.dirname(d)
    if (is.null(meta) || meta$contrast != contrast_keep) next
    out <- c(out, meta$cell_type)
  }
  unique(out)
}

get.cell.order <- function(x, mode) {
  avail <- detect.available.cells(x$base_dir, mode, x$contrast)
  intersect(x$preferred_order, avail)   # strictly respect preferred_order whitelist
}

read.one.homer <- function(fp, motif.order, mode, bg.name,
                           contrast_keep, keep.cells) {
  meta <- parse.dirname(basename(dirname(fp)))
  if (is.null(meta))                          return(NULL)
  if (meta$contrast != contrast_keep)         return(NULL)
  if (!(meta$cell_type %in% keep.cells))      return(NULL)

  valid.comps <- if (mode == "stable")
    c("opening_vs_stable","closing_vs_stable","stable_vs_opening","stable_vs_closing")
  else
    c("opening_vs_NS","closing_vs_NS","NS_vs_opening","NS_vs_closing")

  if (!(meta$comparison %in% valid.comps)) return(NULL)

  tab <- tryCatch(
    read.delim(fp, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || nrow(tab) < 1)           return(NULL)
  if (!("Motif Name" %in% colnames(tab)))       return(NULL)
  if (!("P-value"    %in% colnames(tab)))       return(NULL)

  motif.vec <- as.character(tab[["Motif Name"]])
  p.vec     <- to.num(tab[["P-value"]])
  ok        <- !is.na(motif.vec) & !is.na(p.vec)
  motif.vec <- motif.vec[ok]; p.vec <- p.vec[ok]
  if (length(motif.vec) == 0) return(NULL)

  best.p <- tapply(p.vec, motif.vec, min, na.rm = TRUE)
  score  <- rep(0, length(motif.order))
  names(score) <- motif.order
  hit  <- intersect(motif.order, names(best.p))
  vals <- -log10(pmax(best.p[hit], 1e-300))

  # article-style sign: opening-axis = positive (red), closing-axis = negative (blue)
  if (meta$comparison %in% c(paste0("opening_vs_", bg.name),
                             paste0(bg.name, "_vs_opening"))) {
    vals <-  vals
  } else {
    vals <- -vals
  }
  score[hit] <- vals

  data.frame(
    motif_full = motif.order, score = score,
    cell_type = meta$cell_type, comparison = meta$comparison,
    stringsAsFactors = FALSE
  )
}

build.matrix.and.meta <- function(mode, motif.order, cfg) {
  bg.name    <- ifelse(mode == "stable", "stable", "NS")
  bg.pretty  <- ifelse(mode == "stable", "Stable", "NS")
  comp.order <- c(
    paste0("opening_vs_", bg.name), paste0("closing_vs_", bg.name),
    paste0(bg.name, "_vs_opening"), paste0(bg.name, "_vs_closing")
  )
  comp_pretty <- setNames(
    c(paste0("Opening vs\ncell type specific ", bg.pretty),
      paste0("Closing vs\ncell type specific ", bg.pretty),
      paste0(bg.pretty, " cell type\nspecific vs Opening"),
      paste0(bg.pretty, " cell type\nspecific vs Closing")),
    comp.order
  )

  score.store <- list(); col_meta <- list()

  for (x in cfg) {
    keep.cells <- get.cell.order(x, mode)
    homer.dir  <- file.path(x$base_dir,
                            ifelse(mode == "stable", "HOMER_stable_bg", "HOMER_NS_bg"))
    known.files <- list.files(homer.dir, pattern = "^knownResults\\.txt$",
                              recursive = TRUE, full.names = TRUE)
    if (length(known.files) == 0) next

    for (fp in known.files) {
      tmp <- read.one.homer(fp, motif.order, mode, bg.name, x$contrast, keep.cells)
      if (is.null(tmp)) next
      meta   <- parse.dirname(basename(dirname(fp)))
      col_id <- paste(x$tissue, meta$cell_type, meta$comparison, sep = "|")
      score.store[[col_id]] <- tmp$score
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

  if (length(score.store) == 0)
    stop("No matching HOMER results found for mode = ", mode)

  # fixed column order: comparison block -> tissue -> cell type
  ordered.meta <- list()
  for (comp_use in comp.order) {
    for (x in cfg) {
      for (cell_use in get.cell.order(x, mode)) {
        col_id <- paste(x$tissue, cell_use, comp_use, sep = "|")
        if (!is.null(col_meta[[col_id]])) ordered.meta[[col_id]] <- col_meta[[col_id]]
      }
    }
  }

  meta.df <- do.call(rbind, ordered.meta)
  mat <- matrix(0, nrow = length(motif.order), ncol = nrow(meta.df),
                dimnames = list(motif.order, meta.df$col_id))
  for (cid in meta.df$col_id) mat[, cid] <- score.store[[cid]]

  mat[mat >  score.cap] <-  score.cap
  mat[mat < -score.cap] <- -score.cap
  mat[!is.finite(mat)]  <- 0

  list(mat = mat, meta = meta.df)
}

## Row mark annotation (AP-1 / CTCF); does not reorder rows
make.row.mark.anno <- function(row_labels_vec) {
  ap1.pat  <- "AP-1|BATF|^JUN|JUNB|JUND|^FOS$|FOSL|IRF4|IRF8|NFAT|ETS|bZIP"
  ctcf.pat <- "CTCF|BORIS|YY1|CCCTC"

  ap1.idx  <- grep(ap1.pat,  row_labels_vec, ignore.case = TRUE)
  ctcf.idx <- grep(ctcf.pat, row_labels_vec, ignore.case = TRUE)

  # pick representative index for the label position (middle of each group)
  mark.idx    <- c(ap1.idx[ceiling(length(ap1.idx)/2)],
                   ctcf.idx[ceiling(length(ctcf.idx)/2)])
  mark.labels <- c("AP-1", "CTCF")

  ok          <- !is.na(mark.idx)
  list(at = mark.idx[ok], labels = mark.labels[ok])
}

## Draw and save one heatmap (paper style)
draw.one.heatmap <- function(obj, mode, cfg, out.dir, suffix = "") {
  mat     <- obj$mat
  meta.df <- obj$meta

  tissue.colors <- setNames(
    sapply(cfg, `[[`, "tissue_color"),
    sapply(cfg, `[[`, "tissue")
  )

  col.fun.signed <- colorRamp2(
    c(-score.cap, -10, 0, 10, score.cap),
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
  )

  ## top annotation: tissue color bar
  ha.top <- HeatmapAnnotation(
    Tissue = meta.df$tissue,
    col    = list(Tissue = tissue.colors),
    annotation_name_gp  = gpar(fontsize = 8, fontface = "bold"),
    annotation_name_side = "right",
    simple_anno_size     = unit(3, "mm"),
    show_legend          = TRUE,
    show_annotation_name = TRUE
  )

  ## column split: comparison blocks, labels at bottom
  split.levels <- unique(meta.df$comp_pretty)
  split.factor <- factor(meta.df$comp_pretty, levels = split.levels)

  ## row labels (fixed order from motif_names.txt, never reordered)
  row.labels.short <- short.motif(rownames(mat))

  # column labels: Tissue_CellType, but skip prefix if tissue == cell_type
  col.labels <- ifelse(
    meta.df$tissue == meta.df$cell_type,
    meta.df$tissue,
    paste0(meta.df$tissue, "_", meta.df$cell_type)
  )

  ht <- Heatmap(
    mat,
    name              = "signed\n-log10(P)",
    col               = col.fun.signed,
    top_annotation    = ha.top,

    # rows: strict fixed order from motif_names.txt
    cluster_rows      = FALSE,
    show_row_dend     = FALSE,
    show_row_names    = FALSE,

    # columns
    cluster_columns   = FALSE,
    show_column_dend  = FALSE,
    column_split      = split.factor,
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

  write.csv(mat,
            file.path(out.dir, paste0("Integrated_", tag, "_matrix.csv")))
  write.csv(meta.df,
            file.path(out.dir, paste0("Integrated_", tag, "_column_metadata.csv")),
            row.names = FALSE)

  pdf(file.path(out.dir, paste0("Integrated_", tag, "_heatmap.pdf")),
      width = 20, height = 16)
  draw(ht, merge_legend = TRUE,
       padding = unit(c(5, 120, 5, 2), "mm"))
  dev.off()

  png(file.path(out.dir, paste0("Integrated_", tag, "_heatmap.png")),
      width = 6500, height = 3500, res = 250)
  draw(ht, merge_legend = TRUE,
       padding = unit(c(5, 150, 5, 2), "mm"))
  dev.off()

  message("Saved: ", tag)
}

## Run
motif.order <- readLines(motif.file)
motif.order <- motif.order[nzchar(motif.order)]

message("=== 4-tissue Stable ===")
stable.obj <- build.matrix.and.meta("stable", motif.order, cfg.all)
draw.one.heatmap(stable.obj, "Stable", cfg.all, out.dir, suffix = "")

message("=== 4-tissue NS ===")
ns.obj <- build.matrix.and.meta("NS", motif.order, cfg.all)
draw.one.heatmap(ns.obj, "NS", cfg.all, out.dir, suffix = "")

message("Done. Output dir: ", out.dir)
