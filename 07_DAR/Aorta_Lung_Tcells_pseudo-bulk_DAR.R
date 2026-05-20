# ============================================================
# Pseudo-bulk DAR wrapper for Aorta / Lung / Tcells (DESeq2)
# Uses DATesting.R (apply_DESeq2_test_seurat)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# 0) User paths (EDIT THESE 3)
# -----------------------------
obj_rds_path <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/Integrated/",
  "All_Tissues_Integrated_Annotated_Clean_for_DAR.rds"
)

datesting_r_path <- paste0(
  "/home/s4869245/scripts/DAR/",
  "DATesting.R"
)

base_out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_DESeq2"
)

# -----------------------------
# 1) Global parameters
# -----------------------------
assay_name <- "peaks_universal"
tissue_col <- "Tissue"
celltype_col <- "cell_type"

# For Lung small N: set to 30
min_cells_per_group_default <- 80
min_cells_per_group_by_tissue <- list(
  Aorta = 80,
  Lung  = 50,
  Tcells = 80
)

# Ralph suggested testing exp.thresh
exp_thresh <- 0.05
num_splits <- 10
padj_cutoff <- 0.05

# Which tissues to run
target_tissues <- c("Aorta", "Lung", "Tcells")

# Candidate group columns (auto-pick first valid)
candidate_group_cols <- c(
  "condition", "Group", "group", "deMultliplex2_final_mapped",
  "SampleID", "orig.ident"
)

# If TRUE: for >=3 groups, do reference vs all + pairwise (Kidney-like)
run_pairwise_if_multigroup <- TRUE

# -----------------------------
# 2) Load helper functions
# -----------------------------
if (!file.exists(datesting_r_path)) {
  stop("DATesting.R not found: ", datesting_r_path)
}
source(datesting_r_path)

stopifnot(exists("apply_DESeq2_test_seurat"))
stopifnot(exists("GetExpressedPeaks"))

# -----------------------------
# 3) Load Seurat object
# -----------------------------
if (!file.exists(obj_rds_path)) {
  stop("Object RDS not found: ", obj_rds_path)
}

message("Loading object...")
obj <- readRDS(obj_rds_path)
message("Loaded. Cells = ", ncol(obj))

if (!assay_name %in% Assays(obj)) {
  stop("Assay not found: ", assay_name,
       ". Available: ", paste(Assays(obj), collapse = ", "))
}

if (!tissue_col %in% colnames(obj[[]])) {
  stop("Metadata column not found: ", tissue_col)
}

if (!celltype_col %in% colnames(obj[[]])) {
  stop("Metadata column not found: ", celltype_col)
}

DefaultAssay(obj) <- assay_name

# -----------------------------
# 4) Helpers
# -----------------------------
safe_dir_create <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

peak_to_bed <- function(peaks, out_file) {
  if (length(peaks) == 0) {
    return(invisible(NULL))
  }
  peak_table <- data.frame(
    Chr = sub("(.*)-.*-.*", "\\1", peaks),
    Start = sub(".*-(.*)-.*", "\\1", peaks),
    End = sub(".*-.*-(.*)", "\\1", peaks),
    stringsAsFactors = FALSE
  )
  write.table(
    peak_table,
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

pick_group_col <- function(meta_df, cols) {
  for (cc in cols) {
    if (!cc %in% colnames(meta_df)) next
    v <- meta_df[[cc]]
    v <- as.character(v)
    v <- v[!is.na(v)]
    v <- v[v != ""]
    if (length(v) == 0) next
    if (length(unique(v)) >= 2) return(cc)
  }
  return(NA_character_)
}

choose_reference_level <- function(levels_vec) {
  if ("Control" %in% levels_vec) return("Control")
  if ("Sham" %in% levels_vec) return("Sham")
  if ("Young control" %in% levels_vec) return("Young control")
  levels_vec[1]
}

make_contrasts <- function(levels_vec) {
  lv <- as.character(levels_vec)
  if (length(lv) < 2) return(list())
  if (length(lv) == 2) return(list(c(lv[1], lv[2])))
  ref <- choose_reference_level(lv)
  others <- setdiff(lv, ref)
  out <- lapply(others, function(x) c(x, ref))
  if (run_pairwise_if_multigroup) {
    pairs <- combn(lv, 2, simplify = FALSE)
    for (p in pairs) {
      nm1 <- paste0(p[1], "_vs_", p[2])
      nm2 <- paste0(p[2], "_vs_", p[1])
      already <- vapply(out, function(z) {
        paste0(z[1], "_vs_", z[2]) %in% c(nm1, nm2)
      }, logical(1))
      if (!any(already)) out[[length(out) + 1]] <- p
    }
  }
  out
}

save_hist <- function(x, title_txt, out_png) {
  png(out_png, width = 1600, height = 1200, res = 200)
  hist(x, breaks = 50, main = title_txt, xlab = "p-value")
  dev.off()
}

# -----------------------------
# 5) Run per tissue
# -----------------------------
all_tissue_summaries <- list()

for (tt in target_tissues) {
  message("\n==============================")
  message("Tissue: ", tt)
  message("==============================")
  
  # Output dirs
  out_dir <- file.path(base_out_dir,
                       paste0("DAR_pseudobulk_", tt, "_DESeq2"))
  safe_dir_create(out_dir)
  safe_dir_create(file.path(out_dir, "DAR_tables"))
  safe_dir_create(file.path(out_dir, "DAR_BED"))
  safe_dir_create(file.path(out_dir, "QC"))
  safe_dir_create(file.path(out_dir, "Figures"))
  
  # Tissue subset (cells-based, avoids .data / FetchData issues)
  tissue_vec <- as.character(obj[[tissue_col]][, 1])
  cells_keep <- colnames(obj)[which(tissue_vec == tt)]
  obj_t <- subset(obj, cells = cells_keep)
  DefaultAssay(obj_t) <- assay_name
  
  message("Cells: ", ncol(obj_t))
  
  meta_t <- obj_t[[]]
  
  # pick group column
  group_col <- pick_group_col(meta_t, candidate_group_cols)
  if (is.na(group_col)) {
    message("ERROR: No valid group column found for ", tt)
    message("Tried: ", paste(candidate_group_cols, collapse = ", "))
    next
  }
  
  message("Using group column: ", group_col)
  
  # build dar_group + cell_type_dar
  obj_t$dar_group <- as.character(obj_t[[group_col]][, 1])
  obj_t$cell_type_dar <- as.character(obj_t[[celltype_col]][, 1])
  
  # drop NA groups
  ok <- !is.na(obj_t$dar_group) & obj_t$dar_group != ""
  obj_t <- subset(obj_t, cells = colnames(obj_t)[ok])
  
  # group levels
  group_levels <- sort(unique(obj_t$dar_group))
  obj_t$dar_group <- factor(obj_t$dar_group, levels = group_levels)
  
  # contrasts
  contrasts_list <- make_contrasts(group_levels)
  if (length(contrasts_list) == 0) {
    message("Skip tissue (need >=2 groups): ", tt)
    next
  }
  
  # min cells threshold (Lung=30)
  min_cells_per_group <- min_cells_per_group_default
  if (tt %in% names(min_cells_per_group_by_tissue)) {
    min_cells_per_group <- min_cells_per_group_by_tissue[[tt]]
  }
  message("min_cells_per_group = ", min_cells_per_group)
  
  # QC exports
  write.csv(
    as.data.frame(table(obj_t$dar_group, useNA = "ifany")),
    file.path(out_dir, "QC", paste0(tt, "_dar_group_counts.csv")),
    row.names = FALSE
  )
  
  write.csv(
    as.data.frame(table(obj_t$cell_type_dar, useNA = "ifany")),
    file.path(out_dir, "QC", paste0(tt, "_celltype_counts.csv")),
    row.names = FALSE
  )
  
  write.csv(
    as.data.frame(table(obj_t$cell_type_dar, obj_t$dar_group)),
    file.path(out_dir, "QC",
              paste0(tt, "_celltype_by_group_counts.csv")),
    row.names = FALSE
  )
  
  # eligible cell types: keep those with >= min_cells in >=2 groups
  celltypes_all <- sort(unique(obj_t$cell_type_dar))
  eligible_celltypes <- c()
  for (ct in celltypes_all) {
    ct_obj <- subset(obj_t, subset = cell_type_dar == ct)
    tab <- table(ct_obj$dar_group)
    if (sum(tab >= min_cells_per_group) >= 2) {
      eligible_celltypes <- c(eligible_celltypes, ct)
    }
  }
  
  writeLines(
    eligible_celltypes,
    con = file.path(out_dir, "QC",
                    paste0(tt, "_eligible_celltypes.txt"))
  )
  
  message("Eligible cell types (n=",
          length(eligible_celltypes), "): ",
          paste(head(eligible_celltypes, 20), collapse = ", "))
  
  if (length(eligible_celltypes) == 0) {
    message("Skip tissue (no eligible cell types): ", tt)
    next
  }
  
  # Run loop
  tissue_summary <- list()
  
  for (ct in eligible_celltypes) {
    message("\n-- Cell type: ", ct)
    
    ct_obj <- subset(obj_t, subset = cell_type_dar == ct)
    DefaultAssay(ct_obj) <- assay_name
    Idents(ct_obj) <- "dar_group"
    
    ct_counts <- table(ct_obj$dar_group, useNA = "ifany")
    write.csv(
      as.data.frame(ct_counts),
      file.path(out_dir, "QC",
                paste0(gsub("[ /]", "_", ct), "__group_counts.csv")),
      row.names = FALSE
    )
    
    for (cc in contrasts_list) {
      g1 <- cc[1]
      g2 <- cc[2]
      contrast_name <- paste0(g1, "_vs_", g2)
      
      cells_g1 <- colnames(ct_obj)[Idents(ct_obj) == g1]
      cells_g2 <- colnames(ct_obj)[Idents(ct_obj) == g2]
      
      n1 <- length(cells_g1)
      n2 <- length(cells_g2)
      
      if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
        next
      }
      
      message("   Contrast: ", contrast_name,
              " (", n1, " vs ", n2, ")")
      
      res <- tryCatch(
        {
          apply_DESeq2_test_seurat(
            seurat.object = ct_obj,
            population.1 = cells_g1,
            population.2 = cells_g2,
            exp.thresh = exp_thresh,
            num.splits = num_splits,
            assay.use = assay_name,
            verbose = TRUE
          )
        },
        error = function(e) NULL
      )
      
      if (is.null(res)) next
      
      res_df <- as.data.frame(res)
      res_df$peak <- rownames(res_df)
      
      res_sig <- subset(res_df, !is.na(padj) & padj < padj_cutoff)
      res_open <- subset(res_sig, log2FoldChange > 0)
      res_close <- subset(res_sig, log2FoldChange < 0)
      
      out_label <- paste0(
        gsub("[ /]", "_", ct), "__",
        contrast_name, "__exp",
        gsub("\\.", "", format(exp_thresh))
      )
      
      write.table(
        res_df,
        file = file.path(out_dir, "DAR_tables",
                         paste0(out_label, "_DESeq2_all.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      
      write.table(
        res_sig,
        file = file.path(out_dir, "DAR_tables",
                         paste0(out_label, "_DESeq2_padj005.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      
      write.table(
        res_open,
        file = file.path(out_dir, "DAR_tables",
                         paste0(out_label, "_opening.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      
      write.table(
        res_close,
        file = file.path(out_dir, "DAR_tables",
                         paste0(out_label, "_closing.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      
      peak_to_bed(
        peaks = res_open$peak,
        out_file = file.path(out_dir, "DAR_BED",
                             paste0(gsub("[ /]", "_", ct),
                                    "__", contrast_name, "__opening.bed"))
      )
      
      peak_to_bed(
        peaks = res_close$peak,
        out_file = file.path(out_dir, "DAR_BED",
                             paste0(gsub("[ /]", "_", ct),
                                    "__", contrast_name, "__closing.bed"))
      )
      
      if ("pvalue" %in% colnames(res_df)) {
        save_hist(
          x = res_df$pvalue,
          title_txt = paste(tt, ct, contrast_name, "p-values"),
          out_png = file.path(out_dir, "Figures",
                              paste0(gsub("[ /]", "_", ct),
                                     "__", contrast_name, "__pvalue_hist.png"))
        )
      }
      
      tissue_summary[[paste(ct, contrast_name, sep = "__")]] <-
        data.frame(
          tissue = tt,
          group_col = group_col,
          cell_type = ct,
          contrast = contrast_name,
          group1 = g1,
          group2 = g2,
          n_cells_group1 = n1,
          n_cells_group2 = n2,
          exp_thresh = exp_thresh,
          num_splits = num_splits,
          n_tested_peaks = nrow(res_df),
          n_sig_total = nrow(res_sig),
          n_opening = nrow(res_open),
          n_closing = nrow(res_close),
          stringsAsFactors = FALSE
        )
    }
    
    rm(ct_obj)
    gc()
  }
  
  if (length(tissue_summary) == 0) {
    message("No results produced for tissue: ", tt)
    next
  }
  
  summary_df <- bind_rows(tissue_summary)
  
  write.csv(
    summary_df,
    file.path(out_dir,
              paste0("DAR_pseudobulk_summary_", tt, ".csv")),
    row.names = FALSE
  )
  
  # Long format for plots
  plot_df <- bind_rows(
    summary_df %>% transmute(
      tissue, cell_type, contrast,
      direction = "Opening", n_DAR = n_opening
    ),
    summary_df %>% transmute(
      tissue, cell_type, contrast,
      direction = "Closing", n_DAR = n_closing
    )
  )
  
  write.csv(
    plot_df,
    file.path(out_dir,
              paste0("DAR_pseudobulk_summary_", tt, "_long.csv")),
    row.names = FALSE
  )
  
  # Bar plot (Opening red, Closing blue)
  p_counts <- ggplot(
    plot_df,
    aes(x = cell_type, y = n_DAR, fill = direction)
  ) +
    geom_col(position = "dodge") +
    facet_wrap(~ contrast, nrow = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("Opening" = "#D73027", "Closing" = "#4575B4")
    ) +
    labs(
      title = paste0(tt, " pseudo-bulk DAR counts (DESeq2)"),
      x = "Cell type",
      y = "Number of DARs"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(out_dir, "Figures",
                         paste0(tt, "_pseudobulk_DAR_counts_barplot.pdf")),
    plot = p_counts, width = 14, height = 6
  )
  
  ggsave(
    filename = file.path(out_dir, "Figures",
                         paste0(tt, "_pseudobulk_DAR_counts_barplot.png")),
    plot = p_counts, width = 14, height = 6, dpi = 300
  )
  
  # Tile plots: separate gradients (Opening red, Closing blue)
  p_open <- ggplot(
    subset(plot_df, direction == "Opening"),
    aes(x = contrast, y = cell_type, fill = n_DAR)
  ) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_DAR), size = 3) +
    scale_fill_gradient(low = "white", high = "#D73027") +
    labs(title = paste0(tt, " Opening DARs"),
         x = "Contrast", y = "Cell type") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p_close <- ggplot(
    subset(plot_df, direction == "Closing"),
    aes(x = contrast, y = cell_type, fill = n_DAR)
  ) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_DAR), size = 3) +
    scale_fill_gradient(low = "white", high = "#4575B4") +
    labs(title = paste0(tt, " Closing DARs"),
         x = "Contrast", y = "Cell type") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p_tile <- p_open | p_close
  
  ggsave(
    filename = file.path(out_dir, "Figures",
                         paste0(tt, "_pseudobulk_DAR_counts_tile_AB.pdf")),
    plot = p_tile, width = 14, height = 6.5
  )
  
  ggsave(
    filename = file.path(out_dir, "Figures",
                         paste0(tt, "_pseudobulk_DAR_counts_tile_AB.png")),
    plot = p_tile, width = 14, height = 6.5, dpi = 300
  )
  
  all_tissue_summaries[[tt]] <- summary_df
  message("Finished tissue: ", tt)
}

# Save combined summary across tissues (if any)
if (length(all_tissue_summaries) > 0) {
  all_df <- bind_rows(all_tissue_summaries)
  safe_dir_create(base_out_dir)
  write.csv(
    all_df,
    file.path(base_out_dir, "DAR_pseudobulk_summary_ALL_tissues.csv"),
    row.names = FALSE
  )
  message("\nAll tissues done.")
  message("Combined summary: ",
          file.path(base_out_dir,
                    "DAR_pseudobulk_summary_ALL_tissues.csv"))
} else {
  message("No tissues produced results.")
}
