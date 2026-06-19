# Kidney pseudo-bulk DAR wrapper (DESeq2 via DATesting.R)
# Based on previous workflow

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
})

## User paths (EDIT THESE)
# Per-tissue universal-peak object.
obj_rds_path <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/Kidney/",
  "kidney_merged_universal.rds"
)

datesting_r_path <- paste0(
  "/home/s4869245/scripts/DAR/",
  "DATesting.R"
)

out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_Kidney_DESeq2"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DAR_tables"),
           showWarnings = FALSE)
dir.create(file.path(out_dir, "DAR_BED"),
           showWarnings = FALSE)
dir.create(file.path(out_dir, "QC"),
           showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"),
           showWarnings = FALSE)

## Parameters (EDIT if needed)
assay_name <- "peaks_universal"
# tissue_col removed: kidney_merged_universal.rds is already a per-tissue object,
# no `Tissue == "Kidney"` subset needed.
group_col <- "condition"
celltype_col <- "cell_type"

group_levels <- c("Sham", "Day14", "Day42")

# Ralph suggested testing exp.thresh (5% default)
exp_thresh <- 0.05
num_splits <- 10
padj_cutoff <- 0.05

# Optional stricter minimum cell number per group
min_cells_per_group <- 80

# Cell types to run (Kidney validated set)
eligible_celltypes <- c(
  "DCT", "Endothelial", "IC", "Macrophages",
  "PC", "PT", "Podocytes", "TAL"
)

# Contrasts
contrasts_list <- list(
  c("Day14", "Sham"),
  c("Day42", "Sham"),
  c("Day42", "Day14")
)

## Load helper functions
if (!file.exists(datesting_r_path)) {
  stop("DATesting.R not found: ", datesting_r_path)
}
source(datesting_r_path)

# Check functions exist

## Load Seurat object
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

meta_cols <- colnames(obj[[]])
req_cols <- c(group_col, celltype_col)
missing_cols <- setdiff(req_cols, meta_cols)
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ",
       paste(missing_cols, collapse = ", "))
}

DefaultAssay(obj) <- assay_name

## Metadata prep (no tissue subset - obj is already kidney-only)
obj_kidney <- obj

obj_kidney$dar_group <- as.character(obj_kidney[[group_col]][, 1])
obj_kidney$cell_type_dar <- as.character(obj_kidney[[celltype_col]][, 1])

obj_kidney <- subset(
  obj_kidney,
  subset = dar_group %in% group_levels
)

obj_kidney$dar_group <- factor(
  obj_kidney$dar_group,
  levels = group_levels
)

# Keep only selected cell types
eligible_celltypes <- intersect(
  eligible_celltypes,
  unique(obj_kidney$cell_type_dar)
)

write.csv(
  as.data.frame(table(obj_kidney$dar_group)),
  file.path(out_dir, "QC", "Kidney_dar_group_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(obj_kidney$cell_type_dar)),
  file.path(out_dir, "QC", "Kidney_cell_type_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(obj_kidney$cell_type_dar, obj_kidney$dar_group)),
  file.path(out_dir, "QC", "Kidney_celltype_by_group_counts.csv"),
  row.names = FALSE
)

message("Kidney cells: ", ncol(obj_kidney))
message("Eligible cell types: ",
        paste(eligible_celltypes, collapse = ", "))

## Helpers
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

safe_hist <- function(x, title_txt, out_png) {
  png(out_png, width = 1600, height = 1200, res = 200)
  hist(x, breaks = 50, main = title_txt, xlab = "p-value")
  dev.off()
}

## Run pseudo-bulk DAR loop
all_summary <- list()

for (ct in eligible_celltypes) {
  message("Processing cell type: ", ct)
  ct_obj <- subset(obj_kidney, subset = cell_type_dar == ct)
  DefaultAssay(ct_obj) <- assay_name
  
  # Use dar_group as identity for selecting cells
  Idents(ct_obj) <- "dar_group"
  
  # QC counts for this cell type
  ct_counts <- table(ct_obj$dar_group, useNA = "ifany")
  print(ct_counts)
  
  write.csv(
    as.data.frame(ct_counts),
    file.path(
      out_dir,
      "QC",
      paste0(gsub("[ /]", "_", ct), "__group_counts.csv")
    ),
    row.names = FALSE
  )
  
  for (cc in contrasts_list) {
    g1 <- cc[1]
    g2 <- cc[2]
    contrast_name <- paste0(g1, "_vs_", g2)
    
    message("  -> Contrast: ", contrast_name)
    
    cells_g1 <- colnames(ct_obj)[Idents(ct_obj) == g1]
    cells_g2 <- colnames(ct_obj)[Idents(ct_obj) == g2]
    
    n1 <- length(cells_g1)
    n2 <- length(cells_g2)
    
    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      message("     Skip: too few cells (",
              g1, "=", n1, ", ", g2, "=", n2, ")")
      next
    }
    
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
      error = function(e) {
        message("     ERROR in pseudo-bulk DAR: ",
                conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(res)) {
      next
    }
    
    res_df <- as.data.frame(res)
    res_df$peak <- rownames(res_df)
    
    # Reorder columns
    first_cols <- c("peak", "baseMean", "log2FoldChange",
                    "lfcSE", "stat", "pvalue", "padj")
    keep_first <- intersect(first_cols, colnames(res_df))
    res_df <- res_df[, c(keep_first,
                         setdiff(colnames(res_df), keep_first)),
                     drop = FALSE]
    
    # Significant
    res_sig <- subset(res_df, !is.na(padj) & padj < padj_cutoff)
    
    # Opening = positive log2FC (g1 > g2)
    res_open <- subset(res_sig, log2FoldChange > 0)
    
    # Closing = negative log2FC (g1 < g2)
    res_close <- subset(res_sig, log2FoldChange < 0)
    
    # Prefix label
    out_label <- paste0(
      gsub("[ /]", "_", ct),
      "__", contrast_name,
      "__",
      gsub("\\.", "", sub("^0\\.", "0", format(exp_thresh)))
    )
    
    write.table(
      res_df,
      file = file.path(out_dir, "DAR_tables",
                       paste0(out_label, "_DESeq2_all.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.table(
      res_sig,
      file = file.path(out_dir, "DAR_tables",
                       paste0(out_label, "_DESeq2_padj005.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.table(
      res_open,
      file = file.path(out_dir, "DAR_tables",
                       paste0(out_label, "_opening.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.table(
      res_close,
      file = file.path(out_dir, "DAR_tables",
                       paste0(out_label, "_closing.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # Save BED for HOMER
    peak_to_bed(
      peaks = res_open$peak,
      out_file = file.path(
        out_dir, "DAR_BED",
        paste0(gsub("[ /]", "_", ct),
               "__", contrast_name, "__opening.bed")
      )
    )
    
    peak_to_bed(
      peaks = res_close$peak,
      out_file = file.path(
        out_dir, "DAR_BED",
        paste0(gsub("[ /]", "_", ct),
               "__", contrast_name, "__closing.bed")
      )
    )
    
    # p-value histogram
    if ("pvalue" %in% colnames(res_df)) {
      safe_hist(
        x = res_df$pvalue,
        title_txt = paste(ct, contrast_name, "pseudo-bulk DAR p-values"),
        out_png = file.path(
          out_dir, "Figures",
          paste0(gsub("[ /]", "_", ct),
                 "__", contrast_name, "__pvalue_hist.png")
        )
      )
    }
    
    # Summary row
    all_summary[[paste(ct, contrast_name, sep = "__")]] <- data.frame(
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
    
    message("     Done. tested=", nrow(res_df),
            ", sig=", nrow(res_sig),
            ", opening=", nrow(res_open),
            ", closing=", nrow(res_close))
  }
  
  rm(ct_obj)
  gc()
}

## Save summary + Count plot
if (length(all_summary) == 0) {
  stop("No pseudo-bulk DAR results were generated.")
}

summary_df <- bind_rows(all_summary)

write.csv(
  summary_df,
  file.path(out_dir, "DAR_pseudobulk_summary_Kidney.csv"),
  row.names = FALSE
)

# Long format for plotting opening/closing counts
plot_df <- bind_rows(
  summary_df %>%
    transmute(cell_type, contrast, direction = "Opening",
              n_DAR = n_opening),
  summary_df %>%
    transmute(cell_type, contrast, direction = "Closing",
              n_DAR = n_closing)
)

write.csv(
  plot_df,
  file.path(out_dir, "DAR_pseudobulk_summary_Kidney_long.csv"),
  row.names = FALSE
)

# Bar plot: opening/closing DAR counts by cell type + contrast
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
    title = "Kidney pseudo-bulk DAR counts (DESeq2)",
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
                       "Kidney_pseudobulk_DAR_counts_barplot.pdf"),
  plot = p_counts,
  width = 14,
  height = 6
)

ggsave(
  filename = file.path(out_dir, "Figures",
                       "Kidney_pseudobulk_DAR_counts_barplot.png"),
  plot = p_counts,
  width = 14,
  height = 6,
  dpi = 300
)

# Optional heatmap-like tile of counts
p_tile <- ggplot(
  plot_df,
  aes(x = contrast, y = cell_type, fill = n_DAR)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_DAR), size = 3) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#CC0000") +
  labs(
    title = "Kidney pseudo-bulk DAR counts by contrast and cell type",
    x = "Contrast",
    y = "Cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "Figures",
                       "Kidney_pseudobulk_DAR_counts_tile.pdf"),
  plot = p_tile,
  width = 12,
  height = 6.5
)

ggsave(
  filename = file.path(out_dir, "Figures",
                       "Kidney_pseudobulk_DAR_counts_tile.png"),
  plot = p_tile,
  width = 12,
  height = 6.5,
  dpi = 300
)

message("Pseudo-bulk DAR wrapper finished.")
message("Results saved to: ", out_dir)
