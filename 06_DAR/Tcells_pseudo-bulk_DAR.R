#!/usr/bin/env Rscript

# Tcells pseudo-bulk DAR pipeline (DESeq2 via DATesting.R)
# Lung-style structure:
# - explicit group mapping from metadata column
# - optional fragment mounting (recommended for downstream tracks)
# - FORCED global cell type assignment to "Tcell"

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

## User paths (EDIT if needed)
obj_rds_path <- "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
frag_path <- "/QRISdata/Q8448/Mouse_disease_data/Tcells/atac_fragments.tsv.gz"
datesting_r_path <- "/home/s4869245/scripts/DAR/DATesting.R"

out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_Tcells_DESeq2"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DAR_tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "DAR_BED"), showWarnings = FALSE)
dir.create(file.path(out_dir, "QC"), showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), showWarnings = FALSE)

## Parameters
# Assay name in Tcells object (auto-detect if peaks_universal missing)
assay_name_prefer <- "peaks_universal"

# Tcells group label column
group_col <- "deMultliplex2_final_mapped"

# Group levels
group_levels <- c(
  "Young control",
  "Young acute",
  "Young chronic",
  "Juvenile",
  "Aged"
)

# Contrasts to run
contrasts_list <- list(
  c("Aged", "Juvenile"),
  c("Aged", "Young acute"),
  c("Aged", "Young chronic"),
  c("Aged", "Young control"),
  c("Juvenile", "Young acute"),
  c("Juvenile", "Young chronic"),
  c("Juvenile", "Young control"),
  c("Young acute", "Young control"),
  c("Young chronic", "Young control")
)

# DESeq2 pseudo-bulk settings
exp_thresh <- 0.05
num_splits <- 10
padj_cutoff <- 0.05
min_cells_per_group <- 80

# Lung-style: optionally mount fragments onto assay (for track/footprint later)
mount_fragments <- TRUE

## Helpers
stop_with <- function(...) stop(paste0(...), call. = FALSE)
sanitize_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

peak_to_bed <- function(peaks, out_file) {
  if (length(peaks) == 0) return(invisible(NULL))
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
  if (length(x) == 0) return(invisible(NULL))
  png(out_png, width = 1600, height = 1200, res = 200)
  hist(x, breaks = 50, main = title_txt, xlab = "p-value")
  dev.off()
}

# Best-effort mapping from Seurat cell names -> raw barcodes in fragments file
make_raw_barcodes <- function(cell_names) {
  x <- sub("^Tcells_", "", cell_names)
  if (any(grepl("_", x))) {
    x2 <- sub("^.*_", "", x)
    use_x2 <- grepl("-[0-9]+$", x2)
    x[use_x2] <- x2[use_x2]
  }
  x
}

## Load helper methods + object
if (!file.exists(datesting_r_path)) {
  stop_with("DATesting.R not found: ", datesting_r_path)
}
source(datesting_r_path)
stopifnot(exists("apply_DESeq2_test_seurat"))
stopifnot(exists("GetExpressedPeaks"))

if (!file.exists(obj_rds_path)) stop_with("Object not found: ", obj_rds_path)
message("Loading Tcells object...")
obj <- readRDS(obj_rds_path)
message("Loaded. Cells = ", ncol(obj))

# Assay auto-detect
assays <- Assays(obj)
assay_name <- if (assay_name_prefer %in% assays) {
  assay_name_prefer
} else if ("ATAC" %in% assays) {
  "ATAC"
} else {
  stop_with(
    "Cannot find assay peaks_universal or ATAC. Available: ",
    paste(assays, collapse = ", ")
  )
}
DefaultAssay(obj) <- assay_name
message("Using assay: ", assay_name)

if (!group_col %in% colnames(obj[[]])) {
  stop_with("Missing group column: ", group_col)
}

## Group mapping & Global Cell Type Assignment
obj$dar_group <- as.character(obj[[group_col]][, 1])

# Assign all cells as "Tcell" uniformly
obj$cell_type_dar <- "Tcell"

# QC: show what is inside
write.csv(
  as.data.frame(table(obj$dar_group, useNA = "ifany")),
  file.path(out_dir, "QC", "Tcells_group_counts_raw.csv"),
  row.names = FALSE
)

# Keep only the desired groups
obj <- subset(obj, subset = dar_group %in% group_levels)
obj$dar_group <- factor(obj$dar_group, levels = group_levels)

write.csv(
  as.data.frame(table(obj$dar_group, useNA = "ifany")),
  file.path(out_dir, "QC", "Tcells_group_counts_filtered.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(obj$cell_type_dar, useNA = "ifany")),
  file.path(out_dir, "QC", "Tcells_celltype_counts.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(obj$cell_type_dar, obj$dar_group, useNA = "ifany")),
  file.path(out_dir, "QC", "Tcells_celltype_by_group_counts.csv"),
  row.names = FALSE
)

message("After filtering: cells = ", ncol(obj))

## Optional: mount fragments
if (mount_fragments) {
  if (!file.exists(frag_path)) stop_with("Fragment file missing: ", frag_path)
  
  message("Mounting fragments to assay: ", assay_name)
  cells_seurat <- colnames(obj)
  raw_barcodes <- make_raw_barcodes(cells_seurat)
  names(raw_barcodes) <- cells_seurat
  
  frag_obj <- CreateFragmentObject(
    path = frag_path,
    cells = raw_barcodes,
    validate.fragments = FALSE
  )
  
  # Attach to the active ChromatinAssay
  tryCatch(
    {
      Fragments(obj[[assay_name]]) <- list(frag_obj)
    },
    error = function(e) {
      message("WARNING: Failed to attach fragments: ", conditionMessage(e))
    }
  )
  
  saveRDS(obj, file.path(out_dir, "tcells_with_fragments.rds"))
  message("Saved: ", file.path(out_dir, "tcells_with_fragments.rds"))
}

## Pseudo-bulk DAR loop
DefaultAssay(obj) <- assay_name
eligible_celltypes <- sort(unique(obj$cell_type_dar))
message("Cell types: ", paste(eligible_celltypes, collapse = ", "))

all_summary <- list()

for (ct in eligible_celltypes) {
  message("Processing cell type: ", ct)
  ct_obj <- subset(obj, subset = cell_type_dar == ct)
  DefaultAssay(ct_obj) <- assay_name
  Idents(ct_obj) <- "dar_group"
  
  ct_counts <- table(ct_obj$dar_group, useNA = "ifany")
  write.csv(
    as.data.frame(ct_counts),
    file.path(out_dir, "QC", paste0(sanitize_name(ct), "__group_counts.csv")),
    row.names = FALSE
  )
  
  for (cc in contrasts_list) {
    g1 <- cc[1]
    g2 <- cc[2]
    contrast_name <- paste0(sanitize_name(g1), "_vs_", sanitize_name(g2))
    message("  -> Contrast: ", g1, " vs ", g2)
    
    cells_g1 <- colnames(ct_obj)[Idents(ct_obj) == g1]
    cells_g2 <- colnames(ct_obj)[Idents(ct_obj) == g2]
    n1 <- length(cells_g1)
    n2 <- length(cells_g2)
    
    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      message("     Skip: too few cells (", g1, "=", n1, ", ", g2, "=", n2, ")")
      next
    }
    
    res <- tryCatch(
      apply_DESeq2_test_seurat(
        seurat.object = ct_obj,
        population.1 = cells_g1,
        population.2 = cells_g2,
        exp.thresh = exp_thresh,
        num.splits = num_splits,
        assay.use = assay_name,
        verbose = TRUE
      ),
      error = function(e) {
        message("     ERROR: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res)) next
    
    res_df <- as.data.frame(res)
    res_df$peak <- rownames(res_df)
    
    res_sig <- subset(res_df, !is.na(padj) & padj < padj_cutoff)
    res_open <- subset(res_sig, log2FoldChange > 0)
    res_close <- subset(res_sig, log2FoldChange < 0)
    
    out_label <- paste0(
      sanitize_name(ct), "__", contrast_name, "__",
      gsub("\\.", "", sub("^0\\.", "0", format(exp_thresh)))
    )
    
    write.table(
      res_df,
      file.path(out_dir, "DAR_tables", paste0(out_label, "_DESeq2_all.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      res_sig,
      file.path(out_dir, "DAR_tables", paste0(out_label, "_DESeq2_padj005.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    peak_to_bed(
      res_open$peak,
      file.path(out_dir, "DAR_BED",
                paste0(sanitize_name(ct), "__", contrast_name, "__opening.bed"))
    )
    peak_to_bed(
      res_close$peak,
      file.path(out_dir, "DAR_BED",
                paste0(sanitize_name(ct), "__", contrast_name, "__closing.bed"))
    )
    
    if ("pvalue" %in% colnames(res_df)) {
      safe_hist(
        x = res_df$pvalue,
        title_txt = paste(ct, g1, "vs", g2, "pseudo-bulk p-values"),
        out_png = file.path(out_dir, "Figures",
                            paste0(out_label, "__pvalue_hist.png"))
      )
    }
    
    all_summary[[paste(ct, contrast_name, sep = "__")]] <- data.frame(
      tissue = "Tcells",
      cell_type = ct,
      contrast = paste0(g1, "_vs_", g2),
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

if (length(all_summary) == 0) stop_with("No DAR results generated.")

summary_df <- bind_rows(all_summary)
write.csv(summary_df,
          file.path(out_dir, "DAR_pseudobulk_summary_Tcells.csv"),
          row.names = FALSE)

message("Done. Output folder: ", out_dir)

# 8) Save summaries + plots
if (length(all_summary) == 0) stop_with("No pseudo-bulk DAR results were generated.")

summary_df <- bind_rows(all_summary)
write.csv(summary_df, file.path(out_dir, "DAR_pseudobulk_summary_Tcells.csv"), row.names = FALSE)

plot_df <- bind_rows(
  summary_df %>% transmute(cell_type, contrast, direction = "Opening", n_DAR = n_opening),
  summary_df %>% transmute(cell_type, contrast, direction = "Closing", n_DAR = n_closing)
)
write.csv(plot_df, file.path(out_dir, "DAR_pseudobulk_summary_Tcells_long.csv"), row.names = FALSE)

## Bar plot
p_counts <- ggplot(plot_df, aes(x = cell_type, y = n_DAR, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ contrast, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = c("Opening" = "#D73027", "Closing" = "#4575B4")) +
  labs(
    title = "Tcells pseudo-bulk DAR counts (DESeq2)",
    x = "Cell type",
    y = "Number of DARs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Figures", "Tcells_pseudobulk_DAR_counts_barplot.pdf"), p_counts, width = 14, height = 6)
ggsave(file.path(out_dir, "Figures", "Tcells_pseudobulk_DAR_counts_barplot.png"), p_counts, width = 14, height = 6, dpi = 300)

## Tile plot (Heatmap)
plot_df$n_DAR_signed <- ifelse(plot_df$direction == "Closing", -plot_df$n_DAR, plot_df$n_DAR)

p_tile <- ggplot(plot_df, aes(x = contrast, y = cell_type, fill = n_DAR_signed)) +
  geom_tile(color = "white") +
  # Note: The text labels still use the absolute positive numbers (n_DAR)
  geom_text(aes(label = n_DAR), size = 3) + 
  facet_wrap(~ direction, nrow = 1) +
  # Use gradient2 for the red-white-blue diverging color scale
  scale_fill_gradient2(
    low = "#4575B4",    # Blue represents Closing (mapped to negative values)
    mid = "white",      # White represents 0
    high = "#D73027",   # Red represents Opening (mapped to positive values)
    midpoint = 0,
    labels = abs,       # Crucial: Ensures the legend displays absolute (positive) numbers!
    name = "DAR Count"
  ) +
  labs(
    title = "Tcells pseudo-bulk DAR counts by contrast and cell type",
    x = "Contrast",
    y = "Cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Figures", "Tcells_pseudobulk_DAR_counts_tile.pdf"), p_tile, width = 12, height = 6.5)
ggsave(file.path(out_dir, "Figures", "Tcells_pseudobulk_DAR_counts_tile.png"), p_tile, width = 12, height = 6.5, dpi = 300)

message("Done.")
message("Output folder: ", out_dir)
message("Design: Global Tcell pseudo-bulk DAR complete across 5 contrasts.")
