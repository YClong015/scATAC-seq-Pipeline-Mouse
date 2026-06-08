#!/usr/bin/env Rscript

# Lung pseudo-bulk DAR pipeline (DESeq2 via DATesting.R)
# Same structure as Tcells_pseudo-bulk_DAR.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
})

## Paths
obj_rds_path     <- "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds"
datesting_r_path <- "/home/s4869245/scripts/DAR/DATesting.R"

out_dir <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Lung_DESeq2"

dir.create(file.path(out_dir, "DAR_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DAR_BED"),   showWarnings = FALSE)
dir.create(file.path(out_dir, "QC"),        showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"),   showWarnings = FALSE)

## Parameters
assay_name          <- "peaks_universal"
group_col           <- "Group"            # metadata column for Case/Control
celltype_col        <- "cell_type"        # metadata column for cell type

group_levels        <- c("Control", "Case")
contrasts_list      <- list(c("Case", "Control"))

exp_thresh          <- 0.05
num_splits          <- 10
padj_cutoff         <- 0.05
min_cells_per_group <- 80

## Helpers
sanitize_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

peak_to_bed <- function(peaks, out_file) {
  if (length(peaks) == 0) return(invisible(NULL))
  df <- data.frame(
    Chr   = sub("(.*)-.*-.*",  "\\1", peaks),
    Start = sub(".*-(.*)-.*",  "\\1", peaks),
    End   = sub(".*-.*-(.*)",  "\\1", peaks),
    stringsAsFactors = FALSE
  )
  write.table(df, file = out_file, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

safe_hist <- function(x, title_txt, out_png) {
  if (length(x) == 0) return(invisible(NULL))
  png(out_png, width = 1600, height = 1200, res = 200)
  hist(x, breaks = 50, main = title_txt, xlab = "p-value")
  dev.off()
}

## Load DATesting.R + object
if (!file.exists(datesting_r_path)) stop("DATesting.R not found: ", datesting_r_path)
source(datesting_r_path)

if (!file.exists(obj_rds_path)) stop("Object not found: ", obj_rds_path)
message("Loading Lung object...")
obj <- readRDS(obj_rds_path)
message("Loaded. Cells = ", ncol(obj))

if (!assay_name %in% Assays(obj))
  stop("Assay not found: ", assay_name,
       ". Available: ", paste(Assays(obj), collapse = ", "))
if (!group_col %in% colnames(obj[[]]))
  stop("Missing group column: ", group_col,
       ". Available: ", paste(colnames(obj[[]]), collapse = ", "))
if (!celltype_col %in% colnames(obj[[]]))
  stop("Missing celltype column: ", celltype_col)

DefaultAssay(obj) <- assay_name

## Group + cell type mapping
obj$dar_group     <- as.character(obj[[group_col]][, 1])
obj$cell_type_dar <- as.character(obj[[celltype_col]][, 1])

# QC before filter
write.csv(as.data.frame(table(obj$dar_group, useNA = "ifany")),
          file.path(out_dir, "QC", "Lung_group_counts_raw.csv"), row.names = FALSE)

obj <- subset(obj, subset = dar_group %in% group_levels)
obj$dar_group <- factor(obj$dar_group, levels = group_levels)

message("Cells after group filter: ", ncol(obj))

ct_grp <- as.data.frame(table(CellType = obj$cell_type_dar, Group = obj$dar_group))
message("\n=== Cells per cell type x group ===")
print(ct_grp[order(ct_grp$CellType), ], row.names = FALSE)

write.csv(ct_grp,
          file.path(out_dir, "QC", "Lung_celltype_by_group_counts.csv"), row.names = FALSE)
write.csv(as.data.frame(table(obj$dar_group)),
          file.path(out_dir, "QC", "Lung_group_counts_filtered.csv"), row.names = FALSE)

# Eligible cell types: min cells >= threshold in BOTH groups
eligible_celltypes <- ct_grp %>%
  group_by(CellType) %>%
  summarise(min_cells = min(Freq), .groups = "drop") %>%
  filter(min_cells >= min_cells_per_group) %>%
  pull(CellType) %>% sort()

skipped <- setdiff(unique(obj$cell_type_dar), eligible_celltypes)

message("\nEligible (>= ", min_cells_per_group, " cells/group): ",
        paste(eligible_celltypes, collapse = ", "))
if (length(skipped) > 0)
  message("Skipped (too few cells): ", paste(skipped, collapse = ", "))

## Pseudo-bulk DAR loop
all_summary <- list()
Idents(obj)  <- "dar_group"

for (ct in eligible_celltypes) {
  message("Processing: ", ct)
  ct_obj <- subset(obj, subset = cell_type_dar == ct)
  DefaultAssay(ct_obj) <- assay_name
  Idents(ct_obj) <- "dar_group"

  write.csv(as.data.frame(table(ct_obj$dar_group, useNA = "ifany")),
            file.path(out_dir, "QC", paste0(sanitize_name(ct), "__group_counts.csv")),
            row.names = FALSE)

  for (cc in contrasts_list) {
    g1 <- cc[1]; g2 <- cc[2]
    contrast_name <- paste0(g1, "_vs_", g2)
    message("  -> Contrast: ", contrast_name)

    cells_g1 <- colnames(ct_obj)[Idents(ct_obj) == g1]
    cells_g2 <- colnames(ct_obj)[Idents(ct_obj) == g2]
    n1 <- length(cells_g1); n2 <- length(cells_g2)

    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      message("     Skip: too few cells (", g1, "=", n1, ", ", g2, "=", n2, ")")
      next
    }

    res <- tryCatch(
      apply_DESeq2_test_seurat(
        seurat.object = ct_obj,
        population.1  = cells_g1,
        population.2  = cells_g2,
        exp.thresh    = exp_thresh,
        num.splits    = num_splits,
        assay.use     = assay_name,
        verbose       = TRUE
      ),
      error = function(e) { message("     ERROR: ", e$message); NULL }
    )
    if (is.null(res)) next

    res_df <- as.data.frame(res)
    res_df$peak <- rownames(res_df)
    first_cols  <- c("peak", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
    res_df      <- res_df[, c(intersect(first_cols, colnames(res_df)),
                               setdiff(colnames(res_df), first_cols)), drop = FALSE]

    res_sig   <- subset(res_df, !is.na(padj) & padj < padj_cutoff)
    res_open  <- subset(res_sig, log2FoldChange > 0)
    res_close <- subset(res_sig, log2FoldChange < 0)

    out_label <- paste0(sanitize_name(ct), "__", contrast_name, "__",
                        gsub("\\.", "", sub("^0\\.", "0", format(exp_thresh))))

    write.table(res_df,
                file.path(out_dir, "DAR_tables", paste0(out_label, "_DESeq2_all.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(res_sig,
                file.path(out_dir, "DAR_tables", paste0(out_label, "_DESeq2_padj005.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(res_open,
                file.path(out_dir, "DAR_tables", paste0(out_label, "_opening.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(res_close,
                file.path(out_dir, "DAR_tables", paste0(out_label, "_closing.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    peak_to_bed(res_open$peak,
                file.path(out_dir, "DAR_BED",
                          paste0(sanitize_name(ct), "__", contrast_name, "__opening.bed")))
    peak_to_bed(res_close$peak,
                file.path(out_dir, "DAR_BED",
                          paste0(sanitize_name(ct), "__", contrast_name, "__closing.bed")))

    if ("pvalue" %in% colnames(res_df))
      safe_hist(res_df$pvalue,
                paste(ct, contrast_name, "p-values"),
                file.path(out_dir, "Figures",
                          paste0(out_label, "__pvalue_hist.png")))

    all_summary[[paste(ct, contrast_name, sep = "__")]] <- data.frame(
      tissue         = "Lung",
      cell_type      = ct,
      contrast       = contrast_name,
      group1         = g1,
      group2         = g2,
      n_cells_group1 = n1,
      n_cells_group2 = n2,
      exp_thresh     = exp_thresh,
      num_splits     = num_splits,
      n_tested_peaks = nrow(res_df),
      n_sig_total    = nrow(res_sig),
      n_opening      = nrow(res_open),
      n_closing      = nrow(res_close),
      stringsAsFactors = FALSE
    )

    message("     Done. tested=", nrow(res_df),
            " | opening=", nrow(res_open),
            " | closing=", nrow(res_close))
  }

  rm(ct_obj); gc()
}

## Summary + plots
if (length(all_summary) == 0) stop("No DAR results generated.")

summary_df <- bind_rows(all_summary)
write.csv(summary_df,
          file.path(out_dir, "DAR_pseudobulk_summary_Lung.csv"), row.names = FALSE)

plot_df <- bind_rows(
  summary_df %>% transmute(cell_type, contrast, direction = "Opening", n_DAR = n_opening),
  summary_df %>% transmute(cell_type, contrast, direction = "Closing", n_DAR = n_closing)
)
write.csv(plot_df,
          file.path(out_dir, "DAR_pseudobulk_summary_Lung_long.csv"), row.names = FALSE)

p_counts <- ggplot(plot_df, aes(x = cell_type, y = n_DAR, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ contrast, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(Opening = "#D73027", Closing = "#4575B4")) +
  labs(title = "Lung pseudo-bulk DAR counts (DESeq2, Case vs Control)",
       x = "Cell type", y = "Number of DARs") +
  theme_bw(base_size = 12) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figures", "Lung_DAR_counts_barplot.pdf"),
       p_counts, width = 14, height = 6)
ggsave(file.path(out_dir, "Figures", "Lung_DAR_counts_barplot.png"),
       p_counts, width = 14, height = 6, dpi = 300)

plot_df$n_DAR_signed <- ifelse(plot_df$direction == "Closing", -plot_df$n_DAR, plot_df$n_DAR)

p_tile <- ggplot(plot_df, aes(x = contrast, y = cell_type, fill = n_DAR_signed)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_DAR), size = 3) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, labels = abs, name = "DAR Count") +
  labs(title = "Lung pseudo-bulk DAR counts by cell type",
       x = "Contrast", y = "Cell type") +
  theme_bw(base_size = 12) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figures", "Lung_DAR_counts_tile.pdf"),
       p_tile, width = 12, height = 6.5)
ggsave(file.path(out_dir, "Figures", "Lung_DAR_counts_tile.png"),
       p_tile, width = 12, height = 6.5, dpi = 300)

message("\nDone. Output: ", out_dir)
