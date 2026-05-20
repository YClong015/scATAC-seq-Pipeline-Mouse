#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# Make opening / closing / stable / NS BED files
# for selected tissues and selected contrasts only
# ============================================================

# -----------------------------
# 1) User config
# -----------------------------
cfg <- list(
  list(
    tissue = "Kidney",
    base_dir = paste0(
      "/QRISdata/Q8448/Mouse_disease_data/DAR/",
      "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2"
    ),
    contrast_keep = "Day42_vs_Sham"
  ),
  list(
    tissue = "Aorta",
    base_dir = paste0(
      "/QRISdata/Q8448/Mouse_disease_data/DAR/",
      "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2"
    ),
    contrast_keep = "Challenge_vs_Control"
  ),
  list(
    tissue = "Lung",
    base_dir = paste0(
      "/QRISdata/Q8448/Mouse_disease_data/DAR/",
      "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2"
    ),
    contrast_keep = "Case_vs_Control"
  ),
  list(
    tissue = "Tcell",
    base_dir = paste0(
      "/QRISdata/Q8448/Mouse_disease_data/DAR/",
      "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2"
    ),
    contrast_keep = "Young_chronic_vs_Young_control"
  )
)

stable_padj_cut <- 0.9
stable_lfc_cut <- 0.05
sig_padj_cut <- 0.05

# -----------------------------
# 2) Helpers
# -----------------------------
peak_to_bed <- function(df, outfile) {
  if (nrow(df) == 0) {
    file.create(outfile)
    return(invisible(NULL))
  }
  
  bed <- data.frame(
    chr = sub("(.*)-.*-.*", "\\1", df$peak),
    start = sub(".*-(.*)-.*", "\\1", df$peak),
    end = sub(".*-.*-(.*)", "\\1", df$peak),
    stringsAsFactors = FALSE
  )
  
  write.table(
    bed,
    file = outfile,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# -----------------------------
# 3) Main
# -----------------------------
for (x in cfg) {
  tissue <- x$tissue
  base_dir <- x$base_dir
  contrast_keep <- x$contrast_keep
  
  in_dir <- file.path(base_dir, "DAR_tables")
  out_stable <- file.path(base_dir, "DAR_BED_stable")
  out_ns <- file.path(base_dir, "DAR_BED_NS")
  
  dir.create(out_stable, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_ns, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(
    in_dir,
    pattern = "_DESeq2_all\\.tsv$",
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message("No DESeq2 all files found for ", tissue)
    next
  }
  
  summary_list <- list()
  
  for (f in files) {
    bn <- basename(f)
    
    if (!grepl(contrast_keep, bn)) {
      next
    }
    
    res <- read.delim(f, stringsAsFactors = FALSE)
    
    need_cols <- c("peak", "log2FoldChange", "padj")
    if (!all(need_cols %in% colnames(res))) {
      message("Skipping missing columns: ", bn)
      next
    }
    
    stub <- sub("_DESeq2_all\\.tsv$", "", bn)
    
    opening <- res %>%
      filter(!is.na(padj), padj < sig_padj_cut, log2FoldChange > 0)
    
    closing <- res %>%
      filter(!is.na(padj), padj < sig_padj_cut, log2FoldChange < 0)
    
    stable <- res %>%
      filter(
        !is.na(padj),
        padj > stable_padj_cut,
        abs(log2FoldChange) < stable_lfc_cut
      )
    
    ns <- res %>%
      filter(!is.na(padj), padj > sig_padj_cut)
    
    # stable BED set
    peak_to_bed(opening, file.path(out_stable, paste0(stub, "__opening.bed")))
    peak_to_bed(closing, file.path(out_stable, paste0(stub, "__closing.bed")))
    peak_to_bed(stable, file.path(out_stable, paste0(stub, "__stable.bed")))
    
    # NS BED set
    peak_to_bed(opening, file.path(out_ns, paste0(stub, "__opening.bed")))
    peak_to_bed(closing, file.path(out_ns, paste0(stub, "__closing.bed")))
    peak_to_bed(ns, file.path(out_ns, paste0(stub, "__NS.bed")))
    
    summary_list[[stub]] <- data.frame(
      tissue = tissue,
      comparison = stub,
      n_opening = nrow(opening),
      n_closing = nrow(closing),
      n_stable = nrow(stable),
      n_NS = nrow(ns),
      stringsAsFactors = FALSE
    )
    
    message(
      tissue, " | ", stub,
      " | opening=", nrow(opening),
      " closing=", nrow(closing),
      " stable=", nrow(stable),
      " NS=", nrow(ns)
    )
  }
  
  if (length(summary_list) > 0) {
    summary_df <- bind_rows(summary_list)
    write.csv(
      summary_df,
      file.path(base_dir, paste0(tissue, "_selected_counts.csv")),
      row.names = FALSE
    )
  }
}

message("All done.")

