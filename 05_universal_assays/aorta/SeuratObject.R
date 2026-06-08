#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
})

msg <- function(x) cat("[INFO]", x, "\n")

get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste0("Missing value after ", key))
  args[[hit + 1]]
}

# Barcode cleaner
# specifically built to handle Aorta suffixes (_1) and Lung prefixes (Control_F2_...)
clean_barcodes <- function(x) {
  # 1. Remove Seurat Merge suffixes (e.g., AAACGAA...-1_1 -> AAACGAA...-1)
  x <- sub("_[0-9]+$", "", x)
  # 2. Remove complex Lung prefixes (e.g., Control_F2_CELL771_N1 -> CELL771_N1)
  x <- sub("^.*?(CELL.*)", "\\1", x)
  return(x)
}

load_peaks_bed_gz <- function(path) {
  con <- gzfile(path, "rt")
  df <- read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  close(con)
  
  if (ncol(df) < 3) stop("Universal bed has <3 columns.")
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  
  bad <- is.na(df$start) | is.na(df$end) | df$end <= df$start
  df <- df[!bad, , drop = FALSE]
  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1L, end = df$end)
  )
  gr
}

args <- commandArgs(trailingOnly = TRUE)

obj_path <- get_arg(args, "--obj", "")
obj_type <- get_arg(args, "--objtype", "rds")
up_bed <- get_arg(args, "--up", "")
out_rds <- get_arg(args, "--out", "")
assay_in <- get_arg(args, "--assay_in", "peaks")
assay_out <- get_arg(args, "--assay_out", "peaks_universal")
mode <- get_arg(args, "--mode", "per_sample")
sample_key <- get_arg(args, "--sample_key", "dataset")
frag_file <- get_arg(args, "--frag_file", "")
frag_tpl <- get_arg(args, "--frag_tpl", "")

if (obj_path == "" || up_bed == "" || out_rds == "") {
  stop("Need --obj --up --out")
}

msg("Step 1: Load object")
if (obj_type == "rds") {
  obj <- readRDS(obj_path)
} else {
  stop("objtype must be rds")
}

if (!assay_in %in% names(obj@assays)) {
  stop(paste0("Assay not found: ", assay_in))
}
DefaultAssay(obj) <- assay_in
msg(paste0("Successfully loaded: ", obj_path))
msg(paste0("Total Cells: ", ncol(obj)))

msg("Step 2: Load universal peaks (bed.gz)")
peaks <- load_peaks_bed_gz(up_bed)
msg(paste0("Universal Peaks count: ", length(peaks)))

msg("Step 3: Build FeatureMatrix on universal peaks")
counts_list <- list()

if (mode == "per_sample") {
  if (frag_tpl == "") stop("per_sample mode needs --frag_tpl")
  if (!sample_key %in% colnames(obj@meta.data)) {
    stop(paste0("Missing sample_key in metadata: ", sample_key))
  }

  sample_ids <- sort(unique(obj[[sample_key]][, 1]))
  msg(paste0("Found ", length(sample_ids), " samples to process."))
  
  for (id in sample_ids) {
    frag_path <- gsub("\\{id\\}", id, frag_tpl)
    if (!file.exists(frag_path)) stop(paste0("Missing fragments file: ", frag_path))

    # Get cells for the current sample
    cells_pref <- colnames(obj)[obj[[sample_key]][, 1] == id]
    if (length(cells_pref) == 0) next
    
    # Clean the cell barcodes to strictly match the Fragment file
    cells_raw <- clean_barcodes(cells_pref)
    msg(paste0("  -> Processing sample: ", id, " (", length(cells_pref), " cells)"))

    frag <- CreateFragmentObject(path = frag_path)
    mat <- FeatureMatrix(fragments = frag, features = peaks, cells = cells_raw)

    if (ncol(mat) == 0) {
      msg(paste0("    [WARNING] No cells matched in fragment file for ", id))
      next
    }

    # Restore the original Seurat names so merging works later
    colnames(mat) <- cells_pref
    counts_list[[id]] <- mat
  } # The For loop is now correctly closed here!

  if (length(counts_list) == 0) stop("FATAL ERROR: No matrices generated at all!")
  counts <- do.call(cbind, counts_list)

} else {
  stop("This script is currently optimized for per_sample mode.")
}

common <- intersect(colnames(obj), colnames(counts))
msg(paste0("Cells successfully extracted: ", ncol(counts)))
if (length(common) == 0) stop("No overlap between original Seurat cells and newly extracted counts!")

# Ensure perfect order alignment
counts <- counts[, colnames(obj), drop = FALSE]

msg("Step 4: Add new assay and save")
assay_u <- CreateChromatinAssay(counts = counts, genome = "mm10")

obj[[assay_out]] <- assay_u
DefaultAssay(obj) <- assay_out
saveRDS(obj, out_rds)
msg(paste0("Successfully saved to: ", out_rds))
msg("Pipeline Complete. <")
