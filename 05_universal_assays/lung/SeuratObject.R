#!/usr/bin/env Rscript
# Generic per-tissue universal-peak re-quantification script.
# Reads any per-tissue Seurat .rds with a `peaks` assay, re-quantifies cells
# against a supplied universal peak BED, writes a new object with a
# `peaks_universal` assay.
#
# Used by:
#   05_universal_assays/lung/SeuratObject.slurm  (Lung-specific arguments)
#
# Tissue-specific barcode-cleaning logic is embedded in clean_barcodes().
# The Kidney / Aorta / Tcell sibling scripts have different clean_barcodes().

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(future)
})

plan("multicore", workers = 4)
options(future.globals.maxSize = 50 * 1024^3)

msg <- function(x) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] "), x, "\n")

get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste0("Missing value after ", key))
  args[[hit + 1]]
}

# Lung barcode cleaner: strip "Condition_Sample_" prefix to recover BGI raw barcode.
# Example: Control_F2_CELL771_N1 -> CELL771_N1
clean_barcodes <- function(x) sub("^.*?(CELL.*)", "\\1", x)

load_peaks_bed_gz <- function(path) {
  con <- gzfile(path, "rt")
  df <- read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  close(con)
  if (ncol(df) < 3) stop("Universal bed file has <3 columns.")
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  bad <- is.na(df$start) | is.na(df$end) | df$end <= df$start
  df <- df[!bad, , drop = FALSE]
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1L, end = df$end)
  )
}

args <- commandArgs(trailingOnly = TRUE)

obj_path    <- get_arg(args, "--obj", "")
obj_type    <- get_arg(args, "--objtype", "rds")
up_bed      <- get_arg(args, "--up", "")
out_rds     <- get_arg(args, "--out", "")
assay_in    <- get_arg(args, "--assay_in", "peaks")
assay_out   <- get_arg(args, "--assay_out", "peaks_universal")
mode        <- get_arg(args, "--mode", "per_sample")
sample_key  <- get_arg(args, "--sample_key", "dataset")
frag_tpl    <- get_arg(args, "--frag_tpl", "")

if (obj_path == "" || up_bed == "" || out_rds == "") {
  stop("Error: Missing required arguments: --obj, --up, or --out")
}

msg("=== Step 1: Load Seurat Object ===")
obj <- readRDS(obj_path)
if (!assay_in %in% names(obj@assays)) {
  stop(paste0("Error: Assay not found in object: ", assay_in))
}
DefaultAssay(obj) <- assay_in
msg(paste0("Loaded object: ", obj_path))
msg(paste0("Total cells: ", ncol(obj)))
msg(paste0("Using sample key column: ", sample_key))

if (!sample_key %in% colnames(obj@meta.data)) {
  stop(paste0("Error: The sample_key '", sample_key, "' was not found!"))
}

msg("=== Step 2: Load Universal Peaks (BED) ===")
peaks <- load_peaks_bed_gz(up_bed)
msg(paste0("Loaded universal peaks: ", length(peaks)))

msg("=== Step 3: Build FeatureMatrix ===")
counts_list <- list()

if (mode == "per_sample") {
  if (frag_tpl == "") stop("Error: --mode per_sample requires --frag_tpl")

  sample_ids <- sort(unique(obj[[sample_key]][, 1]))
  msg(paste0("Found ", length(sample_ids), " samples to process."))

  for (id in sample_ids) {
    frag_path <- gsub("\\{id\\}", id, frag_tpl)

    if (!file.exists(frag_path)) {
      msg(paste0("WARNING: Fragment file missing for ID '", id, "': ", frag_path))
      next
    }

    cells_pref <- colnames(obj)[obj[[sample_key]][, 1] == id]
    if (length(cells_pref) == 0) next

    cells_raw <- clean_barcodes(cells_pref)
    msg(paste0("Processing Sample: ", id, " ... (", length(cells_pref), " cells)"))

    frag <- CreateFragmentObject(path = frag_path, validate.fragments = FALSE)
    mat <- FeatureMatrix(fragments = frag, features = peaks, cells = cells_raw)

    if (ncol(mat) == 0) {
      msg(paste0("    WARNING: No cells matched in sample ", id))
      next
    }

    if (ncol(mat) == length(cells_pref)) {
      colnames(mat) <- cells_pref
    } else {
      match_idx <- match(colnames(mat), cells_raw)
      colnames(mat) <- cells_pref[match_idx]
    }

    counts_list[[id]] <- mat
  }

  if (length(counts_list) == 0) stop("FATAL: No matrices generated.")
  counts <- do.call(cbind, counts_list)

} else {
  stop("Error: This script is configured for --mode per_sample only.")
}

common <- intersect(colnames(obj), colnames(counts))
msg(paste0("Total matched cells: ", length(common), " / ", ncol(obj)))
if (length(common) == 0) stop("FATAL: No overlap.")

counts <- counts[, colnames(obj), drop = FALSE]

msg("=== Step 4: Create Assay and Save ===")
assay_u <- CreateChromatinAssay(counts = counts, genome = "mm10")
obj[[assay_out]] <- assay_u
DefaultAssay(obj) <- assay_out

msg(paste0("New assay dimensions: ", nrow(obj[[assay_out]]), " x ", ncol(obj[[assay_out]])))

saveRDS(obj, out_rds)
msg(paste0("Success! Object saved to: ", out_rds))
