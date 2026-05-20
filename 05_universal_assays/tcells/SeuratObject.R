#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(future)
})

# Enable parallel processing (Optional: speeds up FeatureMatrix)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50 * 1024^3) # Set limit to 50GB

# Helper function for logging
msg <- function(x) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] "), x, "\n")

# Helper function to parse command line arguments
get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste0("Missing value after ", key))
  args[[hit + 1]]
}

# Function to clean cell barcodes (remove suffixes or prefixes to match fragment file)
strip_by_id <- function(x, id) {
  # 1. Try standard prefix matching (e.g., SampleID_Barcode)
  pref <- paste0(id, "_")
  if (all(startsWith(x, pref))) {
    return(substring(x, nchar(pref) + 1))
  }
  
  # 2. If prefix doesn't match, try stripping Seurat suffixes (e.g., _1, _2)
  # Regex: Find "_number" at the end of the string and replace with empty string
  cleaned <- sub("_\\d+$", "", x)
  return(cleaned)
}

# Function to load the universal peaks BED file
load_peaks_bed_gz <- function(path) {
  con <- gzfile(path, "rt")
  df <- read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  close(con)

  if (ncol(df) < 3) stop("Universal bed file has <3 columns.")
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")

  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)

  # Filter bad ranges
  bad <- is.na(df$start) | is.na(df$end) | df$end <= df$start
  df <- df[!bad, , drop = FALSE]

  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1L, end = df$end)
  )
  gr
}

# ==============================================================================
# Main Execution Logic
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
obj_path    <- get_arg(args, "--obj", "")
obj_type    <- get_arg(args, "--objtype", "rds")
up_bed      <- get_arg(args, "--up", "")
out_rds     <- get_arg(args, "--out", "")
assay_in    <- get_arg(args, "--assay_in", "peaks")
assay_out   <- get_arg(args, "--assay_out", "peaks_universal")
mode        <- get_arg(args, "--mode", "per_sample")
sample_key  <- get_arg(args, "--sample_key", "dataset") 
frag_file   <- get_arg(args, "--frag_file", "")
frag_tpl    <- get_arg(args, "--frag_tpl", "")

# Validate required arguments
if (obj_path == "" || up_bed == "" || out_rds == "") {
  stop("Error: Missing required arguments: --obj, --up, or --out")
}

msg("===============================================")
msg("Step 1: Load Seurat Object")
msg("===============================================")

if (obj_type == "rds") {
  obj <- readRDS(obj_path)
} else {
  stop("Error: --objtype must be 'rds'")
}

# Check if input assay exists
if (!assay_in %in% names(obj@assays)) {
  stop(paste0("Error: Assay not found in object: ", assay_in))
}
DefaultAssay(obj) <- assay_in

msg(paste0("Loaded object: ", obj_path))
msg(paste0("Total cells: ", ncol(obj)))
msg(paste0("Using sample key column: ", sample_key))

# Validate sample_key in metadata
if (!sample_key %in% colnames(obj@meta.data)) {
    stop(paste0("Error: The sample_key '", sample_key, "' was not found in the object metadata!"))
}

msg("===============================================")
msg("Step 2: Load Universal Peaks (BED)")
msg("===============================================")

peaks <- load_peaks_bed_gz(up_bed)
msg(paste0("Loaded universal peaks: ", length(peaks)))

msg("===============================================")
msg("Step 3: Build FeatureMatrix (Quantification)")
msg("===============================================")

counts_list <- list()

if (mode == "per_sample") {
  if (frag_tpl == "") stop("Error: --mode per_sample requires --frag_tpl")
  
  # Get unique sample IDs from the object
  sample_ids <- sort(unique(obj[[sample_key]][, 1]))
  msg(paste0("Found ", length(sample_ids), " samples to process."))
  print(sample_ids)

  for (id in sample_ids) {
    # Construct fragment path by replacing {id}
    frag_path <- gsub("\\{id\\}", id, frag_tpl)

    if (!file.exists(frag_path)) {
      msg(paste0("WARNING: Fragment file missing for ID '", id, "': ", frag_path))
      next
    }

    # Identify cells belonging to this sample
    cells_pref <- colnames(obj)[obj[[sample_key]][, 1] == id]
    if (length(cells_pref) == 0) next

    # Clean cell names to match the Fragment file format (remove suffixes)
    cells_raw <- strip_by_id(cells_pref, id)

    msg(paste0("Processing Sample: ", id, " ... (", length(cells_pref), " cells)"))
    
    # Create Fragment Object (skip validation for speed)
    frag <- CreateFragmentObject(path = frag_path, validate.fragments = FALSE)

    # Calculate Feature Matrix
    mat <- FeatureMatrix(
      fragments = frag,
      features = peaks,
      cells = cells_raw
    )

    if (ncol(mat) == 0) {
      msg(paste0("    WARNING: No cells matched in sample ", id))
      next
    }

    # === CRITICAL FIX: Restore original Seurat cell names ===
    # This ensures that the new matrix matches the existing object exactly
    if (ncol(mat) == length(cells_pref)) {
        colnames(mat) <- cells_pref
    } else {
        # Fallback: if FeatureMatrix dropped some cells, match them back
        match_idx <- match(colnames(mat), cells_raw)
        colnames(mat) <- cells_pref[match_idx]
    }
    
    counts_list[[id]] <- mat
  }

  if (length(counts_list) == 0) stop("FATAL: No matrices generated. Check your paths and IDs!")
  
  # Combine matrices from all samples
  counts <- do.call(cbind, counts_list)
  
} else {
  stop("Error: This script is configured for --mode per_sample only.")
}

# Final check for cell overlap
common <- intersect(colnames(obj), colnames(counts))
msg(paste0("Total matched cells: ", length(common), " / ", ncol(obj)))

if (length(common) == 0) stop("FATAL: No overlap between object cells and calculated counts.")

# Subset object to matched cells (usually keeps all cells)
obj <- subset(obj, cells = common)
counts <- counts[, common, drop = FALSE]

msg("===============================================")
msg("Step 4: Create Assay and Save")
msg("===============================================")

# Try to copy annotations from the old assay
old_annotations <- Annotation(obj[[assay_in]])

if (is.null(old_annotations)) {
  msg("WARNING: No gene annotations found in old assay. Creating assay without annotations.")
  assay_u <- CreateChromatinAssay(counts = counts, genome = "mm10")
} else {
  assay_u <- CreateChromatinAssay(counts = counts, genome = "mm10", annotation = old_annotations)
}

# Add the new assay to the object
obj[[assay_out]] <- assay_u
DefaultAssay(obj) <- assay_out

# Verify dimensions
msg(paste0("New assay dimensions: ", paste(dim(obj[[assay_out]]), collapse=" x ")))

# Save final object
saveRDS(obj, out_rds)
msg(paste0("Success! Object saved to: ", out_rds))
msg("Done.")
