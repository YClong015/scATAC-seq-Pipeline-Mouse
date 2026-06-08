#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
})

msg <- function(x) cat(x, "\n")

get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste0("Missing value after ", key))
  args[[hit + 1]]
}

strip_by_id <- function(x, id) {
  pref <- paste0(id, "_")
  out <- x
  ok <- startsWith(x, pref)
  out[ok] <- substring(x[ok], nchar(pref) + 1)
  out[!ok] <- sub("^[^_]+_", "", x[!ok])
  out
}

load_peaks_bed_gz <- function(path) {
  con <- gzfile(path, "rt")
  df <- read.table(con, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
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
} else if (obj_type == "rdata") {
  e <- new.env()
  load(obj_path, envir = e)
  nm <- NULL
  for (x in ls(e)) {
    if (inherits(get(x, envir = e), "Seurat")) {
      nm <- x
      break
    }
  }
  if (is.null(nm)) stop("No Seurat object in RData.")
  obj <- get(nm, envir = e)
  rm(e)
} else {
  stop("objtype must be rds or rdata")
}

if (!assay_in %in% names(obj@assays)) {
  stop(paste0("Assay not found: ", assay_in))
}
DefaultAssay(obj) <- assay_in
msg(paste0("OK loaded: ", obj_path))
msg(paste0("Cells: ", ncol(obj)))

msg("Step 2: Load universal peaks (bed.gz)")
peaks <- load_peaks_bed_gz(up_bed)
msg(paste0("Peaks: ", length(peaks)))

msg("Step 3: Build FeatureMatrix on universal peaks")
counts_list <- list()

if (mode == "single") {
  if (frag_file == "") stop("single mode needs --frag_file")

  cells_pref <- colnames(obj)

  frag <- CreateFragmentObject(path = frag_file)

  mat <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    cells = cells_pref
  )

  if (ncol(mat) == 0) {
    cells_raw <- sub("^[^_]+_", "", cells_pref)
    mat <- FeatureMatrix(
      fragments = frag,
      features = peaks,
      cells = cells_raw
    )
    if (ncol(mat) == 0) stop("No cells matched fragments in single mode.")
    colnames(mat) <- cells_pref[match(colnames(mat), cells_raw)]
  }

  counts <- mat
} else if (mode == "per_sample") {
  if (frag_tpl == "") stop("per_sample mode needs --frag_tpl")
  if (!sample_key %in% colnames(obj@meta.data)) {
    stop(paste0("Missing sample_key: ", sample_key))
  }

  sample_ids <- sort(unique(obj[[sample_key]][, 1]))
  msg(paste0("Samples: ", length(sample_ids)))

  for (id in sample_ids) {
    frag_path <- gsub("\\{id\\}", id, frag_tpl)

    if (!file.exists(frag_path)) {
      stop(paste0("Missing fragments: ", frag_path))
    }

    cells_pref <- colnames(obj)[obj[[sample_key]][, 1] == id]
    if (length(cells_pref) == 0) next

    cells_raw <- strip_by_id(cells_pref, id)

    msg(paste0("  - ", id, " cells=", length(cells_pref)))

    frag <- CreateFragmentObject(path = frag_path)

    mat <- FeatureMatrix(
      fragments = frag,
      features = peaks,
      cells = cells_raw
    )

    if (ncol(mat) == 0) {
      msg(paste0("    WARNING: no cells matched in ", id))
      next
    }

    colnames(mat) <- paste0(id, "_", colnames(mat))
    counts_list[[id]] <- mat
  }

  if (length(counts_list) == 0) stop("No matrices generated.")
  counts <- do.call(cbind, counts_list)
} else {
  stop("mode must be single or per_sample")
}

common <- intersect(colnames(obj), colnames(counts))
msg(paste0("Cells in counts: ", ncol(counts)))
msg(paste0("Cells intersect: ", length(common)))
if (length(common) == 0) stop("No overlap between obj cells and counts.")

obj <- subset(obj, cells = common)
counts <- counts[, common, drop = FALSE]

msg("Step 4: Add new assay and save")
msg(paste0("Copying annotations from assay: ", assay_in))
old_annotations <- Annotation(obj[[assay_in]])

if (is.null(old_annotations)) {
  msg("WARNING: No annotations in old assay; creating assay without annotation.")
  assay_u <- CreateChromatinAssay(
    counts = counts,
    genome = "mm10"
  )
} else {
  msg(paste0("Found annotations for ", length(old_annotations), " ranges."))
  assay_u <- CreateChromatinAssay(
    counts = counts,
    genome = "mm10",
    annotation = old_annotations
  )
}

obj[[assay_out]] <- assay_u
DefaultAssay(obj) <- assay_out

saveRDS(obj, out_rds)
msg(paste0("OK saved: ", out_rds))
msg("Done.")

