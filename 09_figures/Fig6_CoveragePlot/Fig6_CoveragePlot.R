#!/usr/bin/env Rscript
# CoveragePlot - loads pre-linked objects from CoveragePlot_ready/
# Run Fig6_link_fragments.R first if ready objects don't exist.

library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)

# Build mm10 annotation once and reuse
message("Loading mm10 gene annotation ...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"   # ensure chr1/chr2 format
genome(annotations) <- "mm10"
message("Annotation loaded: ", length(annotations), " features")

ready.dir <- "/QRISdata/Q8448/Mouse_disease_data/CoveragePlot_ready"
out.dir   <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures/CoveragePlots"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

tissue.configs <- list(
  Tcells = list(
    rds       = file.path(ready.dir, "Tcells_ready.rds"),
    label_col = "cell_type",
    gene      = "Cd3e"
  ),
  Aorta = list(
    rds       = file.path(ready.dir, "Aorta_ready.rds"),
    label_col = "cell_type",
    gene      = "Acta2"
  ),
  Lung = list(
    rds       = file.path(ready.dir, "Lung_ready.rds"),
    label_col = "cell_type",
    gene      = "Sftpc"
  ),
  Kidney = list(
    rds       = file.path(ready.dir, "Kidney_ready.rds"),
    label_col = "cell_type",
    gene      = "Epcam"
  )
)

for (tis in names(tissue.configs)) {
  cfg <- tissue.configs[[tis]]

  if (!file.exists(cfg$rds)) {
    message(tis, ": ready object not found. Run Fig6_link_fragments.R first.")
    next
  }

  message("\n=== ", tis, " | gene: ", cfg$gene, " ===")
  obj <- readRDS(cfg$rds)
  DefaultAssay(obj) <- "peaks_universal"
  if (is.null(Annotation(obj[["peaks_universal"]])))
    Annotation(obj[["peaks_universal"]]) <- annotations
  # Add Tcells cell type annotation from cluster mapping
  if (tis == "Tcells") {
    tcell.annotations <- c(
      "0" = "Effector T cell",   "1"  = "Naive T cell",
      "2" = "B cell",            "3"  = "Treg",
      "4" = "CD8+ T cell",       "5"  = "Activated T cell",
      "6" = "CD8+ effector",     "7"  = "Cycling T cell",
      "8" = "Innate-like T cell","9"  = "NK cell",
      "10" = "Memory T cell",    "11" = "Memory T cell",
      "12" = "Memory T cell",    "13" = "Memory T cell",
      "14" = "CD8+ effector"
    )
    obj$cell_type <- unname(tcell.annotations[as.character(obj$seurat_clusters)])
  }
  Idents(obj) <- cfg$label_col
  message("  Cells = ", ncol(obj))

  # Downsample to max 200 cells per cell type to speed up CoveragePlot on NFS
  set.seed(1234)
  cells.keep <- unlist(lapply(levels(Idents(obj)), function(ct) {
    ct.cells <- WhichCells(obj, idents = ct)
    if (length(ct.cells) > 200) sample(ct.cells, 200) else ct.cells
  }))
  obj <- subset(obj, cells = cells.keep)
  message("  Downsampled to ", ncol(obj), " cells (max 200 per cell type)")

  p <- tryCatch({
    CoveragePlot(
      object            = obj,
      region            = cfg$gene,
      annotation        = TRUE,
      peaks             = TRUE,
      extend.upstream   = 1000,
      extend.downstream = 1000
    )
  }, error = function(e) {
    message("  ERROR: ", conditionMessage(e)); NULL
  })

  if (!is.null(p)) {
    pdf(file.path(out.dir, paste0("CoveragePlot_", tis, "_", cfg$gene, ".pdf")),
        width = 10, height = 7)
    print(p)
    dev.off()

    png(file.path(out.dir, paste0("CoveragePlot_", tis, "_", cfg$gene, ".png")),
        width = 3000, height = 2100, res = 300)
    print(p)
    dev.off()

    message("  Saved.")
  }

  rm(obj); gc()
}

message("\nDone. Outputs -> ", out.dir)
