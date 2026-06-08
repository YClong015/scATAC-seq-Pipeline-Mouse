#!/usr/bin/env Rscript
# Run ONCE: link fragment files to the peaks_universal assay, save per-tissue
# {Tissue}_ready.rds for CoveragePlot. Cell -> fragment mapping differs by tissue
# (Kidney SRR prefix, Lung BGI run id, Aorta merge suffix, Tcells single file).

library(Signac)
library(Seurat)

save.dir <- "/QRISdata/Q8448/Mouse_disease_data/CoveragePlot_ready"
dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)

base.lung  <- "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac"
base.kid   <- "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data"
base.aorta <- "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac"

## Tissue configs
tissue.rds <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
)

# Kidney: SRR sample ID -> fragment file
kidney.map <- c(
  SRR27367347 = file.path(base.kid, "SRR27367347_Kidney_atac/SRR27367347_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367344 = file.path(base.kid, "SRR27367344_Kidney_atac/SRR27367344_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367332 = file.path(base.kid, "SRR27367332_Kidney_atac/SRR27367332_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367351 = file.path(base.kid, "SRR27367351_Kidney_atac/SRR27367351_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367349 = file.path(base.kid, "SRR27367349_Kidney_atac/SRR27367349_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367346 = file.path(base.kid, "SRR27367346_Kidney_atac/SRR27367346_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367331 = file.path(base.kid, "SRR27367331_Kidney_atac/SRR27367331_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367330 = file.path(base.kid, "SRR27367330_Kidney_atac/SRR27367330_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367340 = file.path(base.kid, "SRR27367340_Kidney_atac/SRR27367340_Kidney_atac/outs/fragments.tsv.gz")
)

# Lung: orig.ident -> fragment file (confirmed from CNGB metadata)
lung.map <- c(
  Control_F2 = file.path(base.lung, "CL100168054_L01/outs/fragments.tsv.gz"),
  Control_M1 = file.path(base.lung, "CL100167942_L01/outs/fragments.tsv.gz"),
  Case_F1    = file.path(base.lung, "CL100168078_L02/outs/fragments.tsv.gz"),
  Case_F3    = file.path(base.lung, "CL100168054_L02/outs/fragments.tsv.gz"),
  Case_M2    = file.path(base.lung, "CL100167942_L02/outs/fragments.tsv.gz"),
  Case_M3    = file.path(base.lung, "CL100168078_L01/outs/fragments.tsv.gz")
)

# Aorta: Seurat merge suffix _1/_2 -> fragment file (confirmed by barcode scan)
aorta.map <- c(
  "1" = file.path(base.aorta, "SRR21686724_Aorta_atac/outs/fragments.tsv.gz"),  # Control
  "2" = file.path(base.aorta, "SRR21686722_Aorta_atac/outs/fragments.tsv.gz")   # Challenge
)

tcell.frag <- "/QRISdata/Q8448/Mouse_disease_data/Tcells/atac_fragments.tsv.gz"

## Main
for (tis in names(tissue.rds)) {
  out.rds <- file.path(save.dir, paste0(tis, "_ready.rds"))

  if (file.exists(out.rds)) {
    message(tis, ": already exists, skipping. Delete to re-run.")
    next
  }

  message("\n=== ", tis, " ===")
  obj <- readRDS(tissue.rds[[tis]])
  DefaultAssay(obj) <- "peaks_universal"
  message("  Loaded. Cells = ", ncol(obj))

  # Skip if already linked
  existing <- tryCatch(Fragments(obj[["peaks_universal"]]), error = function(e) list())
  if (length(existing) > 0) {
    message("  Fragments already linked (", length(existing), "). Saving as-is.")
    saveRDS(obj, out.rds)
    message("  Saved -> ", out.rds)
    rm(obj); gc(); next
  }

  all.cells <- colnames(obj)
  frag.list <- list()

  ## Kidney
  if (tis == "Kidney") {
    for (sid in names(kidney.map)) {
      fp       <- kidney.map[sid]
      cell.idx <- grepl(sid, all.cells, fixed = TRUE)
      if (!any(cell.idx)) { message("  [SKIP] no cells: ", sid); next }
      sub.cells <- all.cells[cell.idx]
      # Raw barcode: strip "SRR_ID_" prefix
      raw.bc         <- sub("^[^_]+_", "", sub.cells)
      names(raw.bc)  <- sub.cells
      frag.list[[sid]] <- CreateFragmentObject(
        path = fp, cells = raw.bc, validate.fragments = FALSE
      )
      message("  Linked: ", sid, " (", length(sub.cells), " cells)")
    }
  }

  ## Lung
  if (tis == "Lung") {
    for (sid in names(lung.map)) {
      fp       <- lung.map[sid]
      cell.idx <- obj$orig.ident == sid
      if (!any(cell.idx)) { message("  [SKIP] no cells: ", sid); next }
      sub.cells <- all.cells[cell.idx]
      # Raw barcode: strip "Condition_Sample_" prefix (first 2 fields)
      raw.bc        <- sub("^[^_]+_[^_]+_", "", sub.cells)
      names(raw.bc) <- sub.cells
      frag.list[[sid]] <- CreateFragmentObject(
        path = fp, cells = raw.bc, validate.fragments = FALSE
      )
      message("  Linked: ", sid, " (", length(sub.cells), " cells)")
    }
  }

  ## Aorta
  if (tis == "Aorta") {
    for (sfx in names(aorta.map)) {
      fp       <- aorta.map[sfx]
      cell.idx <- grepl(paste0("_", sfx, "$"), all.cells)
      if (!any(cell.idx)) { message("  [SKIP] no cells with suffix _", sfx); next }
      sub.cells <- all.cells[cell.idx]
      # Raw barcode: strip trailing "_N" Seurat merge suffix
      raw.bc        <- sub("_[0-9]+$", "", sub.cells)
      names(raw.bc) <- sub.cells
      frag.list[[sfx]] <- CreateFragmentObject(
        path = fp, cells = raw.bc, validate.fragments = FALSE
      )
      label <- if (sfx == "1") "SRR21686724 Control" else "SRR21686722 Challenge"
      message("  Linked: ", label, " (", length(sub.cells), " cells)")
    }
  }

  ## Tcells
  if (tis == "Tcells") {
    if (!file.exists(tcell.frag)) {
      message("  ERROR: Tcells fragment file not found"); rm(obj); gc(); next
    }
    # Barcodes are already raw (no prefix/suffix)
    raw.bc        <- all.cells
    names(raw.bc) <- all.cells
    frag.list[["Tcells"]] <- CreateFragmentObject(
      path = tcell.frag, cells = raw.bc, validate.fragments = FALSE
    )
    message("  Linked: Tcells (", length(all.cells), " cells)")
  }

  if (length(frag.list) == 0) {
    message("  ERROR: no fragments linked for ", tis); rm(obj); gc(); next
  }

  Fragments(obj[["peaks_universal"]]) <- frag.list
  saveRDS(obj, out.rds)
  message("  Saved -> ", out.rds)
  rm(obj); gc()
}

message("\nDone. Ready objects in: ", save.dir)
message("Now run Fig6_CoveragePlot.R to generate plots.")
