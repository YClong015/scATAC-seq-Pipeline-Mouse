#!/usr/bin/env Rscript
# ============================================================
# Run ONCE: link fragment files to peaks_universal assay
# and save per-tissue objects ready for CoveragePlot.
#
# Output: /QRISdata/.../CoveragePlot_ready/{Tissue}_ready.rds
#
# Confirmed cell → fragment file mappings:
#   Kidney : cell name contains SRR prefix  → per-sample match
#   Lung   : orig.ident (Control_F2 etc.)   → BGI run ID map
#   Aorta  : Seurat merge suffix _1/_2      → SRR21686724/SRR21686722
#   Tcells : single fragment file           → all cells
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
})

save_dir <- "/QRISdata/Q8448/Mouse_disease_data/CoveragePlot_ready"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

base_lung  <- "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac"
base_kid   <- "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data"
base_aorta <- "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac"

# ---------------------------------------------------------------
# Tissue configs
# ---------------------------------------------------------------
tissue_rds <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
)

# Kidney: SRR sample ID → fragment file
kidney_map <- c(
  SRR27367347 = file.path(base_kid, "SRR27367347_Kidney_atac/SRR27367347_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367344 = file.path(base_kid, "SRR27367344_Kidney_atac/SRR27367344_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367332 = file.path(base_kid, "SRR27367332_Kidney_atac/SRR27367332_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367351 = file.path(base_kid, "SRR27367351_Kidney_atac/SRR27367351_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367349 = file.path(base_kid, "SRR27367349_Kidney_atac/SRR27367349_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367346 = file.path(base_kid, "SRR27367346_Kidney_atac/SRR27367346_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367331 = file.path(base_kid, "SRR27367331_Kidney_atac/SRR27367331_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367330 = file.path(base_kid, "SRR27367330_Kidney_atac/SRR27367330_Kidney_atac/outs/fragments.tsv.gz"),
  SRR27367340 = file.path(base_kid, "SRR27367340_Kidney_atac/SRR27367340_Kidney_atac/outs/fragments.tsv.gz")
)

# Lung: orig.ident → fragment file (confirmed from CNGB metadata)
lung_map <- c(
  Control_F2 = file.path(base_lung, "CL100168054_L01/outs/fragments.tsv.gz"),
  Control_M1 = file.path(base_lung, "CL100167942_L01/outs/fragments.tsv.gz"),
  Case_F1    = file.path(base_lung, "CL100168078_L02/outs/fragments.tsv.gz"),
  Case_F3    = file.path(base_lung, "CL100168054_L02/outs/fragments.tsv.gz"),
  Case_M2    = file.path(base_lung, "CL100167942_L02/outs/fragments.tsv.gz"),
  Case_M3    = file.path(base_lung, "CL100168078_L01/outs/fragments.tsv.gz")
)

# Aorta: Seurat merge suffix _1/_2 → fragment file (confirmed by barcode scan)
aorta_map <- c(
  "1" = file.path(base_aorta, "SRR21686724_Aorta_atac/outs/fragments.tsv.gz"),  # Control
  "2" = file.path(base_aorta, "SRR21686722_Aorta_atac/outs/fragments.tsv.gz")   # Challenge
)

tcell_frag <- "/QRISdata/Q8448/Mouse_disease_data/Tcells/atac_fragments.tsv.gz"

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
for (tis in names(tissue_rds)) {
  out_rds <- file.path(save_dir, paste0(tis, "_ready.rds"))

  if (file.exists(out_rds)) {
    message(tis, ": already exists, skipping. Delete to re-run.")
    next
  }

  message("\n=== ", tis, " ===")
  obj <- readRDS(tissue_rds[[tis]])
  DefaultAssay(obj) <- "peaks_universal"
  message("  Loaded. Cells = ", ncol(obj))

  # Skip if already linked
  existing <- tryCatch(Fragments(obj[["peaks_universal"]]), error = function(e) list())
  if (length(existing) > 0) {
    message("  Fragments already linked (", length(existing), "). Saving as-is.")
    saveRDS(obj, out_rds)
    message("  Saved -> ", out_rds)
    rm(obj); gc(); next
  }

  all_cells <- colnames(obj)
  frag_list <- list()

  # ----- Kidney -----
  if (tis == "Kidney") {
    for (sid in names(kidney_map)) {
      fp       <- kidney_map[sid]
      cell_idx <- grepl(sid, all_cells, fixed = TRUE)
      if (!any(cell_idx)) { message("  [SKIP] no cells: ", sid); next }
      sub_cells <- all_cells[cell_idx]
      # Raw barcode: strip "SRR_ID_" prefix
      raw_bc         <- sub("^[^_]+_", "", sub_cells)
      names(raw_bc)  <- sub_cells
      frag_list[[sid]] <- CreateFragmentObject(
        path = fp, cells = raw_bc, validate.fragments = FALSE
      )
      message("  Linked: ", sid, " (", length(sub_cells), " cells)")
    }
  }

  # ----- Lung -----
  if (tis == "Lung") {
    for (sid in names(lung_map)) {
      fp       <- lung_map[sid]
      cell_idx <- obj$orig.ident == sid
      if (!any(cell_idx)) { message("  [SKIP] no cells: ", sid); next }
      sub_cells <- all_cells[cell_idx]
      # Raw barcode: strip "Condition_Sample_" prefix (first 2 fields)
      raw_bc        <- sub("^[^_]+_[^_]+_", "", sub_cells)
      names(raw_bc) <- sub_cells
      frag_list[[sid]] <- CreateFragmentObject(
        path = fp, cells = raw_bc, validate.fragments = FALSE
      )
      message("  Linked: ", sid, " (", length(sub_cells), " cells)")
    }
  }

  # ----- Aorta -----
  if (tis == "Aorta") {
    for (sfx in names(aorta_map)) {
      fp       <- aorta_map[sfx]
      cell_idx <- grepl(paste0("_", sfx, "$"), all_cells)
      if (!any(cell_idx)) { message("  [SKIP] no cells with suffix _", sfx); next }
      sub_cells <- all_cells[cell_idx]
      # Raw barcode: strip trailing "_N" Seurat merge suffix
      raw_bc        <- sub("_[0-9]+$", "", sub_cells)
      names(raw_bc) <- sub_cells
      frag_list[[sfx]] <- CreateFragmentObject(
        path = fp, cells = raw_bc, validate.fragments = FALSE
      )
      label <- if (sfx == "1") "SRR21686724 Control" else "SRR21686722 Challenge"
      message("  Linked: ", label, " (", length(sub_cells), " cells)")
    }
  }

  # ----- Tcells -----
  if (tis == "Tcells") {
    if (!file.exists(tcell_frag)) {
      message("  ERROR: Tcells fragment file not found"); rm(obj); gc(); next
    }
    # Barcodes are already raw (no prefix/suffix)
    raw_bc        <- all_cells
    names(raw_bc) <- all_cells
    frag_list[["Tcells"]] <- CreateFragmentObject(
      path = tcell_frag, cells = raw_bc, validate.fragments = FALSE
    )
    message("  Linked: Tcells (", length(all_cells), " cells)")
  }

  if (length(frag_list) == 0) {
    message("  ERROR: no fragments linked for ", tis); rm(obj); gc(); next
  }

  Fragments(obj[["peaks_universal"]]) <- frag_list
  saveRDS(obj, out_rds)
  message("  Saved -> ", out_rds)
  rm(obj); gc()
}

message("\nDone. Ready objects in: ", save_dir)
message("Now run Fig_CoveragePlot.R to generate plots.")
