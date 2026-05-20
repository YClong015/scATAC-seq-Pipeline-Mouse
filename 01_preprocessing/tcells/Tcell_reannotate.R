#!/usr/bin/env Rscript
# Re-apply corrected cell type annotations to existing tcells_processed.rds
# Run this instead of re-running the full Tcell_scATAC.R pipeline

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(dplyr)
})

BASE_DIR <- "/QRISdata/Q8448/Mouse_disease_data/Tcells"
RDS_IN   <- file.path(BASE_DIR, "tcells_processed.rds")
RDS_OUT  <- RDS_IN

message("Loading tcells_processed.rds...")
tcells <- readRDS(RDS_IN)
message(sprintf("Loaded: %d cells, clusters: %s",
                ncol(tcells),
                paste(sort(unique(tcells$seurat_clusters)), collapse = ",")))

# ------------------------------------------------------------------------------
# Corrected annotations (based on DotPlot review)
# ------------------------------------------------------------------------------
new.cluster.ids <- c(
  "0"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1
  "1"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1, Il7r
  "2"  = "Naive_T",            # Sell, Ccr7 dominant; no Ebf1/Cd19/Ms4a1
  "3"  = "Tfh_like_T",         # Cxcr5
  "4"  = "Naive_T",            # Sell, Ccr7, S1pr1 dominant
  "5"  = "Effector_CD8_T",     # Cd8b1, Tbx21, Klrg1, Bhlhe40, Cx3cr1
  "6"  = "Naive_CD8_T",        # Sell, S1pr1, Cd8b1
  "7"  = "Cytotoxic_CD8_T",    # Klrg1, Cx3cr1, Klrb1c
  "8"  = "CD8_Eff",            # Zeb2, Cx3cr1
  "9"  = "Treg",               # Foxp3++, Ikzf2++
  "10" = "Naive_CD8_T",        # Sell, Ccr7, Tcf7, Cd8a, Cd8b1
  "11" = "Low_quality",        # Broad signal across all markers
  "12" = "B_cell",             # Ebf1++, Cd19++, Ms4a1++
  "13" = "Memory_CD8_T",       # Itga4++, Cx3cr1, Zeb2
  "14" = "NK",                 # Nkg7++, Klra7++
  "15" = "Memory_CD8_T",       # Zeb2++, Il18r1, Itga4
  "16" = "CD8_Eff",            # Zeb2++ dominant
  "17" = "Naive_T"             # Sell, Ccr7, Tcf7
)

Idents(tcells) <- "seurat_clusters"
tcells <- RenameIdents(tcells, new.cluster.ids)
tcells$cell_type <- as.character(Idents(tcells))

message("\nCell type distribution (before removing contamination):")
print(table(tcells$cell_type))

# Remove contaminating clusters
remove_types <- c("B_cell", "Low_quality")
n_before <- ncol(tcells)
tcells   <- subset(tcells, subset = cell_type %in% remove_types, invert = TRUE)
message(sprintf("\nRemoved %d contaminating cells (%s); %d cells remaining",
                n_before - ncol(tcells),
                paste(remove_types, collapse = ", "),
                ncol(tcells)))

message("\nFinal cell type counts:")
print(table(tcells$cell_type))

write.csv(as.data.frame(table(tcells$cell_type)),
          file.path(BASE_DIR, "Tcells_CellType_Counts.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Split fragments by corrected cell type
# ------------------------------------------------------------------------------
message("\nSplitting fragment files by corrected cell type...")
split_dir <- file.path(BASE_DIR, "fragment_files_split_by_celltype")
dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(tcells) <- "ATAC"
SplitFragments(
  object   = tcells,
  assay    = "ATAC",
  group.by = "cell_type",
  outdir   = split_dir,
  verbose  = TRUE
)
message(sprintf("Fragment files saved to: %s", split_dir))

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
message(sprintf("\nSaving to: %s", RDS_OUT))
saveRDS(tcells, file = RDS_OUT)
message("Done.")
