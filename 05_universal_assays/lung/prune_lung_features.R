library(Seurat)
library(Signac)

## Lung is re-quantified against the universal peak set in SeuratObject.R but comes out
## with 14 extra peak features (667,473 vs the 667,459 in Kidney/Aorta/Tcells), an
## artefact of the chromosome-set difference between the source fragment files. Prune
## Lung's feature space to exactly match the Kidney reference so all four objects are
## peak-aligned for downstream merging and comparison.

ref_rds  <- "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds"
lung_rds <- "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal.rds"
out_rds  <- "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds"

## reference feature set from Kidney
kidney <- readRDS(ref_rds)
master_features <- rownames(kidney)
rm(kidney); gc()

## keep only the peaks present in the Kidney reference
lung <- readRDS(lung_rds)
lung_pruned <- lung[master_features, ]

message("Pruned Lung: ", nrow(lung_pruned), " peaks x ", ncol(lung_pruned), " cells")
if (nrow(lung_pruned) != length(master_features)) {
  stop("Feature counts do not match the Kidney reference.")
}

saveRDS(lung_pruned, out_rds)
message("Saved: ", out_rds)
