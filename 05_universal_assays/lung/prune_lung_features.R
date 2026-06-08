library(Seurat)
library(Signac)

## Lung re-quantification yields 14 extra peaks (667,473 vs 667,459) from a chromosome-set
## difference; prune to the Kidney reference so all four objects are peak-aligned.

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
