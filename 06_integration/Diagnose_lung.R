library(Seurat)

# 1. Define file paths
files <- c(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds", 
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal.rds"
)

message("===============================================")
message("Starting Diagnostic Check for ATAC-seq Data Objects")
message("===============================================")

for (tissue in names(files)) {
  message("\nChecking: [ ", tissue, " ]")
  
  if (!file.exists(files[[tissue]])) {
    message("File does not exist, please check the path!")
    next
  }
  
  obj <- readRDS(files[[tissue]])
  
  # Check 1: Assay exists and is default
  current_assay <- DefaultAssay(obj)
  if (current_assay != "peaks_universal") {
    message("Warning: Default Assay is not peaks_universal, it is ", current_assay)
  }
  
  # Check 2: Dimensions (Number of Peaks and Cells)
  mat <- GetAssayData(obj, assay = "peaks_universal", layer = "counts")
  n_peaks <- nrow(mat)
  n_cells <- ncol(mat)
  message("Dimensions: ", n_peaks, " Peaks x ", n_cells, " Cells")
  
  # Check 3: Real Signal Amount (Prevent Empty Shell Replays)
  total_reads <- sum(mat, na.rm = TRUE)
  avg_reads <- total_reads / n_cells
  if (total_reads == 0) {
    message("Fatal Error: Matrix is empty (0 Reads)!")
  } else {
    message("Signal Health: Total Reads = ", total_reads, ", Average per Cell = ", round(avg_reads, 1))
  }
  
  # Check 4: Naming Convention Format
  first_peak <- head(rownames(mat), 1)
  first_cell <- head(colnames(mat), 1)
  message("Peak Naming Format: ", first_peak)
  message("Cell Naming Format: ", first_cell)
  
  rm(obj, mat)
  gc()
}

message("\n===============================================")
message("<Á Diagnostic Check Complete!")

#############################################################
library(Seurat)
library(Signac)

# 1. Load the Golden Standard (Kidney) to get the correct 667,459 peaks
message("Loading Kidney as reference template...")
kidney <- readRDS("/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds")
master_features <- rownames(kidney)
rm(kidney); gc()

# 2. Load the newly rescued Lung object (the one with 667,473 peaks)
message("Loading Lung for pruning...")
lung <- readRDS("/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal.rds")

# 3. Perform the pruning
message("Pruning 667,473 -> 667,459 peaks...")
# This keeps ONLY the peaks that exist in the Kidney master list
lung_pruned <- lung[master_features, ]

# 4. Final verification
message("=== Final Verification ===")
message("New Peak Count: ", nrow(lung_pruned))
message("Cell Count: ", ncol(lung_pruned))

if (nrow(lung_pruned) == 667459) {
  message("Alignment Successful!")
} else {
  stop("L Alignment Failed. Feature counts still do not match.")
}

# 5. Save as the specific file your Harmony script expects
out_path <- "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds"
message("Saving pruned object to: ", out_path)
saveRDS(lung_pruned, out_path)

message("Done! Lung is now ready for the Harmony pipeline.")
