library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)

# ==============================================================================
# Define file paths
# ==============================================================================
files <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_universal_v5.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_universal_v5.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_v5.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal_v5.rds"
)

# Output path
OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/Integrated"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Read and prepare objects
# ==============================================================================
obj_list <- list()

for (tissue in names(files)) {
  message(paste0("Loading: ", tissue, "..."))
  x <- readRDS(files[[tissue]])

  # Key point: Switch to Universal Assay
  DefaultAssay(x) <- "peaks_universal"

  # Remove old assays to save memory (keep only universal)
  # If you want to keep RNA assay, you can, but it's recommended to keep it light during ATAC integration
  for (assay in names(x@assays)) {
    if (assay != "peaks_universal") {
      x[[assay]] <- NULL
    }
  }

  # Add Tissue metadata for Harmony
  x$Tissue <- tissue
  
  obj_list[[tissue]] <- x
}

# ==============================================================================
# 3. Merge all objects into one Seurat object
# ==============================================================================
message("Merging all objects...")

combined <- merge(
  x = obj_list[["Kidney"]],
  y = list(
    obj_list[["Aorta"]],
    obj_list[["Lung"]],
    obj_list[["Tcells"]]
  ),
  add.cell.ids = c("Kidney", "Aorta", "Lung", "Tcells"),
  project = "Mouse_Atlas"
)


if (!"Tissue" %in% colnames(combined@meta.data)) {
  message("Warning: Tissue column missing after merge, re-assigning...")
  combined$Tissue <- sub("_.*", "", colnames(combined))
}

message("=== Verification After Explicit Merge ===")
print(table(combined$Tissue))
message("Total Cells: ", ncol(combined))

# Directly extract real counts from the underlying matrix, without relying on potentially corrupted metadata
counts_matrix <- GetAssayData(combined, assay = "peaks_universal", layer = "counts")
combined$real_nCount <- colSums(counts_matrix)

# Print the average count for each tissue to confirm there is indeed data
message("=== Average Universal Peak Counts by Tissue ===")
print(tapply(combined$real_nCount, combined$Tissue, mean))

# Perform cleanup (using the real data we just calculated)
message("Cleaning up low count cells (< 100)...")
combined <- subset(combined, subset = real_nCount > 100)

message("=== After Filtering ===")
print(table(combined$Tissue))

# Cleanup: Remove the temporary counts matrix to free memory
rm(counts_matrix)
gc()

# ==============================================================================
# 4. Standard processing workflow (LSI)
# ==============================================================================
message("Starting (TF-IDF + SVD)...")

DefaultAssay(combined) <- "peaks_universal"

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)

# Initialize UMAP based on LSI
message("Calculating Uncorrected UMAP...")
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:30, reduction.name = "umap.lsi", reduction.key = "lsiUMAP_")

gc()
message(sprintf("Memory after LSI UMAP: %.1f GB used",
                sum(gc()[,2]) * 8 / 1024))

# ==============================================================================
# 5. Harmony integration (removing batch effects)
# ==============================================================================
message("Starting Harmony integration (removing Tissue batch effects)...")

# Harmony will mix cells from different tissues based on the 'Tissue' column
combined <- RunHarmony(
  object = combined,
  group.by.vars = "Tissue", 
  reduction.use = "lsi",             
  reduction.save = "harmony",        
  assay.use = "peaks_universal",
  project.dim = FALSE,
  dims.use = 2:30
)
# BASED ON HARMONY EMBEDDINGS, CALCULATE UMAP
message("Calculating Integrated UMAP...")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:29, reduction.name = "umap.harmony", reduction.key = "harmonyUMAP_")

gc()
message(sprintf("Memory after Harmony UMAP: %.1f GB used",
                sum(gc()[,2]) * 8 / 1024))

# Find neighbors and clusters based on Harmony embeddings
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:29)
combined <- FindClusters(combined, algorithm = 3, resolution = 0.5) 

# ==============================================================================
# 6. Save and check results
# ==============================================================================
# Save the final object
save_path <- file.path(OUT_DIR, "All_Tissues_Integrated.rds")
message(paste0("Saving final object to: ", save_path))
saveRDS(combined, save_path)

# Plot comparison (save as PDF)
pdf(file.path(OUT_DIR, "Integration_Check.pdf"), width = 12, height = 6)

p1 <- DimPlot(combined, reduction = "umap.lsi", group.by = "Tissue") + ggtitle("Before Integration (LSI)")
p2 <- DimPlot(combined, reduction = "umap.harmony", group.by = "Tissue") + ggtitle("After Integration (Harmony)")

print(p1 + p2)
dev.off()

message("Done! Check the PDF for integration results and the RDS file for the integrated object.")


