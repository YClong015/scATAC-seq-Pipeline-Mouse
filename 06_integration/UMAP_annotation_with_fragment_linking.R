library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)

combined <- readRDS("/QRISdata/Q8448/Mouse_disease_data/Integrated/All_Tissues_Integrated.rds")

# Define the mapping of SampleID to Fragment file paths
fragment_map <- list(
  # Kidney samples (9)
  "SRR27367347" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367347_Kidney_atac/SRR27367347_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367344" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367344_Kidney_atac/SRR27367344_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367332" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367332_Kidney_atac/SRR27367332_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367351" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367351_Kidney_atac/SRR27367351_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367349" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367349_Kidney_atac/SRR27367349_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367346" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367346_Kidney_atac/SRR27367346_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367331" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367331_Kidney_atac/SRR27367331_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367330" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367330_Kidney_atac/SRR27367330_Kidney_atac/outs/fragments.tsv.gz",
  "SRR27367340" = "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data/SRR27367340_Kidney_atac/SRR27367340_Kidney_atac/outs/fragments.tsv.gz",
  
  # Lung samples (6)
  "CL100168054_L01" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100168054_L01/outs/fragments.tsv.gz",
  "CL100167942_L01" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100167942_L01/outs/fragments.tsv.gz",
  "CL100168054_L02" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100168054_L02/outs/fragments.tsv.gz",
  "CL100168078_L02" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100168078_L02/outs/fragments.tsv.gz",
  "CL100167942_L02" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100167942_L02/outs/fragments.tsv.gz",
  "CL100168078_L01" = "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac/CL100168078_L01/outs/fragments.tsv.gz",
  
  # Aorta samples (2)
  "SRR21686724" = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac/SRR21686724_Aorta_atac/outs/fragments.tsv.gz",
  "SRR21686722" = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac/SRR21686722_Aorta_atac/outs/fragments.tsv.gz",
  
  # Tcell samples (1)
  "Tcell_Single" = "/QRISdata/Q8448/Mouse_disease_data/Tcells/atac_fragments.tsv.gz"
)

#  Create a new list to store Fragment objects
message("Starting to create Fragment objects for each sample...")
new_fragments_list <- list()


fix_barcodes <- function(x) {
  x <- ifelse(grepl("^Aorta_", x), sub("_[0-9]+$", "", sub("^Aorta_", "", x)), x)
  x <- ifelse(grepl("^Kidney_", x), sub("^Kidney_[^_]+_", "", x), x)
  x <- ifelse(grepl("^Tcells_", x), sub("^Tcells_", "", x), x)
  x <- ifelse(grepl("^Lung_", x), sub(".*(CELL[0-9]+_N[0-9]+).*", "\\1", x), x)
  return(x)
}

for (sample_key in names(fragment_map)) {
  cells_in_sample <- colnames(combined)[grepl(sample_key, colnames(combined))]
  if (length(cells_in_sample) > 0) {
    raw_barcodes <- fix_barcodes(cells_in_sample)
    names(raw_barcodes) <- cells_in_sample
    
    fobj <- CreateFragmentObject(
      path = fragment_map[[sample_key]],
      cells = raw_barcodes,
      validate.fragments = FALSE 
    )
    new_fragments_list[[sample_key]] <- fobj
  }
}

Fragments(combined[["peaks_universal"]]) <- new_fragments_list
DefaultAssay(combined) <- 'peaks_universal'
gene.activities <- GeneActivity(combined)

combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
DefaultAssay(combined) <- 'RNA'
tapply(combined$nCount_RNA, combined$Tissue, summary)


message("1. Performing lightweight downsampling (up to 300 cells per cluster)...")
# This step finishes quickly and greatly reduces the matrix size
combined_small <- subset(combined, downsample = 300)

message("2. Finding marker genes with a fast workflow...")
# Use the downsampled subset for calculation
# Set only.pos = TRUE to find only upregulated genes
fast_markers <- FindAllMarkers(
  object = combined_small,
  only.pos = TRUE,
  min.pct = 0.1,          # Slightly relaxed threshold (ATAC-inferred activity
  # can be lower)
  logfc.threshold = 0.25
)

message("3. Extracting the top 10 genes for each cluster...")
# Use dplyr to group by cluster and select the top 10 by log2FC
top10_genes <- fast_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# 4. Save the results as a CSV file for easy viewing in Excel
output_file <- "/QRISdata/Q8448/Mouse_disease_data/Integrated/Top10_Markers.csv"
write.csv(top10_genes, file = output_file, row.names = FALSE)

message("Calculation complete! Please check the file: ", output_file)

# 5. Print the comparison UMAP between harmony or not
p_before <- DimPlot(combined, reduction = "umap.lsi", group.by = "Tissue", pt.size = 0.1, raster = FALSE) +
  ggtitle("Before Harmony (Unintegrated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p_after <- DimPlot(combined, reduction = "umap.harmony", group.by = "Tissue", pt.size = 0.1, raster = FALSE) +
  ggtitle("After Harmony (Integrated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf("/QRISdata/Q8448/Mouse_disease_data/Integrated/UMAP_Harmony_Comparison.pdf", width = 14, height = 6)
print(p_before | p_after)
dev.off()


Idents(combined) <- "seurat_clusters"
# Draw Harmony UMAP with cluster number
p_cluster_num <- DimPlot(
  object = combined,
  reduction = "umap.harmony",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.1,
  raster = FALSE
) +
  ggtitle("Harmony UMAP (Seurat Clusters)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_cluster_num)

# save pdf
pdf(
  "/QRISdata/Q8448/Mouse_disease_data/Integrated/"
  |> paste0("UMAP_Harmony_ClusterNumbers.pdf"),
  width = 10,
  height = 8
)
print(p_cluster_num)
dev.off()

message(" Harmony UMAP with cluster numbers saved.")

# Starr to annotate cluster
cluster_annotations <- c(
  "0" = "PT",
  "1" = "TAL",
  "2" = "Tcell", 
  "3" = "SMC",
  "4" = "DCT",
  "5" = "Urothelium",
  "6" = "PC",
  "7" = "Unassigned (VN/OR-high)",
  "8" = "Macrophages",
  "9" = "PT",
  "10" = "PT",
  "11" = "Lung-derived (Unassigned, OR-high)",
  "12" = "Endothelial",
  "13" = "PT",
  "14" = "PT",
  "15" = "TAL",
  "16" = "Tcell",
  "17" = "Fibroblasts",
  "18" = "Bcell",
  "19" = "Urothelium",
  "20" = "IC",
  "21" = "Pericytes",
  "22" = "Capillary Endothelial",
  "23" = "Lung-derived (Unassigned, OR-high)",
  "24" = "Doublets",
  "25" = "DC",
  "26" = "Endothelial",
  "27" = "Endothelial",
  "28" = "Podocytes",
  "29" = "TAL",
  "30" = "Basal Urothelium"
)

missing_clusters <- setdiff(
  levels(Idents(combined)),
  names(cluster_annotations)
)
extra_clusters <- setdiff(
  names(cluster_annotations),
  levels(Idents(combined))
)

print(missing_clusters)
print(extra_clusters)

combined <- RenameIdents(combined, cluster_annotations)

combined$cell_type <- Idents(combined)

message("Cluster annotation completed and saved to metadata.")
# -----------------------------
# Save annotated Seurat object
# -----------------------------

output_dir <- "/QRISdata/Q8448/Mouse_disease_data/Integrated"

# Keep a copy of the original numeric cluster labels (for safety)
combined$cluster_id <- as.character(combined$seurat_clusters)

# Save the fully annotated Seurat object
annotated_rds <- file.path(
  output_dir,
  "All_Tissues_Integrated_Annotated.rds"
)
saveRDS(combined, file = annotated_rds)

# Save metadata table (easy to inspect in Excel/R)
metadata_csv <- file.path(
  output_dir,
  "All_Tissues_Integrated_Annotated_Metadata.csv"
)
write.csv(combined[[]], file = metadata_csv, row.names = TRUE)

# Save the cluster annotation dictionary used in this run
anno_df <- data.frame(
  seurat_cluster = names(cluster_annotations),
  cell_type = unname(cluster_annotations),
  stringsAsFactors = FALSE
)

anno_csv <- file.path(
  output_dir,
  "Cluster_Annotation_Dictionary.csv"
)
write.csv(anno_df, file = anno_csv, row.names = FALSE)

# Save cell counts per annotated cell type
celltype_counts <- as.data.frame(table(combined$cell_type))
colnames(celltype_counts) <- c("cell_type", "n_cells")

celltype_counts_csv <- file.path(
  output_dir,
  "CellType_Counts.csv"
)
write.csv(celltype_counts, file = celltype_counts_csv, row.names = FALSE)

message("Annotated Seurat object saved: ", annotated_rds)
message("Metadata CSV saved: ", metadata_csv)
message("Annotation dictionary saved: ", anno_csv)
message("Cell type counts saved: ", celltype_counts_csv)

# Keep full annotated object unchanged
combined_annotated <- combined

# Remove suspicious / artifact-like clusters + doublets for DAR
clusters_to_remove <- c("7", "11", "23", "24")

combined_clean <- subset(
  combined_annotated,
  subset = !(seurat_clusters %in% clusters_to_remove)
)

# Quick check
table(combined_clean$seurat_clusters)
table(combined_clean$cell_type, useNA = "ifany")

# Save clean object for downstream DAR/enrichment
saveRDS(
  combined_clean,
  "/QRISdata/Q8448/Mouse_disease_data/Integrated/All_Tissues_Integrated_Annotated_Clean_for_DAR.rds"
)

message("Clean object for DAR/enrichment saved.")
