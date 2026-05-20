# ==============================================================================
# Kidney scATAC-seq Workflow (H5 Optimized Version)
# Paper: Muto et al. (Science Advances 2024)
# Input: .h5 matrix + .tsv.gz fragments
# Samples: 8 Samples (40 excluded)
# ==============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(harmony)
library(Matrix)
library(hdf5r)

set.seed(1234)

# ----------------------------
# 1) Annotation Setup (mm10)
# ----------------------------
print("Setting up genome annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ----------------------------
# 2) Sample Info & Data Loading
# ----------------------------
base_dir <- "/QRISdata/Q8448/Mouse_disease_data/Kidney/cellranger_unpacked_data"

samples <- data.frame(
  sample_id = c(
    "SRR27367347", "SRR27367344", "SRR27367332", # Sham
    "SRR27367351", "SRR27367349", "SRR27367346", # Day14
    "SRR27367331", "SRR27367330", "SRR27367340"  # Day42 (New Sample!)
  ),
  condition = c(
    "Sham", "Sham", "Sham",
    "Day14", "Day14", "Day14",
    "Day42", "Day42", "Day42" # °7,/ Day42
  ),
  group = c(
    "Control", "Control", "Control",
    "Experiment", "Experiment", "Experiment",
    "Experiment", "Experiment", "Experiment"
  ),
  stringsAsFactors = FALSE
)

obj_list <- list()


for (i in 1:nrow(samples)) {
  
  cur_id <- samples$sample_id[i]
  
  folder_name <- paste0(cur_id, "_Kidney_atac")
  
  print(paste0("[", i, "/", nrow(samples), "] Processing: ", cur_id))
  
  outs_dir <- file.path(base_dir, folder_name, folder_name, "outs")
  
  h5_file   <- file.path(outs_dir, "filtered_peak_bc_matrix.h5")
  frag_file <- file.path(outs_dir, "fragments.tsv.gz")
  meta_path <- file.path(outs_dir, "singlecell.csv")
  
  if(file.exists(h5_file)) {
    counts <- Read10X_h5(filename = h5_file)
  } else {
    warning(paste("H5 file not found for", cur_id))
    next
  }
  
  # metadata
  if(file.exists(meta_path)){
    metadata <- read.csv(meta_path, header = TRUE, row.names = 1)
  } else {
    metadata <- NULL
    warning(paste("Metadata missing for", cur_id))
  }
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag_file,
    min.cells = 10,
    min.features = 200,
    annotation = annotations
  )
  
  cur_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata,
    project = cur_id
  )
  
  cur_obj$dataset <- cur_id
  cur_obj$condition <- samples$condition[i]
  cur_obj$group <- samples$group[i]
  
  # ----------------------------
  # 3) QC Calculation
  # ----------------------------
  # print(paste("  - Calculating QC..."))
  cur_obj <- NucleosomeSignal(cur_obj)
  cur_obj <- TSSEnrichment(cur_obj, fast = FALSE)
  
  cur_obj$peak_region_fragments <- cur_obj$nCount_peaks
  
  if ("passed_filters" %in% colnames(cur_obj@meta.data)) {
    cur_obj$pct_reads_in_peaks <- cur_obj$peak_region_fragments / cur_obj$passed_filters * 100
  } else {
    cur_obj$pct_reads_in_peaks <- NA
  }
  
  # blacklist_ratio
  if ("blacklist_region_fragments" %in% colnames(cur_obj@meta.data)) {
    cur_obj$blacklist_ratio <- cur_obj$blacklist_region_fragments / cur_obj$peak_region_fragments
  } else {
    cur_obj$blacklist_ratio <- 0 
  }
  
  # ----------------------------
  # 4) Doublet Finding
  # ----------------------------
  # print(paste("  - scDblFinder..."))
  mtx <- GetAssayData(cur_obj, layer = "counts")
  dbl_results <- suppressMessages(scDblFinder(mtx, verbose = FALSE))
  cur_obj$scDblFinder.class <- dbl_results$scDblFinder.class
  cur_obj$scDblFinder.score <- dbl_results$scDblFinder.score
  
  # ----------------------------
  # 5) Filtering (Per Paper Methods)
  # ----------------------------
  print(paste("  - Filtering", cur_id, "..."))
  print(paste("    Before:", ncol(cur_obj)))
  
  cur_obj <- subset(
    x = cur_obj,
    subset = peak_region_fragments > 2000 &    
      peak_region_fragments < 100000 &  
      pct_reads_in_peaks > 25 &         
      blacklist_ratio < 0.08 &          
      nucleosome_signal < 4 &           
      TSS.enrichment > 2 &              
      scDblFinder.class == "singlet"
  )
  
  print(paste("    After:", ncol(cur_obj)))
  
  obj_list[[i]] <- cur_obj
}
# ==============================================================================
# STEP 1: Calculate Missing QC Metrics for ALL Samples (In Memory)
# ==============================================================================

print("--- Updating QC Metrics for all samples in obj_list ---")

# Ensure list has names
names(obj_list) <- samples$sample_id

# Loop through list and update metrics
for (id in names(obj_list)) {
  
  # Retrieve object
  x <- obj_list[[id]]
  
  # 1. Nucleosome Signal
  if (!"nucleosome_signal" %in% colnames(x@meta.data)) {
    x <- NucleosomeSignal(x, verbose = FALSE)
  }
  
  # 2. TSS Enrichment
  if (!"TSS.enrichment" %in% colnames(x@meta.data)) {
    # This might take a moment, so we print progress
    print(paste("Calculating TSS for:", id))
    x <- TSSEnrichment(x, fast = FALSE, verbose = FALSE)
  }
  
  # 3. Peak Region Fragments (Rename for paper consistency)
  x$peak_region_fragments <- x$nCount_peaks
  
  # 4. pct_reads_in_peaks (FRiP)
  if ("passed_filters" %in% colnames(x@meta.data)) {
    x$pct_reads_in_peaks <- x$peak_region_fragments / x$passed_filters * 100
  } else {
    x$pct_reads_in_peaks <- NA # Handle missing metadata
  }
  
  # 5. Blacklist Ratio
  if ("blacklist_region_fragments" %in% colnames(x@meta.data)) {
    x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
  } else {
    x$blacklist_ratio <- 0
  }
  
  # Update the list with the calculated object
  obj_list[[id]] <- x
}

print("--- All metrics calculated. Ready to view! ---")


# ==============================================================================
# STEP 2: Define Interactive Viewer Function
# ==============================================================================

ViewQC <- function(sample_name) {
  
  # Check if sample exists
  if (!sample_name %in% names(obj_list)) {
    stop(paste("Error: Sample", sample_name, "not found in obj_list!"))
  }
  
  # Get object
  obj <- obj_list[[sample_name]]
  
  # Define the 5 metrics
  feats <- c("peak_region_fragments", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks")
  
  # Generate Plot
  p <- VlnPlot(
    object = obj,
    features = feats,
    ncol = 5,
    pt.size = 0.1,
    cols = "#F8766D"
  ) + 
    NoLegend() +
    plot_annotation(
      title = paste("QC Metrics for Single Sample:", sample_name),
      subtitle = paste("Total Cells:", ncol(obj)),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Force print the plot to the Viewer
  print(p)
}
# ----------------------------
# 6) Merge & Harmony Integration
# ----------------------------
print("Merging all samples...")
obj_list <- obj_list[!sapply(obj_list, is.null)]

kidney_merged <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  add.cell.ids = samples$sample_id
)

rm(obj_list); gc()

print("Running LSI & Harmony...")
DefaultAssay(kidney_merged) <- "peaks"
kidney_merged <- RunTFIDF(kidney_merged)
kidney_merged <- FindTopFeatures(kidney_merged, min.cutoff = 'q0')
kidney_merged <- RunSVD(kidney_merged)

# Run Harmony (Correcting Batch)
kidney_merged <- RunHarmony(
  object = kidney_merged,
  group.by.vars = "dataset",
  reduction.use = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE,
  verbose = TRUE
)

# ----------------------------
# 7) Clustering & UMAP
# ----------------------------
print("Running UMAP & Clustering...")
kidney_merged <- RunUMAP(kidney_merged, reduction = "harmony", dims = 2:30, seed.use = 1234)
kidney_merged <- FindNeighbors(kidney_merged, reduction = "harmony", dims = 2:30)
kidney_merged <- FindClusters(kidney_merged, algorithm = 3, resolution = 0.8, verbose = FALSE, random.seed = 1234)

# ----------------------------
# 8) Visualization
# ----------------------------
print("Calculating Gene Activities...")
gene.activities <- GeneActivity(kidney_merged)
kidney_merged[['RNA']] <- CreateAssayObject(counts = gene.activities)
kidney_merged <- NormalizeData(
  object = kidney_merged,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(kidney_merged$nCount_RNA)
)
# ==============================================================================
# STEP: Save the Integrated Object
# ==============================================================================
print("Saving merged object to disk...")

# Define output path (modify this path if needed)
output_file <- file.path(base_dir, "kidney_merged_harmony.rds")

# Save the object
saveRDS(kidney_merged, file = output_file)

print(paste("Object saved successfully to:", output_file))
print("Workflow Complete!")

DefaultAssay(kidney_merged) <- 'RNA'
p1 <- DimPlot(kidney_merged, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ggtitle("Kidney Clusters")
p2 <- DimPlot(kidney_merged, reduction = "umap", group.by = "condition") + ggtitle("Condition")

print(p1 + p2)

print("Workflow Complete!")

# ==============================================================================
# STEP 9: Cell Type Annotation (Marker Visualization)
# ==============================================================================

print("Generating Marker Gene DotPlot...")

# 1. Define Marker Gene List
marker_list_clean <- list(
  PCT = c("Slc5a12", "Slc13a3", "Acsm1", "Lrp2", "Slc34a1", "Slc7a13"),
  PST = c("Mep1b"), 
  DTL_ATL = c("Epha7", "Slc14a2"),
  TAL = c("Slc12a1", "Umod"),
  DCT_CNT = c("Slc12a3", "Trpm6", "Slc8a1"), 
  PC_URO = c("Scnn1b", "Aqp2", "Upk1b"), 
  IC = c("Slc26a4", "Kit"),
  PODO_PEC = c("Wt1", "Nphs1"),
  EC = c("Flt1", "Pecam1"),
  FIB = c("Pdgfrb"),
  LEUK = c("Ptprc", "Ikzf1", "Cd86"),
  LowQC_Injury = c("Havcr1", "Vcam1", "Slc5a8") # Suggest renaming to LowQC_Injury
)

# 2. Switch to RNA (Gene Activity) Assay
# Note: Must use 'kidney_merged' as this is your integrated object
DefaultAssay(kidney_merged) <- 'RNA' 

# 3. Generate DotPlot
p_dot <- DotPlot(
  object = kidney_merged, 
  features = marker_list_clean,
  cols = c("lightgrey", "red"),
  dot.scale = 8,
  cluster.idents = FALSE # Keep the numeric cluster order (0, 1, 2...)
) + 
  RotatedAxis() +
  ggtitle("Cell Type Markers (Gene Activity)")

# 4. Display and Save Plot
print(p_dot)

# Save the plot for manual annotation reference
ggsave("Kidney_Annotation_DotPlot.png", plot = p_dot, width = 14, height = 6)

print("DotPlot generated! Please check the image to start manual annotation.")

# ==============================================================================
# STEP 10: Rename Clusters Based on Marker Annotation
# ==============================================================================
DefaultAssay(kidney_merged) <- "RNA"
Idents(kidney_merged) <- "seurat_clusters"
clusters_to_check <- c("3", "6", "28")
for (clust in clusters_to_check) {
  print(paste("====== Cluster", clust, "Top 5 Markers ======"))
  cluster_markers <- FindMarkers(kidney_merged, 
                                 ident.1 = clust, 
                                 only.pos = TRUE, 
                                 logfc.threshold = 0.25)
  top5_markers <- head(cluster_markers, 5)
  print(top5_markers)
}
print("Renaming clusters to formal cell types...")
# ==============================================================================
# Find and export Top Markers for all Clusters (for rigorous cell annotation)
# ==============================================================================

print("=
 Calculating specific Marker genes for all Clusters (this may take a while)...")

# 1. Ensure Gene Activity (RNA) Assay is used, and revert to numeric clusters
DefaultAssay(kidney_merged) <- "RNA"
Idents(kidney_merged) <- "seurat_clusters"

# 2. Core function: Find Markers for all Clusters
# only.pos = TRUE: Only find upregulated (highly expressed) genes
# min.pct = 0.25: Gene must be expressed in at least 25% of cells
# logfc.threshold = 0.25: Log fold-change threshold
all_markers <- FindAllMarkers(
  object = kidney_merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 3. Extract Top 5 for each Cluster (can be changed to Top 10)
top5_all <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# 4. Print to console (for quick browsing)
print(top5_all, n = 100)

# 5. [CRITICAL] Save all Markers as a CSV file for manual review in Excel
output_csv <- file.path(base_dir, "All_Clusters_Top_Markers.csv")
write.csv(all_markers, file = output_csv, row.names = FALSE)

print(paste("All Markers have been calculated and saved to:", output_csv))
print("=¡ Tip: Open this CSV in Excel, filter for p_val_adj < 0.05, and cross-reference with databases like CellMarker.")

# Create a mapping vector from cluster IDs to annotated cell types
new.cluster.ids <- c(
  "0"  = "PCT",
  "1"  = "PCT",
  "2"  = "PCT",
  "3"  = "Injured_PT",
  "4"  = "EC",
  "5"  = "TAL",
  "6"  = "PST",
  "7"  = "PCT",
  "8"  = "PC_URO",
  "9"  = "PC_URO",
  "10" = "TAL",
  "11" = "FIB",
  "12" = "LEUK",
  "13" = "DCT_CNT",
  "14" = "TAL",
  "15" = "IC",
  "16" = "DTL_ATL",
  "17" = "PCT",
  "18" = "PC_URO",
  "19" = "TAL",
  "20" = "TAL",
  "21" = "PODO_PEC",
  "22" = "DTL_ATL",
  "23" = "DCT_CNT",
  "24" = "LEUK",
  "25" = "EC",
  "26" = "PC_URO",
  "27" = "FIB",
  "28" = "PCT",
  "29" = "PCT",
  "30" = "Injured_PT",
  "31" = "IC",
  "32" = "PCT"
)

# Ensure the current identities are set to Seurat clusters
Idents(kidney_merged) <- "seurat_clusters"

# Rename identities using the mapping
kidney_merged <- RenameIdents(kidney_merged, new.cluster.ids)

# [Key step] Save the annotated cell type labels into metadata as "cell_type"
kidney_merged$cell_type <- Idents(kidney_merged)

# Plot a UMAP with the new labels for a quick sanity check
p_annotated <- DimPlot(
  kidney_merged,
  reduction = "umap",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) + ggtitle("Annotated Kidney Cell Types")

print(p_annotated)

ggsave(
  "Kidney_UMAP_Annotated.png",
  plot = p_annotated,
  width = 10,
  height = 8
)

# ==============================================================================
# STEP 11: Split Fragment Files by Cell Type
# ==============================================================================
# Now we have the 'cell_type' column, so we can safely split the fragment files.

output_dir <- file.path(base_dir, "fragment_files_split_by_celltype")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

print(paste("Start splitting fragment files by cell type into:", output_dir))

# Make sure we are using the peaks assay fragments
DefaultAssay(kidney_merged) <- "peaks"

# Run fragment splitting
SplitFragments(
  object = kidney_merged,
  assay = "peaks",
  group.by = "cell_type",  # Split by our annotated labels
  outdir = output_dir,
  verbose = TRUE
)

print(" Finished splitting fragment files by cell type.")
print(" Entire annotation and splitting workflow complete!")

# Recommended: save the object again to persist the cell_type metadata
saveRDS(
  kidney_merged,
  file = file.path(base_dir, "kidney_merged_annotated.rds")
)

