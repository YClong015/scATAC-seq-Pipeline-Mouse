#!/usr/bin/env Rscript
# ==============================================================================
# T Cell scATAC-seq Workflow
# Follows the same pipeline as Kidney_scATAC(Combine).R:
#   QC → Doublet removal → Filter → LSI → Harmony → UMAP → Annotate by cell type → Split fragments
# Input: Tcells_Seurat_filtered.RData (single multiplexed experiment)
# Batch correction: Harmony by HTO condition (deMultliplex2_final_mapped)
# ==============================================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Mmusculus.v79)
  library(ggplot2)
  library(patchwork)
  library(scDblFinder)
  library(harmony)
  library(dplyr)
  library(Matrix)
})

set.seed(1234)

BASE_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/Tcells"
RDATA     <- file.path(BASE_DIR, "Tcells_Seurat_filtered.RData")
FRAG_FILE <- file.path(BASE_DIR, "atac_fragments.tsv.gz")
OUT_RDS   <- file.path(BASE_DIR, "tcells_processed.rds")

# ----------------------------
# 1) Genome Annotations (mm10)
# ----------------------------
print("Setting up mm10 genome annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ----------------------------
# 2) Load Data
# ----------------------------
print("Loading T cell Seurat object from RData...")
local_env <- new.env()
load(RDATA, envir = local_env)

# Auto-detect Seurat object name
obj_name <- Filter(function(n) inherits(get(n, envir = local_env), "Seurat"),
                   ls(local_env))[1]
if (is.na(obj_name)) stop("No Seurat object found in RData file.")
tcells <- get(obj_name, envir = local_env)
rm(local_env); gc()
print(paste("Loaded object:", obj_name, "| Total cells:", ncol(tcells)))

# Ensure ATAC assay is active and fragment path is current
DefaultAssay(tcells) <- "ATAC"
Fragments(tcells)    <- NULL
Fragments(tcells)    <- CreateFragmentObject(path = FRAG_FILE,
                                             cells = colnames(tcells),
                                             validate.fragments = FALSE)
Annotation(tcells)   <- annotations
print("Fragment path and annotations updated.")

# ----------------------------
# 3) QC Metrics (visualisation only — object is pre-filtered)
# ----------------------------
print("Calculating QC metrics for visualisation...")
tcells <- ATACqc(tcells, functions = "NucleosomeSignal")
suppressWarnings(tcells <- TSSEnrichment(tcells, fast = FALSE))

tcells$peak_region_fragments <- tcells$nCount_ATAC

if ("blacklist_region_fragments" %in% colnames(tcells@meta.data)) {
  tcells$blacklist_ratio <- tcells$blacklist_region_fragments /
                            tcells$peak_region_fragments
} else {
  tcells$blacklist_ratio <- 0
}

# Only include features with non-NA values
candidate_features <- c("peak_region_fragments", "TSS.enrichment",
                         "nucleosome_signal", "blacklist_ratio")
qc_features <- candidate_features[sapply(candidate_features, function(f) {
  f %in% colnames(tcells@meta.data) && !all(is.na(tcells@meta.data[[f]]))
})]
print(paste("QC features available:", paste(qc_features, collapse = ", ")))

p_qc <- VlnPlot(
  tcells,
  features  = qc_features,
  ncol      = length(qc_features),
  pt.size   = 0,
  group.by  = "deMultliplex2_final_mapped"
) + plot_annotation(title = "T cells — QC metrics")

ggsave(file.path(BASE_DIR, "Tcells_QC.pdf"),
       p_qc, width = 4 * length(qc_features), height = 5)
ggsave(file.path(BASE_DIR, "Tcells_QC.png"),
       p_qc, width = 4 * length(qc_features), height = 5, dpi = 300)
print(paste("QC plot saved. Total cells:", ncol(tcells)))

# ----------------------------
# 4) Doublet Detection & Removal
# ----------------------------
print("Running scDblFinder for doublet detection...")
mtx <- GetAssayData(tcells, assay = "ATAC", layer = "counts")
dbl <- suppressMessages(scDblFinder(mtx, verbose = FALSE))
tcells$scDblFinder.class <- dbl$scDblFinder.class
tcells$scDblFinder.score <- dbl$scDblFinder.score
print(table(tcells$scDblFinder.class))

tcells <- subset(tcells, subset = scDblFinder.class == "singlet")
print(paste("Cells after doublet removal:", ncol(tcells)))

# QC violin plots (after filtering)
p_qc_post <- VlnPlot(
  tcells,
  features = c("peak_region_fragments", "TSS.enrichment",
               "nucleosome_signal", "pct_reads_in_peaks"),
  ncol = 4, pt.size = 0, group.by = "deMultliplex2_final_mapped"
) + plot_annotation(title = "T cells — QC metrics (after filtering)")

ggsave(file.path(BASE_DIR, "Tcells_QC_after_filter.pdf"),
       p_qc_post, width = 16, height = 5)

# ----------------------------
# 6) LSI Dimensionality Reduction
# ----------------------------
print("Running LSI (TF-IDF + SVD)...")
DefaultAssay(tcells) <- "ATAC"
tcells <- RunTFIDF(tcells)
tcells <- FindTopFeatures(tcells, min.cutoff = "q0")
tcells <- RunSVD(tcells)

# Check correlation of each LSI component with sequencing depth (skip dim 1)
p_lsi <- DepthCor(tcells)
ggsave(file.path(BASE_DIR, "Tcells_LSI_DepthCor.pdf"), p_lsi, width = 6, height = 4)
print("LSI depth correlation plot saved — dim 1 should be excluded (high correlation).")

# ----------------------------
# 7) Harmony Batch Correction (by HTO condition)
# ----------------------------
print("Running Harmony batch correction by HTO condition (deMultliplex2_final_mapped)...")
tcells <- RunHarmony(
  object       = tcells,
  group.by.vars = "deMultliplex2_final_mapped",
  reduction.use = "lsi",
  assay.use     = "ATAC",
  project.dim   = FALSE,
  verbose       = TRUE
)

# ----------------------------
# 8) UMAP & Clustering
# ----------------------------
print("Running UMAP and clustering on Harmony-corrected LSI (dims 2:30)...")
tcells <- RunUMAP(tcells, reduction = "harmony", dims = 2:30, seed.use = 1234,
                  reduction.name = "umap")
tcells <- FindNeighbors(tcells, reduction = "harmony", dims = 2:30)
tcells <- FindClusters(tcells, algorithm = 3, resolution = 0.6,
                       verbose = FALSE, random.seed = 1234)

# Quick UMAP previews
p1 <- DimPlot(tcells, reduction = "umap", label = TRUE, pt.size = 0.3) +
      ggtitle("Clusters (Harmony LSI UMAP)")
p2 <- DimPlot(tcells, reduction = "umap",
              group.by = "deMultliplex2_final_mapped", pt.size = 0.3) +
      ggtitle("HTO Condition")

ggsave(file.path(BASE_DIR, "Tcells_UMAP_clusters_vs_condition.pdf"),
       p1 + p2, width = 14, height = 6)
print("UMAP preview saved.")

# ----------------------------
# 9) Gene Activity Matrix
# ----------------------------
print("Calculating Gene Activity matrix...")
gene.activities <- GeneActivity(tcells, extend.upstream = 2000,
                                extend.downstream = 0)
tcells[["RNA"]] <- CreateAssayObject(counts = gene.activities)
tcells <- NormalizeData(tcells, assay = "RNA",
                        normalization.method = "LogNormalize",
                        scale.factor = median(tcells$nCount_RNA))

# ----------------------------
# 10) Find Cluster Markers
# ----------------------------
print("Finding top marker genes per cluster (downsampled for speed)...")
DefaultAssay(tcells) <- "RNA"
tcells_small <- subset(tcells, downsample = 300)

all_markers <- FindAllMarkers(
  tcells_small,
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25
)

top5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(all_markers,
          file.path(BASE_DIR, "Tcells_All_Cluster_Markers.csv"),
          row.names = FALSE)
print("Marker genes saved to: Tcells_All_Cluster_Markers.csv")
print(top5, n = 100)

# DotPlot of T cell marker genes for manual annotation reference
marker_list <- list(
  Naive_T    = c("Sell", "Ccr7", "Tcf7", "S1pr1", "Il7r"),
  CD4_T      = c("Cd4", "Cd40lg", "Cxcr5"),
  CD8_T      = c("Cd8a", "Cd8b1", "Gzmb", "Prf1"),
  Treg       = c("Foxp3", "Ikzf2", "Il2ra"),
  Effector_T = c("Tbx21", "Ifng", "Klrg1", "Bhlhe40"),
  Memory_T   = c("Itga4", "Cxcr3", "Il18r1"),
  CD8_Eff    = c("Zeb2", "Cx3cr1"),
  NK         = c("Klrb1c", "Nkg7", "Klra7"),
  Cycling_T  = c("Mki67", "Top2a", "Rgcc"),
  B_cell     = c("Ebf1", "Cd19", "Ms4a1")
)

p_dot <- DotPlot(tcells, features = marker_list,
                 cols = c("lightgrey", "red"), dot.scale = 6) +
  RotatedAxis() +
  ggtitle("T cell marker gene activity per cluster")

ggsave(file.path(BASE_DIR, "Tcells_Annotation_DotPlot.pdf"),
       p_dot, width = 16, height = 6)
print("DotPlot saved — use this to assign cell type labels below.")

# ==============================================================================
# STEP 10: Manual Cell Type Annotation
# Fill in the mapping below AFTER reviewing the DotPlot and marker CSV.
# Current placeholders based on marker gene patterns from prior analysis.
# ==============================================================================
# Key markers per cluster (from prior analysis in Fig5_Tcells_annotation.R):
#   S1pr1        → Naïve T cell
#   Ikzf2/Helios → Treg
#   Cd8a/Zeb2    → CD8+ effector
#   Hsph1        → Activated T cell
#   Rgcc         → Cycling T cell
#   Klra7        → NK cell
#   Itga4        → Memory T cell
#   Ebf1         → B cell (contaminant — will be removed)

new.cluster.ids <- c(
  "0"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1 (DotPlot)
  "1"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1, Il7r (DotPlot)
  "2"  = "Naive_T",            # Sell, Ccr7 dominant; no Ebf1/Cd19/Ms4a1
  "3"  = "Tfh_like_T",         # Cxcr5 (follicular helper marker)
  "4"  = "Naive_T",            # Sell, Ccr7, S1pr1 dominant (DotPlot)
  "5"  = "Effector_CD8_T",     # Cd8b1, Tbx21, Klrg1, Bhlhe40, Cx3cr1
  "6"  = "Naive_CD8_T",        # Sell, S1pr1, Cd8b1 (DotPlot)
  "7"  = "Cytotoxic_CD8_T",    # Klrg1, Cx3cr1, Klrb1c (terminally differentiated)
  "8"  = "CD8_Eff",            # Zeb2, Cx3cr1
  "9"  = "Treg",               # Foxp3++, Ikzf2++ (DotPlot — confirmed Treg)
  "10" = "Naive_CD8_T",        # Sell, Ccr7, Tcf7, Cd8a, Cd8b1 (DotPlot)
  "11" = "Low_quality",        # Broad signal across all markers — likely doublets
  "12" = "B_cell",             # Ebf1++, Cd19++, Ms4a1++ — contamination
  "13" = "Memory_CD8_T",       # Itga4++, Cx3cr1, Zeb2
  "14" = "NK",                 # Nkg7++, Klra7++ (DotPlot)
  "15" = "Memory_CD8_T",       # Zeb2++, Il18r1, Itga4
  "16" = "CD8_Eff",            # Zeb2++ dominant
  "17" = "Naive_T"             # Sell, Ccr7, Tcf7 (DotPlot)
)

Idents(tcells) <- "seurat_clusters"
tcells <- RenameIdents(tcells, new.cluster.ids)
tcells$cell_type <- as.character(Idents(tcells))

# Remove contaminating clusters
remove_types <- c("B_cell", "Low_quality")
n_before <- ncol(tcells)
tcells   <- subset(tcells, subset = cell_type %in% remove_types, invert = TRUE)
message(sprintf("Removed %d contaminating cells (%s); %d T cells remaining",
                n_before - ncol(tcells),
                paste(remove_types, collapse = ", "),
                ncol(tcells)))

# Final annotated UMAP
p_ann <- DimPlot(tcells, reduction = "umap", group.by = "cell_type",
                 label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.3) +
         ggtitle("T cells — Harmony LSI UMAP (cell type annotation)") +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(BASE_DIR, "Tcells_UMAP_Annotated.pdf"),
       p_ann, width = 10, height = 8)
ggsave(file.path(BASE_DIR, "Tcells_UMAP_Annotated.png"),
       p_ann, width = 10, height = 8, dpi = 300)
print("Annotated UMAP saved.")

# Cell type counts
print(table(tcells$cell_type))
write.csv(as.data.frame(table(tcells$cell_type)),
          file.path(BASE_DIR, "Tcells_CellType_Counts.csv"),
          row.names = FALSE)

# ----------------------------
# 11) Save Annotated Object
# ----------------------------
print(paste("Saving annotated object to:", OUT_RDS))
saveRDS(tcells, file = OUT_RDS)

# ----------------------------
# 12) Split Fragments by Cell Type
# ----------------------------
print("Splitting fragment files by cell type...")
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
print(paste("Fragment files split by cell type saved to:", split_dir))
print("Pipeline complete.")
