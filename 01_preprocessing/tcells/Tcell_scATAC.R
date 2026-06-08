#!/usr/bin/env Rscript
# T-cell scATAC-seq preprocessing: QC, doublet removal, LSI, Harmony, UMAP,
# initial annotation. Same workflow as Kidney_scATAC_Combine.R.


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

## mm10 annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

## load (auto-detect the Seurat object inside the RData)
local_env <- new.env()
load(RDATA, envir = local_env)
obj_name <- Filter(function(n) inherits(get(n, envir = local_env), "Seurat"),
                   ls(local_env))[1]
if (is.na(obj_name)) stop("No Seurat object found in RData file.")
tcells <- get(obj_name, envir = local_env)
rm(local_env); gc()

DefaultAssay(tcells) <- "ATAC"
Fragments(tcells)    <- NULL
Fragments(tcells)    <- CreateFragmentObject(path = FRAG_FILE,
                                             cells = colnames(tcells),
                                             validate.fragments = FALSE)
Annotation(tcells)   <- annotations

## QC metrics (object is pre-filtered; these are for visualisation only)
tcells <- ATACqc(tcells, functions = "NucleosomeSignal")
suppressWarnings(tcells <- TSSEnrichment(tcells, fast = FALSE))
tcells$peak_region_fragments <- tcells$nCount_ATAC
if ("blacklist_region_fragments" %in% colnames(tcells@meta.data)) {
  tcells$blacklist_ratio <- tcells$blacklist_region_fragments /
                            tcells$peak_region_fragments
} else {
  tcells$blacklist_ratio <- 0
}
candidate_features <- c("peak_region_fragments", "TSS.enrichment",
                         "nucleosome_signal", "blacklist_ratio")
qc_features <- candidate_features[sapply(candidate_features, function(f) {
  f %in% colnames(tcells@meta.data) && !all(is.na(tcells@meta.data[[f]]))
})]
p_qc <- VlnPlot(tcells, features = qc_features, ncol = length(qc_features),
                pt.size = 0, group.by = "deMultliplex2_final_mapped") +
  plot_annotation(title = "T cells - QC metrics")
ggsave(file.path(BASE_DIR, "Tcells_QC.pdf"),
       p_qc, width = 4 * length(qc_features), height = 5)
ggsave(file.path(BASE_DIR, "Tcells_QC.png"),
       p_qc, width = 4 * length(qc_features), height = 5, dpi = 300)

## doublet removal
mtx <- GetAssayData(tcells, assay = "ATAC", layer = "counts")
dbl <- suppressMessages(scDblFinder(mtx, verbose = FALSE))
tcells$scDblFinder.class <- dbl$scDblFinder.class
tcells$scDblFinder.score <- dbl$scDblFinder.score
tcells <- subset(tcells, subset = scDblFinder.class == "singlet")
message("Cells after doublet removal: ", ncol(tcells))

p_qc_post <- VlnPlot(tcells,
  features = c("peak_region_fragments", "TSS.enrichment",
               "nucleosome_signal", "pct_reads_in_peaks"),
  ncol = 4, pt.size = 0, group.by = "deMultliplex2_final_mapped") +
  plot_annotation(title = "T cells - QC metrics (after filtering)")
ggsave(file.path(BASE_DIR, "Tcells_QC_after_filter.pdf"),
       p_qc_post, width = 16, height = 5)

## LSI (TF-IDF + SVD); dim 1 correlates with depth and is excluded downstream
DefaultAssay(tcells) <- "ATAC"
tcells <- RunTFIDF(tcells)
tcells <- FindTopFeatures(tcells, min.cutoff = "q0")
tcells <- RunSVD(tcells)
ggsave(file.path(BASE_DIR, "Tcells_LSI_DepthCor.pdf"), DepthCor(tcells),
       width = 6, height = 4)

## Harmony batch correction by HTO condition
tcells <- RunHarmony(tcells, group.by.vars = "deMultliplex2_final_mapped",
                     reduction.use = "lsi", assay.use = "ATAC",
                     project.dim = FALSE, verbose = TRUE)

## UMAP + clustering on Harmony-corrected LSI (dims 2:30)
tcells <- RunUMAP(tcells, reduction = "harmony", dims = 2:30, seed.use = 1234,
                  reduction.name = "umap")
tcells <- FindNeighbors(tcells, reduction = "harmony", dims = 2:30)
tcells <- FindClusters(tcells, algorithm = 3, resolution = 0.6,
                       verbose = FALSE, random.seed = 1234)
p1 <- DimPlot(tcells, reduction = "umap", label = TRUE, pt.size = 0.3) +
      ggtitle("Clusters (Harmony LSI UMAP)")
p2 <- DimPlot(tcells, reduction = "umap",
              group.by = "deMultliplex2_final_mapped", pt.size = 0.3) +
      ggtitle("HTO Condition")
ggsave(file.path(BASE_DIR, "Tcells_UMAP_clusters_vs_condition.pdf"),
       p1 + p2, width = 14, height = 6)

## gene activity (for marker-based annotation)
gene.activities <- GeneActivity(tcells, extend.upstream = 2000, extend.downstream = 0)
tcells[["RNA"]] <- CreateAssayObject(counts = gene.activities)
tcells <- NormalizeData(tcells, assay = "RNA",
                        normalization.method = "LogNormalize",
                        scale.factor = median(tcells$nCount_RNA))

## per-cluster markers (downsampled for speed) + marker DotPlot for annotation
DefaultAssay(tcells) <- "RNA"
tcells_small <- subset(tcells, downsample = 300)
all_markers <- FindAllMarkers(tcells_small, only.pos = TRUE,
                              min.pct = 0.1, logfc.threshold = 0.25)
write.csv(all_markers, file.path(BASE_DIR, "Tcells_All_Cluster_Markers.csv"),
          row.names = FALSE)

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
  RotatedAxis() + ggtitle("T cell marker gene activity per cluster")
ggsave(file.path(BASE_DIR, "Tcells_Annotation_DotPlot.pdf"),
       p_dot, width = 16, height = 6)

## initial annotation from the DotPlot (provisional; revised downstream)
new.cluster.ids <- c(
  "0"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1
  "1"  = "Naive_T",            # Sell, Ccr7, Tcf7, S1pr1, Il7r
  "2"  = "Naive_T",            # Sell, Ccr7 dominant
  "3"  = "Tfh_like_T",         # Cxcr5
  "4"  = "Naive_T",            # Sell, Ccr7, S1pr1
  "5"  = "Effector_CD8_T",     # Cd8b1, Tbx21, Klrg1, Bhlhe40, Cx3cr1
  "6"  = "Naive_CD8_T",        # Sell, S1pr1, Cd8b1
  "7"  = "Cytotoxic_CD8_T",    # Klrg1, Cx3cr1, Klrb1c
  "8"  = "CD8_Eff",            # Zeb2, Cx3cr1
  "9"  = "Treg",               # Foxp3, Ikzf2
  "10" = "Naive_CD8_T",        # Sell, Ccr7, Tcf7, Cd8a, Cd8b1
  "11" = "Low_quality",        # broad signal across markers (likely doublets)
  "12" = "B_cell",             # Ebf1, Cd19, Ms4a1 (contamination)
  "13" = "Memory_CD8_T",       # Itga4, Cx3cr1, Zeb2
  "14" = "NK",                 # Nkg7, Klra7
  "15" = "Memory_CD8_T",       # Zeb2, Il18r1, Itga4
  "16" = "CD8_Eff",            # Zeb2
  "17" = "Naive_T"             # Sell, Ccr7, Tcf7
)

Idents(tcells) <- "seurat_clusters"
tcells <- RenameIdents(tcells, new.cluster.ids)
tcells$cell_type <- as.character(Idents(tcells))

## drop contaminating clusters
remove_types <- c("B_cell", "Low_quality")
n_before <- ncol(tcells)
tcells   <- subset(tcells, subset = cell_type %in% remove_types, invert = TRUE)
message(sprintf("Removed %d contaminating cells (%s); %d remaining",
                n_before - ncol(tcells), paste(remove_types, collapse = ", "),
                ncol(tcells)))

p_ann <- DimPlot(tcells, reduction = "umap", group.by = "cell_type",
                 label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.3) +
         ggtitle("T cells - Harmony LSI UMAP (cell type annotation)") +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(BASE_DIR, "Tcells_UMAP_Annotated.pdf"), p_ann, width = 10, height = 8)
ggsave(file.path(BASE_DIR, "Tcells_UMAP_Annotated.png"), p_ann, width = 10, height = 8, dpi = 300)

write.csv(as.data.frame(table(tcells$cell_type)),
          file.path(BASE_DIR, "Tcells_CellType_Counts.csv"), row.names = FALSE)

saveRDS(tcells, file = OUT_RDS)
