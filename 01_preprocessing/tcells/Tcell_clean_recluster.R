#!/usr/bin/env Rscript
# ==============================================================================
# Tcell_clean_recluster.R
# Stage 1 of T-cell re-annotation:
#   - Remove contaminating non-T clusters confirmed by FindAllMarkers
#     (B cells: Tfh_like_T, CD8_Eff;  plasma: NK;  myeloid: Cytotoxic_CD8_T, Memory_CD8_T)
#   - Optional marker-score safety gate to catch residual stragglers
#   - Re-run the SAME LSI/Harmony/UMAP/clustering pipeline as Tcell_scATAC.R
#   - Output FindAllMarkers + canonical-marker DotPlot for re-annotation
#
# After running, INSPECT the outputs:
#   Tcells_clean_cluster_topmarkers.csv   -> top genes per new cluster
#   Tcells_clean_recheck_dotplot.pdf      -> confirm no residual B/plasma/myeloid
#   Tcells_clean_UMAP_clusters.pdf        -> new UMAP
# If any NEW cluster still shows Pax5/Cd79a/Sdc1/Jchain/Adgre1 -> add to
# `extra_drop_clusters` below and re-run (iterative cleanup).
# ==============================================================================

suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(harmony)
  library(ggplot2); library(dplyr)
})

set.seed(1234)
BASE_DIR <- "/QRISdata/Q8448/Mouse_disease_data/Tcells"
IN_RDS   <- file.path(BASE_DIR, "tcells_processed.rds")
OUT_RDS  <- file.path(BASE_DIR, "tcells_clean_reclustered.rds")

# Contaminating clusters confirmed by FindAllMarkers (by current cell_type label)
CONTAM_TYPES <- c("Tfh_like_T", "CD8_Eff", "NK", "Cytotoxic_CD8_T", "Memory_CD8_T")
# If a NEW cluster after re-clustering is still contamination, add its
# seurat_clusters id here and re-run (leave empty on first pass).
# Pass 2: 9 = high-mito low-quality, 10 = residual B cells, 12 = endothelial/stromal.
extra_drop_clusters <- c("9", "10", "12")

# Marker-score safety gate (remove individual straggler non-T cells)
USE_SAFETY_GATE <- TRUE
GATE_THRESH     <- 0.25   # module-score cutoff; raise to be more lenient

tcells <- readRDS(IN_RDS)
message("Loaded: ", ncol(tcells), " cells")
cat("\nOriginal cell_type counts:\n"); print(table(tcells$cell_type))

# ---- 1) Remove labelled contaminating clusters ------------------------------
tcells <- subset(tcells, subset = cell_type %in% CONTAM_TYPES, invert = TRUE)
message("\nAfter removing labelled contamination: ", ncol(tcells), " cells")

# ---- 2) Optional safety gate by contamination marker score ------------------
if (USE_SAFETY_GATE) {
  DefaultAssay(tcells) <- "SCT"
  contam_markers <- list(
    Bcell   = c("Cd79a","Cd79b","Ebf1","Pax5","Bank1","Cd19","Ms4a1","Cd180","Lyn"),
    Plasma  = c("Sdc1","Jchain","Derl3","Xbp1","Prdm1"),
    Myeloid = c("Adgre1","Adgre4","Pparg","Lyz2","Itgam","Csf1r"),
    Endoth  = c("Vwf","Pecam1","Cdh5","Hspg2","Adamts18")
  )
  for (nm in names(contam_markers)) {
    g <- contam_markers[[nm]][contam_markers[[nm]] %in% rownames(tcells)]
    if (length(g) >= 2)
      tcells <- AddModuleScore(tcells, features = list(g), name = paste0(nm, "_score"))
  }
  sc_cols <- grep("_score1$", colnames(tcells@meta.data), value = TRUE)
  flag <- rep(FALSE, ncol(tcells))
  for (sc in sc_cols) flag <- flag | (tcells@meta.data[[sc]] > GATE_THRESH)
  message("Residual contamination flagged by marker score (>", GATE_THRESH, "): ", sum(flag))
  keep_cells <- colnames(tcells)[!flag]
  tcells <- subset(tcells, cells = keep_cells)
  message("After safety gate: ", ncol(tcells), " cells")
}

# ---- 3) Re-run LSI / Harmony / UMAP / clustering (same params as original) ---
DefaultAssay(tcells) <- "ATAC"
message("\nRe-running LSI + Harmony + UMAP + clustering on clean T cells...")
tcells <- RunTFIDF(tcells)
tcells <- FindTopFeatures(tcells, min.cutoff = "q0")
tcells <- RunSVD(tcells)
tcells <- RunHarmony(
  object        = tcells,
  group.by.vars = "deMultliplex2_final_mapped",
  reduction.use = "lsi",
  assay.use     = "ATAC",
  project.dim   = FALSE,
  verbose       = FALSE
)
tcells <- RunUMAP(tcells, reduction = "harmony", dims = 2:30,
                  seed.use = 1234, reduction.name = "umap")
tcells <- FindNeighbors(tcells, reduction = "harmony", dims = 2:30)
tcells <- FindClusters(tcells, algorithm = 3, resolution = 0.6,
                       verbose = FALSE, random.seed = 1234)

# Second-pass: drop residual contaminated/low-quality clusters, then RE-EMBED
# (so the final UMAP/clusters are computed without them).
if (length(extra_drop_clusters) > 0) {
  tcells <- subset(tcells, subset = seurat_clusters %in% extra_drop_clusters, invert = TRUE)
  message("Dropped extra clusters [", paste(extra_drop_clusters, collapse=","), "]; ",
          ncol(tcells), " cells remain. Re-embedding...")
  DefaultAssay(tcells) <- "ATAC"
  tcells <- RunTFIDF(tcells)
  tcells <- FindTopFeatures(tcells, min.cutoff = "q0")
  tcells <- RunSVD(tcells)
  tcells <- RunHarmony(tcells, group.by.vars = "deMultliplex2_final_mapped",
                       reduction.use = "lsi", assay.use = "ATAC",
                       project.dim = FALSE, verbose = FALSE)
  tcells <- RunUMAP(tcells, reduction = "harmony", dims = 2:30,
                    seed.use = 1234, reduction.name = "umap")
  tcells <- FindNeighbors(tcells, reduction = "harmony", dims = 2:30)
  tcells <- FindClusters(tcells, algorithm = 3, resolution = 0.6,
                         verbose = FALSE, random.seed = 1234)
}

# ---- 4) UMAP preview --------------------------------------------------------
p1 <- DimPlot(tcells, reduction = "umap", label = TRUE, pt.size = 0.3) +
      ggtitle(paste0("Re-clustered clean T cells (n=", ncol(tcells), ")"))
ggsave(file.path(BASE_DIR, "Tcells_clean_UMAP_clusters.pdf"), p1, width = 8, height = 6)

# ---- 5) FindAllMarkers + canonical-marker DotPlot for re-annotation ---------
DefaultAssay(tcells) <- "SCT"
Idents(tcells) <- "seurat_clusters"
fam <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top <- fam %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15)
write.csv(top, file.path(BASE_DIR, "Tcells_clean_cluster_topmarkers.csv"), row.names = FALSE)

# Canonical T markers (Masopust 2026) + contamination re-check
recheck <- c("Cd3e","Cd3d","Cd8a","Cd8b1","Cd4",                 # lineage
             "Sell","Ccr7","Tcf7","Lef1","Il7r",                 # naive/memory
             "Klrg1","Tbx21","Gzmb","Ifng","Cx3cr1","Zeb2",      # effector
             "Foxp3","Ikzf2","Ctla4",                            # Treg
             "Cxcr5","Pdcd1","Tox","Havcr2",                     # Tfh/exhaustion
             "Cd79a","Pax5","Ebf1","Sdc1","Jchain","Adgre1","Pparg")  # contamination check
recheck <- recheck[recheck %in% rownames(tcells)]
pd <- DotPlot(tcells, features = recheck) + RotatedAxis() +
      theme(axis.text.x = element_text(size = 7))
ggsave(file.path(BASE_DIR, "Tcells_clean_recheck_dotplot.pdf"), pd, width = 14, height = 6)

saveRDS(tcells, OUT_RDS)
message("\nSaved cleaned re-clustered object: ", OUT_RDS)
message("Now inspect topmarkers.csv + recheck_dotplot.pdf, then annotate.")
