#!/usr/bin/env Rscript
# ==============================================================================
# Tcell_final_annotate.R  — FINAL T-cell annotation
#
# Input : tcells_clean_reclustered.rds  (11 clusters 0-10, contamination-cleaned)
# Steps : 1) drop cluster 10 (endothelial residual, 7 cells) — NO re-clustering
#         2) assign cell_type_final by canonical markers
#         3) save annotated object + final UMAP (Figure 7)
#
# Annotation basis (B approach: Masopust backbone + additional sources):
#   Naive_central_T   (g0,2,4)  Masopust 2026: Lef1/Sell/Ccr7/Tcf7
#   Activated_Cd69_T  (g1)      Masopust 2026: Cd69 (+Tbx21)
#   Effector_SLEC_CD8 (g3)      Masopust 2026: Klrg1/Cx3cr1/Tbx21/Gzmb/Ifng
#                               +Zeb2  -> Omilusik et al. 2015 J Exp Med (10.1084/jem.20150194)
#   Exhausted_Tex_CD8 (g5)      Masopust 2026: Tox/Havcr2(TIM3)/Pdcd1/Entpd1(CD39)
#                               +Maf,Nr4a2 -> Giordano et al. 2015 EMBO J (10.15252/embj.201490786)
#   Treg              (g6)      Masopust 2026: Foxp3/Il2ra/Ctla4/Ikzf2/Tnfrsf18  [pure Masopust]
#   TRM_CD8_CD103     (g7)      Masopust 2026: Itgae(CD103)/Cxcr6
#   TRM_CD8_CD49a     (g8)      Masopust 2026: Itga1(CD49a)/Cxcr6/Cd44  [Ly49 noted, not cited]
#   Tpex_Tfh_like     (g9)      Masopust 2026: Cxcr5/Pdcd1  (small; LCMV-chronic context)
#   [dropped]         (g10)     endothelial residual (Vwf)
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

BASE_DIR <- "/QRISdata/Q8448/Mouse_disease_data/Tcells"
IN_RDS   <- file.path(BASE_DIR, "tcells_clean_reclustered.rds")
OUT_RDS  <- file.path(BASE_DIR, "tcells_final_annotated.rds")

tc <- readRDS(IN_RDS)
message("Loaded: ", ncol(tc), " cells, ",
        length(unique(tc$seurat_clusters)), " clusters")

# ---- 1) Drop cluster 10 (endothelial residual) — keep existing UMAP/clusters ----
tc <- subset(tc, subset = seurat_clusters %in% c("10"), invert = TRUE)
message("After dropping cluster 10: ", ncol(tc), " cells")

# ---- 2) Assign final cell types ---------------------------------------------
anno <- c(
  "0" = "Naive_central_T",
  "1" = "Activated_Cd69_T",
  "2" = "Naive_central_T",
  "3" = "Effector_SLEC_CD8",
  "4" = "Naive_central_T",
  "5" = "Exhausted_Tex_CD8",
  "6" = "Treg",
  "7" = "TRM_CD8_CD103",
  "8" = "TRM_CD8_CD49a",
  "9" = "Tpex_Tfh_like"
)
# Safety check: every remaining cluster must be in the map
present <- as.character(sort(unique(as.integer(as.character(tc$seurat_clusters)))))
missing <- setdiff(present, names(anno))
if (length(missing) > 0) stop("Unmapped clusters: ", paste(missing, collapse=","))

Idents(tc) <- "seurat_clusters"
tc <- RenameIdents(tc, anno)
tc$cell_type_final <- as.character(Idents(tc))

# Order levels for a tidy legend (functional grouping)
lvl <- c("Naive_central_T","Activated_Cd69_T","Effector_SLEC_CD8",
         "Exhausted_Tex_CD8","Tpex_Tfh_like","Treg",
         "TRM_CD8_CD103","TRM_CD8_CD49a")
tc$cell_type_final <- factor(tc$cell_type_final, levels = lvl)
Idents(tc) <- "cell_type_final"

cat("\nFinal cell_type_final counts:\n")
print(table(tc$cell_type_final))

# ---- 3) Final annotated UMAP (Figure 7) -------------------------------------
pal <- c("Naive_central_T"="#4C72B0","Activated_Cd69_T"="#DD8452",
         "Effector_SLEC_CD8"="#C44E52","Exhausted_Tex_CD8"="#8172B3",
         "Tpex_Tfh_like"="#937860","Treg"="#DA8BC3",
         "TRM_CD8_CD103"="#55A868","TRM_CD8_CD49a"="#8C8C00")

p <- DimPlot(tc, reduction = "umap", group.by = "cell_type_final",
             cols = pal, label = TRUE, repel = TRUE, label.size = 4,
             pt.size = 0.3) +
     ggtitle(paste0("T cells — final annotation (n=", ncol(tc), ")")) +
     theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(BASE_DIR, "Tcells_final_annotated_UMAP.pdf"), p, width = 9, height = 7)
ggsave(file.path(BASE_DIR, "Tcells_final_annotated_UMAP.png"), p, width = 9, height = 7, dpi = 300)

saveRDS(tc, OUT_RDS)
message("\nSaved: ", OUT_RDS)
message("Figure: Tcells_final_annotated_UMAP.pdf / .png")
